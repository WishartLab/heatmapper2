#
# Heatmapper
# 3D
#
# This file contains the Shiny application for 3D Heatmapper.
# It can be run with the following command within this directory:
#		shiny run
#
# Exporting via ShinyLive is not currently supported, as pyvista
# is not yet available in the Pyodide environment. Required libraries
# include: openmpi, verdict, glew, alongside python libraries in requirements.txt
# WebGL is required for this application.
#

from numpy.core.multiarray import MAXDIMS
from shiny import App, reactive, render, ui
from pandas import DataFrame
from Bio.PDB import PDBParser, PDBIO
from io import StringIO
from numpy import mean
from numpy.linalg import norm

# Shared functions
from shared import Cache, MainTab, NavBar, FileSelection, Filter, ColumnType, TableOptions, InitializeConfig, ColorMaps, Update, Pyodide, Error, Msg, File

if not Pyodide: from pyvista import Plotter, plotting, read_texture, read as VistaRead

from py3Dmol import view

try:
	from user import config
except ImportError:
	from config import config


def server(input, output, session):

	# Information regarding example files.
	Info = {
		"example1.csv": {
			"Object": "bunny.obj",
			"Description": "Input type: csv, obj<br>Contents: A bunny, mapped with random data."
		},
		"texture.jpg": {
			"Object": "FinalBaseMesh.obj",
			"Description": "Input type: jpg, obj<br>Contents: A human model with a sample heatmap texture applied.<br>Source: https://free3d.com/3d-model/male-base-mesh-6682.html"
		},
		"4K8X.pdb": {
			"Object": None,
			"Description": "Input type: pdb<br>Contents: An example protein PDB.<br>Source: https://dash.plotly.com/dash-bio/molecule3dviewer"
		}
	}

	Schemes = ["spectrum", "b-factor", "b-factor (norm)", "RMSF", "rainbow", "ssPyMol", "ssJmol", "Jmol", "amino", "shapely", "nucleic", "chain", "rasmol", "default", "greenCarbon", "cyanCarbon", "magentaCarbon", "purpleCarbon", "whiteCarbon", "orangeCarbon", "yellowCarbon", "blueCarbon", "chainHetatm"]


	def HandleData(path, p=None):
		"""
		@brief A custom Data Handler for the Cache.
		@param n: The Path object to the file.
		@returns A data object from the cache.
		@info This Data Handler supports object files, and images as textures.
		"""

		suffix = path.suffix
		if suffix == ".obj": return VistaRead(path.resolve())
		if suffix == ".png" or suffix == ".jpg": return read_texture(path.resolve())
		else: return DataCache.DefaultHandler(path)
	DataCache = Cache("3d", DataHandler=HandleData)

	Data = reactive.value(None)
	Valid = reactive.value(False)
	Object = reactive.value(None)

	InitializeConfig(config, input)


	@reactive.effect
	@reactive.event(input.SourceFile, input.File, input.Example, input.Reset, input.ID)
	async def UpdateData():
		"""
		@brief Updates the Data variables when input changes.
		@info When any relevant reactive input changes, this function requests data from the Cache,
		and then invalidates the current data in the table.
		"""
		if input.SourceFile() == "ID":
			data = await DataCache.Download(f"https://files.rcsb.org/view/{input.ID()}.pdb")
			Data.set(data)
		else:
			Data.set((await DataCache.Load(input, default=None, p=ui.Progress(), wasm_blacklist=(".csv", ".txt", ".dat", ".tsv", ".tab", ".xlsx", ".xls", ".odf", ".png", ".jpg"))))
		Valid.set(False)
		DataCache.Invalidate(File(input))


	@reactive.effect
	@reactive.event(input.SourceFile, input.Object, input.Example)
	async def UpdateObject():
		"""
		@brief Updates the Object variable when input changes.
		@info When the Object's source changes, this function fetches the most up to date information from the Cache.
		@info This function will not do anything in WebAssembly.
		"""
		if type(GetData()) != str:
			Object.set(await DataCache.Load(input,
				source_file=input.Object(),
				example_file=Info[input.Example()]["Object"],
				default=None,
				p=ui.Progress(),
				p_name="object",
				wasm=False
			))


	def GetData(): return Table.data_view() if Valid() else Data()


	@output
	@render.data_frame
	def Table():
		df = Data()
		try:
			grid = render.DataGrid(Data(), editable=True)
			Valid.set(True)
			return grid
		except TypeError:
			Error("The provided input format cannot be rendered")


	@Table.set_patch_fn
	def UpdateTable(*, patch: render.CellPatch) -> render.CellValue:
		if input.Type() == "Integer": value = int(patch["value"])
		elif input.Type() == "Float": value = float(patch["value"])
		else: value = patch["value"]
		DataCache.Invalidate(File(input))
		return value


	def PDBViewer(source, p):
		"""
		@brief Returns an HTML string containing the Py3DMol.js viewer of the source.
		@param source: The source containing PDB data as a string.
		@param p: The progress bar.
		@returns An HTML string that should be wrapped with ui.HTML
		"""

		# Used for caching.
		global_inputs = [
			input.File() if input.SourceFile() == "Upload" else input.ID() if input.SourceFile() == "ID" else input.Example(),
			config.Size(),
			config.ColorScheme(),
			config.PStyle(),
			config.PFeatures(),
			config.Thickness(),
			config.Width(),
			config.Opacity(),
			config.Radius(),
			config.Scale(),
			config.SurfaceOpacity(),
			config.SurfaceType(),
			config.SurfaceScheme(),
			config.Model(),
		]

		if not DataCache.In(global_inputs):

			parser = PDBParser()
			structure = parser.get_structure("protein", StringIO(source))
			model = config.Model()

			def GenerateScheme(source, initial_scheme, function_declared=False, model=0):
				"""
				@brief Py3DMol has a color, colorscheme, and colorfunc attribute. This function puts the right one in
				without cluttering the interface with three different options.
				@param initial_scheme: The value of the scheme. Spectrum is a color, B-Color requires a function, all others
				use a colorscheme
				@param function_declared: Since there is a colorscheme for both the heatmap and the structure, we can accidentally
				redefine the same JavaScript function twice if they both use the same B-Color scheme. This avoids that.
				"""
				prop = "color"
				scheme = initial_scheme

				# B Color requires a custom function.
				if scheme == "b-factor" or scheme == "b-factor (norm)":

					# Initial weights
					darkblue, blue, lightblue, white, orange, red = 5, 10, 15, 20, 40, 50

					# If we're normalizing, get the average B-Factor, then assign that as white.
					if "norm" in scheme:

						# For caching.
						entry = [input.File() if input.SourceFile() == "Upload" else input.ID() if input.SourceFile() == "ID" else input.Example()]
						a = 0.0
						l = 0
						if not DataCache.In(entry):
							p.inc(message="Normalizing...")
							l = 0
							# Get all the Atoms from the string, split by spaces, and remove empty entries.
							for atom in [atom for atom in source.split("\n") if atom.startswith("ATOM")]:
								entries = list(filter(None, atom.split(" ")))

								# This should be B-Factor, but sometimes the B-Factor is absent, in which case its an element.
								if not entries[10].isalpha():
									a += float(entries[10])
									l += 1
							a /= l

							# White is the average
							DataCache.Store((a * 0.25, a * 0.50, a * 0.75, a, a * 1.25, a * 1.5), entry)
						darkblue, blue, lightblue, white, orange, red = DataCache.Get(entry)
						Msg(f"Using normalized blue/white/red cutoffs at {lightblue:.2f}/{white:.2f}/{orange:.2f}")
						scheme = "NormalizedScheme"
					else: scheme = "Scheme"

					# Declare the function.
					if not function_declared:
						viewer.startjs += f"""\n
							let {scheme} = function(atom) {{
								if (atom.b < {darkblue}) return "darkblue"
								if (atom.b < {blue}) return "blue"
								else if (atom.b < {lightblue}) return "lightblue"
								else if (atom.b < {white}) return "white"
								else if (atom.b < {orange}) return "orange"
								else if (atom.b < {red}) return "red"
								else return "darkred"
							}}\n"""
					prop = "colorfunc"
				elif scheme == "RMSF":

					if len(structure) == 1:
						Error("RMSF requires a PDB with more than one model to compute difference!")
						return source, prop, scheme


					# List of atom names of interest
					atom_names_of_interest = ["C", "CA", "N"]

					entry = [input.File() if input.SourceFile() == "Upload" else input.ID() if input.SourceFile() == "ID" else input.Example(), model]
					if not DataCache.In(entry):
						main_model = structure[model]
						for chain in main_model:
							for residue in chain:
								for atom in residue:
									if atom.get_id() in atom_names_of_interest:
										distances = []
										for model in structure:
											if model != main_model:
												try:
													corresponding_atom = model[chain.id][residue.id][atom.get_id()]
													distance = norm(atom.coord - corresponding_atom.coord)
													distances.append(distance)
												except KeyError: continue
										# Calculate mean distance
										if distances: atom.set_bfactor(mean(distances))
										output = StringIO()

						output = StringIO()
						io = PDBIO()
						io.set_structure(main_model)
						io.save(output)
						DataCache.Store(output.getvalue(), entry)
						output.close()

					scheme = "RMSD"
					source = DataCache.Get(entry)

					darkblue, blue, lightblue, white, orange, red = 0.5, 1.0, 1.5, 2, 3, 4
					if not function_declared:
						viewer.startjs += f"""\n
							let {scheme} = function(atom) {{
								if (atom.b == 0) return "grey"
								else if (atom.b < {darkblue}) return "darkblue"
								else if (atom.b < {blue}) return "blue"
								else if (atom.b < {lightblue}) return "lightblue"
								else if (atom.b < {white}) return "white"
								else if (atom.b < {orange}) return "orange"
								else if (atom.b < {red}) return "red"
								else return "darkred"
							}}\n"""
					prop = "colorfunc"

				elif scheme == "rainbow":
					max = len([atom for atom in source.split("\n") if atom.startswith("ATOM")])
					prop = "colorscheme"
					scheme = {"prop": "index", "gradient": "ROYGB", "min": 0, "max": max}
				elif scheme != "spectrum": prop = "colorscheme"
				return source, prop, scheme

			viewer = view(width=f"{input.Size()}vw", height=f"{input.Size()}vh")
			source, heatmap_property, heatname_name = GenerateScheme(source, config.ColorScheme(), model=model)
			source, surface_property, surface_name = GenerateScheme(source, config.SurfaceScheme(), function_declared=config.SurfaceScheme() == config.ColorScheme(), model=model)

			viewer.addModelsAsFrames(source)
			viewer.zoomTo()


			p.inc(message="Styling...")
			# handle: if both radius and scale are included, radius overrides scale
			if config.Radius() > 0:
				var = "radius"
				val = config.Radius()
			else:
				var = "scale"
				if config.Scale() > 0:
					val = config.Scale()
				else:
					val = 1
			viewer.setStyle({config.PStyle().lower(): {
				heatmap_property: heatname_name,
				"style": "trace" if "Trace" in config.PFeatures() else "rectangle",
				"thickness": config.Thickness(),
				"tubes": "Tubes" in config.PFeatures(),
				"width": config.Width(),
				"opacity": config.Opacity(),
				"dashedBonds": "Dashed Bonds" in config.PFeatures(),
				"showNonBonded": "Show Non-Bonded" in config.PFeatures(),
				"singleBonds": "Single Bonds" in config.PFeatures(),
				var: val,
				# "scale": config.Scale(),
				# "radius": config.Radius(),
			}})
			viewer.addSurface(config.SurfaceType(), {"opacity": config.SurfaceOpacity(), surface_property: surface_name})

			if heatmap_property == "colorfunc": viewer.startjs = viewer.startjs.replace(f'"{heatname_name}"', f'{heatname_name}')
			if surface_property == "colorfunc": viewer.startjs = viewer.startjs.replace(f'"{surface_name}"', f'{surface_name}')


			p.inc(message="Exporting...")
			DataCache.Store(viewer.write_html(), global_inputs)

		return DataCache.Get(global_inputs)


	def ModelViewer(source, p):
		"""
		@brief Generates an HTML string of the PyVista Model viewer.
		@param source: The data to be applied to the object.
		@param p: The progress bar.
		@returns An HTML string that should be wrapped with ui.HTML
		@info Object will also need to be defined.
		"""
		if Pyodide:
			Error(f"Cannot render objects in WebAssembly! Please use the Server version for this functionality.")
			return

		# For Caching.
		inputs = [
			File(input),
			config.Style(),
			config.Opacity(),
			config.Features(),
			config.ColorMap(),
			config.Colors(),
		]

		if not DataCache.In(inputs):
			model = Object()
			if model is None: return

			p.inc(message="Plotting...")
			pl = Plotter()

			style = config.Style().lower()
			opacity = config.Opacity()
			features = config.Features()
			cmap = config.ColorMap().lower()
			colors = config.Colors()

			# If there's no source, just render the model
			if source is None:
				pl.add_mesh(
					model,
					style=style,
					opacity=opacity,
					show_edges="Edges" in features,
					lighting="Lighting" in features,
					smooth_shading="Smooth Shading" in features,
				)

			# If are data source is a table, render it as a heatmap.
			elif type(source) is DataFrame:
				values = source[Filter(source.columns, ColumnType.Name)]
				pl.add_mesh(
					model,
					scalars=values,
					style=style,
					cmap=cmap,
					opacity=opacity,
					n_colors=colors,
					show_edges="Edges" in features,
					lighting="Lighting" in features,
					smooth_shading="Smooth Shading" in features,
				)

			# If we have a texture, map it.
			elif type(source) is plotting.texture.Texture:
				mesh = model.texture_map_to_plane()
				pl.add_mesh(mesh, texture=source)

			# Exporting as None returns the HTML as a file handle, which we read.
			p.inc(message="Exporting...")
			DataCache.Store(pl.export_html(filename=None).read(), inputs)
		return DataCache.Get(inputs)


	def GenerateHeatmap():
		"""
		@brief Generates the Heatmap based on input, returns HTML
		@returns HTML of the heatmap.
		"""
		with ui.Progress() as p:

			# Get the model and data.
			p.inc(message="Loading input...")

			source = GetData()
			if source is None: return

			if type(source) == str:
				return PDBViewer(source, p)
			return ModelViewer(source, p)


	@output
	@render.ui
	def HeatmapReactive(): return ui.HTML(GenerateHeatmap())


	@reactive.effect
	@reactive.event(input.ExampleInfo)
	def ExampleInfo():
		Msg(ui.HTML(Info[input.Example()]["Description"]))


	@render.download(filename="table.csv")
	def DownloadTable():
		df = GetData()
		if df is not None:
			yield df.to_string()


	@render.download(filename="heatmap.html")
	def DownloadHeatmap():
		html = GenerateHeatmap()
		if html is not None:
			yield html


	@output
	@render.ui
	def ConditionalElements():
		elements = []
		data = GetData()

		if data is None: return

		if type(data) == str or input.SourceFile() == "ID":
			elements += [
				ui.HTML("<b>Model</b>"),
				config.ColorScheme.UI(ui.input_select, id="ColorScheme", label="Color Scheme", choices=Schemes, tooltip=ui.HTML('Define the coloring of the model. The default option `spectrum` applies a reversed gradient based on residue number. Read about other options <a href="https://3dmol.org/doc/global.html#builtinColorSchemes"; target="_blank">here</a>.')),
				config.Opacity.UI(ui.input_numeric, id="Opacity", label="Model Opacity", min=0.0, max=1.0, step=0.1, tooltip="Specify the opacity of the <i>model</i>. 1.0 indicates full opacity, while lower values make the model more transparent."),
				ui.HTML("<b>Surface</b>"),
				config.SurfaceScheme.UI(ui.input_select, id="SurfaceScheme", label="Color Scheme", choices=Schemes, tooltip='Define the coloring of the surface drawn on top of the model. The default option `spectrum` applies a reversed gradient based on residue number. Read about other options <a href="https://3dmol.org/doc/global.html#builtinColorSchemes"; target="_blank">here</a>.'),
				config.SurfaceOpacity.UI(ui.input_numeric, id="SurfaceOpacity", label="Surface Opacity", min=0.0, max=1.0, step=0.1, tooltip="Specify the opacity of the <i>surface</i> drawn on top of the model. 1.0 indicates full opacity, while lower values make the surface more transparent."),
				ui.HTML("<b>Customization</b>"),
				config.Model.UI(ui.input_numeric, id="Model", label="Model #", min=0, tooltip="Select which model to use from the PDB file, if the PDB contains multiple models."),
				config.PStyle.UI(ui.input_select, id="PStyle", label="Style", choices=["Cartoon", "Stick", "Sphere", "Line", "Cross"], tooltip=ui.HTML("Specify the rendering style of the model. <br>Cartoon visualizes secondary structures as ribbons, cylinders, arrows, and lines. <br>Stick depicts atoms as colored nodes and bonds as sticks. <br>Sphere depicts atoms as spheres. <br>Line depicts atoms and bonds as lines. <br>Cross depicts atoms as crosses.")),
				config.SurfaceType.UI(ui.input_select, id="SurfaceType", label="Surface Type", choices=["VDW", "MS", "SAS", "SES"], tooltip="Specify a surface to draw on top of the model. To see the surface, change the Surface Opacity above to be greater than 0. <br>VDW: van der Waals surface, each atom is surrounded by a sphere whose size is proportional to the van der Waals radius of that atom. <br>MS: Molecular surface, the outer boundary of the molecule, accessible by a probe sphere rolling over the VDW surface. <br>SAS: Solvent accessible surface, the boundary traced by the <i>centre</i> of a probe sphere rolling over the VDW surface. <br>SES: Solvent exposed surface, the outer boundary of a molecule where a solvent can come into contact with the molecule. This excludes parts of the surface that are shielded by other atoms."),
				config.Thickness.UI(ui.input_numeric, id="Thickness", label="Thickness", min=0, max=10, step=0.1, tooltip="Cartoon style only - specify the thickness of the visualized components. Lower values make components thinner, while higher values (to a maximum of 10) make components thicker."),
				config.Width.UI(ui.input_numeric, id="Width", label="Width", min=0, max=10, step=0.1, tooltip="Cartoon style only - specify the width of the visualized components. Lower values make components narrower, while higher values (to a maximum of 10) make components wider."),
				config.Radius.UI(ui.input_numeric, id="Radius", label="Radius", min=0, max=5, step=0.05, tooltip="Stick, Sphere, or Cross style only - Specify a fixed radius (in Angstroms) for atoms in the <i>model</i>. Lower values make atoms smaller, while higher values (to a maximum of 5) make atoms larger. This value overrides 'Scale'."),
				config.Scale.UI(ui.input_numeric, id="Scale", label="Scale", min=0, max=10, step=1, tooltip="Stick, Sphere, or Cross style only - specify a scalar to modify the van der Waals radius of atoms in the <i>model</i>. If a 'Radius' is specified above, this value is ignored. Set 'Radius' to 0 to visualize 'Scale'."),
				config.Size.UI(ui.input_numeric, id="Size", label="View Size", min=1, max=100, step=1, tooltip="Change the size of the viewer in your browser."),
				ui.HTML("<b>Features</b>"),
				config.PFeatures.UI(ui.input_checkbox_group, make_inline=False, id="PFeatures", label=None, choices=["Dashed Bonds", "Show Non-Bonded", "Single Bonds", "Tubes", "Trace"], tooltip="Dashed Bonds (Stick Style) - draw bonds as dashed lines. <br>Show Non-Bonded (Stick Style) - display non-bonded atoms as spheres. <br>Single Bonds (Stick Style) - display all bonds as single bonds. <br>Tubes (Cartoon Style) - display alpha helices as simple cylinders. <br>Trace (Cartoon Style) - draw the model as a simple outline. This overrides the 'Tubes' feature."),
			]

		else:
			elements.append(ui.panel_conditional("input.SourceFile === 'Upload'", ui.input_file("Object", "Choose an Object File", accept=[".obj"], multiple=False)))
			if type(data) == DataFrame:
				elements += [
					ui.HTML("<b>Heatmap</b>"),
					config.Opacity.UI(ui.input_numeric, id="Opacity", label="Opacity", min=0.0, max=1.0, step=0.1, tooltip="Specify the opacity of the heatmap. 1.0 indicates full opacity, while lower values make the heatmap more transparent."),
					config.Style.UI(ui.input_select, id="Style", label="Style", choices=["Surface", "Wireframe", "Points"], tooltip="Specify how to render the model. Surface visualizes data as triangles on a surface. Wireframe displays a wireframe of the outer geometry. Points displays values as a collection of points."),
					ui.HTML("<b>Colors</b>"),
					config.Colors.UI(ui.input_numeric, id="Colors", label="Number", value=256, min=1, step=1, tooltip="Specify the number of colors to use to visualize data."),
					config.ColorMap.UI(ui.input_select, id="ColorMap", label="Map", choices=ColorMaps, tooltip="Specify a color scheme to use."),
					ui.HTML("<b>Features</b>"),
					config.Features.UI(ui.input_checkbox_group, make_inline=False, id="Features", label=None, choices=["Edges", "Lighting", "Smooth Shading",], tooltip="Edges - display the wireframe edges on top of the surface visualization. Does not apply to Wireframe or Points. <br>Lighting - visualize the model with an external light source. Lighting may affect color accuracy. Lighting must be enabled for smooth shading to be applied. <br>Smooth Shading - smooth shadows on the surface of the model."),
			]

		return elements



app_ui = ui.page_fluid(

	ui.tags.style("""
		.navbar {
			position: fixed;  /* prevent navbar from scrolling */
			top: 0;
			height: 10vh;
			width: 100%;
			z-index: 1001;
			overflow-x: auto;
        }
		.navbar-nav {
			flex-wrap: nowrap !important;
		}
			   
		.bslib-sidebar-layout {
			margin-top: 10vh;  /* prevent content from being hidden under navbar */
		}
		.bslib-grid {
			display: flex;
			width: 100%;
		    justify-content: space-between;
		}	   

		#MainTab {
			position: sticky;  /* prevent tabs from scrolling */
			top: 0;
			width: 100%;
			z-index: 1000;
			background: rgba(255, 255, 255, 0.25);
		}
	"""),

	ui.panel_title(title=None, window_title="3D"),
	NavBar(),

	ui.layout_sidebar(
		ui.sidebar(
			FileSelection(
				examples={"4K8X.pdb": "Example 1", "example1.csv": "Example 2", "texture.jpg": "Example 3"},
				types=[".csv", ".txt", ".dat", ".tsv", ".tab", ".xlsx", ".xls", ".odf", ".png", ".jpg", ".pdb"],
				project="3D",
				extras=["ID"]),

			ui.panel_conditional(
				"input.SourceFile === 'ID'",
				ui.input_text(id="ID", value="1upp", label="PDB ID"),
			),

			TableOptions(config),

			ui.panel_conditional(
				"input.MainTab === 'HeatmapTab'",

				ui.output_ui(id="ConditionalElements"),

				ui.download_button(id="DownloadHeatmap", label="Download HTML"),
			),
			padding="10px",
			gap="20px",
			width="300px",
		),

		# Add the main interface tabs.
		MainTab(m_type=ui.output_ui),
		height="86vh",
	)
)

app = App(app_ui, server)
