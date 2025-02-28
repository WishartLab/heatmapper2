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
from pandas import DataFrame, read_table
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
	#Schemes = ["spectrum", "b-factor", "b-factor (norm)", "RMSF", "RMSD", "ssJmol", "amino", "shapely", "nucleic", "chain", "rasmol"]
	Schemes = ["Residue #", "Reverse Residue #", "B-factor", "RMSF", "RMSD", "2ndary Structure", "pLDDT"]


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
		if input.SourceFile() == "PDB-ID":
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
		data = Data()
		
		# warning message if input is an image, not a table
		if data is None:
			return DataFrame({"Note": ["No data to display! Please upload your data or select an example data set in the sidebar."]})
		
		if isinstance(data, plotting.texture.Texture):
			df = DataFrame({"Note": ["This heatmap is mapping an image file (.png or .jpg) onto the 3D surface. There is no table data to display."]})
			return df
		# display residue numbers, B-factor data from PDB files
			'''
			Dr. Wishart notes:
			The Table should display the residue numbers (column 1) and the B-factor or RMSD or RMSF values (column 2)
			'''
		elif isinstance(data, str):
			output_data = []
			selection = config.ColorScheme()
			col_name = "B-factor"
			if selection == "RMSF":
				col_name = "RMSF"
			elif selection == "RMSD":
				col_name = "RMSD"
			elif selection == "pLDDT":
				col_name = "pLDDT"

			for line in data.splitlines():
				if line.startswith("ATOM"):
					residue_num = int(line[22:26].strip())
					b_factor = float(line[60:66].strip())
					output_data.append((residue_num, b_factor))
			df = DataFrame(output_data, columns=["Residue Number", col_name])
			return df
		# display table data
		else:
			try:
				grid = render.DataGrid(Data(), editable=True)
				Valid.set(True)
				print(grid)
				return grid
			except TypeError:
				#Error("Please ensure your uploaded file is properly formatted. The provided input format cannot be rendered.")
				return DataFrame({"Note": ["Please ensure your uploaded file is properly formatted. The provided input format cannot be rendered."]})


	@Table.set_patch_fn
	def UpdateTable(*, patch: render.CellPatch) -> render.CellValue:
		if input.Type() == "Integer": value = int(patch["value"])
		elif input.Type() == "Float": value = float(patch["value"])
		else: value = patch["value"]
		DataCache.Invalidate(File(input))
		return value
	

	# Info text in welcome tab
	@render.ui
	def Welcome():
		return ui.HTML("""
			<h1>3D Heatmaps</h1>
			3D heatmaps render PDB models, or display data or textures on a 3D model. 
			<br><br>
			Upload a data file or specify a PDB ID in the sidebar to get started, or select 'Example' to check out a pre-loaded example.
			<br><br>	 
			Navigate to the 'Heatmap' tab to see the heatmap, or 'Table' to look at the input data.
				 
			<br><br>
			<img src="https://github.com/WishartLab/heatmapper2/wiki/assets/3D.png" alt="Image"; style="max-width:500px;">
				 
			<br><br><h3>Format</h3>
			<i>3D heatmaps can be created in two different ways:</i><br><br>
			<b>1 - PDB</b><br>
			Upload a <b>.pdb</b> file, or select 'ID' in the sidebar and enter a <b>PDB ID</b> (see example PDB 4K8X).
			<br><br>
			<b>2 - Object Files</b><br>
			Input an .obj file and either a table file or an image. If an image is used, it will be mapped onto the surface of the object (see Example 3).<br>
			If a table file is used, values in a 'Value' column will be mapped to each face of the model. If a 'Name' column is provided, it should contain the numerical values of the faces, otherwise the values will be applied linearly (see Example 2).
			<br><br>
			<table style="border-spacing: 100px";>
			<tr>
				<th>Table Files</th>
				<th>Image Files</th>
			</tr>
			<tr>
				<td style="padding-right:75px;">
					<li>.csv</li>
					<li>.dat</li>
					<li>.odf</li>
					<li>.tab</li>
					<li>.tsv</li>
					<li>.txt</li>
					<li>.xls</li>
					<li>.xlsx</li>
				</td>
				<td style="vertical-align:top;">
					<li>.jpg</li>
					<li>.png</li>
				</td>
			</tr>
			</table> 
			<br><br>
			<br><h3>Interface</h3>
			Click on the '?' icon beside sidebar options to read more about them.
			
		""")


	def PDBViewer(source, p):
		"""
		@brief Returns an HTML string containing the Py3DMol.js viewer of the source.
		@param source: The source containing PDB data as a string.
		@param p: The progress bar.
		@returns An HTML string that should be wrapped with ui.HTML
		"""

		# Used for caching.
		global_inputs = [
			config.ModelType(),
			input.File() if input.SourceFile() == "Upload" else input.ID() if input.SourceFile() == "PDB-ID" else input.Example(),
			config.Size(),
			config.ColorScheme(),
			config.PFeatures(),
			config.Thickness(),
			config.Width(),
			config.Opacity(),
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
				# dict mapping 
				scheme_dict = {
					"Residue #": "residue",
					"Reverse Residue #": "reverse",
					"B-factor": "b-factor",
					"2ndary Structure": "ssJmol",
					"RMSF": "rmsf",
					"RMSD": "rmsd",
					"pLDDT": "plddt",
				}
				
				# B-factor, RMSD, RMSF, residue number, reverse residue number, secondary structure
				prop = "color"
				scheme = scheme_dict[initial_scheme]

				# B Color requires a custom function.
				#if scheme == "b-factor" or scheme == "b-factor (norm)":
				if scheme == "b-factor":

					# Initial weights
					darkblue, blue, lightblue, white, orange, red = 5, 10, 15, 20, 40, 50

					# If we're normalizing, get the average B-Factor, then assign that as white.
					# if "norm" in scheme:

					# 	# For caching.
					# 	entry = [input.File() if input.SourceFile() == "Upload" else input.ID() if input.SourceFile() == "PDB-ID" else input.Example()]
					# 	a = 0.0
					# 	l = 0
					# 	if not DataCache.In(entry):
					# 		p.inc(message="Normalizing...")
					# 		l = 0
					# 		# Get all the Atoms from the string, split by spaces, and remove empty entries.
					# 		for atom in [atom for atom in source.split("\n") if atom.startswith("ATOM")]:
					# 			entries = list(filter(None, atom.split(" ")))

					# 			# This should be B-Factor, but sometimes the B-Factor is absent, in which case its an element.
					# 			if not entries[10].isalpha():
					# 				a += float(entries[10])
					# 				l += 1
					# 		a /= l

					# 		# White is the average
					# 		DataCache.Store((a * 0.25, a * 0.50, a * 0.75, a, a * 1.25, a * 1.5), entry)
					# 	darkblue, blue, lightblue, white, orange, red = DataCache.Get(entry)
					# 	Msg(f"Using normalized blue/white/red cutoffs at {lightblue:.2f}/{white:.2f}/{orange:.2f}")
					# 	scheme = "NormalizedScheme"
					# else: scheme = "Scheme"
					scheme = "Scheme"

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
				
				# TODO: implement colour by pLDDT
				elif scheme == "plddt":
					pass
				
				# TODO: implement colour using RMSD data
				elif scheme == "rmsd":
					
					if len(structure) == 1:
						Error("RMSD requires a PDB with more than one model to compute difference!")
						return source, prop, scheme

				# colour using RMSF data
				elif scheme == "rmsf":

					if len(structure) == 1:
						Error("RMSF requires a PDB with more than one model to compute difference!")
						return source, prop, scheme

					# List of atom names of interest
					atom_names_of_interest = ["C", "CA", "N"]

					entry = [input.File() if input.SourceFile() == "Upload" else input.ID() if input.SourceFile() == "PDB-ID" else input.Example(), model]
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

					scheme = "RMSD"  # ?????
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

				# elif scheme == "rainbow":
				# 	max = len([atom for atom in source.split("\n") if atom.startswith("ATOM")])
				# 	prop = "colorscheme"
				# 	scheme = {"prop": "index", "gradient": "ROYGB", "min": 0, "max": max}
				
				elif scheme == "residue":
					# get number of residues (protein length)
					max = len([atom for atom in source.split("\n") if atom.startswith("ATOM")])
					prop = "colorscheme"
					scheme = {
						"prop": "index",  # index, b, resi
			   			"gradient": "linear", 
						"colors": ["red", "orange", "yellow", "green", "blue", "purple"],
						"min": 0, 
						"max": max
					}
				
				elif scheme == "reverse":
					# get number of residues (protein length)
					max = len([atom for atom in source.split("\n") if atom.startswith("ATOM")])
					prop = "colorscheme"
					scheme = {
						"prop": "index",  # index, b, resi
			   			"gradient": "linear", 
						"colors": ["purple", "blue", "green", "yellow", "orange", "red"],
						"min": 0, 
						"max": max
					}

				elif scheme == "ssJmol": 
					prop = "colorscheme"
					
				return source, prop, scheme

			viewer = view(width=f"{input.Size()}vw", height=f"{input.Size()}vh")
			source, heatmap_property, heatname_name = GenerateScheme(source, config.ColorScheme(), model=model)

			viewer.addModelsAsFrames(source)
			viewer.zoomTo()


			p.inc(message="Styling...")
			viewer.setStyle({"cartoon": {
				heatmap_property: heatname_name,
				"style": "trace" if "Simplified View" in config.PFeatures() else "rectangle",
				"thickness": config.Thickness(),
				"tubes": "Helices as Tubes" in config.PFeatures(),
				"width": config.Width(),
				"opacity": config.Opacity(),
				"scale": 1,
			}})

			if heatmap_property == "colorfunc": viewer.startjs = viewer.startjs.replace(f'"{heatname_name}"', f'{heatname_name}')


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
			# Error message pop-up
			Error(f"The WebAssembly version of Heatmapper2 does not support object rendering! Please use the Server version (server.heatmapper2.ca/3d) for this functionality.")
			# Error message in Heatmap tab
			return "The WebAssembly version of Heatmapper2 does not support object rendering! Please use the Server version (server.heatmapper2.ca/3d) for this functionality."

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
				try:
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
				except NotImplementedError:
					return "Make sure you have uploaded a Table or Image file, as well as an Object file!"
				except ValueError:
					return "The number of rows in the Table file must match either the number of cells in the Object file, or the number of points in the Object file!"

			# If we have a texture, map it.
			elif type(source) is plotting.texture.Texture:
				try:
					mesh = model.texture_map_to_plane()
					pl.add_mesh(mesh, texture=source)
				except:
					return "Make sure you have uploaded a Table or Image file, as well as an Object file!"

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
			if source is None: return "No data to display! <br>Please upload your data or select an example data set in the sidebar."

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
		# don't download if there is no data
		if df is None:
			Error("The downloaded table is empty! Please upload your data or select an example data set in the sidebar.")
		# if data is already a string
		elif isinstance(df, str):
			yield df
		else:
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

		if type(data) == str or input.SourceFile() == "PDB-ID":
			# add tooltip: two column *.csv file containing the protein residue numbers and the corresponding B-factor or RMSD or RMSF values
			elements.append(
				ui.panel_conditional("input.SourceFile === 'Upload'", ui.input_file("OptFile", "Add Optional B-factor, RMSD, or RMSF Data", accept=[".csv", ".txt", ".dat", ".tsv", ".tab", ".xlsx", ".xls", ".odf"], multiple=False)))
			elements += [
				ui.HTML("<b>Customization</b>"),
				config.ColorScheme.UI(ui.input_select, id="ColorScheme", label="Color Scheme", choices=Schemes, tooltip=ui.HTML('Define the coloring of the model. The default option `spectrum` applies a reversed gradient based on residue number. Read about other options <a href="https://3dmol.org/doc/global.html#builtinColorSchemes"; target="_blank">here</a>.')),			
				
				config.Model.UI(ui.input_numeric, id="Model", label="PDB Model", min=0, tooltip="Select which model to use from the PDB file, if the PDB contains multiple models."),
				
				config.Opacity.UI(ui.input_slider, id="Opacity", label="Model Opacity", min=0.0, max=1.0, step=0.1, tooltip=ui.HTML('Specify the opacity of the <i>model</i>. 1.0 indicates full opacity, while lower values make the model more transparent.')),	
				config.Thickness.UI(ui.input_slider, id="Thickness", label="Ribbon Thickness", min=0, max=10, step=0.1, tooltip="Specify the thickness of the visualized components. Lower values make components thinner, while higher values (to a maximum of 10) make components thicker."),
				config.Width.UI(ui.input_slider, id="Width", label="Ribbon Width", min=0, max=10, step=0.1, tooltip=ui.HTML("Specify the width of the visualized components. Lower values make components narrower, while higher values (to a maximum of 10) make components wider. <br>If 'Simplified View' is selected below, Ribbon Width is ignored.")),
				config.PFeatures.UI(ui.input_checkbox_group, make_inline=False, id="PFeatures", label=None, choices=["Helices as Tubes", "Simplified View"], tooltip=ui.HTML('''Helices as Tubes - display alpha helices as simple cylinders. <br><br>Simplified View - draw the model as a simple outline. This overrides the 'Helices as Tubes' feature.''')),

		
				config.Size.UI(ui.input_numeric, id="Size", label="View Size", min=1, max=100, step=1, tooltip="Change the size of the viewer in your browser."),
			]

		else:
			elements.append(ui.panel_conditional("input.SourceFile === 'Upload'", ui.input_file("Object", "Choose an Object File", accept=[".obj"], multiple=False)))
			if type(data) == DataFrame:
				elements += [
					ui.HTML("<b>Heatmap</b>"),
					config.Opacity.UI(ui.input_slider, id="Opacity", label="Opacity", min=0.0, max=1.0, step=0.1, tooltip="Specify the opacity of the heatmap. 1.0 indicates full opacity, while lower values make the heatmap more transparent."),
					config.Style.UI(ui.input_select, id="Style", label="Style", choices=["Surface", "Wireframe", "Points"], tooltip="Specify how to render the model. Surface visualizes data as triangles on a surface. Wireframe displays a wireframe of the outer geometry. Points displays values as a collection of points."),
					ui.HTML("<b>Colors</b>"),
					config.Colors.UI(ui.input_numeric, id="Colors", label="Number", value=256, min=1, step=1, tooltip="Specify the number of colors to use to visualize data."),
					config.ColorMap.UI(ui.input_select, id="ColorMap", label="Map", choices=ColorMaps, tooltip="Specify a color scheme to use."),
					ui.HTML("<b>Features</b>"),
					config.Features.UI(ui.input_checkbox_group, make_inline=False, id="Features", label=None, choices=["Edges", "Lighting", "Smooth Shading",], tooltip=ui.HTML("Edges - display the wireframe edges on top of the surface visualization. Does not apply to Wireframe or Points. <br>Lighting - visualize the model with an external light source. Lighting may affect color accuracy. Lighting must be enabled for smooth shading to be applied. <br>Smooth Shading - smooth shadows on the surface of the model.")),
			]

		return elements

	@output
	@render.ui
	def GetInputTypes():
		if config.ModelType() =="Object":
			return FileSelection(
				examples={"example1.csv": "Example 2", "texture.jpg": "Example 3"},
				types=[".csv", ".txt", ".dat", ".tsv", ".tab", ".xlsx", ".xls", ".odf", ".png", ".jpg"],
				project="3D",
			)
		else:
			return FileSelection(
				examples={"4K8X.pdb": "PDB 4K8X"},
				types=[".pdb"],
				project="3D",
				extras=["PDB-ID"],
			)


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
		
	"""),

	ui.panel_title(title=None, window_title="3D"),
	NavBar(),

	ui.layout_sidebar(
		ui.sidebar(
			ui.HTML("<br>"),
			#Update(),

			config.ModelType.UI(ui.input_select, id="ModelType", label="Choose 3D Model Format", choices=["Protein", "Object"], tooltip=ui.HTML('Create a Protein heat map using a .pdb file and optional Table file, or an Object heat map using an .obj file and an Image or Table file.')),
			
			ui.output_ui(id="GetInputTypes"),

			ui.panel_conditional(
				"input.SourceFile === 'PDB-ID'",
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
