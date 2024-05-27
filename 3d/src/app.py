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

from shiny import App, reactive, render, ui
from pandas import DataFrame
from os.path import basename

# Shared functions
from shared import Cache, MainTab, NavBar, FileSelection, Filter, ColumnType, TableOptions, InitializeConfig, ColorMaps, Update, Pyodide, Error

if not Pyodide:
	from pyvista import Plotter, plotting, read_texture, read as VistaRead

from py3Dmol import view
from Bio.PDB import PDBParser, Structure, PDBIO
from io import StringIO
from matplotlib.cm import get_cmap

try:
	from user import config
except ImportError:
	from config import config


def server(input, output, session):

	# Information regarding example files.
	Info = {
		"example1.csv": {
			"Object": "bunny.obj",
			"Description": "A bunny, mapped with random data."
		},
		"texture.jpg": {
			"Object": "FinalBaseMesh.obj",
			"Description": "A human model with a sample heatmap texture applied. Sourced from https://free3d.com/3d-model/male-base-mesh-6682.html"
		},
		"4K8X.pdb": {
			"Description": "An example protein PDB from dash-bio at https://dash.plotly.com/dash-bio/molecule3dviewer"
		}
	}


	def HandleData(path):
		"""
		@brief A custom Data Handler for the Cache.
		@param n: The Path object to the file.
		@returns A data object from the cache.
		@info This Data Handler supports object files, and images as textures.
		"""

		suffix = path.suffix
		if suffix == ".obj": return None if Pyodide else VistaRead(path.resolve())
		if suffix == ".png" or suffix == ".jpg": return None if Pyodide else read_texture(path.resolve())
		if suffix == ".pdb": return PDBParser().get_structure("protein", str(path))
		else: return DataCache.DefaultHandler(path)
	DataCache = Cache("3d", DataHandler=HandleData)

	Data = reactive.value(None)
	Valid = reactive.value(False)
	Object = reactive.value(None)

	InitializeConfig(config, input)


	@reactive.effect
	@reactive.event(input.SourceFile, input.File, input.Example, input.Reset, input.ID)
	async def UpdateData(): 
		if input.SourceFile() == "ID": 
			Data.set(input.ID())
		else: 
			Data.set((await DataCache.Load(input, default=None)))
		Valid.set(False)

	@reactive.effect
	@reactive.event(input.SourceFile, input.Object, input.Example)
	async def UpdateObject():
		example = input.Example()
		if not example.endswith(".pdb"):
			Object.set(await DataCache.Load(input,
				source_file=input.Object(),
				example_file=Info[input.Example()]["Object"],
				default=None
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
		return value


	def PDBViewer(source, p):

		# Store computations
		if input.SourceFile() == "ID":
			viewer = view(query=f"pdb:{source}", width=input.Size(), height=input.Size())
		else:
			inputs = [input.File() if input.SourceFile() == "Upload" else input.Example()]
			if not DataCache.In(inputs):
				pio = PDBIO()
				pio.set_structure(source)

				pdb_data = StringIO()
				pio.save(pdb_data)
				DataCache.Store(pdb_data.getvalue(), inputs)
			data = DataCache.Get(inputs)
			viewer = view(data=data, width=input.Size(), height=input.Size())

		if config.HeatmapType() == "Flexibility" and input.SourceFile() != "ID":
			atoms = [list(filter(None, entry.split(" "))) for entry in data.split("\n") if entry.startswith("ATOM")]
			b_factors = [(atom[1], float(atom[10])) for atom in atoms if not atom[10].isalpha()]
			m, M = min(b_factors, key=lambda item: item[1])[1], max(b_factors, key=lambda item: item[1])[1]

			cmap = get_cmap("RdBu_r")

			for serial, b_factor in b_factors:
				color = cmap((b_factor - m) / (M - m))
				hex_color = "#{:02x}{:02x}{:02x}".format(int(color[0] * 255), int(color[1] * 255), int(color[2] * 255))
				viewer.setStyle({"serial": serial}, {config.PStyle().lower(): {
					"color": hex_color,
					"arrows": "Arrows" in config.PFeatures(), 
					"style": config.PCStyle().lower(),
					"thickness": config.Thickness(),
					"tubes": "Tubes" in config.PFeatures(),
					"width": config.Width(),
					"opacity": config.Opacity(),
					"dashedBonds": "Dashed Bonds" in config.PFeatures(),
					"showNonBonded": "Show Non-Bonded" in config.PFeatures(),
					"singleBonds": "Single Bonds" in config.PFeatures(),
					"radius": config.Radius(),
					"scale": config.Scale(),
				}})

		else:
			if config.ColorScheme() == "spectrum":
				color_property = "color"
				color_name = config.ColorScheme()
			else: 
				color_property = "colorscheme"
				color_name = config.ColorScheme()

			color = config.ColorScheme() == "spectrum"
			viewer.setStyle({config.PStyle().lower(): {
				color_property: color_name,
				"arrows": "Arrows" in config.PFeatures(), 
				"style": config.PCStyle().lower(),
				"thickness": config.Thickness(),
				"tubes": "Tubes" in config.PFeatures(),
				"width": config.Width(),
				"opacity": config.Opacity(),
				"dashedBonds": "Dashed Bonds" in config.PFeatures(),
				"showNonBonded": "Show Non-Bonded" in config.PFeatures(),
				"singleBonds": "Single Bonds" in config.PFeatures(),
				"radius": config.Radius(),
				"scale": config.Scale(),
			}})


		return ui.HTML(viewer.write_html())


	def ModelViewer(source, p):
		if Pyodide:
			Error("The Object Viewer is not supported in WebAssembly Mode!")
			return

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
				interpolate_before_map="Interpolation" in features,
				smooth_shading="Smooth Shading" in features,
			)

		# If are data source is a table, render it as a heatmap.
		elif type(source) is DataFrame:
			values = source[Filter(source.columns, ColumnType.Name, only_one=True)]
			pl.add_mesh(
				model,
				scalars=values,
				style=style,
				cmap=cmap,
				opacity=opacity,
				n_colors=colors,
				show_edges="Edges" in features,
				lighting="Lighting" in features,
				interpolate_before_map="Interpolation" in features,
				smooth_shading="Smooth Shading" in features,
			)

		# If we have a texture, map it.
		elif type(source) is plotting.texture.Texture:
			mesh = model.texture_map_to_plane()
			pl.add_mesh(mesh, texture=source)

		# Exporting as None returns the HTML as a file handle, which we read.
		p.inc(message="Exporting...")
		return ui.HTML(pl.export_html(filename=None).read())


	def GenerateHeatmap():
		with ui.Progress() as p:

			# Get the model and data. 
			p.inc(message="Loading input...")

			source = GetData()

			if type(source) == Structure.Structure or input.SourceFile() == "ID":
				return PDBViewer(source, p)
			elif source is None: return
			else: return ModelViewer(source, p)


	@output
	@render.ui
	def HeatmapReactive(): return GenerateHeatmap()


	@output
	@render.text
	@reactive.event(input.SourceFile, input.Example)
	def ExampleInfo(): return Info[input.Example()]["Description"]


	@render.download(filename="table.csv")
	def DownloadTable():
		df = GetData()
		if df is not None:
			yield df.to_string()


	@output
	@render.ui
	def ConditionalElements():
		elements = []
		data = GetData()

		if data is None: return
		elements.append(config.Opacity.UI(ui.input_slider, id="Opacity", label="Heatmap Opacity", min=0.0, max=1.0, step=0.1))

		if type(data) == Structure.Structure or input.SourceFile() == "ID":
			if input.SourceFile() != "ID":
				elements.append(config.HeatmapType.UI(ui.input_radio_buttons, id="HeatmapType", label="Heatmap Type", choices=["Builtin", "Flexibility"], inline=True))
			elements += [
				config.PStyle.UI(ui.input_select, id="PStyle", label="Style", choices=["Cartoon", "Stick", "Sphere", "Line", "Cross"]),
				config.ColorScheme.UI(ui.input_select, id="ColorScheme", label="Color Scheme", choices=["spectrum", "ssPyMol", "ssJmol", "Jmol", "amino", "shapely", "nucleic", "chain", "rasmol", "default", "greenCarbon", "cyanCarbon", "magentaCarbon", "purpleCarbon", "whiteCarbon", "orangeCarbon", "yellowCarbon", "blueCarbon", "chainHetatm"]),
				config.PCStyle.UI(ui.input_select, id="PCStyle", label="Cartoon Style", choices=["Trace", "Oval", "Rectangle", "Parabola", "Edged"]),
				config.Thickness.UI(ui.input_numeric, id="Thickness", label="Strand Thickness", min=0, max=10, step=0.1),
				config.Width.UI(ui.input_numeric, id="Width", label="Strand Width", min=0, max=10, step=1),
				config.Radius.UI(ui.input_numeric, id="Radius", label="Radius", min=0, max=5, step=0.05),
				config.Scale.UI(ui.input_numeric, id="Scale", label="Scale", min=0, max=10, step=1),
				config.Size.UI(ui.input_numeric, id="Size", label="Viewer Size", min=1, step=1),
				config.PFeatures.UI(ui.input_checkbox_group, id="PFeatures", label="PDB Features", choices=["Dashed Bonds", "Show Non-Bonded", "Single Bonds", "Tubes"]),
			]	

		else:
			elements.append(ui.panel_conditional("input.SourceFile === 'Upload'", ui.input_file("Object", "Choose an Object File", accept=[".obj"], multiple=False)))
			if type(data) == DataFrame:
				elements += [
					config.Colors.UI(ui.input_slider, id="Colors", label="Number of Colors", value=256, min=1, max=256, step=1),
					config.ColorMap.UI(ui.input_select, id="ColorMap", label="Color Map", choices=ColorMaps),
					config.Style.UI(ui.input_select, id="Style", label="Style", choices=["Surface", "Wireframe", "Points"]),
					config.Features.UI(ui.input_checkbox_group, id="Features", label="Heatmap Features", choices=["Edges", "Lighting", "Interpolation", "Smooth Shading"]),
			]

		return elements



app_ui = ui.page_fluid(

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
				ui.input_text(id="ID", label="PDB ID", value="1ubq"),
			),

			TableOptions(config),

			ui.panel_conditional(
				"input.MainTab === 'HeatmapTab'",

				ui.output_ui(id="ConditionalElements")
			)
		),

		# Add the main interface tabs.
		MainTab(m_type=ui.output_ui),
		height="90vh",
	)
)

app = App(app_ui, server)
