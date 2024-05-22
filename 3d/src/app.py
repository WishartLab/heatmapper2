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
from pyvista import Plotter, plotting, read_texture, read as VistaRead

# Shared functions
from shared import Cache, MainTab, NavBar, FileSelection, Filter, ColumnType, TableOptions, InitializeConfig, ColorMaps, Update

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
		if suffix == ".obj": return VistaRead(path.resolve())
		if suffix == ".png" or suffix == ".jpg": return read_texture(path.resolve())
		else: return DataCache.DefaultHandler(path)
	DataCache = Cache("3d", DataHandler=HandleData)
	Data = reactive.value(None)
	Valid = reactive.value(False)
	Object = reactive.value(None)

	InitializeConfig(config, input)


	@reactive.effect
	@reactive.event(input.SourceFile, input.File, input.Example, input.Reset)
	async def UpdateData(): Data.set((await DataCache.Load(input))); Valid.set(False)

	@reactive.effect
	@reactive.event(input.SourceFile, input.Object, input.Example)
	async def UpdateObject():Object.set(await DataCache.Load(input,
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
			ui.modal_show(ui.modal("The provided input format cannot be rendered",
			title="Table cannot be rendered", easy_close=True, footer=None))


	@Table.set_patch_fn
	def UpdateTable(*, patch: render.CellPatch) -> render.CellValue:
		if input.Type() == "Integer": value = int(patch["value"])
		elif input.Type() == "Float": value = float(patch["value"])
		else: value = patch["value"]
		return value


	def GenerateHeatmap():
		with ui.Progress() as p:

			# Get the model and data. 
			p.inc(message="Loading input...")
			model = Object()
			source = GetData()

			# We support just rendering a model without data, but we need the model
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


	@output
	@render.ui
	def Heatmap(): return GenerateHeatmap()


	@output
	@render.ui
	@reactive.event(input.Update)
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


app_ui = ui.page_fluid(

	NavBar(),

	ui.layout_sidebar(
		ui.sidebar(

			FileSelection(examples={"example1.csv": "Example 1", "texture.jpg": "Example 2"}, types=[".csv", ".txt", ".dat", ".tsv", ".tab", ".xlsx", ".xls", ".odf", ".png", ".jpg"], project="3D"),

			Update(),

			TableOptions(config),

			ui.panel_conditional("input.SourceFile === 'Upload'", ui.input_file("Object", "Choose an Object File", accept=[".obj"], multiple=False)),

			config.Opacity.UI(ui.input_slider, id="Opacity", label="Heatmap Opacity", min=0.0, max=1.0, step=0.1),
			
			config.Colors.UI(ui.input_slider, id="Colors", label="Number of Colors", value=256, min=1, max=256, step=1),

			# Set the ColorMap used.
			config.ColorMap.UI(ui.input_select, id="ColorMap", label="Color Map", choices=ColorMaps),

			config.Style.UI(ui.input_select, id="Style", label="Style", choices=["Surface", "Wireframe", "Points"]),

			config.Features.UI(ui.input_checkbox_group, id="Features", label="Heatmap Features", choices=["Edges", "Lighting", "Interpolation", "Smooth Shading"]),
		),

		# Add the main interface tabs.
		MainTab(m_type=ui.output_ui),
	)
)

app = App(app_ui, server)
