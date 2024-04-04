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
from numpy import zeros, array, linspace
from matplotlib.pyplot import get_cmap
from matplotlib.colors import Normalize
from scipy.interpolate import interp1d
from io import BytesIO
from pathlib import Path
from tempfile import NamedTemporaryFile
from pandas import DataFrame

from pyvista import Plotter, MultiBlock, PolyData, plotting, read_texture, read as VistaRead
# Requires trame-vtk, trame-vuetify, nest_asyncio

# Shared functions
from shared import Cache, MainTab, NavBar, FileSelection, Filter, ColumnType, TableValueUpdate


def server(input, output, session):

	# Information regarding example files.
	Info = {
		"example1.csv": {
			"Object": "bunny.obj",
			"Description": "A bunny, mapped with random data."
		}
	}


	def HandleData(n, i):
		"""
		@brief A custom Data Handler for the Cache.
		@param n: The name of the file
		@param i: The source of the file. It can be a path to a file (string) or a BytesIO object.
		@returns A data object from the cache.
		@info This Data Handler supports object files, and images as textures.
		"""

		suffix = Path(n).suffix
		if suffix == ".obj":
				# If the file is sourced from a BytesIO, we need to store it in a file temporarily.
				if type(i) is BytesIO:
					temp = NamedTemporaryFile(suffix=suffix); temp.write(i.read())
					i = temp.name
				return VistaRead(i)
		if suffix == ".png" or suffix == ".jpg":
				# If the file is sourced from a BytesIO, we need to store it in a file temporarily.
				if type(i) is BytesIO:
					temp = NamedTemporaryFile(suffix=suffix); temp.write(i.read())
					i = temp.name
				return read_texture(i)
		else: return DataCache.DefaultHandler(n, i)
	DataCache = Cache("3d", DataHandler=HandleData)


	@output
	@render.data_frame
	@reactive.event(input.SourceFile, input.File, input.Example, input.File, input.Update, input.Reset)
	async def LoadedTable(): return await DataCache.Load(input)


	@output
	@render.ui
	@reactive.event(input.SourceFile, input.File, input.Example, input.Object, input.Update, input.Reset, input.ColorMap, input.Style, input.Colors, input.Opacity, input.Features)
	async def Heatmap():

		with ui.Progress() as p:

			p.inc(message="Loading input...")
			model = await DataCache.Load(
				input,
				source_file=input.Object(),
				example_file=Info[input.Example()]["Object"],
				default=None
			)
			source = await DataCache.Load(input, default=None)
			if model is None: return

			p.inc(message="Plotting...")
			pl = Plotter()

			if source is None:
				pl.add_mesh(
					model,
					style=input.Style(),
					opacity=input.Opacity(),
					show_edges="Edges" in input.Features(),
					lighting="Lighting" in input.Features(),
					interpolate_before_map="Interpolation" in input.Features(),
					smooth_shading="Smooth Shading" in input.Features(),
				)

			elif type(source) is DataFrame:
				values = source[Filter(source.columns, ColumnType.Name, only_one=True)]
				pl.add_mesh(
					model,
					scalars=values,
					style=input.Style(),
					cmap=input.ColorMap().lower(),
					opacity=input.Opacity(),
					n_colors=input.Colors(),
					show_edges="Edges" in input.Features(),
					lighting="Lighting" in input.Features(),
					interpolate_before_map="Interpolation" in input.Features(),
					smooth_shading="Smooth Shading" in input.Features(),
				)

			elif type(source) is plotting.texture.Texture:
				mesh = model.texture_map_to_plane()
				pl.add_mesh(mesh, texture=source)

			p.inc(message="Exporting...")
			return ui.HTML(pl.export_html(filename=None).read())


	@output
	@render.text
	@reactive.event(input.SourceFile, input.Example)
	def ExampleInfo(): return Info[input.Example()]["Description"]


	@render.download(filename="table.csv")
	async def DownloadTable():
		df = await DataCache.Load(input);
		if df is not None:
			yield df.to_string()


	@reactive.Effect
	@reactive.event(input.Update)
	async def Update(): await DataCache.Update(input)


	@reactive.Effect
	@reactive.event(input.Reset)
	async def Reset(): await DataCache.Purge(input)


app_ui = ui.page_fluid(

	NavBar("3D"),

	ui.layout_sidebar(
		ui.sidebar(

			FileSelection(examples={"example1.csv": "Example 1"}, types=[".csv", ".txt", ".dat", ".tsv", ".tab", ".xlsx", ".xls", ".odf", ".png", ".jpg"]),

			ui.panel_conditional("input.SourceFile === 'Upload'", ui.input_file("Object", "Choose an Object File", accept=[".obj"], multiple=False)),

			ui.input_slider(id="Opacity", label="Heatmap Opacity", value=1.0, min=0.0, max=1.0, step=0.1),
			ui.input_slider(id="Colors", label="Number of Colors", value=256, min=1, max=256, step=1),

			# Set the ColorMap used.
			ui.input_select(id="ColorMap", label="Color Map", choices=["Viridis", "Plasma", "Inferno", "Magma", "Cividis"], selected="Viridis"),

			ui.input_select(id="Style", label="Style", choices={"surface": "Surface", "wireframe": "Wireframe", "points": "Points", "points_gaussian": "Gaussian"}, selected="surface"),

			ui.input_checkbox_group(id="Features", label="Heatmap Features",
				choices=["Edges", "Lighting", "Interpolation", "Smooth Shading"],
				selected=["Lighting", "Interpolation", "Smooth Shading"]
			),

			# Add the download buttons.
			ui.download_button("DownloadTable", "Download Table"),
		),

		# Add the main interface tabs.
		MainTab(m_type=ui.output_ui),
	)
)

app = App(app_ui, server)
