#
# Heatmapper
# Image
#
# This file contains the ShinyLive application for Image Heatmapper.
# It can be run with the following command within this directory:
#		shinylive export . [site]
# Where [site] is the destination of the site folder.
#
# If you would rather deploy the application as a PyShiny application,
# run the following command within this directory:
#		shiny run
#
#

from shiny import App, reactive, render, ui
from matplotlib.pyplot import subplots, colorbar
from matplotlib.tri import Triangulation
from matplotlib.cm import get_cmap
from PIL import Image
from tempfile import NamedTemporaryFile
from io import BytesIO

from shared import Cache, MainTab, NavBar, FileSelection, Filter, ColumnType, TableOptions, InitializeConfig, ColorMaps, Update, Msg

try:
	from user import config
except ImportError:
	from config import config


def server(input, output, session):

	# Information regarding example files.
	Info = {
		"example1.txt": {
			"Image": "example1.jpg",
			"Description": "Hypothetical example illustrating data overlaid on a satellite image. Input data are count or magnitude values within the overlaid grid sections."
		}
	}

	def HandleData(path):
		"""
		@brief A custom Data Handler for the Cache.
		@param path: The path to the file
		@returns A data object from the cache.
		@info This Data Handler supports png and jpg images as PIL.Image objects
		"""
		if path.suffix in [".bmp", ".gif", ".h5", ".hdf", ".ico", ".jpeg", ".jpg", ".tif", ".tiff", ".webp", ".png"]: return Image.open(path.resolve())
		else: return DataCache.DefaultHandler(path)
	DataCache = Cache("image", DataHandler=HandleData)
	Data = reactive.value(None)
	Valid = reactive.value(False)
	IMG = reactive.value(None)

	InitializeConfig(config, input)


	@reactive.effect
	@reactive.event(input.SourceFile, input.File, input.Example, input.Reset)
	async def UpdateData(): Data.set((await DataCache.Load(input, p=ui.Progress()))); Valid.set(False)


	@reactive.effect
	@reactive.event(input.SourceFile, input.Example, input.Image)
	async def UpdateIMG(): IMG.set(await DataCache.Load(input,
			source_file=input.Image(),
			example_file=Info[input.Example()]["Image"]
		))


	def GetData(): return Table.data_view() if Valid() else Data()


	@output
	@render.data_frame
	def Table(): Valid.set(True); return render.DataGrid(Data(), editable=True)


	@Table.set_patch_fn
	def UpdateTable(*, patch: render.CellPatch) -> render.CellValue:
		if config.Type() == "Integer": value = int(patch["value"])
		elif config.Type() == "Float": value = float(patch["value"])
		else: value = patch["value"]
		return value


	def GenerateHeatmap():
		inputs = [
			input.File() if input.SourceFile() == "Upload" else input.Example(),
			input.Image(),
			config.ColorMap(),
			config.Opacity(),
			config.Algorithm(),
			config.Levels(),
			config.Features(),
			config.TextSize(),
			config.DPI(),
		]

		if not DataCache.In(inputs):
			with ui.Progress() as p:

				p.inc(message="Loading input...")
				df = GetData()

				# We need to reflect the DataFrame so that the first row is plotted at the top
				df = df.iloc[::-1]

				p.inc(message="Loading image...")
				img = IMG()
				if img is None or df.empty: return None

				# Wrangle into an acceptable format.
				p.inc(message="Formatting...")
				v_col = Filter(df.columns, ColumnType.Value)
				x_col = Filter(df.columns, ColumnType.X)
				y_col = Filter(df.columns, ColumnType.Y)

				if v_col != x_col and x_col != y_col and {v_col, x_col, y_col}.issubset(df.columns):
					df = df.pivot(index=x_col, columns=y_col, values=v_col).reset_index(drop=True)

				p.inc(message="Plotting...")
				fig, ax = subplots()

				# Add the image as an overlay, if we have one.
				if img is not None: ax.imshow(img, extent=[0, 1, 0, 1], aspect="auto",zorder=0)

				cmap = config.ColorMap().lower()
				alpha = config.Opacity()
				algorithm = config.Algorithm().lower()
				levels = config.Levels()

				# Make 0 values transparent.
				df = df.replace(0, -1)
				cmap = get_cmap(cmap)
				cmap.set_under(alpha=0)

				im = ax.contourf(
					df,
					cmap=cmap,
					extent=[0, 1, 0, 1],
					zorder=1,
					alpha=alpha,
					algorithm=algorithm,
					levels=levels,
				)

				# Visibility of features
				if "legend" in config.Features():
					cbar = colorbar(im, ax=ax, label="Value")
					cbar.ax.tick_params(labelsize=config.TextSize())

				if "y" in config.Features(): ax.tick_params(axis="y", labelsize=config.TextSize())
				else: ax.set_yticklabels([])

				if "x" in config.Features(): ax.tick_params(axis="x", labelsize=config.TextSize())
				else: ax.set_xticklabels([])

				b = BytesIO()
				fig.savefig(b, format="png", dpi=config.DPI())
				b.seek(0)
				DataCache.Store(b.read(), inputs)

		b = DataCache.Get(inputs)
		with NamedTemporaryFile(delete=False, suffix=".png") as temp:
			temp.write(b)
			temp.close()
			img: types.ImgData = {"src": temp.name, "height": f"{config.Size()}vh"}
			return img


	@output
	@render.image(delete_file=True)
	def Heatmap(): return GenerateHeatmap()


	@output
	@render.image(delete_file=True)
	@reactive.event(input.Update)
	def HeatmapReactive(): return GenerateHeatmap()


	@reactive.effect
	@reactive.event(input.ExampleInfo)
	def ExampleInfo():
		Msg(ui.HTML(Info[input.Example()]["Description"]))


	@render.download(filename="table.csv")
	def DownloadTable(): yield GetData().to_string()

	@render.download(filename="heatmap.png")
	def DownloadHeatmap():
		yield DataCache.Get([
			input.File() if input.SourceFile() == "Upload" else input.Example(),
			input.Image(),
			config.ColorMap(),
			config.Opacity(),
			config.Algorithm(),
			config.Levels(),
			config.Features(),
			config.TextSize(),
			config.DPI(),
		])


app_ui = ui.page_fluid(

	NavBar(),

	ui.layout_sidebar(
		ui.sidebar(

			FileSelection(examples={"example1.txt": "Example 1"}, types=[".csv", ".txt", ".dat", ".tsv", ".tab", ".xlsx", ".xls", ".odf"], project="Image"),

			ui.panel_conditional("input.SourceFile === 'Upload'", ui.input_file("Image", "Choose your Image File",
				accept=[".bmp", ".gif", ".h5", ".hdf", ".ico", ".jpeg", ".jpg", ".tif", ".tiff", ".webp", ".png"],
				multiple=False)),

			TableOptions(config),

			ui.panel_conditional(
				"input.MainTab === 'HeatmapTab'",

				Update(),

				ui.HTML("<b>Heatmap</b>"),
				config.TextSize.UI(ui.input_numeric, id="TextSize", label="Text", min=1, max=50, step=1),
				config.ColorMap.UI(ui.input_select, id="ColorMap", label="Map", choices=ColorMaps),
				config.Algorithm.UI(ui.input_select, id="Algorithm", label="Contour", choices=["MPL2005", "MPL2014", "Serial", "Threaded"]),
				config.Levels.UI(ui.input_numeric, id="Levels", label="Levels", min=1, step=1),
				config.Opacity.UI(ui.input_numeric, id="Opacity", label="Opacity", min=0.0, max=1.0, step=0.1),

				ui.HTML("<b>Image Settings</b>"),
				config.Size.UI(ui.input_numeric, id="Size", label="Size", value=800, min=1),
				config.DPI.UI(ui.input_numeric, id="DPI", label="DPI", value=100, min=1),

				# Customize what aspects of the heatmap are visible
				ui.HTML("<b>Features</b>"),
				config.Features.UI(ui.input_checkbox_group, make_inline=False, id="Features", label=None,
						choices={"x": "X Labels", "y": "Y Labels", "legend": "Legend"}
				),

				ui.download_button(id="DownloadHeatmap", label="Download"),
			),
			padding="10px",
			gap="20px",
			width="250px",
		),

		# Add the main interface tabs.
		MainTab(),
	)
)

app = App(app_ui, server)
