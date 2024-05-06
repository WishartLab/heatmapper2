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
from PIL import Image

from shared import Cache, MainTab, NavBar, FileSelection, Filter, ColumnType, TableOptions


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

	@reactive.effect
	@reactive.event(input.SourceFile, input.File, input.Example, input.Reset)
	async def UpdateData(): Data.set((await DataCache.Load(input))); Valid.set(False)


	def GetData(): return Table.data_view() if Valid() else Data()


	@output
	@render.data_frame
	def Table(): Valid.set(True); return render.DataGrid(Data(), editable=True)


	@Table.set_patch_fn
	def UpdateTable(*, patch: render.CellPatch) -> render.CellValue:
		if input.Type() == "Integer": value = int(patch["value"])
		elif input.Type() == "Float": value = float(patch["value"])
		else: value = patch["value"]
		return value


	@output
	@render.plot
	async def Heatmap():
		with ui.Progress() as p:

			p.inc(message="Loading input...")
			df = GetData()

			p.inc(message="Loading image...")
			img = await DataCache.Load(input, source_file=input.Image(), example_file=Info[input.Example()]["Image"])
			if img is None or df.empty: return None

			# Wrangle into an acceptable format.
			p.inc(message="Formatting...")
			v_col = Filter(df.columns, ColumnType.Value, only_one=True, reject_unknown=True)
			x_col = Filter(df.columns, ColumnType.X, only_one=True, bad = [v_col], reject_unknown=True)
			y_col = Filter(df.columns, ColumnType.Y, only_one=True, bad = [v_col, x_col], reject_unknown=True)

			if {v_col, x_col, y_col}.issubset(df.columns):
				df = df.pivot(index=x_col, columns=y_col, values=v_col).reset_index(drop=True)

			p.inc(message="Plotting...")
			fig, ax = subplots()

			# Add the image as an overlay, if we have one.
			if img is not None: ax.imshow(img, extent=[0, 1, 0, 1], aspect="auto",zorder=0)

			im = ax.contourf(
				df,
				cmap=input.ColorMap().lower(),
				extent=[0, 1, 0, 1],
				zorder=1,
				alpha=input.Opacity(),
				algorithm=input.Algorithm().lower(),
				linestyles=input.Style().lower(),
				levels=input.Levels(),
			)

		# Visibility of features
		if "legend" in input.Features(): 
			cbar = colorbar(im, ax=ax, label="Value")
			cbar.ax.tick_params(labelsize=input.TextSize())

		if "y" in input.Features(): ax.tick_params(axis="y", labelsize=input.TextSize())
		else: ax.set_yticklabels([])

		if "x" in input.Features(): ax.tick_params(axis="x", labelsize=input.TextSize())
		else: ax.set_xticklabels([])

		return ax


	@output
	@render.text
	@reactive.event(input.SourceFile, input.Example)
	def ExampleInfo(): return Info[input.Example()]["Description"]


	@render.download(filename="table.csv")
	def DownloadTable(): yield GetData().to_string()


app_ui = ui.page_fluid(

	NavBar(),

	ui.layout_sidebar(
		ui.sidebar(

			FileSelection(examples={"example1.txt": "Example 1"}, types=[".csv", ".txt", ".dat", ".tsv", ".tab", ".xlsx", ".xls", ".odf"], project="Image"),

			TableOptions(),

			ui.panel_conditional("input.SourceFile === 'Upload'", ui.input_file("Image", "Choose your Image File", 
				accept=[".bmp", ".gif", ".h5", ".hdf", ".ico", ".jpeg", ".jpg", ".tif", ".tiff", ".webp", ".png"], 
				multiple=False)),

			# Customize the text size of the axes.
			ui.input_numeric(id="TextSize", label="Text Size", value=8, min=1, max=50, step=1),

			# Customize the opacity of the heatmap, making the background image more visible.
			ui.input_slider(id="Opacity", label="Heatmap Opacity", value=0.5, min=0.0, max=1.0, step=0.1),

			# Set the ColorMap used.
			ui.input_select(id="ColorMap", label="Color Map", choices=["Viridis", "Plasma", "Inferno", "Magma", "Cividis"], selected="Viridis"),

			# https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.contourf.html#matplotlib.pyplot.contourf
			ui.input_select(id="Algorithm", label="Contour Algorithm", choices=["MPL2005", "MPL2014", "Serial", "Threaded"], selected="MPL2014"),

			# https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp2d.html
			ui.input_select(id="Interpolation", label="Spline Interpolation", choices=["Linear", "Cubic", "Quintic"], selected="Linear"),

			# https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.contourf.html#matplotlib.pyplot.contourf
			ui.input_select(id="Style", label=" Line Style", choices=["Solid", "Dashed", "Dashdot", "Dotted"], selected="Solid"),

			ui.input_slider(id="Levels", label="Number of Levels", value=20, min=1, max=100, step=1),
			ui.input_slider(id="Smoothing", label="Smoothing", value=25, min=10, max=100, step=1),


			# Customize what aspects of the heatmap are visible
			ui.input_checkbox_group(id="Features", label="Heatmap Features",
					choices={"x": "X Labels", "y": "Y Labels", "legend": "Legend"},
					selected=["legend"]),

			# Add the download buttons.
			ui.download_button("DownloadTable", "Download Table"),
		),

		# Add the main interface tabs.
		MainTab(),
	)
)

app = App(app_ui, server)
