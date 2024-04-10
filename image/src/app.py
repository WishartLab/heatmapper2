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
from scipy.interpolate import interp2d
from numpy import arange, meshgrid, linspace
from PIL import Image
from pathlib import Path

from shared import Cache, MainTab, NavBar, FileSelection, Filter, ColumnType, TableValueUpdate


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


	async def GenerateHeatmap():
		"""
		@brief Generates the heatmap, overlaying the Image with the DataFrame
		@returns The Plot's axis, for downloading purposes.
		"""

		with ui.Progress() as p:

			p.inc(message="Loading input...")
			df = await DataCache.Load(input)

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

			# Expand the data for more refined points
			x = arange(df.shape[1])
			y = arange(df.shape[0])
			x_new = linspace(0, df.shape[1] - 1, input.Smoothing())
			y_new = linspace(0, df.shape[0] - 1, input.Smoothing())
			interp_func = interp2d(x, y, df, kind=input.Interpolation().lower())
			data_interp = interp_func(x_new, y_new)
			X_new, Y_new = meshgrid(x_new, y_new)

			p.inc(message="Plotting...")
			fig, ax = subplots()

			# Add the image as an overlay, if we have one.
			if img is not None: ax.imshow(img, extent=[x_new.min(), x_new.max(), y_new.min(), y_new.max()], aspect="auto",zorder=0)

			im = ax.contourf(
				X_new, Y_new, data_interp,
				cmap=input.ColorMap().lower(),
				extent=[0, 1, 0, 1],
				zorder=1,
				alpha=input.Opacity(),
				algorithm=input.Algorithm().lower(),
				linestyles=input.Style().lower(),
				levels=input.Levels(),
			)

		# Visibility of features
		if "legend" in input.Features(): colorbar(im, ax=ax, label="Value")

		if "y" in input.Features(): ax.tick_params(axis="y", labelsize=input.TextSize())
		else: ax.set_yticklabels([])

		if "x" in input.Features(): ax.tick_params(axis="x", labelsize=input.TextSize())
		else: ax.set_xticklabels([])

		return ax


	@output
	@render.data_frame
	@reactive.event(input.SourceFile, input.File, input.Example, input.Update, input.Reset)
	async def LoadedTable(): return await DataCache.Load(input)


	@output
	@render.plot
	@reactive.event(input.SourceFile, input.File, input.Example, input.Image, input.Update, input.Reset, input.TextSize, input.Opacity, input.ColorMap, input.Algorithm, input.Interpolation, input.Style, input.Levels, input.Smoothing, input.Features)
	async def Heatmap(): return await GenerateHeatmap()


	@output
	@render.text
	@reactive.event(input.SourceFile, input.Example)
	def ExampleInfo(): return Info[input.Example()]["Description"]


	@render.download(filename="table.csv")
	async def DownloadTable(): yield (await DataCache.Load(input)).to_string()


	@reactive.Effect
	@reactive.event(input.Update)
	async def Update(): await DataCache.Update(input)


	@reactive.Effect
	@reactive.event(input.Reset)
	async def Reset(): await DataCache.Purge(input)


	@reactive.Effect
	@reactive.event(input.SourceFile, input.File, input.Example, input.TableRow, input.TableCol, input.Update, input.Reset)
	async def UpdateTableValue(): TableValueUpdate(await DataCache.Load(input), input)


app_ui = ui.page_fluid(

	NavBar("Image"),

	ui.layout_sidebar(
		ui.sidebar(

			FileSelection(examples={"example1.txt": "Example 1"}, types=[".csv", ".txt", ".dat", ".tsv", ".tab", ".xlsx", ".xls", ".odf"], project="Image"),

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
