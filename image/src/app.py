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

from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from matplotlib.pyplot import subplots, colorbar
from scipy.interpolate import interp2d
from pandas import DataFrame
from io import BytesIO
from numpy import arange, meshgrid, linspace
from PIL import Image

from shared import Cache, MainTab, NavBar, FileSelection, Filter, ColumnType, FillColumnSelection, TableValueUpdate


def server(input: Inputs, output: Outputs, session: Session):

	# Information regarding example files.
	Info = {
		"example1.txt": {
			"Image": "example1.jpg",
			"Description": "Hypothetical example illustrating data overlaid on a satellite image. Input data are count or magnitude values within the overlaid grid sections."
		}
	}

	DataCache = Cache("image")

	async def GenerateHeatmap():
		"""
		@brief Generates the heatmap, overlaying the Image with the DataFrame
		@returns The Plot's axis, for downloading purposes.
		"""

		df = await DataCache.Load(input)
		img = await DataCache.Load(input, source_file=input.Image(), example_file=Info[input.Example()]["Image"])
		if df is None or img is None: return None

		# Wrangle into an acceptable format.
		v_col = Filter(df.columns, ColumnType.Value, only_one=True)
		x_col = Filter(df.columns, ColumnType.X, only_one=True, bad = [v_col])
		y_col = Filter(df.columns, ColumnType.Y, only_one=True, bad = [v_col, x_col])

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

		if "y" in input.Features():
			ax.tick_params(axis="y", labelsize=input.TextSize())
		else:
			ax.set_yticklabels([])

		if "x" in input.Features():
			ax.tick_params(axis="x", labelsize=input.TextSize())
		else:
			ax.set_xticklabels([])

		return ax


	@output
	@render.data_frame
	@reactive.event(input.Update, input.Reset, input.Example, input.File, ignore_none=False, ignore_init=False)
	async def LoadedTable(): return await DataCache.Load(input)


	@output
	@render.plot
	@reactive.event(input.Update, input.Reset, input.Example, input.File, input.TextSize, input.Opacity, input.ColorMap, input.Algorithm, input.Style, input.Levels, input.Features, input.Smoothing, input.Interpolation, input.Image, ignore_none=False, ignore_init=False)
	async def Heatmap(): return await GenerateHeatmap()


	@output
	@render.text
	def ExampleInfo(): return Info[input.Example()]["Description"]


	@render.download(filename="table.csv")
	async def DownloadTable(): df = await DataCache.Load(input); yield df.to_string()


	@reactive.Effect
	@reactive.event(input.Update)
	async def Update(): await DataCache.Update(input)


	@reactive.Effect
	@reactive.event(input.Reset)
	async def Reset(): await DataCache.Purge(input)


	@reactive.Effect
	@reactive.event(input.TableRow, input.TableCol, input.Example, input.File, input.Reset, input.Update)
	async def UpdateTableValue(): TableValueUpdate(await DataCache.Load(input), input)


app_ui = ui.page_fluid(

	NavBar("Image"),

	ui.layout_sidebar(
		ui.sidebar(

			FileSelection(examples={"example1.txt": "Example 1"}, types=[".csv", ".txt", ".xlsx"]),

			ui.panel_conditional("input.SourceFile === 'Upload'", ui.input_file("Image", "Choose your Image File", accept=[".png", ".jpg"], multiple=False)),

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
