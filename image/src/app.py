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
from pandas import DataFrame
from io import BytesIO
from PIL import Image

from shared import Table, Cache, NavBar, FileSelection


def server(input: Inputs, output: Outputs, session: Session):

	# Information regarding example files.
	Info = {
		"example1.txt": {
			"Image": "example1.jpg",
			"Description": "Hypothetical example illustrating data overlaid on a satellite image. Input data are count or magnitude values within the overlaid grid sections."
		}
	}

	DataCache = Cache("image")

	async def LoadImage():
		"""
		@brief Loads the image to render behind the heatmap.
		@returns an Image object, if an image is specified, otherwise None.
		"""

		# Grab an uploaded file, if its done, or grab an example (Using a cache to prevent redownload)
		if input.SourceFile() == "Upload":
			file: list[FileInfo] | None = input.Image()
			return None if file is None else Image.open(file[0]["datapath"])
		else:
			n = Info[input.Example()]["Image"]
			cache = DataCache.Cache()
			if n not in cache:
				cache[n] = Image.open(BytesIO(await DataCache.Download(DataCache.Source + n)))
			return cache[n]


	async def GenerateHeatmap():
		"""
		@brief Generates the heatmap, overlaying the Image with the DataFrame
		@returns The Plot's axis, for downloading purposes.
		"""

		df = await DataCache.Load(input)
		img = await LoadImage()

		if df.empty: return None

		# Wrangle into an acceptable format.
		if {"x", "y", "value"}.issubset(df.columns):
			df = df.pivot(index="y", columns="x", values="value")

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
	@reactive.event(input.Update, input.Reset, input.Example, input.File, input.TextSize, input.Opacity, input.ColorMap, input.Algorithm, input.Style, input.Levels, input.Features, ignore_none=False, ignore_init=False)
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
	async def UpdateTableValue():
		"""
		@brief Updates the label for the Value input to display the current value.
		"""
		df = await DataCache.Load(input)

		rows, columns = df.shape
		row, column = int(input.TableRow()), int(input.TableCol())

		if 0 <= row <= rows and 0 <= column <= columns:
			ui.update_text(id="TableVal", label="Value (" + str(df.iloc[row, column]) + ")"),


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

			# https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.contourf.html#matplotlib.pyplot.contourf
			ui.input_select(id="Style", label=" Line Style", choices=["Solid", "Dashed", "Dashdot", "Dotted"], selected="Solid"),

			ui.input_slider(id="Levels", label="Number of Levels", value=20, min=1, max=100, step=1),

			# Customize what aspects of the heatmap are visible
			ui.input_checkbox_group(id="Features", label="Heatmap Features",
					choices={"x": "X Labels", "y": "Y Labels", "legend": "Legend"},
					selected=["legend"]),

			# Add the download buttons.
			ui.download_button("DownloadTable", "Download Table"),
		),

		# Add the main interface tabs.
		ui.navset_tab(
				ui.nav_panel("Interactive", ui.output_plot("Heatmap", height="90vh")),
				Table
		),
	)
)

app = App(app_ui, server)