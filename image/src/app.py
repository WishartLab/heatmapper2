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
from matplotlib.pyplot import subplots, colorbar, style
from matplotlib.tri import Triangulation
from PIL import Image
from tempfile import NamedTemporaryFile
from io import BytesIO
from numpy import meshgrid, arange, zeros_like, array, zeros, linspace, column_stack
from scipy.interpolate import griddata

from shared import Cache, MainTab, NavBar, FileSelection, Filter, ColumnType, TableOptions, InitializeConfig, ColorMaps, Update, Msg, File

try:
	from user import config
except ImportError:
	from config import config


def server(input, output, session):

	# Information regarding example files.
	Info = {
		"example1.txt": {
			"Image": "example1.jpg",
			"Description": "Input type: txt, jpg\nContents: Hypothetical example illustrating data overlaid on a satellite image. Input data are count or magnitude values within the overlaid grid sections."
		}
	}

	def HandleData(path, p=None):
		"""
		@brief A custom Data Handler for the Cache.
		@param path: The path to the file
		@returns A data object from the cache.
		@info This Data Handler supports png and jpg images as PIL.Image objects
		"""
		if path.suffix in [".bmp", ".gif", ".h5", ".hdf", ".ico", ".jpeg", ".jpg", ".tif", ".tiff", ".webp", ".png"]:
			return Image.open(path.resolve())
		else: return DataCache.DefaultHandler(path)
	DataCache = Cache("image", DataHandler=HandleData)
	Data = reactive.value(None)
	Valid = reactive.value(False)
	IMG = reactive.value(None)

	InitializeConfig(config, input)


	@reactive.effect
	@reactive.event(input.SourceFile, input.File, input.Example, input.Reset)
	async def UpdateData():
		Data.set((await DataCache.Load(input, p=ui.Progress())))
		Valid.set(False)
		DataCache.Invalidate(File(input))


	@reactive.effect
	@reactive.event(input.SourceFile, input.Example, input.Image)
	async def UpdateIMG(): IMG.set(await DataCache.Load(input,
			source_file=input.Image(),
			example_file=Info[input.Example()]["Image"]
		))


	def GetData(): return Table.data_view() if Valid() else Data()


	def HashString():
		inputs = [
			File(input),
			input.Image(),
			config.ColorMap(),
			config.Opacity(),
			config.Algorithm(),
			config.Levels(),
			config.Features(),
			config.TextSize(),
			config.DPI(),
			config.Quality(),
			input.mode()
		]
		if config.Elevation() != 90: inputs += [config.Elevation(), config.Rotation(), config.Zoom(), config.Slices()]
		return inputs


	@output
	@render.data_frame
	def Table(): Valid.set(True); return render.DataGrid(Data(), editable=True)


	@Table.set_patch_fn
	def UpdateTable(*, patch: render.CellPatch) -> render.CellValue:
		if config.Type() == "Integer": value = int(patch["value"])
		elif config.Type() == "Float": value = float(patch["value"])
		else: value = patch["value"]
		DataCache.Invalidate(File(input))
		return value


	def GenerateHeatmap():
		inputs = HashString()

		if not DataCache.In(inputs):
			with ui.Progress() as p:
				p.inc(message="Loading input...")
				df = GetData()

				p.inc(message="Loading image...")
				img = IMG()
				if img is None or df.empty: return None

				if img is not None:
					try:
						w, h = img.size
						img = img.resize((round(w * config.Quality()), round(h * config.Quality())))
					except TypeError: img = None

				# Wrangle into an acceptable format.
				p.inc(message="Formatting...")
				v_col = Filter(df.columns, ColumnType.Value)
				x_col = Filter(df.columns, ColumnType.X)
				y_col = Filter(df.columns, ColumnType.Y)

				if v_col != x_col and x_col != y_col and {v_col, x_col, y_col}.issubset(df.columns):
					df = df.pivot(index=x_col, columns=y_col, values=v_col).reset_index(drop=True)

				p.inc(message="Plotting...")
				color = input.mode()
				with  style.context('dark_background' if color == "dark" else "default"):

					cmap = config.ColorMap().lower()
					alpha = config.Opacity()
					algorithm = config.Algorithm().lower()
					levels = config.Levels()

					if config.Elevation() == 90:
						fig, ax = subplots()
						# Add the image as an overlay, if we have one.
						if img is not None:
							img = img.transpose(method=Image.FLIP_TOP_BOTTOM)
							ax.imshow(img, extent=[0, 1, 0, 1], aspect="auto",zorder=0)
						im = ax.contourf(df, cmap=cmap, extent=[0, 1, 0, 1], zorder=1, alpha=alpha, algorithm=algorithm, levels=levels)
						ax.invert_yaxis()

					else:
						fig, ax = subplots(subplot_kw={"projection": "3d"})

						z = df.values
						ax.set_zlim([0, z.max()])

						x, y = arange(df.shape[0]), arange(df.shape[1])
						x, y = meshgrid(x, y)

						if img is not None:
							arr = array(img) / 255.0
							ix, iy, _ = arr.shape

							x_new, y_new = meshgrid(linspace(0, df.shape[0]-1, ix), linspace(0, df.shape[1]-1, iy))

							points = column_stack((x.ravel(), y.ravel()))
							points_new = column_stack((x_new.ravel(), y_new.ravel()))
							z = griddata(points, z.flatten(), points_new, method='cubic').reshape(ix, iy)
							x, y = meshgrid(arange(ix), arange(iy))

							ax.plot_surface(x, y, zeros_like(x), rstride=1, cstride=1, facecolors=arr)

						ax.view_init(elev=config.Elevation(), azim=config.Rotation())
						ax.set_box_aspect(None, zoom=config.Zoom())
						im = ax.plot_surface(y, x, z, alpha=alpha, cmap=cmap)
						if config.Slices():
							ax.contour(y, x, z, zdir='z', offset=-z.min(), cmap=cmap)
							ax.contour(y, x, z, zdir='x', offset=x.min(), cmap=cmap)
							ax.contour(y, x, z, zdir='y', offset=y.min(), cmap=cmap)


					# Visibility of features
					if "legend" in input.Features():
						cbar = colorbar(im, ax=ax, label="Value")
						cbar.ax.tick_params(labelsize=config.TextSize())

					if "y" in config.Features(): ax.tick_params(axis="y", labelsize=config.TextSize())
					else: ax.set_yticklabels([])

					if "x" in config.Features(): ax.tick_params(axis="x", labelsize=config.TextSize())
					else: ax.set_xticklabels([])

					if config.Elevation() != 90:
						if "z" in config.Features(): ax.tick_params(axis="z", labelsize=config.TextSize())
						else: ax.set_zticklabels([])

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
	def DownloadHeatmap(): yield DataCache.Get(HashString())


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

	ui.panel_title(title=None, window_title="Image"),
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
				config.TextSize.UI(ui.input_numeric, id="TextSize", label="Text", min=1, max=50, step=1, tooltip="Change the text size of all axis labels. Axis labels can be toggled on and off in the 'Features' section at the bottom of this sidebar."),
				config.ColorMap.UI(ui.input_select, id="ColorMap", label="Color Map", choices=ColorMaps + ["Spring", "Summer", "Autumn", "Winter"], tooltip="Select a color scheme to use for the heatmap."),
				config.Algorithm.UI(ui.input_select, id="Algorithm", label="Contour", choices=["MPL2005", "MPL2014", "Serial", "Threaded"], tooltip="Select a algorithm used to generate the contours of the heatmap (convert the 2D data grid into smooth shapes). Default is MPL2014, while Threaded is best for large datasets."),
				config.Levels.UI(ui.input_numeric, id="Levels", label="Levels", min=1, step=1, tooltip="Specify the number of contour levels. A higher number of levels results in smoother transitions between values, but is more computationally expensive. "),
				config.Opacity.UI(ui.input_numeric, id="Opacity", label="Opacity", min=0.0, max=1.0, step=0.1, tooltip="Specify the opacity of the heatmap. 1.0 indicates full opacity, while lower values make the background image more visible."),

				ui.HTML("<b>3D</b>"),
				config.Elevation.UI(ui.input_numeric, id="Elevation", label="Elevation", tooltip="Control whether the plot is 2D or 3D. Any value other than 90 will display the plot in 3D, with the value specifying the elevation angle of the viewer in respect to the model. Change the angle back to 90 to display the plot in 2D."),
				config.Rotation.UI(ui.input_numeric, id="Rotation",	label="Rotation", conditional="input.Elevation != 90", tooltip="Change the angle of rotation of the viewer in respect to the model. For 3D plots only."),
				config.Zoom.UI(ui.input_numeric, id="Zoom",	label="Zoom", conditional="input.Elevation != 90", step=0.1, tooltip="Crop the view. For 3D plots only."),
				config.Slices.UI(ui.input_switch, id="Slices",	label="Slices", conditional="input.Elevation != 90", tooltip="Toggle on to display 2D projections of the heatmap on the XY, XZ, and YZ planes"),


				ui.HTML("<b>Image Settings</b>"),
				config.Quality.UI(ui.input_numeric, id="Quality", label="Quality", min=0.1, max=1.0, step=0.1, tooltip="Specify a multiplier to downscale the background image. Lower values decrease image quality and improve rendering speed. Set the value to 1.0 to use the original image with no downscaling."),
				config.Size.UI(ui.input_numeric, id="Size", label="Size", min=1, tooltip="Change the width (in pixels) of the heatmap on your screen."),
				config.DPI.UI(ui.input_numeric, id="DPI", label="DPI", min=1, tooltip="Specify the resolution of the image in pixels per inch. Higher DPI values result in higher quality images, but larger file sizes. This setting affects the heatmap on screen as well as the downloaded plot."),

				# Customize what aspects of the heatmap are visible
				ui.HTML("<b>Features</b>"),
				config.Features.UI(
					ui.input_checkbox_group, 
					make_inline=False, 
					id="Features", 
					label=None,
					choices={"x": "X Labels", "y": "Y Labels", "z": "Z Labels", "legend": "Legend"},
					tooltip="X and Y labels toggle the data labels along their respective axes. Z labels toggles the data labels along the Z axis if rendering as a 3D plot. Legend displays a colorbar legend on the heatmap."
				),

				ui.download_button(id="DownloadHeatmap", label="Download PNG"),
			),
			padding="10px",
			gap="20px",
			width="300px",
		),

		# Add the main interface tabs.
		MainTab(),
		#height="86vh",  # removes scroll bar but compresses image
	)
)

app = App(app_ui, server)
