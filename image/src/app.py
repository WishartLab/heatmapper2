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
from matplotlib.pyplot import subplots, colorbar, style, close as fig_close
from matplotlib.tri import Triangulation
from PIL import Image
from tempfile import NamedTemporaryFile
from io import BytesIO
from pandas import DataFrame
from numpy import meshgrid, arange, zeros_like, array, zeros, linspace, column_stack
from scipy.interpolate import griddata

from shared import Cache, Error, MainTab, NavBar, FileSelection, Filter, ColumnType, TableOptions, InitializeConfig, ColorMaps, Update, Msg, File

try:
	from user import config
except ImportError:
	from config import config


def server(input, output, session):

	# Information regarding example files.
	Info = {
		"example1.txt": {
			"Image": "example1.jpg",
			"Description": "<u>Input type:</u> .txt Data, .jpg Image <br><u>Contents:</u> Hypothetical example illustrating data overlaid on a satellite image. Input data are count or magnitude values within the overlaid grid sections."
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
		# catch error here
		p = ui.Progress()
		try:
			Data.set((await DataCache.Load(input, p=p)))
			Valid.set(False)
			DataCache.Invalidate(File(input))
		except:
			p.close()
			Error("File could not be loaded!\nData can be uploaded as a .csv, .tsv, .txt, .xslx, .dat, .tab, or .odf file. \nImages can be uploaded as a .bmp, .gif, .ico, .jpg, .tif, .webp, or .png file.")


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
			config.Legend(),
			input.LegendOrientation(),
			input.LegendSize(),
			input.LegendPadding(),
			config.TextSize(),
			config.DPI(),
			config.Quality(),
			input.mode(),
		]
		if config.Elevation() != 90: inputs += [config.Elevation(), config.Rotation(), config.Zoom(), config.Slices()]
		return inputs


	def CreateErrorImg(text, color, inputs):
		"""
		@brief Generates an image of the provided text
		@param text: The text to display as an error
		@param color: Hex color code for the text
		@param inputs: A list of all the inputs for caching (from HashString())
		@returns 
		"""
		# create image with error text
		fig, ax = subplots()
		ax.text(0, 50, text, color=color, fontsize=32)
		ax.set_xlim(0, 200)
		ax.set_ylim(0, 100)
		# make axes transparent
		[ax.spines[side].set_alpha(0.0) for side in ["top", "bottom", "left", "right"]]
		ax.tick_params(axis='both', which='both', reset=False, color=[0,0,0,0], labelcolor=[0,0,0,0])
		# save image to cache
		b = BytesIO()
		fig.savefig(b, format="png", dpi=100, bbox_inches="tight")
		b.seek(0)
		DataCache.Store(b.read(), inputs)
		fig_close(fig)
		# get image for display
		b = DataCache.Get(inputs)
		with NamedTemporaryFile(delete=False, suffix=".png") as temp:
			temp.write(b)
			temp.close()
			img: types.ImgData = {"src": temp.name, "width": "400px"}
			return img
		

	@output
	@render.data_frame
	def Table():
		df = Data()
		print(f"df.col: {df.columns}")
		print(type(df))
		if len(df.columns) == 0 or df is None:
			df = DataFrame({"_": ["No data to display! Please upload your data or select an example data set in the sidebar."]})
			return df

		# render data as editable table
		try:
			grid = render.DataGrid(df, editable=True)
			Valid.set(True)
			return grid	
		except Exception:
			pass
			return DataFrame({"Error": ["The provided input format cannot be rendered."]})


	@Table.set_patch_fn
	def UpdateTable(*, patch: render.CellPatch) -> render.CellValue:
		if config.Type() == "Integer": value = int(patch["value"])
		elif config.Type() == "Float": value = float(patch["value"])
		else: value = patch["value"]
		DataCache.Invalidate(File(input))
		return value
	
	# Info text in welcome tab
	@render.ui
	def Welcome():
		return ui.HTML("""
			<h1>Image Heatmaps</h1>
			Image heatmaps visualize data from an input table over a user-supplied image. <br>
			Upload a data file and an image file in the sidebar to get started, or select 'Example' to check out a pre-loaded example. <br>
			Navigate to the 'Heatmap' tab to see the heatmap, or 'Table' to look at the input data.
			
			<br><br>
			<img src="https://github.com/WishartLab/heatmapper2/wiki/assets/Image.png" alt="Image"; style="max-width:500px;">
				 
			<br><br><h3>Format</h3>
			<i>Input data can be formatted as follows:</i>
				 <ul>
				 <li>A precomputed 2D matrix of values, which will be displayed as-is as a heatmap. (See Ex 1: Map)</li>
				 <li>X and Y columns indicating coordinates, with an associated value column.</li>
				 </ul>
			<i>Image heatmaps can be generated from the following file formats:</i>
			<table style="border-spacing: 100px";>
			<tr>
				<th>Image Files</th>
				<th>Table Files</th>
			</tr>
			<tr>
				<td style="padding-right:75px;">
					<li>.bmp</li>
					<li>.gif</li>
					<li>.h5</li>
					<li>.hdf</li>
					<li>.ico</li>
					<li>.jpeg</li>
					<li>.tif</li>
					<li>.tiff</li>
					<li>.webp</li>
					<li>.png</li>
				</td>
				<td style="vertical-align:top;">
					<li>.csv</li>
					<li>.dat</li>
					<li>.odf</li>
					<li>.tab</li>
					<li>.tsv</li>
					<li>.txt</li>
					<li>.xls</li>
					<li>.xlsx</li>
				</td>
			</tr>
			</table>
				
			<br><h3>Interface</h3>
			Click on the '?' icon beside sidebar options to read more about them.
		""")


	def GenerateHeatmap():
		inputs = HashString()

		if not DataCache.In(inputs):
			with ui.Progress() as p:
				p.inc(message="Loading input...")
				df = GetData()

				p.inc(message="Loading image...")
				img = IMG()
				if img is None or df.empty: 
					return CreateErrorImg("No data to display!\n\nPlease upload your data or select an example data set in the sidebar.", "#027bc2", inputs)

				if img is not None:
					try:
						w, h = img.size
						img = img.resize((round(w * config.Quality()), round(h * config.Quality())))
					except TypeError: img = None

				if img is None:
					return CreateErrorImg("Please upload a background image.", "#027bc2", inputs)

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
						print(f"img.size: {img.size}")
						# calculate subplot dimensions
						w, h = img.size
						m = float(max(w, h, 15))
						w_new = float(w) * (15.0 / m)
						h_new = float(h) * (15.0 / m)
						if input.LegendOrientation() in ["Right", "Left"]:
							pad = float(input.LegendPadding() / 100) * w_new
							figsize = (w_new + (3.0 * pad), h_new)
						else:
							pad = float(input.LegendPadding() / 100) * h_new
							figsize = (w_new, h_new + (3.0 * pad))
						print(figsize)
						
						fig, ax = subplots(figsize=figsize)
						# Add the image as an overlay, if we have one.
						if img is not None:
							img = img.transpose(method=Image.FLIP_TOP_BOTTOM)
							ax.imshow(img, extent=[0, 1, 0, 1], aspect="auto",zorder=0)
						try:
							im = ax.contourf(df, cmap=cmap, extent=[0, 1, 0, 1], zorder=1, alpha=alpha, algorithm=algorithm, levels=levels)
						except:
							return CreateErrorImg("Data could not be parsed, please check formatting.\nData can be uploaded as a .csv, .tsv, .txt, .xslx, .dat, .tab, or .odf file.", "#027bc2", inputs)
						ax.invert_yaxis()

					else:
						fig, ax = subplots(subplot_kw={"projection": "3d"}, constrained_layout=True)

						z = df.values
						ax.set_zlim([0, z.max()])

						x, y = arange(df.shape[0]), arange(df.shape[1])
						x, y = meshgrid(x, y)

						if img is not None:
							# normalize image
							arr = array(img) / 255.0
							ix, iy, _ = arr.shape
							print(f"arr:\t{arr}")
							print(f"arr.shape:\t{arr.shape}")

							x_new, y_new = meshgrid(linspace(0, df.shape[0]-1, ix), linspace(0, df.shape[1]-1, iy))

							points = column_stack((x.ravel(), y.ravel()))
							points_new = column_stack((x_new.ravel(), y_new.ravel()))
							z = griddata(points, z.flatten(), points_new, method='cubic').reshape(ix, iy)
							x, y = meshgrid(arange(iy), arange(ix))

							ax.plot_surface(y, x, zeros_like(x), rstride=1, cstride=1, facecolors=arr)

						ax.view_init(elev=config.Elevation(), azim=config.Rotation())
						ax.set_box_aspect(None, zoom=config.Zoom())
						im = ax.plot_surface(y, x, z, alpha=alpha, cmap=cmap)
						if config.Slices():
							ax.contour(y, x, z, zdir='z', offset=-z.min(), cmap=cmap)
							ax.contour(y, x, z, zdir='x', offset=x.min(), cmap=cmap)
							ax.contour(y, x, z, zdir='y', offset=y.min(), cmap=cmap)


					# Visibility of features
					if "legend" in input.Features():
						cbar = colorbar(
							im, 
							ax=ax, 
							label=config.Legend(), 
							location=input.LegendOrientation().lower(),
							shrink=input.LegendSize()/100,
							pad=(input.LegendPadding()/100), 
						)
						cbar.ax.tick_params(labelsize=config.TextSize())
						cbar.ax.set_ylabel(cbar.ax.get_ylabel(), fontsize=config.TextSize())
						cbar.ax.set_xlabel(cbar.ax.get_xlabel(), fontsize=config.TextSize())

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


	@render.download(filename=lambda: f"table{config.TableType()}")
	def DownloadTable(): 
		data = GetData()
		
		# return error if no data to download
		if data.empty:
			Error("The downloaded table is empty! Please upload your data or select an example data set in the sidebar.")
		
		file_contents = data.to_string()
		yield file_contents


	@render.download(filename=lambda: f"heatmap{config.HeatmapType()}")
	def DownloadHeatmap(): yield DataCache.Get(HashString())


	@render.download(filename=lambda: f"settings{config.SettingType()}")
	def DownloadSettings(): 
		'''
		Download a table file containing current config settings
		'''
		yield f"Data Filename:\t{File(input)}\nImage Filename:\t{input.Image()}\nText Size:\t{config.TextSize()}\nColor Map:\t{config.ColorMap()}\nContour Algorithm:\t{config.Algorithm()}\nContour Levels:\t{config.Levels()}\nHeatmap Opacity:\t{config.Opacity()}\nImage Quality:\t{config.Quality()}\nResolution(DPI):\t{config.DPI()}\nSelected Features:\t{config.Features()}\nLegend Title:\t{config.Legend()}"


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

	ui.panel_title(title=None, window_title="Image"),
	NavBar(),

	ui.layout_sidebar(
		ui.sidebar(

			FileSelection(
				examples={
					"example1.txt": "Ex 1: Map"
				}, 
				types=[".csv", ".txt", ".dat", ".tsv", ".tab", ".xlsx", ".xls", ".odf"], 
				project="Image"),

			ui.panel_conditional("input.SourceFile === 'Upload'", ui.input_file("Image", "Choose your Image File",
				accept=[".bmp", ".gif", ".h5", ".hdf", ".ico", ".jpeg", ".jpg", ".tif", ".tiff", ".webp", ".png"],
				multiple=False)),

			TableOptions(config),

			ui.panel_conditional(
				"input.MainTab === 'HeatmapTab'",

				Update(),

				ui.HTML("<b>Heatmap</b>"),
				config.TextSize.UI(ui.input_numeric, id="TextSize", label="Text Size", min=1, max=50, step=1, tooltip="Change the text size of all axis labels. Axis labels can be toggled on and off in the 'Features' section at the bottom of this sidebar."),
				config.ColorMap.UI(ui.input_select, id="ColorMap", label="Color Map", choices=ColorMaps + ["Spring", "Summer", "Autumn", "Winter"], tooltip="Select a color scheme to use for the heatmap."),
				config.Algorithm.UI(ui.input_select, id="Algorithm", label="Contour Algorithm", choices=["MPL2005", "MPL2014", "Serial", "Threaded"], tooltip="Select a algorithm used to generate the contours of the heatmap (convert the 2D data grid into smooth shapes). Default is MPL2014, while Threaded is best for large datasets."),
				config.Levels.UI(ui.input_numeric, id="Levels", label="Contour Levels", min=1, step=1, tooltip="Specify the number of contour levels. A higher number of levels results in smoother transitions between values, but is more computationally expensive. "),
				config.Opacity.UI(ui.input_slider, id="Opacity", label="Heatmap Opacity", min=0.0, max=1.0, step=0.1, tooltip="Specify the opacity of the heatmap. 1.0 indicates full opacity, while lower values make the background image more visible."),

				ui.HTML("<b>3D</b>"),
				config.Elevation.UI(ui.input_numeric, id="Elevation", label="View Elevation", tooltip="Control whether the plot is 2D or 3D. Any value other than 90 will display the plot in 3D, with the value specifying the elevation angle of the viewer in respect to the model. Change the angle back to 90 to display the plot in 2D."),
				config.Rotation.UI(ui.input_numeric, id="Rotation",	label="Rotation", conditional="input.Elevation != 90", tooltip="Change the angle of rotation of the viewer in respect to the model. For 3D plots only."),
				config.Zoom.UI(ui.input_numeric, id="Zoom",	label="Zoom", conditional="input.Elevation != 90", step=0.1, tooltip="Crop the view. For 3D plots only."),
				config.Slices.UI(ui.input_switch, id="Slices",	label="Slices", conditional="input.Elevation != 90", tooltip="Toggle on to display 2D projections of the heatmap on the XY, XZ, and YZ planes"),


				ui.HTML("<b>Image Settings</b>"),
				config.Quality.UI(ui.input_slider, id="Quality", label="Image Quality", min=0.1, max=1.0, step=0.1, tooltip="Specify a multiplier to downscale the background image. Lower values decrease image quality and improve rendering speed. Set the value to 1.0 to use the original image with no downscaling."),
				config.Size.UI(ui.input_numeric, id="Size", label="Heatmap Size", min=1, tooltip="Change the width (in pixels) of the heatmap on your screen."),
				config.DPI.UI(ui.input_numeric, id="DPI", label="Resolution (DPI)", min=1, tooltip="Specify the resolution of the image in pixels per inch. Higher DPI values result in higher quality images, but larger file sizes. This setting affects the heatmap on screen as well as the downloaded plot."),

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
				# legend options
				ui.panel_conditional("input.Features.includes('legend')", 
					config.Legend.UI(
						ui.input_text, 
						id="Legend", 
						label="Legend Title", 
						tooltip="Provide a title for the colorbar legend."
					),
					ui.input_slider(
						id="LegendSize",
						label="Legend Size",
						min=10,
						max=100,
						step=1,
						value=100,
					),
					ui.input_select(
						id="LegendOrientation",
						label="Legend Orientation",
						choices=["Left", "Right", "Top", "Bottom"],
						selected="Right",
					),
					ui.input_slider(
						id="LegendPadding",
						label="Legend Padding",
						min=0,
						max=99,
						step=1,
						value=5,
					),
				),
				
				ui.HTML("<b>Downloads</b>"),
				config.HeatmapType.UI(ui.input_radio_buttons, make_inline=False, id="HeatmapType", label="Heatmap File Type", choices=[".png", ".jpg"], inline=True),
				ui.download_button(id="DownloadHeatmap", label="Download Heatmap"),

				config.SettingType.UI(ui.input_radio_buttons, make_inline=False, 
				id="SettingType", label="Settings File Type", choices=[".txt", ".csv", ".tsv", ".xlsx"], inline=True),
				ui.download_button(id="DownloadSettings", label="Download Current Settings"),
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
