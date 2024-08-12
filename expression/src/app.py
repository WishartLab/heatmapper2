#
# Heatmapper
# Expression
#
# This file contains the ShinyLive application for Expression Heatmapper.
# It can be run with the following command within this directory:
#		shinylive export . [site]
# Where [site] is the destination of the site folder.
#
# If you would rather deploy the application as a PyShiny application,
# run the following command within this directory:
#		shiny run
#
#


from shiny import App, reactive, render, ui, types
from matplotlib.pyplot import figure, style, subplots, close as fig_close
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from scipy.cluster import hierarchy
from scipy.stats import zscore
from scipy.interpolate import griddata
from tempfile import NamedTemporaryFile
from io import BytesIO
from numpy import arange, zeros_like, meshgrid, array, column_stack, linspace, min as n_min

from shared import Cache, NavBar, MainTab, FileSelection, Filter, ColumnType, TableOptions, Colors, InterpolationMethods, ClusteringMethods, DistanceMethods, InitializeConfig, Update, Msg, File

try:
	from user import config
except ImportError:
	from config import config


def server(input, output, session):
	# Information about the Examples
	Info = {
		"example1.txt": "This example dataset is sample input retrieved from the website for the Ashley Lab Heatmap Builder.",
		"example2.txt": "This example dataset is sample input retrieved from an online tutorial by Yan Cui (ycui2@uthsc.edu).",
		"example3.txt": "This example dataset is retrieved from the online supplement to Eisen et al. (1998), which is a very well known paper about cluster analysis and visualization. The details of how the data was collected are outlined in the paper."
	}

	DataCache = Cache("expression")
	Data = reactive.value(None)
	Valid = reactive.value(False)

	InitializeConfig(config, input)

	@reactive.effect
	@reactive.event(input.SourceFile, input.File, input.Example, input.Reset)
	async def UpdateData():
		Data.set((await DataCache.Load(input, p=ui.Progress())))
		Valid.set(False)
		Filter(Data().columns, ColumnType.Name, id="NameColumn")
		DataCache.Invalidate(File(input))


	def GetData(): return Table.data_view() if Valid() else Data()


	def HashString():
		"""
		@brief Returns the hash string used for the data cache.
		"""
		inputs = [
			File(input),
			config.NameColumn(),
			config.Features(),
			config.ScaleType(),
			input.CustomColors() if config.Custom() else config.ColorMap().split(),
			config.Interpolation(),
			config.Bins(),
			config.TextSize(),
			config.ClusterMethod(),
			config.DistanceMethod(),
			config.DPI(),
			config.Elevation(),
			input.mode(),
		]
		if config.Elevation() != 90: inputs.extend([config.Rotation(), config.Zoom(), config.InterpolationLevels(), config.MinScale(), config.Opacity()])
		return inputs


	def ProcessData(df):
		"""
		@brief Extracts the labels for each axis, and returns it alongside a DataFrame containing only the relevant data.
		@returns	A list containing the labels for the y axis, a list containing the labels for the x axis, and a
							DataFrame containing the loaded data without those two columns.
		"""

		name = config.NameColumn()
		if name not in df: return None, None , None

		# Drop the naming columns before linkage.
		data = df.drop(columns=Filter(df.columns, ColumnType.Name, all=True))
		x_labels = ["X" + name if list(data.columns).count(name) == 1 else "X" + name + f".{i+1}" for i, name in enumerate(data.columns)]
		return list(df[name]), x_labels, data


	def GenerateDendrogram(data, ax, orientation, progress, labels = [], invert=False):
		"""
		@brief General dendrogram generator.
		@param data: The DataFrame that contains the data to generate the dendrogram from.
		@param ax: The MatPlotLib Axis to assign tick marks to
		@param orientation: What orientation we should set the dendrogram to be. Can be "Left", "Right", "Top", or "Bottom"
		@param labels: An optional list of labels to add the dendrogram, labelling the X axis on Left/Right, and the Y on Top/Bottom
		@param invert: Whether to invert the DataFrame to generate Columnar dendrograms.
		@returns The dendrogram, mostly useful to aligning the Heatmap to the new ordering.
		"""

		if progress is not None: progress.inc(message="Creating linkage matrix...")
		method = config.ClusterMethod().lower()
		metric = config.DistanceMethod().lower()

		matrix = hierarchy.linkage(data.values.T if invert else data.values, method=method, metric=metric)

		if progress is not None: progress.inc(message="Creating dendrogram...")
		dendrogram = hierarchy.dendrogram(matrix, ax=ax, orientation=orientation.lower())

		# If there are labels, sort them according to the dendrogram.
		if labels: labels = [labels[i] for i in dendrogram["leaves"]]

		text_size = config.TextSize()

		# Add ticks depending on the orientation.
		if orientation == "Left" or orientation == "Right":
			ax.set_xticks([])
			ax.set_yticklabels(labels, fontsize=text_size)
		else:
			ax.set_yticks([])
			ax.set_xticklabels(labels, fontsize=text_size)

		return dendrogram


	def RenderDendrogram(data, labels, invert, progress):
		"""
		@brief Renders a Dendrogram
		@param data: The DataFrame
		@param labels: The labels for the Dendrogram
		@param invert: Whether to invert (Use for Column Dendrograms)
		@pararm progress: The progress bar to update when generating the Dendrogram.
		"""
		if data is None: return

		fig = figure(figsize=(12, 10))
		ax = fig.add_subplot(111)

		ax.spines["top"].set_visible(False)
		ax.spines["right"].set_visible(False)
		ax.spines["bottom"].set_visible(False)
		ax.spines["left"].set_visible(False)

		GenerateDendrogram(data, ax, config.Orientation(), progress, labels, invert=invert)
		return fig


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


	def Heatmap2D(df, ax_heatmap, p):
		colors = input.CustomColors() if config.Custom() else config.ColorMap().split()
		interpolation = config.Interpolation().lower()
		bins = config.Bins()

		return ax_heatmap.imshow(
			df,
			cmap=LinearSegmentedColormap.from_list("ColorMap", colors, N=bins),
			interpolation=interpolation,
			aspect="auto",
		)


	def Heatmap3D(df, ax_heatmap, p):

		cached = [
			File(input),
			input.CustomColors() if config.Custom() else config.ColorMap().split(),
			config.Bins(),
			config.InterpolationLevels(),
			config.Features(),
			config.MinScale(),
		]

		if not DataCache.In(cached):
			p.inc(message="Generating...")

			z = df.T.values.flatten()
			if config.MinScale():
				z += abs(n_min(z))

			color = input.mode()
			colors = input.CustomColors() if config.Custom() else config.ColorMap().split()
			cmap = LinearSegmentedColormap.from_list("ColorMap", colors, N=config.Bins())
			norm = Normalize(vmin=z.min(), vmax=z.max())
			c = norm(z)

			x_length, y_length = df.shape[0], df.shape[1]
			x_range, y_range = arange(x_length), arange(y_length)
			x, y = meshgrid(x_range, y_range)
			x, y = x.ravel(), y.ravel()
			width = depth = 1

			if config.InterpolationLevels() != 1:
				p.inc(message="Interpolating...")
				level = config.InterpolationLevels()
				x_grid, y_grid = meshgrid(linspace(0, x_length-1, x_length*2), linspace(0, y_length-1, y_length*2))
				x_new = x_grid.ravel()
				y_new = y_grid.ravel()

				points = column_stack((x, y))
				points_new = column_stack((x_new, y_new))

				z = griddata(points, z, points_new, method='cubic')
				c = griddata(points, c, points_new, method='cubic')

				x, y = x_new, y_new
				width /= level
				depth /= level

			if "legend" in config.Features():
				mappable = ScalarMappable(cmap=cmap, norm=norm)
				mappable.set_array(z)
			else: mappable = None
			DataCache.Store((x, y, width, depth, z, c, cmap, mappable), cached)

		x, y, width, depth, z, c, cmap, mappable = DataCache.Get(cached)

		ax_heatmap.view_init(elev=config.Elevation(), azim=config.Rotation())
		ax_heatmap.set_box_aspect(None, zoom=config.Zoom())
		return ax_heatmap.bar3d(x, y, zeros_like(x), width, depth, z, color=cmap(c), alpha=config.Opacity()), mappable


	def GenerateHeatmap():
		"""
		@brief Generates the Heatmap
		@returns The heatmap
		"""

		# A list of all the inputs for caching.
		inputs = HashString()

		# If we're rendering as images, fetch from the cache if we can
		if not DataCache.In(inputs):
			print(File(input), "NO")
			with ui.Progress() as p:
				p.inc(message="Reading input...")
				index_labels, x_labels, data = ProcessData(GetData())
				if data is None: return

				# Create a figure with a heatmap and associated dendrograms
				p.inc(message="Plotting...")
				color = input.mode()
				with  style.context('dark_background' if color == "dark" else "default"):
					fig = figure(figsize=(12, 10))
					gs = fig.add_gridspec(4, 2, height_ratios=[2, 8, 1, 1], width_ratios=[2, 8], hspace=0, wspace=0)

					# If we render the row dendrogram, we change the order of the index labels to match the dendrogram.
					# However, if we aren't rendering it, and thus row_dendrogram isn't defined, we simply assign df
					# To data, so the order changes when turning the toggle.
					if "row" in config.Features():
						ax_row = fig.add_subplot(gs[1, 0])
						row_dendrogram = GenerateDendrogram(data, ax_row, "Left", progress=p)
						ax_row.axis("off")
						leaves = row_dendrogram["leaves"]
						leaves.reverse()

						index_labels = [index_labels[i] for i in leaves]
						df = data.iloc[leaves]
					else:
						df = data

					# If we render the column dendrogram.
					if "col" in config.Features():
						ax_col = fig.add_subplot(gs[0, 1])
						col_dendrogram = GenerateDendrogram(data, ax_col, "Top", invert=True, progress=p)
						ax_col.axis("off")

					# Handle scaling
					if config.ScaleType() != "None": df = zscore(df, axis=1 if config.ScaleType() == "Row" else 0)

					if config.Elevation() == 90:
						ax_heatmap = fig.add_subplot(gs[1, 1])
						heatmap = Heatmap2D(df, ax_heatmap, p)
					else:
						ax_heatmap = fig.add_subplot(gs[1, 1], projection="3d")
						heatmap, mappable = Heatmap3D(df, ax_heatmap, p)
					text_size = config.TextSize()

					# If we render the Y axis.
					if "y" in config.Features():
						if config.Elevation() == 90: ax_heatmap.set_yticks(range(len(index_labels)))
						ax_heatmap.set_yticklabels(index_labels, fontsize=text_size)
						if config.Elevation() == 90: ax_heatmap.yaxis.tick_right()
					else:
						ax_heatmap.set_yticklabels([])

					# If we render the X axis.
					if "x" in config.Features():
						if config.Elevation() == 90: ax_heatmap.set_xticks(range(len(x_labels)))
						ax_heatmap.set_xticklabels(x_labels, rotation=90, fontsize=text_size)
					else:
						ax_heatmap.set_xticklabels([])

					if config.Elevation() != 90:
						if "z" in config.Features():
							ax_heatmap.tick_params(axis="z", labelsize=text_size)
						else:
							ax_heatmap.set_zticklabels([])

					# If we render the legend.
					if "legend" in config.Features():
						ax_cbar = fig.add_subplot(gs[3, 1])
						cbar = fig.colorbar(heatmap if input.Elevation() == 90 else mappable, cax=ax_cbar, orientation="horizontal")
						cbar.ax.tick_params(labelsize=text_size)


					b = BytesIO()
					fig.savefig(b, format="png", dpi=config.DPI())
					b.seek(0)
					DataCache.Store(b.read(), inputs)
					fig_close(fig)

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


	@output
	@render.plot
	def RowDendrogram():
		index_labels, _, data = ProcessData(GetData());
		with ui.Progress() as p:
			return RenderDendrogram(data=data, labels=index_labels, invert=False, progress=p)


	@output
	@render.plot
	def ColumnDendrogram():
		_, x_labels, data = ProcessData(GetData());
		with ui.Progress() as p:
			return RenderDendrogram(data=data, labels=x_labels, invert=True, progress=p)

	@reactive.effect
	@reactive.event(input.ExampleInfo)
	def ExampleInfo():
		Msg(ui.HTML(Info[input.Example()]))


	@render.download(filename="table.csv")
	def DownloadTable(): yield GetData().to_string()


	@render.download(filename="heatmap.png")
	def DownloadHeatmap(): yield DataCache.Get(HashString())


	@render.ui
	def Color():
		if config.Custom():
			return ui.input_select(id="CustomColors", label=None, choices=Colors, multiple=True, selectize=True, selected=["Blue", "White", "Yellow"])
		else:
			return config.ColorMap.UI(ui.input_select,
				make_inline=False, id="ColorMap", label=None,
				choices={
					"Blue White Yellow": "Blue/Yellow",
					"Red Black Green": "Red/Green",
					"Pink White Green": "Pink/Green",
					"Blue Green Yellow": "Blue/Green/Yellow",
					"Black Gray White": "Grayscale",
					"Red Orange Yellow Green Blue Indigo Violet": "Rainbow",
				}
			)


app_ui = ui.page_fluid(

	NavBar(),

	ui.layout_sidebar(
		ui.sidebar(

			FileSelection(
				examples={"example1.txt": "Example 1", "example2.txt": "Example 2", "example3.txt": "Example 3"},
				types=[".csv", ".txt", ".dat", ".tsv", ".tab", ".xlsx", ".xls", ".odf"],
				project="Expression"
			),

			TableOptions(config),

			# Settings pertaining to the Heatmap view.
			ui.panel_conditional(
				"input.MainTab === 'HeatmapTab'",
				Update(),

				ui.HTML("<b>Heatmap</b>"),

				# The column that holds names for the data.
				config.NameColumn.UI(ui.input_select, id="NameColumn", label="Names", choices=[], multiple=False),

				# https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
				config.ClusterMethod.UI(ui.input_select, id="ClusterMethod", label="Clustering", choices=ClusteringMethods),

				# https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html#scipy.spatial.distance.pdist
				config.DistanceMethod.UI(ui.input_select, id="DistanceMethod", label="Distance", choices=DistanceMethods, selected="Euclidean"),

				# Customize the text size of the axes.
				config.TextSize.UI(ui.input_numeric,id="TextSize", label="Text", min=1, max=50, step=1),

				# Define how the colors are scaled.
				config.ScaleType.UI(ui.input_select, id="ScaleType", label="Scale", choices=["Row", "Column", "None"], selected="Row"),

				# https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.imshow.html
				config.Interpolation.UI(ui.input_select, id="Interpolation", label="Inter", choices=InterpolationMethods, conditional="input.Elevation != 90"),

					ui.HTML("<b>3D</b>"),
					config.Elevation.UI(ui.input_numeric, id="Elevation",	label="Elevation"),
					config.Rotation.UI(ui.input_numeric, id="Rotation",	label="Rotation", conditional="input.Elevation != 90"),
					config.Zoom.UI(ui.input_numeric, id="Zoom",	label="Zoom", conditional="input.Elevation != 90", step=0.1),
					config.InterpolationLevels.UI(ui.input_numeric, id="InterpolationLevels",	label="Inter", conditional="input.Elevation != 90", step=1, min=1),
					config.MinScale.UI(ui.input_switch, id="MinScale",	label="Scaling", conditional="input.Elevation != 90"),
					config.Opacity.UI(ui.input_numeric, id="Opacity",	label="Opacity", conditional="input.Elevation != 90", min=0.0, max=1.0, step=0.1),

				ui.layout_columns(
					ui.HTML("<b>Colors</b>"),
					config.Custom.UI(ui.input_switch, make_inline=False, id="Custom", label="Custom"),
					col_widths=[4,8]
				),
				ui.output_ui("Color"),
				config.Bins.UI(ui.input_numeric, id="Bins", label="Number", min=3, step=1),

				ui.HTML("<b>Image Settings</b>"),
				config.Size.UI(ui.input_numeric, id="Size", label="Size", min=1),
				config.DPI.UI(ui.input_numeric, id="DPI", label="DPI", min=1),

				ui.HTML("<b>Features</b>"),
				config.Features.UI(ui.input_checkbox_group, make_inline=False, id="Features", label=None,
					choices={"row": "Row Dendrogram", "col": "Column Dendrogram", "x": "X Labels", "y": "Y Labels", "z": "Z Labels", "legend": "Legend"}
				),

				ui.download_button(id="DownloadHeatmap", label="Download"),
			),

			# Settings pertaining to the dendrogram view.
			ui.panel_conditional(
				"input.MainTab === 'RowTab' || input.MainTab === 'ColumnTab'",
				# Define the Orientation of the dendrogram in the Tab
				config.Orientation.UI(ui.input_select,id="Orientation", label="Orientation", choices=["Top", "Bottom", "Left", "Right"]),
			),
			padding="10px",
			gap="20px",
			width="250px",
		),

		# Add the main interface tabs.
		MainTab(
			ui.nav_panel("Row Dendrogram", ui.output_plot("RowDendrogram", height="90vh"), value="RowTab"),
			ui.nav_panel("Column Dendrogram", ui.output_plot("ColumnDendrogram", height="90vh"), value="ColumnTab"),
			m_type=ui.output_image
		),
		height="90vh",
	)
)

app = App(app_ui, server)
