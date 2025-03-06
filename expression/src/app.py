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


from shiny import App, reactive, render, ui, types
from matplotlib.pyplot import figure, style, subplots, close as fig_close
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from pandas import DataFrame
from scipy.cluster import hierarchy
from scipy.stats import zscore
from scipy.interpolate import griddata
from tempfile import NamedTemporaryFile
from io import BytesIO
from numpy import arange, zeros_like, meshgrid, array, column_stack, linspace, min as n_min

from shared import Error, Cache, NavBar, MainTab, FileSelection, Filter, ColumnType, TableOptions, Colors, InterpolationMethods, ClusteringMethods, DistanceMethods, InitializeConfig, Update, Msg, File

try:
	from user import config
except ImportError:
	from config import config


# global variable :(
# saves the largest size of "Expand" heat map encountered so far
EXPANDED_SIZE = 0


def server(input, output, session):
	# Information about the Examples
	Info = {
		"mouse_leukemia.tsv": '<u>Input type:</u> .tsv Data <br><u>Contents:</u> A subset of gene expression data from a murine study on acute myeloid leukemia (AML). Columns represent samples, while rows are genes. The samples have various mutations related to altered DNA methylation. <br><u>Source:</u> Shih, A. H. et al. (2017). PMCID: <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC5413413/"; "target="_blank">PMC5413413</a>',
		"human_liver_subset.tsv": '<u>Input type:</u> .tsv Data <br><u>Contents:</u> Subset of human liver RNA-Seq data from ARCHS4. Columns represent samples, while rows are genes. The dataset is composed of samples from multiple different experiments. <br><u>Source:</u> <a href="https://www.kaggle.com/datasets/lachmann12/human-liver-rnaseq-gene-expression-903-samples?resource=download"; target="_blank">kaggle.com</a>.',
		"example3.txt": '<u>Input type:</u> .txt Data <br><u>Contents:</u> Large dataset of gene expression in <i>S. cerevisiae</i> under varying conditions. Data was collected during the cell division cycle, alpha factor arrest, centrifugal elutriation, sporulation, high temperature shock, low temperature shock, and diauxic shift. Data was collected using DNA microarrays.<br><u>Source:</u> Eisen et al. (1998). PMCID: <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC24541/#F2"; target="_blank">PMC24541</a> DOI: 10.1073/pnas.95.25.14863'
	}

	DataCache = Cache("expression")
	Data = reactive.value(None)
	Valid = reactive.value(False)

	InitializeConfig(config, input)

	@reactive.effect
	@reactive.event(input.SourceFile, input.File, input.Example, input.Reset)
	async def UpdateData():
		p = ui.Progress()
		try:
			Data.set((await DataCache.Load(input, p=p)))
			Valid.set(False)
			Filter(Data().columns, ColumnType.Name, id="NameColumn")
			DataCache.Invalidate(File(input))
		except:
			p.close()
			Error("File could not be loaded!\nData can be uploaded as a .csv, .tsv, .txt, .xslx, .dat, .tab, or .odf file.")


	def GetData(): return Table.data_view() if Valid() else Data()


	def HashString():
		"""
		@brief Returns the hash string used for the data cache.
		"""
		inputs = [
			File(input),
			config.NameColumn(),
			config.Features(),
			#config.N(),
			config.ScaleType(),
			input.CustomColors() if config.Custom() else config.ColorMap().split(),
			config.Interpolation(),
			config.Bins(),
			config.TextSize(),
			config.ClusterMethod(),
			config.DistanceMethod(),
			config.DPI(),
			config.AutoSize(),
			config.Elevation(),
			input.mode(),
		]
		if config.Elevation() != 90: inputs.extend([config.Rotation(), config.Zoom(), config.InterpolationLevels(), config.MinScale(), config.Opacity()])
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
		ax.text(0, 50, text, color=color, fontsize=12)
		ax.set_xlim(0, 200)
		ax.set_ylim(0, 100)
		# make axes transparent
		[ax.spines[side].set_alpha(0.0) for side in ["top", "bottom", "left", "right"]]
		ax.tick_params(axis='both', which='both', reset=False, color=[0,0,0,0], labelcolor=[0,0,0,0])
		if inputs == "fig":
			fig_close(fig)
			return fig
		else:
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


	def ProcessData(df):
		"""
		@brief Extracts the labels for each axis, and returns it alongside a DataFrame containing only the relevant data.
		@returns	A list containing the labels for the y axis, a list containing the labels for the x axis, and a
							DataFrame containing the loaded data without those two columns.
		"""

		name = config.NameColumn()
		if name not in df: return None, None , None

		# check if 

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

		try:
			matrix = hierarchy.linkage(data.values.T if invert else data.values, method=method, metric=metric)
		except:
			Error("Could not generate heatmap from the uploaded data. Please check your data format.")
			return None

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
	def Table(): 
		df = Data()

		# instruct users to load data if table is empty
		if len(df.columns) == 0 or df is None:
			df = DataFrame({"_": ["No data to display! Please upload your data or select an example data set in the sidebar."]})
			return df

		# render data as editable table
		try:
			grid = render.DataGrid(df, editable=True)
			Valid.set(True)
			return grid
		except Exception:
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
			<h1>Expression Heatmaps</h1>
			Expression heatmaps display data from data from transcriptomic (microarray or RNAseq), proteomic or metabolomic experiments. Data can be unclustered, or pre-clustered. Clustering can be performed by toggling on column or row dendrograms in the sidebar. <br>
			Upload a data file in the sidebar to get started, or select 'Example' to check out a pre-loaded example. <br>
			Navigate to the 'Heatmap' tab to see the heatmap, or 'Table' to look at the input data.
			
			<br><br>
			<img src="https://github.com/WishartLab/heatmapper2/wiki/assets/Expression.png" alt="Image"; style="max-width:500px;">
				 
			<br><br><h3>Format</h3>
			A Name column is required in the input table, which is used to generate axis labels. All other columns are plotted as data. If columns do not have names - or have repeating names - a unique namming scheme will be applied.
				 
			<i>Expression heatmaps can be generated from the following file formats:</i>
				<li>.csv</li>
				<li>.dat</li>
				<li>.odf</li>
				<li>.tab</li>
				<li>.tsv</li>
				<li>.txt</li>
				<li>.xls</li>
				<li>.xlsx</li>
			
			<br><h3>Interface</h3>
			Click on the '?' icon beside sidebar options to read more about them.
		""")


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

			try:
				z = df.T.values.flatten()
			except:
				z = df.T.flatten()
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
		global EXPANDED_SIZE
		size = 0

		# A list of all the inputs for caching.
		inputs = HashString()

		# If we're rendering as images, fetch from the cache if we can
		if not DataCache.In(inputs):
			with ui.Progress() as p:
				p.inc(message="Reading input...")
				index_labels, x_labels, data = ProcessData(GetData())
				if data is None or len(data.index) == 0: 
					return CreateErrorImg("No data to display!\n\nPlease upload your data or select an example data set in the sidebar.", "#027bc2", inputs)

				# Create a figure with a heatmap and associated dendrograms
				p.inc(message="Plotting...")
				color = input.mode()
				with  style.context('dark_background' if color == "dark" else "default"):
					fig = figure(figsize=(12, 10))
					gs = fig.add_gridspec(4, 2, height_ratios=[2, 8, 1, 1], width_ratios=[2, 8], hspace=0, wspace=0)

					# If we render the row dendrogram, we change the order of the index labels to match the dendrogram.
					# However, if we aren't rendering it, and thus row_dendrogram isn't defined, we simply assign df
					# To data, so the order changes when turning the toggle.
					if "row" in config.Features() and config.Elevation() == 90:
						ax_row = fig.add_subplot(gs[1, 0])
						row_dendrogram = GenerateDendrogram(data, ax_row, "Left", progress=p)
						if row_dendrogram is None:
							return CreateErrorImg("Error generating heat map. Please check your data format.", "#027bc2", inputs)
						ax_row.axis("off")
						leaves = row_dendrogram["leaves"]
						leaves.reverse()

						index_labels = [index_labels[i] for i in leaves]
						df = data.iloc[leaves]
					else:
						df = data
					# get # of columns to calculate expanded view size
					num_col = max(len(df.columns), df.shape[0])
					
					# If we render the column dendrogram.
					if "col" in config.Features() and config.Elevation() == 90:
						ax_col = fig.add_subplot(gs[0, 1])
						col_dendrogram = GenerateDendrogram(data, ax_col, "Top", invert=True, progress=p)
						if col_dendrogram is None:
							return CreateErrorImg("Error generating heat map. Please check your data format.", "#027bc2", inputs)
						ax_col.axis("off")

					# Handle scaling
					if config.ScaleType() != "None": df = zscore(df, axis=1 if config.ScaleType() == "Row" else 0)

					if config.Elevation() == 90:
						ax_heatmap = fig.add_subplot(gs[1, 1])
						heatmap = Heatmap2D(df, ax_heatmap, p)
					else:
						ax_heatmap = fig.add_subplot(gs[1, 1], projection="3d")
						heatmap, mappable = Heatmap3D(df, ax_heatmap, p)
					
					# if "expand" is selected, set text to 8
					# if config.AutoSize() == "expand":
					# 	text_size = 8
					# else:
					# 	text_size = config.TextSize()
					text_size = config.TextSize()

					# If we render the Y axis.
					# TODO: show every N-th (n = config.N()...)
					if "y" in config.Features():
						if config.Elevation() == 90: ax_heatmap.set_yticks(range(len(index_labels)))
						ax_heatmap.set_yticklabels(index_labels, fontsize=text_size)
						if config.Elevation() == 90: ax_heatmap.yaxis.tick_right()
					else:
						ax_heatmap.set_yticklabels([])

					# If we render the X axis.
					# TODO: show every N-th (n = config.N()...)
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


					# set image size based on config
					# if config.AutoSize() == "expand":
					# 	print(f"NUM COL: {num_col}")
					# 	size = num_col * (1/3) * num_col
					# 	if size < 1000:
					# 		size = 1000	
					# 	print(f"size: {size}")					
					# 	# save size to global variable to be used when loading from cache
					# 	if size > EXPANDED_SIZE:
					# 		EXPANDED_SIZE = size
					
					# calculate dpi for auto expand
					# if config.AutoSize() == "expand":
					# 	dpi = size * 0.15
					# 	if config.DPI() > dpi:
					# 		dpi = config.DPI()
					# else:
					# 	dpi = config.DPI()
					dpi = config.DPI()
					
					# catch invalid dpi values
					if dpi > 1000:
						dpi = 1000
					elif dpi < 5:
						dpi = 5

					b = BytesIO()
					fig.savefig(b, format="png", dpi=config.DPI())
					b.seek(0)
					DataCache.Store(b.read(), inputs)
					fig_close(fig)

		# get image size		
		if size == 0:  # loading from cache
			if config.AutoSize() == "fit":
				size = 500
		# 	elif config.AutoSize() == "expand":
		# 		size = EXPANDED_SIZE
		# 		print(f"expand size: {size}")					
			else:
				size = config.Size()
		
		b = DataCache.Get(inputs)
		#with NamedTemporaryFile(delete=False, suffix=".png") as temp:
		with NamedTemporaryFile(delete=False, suffix=config.HeatmapType()) as temp:
			temp.write(b)
			temp.close()
			img: types.ImgData = {"src": temp.name, "height": f"{size}px"}
			return img


	@output
	@render.image(delete_file=True)
	def Heatmap(): 
		return GenerateHeatmap()


	@output
	@render.image(delete_file=True)
	@reactive.event(input.Update)
	def HeatmapReactive(): return GenerateHeatmap()


	@output
	@render.plot
	def RowDendrogram():
		index_labels, _, data = ProcessData(GetData())
		# instruct user to upload files if empty data
		if data is None:
			return CreateErrorImg("No data to display!\n\nPlease upload your data or select an example data set in the sidebar.", "#027bc2", inputs="fig")
		with ui.Progress() as p:
			return RenderDendrogram(data=data, labels=index_labels, invert=False, progress=p)


	@output
	@render.plot
	def ColumnDendrogram():
		_, x_labels, data = ProcessData(GetData())
		# instruct user to upload files if empty data
		if data is None:
			return CreateErrorImg("No data to display!\n\nPlease upload your data or select an example data set in the sidebar.", "#027bc2", inputs="fig")
		with ui.Progress() as p:
			return RenderDendrogram(data=data, labels=x_labels, invert=True, progress=p)

	@reactive.effect
	@reactive.event(input.ExampleInfo)
	def ExampleInfo():
		Msg(ui.HTML(Info[input.Example()]))


	@render.download(filename=lambda: f"table{config.TableType()}")
	def DownloadTable(): 
		data = GetData()
		
		# return error if no data to download
		if data.empty:
			Error("The downloaded table is empty! Please upload your data or select an example data set in the sidebar.")
		
		file_contents = data.to_string()
		yield file_contents


	@render.download(filename=lambda: f"heatmap{config.HeatmapType()}")
	def DownloadHeatmap(): 
		yield DataCache.Get(HashString())


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

	ui.panel_title(title=None, window_title="Expression"),
	NavBar(),

	ui.layout_sidebar(
		ui.sidebar(

			FileSelection(
				examples={
					"mouse_leukemia.tsv": "Ex1: Mouse AML", 
			  		"human_liver_subset.tsv": "Ex2: Human Liver", 
					"example3.txt": "Ex3: S. cerevisiae"},
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
				config.NameColumn.UI(ui.input_select, id="NameColumn", label="Y-axis Label", choices=[], multiple=False, tooltip="Select a column in your data to use for Y-axis labels. Common columns are 'NAME' or 'UNIQID'. Labels can be toggled on/off in the 'Features' section at the bottom of this sidebar."),

				# https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
				config.ClusterMethod.UI(ui.input_select, id="ClusterMethod", label="Clustering Method", choices=ClusteringMethods, tooltip=ui.HTML('''Specify the clustering method used to group input data. To enable clustering by row, toggle on 'Row Dendrogram' in the 'Features' section at the bottom of this sidebar. To enable clustering by column, toggle on 'Column Dendrogram'. Disabling the dendrograms displays input data unclustered. <br><br>Read more <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage" target="_blank">here</a>''')),

				# https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html#scipy.spatial.distance.pdist
				config.DistanceMethod.UI(ui.input_select, id="DistanceMethod", label="Distance Method", choices=DistanceMethods, selected="Euclidean",
				tooltip=ui.HTML('Specify a method to calculate the distance between data points. <br>Read more <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage" target="_blank">here</a>',)),

				# Customize the text size of the axes.
				config.TextSize.UI(ui.input_numeric,id="TextSize", label="Text Size", min=1, max=50, step=1, tooltip="Change the text size of all axis labels. Axis labels can be toggled on and off in the 'Features' section at the bottom of this sidebar."),

				# Define how the colors are scaled.
				config.ScaleType.UI(ui.input_select, id="ScaleType", label="Scale by:", choices=["Row", "Column", "None"], selected="Row", tooltip="Normalize cell values to either row or column using z-scores. For each row or column, the elements are transformed to have a mean of 0 and a standard deviation of 1. This makes data comparable across rows or columns, which is useful for clustering. Select 'None' to display data without normalization."),

				# https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.imshow.html
				config.Interpolation.UI(ui.input_select, id="Interpolation", label="Intrpl Method", choices=InterpolationMethods, conditional="input.Elevation == 90", tooltip=ui.HTML('Specify an interpolation algorithm to apply to the figure. This can cause values to bleed together and appear smoother. <br>Read more <a href="https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.imshow.html" target="_blank">here</a>.')),

				ui.HTML("<b>3D</b>"),
				config.Elevation.UI(ui.input_numeric, id="Elevation", label="View Elevation", tooltip="Control whether the plot is 2D or 3D. Any value other than 90 will display the plot in 3D, with the value specifying the elevation angle of the viewer in respect to the model. Change the angle back to 90 to display the plot in 2D."),
				config.Rotation.UI(ui.input_numeric, id="Rotation",	label="View Rotation", conditional="input.Elevation != 90", tooltip="Change the angle of rotation of the viewer in respect to the model. For 3D plots only."),
				config.Zoom.UI(ui.input_numeric, id="Zoom",	label="Zoom", conditional="input.Elevation != 90", step=0.1,tooltip="Crop the view. For 3D plots only."),
				config.InterpolationLevels.UI(ui.input_numeric, id="InterpolationLevels",	label="Intrpl Level", conditional="input.Elevation != 90", step=1, min=1, tooltip="Specify a multiplier for the resolution of the 3D plot. For example, a value of 2 will interpolate the data from an NxM to a 2Nx2M, effectively quadrupling the resolution of each data point by interpolating it into 4. This results in a smoother looking plot, but can be computationally expensive for large datasets."),
				config.MinScale.UI(ui.input_switch, id="MinScale", label="Scaling", conditional="input.Elevation != 90", tooltip="Scale the height of all points by the minimum value. This removes negative values and prevents data from extending below the XY plane. For 3D plots only."),
				config.Opacity.UI(ui.input_slider, id="Opacity", label="Opacity", conditional="input.Elevation != 90", min=0.0, max=1.0, step=0.1, tooltip="Change the opacity of the height bars in 3D plots."),

				ui.layout_columns(
					ui.HTML("<b>Colors</b>"),
					config.Custom.UI(ui.input_switch, make_inline=False, id="Custom", label="Custom", tooltip="Select colors to use in the heatmap, with the first color representing low values, and the last color representing high values. If the 'Custom' checkbox is enabled you may specify up to 12 colors. A minimum of two colors are needed."),
					col_widths=[4,8]
				),
				ui.output_ui("Color"),
				config.Bins.UI(ui.input_numeric, id="Bins", label="# of Color Bins", min=3, step=1, tooltip="Specify the number of color bins to use. A higher number of color bins results in a smoother gradient between neighbouring values. Fewer bins results in more distinct colors."),

				ui.HTML("<b>Features</b>"),
				config.Features.UI(ui.input_checkbox_group, make_inline=False, id="Features", label=None, choices={"row": "Row Dendrogram", "col": "Column Dendrogram", "x": "X Labels", "y": "Y Labels", "z": "Z Labels", "legend": "Legend"}, tooltip="Row Dendrogram enables clustering of rows. Column Dendrogram enables clustering of columns. X and Y labels toggle the data labels along their respective axes. Z labels toggles the data labels along the Z axis if rendering as a 3D plot. Legend displays a colorbar legend on the heatmap.",
				),
				#config.N.UI(ui.input_slider, id="N", label="Show N-th Label", min=1, max=25, step=1, tooltip=ui.HTML("Display every N-th label. <br>For example, a value of 2 will display only every second label on visualized axes. <br>Set to 1 to display every label.")),

				ui.HTML("<b>Image Settings</b>"),
				ui.div(
					config.DPI.UI(ui.input_numeric, id="DPI", label="Resolution (DPI)", min=5, tooltip="Specify the resolution of the image in pixels per inch. Higher DPI values result in higher quality images, but larger file sizes. This setting affects the heatmap on screen as well as the downloaded plot."),
					ui.HTML("<u>Image Size</u><br><br>"),
					ui.div(
						ui.div(
							config.AutoSize.UI(ui.input_radio_buttons,
						  		make_inline=False, id="AutoSize", label=None, choices={"custom": "Custom Width", "fit": "Fit to Screen"}, 
							),
							style="flex: 1; padding-top: 15px;",
						),
						ui.div(
							config.Size.UI(ui.input_numeric, gap="0px", id="Size", label=None, min=1,
					  		tooltip=ui.HTML("Select <b>Custom Width</b> to specify a custom width (in pixels) for the heat map on your screen. <br><br>Select <b>Fit to Screen</b> to have the entire heat map visible in your browser window. <br><br>Select <b>Expand</b> to expand the heat map so that axis labels for all rows and columns are legible. You may have to scroll to see the entire heat map. 'Expand' can be computationally expensive for large datasets, and overrides the 'Text Size' setting."),
							),
							style="flex: 1;",
						),
						style="display: flex; gap: 0px; margin: 0px; align-items: flex-start;"
					),
					style="margin: 0px;"
				),

				config.HeatmapType.UI(ui.input_radio_buttons, make_inline=False, id="HeatmapType", label="Download File Type", choices=[".png", ".jpg"], inline=True),
				ui.download_button(id="DownloadHeatmap", label="Download Heatmap"),
			),

			# Settings pertaining to the dendrogram view.
			ui.panel_conditional(
				"input.MainTab === 'RowTab' || input.MainTab === 'ColumnTab'",
				# Define the Orientation of the dendrogram in the Tab
				config.Orientation.UI(ui.input_select,id="Orientation", label="Orientation", choices=["Top", "Bottom", "Left", "Right"]),
			),
			padding="10px",
			gap="20px",
			width="300px",
		),

		# Add the main interface tabs.
		MainTab(
			ui.nav_panel("Row Dendrogram", ui.output_plot("RowDendrogram", height="90vh"), value="RowTab"),
			ui.nav_panel("Column Dendrogram", ui.output_plot("ColumnDendrogram", height="90vh"), value="ColumnTab"),
			m_type=ui.output_image
		),
		height="86vh",
	)
)

app = App(app_ui, server)
