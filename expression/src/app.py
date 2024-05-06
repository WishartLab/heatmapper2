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
from matplotlib.pyplot import figure, subplots, colorbar
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster import hierarchy
from scipy.stats import zscore
from pandas import DataFrame

from shared import Cache, NavBar, MainTab, FileSelection, Filter, ColumnType, TableOptions, ColorMap


def server(input, output, session):
	# Information about the Examples
	Info = {
		"example1.txt": "This example dataset is sample input retrieved from the website for the Ashley Lab Heatmap Builder.",
		"example2.txt": "This example dataset is sample input retrieved from an online tutorial by Yan Cui (ycui2@uthsc.edu).",
		"example3.txt": "This example dataset is retrieved from the online supplement to Eisen et al. (1998), which is a very well known paper about cluster analysis and visualization. The details of how the data was collected are outlined in the paper."
	}

	DataCache = Cache("expression")
	Data = reactive.value(DataFrame())
	Valid = reactive.value(False)

	@reactive.effect
	@reactive.event(input.SourceFile, input.File, input.Example, input.Reset)
	async def UpdateData(): Data.set((await DataCache.Load(input))); Valid.set(False)


	def GetData(): return Table.data_view() if Valid() else Data()


	def ProcessData(df):
		"""
		@brief Extracts the labels for each axis, and returns it alongside a DataFrame containing only the relevant data.
		@returns	A list containing the labels for the y axis, a list containing the labels for the x axis, and a
							DataFrame containing the loaded data without those two columns.
		"""

		name = input.NameColumn()
		if name not in df: return None, None , None

		# Drop the naming columns before linkage.
		data = df.drop(columns=Filter(df.columns, ColumnType.Name))
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
		matrix = hierarchy.linkage(data.values.T if invert else data.values, method=input.ClusterMethod().lower(), metric=input.DistanceMethod().lower())

		if progress is not None: progress.inc(message="Creating dendrogram...")
		dendrogram = hierarchy.dendrogram(matrix, ax=ax, orientation=orientation.lower())

		# If there are labels, sort them according to the dendrogram.
		if labels: labels = [labels[i] for i in dendrogram["leaves"]]

		# Add ticks depending on the orientation.
		if orientation == "Left" or orientation == "Right":
			ax.set_xticks([])
			ax.set_yticklabels(labels, fontsize=input.TextSize())
		else:
			ax.set_yticks([])
			ax.set_xticklabels(labels, fontsize=input.TextSize())

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

		GenerateDendrogram(data, ax, input.Orientation(), progress, labels, invert=invert)
		return fig


	@render.data_frame
	def Table(): Valid.set(True); return render.DataGrid(Data(), editable=True)


	@Table.set_patch_fn
	def UpdateTable(*, patch: render.CellPatch) -> render.CellValue:
		if input.Type() == "Integer": value = int(patch["value"])
		elif input.Type() == "Float": value = float(patch["value"])
		else: value = patch["value"]
		return value


	@render.plot
	def Heatmap(): 
		"""
		@brief Generates the Heatmap
		@returns The heatmap
		"""

		with ui.Progress() as p:
			p.inc(message="Reading input...")
			index_labels, x_labels, data = ProcessData(GetData())
			if data is None: return

			# Create a figure with a heatmap and associated dendrograms
			p.inc(message="Plotting...")
			fig = figure(figsize=(12, 10))

			gs = fig.add_gridspec(4, 2, height_ratios=[2, 8, 1, 1], width_ratios=[2, 8], hspace=0, wspace=0)

			# If we render the row dendrogram, we change the order of the index labels to match the dendrogram.
			# However, if we aren't rendering it, and thus row_dendrogram isn't defined, we simply assign df
			# To data, so the order changes when turning the toggle.
			if "row" in input.Features():
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
			if "col" in input.Features():
				ax_col = fig.add_subplot(gs[0, 1])
				col_dendrogram = GenerateDendrogram(data, ax_col, "Top", invert=True, progress=p)
				ax_col.axis("off")

			# Handle scaling
			if input.ScaleType() != "None": df = zscore(df, axis=1 if input.ScaleType() == "Row" else 0)
			
			# Render the heatmap.
			ax_heatmap = fig.add_subplot(gs[1, 1])

			colors = input.CustomColors() if input.Custom() else input.ColorMap().split()
			heatmap = ax_heatmap.imshow(
				df,
				cmap=LinearSegmentedColormap.from_list("ColorMap", colors, N=input.Bins()),
				interpolation=input.Interpolation().lower(),
				aspect="auto",
			)

			# If we render the Y axis.
			if "y" in input.Features():
				ax_heatmap.set_yticks(range(len(index_labels)))
				ax_heatmap.set_yticklabels(index_labels, fontsize=input.TextSize())
				ax_heatmap.yaxis.tick_right()
			else:
				ax_heatmap.set_yticklabels([])

			# If we render the X axis.
			if "x" in input.Features():
				ax_heatmap.set_xticks(range(len(x_labels)))
				ax_heatmap.set_xticklabels(x_labels, rotation=90, fontsize=input.TextSize())
			else:
				ax_heatmap.set_xticklabels([])

			# If we render the legend.
			if "legend" in input.Features():
				ax_cbar = fig.add_subplot(gs[3, 1])
				cbar = fig.colorbar(heatmap, cax=ax_cbar, orientation="horizontal")
				cbar.ax.tick_params(labelsize=input.TextSize())
			return fig


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


	@output
	@render.text
	@reactive.event(input.SourceFile, input.Example)
	def ExampleInfo(): return Info[input.Example()]


	@render.download(filename="table.csv")
	def DownloadTable(): yield GetData().to_string()


	@reactive.Effect
	def UpdateColumnSelection(): Filter(GetData().columns, ColumnType.Name, ui_element="NameColumn")


app_ui = ui.page_fluid(

	NavBar(),

	ui.layout_sidebar(
		ui.sidebar(

			FileSelection(
				examples={"example1.txt": "Example 1", "example2.txt": "Example 2", "example3.txt": "Example 3"},
				types=[".csv", ".txt", ".dat", ".tsv", ".tab", ".xlsx", ".xls", ".odf"],
				project="Expression"
			),

			TableOptions(),

			# Shared, Non-Table settings.
			ui.panel_conditional(
				"input.MainTab != 'TableTab'",
				# The column that holds names for the data.
				ui.input_select(id="NameColumn", label="Names", choices=[], multiple=False),

				# https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
				ui.input_select(id="ClusterMethod", label="Clustering Method", choices=["Single", "Complete", "Average", "Weighted", "Centroid", "Median", "Ward"], selected="Average"),

				# https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html#scipy.spatial.distance.pdist
				ui.input_select(id="DistanceMethod", label="Distance Method", choices=["Braycurtis", "Canberra", "Chebyshev", "Cityblock", "Correlation", "Cosine", "Dice", "Euclidean", "Hamming", "Jaccard", "Jensenshannon", "Kulczynski1", "Mahalanobis", "Matching", "Minkowski", "Rogerstanimoto", "Russellrao", "Seuclidean", "Sokalmichener", "Sokalsneath", "Sqeuclidean", "Yule"], selected="Euclidean"),

				# Customize the text size of the axes.
				ui.input_numeric(id="TextSize", label="Text Size", value=8, min=1, max=50, step=1),
			),

			# Settings pertaining to the Heatmap view.
			ui.panel_conditional(
				"input.MainTab === 'HeatmapTab'",
				# Define how the colors are scaled.
				ui.input_select(id="ScaleType", label="Scale Type", choices=["Row", "Column", "None"], selected="Row"),

				# https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.imshow.html
				ui.input_select(
					id="Interpolation", 
					label="Interpolation", 
					choices=["None", "Antialiased", "Nearest", "Bilinear", "Bicubic", "Spline16", "Spline36", "Hanning", "Hamming", "Hermite", "Kaiser", "Quadric", "Catrom", "Gaussian", "Bessel", "Mitchell", "Sinc", "Lanczos", "Blackman"], 
					selected="Nearest",
				),

				ColorMap(),
				
				ui.input_slider(id="Bins", label="Number of Bins", value=50, min=3, max=100, step=1),

				# Toggle rendering features. All are on by default.
				ui.input_checkbox_group(id="Features", label="Visibility",
					choices={"row": "Row Dendrogram", "col": "Column Dendrogram", "x": "X Labels", "y": "Y Labels", "legend": "Legend"},
					selected=["row", "col", "x", "y", "legend"])
			),

			# Settings pertaining to the dendrogram view.
			ui.panel_conditional(
				"input.MainTab === 'RowTab' || input.MainTab === 'ColumnTab'",
				# Define the Orientation of the dendrogram in the Tab
				ui.input_select(id="Orientation", label="Dendrogram Orientation", choices=["Top", "Bottom", "Left", "Right"], selected="Left"),
			),

			# Add the download buttons. You can download the heatmap by right clicking it :)
			ui.download_button("DownloadTable", "Download Table"),

			id="SidebarPanel",
		),

		# Add the main interface tabs.
		MainTab(
			ui.nav_panel("Row Dendrogram", ui.output_plot("RowDendrogram", height="90vh"), value="RowTab"),
			ui.nav_panel("Column Dendrogram", ui.output_plot("ColumnDendrogram", height="90vh"), value="ColumnTab"),
		),
	)
)

app = App(app_ui, server)
