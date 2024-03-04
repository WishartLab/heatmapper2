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


from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from matplotlib.pyplot import figure, subplots, colorbar
from matplotlib.colors import Normalize
from scipy.cluster import hierarchy
from pandas import DataFrame

from shared import Table, Cache, NavBar, FileSelection

def server(input: Inputs, output: Outputs, session: Session):
	# Information about the Examples
	Info = {
		"example1.txt": "This example dataset is sample input retrieved from the website for the Ashley Lab Heatmap Builder.",
		"example2.txt": "This example dataset is sample input retrieved from an online tutorial by Yan Cui (ycui2@uthsc.edu).",
		"example3.txt": "This example dataset is retrieved from the online supplement to Eisen et al. (1998), which is a very well known paper about cluster analysis and visualization. The details of how the data was collected are outlined in the paper."
	}

	DataCache = Cache("expression")


	async def ProcessData():
		"""
		@brief Extracts the labels for each axis, and returns it alongside a DataFrame containing only the relevant data.
		@returns	A list containing the labels for the y axis, a list containing the labels for the x axis, and a
							DataFrame containing the loaded data without those two columns.
		"""
		df = await DataCache.Load(input)

		names = ["NAME", "ORF", "UNIQID"]

		for name in names:
			if name in df.columns:
				index_labels = df[name]
				break

		# Drop the naming columns before linkage.
		data = df.drop(columns=[col for col in names if col in df.columns])
		x_labels = ["X" + name if list(data.columns).count(name) == 1 else "X" + name + f".{i+1}" for i, name in enumerate(data.columns)]

		return list(index_labels), x_labels, data


	def GenerateDendrogram(data, ax, orientation, labels = [], invert=False):
		"""
		@brief General dendrogram generator.
		@param data: The DataFrame that contains the data to generate the dendrogram from.
		@param ax: The MatPlotLib Axis to assign tick marks to
		@param orientation: What orientation we should set the dendrogram to be. Can be "Left", "Right", "Top", or "Bottom"
		@param labels: An optional list of labels to add the dendrogram, labelling the X axis on Left/Right, and the Y on Top/Bottom
		@param invert: Whether to invert the DataFrame to generate Columnar dendrograms.
		@returns The dendrogram, mostly useful to aligning the Heatmap to the new ordering.
		"""

		matrix = hierarchy.linkage(data.values.T if invert else data.values, method=input.ClusterMethod().lower(), metric=input.DistanceMethod().lower())
		dendrogram = hierarchy.dendrogram(matrix, ax=ax, orientation=orientation.lower())

		# If there are labels, sort them according to the dendrogram.
		if labels: labels = [labels[i] for i in dendrogram["leaves"]]

		# Add ticks depending on the orientation.
		match orientation:
			case "Left" | "Right":
				ax.set_xticks([])
				ax.set_yticklabels(labels, fontsize=input.TextSize())
			case "Top" | "Bottom":
				ax.set_yticks([])
				ax.set_xticklabels(labels, fontsize=input.TextSize())

		return dendrogram


	async def GenerateHeatmap():
		"""
		@brief Generates the Heatmap
		@returns The heatmap
		"""

		index_labels, x_labels, data = await ProcessData()

		# Create a figure with a heatmap and associated dendrograms
		fig = figure(figsize=(12, 10))
		gs = fig.add_gridspec(4, 2, height_ratios=[2, 8, 1, 1], width_ratios=[2, 8], hspace=0, wspace=0)

		# If we render the row dendrogram, we change the order of the index labels to match the dendrogram.
		# However, if we aren't rendering it, and thus row_dendrogram isn't defined, we simply assign df
		# To data, so the order changes when turning the toggle.
		if "row" in input.Features():
			ax_row = fig.add_subplot(gs[1, 0])
			row_dendrogram = GenerateDendrogram(data, ax_row, "Left")
			ax_row.axis("off")
			index_labels = [index_labels[i] for i in row_dendrogram["leaves"]]
			df = data.iloc[row_dendrogram["leaves"]]
		else:
			df = data

		# If we render the column dendrogram.
		if "col" in input.Features():
			ax_col = fig.add_subplot(gs[0, 1])
			col_dendrogram = GenerateDendrogram(data, ax_col, "Top", invert=True)
			ax_col.axis("off")

		# Handle normalization
		match input.ScaleType():
			case "Row": df = df.div(df.max(axis=1), axis=0)
			case "Column": df = df.div(df.max(axis=0), axis=1)

		# Render the heatmap.
		ax_heatmap = fig.add_subplot(gs[1, 1])
		heatmap = ax_heatmap.imshow(
			df,
			cmap=input.ColorMap().lower(),
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

		return fig


	@output
	@render.data_frame
	@reactive.event(input.Update, input.Reset, input.Example, input.File, ignore_none=False, ignore_init=False)
	async def LoadedTable(): return await DataCache.Load(input)

	@output
	@render.plot
	@reactive.event(input.Update, input.Reset, input.Example, input.File, input.ClusterMethod, input.DistanceMethod, input.TextSize, input.ScaleType, input.Interpolation, input.ColorMap, input.Features, ignore_none=False, ignore_init=False)
	async def Heatmap(): return await GenerateHeatmap()


	@output
	@render.plot
	@reactive.event(input.Update, input.Reset, input.Example, input.File, input.Orientation, input.ClusterMethod, input.DistanceMethod, input.TextSize, ignore_none=False, ignore_init=False)
	async def RowDendrogram():
		index_labels, _, data = await ProcessData()

		fig = figure(figsize=(12, 10))
		ax = fig.add_subplot(111)

		ax.spines["top"].set_visible(False)
		ax.spines["right"].set_visible(False)
		ax.spines["bottom"].set_visible(False)
		ax.spines["left"].set_visible(False)

		GenerateDendrogram(data, ax, input.Orientation(), index_labels)
		return fig


	@output
	@render.plot
	@reactive.event(input.Update, input.Reset, input.Example, input.File, input.Orientation, input.ClusterMethod, input.DistanceMethod, input.TextSize, ignore_none=False, ignore_init=False)
	async def ColumnDendrogram():
		_, x_labels, data = await ProcessData()

		fig = figure(figsize=(12, 10))
		ax = fig.add_subplot(111)

		ax.spines["top"].set_visible(False)
		ax.spines["right"].set_visible(False)
		ax.spines["bottom"].set_visible(False)
		ax.spines["left"].set_visible(False)

		GenerateDendrogram(data, ax, input.Orientation(), x_labels, invert=True)
		return fig


	@output
	@render.text
	def ExampleInfo(): return Info[input.Example()]


	@render.download(filename="table.csv")
	async def DownloadTable(): df = await LoadData(); yield df.to_string()

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

	NavBar("Expression"),

	ui.layout_sidebar(
		ui.sidebar(

			FileSelection(
				examples={"example1.txt": "Example 1", "example2.txt": "Example 2", "example3.txt": "Example 3"},
				types=[".csv", ".txt", ".xlsx", ".pdb", ".dat"]
			),

			# https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
			ui.input_select(id="ClusterMethod", label="Clustering Method", choices=["Single", "Complete", "Average", "Weighted", "Centroid", "Median", "Ward"], selected="Average"),

			# https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html#scipy.spatial.distance.pdist
			ui.input_select(id="DistanceMethod", label="Distance Method", choices=["Braycurtis", "Canberra", "Chebyshev", "Cityblock", "Correlation", "Cosine", "Dice", "Euclidean", "Hamming", "Jaccard", "Jensenshannon", "Kulczynski1", "Mahalanobis", "Matching", "Minkowski", "Rogerstanimoto", "Russellrao", "Seuclidean", "Sokalmichener", "Sokalsneath", "Sqeuclidean", "Yule"], selected="Euclidean"),

			# Customize the text size of the axes.
			ui.input_numeric(id="TextSize", label="Text Size", value=8, min=1, max=50, step=1),

			# Settings pertaining to the Heatmap view.
			ui.panel_conditional(
				"input.MainTab === 'Interactive'",
				ui.br(),

				# Define how the colors are scaled.
				ui.input_select(id="ScaleType", label="Scale Type", choices=["Row", "Column", "None"], selected="None"),

				# https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.imshow.html
				ui.input_select(id="Interpolation", label="Interpolation", choices=["None", "Antialiased", "Nearest", "Bilinear", "Bicubic", "Spline16", "Spline36", "Hanning", "Hamming", "Hermite", "Kaiser", "Quadric", "Catrom", "Gaussian", "Bessel", "Mitchell", "Sinc", "Lanczos", "Blackman"], selected="Nearest"),

				# Set the ColorMap used.
				ui.input_select(id="ColorMap", label="Color Map", choices=["Viridis", "Plasma", "Inferno", "Magma", "Cividis"], selected="Viridis"),

				# Toggle rendering features. All are on by default.
				ui.input_checkbox_group(id="Features", label="Heatmap Features",
					choices={"row": "Row Dendrogram", "col": "Column Dendrogram", "x": "X Labels", "y": "Y Labels", "legend": "Legend"},
					selected=["row", "col", "x", "y", "legend"])
			),

			# Settings pertaining to the dendrogram view.
			ui.panel_conditional(
				"input.MainTab === 'Row' || input.MainTab === 'Column'",
				ui.br(),

				# Define the Orientation of the dendrogram in the Tab
				ui.input_select(id="Orientation", label="Dendrogram Orientation", choices=["Top", "Bottom", "Left", "Right"], selected="Left"),
			),

			# Add the download buttons. You can download the heatmap by right clicking it :)
			ui.download_button("DownloadTable", "Download Table"),

			id="SidebarPanel",
		),

		# Add the main interface tabs.
		ui.navset_tab(
				ui.nav_panel("Interactive", ui.output_plot("Heatmap", height="90vh"), value="Interactive"),
				ui.nav_panel("Row Dendrogram", ui.output_plot("RowDendrogram", height="90vh"), value="Row"),
				ui.nav_panel("Column Dendrogram", ui.output_plot("ColumnDendrogram", height="90vh"), value="Column"),
				Table,
				id="MainTab"
		),
	)
)

app = App(app_ui, server)