#
# Heatmapper
# Pairwise
#
# This file contains the ShinyLive application for Pairwise Heatmapper.
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
from scipy.spatial.distance import pdist, squareform
from Bio.PDB import PDBParser
from Bio import SeqIO
from pandas import DataFrame, read_csv, read_excel, read_table
from pathlib import Path

from shared import Table, Cache, NavBar, FileSelection


def server(input: Inputs, output: Outputs, session: Session):

	# Information about the Examples
	Info = {
		"example1.txt": "This example dataset represents pairwise distances between C-alpha atoms in ubiquitin (1ubq).",
		"example2.txt": "This example dataset was generated randomly.",
		"example3.txt": "This example dataset was generated randomly."
	}

	def HandleData(n, i):
		match Path(n).suffix:
			case ".csv": return read_csv(i)
			case ".xlsx": return read_excel(i)
			case ".pdb": return PDBMatrix(i)
			case ".fasta": return FASTAMatrix(i)
			case _: return read_table(i)
	DataCache = Cache("pairwise", HandleData)


	async def ParseData():
		"""
		@brief Returns a table containing the pairwise matrix.
		@returns	A DataFrame containing the data requested, formatted as a pairwise matrix, or
							an empty DataFrame if we're on Upload, but the user has not supplied a file.
		"""
		n = await DataCache.N(input)
		df = await DataCache.Load(input)

		if n is None: return DataFrame()
		match Path(n).suffix:
			case ".csv": df = ChartMatrix(df)
			case ".xlsx": df = ChartMatrix(df)
			case ".pdb": pass
			case ".fasta": pass
			case _: df = ChartMatrix(df)

		# Fix garbage data and return the resultant DataFrame.
		return df.fillna(0)


	def FASTAMatrix(file):
		"""
		@brief Computes the pairwise matrix from a FASTA file.
		@param file: The path to the FASTA File
		@returns a pairwise matrix.
		"""

		# Get information from the file
		records = list(SeqIO.parse(open(file), "fasta"))
		sequences = [str(record.seq) for record in records]
		column_names = [record.id for record in records]

		# Get our K-Mer value
		k = input.K()

		# Generate the value
		dictionary = {}
		for x, seq in enumerate(sequences):
			kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
			increment = 1 / len(kmers)
			for kmer in kmers:
					if kmer not in dictionary:
							dictionary[kmer] = [0.0] * len(sequences)
					dictionary[kmer][x] += increment
		frequencies = DataFrame.from_dict(dictionary, orient='index')

		# Calculate matrix
		if input.MatrixType() == "Distance":
			distances = pdist(frequencies.T, metric=input.DistanceMethod().lower())
			return DataFrame(squareform(distances), index=column_names, columns=column_names)
		else:
			return frequencies.corr(method=input.CorrelationMethod().lower())


	def PDBMatrix(file):
		"""
		@brief Generates a pairwise matrix from a PDB file
		@param file: The path to a PDB file (Or BytesIO file if applicable)
		@returns The pairwise matrix.
		"""

		parser = PDBParser()
		structure = parser.get_structure("protein", file)

		# Extract atomic coordinates
		coordinates = []
		for model in structure:
			for chain in model:
					if chain.id == input.Chain():
							for residue in chain:
									for atom in residue:
											coordinates.append(atom.coord)

		# Calculate matrix
		if input.MatrixType() == "Distance":
			distances = pdist(coordinates, metric=input.DistanceMethod().lower())
			return DataFrame(squareform(distances))
		else:
			return DataFrame(coordinates).corr(method=input.CorrelationMethod().lower())


	def ChartMatrix(df):
		"""
		@brief Generates a pairwise matrix from charts
		@param df:	The DataFrame containing the data. This can either be a chart
								containing {x,y,z} columns outlining each point on a row, with
								an optional name column (Any fourth column), a chart to which
								an explicit "Name" column is provided, to which the first row
								and column are assumed variable names for an existing matrix,
								or the default, where it is assumed that the chart is an
								unlabeled collection either of points, or an existing matrix.
		@returns A DataFrame containing the provided data as a pairwise matrix
		"""

		# If "Name" is found, its assumed to be the label for the points.
		if "Name" in df:
			point_names = df["Name"]

		# If explicit coordinates ar eprovided, use them, with the final column used as labels.
		if "x" in df.columns and "y" in df.columns and "z" in df.columns:
			coordinates = df[["x", "y", "z"]].values
			point_names = df[list(set(df.columns) - set(["x", "y", "z"]))[0]].values

		else:

			# If the first value is an integer, this is a distance matrix.
			try:
				int(df.iloc[0,0])
				coordinates = df.values

			# Otherwise, we assume the first row/column define the axis names.
			except ValueError:
				coordinates = df.iloc[:, 1:].values

			point_names = None

		# Calculate a distant matrix, and return it
		if input.MatrixType() == "Distance":
			distances = pdist(coordinates, metric=input.DistanceMethod().lower())
			return DataFrame(squareform(distances), index=point_names, columns=point_names)
		else:
			return DataFrame(coordinates, index=point_names, columns=point_names).corr(method=input.CorrelationMethod().lower())


	async def GenerateHeatmap():
		"""
		@brief Generates the Heatmap
		@returns The heatmap
		"""

		df = await ParseData()
		fig, ax = subplots()

		im = ax.imshow(df, cmap=input.ColorMap().lower(), interpolation=input.Interpolation().lower())

		# Visibility of features
		if "legend" in input.Features(): colorbar(im, ax=ax, label="Distance")

		if "y" in input.Features():
			ax.tick_params(axis="y", labelsize=input.TextSize())
			ax.set_yticks(range(len(df.columns)))
			ax.set_yticklabels(df.columns)
		else:
			ax.set_yticklabels([])

		if "x" in input.Features():
			ax.tick_params(axis="x", labelsize=input.TextSize())
			ax.set_xticks(range(len(df.columns)))
			ax.set_xticklabels(df.columns, rotation=90)
		else:
			ax.set_xticklabels([])

		# Annotate each cell with its value
		if "label" in input.Features():
			for i in range(df.shape[0]):
					for j in range(df.shape[1]):
							ax.text(j, i, '{:.2f}'.format(df.iloc[i, j]), ha='center', va='center', color='white')

		return ax


	@output
	@render.data_frame
	@reactive.event(input.Update, input.Reset, input.Example, input.File, ignore_none=False, ignore_init=False)
	async def LoadedTable(): return await DataCache.Load(input)


	@output
	@render.plot
	@reactive.event(input.Update, input.Reset, input.Example, input.File, input.MatrixType, input.TextSize, input.DistanceMethod, input.CorrelationMethod, input.Interpolation, input.ColorMap, input.Features, input.Chain, input.K, ignore_none=False, ignore_init=False)
	async def Heatmap(): return await GenerateHeatmap()

	@output
	@render.text
	def ExampleInfo(): return Info[input.Example()]


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

	NavBar("Pairwise"),

	ui.layout_sidebar(
		ui.sidebar(

			FileSelection(
				examples={"example1.txt": "Example 1", "example2.txt": "Example 2", "example3.txt": "Example 3"},
				types=[".csv", ".txt", ".xlsx", ".pdb", ".dat", ".fasta"]
			),

			# Specify Matrix Type
			ui.input_radio_buttons(id="MatrixType", label="Matrix Type", choices=["Distance", "Correlation"], selected="Distance", inline=True),

			# Customize the text size of the axes.
			ui.input_numeric(id="TextSize", label="Text Size", value=8, min=1, max=50, step=1),

			# https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html
			ui.panel_conditional(
				"input.MatrixType === 'Distance'",
				ui.input_select(id="DistanceMethod", label="Distance Method", choices=[
					"Braycurtis", "Canberra", "Chebyshev", "Cityblock", "Correlation", "Cosine", "Dice", "Euclidean", "Hamming", "Jaccard", "Jensenshannon", "Kulczynski1", "Mahalanobis", "Matching", "Minkowski", "Rogerstanimoto", "Russellrao", "Seuclidean", "Sokalmichener", "Sokalsneath", "Sqeuclidean", "Yule"], selected="Euclidean"),
			),

			# https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.corr.html
			ui.panel_conditional(
				"input.MatrixType === 'Correlation'",
				ui.input_select(id="CorrelationMethod", label="Correlation Method", choices=["Pearson", "Kendall", "Spearman"], selected="Pearson"),
			),

			# https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.imshow.html
			ui.input_select(id="Interpolation", label="Interpolation", choices=["None", "Antialiased", "Nearest", "Bilinear", "Bicubic", "Spline16", "Spline36", "Hanning", "Hamming", "Hermite", "Kaiser", "Quadric", "Catrom", "Gaussian", "Bessel", "Mitchell", "Sinc", "Lanczos", "Blackman"], selected="Nearest"),

			# Set the ColorMap used.
			ui.input_select(id="ColorMap", label="Color Map", choices=["Viridis", "Plasma", "Inferno", "Magma", "Cividis"], selected="Viridis"),

			# Customize what aspects of the heatmap are visible
			ui.input_checkbox_group(id="Features", label="Heatmap Features",
					choices={"x": "X Labels", "y": "Y Labels", "label": "Data Labels", "legend": "Legend"},
					selected=["legend"]),

			# Specify the PDB Chain
			ui.input_text("Chain", "PDB Chain", "A"),

			# Customize the K-mer to compute for FASTA sequences
			ui.input_numeric(id="K", label="K-Mer Length", value=3, min=3, max=5, step=1),

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
