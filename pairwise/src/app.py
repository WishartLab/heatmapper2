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


from shiny import App, reactive, render, ui
from matplotlib.pyplot import subplots, colorbar
from scipy.spatial.distance import pdist, squareform
from Bio.PDB import PDBParser
from Bio import SeqIO
from pandas import DataFrame
from pathlib import Path

from shared import Cache, NavBar, MainTab, Filter, ColumnType, FileSelection, TableValueUpdate


def server(input, output, session):

	# Information about the Examples
	Info = {
		"example1.txt": "This example dataset represents pairwise distances between C-alpha atoms in ubiquitin (1ubq).",
		"example2.txt": "This example dataset was generated randomly.",
		"example3.txt": "This example dataset was generated randomly.",
		"example4.fasta": "An example FASTA file.",
		"ala_phe_ala.pdb": "An example PDB file.",
	}

	def HandleData(n, i):
		suffix = Path(n).suffix
		if suffix == ".pdb": return PDBMatrix(n)
		elif suffix == ".fasta": return FASTAMatrix(n)
		else: return DataCache.DefaultHandler(n, i)
	DataCache = Cache("pairwise", HandleData)


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
		name_col = Filter(df.columns, ColumnType.Name, only_one=True, reject_unknown=True)
		if name_col: point_names = df[name_col]

		# If explicit coordinates ar eprovided, use them, with the final column used as labels.
		x_col = Filter(df.columns, ColumnType.X, only_one=True, reject_unknown=True)
		y_col = Filter(df.columns, ColumnType.Y, only_one=True, reject_unknown=True)
		z_col = Filter(df.columns, ColumnType.Z, only_one=True, reject_unknown=True)

		if x_col and y_col and z_col:
			coordinates = df[[x_col, y_col, z_col]].values
			point_names = df[list(set(df.columns) - set([x_col, y_col, z_col]))[0]].values

		else:

			# If the first value is an integer, this is a distance matrix.
			try:
				float(df.iloc[0,0])
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

		with ui.Progress() as p:

			p.inc(message="Reading input...")
			n, data = await DataCache.Load(input, return_n=True)
			if data.empty: return
			if Path(n).suffix not in [".pdb", ".fasta"]: df = ChartMatrix(data)
			else: df = data

			p.inc(message="Plotting...")
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
	@reactive.event(input.SourceFile, input.File, input.Example, input.Update, input.Reset)
	async def LoadedTable(): return await DataCache.Load(input)


	@output
	@render.plot
	@reactive.event(input.SourceFile, input.File, input.Example, input.Update, input.Reset, input.MatrixType, input.TextSize, input.DistanceMethod, input.CorrelationMethod, input.Interpolation, input.ColorMap, input.Features, input.Chain, input.K)
	async def Heatmap(): return await GenerateHeatmap()

	@output
	@render.text
	@reactive.event(input.SourceFile, input.Example)
	def ExampleInfo(): return Info[input.Example()]


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

	NavBar("Pairwise"),

	ui.layout_sidebar(
		ui.sidebar(

			FileSelection(
				examples={
				"example1.txt": "Example 1", 
				"example2.txt": "Example 2", 
				"example3.txt": "Example 3",
				"example4.fasta": "Example 4",
				"ala_phe_ala.pdb": "Example 5",
				},
				types=[".csv", ".txt", ".xlsx", ".xls", ".odf", ".pdb", ".dat", ".fasta"]
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
		MainTab(),
	)
)

app = App(app_ui, server)
