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


from shiny import App, reactive, render, ui, types
from matplotlib.pyplot import subplots, colorbar
from scipy.spatial.distance import pdist, squareform
from matplotlib.colors import LinearSegmentedColormap
from Bio.PDB import PDBParser
from Bio import SeqIO
from pandas import DataFrame
from tempfile import NamedTemporaryFile
from io import BytesIO

from shared import Cache, NavBar, MainTab, Filter, ColumnType, FileSelection, TableOptions, Colors, DistanceMethods, InterpolationMethods, InitializeConfig, Error, Update

try:
	from user import config
except ImportError:
	from config import config


def server(input, output, session):

	# Information about the Examples
	Info = {
		"example1.txt": "This example dataset represents pairwise distances between C-alpha atoms in ubiquitin (1ubq).",
		"example2.txt": "This example dataset was generated randomly.",
		"example3.txt": "This example dataset was generated randomly.",
		"example4.fasta": "An example FASTA file.",
		"ala_phe_ala.pdb": "An example PDB file.",
	}

	def HandleData(path):
		suffix = path.suffix
		if suffix == ".pdb": return PDBMatrix(path.resolve())
		elif suffix == ".fasta": return FASTAMatrix(path.resolve())
		else: return ChartMatrix(DataCache.DefaultHandler(path))


	DataCache = Cache("pairwise", HandleData)
	Data = reactive.value(None)
	Valid = reactive.value(False)

	InitializeConfig(config, input)


	# We add Matrix and Method as they are calculated in the Matrix handlers.
	@reactive.effect
	@reactive.event(input.SourceFile, input.File, input.Example, input.Reset)
	async def UpdateData(): Data.set((await DataCache.Load(input, p=ui.Progress()))); Valid.set(False)


	def GetData(): return Table.data_view() if Valid() else Data()


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
		k = config.K()

		# Generate the value
		dictionary = {}
		for x, seq in enumerate(sequences):
			kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
			increment = 1 / len(kmers)
			for kmer in kmers:
					if kmer not in dictionary:
							dictionary[kmer] = [0.0] * len(sequences)
					dictionary[kmer][x] += increment
		return DataFrame.from_dict(dictionary, orient='index', columns=column_names)


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
					if chain.id == config.Chain():
							for residue in chain:
									for atom in residue:
											coordinates.append(atom.coord)
		return DataFrame(coordinates)


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
		name_col = Filter(df.columns, ColumnType.Name)
		if name_col: point_names = df[name_col]

		# If explicit coordinates are provided, use them, with the final column used as labels.
		x_col = Filter(df.columns, ColumnType.X)
		y_col = Filter(df.columns, ColumnType.Y)
		z_col = Filter(df.columns, ColumnType.Z)

		if x_col != y_col and y_col != z_col:
			coordinates = df[[x_col, y_col, z_col]].values
			point_names = df[list(set(df.columns) - set([x_col, y_col, z_col]))[0]].values

		else:

			# If the first value is an integer, this is a distance matrix.
			try:
				float(df.iloc[0,0])
				coordinates = df.values
				point_names = df.columns

			# Otherwise, we assume the first row/column define the axis names.
			except ValueError:
				coordinates = df.iloc[:, 1:].values
				point_names = df.columns[1:]
		return DataFrame(coordinates, columns=point_names, index=point_names).T


	@output
	@render.data_frame
	def Table():
		df = Data()
		if df is None or len(df.columns) == 0: return

		if df.columns[0] == 0:
			ui.modal_show(ui.modal("The provided input format cannot be rendered",
			title="Table cannot be rendered", easy_close=True, footer=None))
		else:
			Valid.set(True)
			return render.DataGrid(Data(), editable=True)


	@Table.set_patch_fn
	def UpdateTable(*, patch: render.CellPatch) -> render.CellValue:
		if config.Type() == "Integer": value = int(patch["value"])
		elif config.Type() == "Float": value = float(patch["value"])
		else: value = patch["value"]
		return value


	def GenerateHeatmap():

		inputs = [
			input.File() if input.SourceFile() == "Upload" else input.Example(),
			config.DistanceMethod() if config.MatrixType() == "Distance" else config.CorrelationMethod(),
			input.CustomColors() if config.Custom() else config.ColorMap().split(),
			config.Interpolation(),
			config.Bins(),
			config.TextSize(),
			config.Features(),
			config.DPI(),
		]

		if not DataCache.In(inputs):
			with ui.Progress() as p:
				p.inc(message="Reading input...")
				data = GetData()
				if data is None or len(data.index) == 0: return

				p.inc(message="Calculating...")
				try:
					# Calculate matrix
					if config.MatrixType() == "Distance":
						metric = config.DistanceMethod().lower()
						distances = pdist(data, metric=metric)
						df = DataFrame(squareform(distances), columns=data.index, index=data.index)
					else:
						method = config.CorrelationMethod().lower()
						df = data.T.corr(method=method)
				except Exception:
					Error("Could not compute matrix. Ensure your input data is correct!")
					return

				p.inc(message="Plotting...")
				fig, ax = subplots()

				colors = input.CustomColors() if config.Custom() else config.ColorMap().split()
				interpolation = config.Interpolation().lower()
				im = ax.imshow(
					df,
					cmap=LinearSegmentedColormap.from_list("ColorMap", colors, N=config.Bins()),
					interpolation=interpolation,
					aspect="equal",
				)

				text_size = config.TextSize()

				# Visibility of features
				if "legend" in config.Features():
					cbar = colorbar(im, ax=ax, label="Distance")
					cbar.ax.tick_params(labelsize=text_size)

				if "y" in config.Features():
					ax.tick_params(axis="y", labelsize=text_size)
					ax.set_yticks(range(len(df.columns)))
					ax.set_yticklabels(df.columns)
				else:
					ax.set_yticklabels([])

				if "x" in config.Features():
					ax.tick_params(axis="x", labelsize=text_size)
					ax.set_xticks(range(len(df.columns)))
					ax.set_xticklabels(df.columns, rotation=90)
				else:
					ax.set_xticklabels([])

				# Annotate each cell with its value
				if "label" in config.Features():
					for i in range(df.shape[0]):
							for j in range(df.shape[1]):
									ax.text(j, i, '{:.2f}'.format(df.iloc[i, j]), ha='center', va='center', color='white')

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


	@output
	@render.text
	@reactive.event(input.SourceFile, input.Example)
	def ExampleInfo(): return Info[input.Example()]


	@render.download(filename="table.csv")
	def DownloadTable(): yield GetData().to_string()


	@render.download(filename="heatmap.png")
	def DownloadHeatmap():
		yield DataCache.Get([
			input.File() if input.SourceFile() == "Upload" else input.Example(),
			config.DistanceMethod() if config.MatrixType() == "Distance" else config.CorrelationMethod(),
			input.CustomColors() if config.Custom() else config.ColorMap().split(),
			config.Interpolation(),
			config.Bins(),
			config.TextSize(),
			config.Features(),
			config.DPI(),
		])


	@render.ui
	def Method():
		if config.MatrixType() == "Distance":
			return config.DistanceMethod.UI(ui.input_select, id="DistanceMethod", label="Distance", choices=DistanceMethods)
		elif config.MatrixType() == "Correlation":
			return config.CorrelationMethod.UI(ui.input_select, id="CorrelationMethod", label="Correlation", choices=["Pearson", "Kendall", "Spearman"])


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
				examples={
				"example1.txt": "Example 1",
				"example2.txt": "Example 2",
				"example3.txt": "Example 3",
				"example4.fasta": "Example 4",
				"ala_phe_ala.pdb": "Example 5",
				},
				types=[".csv", ".txt", ".dat", ".tsv", ".tab", ".xlsx", ".xls", ".odf", ".pdb", ".dat", ".fasta"],
				project="Pairwise"
			),

			TableOptions(config),

			ui.panel_conditional(
				"input.MainTab != 'TableTab'",

				Update(),

				ui.HTML("<b>Heatmap</b>"),
				config.MatrixType.UI(ui.input_select, id="MatrixType",  label="Matrix",  choices=["Distance", "Correlation"]),
				config.TextSize.UI(ui.input_numeric, id="TextSize", label="Text Size", min=1, max=20, step=1),
				config.Interpolation.UI(ui.input_select, id="Interpolation", label="Inter", choices=InterpolationMethods),
				config.Chain.UI(ui.input_text, id="Chain", label="Chain"),
				config.K.UI(ui.input_numeric, id="K", label="K-Mer", min=3, max=5, step=1),
				ui.output_ui("Method"),

				ui.HTML("<b>Colors</b>"),
				ui.output_ui("Color"),
				config.Custom.UI(ui.input_checkbox, id="Custom", label="Custom"),
				config.Bins.UI(ui.input_slider, id="Bins", label="Number", min=3, max=100, step=1),

				ui.HTML("<b>Image Settings</b>"),
				config.Size.UI(ui.input_numeric, id="Size", label="Size", min=1),
				config.DPI.UI(ui.input_numeric, id="DPI", label="DPI", min=1),

				ui.HTML("<b>Features</b>"),
				config.Features.UI(ui.input_checkbox_group,
					make_inline=False, id="Features", label=None,
					choices={"x": "X Labels", "y": "Y Labels", "label": "Data Labels", "legend": "Legend"},
				),

				ui.download_button(id="DownloadHeatmap", label="Download"),
			),
			padding="1rem",
			width="255px",
		),

		# Add the main interface tabs.
		MainTab(m_type=ui.output_image),
		height="90vh",
	)
)

app = App(app_ui, server)
