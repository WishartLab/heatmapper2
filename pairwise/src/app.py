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
from matplotlib.pyplot import subplots, colorbar, style
from scipy.spatial.distance import pdist, squareform
from scipy.interpolate import griddata
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from Bio.PDB import PDBParser
from Bio import SeqIO
#from Bio.SeqIO.FastaIO import SimpleFastaParser
from pandas import DataFrame
from tempfile import NamedTemporaryFile
from io import BytesIO
from numpy import arange, zeros_like, meshgrid, array, column_stack, linspace, min as n_min, concatenate

from shared import Cache, NavBar, MainTab, Filter, ColumnType, FileSelection, TableOptions, Colors, DistanceMethods, InterpolationMethods, InitializeConfig, Error, Update, Msg, File

try:
	from user import config
except ImportError:
	from config import config


def server(input, output, session):

	# Information about the Examples
	Info = {
		"example1.txt": "Input type: txt<br>Contents: Pairwise distances between C-alpha atoms in ubiquitin (1ubq).",
		"example2.txt": "Input type: txt<br>Contents: This example dataset was generated randomly.",
		"example3.txt": "Input type: txt<br>Contents: This example dataset was generated randomly.",
		"example4.fasta": "Input type: FASTA<br>Contents: ",
		"ala_phe_ala.pdb": "Input type: PDB<br>Contents: ",
		"example6.txt": "Input type: txt<Br>Contents: Randomly generated data with x, y, and z columns."
	}

	def HandleData(path, p=None):
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
	async def UpdateData():
		Data.set((await DataCache.Load(input, p=ui.Progress())));
		Valid.set(False)
		DataCache.Invalidate(File(input))


	def GetData(): return Table.data_view() if Valid() else Data()


	def HashString():
		inputs = [
			File(input),
			config.DistanceMethod() if config.MatrixType() == "Distance" else config.CorrelationMethod(),
			input.CustomColors() if config.Custom() else config.ColorMap().split(),
			config.Interpolation(),
			config.Bins(),
			config.TextSize(),
			config.Features(),
			config.DPI(),
			config.Elevation(),
			input.mode(),
		]
		if config.Elevation() != 90: inputs.extend([config.Rotation(), config.HeightMatrix(), config.Zoom(), config.InterpolationLevels(), config.MinScale(), config.Opacity()])
		return inputs


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

		point_names = None

		# If "Name" is found, its assumed to be the label for the points.
		name_col = Filter(df.columns, ColumnType.Name)
		if name_col:
			point_names = df[name_col]
			df = df.drop(name_col, inplace=False, axis=1)

		#x_col = Filter(df.columns, ColumnType.X)
		#y_col = Filter(df.columns, ColumnType.Y)
		#z_col = Filter(df.columns, ColumnType.Z)

		#if x_col != y_col and y_col != z_col:
		#	if name_col is not None:
		#		return df[[name_col, x_col, y_col, z_col]]
		#	return df[[x_col, y_col, z_col]]

		# If the first value is an integer, this is a distance matrix.
		try:
			float(df.iloc[0,0])
			coordinates = df.values
			if point_names is None: point_names = df.columns
			columns = df.columns

		# Otherwise, we assume the first row/column define the axis names.
		except ValueError:
			coordinates = df.iloc[:, 1:].values
			if point_names is None: point_names = df.columns[1:]
			columns = df.columns[1:]
		try:
			return DataFrame(coordinates, index=point_names, columns=columns)
		except Exception:
			return DataFrame(coordinates)


	@output
	@render.data_frame
	def Table():
		df = Data()
		if df is None or len(df.columns) == 0: return
		if df.columns[0] == 0:
			Error("The provided input format cannot be rendered")
		else:
			Valid.set(True)
			return render.DataGrid(df, editable=True)


	@Table.set_patch_fn
	def UpdateTable(*, patch: render.CellPatch) -> render.CellValue:
		if config.Type() == "Integer": value = int(patch["value"])
		elif config.Type() == "Float": value = float(patch["value"])
		else: value = patch["value"]

		# If changes are made, invalidate all the cached objects relying on it.
		DataCache.Invalidate(File(input))

		return value


	def GenerateMatrix(data, value):
		'''
		@param data: Pandas df
		'''
		name_col = Filter(data.columns, ColumnType.Name)
		if name_col is not None:
			names = data[name_col]
			data = data.drop(name_col, inplace=False, axis=1)
		else:
			names = data.index

		try:
			# Calculate matrix
			if value == "Distance":
				metric = config.DistanceMethod().lower()
				distances = pdist(data, metric=metric)
				return DataFrame(squareform(distances), columns=names, index=names)
			else:
				method = config.CorrelationMethod().lower()
				return data.T.corr(method=method)
		except Exception:
			Error("Could not compute matrix. Ensure your input data is correct!")
			return None
		

	def HeatmapCube(df, cmap, p):
		fig, ax = subplots(subplot_kw={"projection": "3d"})

		x, y, z = Filter(df.columns, ColumnType.X), Filter(df.columns, ColumnType.Y), Filter(df.columns, ColumnType.Z)
		if not x or not y or not z:
			Error("An X, Y, and Z column are needed to compute a Cube Visualization!")

		name_col = Filter(df.columns, ColumnType.Name)
		if name_col is not None:
			dxy = GenerateMatrix(df[[name_col, x, y]], config.MatrixType())
			dxz = GenerateMatrix(df[[name_col, x, z]], config.MatrixType())
			dyz = GenerateMatrix(df[[name_col, y, z]], config.MatrixType())
		else:
			dxy = GenerateMatrix(df[[x, y]], config.MatrixType())
			dxz = GenerateMatrix(df[[x, z]], config.MatrixType())
			dyz = GenerateMatrix(df[[y, z]], config.MatrixType())

		d = concatenate([dxy.values.flatten(), dxz.values.flatten(), dyz.values.flatten()])
		norm = Normalize(vmin=d.min(), vmax=d.max())

		# Create mesh grids
		ix, iy = dxy.shape
		x, y = meshgrid(arange(ix), arange(iy))

		ix, iz = dxz.shape
		xz, zz = meshgrid(arange(ix), arange(iz))

		iy, iz = dyz.shape
		yz, zz2 = meshgrid(arange(iy), arange(iz))

		# Plot surfaces on the respective planes
		ax.view_init(elev=config.Elevation(), azim=config.Rotation())
		#ax.tick_params(axis='x', labelrotation=90)
		ax.set_box_aspect(None, zoom=config.Zoom())

		im = ax.plot_surface(x, y, zeros_like(x), facecolors=cmap(norm(dxy.values)), shade=False)
		ax.plot_surface(xz, zeros_like(xz), zz, facecolors=cmap(norm(dxz.values)), shade=False)
		ax.plot_surface(zeros_like(yz), yz, zz2, facecolors=cmap(norm(dyz.values)), shade=False)

		return fig, ax, im, norm, d, dxy


	def Heatmap2D(df, cmap, p):
		"""
		@brief Generate a 2D heatmap
		@param df: The dataframe to render
		@param cmap: The LinearSegmentedColorMap to use
		@param p: The progress indicator.
		"""
		p.inc(message="Generating...")
		fig, ax = subplots()
		interpolation = config.Interpolation().lower()
		im = ax.imshow(df, cmap=cmap, interpolation=interpolation, aspect="equal")
		return fig, ax, im


	def Heatmap3D(df, data, cmap, p):
		"""
		@brief Generate a 3D heatmap
		@param df: The dataframe to render
		@param cmap: The LinearSegmentedColorMap to use
		@param p: The progress indicator.
		"""
		fig, ax = subplots(subplot_kw={"projection": "3d"})

		cached = [
			File(input),
			config.DistanceMethod() if config.MatrixType() == "Distance" else config.CorrelationMethod(),
			input.CustomColors() if config.Custom() else config.ColorMap().split(),
			config.Bins(),
			config.HeightMatrix(),
			config.InterpolationLevels(),
			config.MinScale()
		]

		if not DataCache.In(cached):
			p.inc(message="Generating...")
			if config.HeightMatrix()  != config.MatrixType():
				df_height = GenerateMatrix(data, config.HeightMatrix())
			else: df_height = df

			z = df_height.values.flatten()
			if config.MinScale():
				z += abs(n_min(z))

			color_array = df.values.flatten()
			norm = Normalize(vmin=color_array.min(), vmax=color_array.max())
			c = norm(color_array)

			length = df_height.shape[0]
			x, y = meshgrid(arange(length), arange(length))
			x = x.ravel()
			y = y.ravel()
			width = depth = 1

			level = config.InterpolationLevels()
			if level < 1:
				level = 1
			elif level > 10:
				level = 10
			if level != 1:
				p.inc(message="Interpolating...")
				length_new = df_height.shape[0] * level
				space = linspace(0, length-1, length_new)
				x_grid, y_grid = meshgrid(space, space)
				x_new = x_grid.ravel()
				y_new = y_grid.ravel()

				points = column_stack((x, y))
				points_new = column_stack((x_new, y_new))

				z = griddata(points, z, points_new, method='cubic')
				c = griddata(points, c, points_new, method='cubic')

				x, y = x_new, y_new
				width /= level
				depth /= level
			DataCache.Store((x, y, width, depth, z, c, norm), cached)
		x, y, width, depth, z, c, norm = DataCache.Get(cached)

		ax.view_init(elev=config.Elevation(), azim=config.Rotation())
		ax.set_box_aspect(None, zoom=config.Zoom())
		im = ax.bar3d(x, y, zeros_like(x), width, depth, z, color=cmap(c), alpha=config.Opacity())

		return fig, ax, im, norm, z



	def GenerateHeatmap():
		inputs = HashString()

		if not DataCache.In(inputs):
			with ui.Progress() as p:
				p.inc(message="Reading input...")
				data = GetData()
				if data is None or len(data.index) == 0: return

				p.inc(message="Calculating...")
				if config.HeightMatrix() == "Cube":
					df = data
				else:
					df = GenerateMatrix(data, config.MatrixType())
				if df is None: return

				color = input.mode()
				colors = input.CustomColors() if config.Custom() else config.ColorMap().split()
				cmap = LinearSegmentedColormap.from_list("ColorMap", colors, N=config.Bins())

				with style.context('dark_background' if color == "dark" else "default"):
					rotation = config.Rotation()
					elevation = config.Elevation()

					if config.HeightMatrix() == "Cube":
						fig, ax, im, norm, z, df = HeatmapCube(df, cmap, p)
						d3 = True
					elif elevation != 90:
						fig, ax, im, norm, z = Heatmap3D(df, data, cmap, p)
						d3 = True
					else:
						fig, ax, im = Heatmap2D(df, cmap, p)
						d3 = False


					p.inc(message="Plotting...")
					text_size = config.TextSize()

					# Visibility of features
					if "legend" in config.Features():
						if not d3:
							cbar = colorbar(im, ax=ax, label=config.MatrixType())
						else:
							mappable = ScalarMappable(cmap=cmap, norm=norm)
							mappable.set_array(z)
							cbar = colorbar(mappable, ax=ax, label='Value', orientation='vertical')
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

					if d3:
						if "z" in config.Features():
							ax.tick_params(axis="z", labelsize=text_size)
							ax.set_zticks(range(len(df.columns)))
							ax.set_zticklabels(df.columns)
						else:
							ax.set_zticklabels([])

					# Annotate each cell with its value
					if "label" in config.Features():
						for i in range(df.shape[0]):
								for j in range(df.shape[1]):
									if not d3:
										ax.text(j, i, '{:.2f}'.format(df.iloc[i, j]), ha='center', va='center', color='white')
									else:
										ax.text(j, i, z[i * df.shape[1] + j], '{:.2f}'.format(df.iloc[i, j]), ha='center', va='center', color='black')

					# catch invalid DPI values
					if config.DPI() < 5:
						dpi = 5
					else:
						dpi = config.DPI()

					b = BytesIO()
					fig.savefig(b, format="png", dpi=dpi, bbox_inches="tight")
					b.seek(0)
					DataCache.Store(b.read(), inputs)

		b = DataCache.Get(inputs)
		with NamedTemporaryFile(delete=False, suffix=".png") as temp:
			temp.write(b)
			temp.close()
			img: types.ImgData = {"src": temp.name, "width": f"{config.Size()}px"}
			return img


	@output
	@render.image(delete_file=True)
	def Heatmap(): return GenerateHeatmap()


	@output
	@render.image(delete_file=True)
	@reactive.event(input.Update)
	def HeatmapReactive():
		return GenerateHeatmap()


	@reactive.effect
	@reactive.event(input.ExampleInfo)
	def ExampleInfo():
		Msg(ui.HTML(Info[input.Example()]))


	@render.download(filename="table.csv")
	def DownloadTable(): yield GetData().to_string()


	@render.download(filename="heatmap.png")
	def DownloadHeatmap(): yield DataCache.Get(HashString())


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

	ui.panel_title(title=None, window_title="Pairwise"),
	NavBar(),

	ui.layout_sidebar(
		ui.sidebar(

			FileSelection(
				examples={
				"example1.txt": "Ex1: Matrix",
				"example2.txt": "Ex2: RandomData",
				"example3.txt": "Ex3: Matrix",
				"example4.fasta": "Ex4: FASTA",
				"ala_phe_ala.pdb": "Ex5: PDB",
				"example6.txt": "Ex6: CubeMatrix",
				},
				types=[".csv", ".txt", ".dat", ".tsv", ".tab", ".xlsx", ".xls", ".odf", ".pdb", ".dat", ".fasta"],
				project="Pairwise"
			),

			TableOptions(config),

			ui.panel_conditional(
				"input.MainTab != 'TableTab'",

				Update(),

				ui.HTML("<b>Heatmap</b>"),
				config.MatrixType.UI(ui.input_select, id="MatrixType",	label="Matrix",	choices=["Distance", "Correlation"], tooltip="Visualize either the distance or correlation between values. Based on the matrix type, you can further select a distance calculation method or correlation calculation method below."),
				ui.output_ui("Method"),

				config.TextSize.UI(ui.input_numeric, id="TextSize", label="Text", min=1, max=20, step=1, tooltip="Change the text size of all axis labels. Axis labels can be toggled on and off in the 'Features' section at the bottom of this sidebar."),
				config.Interpolation.UI(ui.input_select, id="Interpolation", label="Inter", choices=InterpolationMethods, conditional="input.Elevation === 90", tooltip="Specify an interpolation algorithm to apply to the figure. This can cause values to bleed together and appear smoother."),
				config.Chain.UI(ui.input_text, id="Chain", label="Chain", tooltip="This setting only applies if a PDB file is used. Select a chain within the PDB file to display."),
				config.K.UI(ui.input_numeric, id="K", label="K-Mer", min=3, max=5, step=1, tooltip="This setting only applies if a FASTA file is used. Specify the length of K-Mer (3, 4, or 5) to use for alignment-free sequence comparison. The file is partitioned into K-Mers and a distance or correlation matrix is generated based on the counts of each K-Mer."),
				
				ui.HTML("<b>3D</b>"),
				config.HeightMatrix.UI(ui.input_select, id="HeightMatrix",	label="Height",	choices=["Distance", "Correlation", "Cube"], conditional="input.Elevation != 90", tooltip="Specify a metric to use for the height of the bars in a 3D plot. This metric can be different from the metric used for the color of the bars. The option 'Cube' displays matrices on the XY, XZ, and YZ planes, and requires an input data file with an X, Y, and Z column."),

				config.Elevation.UI(ui.input_numeric, id="Elevation", label="Elevation", tooltip="Control whether the plot is 2D or 3D. Any value other than 90 will display the plot in 3D, with the value specifying the elevation angle of the viewer in respect to the model. Change the angle back to 90 to display the plot in 2D."),
				config.Rotation.UI(ui.input_numeric, id="Rotation",	label="Rotation", conditional="input.Elevation != 90", step=1, min=1, tooltip="Change the angle of rotation of the viewer in respect to the model. For 3D plots only."),
				config.Zoom.UI(ui.input_numeric, id="Zoom",	label="Zoom", conditional="input.Elevation != 90", step=1, min=1, tooltip="Crop the view. For 3D plots only."),
				config.InterpolationLevels.UI(ui.input_numeric, id="InterpolationLevels",	label="Inter", conditional="input.Elevation != 90", step=1, min=1, max=10, tooltip="Specify a multiplier for the resolution of the 3D plot. For example, a value of 2 will interpolate the data from an NxM to a 2Nx2M, effectively quadrupling the resolution of each data point by interpolating it into 4. This results in a smoother looking plot, but can be computationally expensive for large datasets. The minimum value is 1 (default) and the maximum value is 10."),
				config.MinScale.UI(ui.input_switch, id="MinScale",	label="Scaling", conditional="input.Elevation != 90", tooltip="Scale the height of all points by the minimum value. This removes negative values and prevents data from extending below the XY plane. For 3D plots only."),
				config.Opacity.UI(ui.input_numeric, id="Opacity",	label="Opacity", conditional="input.Elevation != 90", min=0.0, max=1.0, step=0.1, tooltip="Change the opacity of the height bars in 3D plots."),

				ui.layout_columns(
					ui.HTML("<b>Colors</b>"),
					config.Custom.UI(ui.input_switch, make_inline=False, id="Custom", label="Custom", tooltip="Select colors to use in the heatmap, with the first color representing low values, and the last color representing high values. If the 'Custom' checkbox is enabled you may specify up to 12 colors. A minimum of two colors are needed."),
					col_widths=[4,8]
				),
				ui.output_ui("Color"),
				config.Bins.UI(ui.input_numeric, id="Bins", label="Number", min=3, step=1, tooltip="Specify the number of color bins to use. A higher number of color bins results in a smoother gradient between neighbouring values. Fewer bins results in more distinct colors."),

				ui.HTML("<b>Image Settings</b>"),
				config.Size.UI(ui.input_numeric, id="Size", label="Size", min=1, tooltip="Change the width (in pixels) of the heatmap on your screen."),
				config.DPI.UI(ui.input_numeric, id="DPI", label="DPI", min=5, tooltip="Specify the resolution of the image in pixels per inch. Higher DPI values result in higher quality images, but larger file sizes. This setting affects the heatmap on screen as well as the downloaded plot."),

				ui.HTML("<b>Features</b>"),
				config.Features.UI(ui.input_checkbox_group,
					make_inline=False, id="Features", label=None,
					choices={"x": "X Labels", "y": "Y Labels", "z": "Z-Labels", "label": "Data Labels", "legend": "Legend"},
					tooltip="X and Y labels toggle the data labels along their respective axes. Z labels toggles the data labels along the Z axis if rendering as a 3D plot. Data labels displays the associated value for every point on the heatmap - this can be illegible for large datasets. Legend displays a colorbar legend on the heatmap.",
				),

				ui.download_button(id="DownloadHeatmap", label="Download PNG"),
			),
			padding="10px",
			gap="20px",
			width="300px",
		),

		# Add the main interface tabs.
		MainTab(m_type=ui.output_image),
		height="86vh",
	)
)

app = App(app_ui, server)
