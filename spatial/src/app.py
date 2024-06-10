#
# Heatmapper
# Spatial
#
# This file contains the Shiny application for Spatial Heatmapper.
# It can be run with the following command within this directory:
#		shiny run
#
# Exporting via ShinyLive is not currently supported, as squidpy
# is not yet available in the Pyodide environment. Required libraries
# include: openmpi, verdict, glew, alongside python libraries in requirements.txt
# WebGL is required for this application.
#

from pandas.core.arrays.arrow.array import pa
from shiny import App, reactive, render, ui
from matplotlib.pyplot import get_cmap
from tempfile import TemporaryDirectory, NamedTemporaryFile
from anndata import read_h5ad
from squidpy import gr, pl, read
from scanpy import pp, tl
from pathlib import Path

# Shared functions
from shared import Cache, MainTab, NavBar, FileSelection, Filter, ColumnType, InitializeConfig, ColorMaps, DistanceMethods, Update, Msg, Error, Inlineify

try:
	from user import config
except ImportError:
	from config import config


def server(input, output, session):

	# Information regarding example files.
	Info = {
		"visium_hne_adata.h5ad": "Pre-processed example files provided by SquidPy",
		"seqfish.h5ad": "Pre-processed example files provided by SquidPy",
		"imc.h5ad": "Pre-processed example files provided by SquidPy",
	}

	InitializeConfig(config, input)


	def HandleData(path):
		"""
		@brief A custom Data Handler for the Cache.
		@param path: the Path to the file
		@returns A data object from the cache.
		@info This Data Handler supports h5ad files via scanpy.
		"""

		suffix = path.suffix
		if suffix == ".h5ad":
			adata = read_h5ad(path.resolve())
			gr.spatial_neighbors(adata)
			return adata

		# Pass loading .h5, we deal with them in UpdateData()
		elif suffix == ".h5": return None
		else: return DataCache.DefaultHandler(path)
	DataCache = Cache("spatial", DataHandler=HandleData)
	Data = reactive.value(None)


	def ColumnNames(adata, p):
		p.inc(message="Generating Annotation Keys...")
		Filter(adata.obs["cluster"].cat.categories.to_list(), ColumnType.Free, id="CoCluster")
		Filter(adata.obs.columns.to_list(), ColumnType.Count, id="Count")

		choices = []
		if input.UploadType() == "Visium" or input.SourceFile() == "Example":
			choices = adata.var.gene_ids.index.drop_duplicates().to_list()
		elif input.UploadType() == "NanoString":
			choices = adata.obs["fov"].drop_duplicates().to_list()
			ui.update_select(id="Count", choices=adata.obs.columns.to_list())
		if choices: ui.update_select(id="Keys", choices=choices, selected=choices[0])


	async def VisiumReader(temp, p):
		Path(f"{temp.name}/spatial").mkdir()
		counts = None
		for file in input.File():
			n = file["datapath"]
			base = file["name"]

			# These files are located in the spatial subdir
			suffix = Path(base).suffix

			# These files are considered spatial.
			if suffix in [".png", ".json", ".csv"]: base = f"spatial/{base}"
			elif suffix == ".h5": counts = base

			# If the user uploaded a .h5ad, we already have that information Cached, so just return it.
			elif suffix == ".h5ad":
				adata = await DataCache.Load(input, default=None, p=p)
				ColumnNames(adata, p)
				p.close()
				return adata

			path = f"{temp.name}/{base}"
			open(path, "wb").write(open(n, "rb").read())

		# Make SquidPy generate an object from the folder.
		p.inc(message="Consolidating Data...")
		return read.visium(temp.name, counts_file=counts)


	async def NanoStringReader(temp, p):
		counts = None
		meta = None
		fov = None

		Path(temp.name, "CellComposite").mkdir()
		Path(temp.name, "CellLabels").mkdir()

		for file in input.File():
			n = file["datapath"]
			base = file["name"]

			# These files are located in the spatial subdir
			suffix = Path(base).suffix

			# If the user uploaded a .h5ad, we already have that information Cached, so just return it.
			if suffix == ".h5ad":
				adata = await DataCache.Load(input, default=None, p=p)
				ColumnNames(adata, p)
				p.close()
				return adata
			elif suffix == ".csv":
				if base.endswith("_exprMat_file.csv"): count = base
				if base.endswith("_metadata_file.csv"): meta = base
				if "fov" in base: fov = base
			elif suffix == ".tif": base = f"CellLabels/{base}"
			elif suffix in (".png", ".jpg"): base = f"CellComposite/{base}"

			open(f"{temp.name}/{base}", "wb").write(open(n, "rb").read())

		# Make SquidPy generate an object from the folder.
		p.inc(message="Consolidating Data...")
		adata = read.nanostring(
			temp.name,
			counts_file=count,
			meta_file=meta,
			fov_file=fov
		)
		return adata

	@reactive.effect
	@reactive.event(input.SourceFile, input.File, input.Example, input.CellCount, input.GeneCount, input.UploadType)
	async def UpdateData():
		"""
		@brief Returns AnnData objects with data for Spatial Mapping.
		@info SquidPy's Visium Reader expect a directory, so Spatial will accept multiple files
			and then parse them into the correct structure.
		"""

		with ui.Progress() as p:
			p.inc(message="Loading Data...")
			if input.SourceFile() == "Upload":

				# Get all the files, to generate a name.
				if input.File() is None: return
				name = [f["datapath"] for f in input.File()]

				# If the name hasn't been cached, we need to construct the object.
				if not DataCache.In(name):
					p.inc(message="Organizing Data...")
					temp = TemporaryDirectory()
					try:
						if input.UploadType() == "Visium": adata = await VisiumReader(temp, p)
						elif input.UploadType() == "NanoString": adata = await NanoStringReader(temp, p)
						else: return None
					except Exception:
						Error("Couldn't parse the provided input! Make sure all files needed files are uploaded, and the right Upload Type is selected!")
						return None

					if adata is None: return

					# Throw it into the Cache.
					DataCache.Store(adata, name)

					# Now that it's cached, remove the origin
					for file in name:
							Path(file).unlink()

				adata = DataCache.Get(name)
				cell, gene = input.CellCount(), input.GeneCount()
				if cell is None or gene is None: return

				filtered = name + [cell, gene]
				if not DataCache.In(filtered):

					bdata = adata.copy()
					pp.filter_cells(bdata, min_counts=cell)
					pp.filter_genes(bdata , min_cells=gene)

					if input.UploadType() == "Visium":
						adata.var_names_make_unique()
						p.inc(message="Normalizing...")
						pp.normalize_total(bdata, inplace=True)
						pp.log1p(bdata)

						p.inc(message="Calculating Neighbors...")
						pp.neighbors(bdata)
						tl.umap(bdata)
						gr.spatial_neighbors(bdata)

						p.inc(message="Calculating QC Metrics...")
						pp.calculate_qc_metrics(bdata, inplace=True)

						p.inc(message="Clustering...")
						tl.leiden(bdata, key_added="cluster", neighbors_key="spatial_neighbors", resolution=input.Resolution())

						p.inc(message="Finding Highly Variable Genes...")
						pp.highly_variable_genes(bdata, inplace=True, n_top_genes=100, flavor="seurat_v3")


					elif input.UploadType() == "NanoString":
						p.inc(message="Obtaining Control Probes...")
						bdata.var["NegPrb"] = bdata.var_names.str.startswith("NegPrb")
						pp.calculate_qc_metrics(bdata, qc_vars=["NegPrb"], inplace=True)

						p.inc(message="Normalizing...")
						bdata.layers["counts"] = bdata.X.copy()
						pp.normalize_total(bdata, inplace=True)
						pp.log1p(bdata)

						p.inc(message="Calculating Neighbors...")
						pp.pca(bdata)
						pp.neighbors(bdata)
						tl.umap(bdata)
						gr.spatial_neighbors(bdata, coord_type="generic", delaunay=True)

						p.inc(message="Clustering...")
						tl.leiden(bdata, key_added="cluster")

					ColumnNames(bdata ,p)
					DataCache.Store(bdata, filtered)
				Data.set(DataCache.Get(filtered))
				p.close()

			# With an example, just return it.
			else:
				Data.set(await DataCache.Load(input, default=None))
				ColumnNames(Data(), p)
				p.close()


	@output
	@render.data_frame
	def Table():
		state = config.TableType()
		df = Data()
		if df is None: return
		if state == "obs": return render.DataGrid(df.obs)
		elif state == "var": return render.DataGrid(df.var)


	def GenerateNanoString(adata, file, p):
		id = config.Keys()
		count = config.Count()
		if adata is None or id is None or count is None: return

		shape = config.Shape().lower()
		features = config.Features()
		img_alpha = config.ImgOpacity()
		cmap = config.ColorMap().lower()
		alpha = config.Opacity()
		columns = config.Columns()
		spacing = config.Spacing()
		dpi = config.DPI()

		p.inc(message="Plotting...")
		pl.spatial_segment(
			adata,
			color=count,
			library_key="fov",
			seg_cell_id="cell_ID",
			library_id=id,
			shape=shape,
			img="Image" in features,
			img_alpha=img_alpha,
			cmap=get_cmap(cmap),
			alpha=alpha,
			colorbar="Legend" in features,
			frameon="Frame" in features,
			ncols=columns,
			wspace=spacing,
			hspace=spacing,
			save=file.name,
			dpi=dpi
		)
		img: types.ImgData = {"src": file.name, "height": f"{config.Size()}vh"}
		return img

	def GenerateVisium(adata, file, p):
		colors = config.Keys()
		if adata is None or colors is None: return
		try:
			genes = adata[:, adata.var.highly_variable].var_names.values[:100]
		except AttributeError:
			genes=None

		p.inc(message="Computing statistic...")
		stat = config.Statistic()
		if stat == "sepal" and "sepal_score" not in adata.uns:
			gr.sepal(
				adata,
				genes=genes,
				max_neighs=6,
				show_progress_bar=False,
			)
		elif (stat == "moran" and "moranI" not in adata.uns) or (stat == "geary" and "gearyC" not in adata.uns):
			gr.spatial_autocorr(
				adata,
				genes=genes,
				mode=stat,
			)

		p.inc(message="Plotting...")
		shape = config.Shape().lower()
		features = config.Features()
		img_alpha = config.ImgOpacity()
		cmap = config.ColorMap().lower()
		alpha = config.Opacity()
		columns = config.Columns()
		spacing = config.Spacing()
		dpi = config.DPI()

		pl.spatial_scatter(
			adata,
			color=colors,
			shape=shape,
			img="Image" in features,
			img_alpha=img_alpha,
			cmap=get_cmap(cmap),
			alpha=alpha,
			colorbar=len(colors) > 1 and "Legend" in features,
			frameon="Frame" in features,
			ncols=columns,
			wspace=spacing,
			hspace=spacing,
			save=file.name,
			dpi=dpi
		)
		img: types.ImgData = {"src": file.name, "height": f"{config.Size()}vh"}
		return img


	def GenerateHeatmap(file=None):
		"""
		@brief Generates the Annotation Key Spatial Scatter Heatmap, returning an ImgData dictionary
		@param file: An optional file to use.
		@returns: An ImgData dictionary

		@info By default, this function will generate a non-deleting TemporaryFile, which means the calling function has the
		responsibility of deleting it (Which Shiny handles via render.image). However, for applications like DownloadHeatmap,
		it has no way to tell Shiny to delete the file, so rather than manually dealing with cleanup, we can just pass it a
		NamedTemporaryFile within the scope if it's function, and let the Operating System delete it after the call.
		"""
		with ui.Progress() as p:

			p.inc(message="Loading input...")
			adata = Data()
			if file is None: file = NamedTemporaryFile(delete=False, suffix=".png")
			if input.SourceFile() == "Example" or input.UploadType() == "Visium":
				return GenerateVisium(adata, file, p)
			elif input.UploadType() == "NanoString":
				return GenerateNanoString(adata, file, p)


	@output
	@render.image(delete_file=True)
	def Heatmap(): return GenerateHeatmap()


	@output
	@render.image(delete_file=True)
	def HeatmapReactive(): return GenerateHeatmap()


	@output
	@render.plot
	def Centrality():
		with ui.Progress() as p:

			p.inc(message="Loading input...")
			adata = Data()
			score = config.Score()

			if adata is None: return

			key = "cluster"
			location = f"{key}_centrality_scores"

			p.inc(message="Computing score...")
			if location not in adata.uns or score not in adata.uns[location]:
				gr.centrality_scores(
					adata,
					cluster_key=key,
					show_progress_bar=False,
					score=score,
				)

			p.inc(message="Plotting...")
			pl.centrality_scores(adata, key, score=score)


	@output
	@render.plot
	def Ripley():
		with ui.Progress() as p:

			p.inc(message="Loading input...")
			adata = Data()
			if adata is None: return

			function = config.Function()
			metric = config.Distance().lower()

			# Because the metric is not uniquely identified within the adata, we cache it
			# and check if the user has changed it. If it has changed, we need to recompute.
			# However, we don't Cache the actual calculation, just the metric, as we would
			# be caching the information twice.
			hash_list = [config.Function(), input.SourceFile(), input.Example(), input.File()]
			old_metric = DataCache.Get(hash_list)

			key = "cluster"

			p.inc(message="Generating function...")
			if f"{key}_ripley_{function}" not in adata.uns or metric != old_metric:
				gr.ripley(
					adata,
					cluster_key=key,
					mode=function,
					metric=metric,
				)
				DataCache.Store(metric, hash_list)

			p.inc(message="Plotting...")
			pl.ripley(adata, cluster_key=key, mode=function)


	@output
	@render.plot
	def Occurrence():
		with ui.Progress() as p:

			p.inc(message="Loading input...")
			adata = Data()
			if adata is None: return

			if input.UploadType() == "NanoString":
				adata = adata[adata.obs.fov.isin(input.Keys())].copy()

			key = "cluster"
			p.inc(message="Calculating...")

			interval = config.Interval()
			splits = None if input.Splits() == 0 else input.Splits()

			if f"{key}_co_occurrence" not in adata.uns:
				gr.co_occurrence(
					adata,
					cluster_key=key,
					interval=interval,
					n_splits=splits,
					show_progress_bar=False,
				)

			p.inc(message="Plotting...")
			if config.OccurrenceGraph() == "Line" and config.CoCluster() is not None:
				pl.co_occurrence(adata, cluster_key=key, clusters=config.CoCluster())
			else:
				pl.spatial_scatter(adata, color=key, size=10, shape=None)


	@reactive.effect
	@reactive.event(input.ExampleInfo)
	def ExampleInfo():
		Msg(ui.HTML(Info[input.Example()]))


	@render.download(filename="adata.h5ad")
	def DownloadTable():
		adata = Data()
		if adata is None: return
		temp = NamedTemporaryFile();
		adata.write(temp.name)
		yield open(temp.name, "rb").read()


	@render.download(filename="heatmap.png")
	def DownloadHeatmap():
		# Generate a TemporaryFile that deletes itself after scope.
		with NamedTemporaryFile(suffix=".png") as file:
			data = GenerateHeatmap(file)
			if data is None: return
			yield file.read()


app_ui = ui.page_fluid(

	NavBar(),

	ui.layout_sidebar(
		ui.sidebar(
			FileSelection(
				examples={
					"visium_hne_adata.h5ad": "Example 1",
					"seqfish.h5ad": "Example 2",
					"imc.h5ad": "Example 3",
				},
				types=[".h5", ".png", ".csv", ".json", ".h5ad", ".jpg", ".tif"],
				multiple=True,
				default="Upload",
				project="Spatial"
			),

			ui.input_select(id="UploadType", label=None, choices=["Visium", "VizGen", "NanoString"]),

			ui.panel_conditional(
				"input.MainTab === 'TableTab'",
				config.TableType.UI(ui.input_select, id="TableType", label="Table", choices={"obs": "Observations", "var": "Variable"})
			),

			ui.panel_conditional("input.MainTab != 'TableTab'",
				Update(),

				ui.HTML("<b>Minium Count Filtering</b>"),
				Inlineify(ui.input_numeric, id="GeneCount", label="Gene", min=0, value=400),
				Inlineify(ui.input_numeric, id="CellCount", label="Cell", min=0, value=100),

				ui.HTML("<b>Keys</b>"),
				config.Keys.UI(ui.input_select, id="Keys", label="Keys", choices=[], selectize=True, multiple=True),
				config.Count.UI(ui.input_select, id="Count", label="Count", choices=[]),
			),

			ui.panel_conditional(
				"input.MainTab === 'HeatmapTab'",
				ui.HTML("<b>Heatmap</b>"),
				config.Statistic.UI(ui.input_select, id="Statistic", label="Statistic", choices={"moran": "Moran's I", "sepal": "Sepal", "geary": "Geary's C"}),
				config.ColorMap.UI(ui.input_select, id="ColorMap", label="Map", choices=ColorMaps),
				config.Shape.UI(ui.input_select, id="Shape", label="Shape", choices=["Circle", "Square", "Hex"]),
				config.Columns.UI(ui.input_slider, id="Columns", label="Columns", min=1, max=10, step=1),
				config.Spacing.UI(ui.input_slider, id="Spacing", label="Spacing", min=0.0, max=1.0, step=0.1),

				ui.HTML("<b>Opacity</b>"),
				config.ImgOpacity.UI(ui.input_slider, id="ImgOpacity", label="Image", min=0.0, max=1.0, step=0.1),
				config.Opacity.UI(ui.input_slider, id="Opacity", label="Data", min=0.0, max=1.0, step=0.1),

				ui.HTML("<b>Image Settings</b>"),
				config.Size.UI(ui.input_numeric, id="Size", label="Size", min=1),
				config.DPI.UI(ui.input_numeric, id="DPI", label="DPI", min=1),

				ui.HTML("<b>Features</b>"),
				config.Features.UI(ui.input_checkbox_group, make_inline=False, id="Features", label=None, 	choices=["Image", "Legend", "Frame"]),

				ui.download_button(id="DownloadHeatmap", label="Download"),
			),

			ui.panel_conditional(
				"input.MainTab === 'Centrality'",
				ui.HTML("<b>Centrality</b>"),
				config.Score.UI(ui.input_select,
					id="Score",
					label="Score",
					choices={
						"closeness_centrality": "Closeness Centrality",
						"average_clustering": "Average Clustering",
						"degree_centrality": "Degree Centrality"
					},
				),
			),

			ui.panel_conditional(
				"input.MainTab === 'Ripley'",
				ui.HTML("<b>Ripley</b>"),
				config.Function.UI(ui.input_select, id="Function", label="Function", choices=["L", "F", "G"]),
				config.Distance.UI(ui.input_select, id="Distance", label="Distance", choices=DistanceMethods),
			),


			ui.panel_conditional(
				"input.MainTab === 'Occurrence'",
				ui.HTML("<b>Co-Occurrence</b>"),
				config.CoCluster.UI(ui.input_select, id="CoCluster", label="Group", choices=[]),
				config.OccurrenceGraph.UI(ui.input_select, id="OccurrenceGraph", label="Graph", choices=["Scatter", "Line"]),
				config.Interval.UI(ui.input_slider, id="Interval", label="Interval", min=1, max=100, step=1),
				config.Splits.UI(ui.input_slider, id="Splits", label="Splits", min=0, max=10, step=0),
			),


			ui.panel_conditional(
				"input.MainTab === 'TableTab'",

				# Add the download buttons.
				ui.download_button("DownloadTable", "Download Table"),
			),
			padding="10px",
			gap="20px",
			width="250px",
		),

		# Add the main interface tabs.
		MainTab(
			ui.nav_panel("Centrality Scores", ui.output_plot("Centrality", height="90vh"), value="Centrality"),
			ui.nav_panel("Ripley's Function", ui.output_plot("Ripley", height="90vh"), value="Ripley"),
			ui.nav_panel("Co-occurrence", ui.output_plot("Occurrence", height="90vh"), value="Occurrence")
		),
		height="90vh",
	)
)

app = App(app_ui, server)
