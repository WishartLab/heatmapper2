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

from shiny import App, reactive, render, ui
from matplotlib.pyplot import get_cmap
from tempfile import TemporaryDirectory, NamedTemporaryFile
from anndata import read_h5ad
from squidpy import gr, pl, read
from scanpy import pp, tl

# Shared functions
from shared import Cache, MainTab, NavBar, FileSelection, Filter, ColumnType


def server(input, output, session):

	# Information regarding example files.
	Info = {
		"visium_hne_adata.h5ad": "Pre-processed example files provided by SquidPy",
		"seqfish.h5ad": "Pre-processed example files provided by SquidPy",
		"imc.h5ad": "Pre-processed example files provided by SquidPy",
	}

	# Concurrent jobs to run for calculations.
	Jobs = 8


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


	@reactive.effect
	@reactive.event(input.SourceFile, input.File, input.Example)
	async def UpdateData():
		"""
		@brief Returns AnnData objects with data for Spatial Mapping.
		@info SquidPy's Visium Reader expect a directory, so Spatial will accept multiple files
			and then parse them into the correct structure.
		"""

		if input.SourceFile() == "Upload":

			# Get all the files, to generate a name.
			if input.File() is None: return
			name = [f["datapath"] for f in input.File()]
			
	
			# If the name hasn't been cached, we need to construct the object.
			if not DataCache.In(name):

				# For each file uploaded, dump it into a temporary directory
				temp = TemporaryDirectory()
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
					elif suffix == ".h5ad": return await DataCache.Load(input, default=None)

					open(f"{temp.name}/{base}", "wb").write(open(n, "rb").read())

				# Make SquidPy generate an object from the folder.
				adata = read.visium(temp.name, counts_file=counts)

				# Post-processing so all the needed statistics are present.
				adata.var_names_make_unique()
				pp.filter_genes(adata, inplace=True, min_counts=100)
				gr.spatial_neighbors(adata)
				pp.calculate_qc_metrics(adata, inplace=True)
				tl.leiden(adata, key_added="cluster", neighbors_key="spatial_neighbors")
				pp.highly_variable_genes(adata, inplace=True, n_top_genes=100, flavor="seurat_v3")

				# Throw it into the Cache.
				DataCache.Store(adata, name)

			# Return the cached information.
			Data.set(DataCache.Get(name))

		# With an example, just return it.
		else: Data.set(await DataCache.Load(input, default=None))


	@output
	@render.data_frame
	def Table(): 
		if input.TableType() == "obs":
			return render.DataGrid(Data().obs)
		elif input.TableType() == "var":
			return render.DataGrid(Data().var)


	@output
	@render.plot
	def Heatmap():

		with ui.Progress() as p:

			p.inc(message="Loading input...")
			adata = Data()
			colors = input.Keys()

			if adata is None or colors is None: return

			try:
				genes = adata[:, adata.var.highly_variable].var_names.values[:100]
			except AttributeError:
				genes=None


			p.inc(message="Computing statistic...")
			stat = input.Statistic()
			if stat == "sepal" and "sepal_score" not in adata.uns:
				gr.sepal(
					adata,
					genes=genes,
					max_neighs=6,
					n_jobs=Jobs,
					show_progress_bar=False,

				)
			elif (stat == "moran" and "moranI" not in adata.uns) or (stat == "geary" and "gearyC" not in adata.uns):
				gr.spatial_autocorr(
					adata,
					genes=genes,
					mode=stat,
				)

			p.inc(message="Plotting...")

			# Shiny does not seem to properly setup dependencies to reactive values if
			# They are used directly in function invocations. Therefore, we need
			# To create variables using these reactives, that way changes properly 
			# Recall this function.
			shape = input.Shape().lower()
			features = input.Features()
			img_alpha = input.ImgOpacity()
			cmap = input.ColorMap().lower()
			alpha = input.Opacity()
			columns = input.Columns()
			spacing = input.Spacing()

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
			)


	@output
	@render.plot
	def Centrality():
		with ui.Progress() as p:

			p.inc(message="Loading input...")
			adata = Data()
			score = input.Score()

			if adata is None: return

			key = Filter(adata.obs.columns, ColumnType.Cluster, only_one=True)
			location = f"{key}_centrality_scores"

			p.inc(message="Computing score...")
			if location not in adata.uns or score not in adata.uns[location]:
				gr.centrality_scores(
					adata,
					cluster_key=key,
					n_jobs=Jobs,
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

			function = input.Function()
			metric = input.Distance().lower()

			# Because the metric is not uniquely identified within the adata, we cache it
			# and check if the user has changed it. If it has changed, we need to recompute.
			# However, we don't Cache the actual calculation, just the metric, as we would
			# be caching the information twice.
			hash_list = [input.Function(), input.SourceFile(), input.Example(), input.File()]
			old_metric = DataCache.Get(hash_list)

			key = Filter(adata.obs.columns, ColumnType.Cluster, only_one=True)

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

			key = Filter(adata.obs.columns, ColumnType.Cluster, only_one=True)

			p.inc(message="Calculating...")

			if f"{key}_co_occurrence" not in adata.uns:
				gr.co_occurrence(
					adata,
					cluster_key=key,
					interval=input.Interval(),
					n_splits=None if input.Splits() == 0 else input.Splits(),
					show_progress_bar=False,
					n_jobs=Jobs,
				)

			p.inc(message="Plotting...")
			if input.OccurrenceGraph() == "Line" and input.Cluster() is not None:
				pl.co_occurrence(adata, cluster_key=key, clusters=input.Cluster())
			else:
				pl.spatial_scatter(adata, color=key, size=10, shape=None)


	@output
	@render.text
	@reactive.event(input.SourceFile, input.Example)
	def ExampleInfo(): return Info[input.Example()]


	@render.download(filename="adata.h5ad")
	def DownloadTable():
		adata = Data()
		temp = NamedTemporaryFile(); 
		adata.write(temp.name)
		yield open(temp.name, "rb").read()


	@reactive.Effect
	def UpdateColumnSelection():
		with ui.Progress() as p:
			p.inc(message="Loading data...")
			adata = Data()

			p.inc(message="Generating Annotation Keys...")
			if adata is None: return
			if input.MainTab() == "Occurrence":
				key = Filter(adata.obs.columns, ColumnType.Cluster, only_one=True)
				if key is not None:
					Filter(adata.obs[key].cat.categories.tolist(), ColumnType.Free, ui_element="Cluster")
			if not input.Keys():
				try:
					choices = adata.var.gene_ids.index.drop_duplicates().to_list()
					ui.update_select(id="Keys", label="Annotation Keys", choices=choices, selected=choices[0])
				except AttributeError:
					pass


app_ui = ui.page_fluid(

	NavBar(),

	ui.layout_sidebar(
		ui.sidebar(
			FileSelection(
				upload_label="Upload Files",
				examples={
					"visium_hne_adata.h5ad": "Example 1",
					"seqfish.h5ad": "Example 2",
					"imc.h5ad": "Example 3",
				},
				types=[".h5", ".png", ".csv", ".json", ".h5ad"],
				multiple=True,
				default="Upload",
				project="Spatial"
			),

			ui.panel_conditional(
				"input.MainTab === 'TableTab'",
				ui.input_radio_buttons(id="TableType", label="Table", choices={"obs": "Observations", "var": "Variable"})
			),

			ui.panel_conditional(
				"input.MainTab === 'HeatmapTab'",
				ui.HTML("<b>Heatmap Settings</b>"),
				ui.input_select(id="Statistic", label="Statistic", choices={"moran": "Moran's I", "sepal": "Sepal", "geary": "Geary's C"}),
				ui.input_select(id="Keys", label="Annotation Keys", choices=[], selected=None, selectize=True, multiple=True),

				ui.HTML("<b>Visual Settings</b>"),
				ui.input_select(id="ColorMap", label="Color Map", choices=["Viridis", "Plasma", "Inferno", "Magma", "Cividis"], selected="Viridis"),
				ui.input_select(id="Shape", label="Shape", choices=["Circle", "Square", "Hex"], selected="Hex"),

				ui.input_slider(id="ImgOpacity", label="Image Opacity", value=1, min=0.0, max=1.0, step=0.1),
				ui.input_slider(id="Opacity", label="Data Opacity", value=1, min=0.0, max=1.0, step=0.1),

				ui.input_slider(id="Columns", label="Columns", value=2, min=1, max=10, step=1),
				ui.input_slider(id="Spacing", label="Spacing", value=0.3, min=0.0, max=1.0, step=0.1),

				ui.input_checkbox_group(id="Features", label="Heatmap Features", choices=
					["Image", "Legend", "Frame"], 
					selected=["Image", "Legend"]
				),
			),

			ui.panel_conditional(
				"input.MainTab === 'Centrality'",
				ui.HTML("<b>Centrality Settings</b>"),
				ui.input_select(
				id="Score",
				label="Score",
				choices={
					"closeness_centrality": "Closeness Centrality", 
					"average_clustering": "Average Clustering",
					"degree_centrality": "Degree Centrality"
				},
				selected=["Closeness Centrality"],
				),
			),

			ui.panel_conditional(
				"input.MainTab === 'Ripley'",
				ui.HTML("<b>Ripley Settings</b>"),
				ui.input_select(id="Function", label="Function", choices=["L", "F", "G"]),
				ui.input_select(id="Distance", label="Distance Method", choices=["Euclidean", "Manhattan", "Chebyshev", "Minkowski"]),
			),

			ui.panel_conditional(
				"input.MainTab === 'Occurrence'",
				ui.HTML("<b>Co-occurrence Settings</b>"),
				ui.input_radio_buttons(id="OccurrenceGraph", label="Graph Type", choices=["Scatter", "Line"], inline=True),

				ui.panel_conditional(
					"input.OccurrenceGraph === 'Line'",
					ui.input_select(id="Cluster", label="Cluster", choices=[], selected=None),
				),

				ui.input_slider(id="Interval", label="Distance Interval", value=50, min=1, max=100, step=1),
				ui.input_slider(id="Splits", label="Splits (0 = auto)", value=0, min=0, max=10, step=0),
			),


			# Add the download buttons.
			ui.download_button("DownloadTable", "Download Table"),
		),

		# Add the main interface tabs.
		MainTab(
			ui.nav_panel("Centrality Scores", ui.output_plot("Centrality", height="90vh"), value="Centrality"),
			ui.nav_panel("Ripley's Function", ui.output_plot("Ripley", height="90vh"), value="Ripley"),
			ui.nav_panel("Co-occurrence", ui.output_plot("Occurrence", height="90vh"), value="Occurrence")
		),
	)
)

app = App(app_ui, server)
