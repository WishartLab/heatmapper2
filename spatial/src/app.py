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
from io import BytesIO
from matplotlib.pyplot import get_cmap
from pathlib import Path
from tempfile import NamedTemporaryFile, TemporaryDirectory
from pandas import DataFrame
from numpy import append
from anndata import read_h5ad
from squidpy import gr, pl, read
from shiny.types import FileInfo
from anndata import AnnData
from scanpy import pp, tl

# scikit-misc needed for highly_variable_genes

# Shared functions
from shared import Cache, MainTab, NavBar, FileSelection, Filter, ColumnType, TableValueUpdate


def server(input, output, session):

	# Information regarding example files.
	Info = {
		"example1.h5ad": "Adult Mouse Brain Section 2, from https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Adult_Mouse_Brain_Coronal_Section_2?",
	}

	# The DataCache holds our H5AD files, but DataCache serves as a Cache to conserve
	# bandwidth, not computation. The SpatialCache caches objects computed via
	# squidpy.gr functions, using a hash of the input parameters
	SpatialCache = {}

	# Concurrent jobs to run for calculations.
	Jobs = 8


	def HandleData(n, i):
		"""
		@brief A custom Data Handler for the Cache.
		@param n: The name of the file
		@param i: The source of the file. It can be a path to a file (string) or a BytesIO object.
		@returns A data object from the cache.
		@info This Data Handler supports h5ad files via scanpy.
		"""

		suffix = Path(n).suffix
		if suffix == ".h5ad":
			if type(i) is BytesIO:
				temp = NamedTemporaryFile(suffix=suffix); temp.write(i.read())
				i = temp.name
			adata = read_h5ad(i)
			gr.spatial_neighbors(adata)
			return adata

		# Pass loading .h5, we deal with them in Load()
		elif suffix == ".h5": return None
		else: return DataCache.DefaultHandler(n, i)
	DataCache = Cache("spatial", DataHandler=HandleData)


	def Load():
		"""
		@brief Returns AnnData objects with data for Spatial Mapping.
		"""

		if input.SourceFile() == "Upload":

			# Get all the files, to generate a name.
			if input.File() is None: return None
			name = ""
			for file in input.File():
				name += file["datapath"]

			# If the name hasn't been cached, we need to construct the object.
			if name not in SpatialCache:

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
					elif suffix == ".h5ad": return DataCache.SyncLoad(input, default=None)

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
				SpatialCache[name] = adata

			# Return the cached information.
			return SpatialCache[name]

		# With an example, just return it.
		else: return DataCache.SyncLoad(input, default=None)


	@output
	@render.data_frame
	async def LoadedTable(): return Load().to_df()


	@output
	@render.plot
	def Heatmap():
		adata = Load()
		if adata is None: return
		genes = adata[:, adata.var.highly_variable].var_names.values[:100]


		stat = input.Statistic()
		if stat == "sepal":
			gr.sepal(
				adata,
				genes=genes,
				max_neighs=6,
				n_jobs=Jobs,
				show_progress_bar=False,

			)
		else:
			gr.spatial_autocorr(
				adata,
				genes=genes,
				mode=stat,
				n_perms=input.Perm(),
				two_tailed=input.P() == "Two-Tailed",
				corr_method=input.Correlation()
			)

		pl.spatial_scatter(
			adata,
			color=input.Keys(),
			shape=input.Shape().lower(),
			img="Image" in input.Features(),
			img_alpha=input.ImgOpacity(),
			cmap=get_cmap(input.ColorMap().lower()),
			alpha=input.Opacity(),
			colorbar=False
		)


	@output
	@render.plot
	def Centrality():
		adata = Load()
		if adata is None: return

		gr.centrality_scores(adata, input.Key())
		pl.centrality_scores(adata, input.Key())


	@output
	@render.plot
	def Ripley():
		adata = Load()
		if adata is None: return
		gr.ripley(adata, cluster_key=input.Key(), mode=input.Function())
		pl.ripley(adata, cluster_key=input.Key(), mode=input.Function())


	@output
	@render.plot
	def Occurrence():
		adata = Load()
		if adata is None: return

		gr.co_occurrence(
			adata,
			cluster_key=input.Key(),
			interval=input.Interval(),
			n_splits=None if input.Splits() == 0 else input.Splits(),
			show_progress_bar=False
		)

		if input.OccurrenceGraph() == "Line" and input.Cluster() is not None:
			pl.co_occurrence(adata, cluster_key=input.Key(), clusters=input.Cluster())
		else:
			pl.spatial_scatter(adata, color=input.Key(), size=10, shape=None)


	@output
	@render.text
	@reactive.event(input.SourceFile, input.Example)
	def ExampleInfo(): return Info[input.Example()]


	@render.download(filename="table.csv")
	def DownloadTable():
		df = Load().to_df();
		if df is not None:
			yield df.to_string()


	@reactive.Effect
	async def UpdateColumnSelection():
		adata = Load()
		if adata is None: return
		key = Filter(adata.obs.columns, ColumnType.Cluster, only_one=True, ui_element="Key")
		if key is not None:
			Filter(adata.obs[key].cat.categories.tolist(), ColumnType.Free, ui_element="Cluster")
		if not input.Keys():
			ui.update_select(id="Keys", label="Annotation Keys", choices=adata.var.gene_ids.index.drop_duplicates().to_list())


app_ui = ui.page_fluid(

	NavBar("Spatial"),

	ui.layout_sidebar(
		ui.sidebar(
			FileSelection(
				upload_label="Upload Files",
				examples={
					"seqfish.h5ad": "Example 1",
					"imc.h5ad": "Example 2",
				},
				types=[".h5", ".png", ".csv", ".json", ".h5ad"],
				multiple=True,
				default="Upload"
			),

			ui.input_select(id="Key", label="Key", choices=[], selected=None),
			ui.input_select(id="Cluster", label="Cluster", choices=[], selected=None),


			ui.panel_conditional(
				"input.MainTab === 'Interactive'",
				ui.HTML("<b>Interactive Settings</b>"),
				ui.input_select(id="Statistic", label="Statistic", choices={"moran": "Moran's I", "sepal": "Sepal", "geary": "Geary's C"}),

				ui.input_slider(id="Perm", label="Permutations", value=100, min=1, max=200, step=1),

				ui.input_radio_buttons(id="P", label="P-Value Distribution", choices=["Two-Tailed", "One-Tailed"], selected="One-Tailed", inline=True),

				# https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html#statsmodels.stats.multitest.multipletests
				ui.input_select(id="Correlation", label="Correlation Method", choices=["bonferroni" "sidak", "holm-sidak", "holm", "simes-hochberg", "hommel", "fdr_bh", "fdr_by", "fdr_tsbh", "fdr_tsbky"], selected="holm-sidak"),

				ui.input_select(id="Keys", label="Annotation Keys", choices=[], selected=None),

				ui.HTML("<b>Visual Settings</b>"),

				ui.input_select(id="ColorMap", label="Color Map", choices=["Viridis", "Plasma", "Inferno", "Magma", "Cividis"], selected="Viridis"),
				ui.input_select(id="Shape", label="Shape", choices=["Circle", "Square", "Hex"], selected="Hex"),

				ui.input_slider(id="ImgOpacity", label="Image Opacity", value=1, min=0.0, max=1.0, step=0.1),
				ui.input_slider(id="Opacity", label="Data Opacity", value=1, min=0.0, max=1.0, step=0.1),

				ui.input_checkbox_group(id="Features", label="Heatmap Features", choices=["Image"], selected=["Image"]),
			),

			ui.panel_conditional(
				"input.MainTab === 'Neighbors'",
				ui.HTML("<b>Spatial Neighbors Settings</b>"),
			),


			ui.panel_conditional(
				"input.MainTab === 'Ripley'",
				ui.HTML("<b>Ripley Settings</b>"),
				ui.input_select(id="Function", label="Function", choices=["L", "F", "G"]),
			),

			ui.panel_conditional(
				"input.MainTab === 'Occurrence'",
				ui.HTML("<b>Co-occurrence Settings</b>"),
				ui.input_radio_buttons(id="OccurrenceGraph", label="Graph Type", choices=["Scatter", "Line"], inline=True),

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
