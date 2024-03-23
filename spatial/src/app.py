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
from tempfile import NamedTemporaryFile
from pandas import DataFrame
from numpy import append
from anndata import read_h5ad
from squidpy import gr, pl

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


	def HandleData(n, i):
		"""
		@brief A custom Data Handler for the Cache.
		@param n: The name of the file
		@param i: The source of the file. It can be a path to a file (string) or a BytesIO object.
		@returns A data object from the cache.
		@info This Data Handler supports h5ad files via scanpy.
		"""

		suffix = Path(n).suffix
		match suffix:
			case ".h5ad" | ".h5":
				if type(i) is BytesIO:
					temp = NamedTemporaryFile(suffix=suffix); temp.write(i.read())
					i = temp.name
				return read_h5ad(i)
			case _: return DataCache.DefaultHandler(n, i)
	DataCache = Cache("spatial", DataHandler=HandleData)


	def GenerateSpatial(n, adata):
		h = hash(f"Spatial{n}{input.Coordinate()}{input.Neighs()}{input.Rings()}{input.Radius()}{input.Percentile()}{input.SpatialFeatures()}{input.Transform()}")

		if h not in SpatialCache:
			copy = adata.copy()
			gr.spatial_neighbors(
				copy,
				coord_type=input.Coordinate().lower(),
				n_neighs=input.Neighs(),
				n_rings=input.Rings(),
				radius=None if input.Radius() == 0 else input.Radius(),
				percentile=None if input.Percentile() == 0.0 else input.Percentile(),
				delaunay="Delaunay" in input.SpatialFeatures(),
				set_diag="Diagonal = 1.0" in input.SpatialFeatures(),
				transform=None if input.Transform() == "None" else input.Transform(),
			)
			SpatialCache[h] = copy
		return SpatialCache[h]


	def GenerateOccurrence(n, adata):
		h = hash(f"Occurrence{n}{input.Key()}{input.Interval()}{input.Clusters()}{input.Splits()}")
		if h not in SpatialCache:
			copy = adata.copy()
			gr.co_occurrence(
				adata,
				cluster_key=input.Key(),
				interval=input.Interval(),
				n_splits=None if input.Splits() == 0 else input.Splits(),
				show_progress_bar=False
			)
			SpatialCache[h] = copy
		return SpatialCache[h]


	@output
	@render.data_frame
	@reactive.event(input.SourceFile, input.File, input.Example, input.File, input.Update, input.Reset)
	async def LoadedTable(): return (await DataCache.Load(input)).to_df()


	@output
	@render.plot
	@reactive.event(input.SourceFile, input.File, input.Example, input.Update, input.Reset, input.Key, input.Neighs, input.Radius, input.Rings, input.Percentile, input.Transform, input.SpatialFeatures, input.Statistic, input.Shape, input.Features, input.Opacity, input.ColorMap, input.Keys)
	def Heatmap():
		n, adata = DataCache.SyncLoad(input, default=None, return_n=True)
		if adata is None: return

		spatial = GenerateSpatial(n, adata)

		match input.Statistic():
			case "Moran's I":
				genes = adata[:, adata.var.highly_variable].var_names.values[:100]
				gr.spatial_autocorr(
					spatial,
					mode="moran",
					genes=genes,
					n_perms=100,
				)
			case "Sepal":
				genes = adata.var_names[(adata.var.n_cells > 100) & adata.var.highly_variable][0:100]
				gr.sepal(
					spatial,
					max_neighs=6,
					genes=genes,
				)

		pl.spatial_scatter(
			adata,
			color=input.Keys(),
			shape=input.Shape().lower(),
			img="Image" in input.Features(),
			img_alpha=input.ImgOpacity(),
			cmap=get_cmap(input.ColorMap().lower()),
			alpha=input.Opacity(),
			colorbar="Legend" in input.Features() and len(input.Keys()) > 1,
		)


	@output
	@render.plot
	@reactive.event(input.SourceFile, input.File, input.Example, input.Update, input.Reset, input.Key)
	def Centrality():
		n, adata = DataCache.SyncLoad(input, default=None, return_n=True)
		if adata is None: return

		spatial = GenerateSpatial(n, adata)

		gr.centrality_scores(spatial, input.Key())
		pl.centrality_scores(spatial, input.Key())


	@output
	@render.plot
	@reactive.event(input.SourceFile, input.File, input.Example, input.Update, input.Reset, input.Key)
	def Neighbors():
		n, adata = DataCache.SyncLoad(input, default=None, return_n=True)
		if adata is None: return

		spatial = GenerateSpatial(n, adata)
		_, idx = spatial.obsp["spatial_connectivities"][420, :].nonzero()
		idx = append(idx, 420)

		pl.spatial_scatter(
				spatial[idx, :],
				shape=None,
				color=input.Key(),
				connectivity_key="spatial_connectivities",
				size=100,
		)


	@output
	@render.plot
	@reactive.event(input.SourceFile, input.File, input.Example, input.Update, input.Reset, input.Function)
	def Ripley():
		adata = DataCache.SyncLoad(input, default=None)
		if adata is None: return

		gr.ripley(adata, cluster_key=input.Key(), mode=input.Function())
		pl.ripley(adata, cluster_key=input.Key(), mode=input.Function())


	@output
	@render.plot
	@reactive.event(input.SourceFile, input.File, input.Example, input.Update, input.Reset)
	def ReceptorLigand():
		adata = DataCache.SyncLoad(input, default=None)
		if adata is None: return

		res = gr.ligrec(
			adata,
			n_perms=1000,
			cluster_key=input.Key(),
			copy=True,
			use_raw=False,
			transmitter_params={"categories": "ligand"},
			receiver_params={"categories": "receptor"},
		)
		pl.ligrec(res, source_groups="Erythroid", alpha=0.005)


	@output
	@render.plot
	@reactive.event(input.SourceFile, input.File, input.Example, input.Update, input.Reset, input.OccurrenceGraph, input.Key, input.Clusters, input.Coordinate, input.Interval, input.Neighs, input.Rings, input.Radius, input.SpatialFeatures, input.Percentile, input.Transform, input.Splits)
	def Occurrence():
		n, adata = DataCache.SyncLoad(input, default=None, return_n=True)
		if adata is None: return

		spatial = GenerateSpatial(n, adata)
		occurrence = GenerateOccurrence(n, spatial)

		if input.OccurrenceGraph() == "Line" and input.Clusters() is not None:
			pl.co_occurrence(
				occurrence,
				cluster_key=input.Key(),
				clusters=input.Clusters(),

			)
		else:
			pl.spatial_scatter(occurrence, color=input.Key(), size=10, shape=None)


	@output
	@render.text
	@reactive.event(input.SourceFile, input.Example)
	def ExampleInfo(): return Info[input.Example()]


	@render.download(filename="table.csv")
	async def DownloadTable():
		df = (await DataCache.Load(input)).to_df();
		if df is not None:
			yield df.to_string()


	@reactive.Effect
	@reactive.event(input.Update)
	async def Update(): await DataCache.Update(input)


	@reactive.Effect
	@reactive.event(input.Reset)
	async def Reset(): await DataCache.Purge(input)


	@reactive.Effect
	@reactive.event(input.SourceFile, input.File, input.Example, input.Update, input.Reset, input.MainTab, input.OccurrenceGraph)
	async def UpdateColumnSelection():
		adata = DataCache.SyncLoad(input, default=None)
		if adata is None: return

		Filter(adata.obs.columns, ColumnType.Cluster, ui_element="Key")

		if input.MainTab() == "Occurrence" and input.OccurrenceGraph() == "Line":
			Filter(adata.obs[input.Key()].cat.categories.tolist(), ColumnType.Free, ui_element="Clusters")


app_ui = ui.page_fluid(

	NavBar("Spatial"),

	ui.layout_sidebar(
		ui.sidebar(

			FileSelection(
				examples={
					"seqfish.h5ad": "Example 1",
					"imc.h5ad": "Example 2",
				}, types=[".h5ad", ".h5"]
			),

			# The column that holds names for the data.
			ui.HTML("<b>Spatial Settings</b>"),
			ui.input_select(id="Key", label="Key", choices=[], multiple=False),
			ui.input_radio_buttons(id="Coordinate", label="Coordinate System", choices=["Grid", "Generic"], inline=True),
			ui.input_slider(id="Neighs", label="Neighboring Tiles", value=6, min=1, max=10, step=1),
			ui.input_slider(id="Radius", label="Radius (0 = auto)", value=0, min=0, max=10, step=0.1),
			ui.input_slider(id="Rings", label="Neighbor Rings", value=1, min=1, max=5, step=1),
			ui.input_slider(id="Percentile", label="Distance Percentile (0 = auto)", value=0.0, min=0.0, max=1.0, step=0.01),

			ui.input_select(id="Transform", label="Matrix Transform", choices={"None": "None", "spectral": "Spectral", "cosine": "Cosine"}),

			ui.input_checkbox_group(id="SpatialFeatures", label="Spatial Features", choices=["Delaunay", "Diagonal = 1.0"], selected=[]),


			ui.panel_conditional(
				"input.MainTab === 'Interactive'",
				ui.HTML("<b>Interactive Settings</b>"),
				ui.input_select(id="Statistic", label="Statistic", choices=["Moran's I", "Sepal"]),

				ui.input_select(id="ColorMap", label="Color Map", choices=["Viridis", "Plasma", "Inferno", "Magma", "Cividis"], selected="Viridis"),

				ui.input_select(id="Shape", label="Shape", choices=["Circle", "Square", "Hex"], selected="Hex"),

				ui.input_slider(id="ImgOpacity", label="Image Opacity", value=1, min=0.0, max=1.0, step=0.1),
				ui.input_slider(id="Opacity", label="Data Opacity", value=1, min=0.0, max=1.0, step=0.1),

				ui.input_checkbox_group(id="Features", label="Heatmap Features", choices=["Image", "Legend"], selected=["Image", "Legend"]),
				ui.input_checkbox_group(id="Keys", label="Annotation Keys", choices=["Lct", "Ecel1", "Cfap65", "Resp18", "Tuba4a"], selected=["Lct"]),
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

				ui.panel_conditional(
					"input.OccurrenceGraph === 'Line'",
					ui.input_select(id="Clusters", label="Clusters", choices=[], multiple=False),
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
			ui.nav_panel("Spatial Neighbors", ui.output_plot("Neighbors", height="90vh"), value="Neighbors"),
			ui.nav_panel("Ripley's Function", ui.output_plot("Ripley", height="90vh"), value="Ripley"),
			ui.nav_panel("Receptor-Ligand", ui.output_plot("ReceptorLigand"), value="ReceptorLigand"),
			ui.nav_panel("Co-occurrence", ui.output_plot("Occurrence", height="90vh"), value="Occurrence")
		),
	)
)

app = App(app_ui, server)
