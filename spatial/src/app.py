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

import squidpy as sq
import scanpy as sp

# Shared functions
from shared import Cache, MainTab, NavBar, FileSelection, Filter, ColumnType, TableValueUpdate


def server(input, output, session):

	# Information regarding example files.
	Info = {
		"example1.h5ad": "Adult Mouse Brain Section 2, from https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Adult_Mouse_Brain_Coronal_Section_2?",
	}


	def HandleData(n, i):
		"""
		@brief A custom Data Handler for the Cache.
		@param n: The name of the file
		@param i: The source of the file. It can be a path to a file (string) or a BytesIO object.
		@returns A data object from the cache.
		@info This Data Handler supports object files, and images as textures.
		"""

		suffix = Path(n).suffix
		match suffix:
			case ".h5ad":
				if type(i) is BytesIO:
					temp = NamedTemporaryFile(suffix=suffix); temp.write(i.read())
					i = temp.name
				return sp.read_h5ad(i)
			case _: return DataCache.DefaultHandler(n, i)
	DataCache = Cache("spatial", DataHandler=HandleData)


	def GenerateMoran(adata):
		genes = adata[:, adata.var.highly_variable].var_names.values[:100]
		sq.gr.spatial_neighbors(adata)
		sq.gr.spatial_autocorr(
				adata,
				mode="moran",
				genes=genes,
				n_perms=100,
				n_jobs=1,
		)
		return adata


	def GenerateSepal(adata):
		sq.gr.spatial_neighbors(adata)
		genes = adata.var_names[(adata.var.n_cells > 100) & adata.var.highly_variable][0:100]
		sq.gr.sepal(adata, max_neighs=6, genes=genes, n_jobs=1)
		return adata


	@output
	@render.data_frame
	@reactive.event(input.SourceFile, input.File, input.Example, input.File, input.Update, input.Reset)
	async def LoadedTable(): return (await DataCache.Load(input)).to_df()


	@output
	@render.plot
	@reactive.event(input.SourceFile, input.File, input.Example, input.Update, input.Reset, input.Statistic, input.Shape, input.Features, input.Opacity, input.ColorMap, input.Keys)
	def Heatmap():
		data = DataCache.SyncLoad(input, default=None)
		if data is None: return

		match input.Statistic():
			case "Moran's I": adata = GenerateMoran(data)
			case "Sepal": adata = GenerateSepal(data)

		sq.pl.spatial_scatter(
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
	@reactive.event(input.SourceFile, input.File, input.Example, input.Update, input.Reset, input.ClusterKey)
	def Centrality():
		adata = DataCache.SyncLoad(input, default=None)
		if adata is None: return
		sq.gr.spatial_neighbors(adata)
		sq.gr.centrality_scores(adata, input.ClusterKey())
		sq.pl.centrality_scores(adata, input.ClusterKey())


	@output
	@render.plot
	@reactive.event(input.SourceFile, input.File, input.Example, input.Update, input.Reset, input.Key, input.ConnectKey)
	def Neighbors():
		adata = DataCache.SyncLoad(input, default=None)
		if adata is None: return

		sq.gr.spatial_neighbors(adata, n_rings=2, coord_type="grid", n_neighs=6)
		_, idx = adata.obsp[input.ConnectKey()][420, :].nonzero()
		idx = append(idx, 420)
		sq.pl.spatial_scatter(
				adata[idx, :],
				shape=None,
				color=input.Key(),
				connectivity_key=input.ConnectKey(),
				size=100,
		)


	@output
	@render.plot
	@reactive.event(input.SourceFile, input.File, input.Example, input.Update, input.Reset, input.Function)
	def Ripley():
		adata = DataCache.SyncLoad(input, default=None)
		if adata is None: return

		sq.gr.ripley(adata, cluster_key=input.ClusterKey(), mode=input.Function())
		sq.pl.ripley(adata, cluster_key=input.ClusterKey(), mode=input.Function())


	@output
	@render.plot
	@reactive.event(input.SourceFile, input.File, input.Example, input.Update, input.Reset)
	def ReceptorLigand():
		adata = DataCache.SyncLoad(input, default=None)
		if adata is None: return

		res = sq.gr.ligrec(
			adata,
			n_perms=1000,
			cluster_key="celltype_mapped_refined",
			copy=True,
			use_raw=False,
			transmitter_params={"categories": "ligand"},
			receiver_params={"categories": "receptor"},
		)
		sq.pl.ligrec(res, source_groups="Erythroid", alpha=0.005)


	@output
	@render.plot
	@reactive.event(input.SourceFile, input.File, input.Example, input.Update, input.Reset, input.OccurrenceGraph)
	def Occurrence():
		adata = DataCache.SyncLoad(input, default=None)
		if adata is None: return

		sq.gr.spatial_neighbors(adata)
		sq.gr.co_occurrence(adata, cluster_key=input.ClusterKey())

		if input.OccurrenceGraph() == "Line":
			sq.pl.co_occurrence(adata, cluster_key=input.ClusterKey(), clusters="basal CK tumor cell")
		else:
			sq.pl.spatial_scatter(adata, color=input.ClusterKey(), size=10, shape=None)


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


app_ui = ui.page_fluid(

	NavBar("Spatial"),

	ui.layout_sidebar(
		ui.sidebar(

			FileSelection(
				examples={
					"seqfish.h5ad": "Example 1",
					"imc.h5ad": "Example 2",
				}, types=[".h5ad"]
			),

			ui.input_text("ClusterKey", "Cluster Key", value="cell type"),

			ui.panel_conditional(
				"input.MainTab === 'Interactive'",
				ui.input_select(id="Statistic", label="Statistic", choices=["Moran's I", "Sepal"]),

				ui.input_select(id="ColorMap", label="Color Map", choices=["Viridis", "Plasma", "Inferno", "Magma", "Cividis"], selected="Viridis"),

				ui.input_select(id="Shape", label="Shape", choices=["Circle", "Square", "Hex"], selected="Hex"),

				ui.input_slider(id="ImgOpacity", label="Image Opacity", value=1, min=0.0, max=1.0, step=0.1),
				ui.input_slider(id="Opacity", label="Data Opacity", value=1, min=0.0, max=1.0, step=0.1),

				ui.input_checkbox_group(id="Features", label="Heatmap Features", choices=["Image", "Legend"], selected=["Image", "Legend"]),
				ui.input_checkbox_group(id="Keys", label="Annotation Keys", choices=["Lct", "Ecel1", "Cfap65", "Resp18", "Tuba4a"], selected=["Lct"]),
			),

			ui.panel_conditional(
				"input.MainTab === 'Centrality'",
			),

			ui.panel_conditional(
				"input.MainTab === 'Neighbors'",
				ui.input_text("Key", "Key", value="cell type"),
				ui.input_text("ConnectKey", "Connectivity Key", value="spatial_connectivities"),
			),


			ui.panel_conditional(
				"input.MainTab === 'Ripley'",
				ui.input_select(id="Function", label="Function", choices=["L", "F", "G"]),
			),

			ui.panel_conditional(
				"input.MainTab === 'Occurrence'",
				ui.input_radio_buttons(id="OccurrenceGraph", label="Graph Type", choices=["Line", "Scatter"], inline=True),
			),


			# Add the download buttons.
			ui.download_button("DownloadTable", "Download Table"),
		),

		# Add the main interface tabs.
		MainTab(
			ui.nav_panel("Centrality Scores", ui.output_plot("Centrality", height="75vh"), value="Centrality"),
			ui.nav_panel("Spatial Neighbors", ui.output_plot("Neighbors", height="75vh"), value="Neighbors"),
			ui.nav_panel("Ripley's Function", ui.output_plot("Ripley", height="75vh"), value="Ripley"),
			ui.nav_panel("Receptor-Ligand", ui.output_plot("ReceptorLigand", height="75vh"), value="ReceptorLigand"),
			ui.nav_panel("Co-occurrence", ui.output_plot("Occurrence", height="75vh"), value="Occurrence")
		),
	)
)

app = App(app_ui, server)
