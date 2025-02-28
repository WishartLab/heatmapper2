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

'''
I don’t see any legend in the example file. The Auto Update slider looks different than the other Auto Update slides.  These need to be consistent.
'''

from pandas import DataFrame
from shiny import App, reactive, render, ui
from matplotlib.pyplot import get_cmap, subplots, close as fig_close
from tempfile import TemporaryDirectory, NamedTemporaryFile
from anndata import read_h5ad
from squidpy import gr, pl, read
from scanpy import pp, tl
from pathlib import Path

# Shared functions
from shared import Cache, MainTab, NavBar, File, FileSelection, Filter, ColumnType, InitializeConfig, ColorMaps, DistanceMethods, Update, Msg, Error, Inlineify, TableOptions, TooltipIcon

try:
	from user import config
except ImportError:
	from config import config


def server(input, output, session):

	# Information regarding example files.
	Info = {
		"visium_hne_adata.h5ad": "Input type: h5ad\nContents: Pre-processed example file.\nSource: SquidPy",
	}

	InitializeConfig(config, input)


	def CreateErrorImg(text, color, file):
		"""
		@brief Generates an image of the provided text
		@param text: The text to display as an error
		@param color: Hex color code for the text
		@param inputs: A list of all the inputs for caching (from HashString())
		@returns 
		"""
		# create image with error text
		fig, ax = subplots()
		ax.text(0, 50, text, color=color, fontsize=32)
		ax.set_xlim(0, 200)
		ax.set_ylim(0, 100)
		# make axes transparent
		[ax.spines[side].set_alpha(0.0) for side in ["top", "bottom", "left", "right"]]
		ax.tick_params(axis='both', which='both', reset=False, color=[0,0,0,0], labelcolor=[0,0,0,0])
		# get image for display
		fig.savefig(file.name, format="png", dpi=100, bbox_inches="tight")
		fig_close(fig)
		img: types.ImgData = {"src": file.name, "width": "400px"}
		return img


	def HandleData(path, p=None):
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
		print(f"\nCOLUMNS TO LIST: \n{adata.obs.columns.to_list()}")
		print(f"\nadata.obs: \n{adata.obs}")
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
		'''
		Returns AnnData object with the following structure:
		obsm:	[spatial] - spatial spot coordinate matrix
		uns:	[spatial][LIBRARY ID][images] - paths to hires and lowres images
		uns:	[spatial][LIBRARY ID][scalefactors] - scalefactors for the spots
		uns:	[spatial][LIBARARY ID][metadata] - misc metadata
		'''
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
	@reactive.event(input.SourceFile, input.File, input.Example, input.Reset, input.CellCount, input.GeneCount, input.UploadType)
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

					if adata is None: return

					# Throw it into the Cache.
					DataCache.Store(adata, name)

					# Now that it's cached, remove the origin
					for file in name:
							Path(file).unlink()

				adata = DataCache.Get(name)

				if input.File()[0]["name"].endswith(".h5ad"):
					# check for cluster obs, counts obsm, ....,
					# TODO: !!! 
					print("Returning")
					Data.set(adata)
					return

				cell, gene = input.CellCount(), input.GeneCount()
				if cell is None or gene is None: return

				filtered = name + [cell, gene]
				if not DataCache.In(filtered):

					bdata = adata.copy()
					# keep cells that have at least min_counts RNA counts
					pp.filter_cells(bdata, min_counts=cell)
					# keep genes that are expressed in at least min_cells
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
				# check for cluster obs, counts obsm, ....,
				# TODO: !!! 
				Data.set(await DataCache.Load(input, default=None))
				ColumnNames(Data(), p)
				p.close()

	@output
	@render.data_frame
	def Table():
		state = config.State()
		df = Data()
		# add placeholder message if no data is uploaded
		if df is None: 
			return DataFrame({"Note": ["No data to display! Please upload your data or select an example data set in the sidebar."]})
		if state == "obs": 
			try:
				return render.DataGrid(df.obs, editable=True)
			except:
				return DataFrame({"Error": ["Observation table could not be rendered."]})
		elif state == "var": 
			try:
				return render.DataGrid(df.var, editable=True)
			except:
				return DataFrame({"Error": ["Variable table could not be rendered."]})


	@Table.set_patch_fn
	def UpdateTable(*, patch: render.CellPatch) -> render.CellValue:
		if config.Type() == "Integer": value = int(patch["value"])
		elif config.Type() == "Float": value = float(patch["value"])
		else: value = patch["value"]

		row = patch["row_index"]
		col = patch["column_index"]

		df = Data()
		table = df.obs if config.State() == "obs" else df.var
		table.iloc[row, col] = value
		Data.set(df)

		DataCache.Invalidate(File(input))

		return value


	# Info text in welcome tab
	@render.ui
	def Welcome():
		return ui.HTML("""
			<h1>Spatial heatmaps</h1>
			Spatial heatmaps display spatial molecular data, and visualize various related metrics. <br><br>
			Upload your data files in the sidebar to get started, or select 'Example' to check out a pre-loaded example. <br><br>
			Navigate to the 'Heatmap' tab to see the heatmap, 'Table' to look at the input data, or one of the statistics tabs to see a visualization of statistics.
				 
			<br><br>
			<img src="https://github.com/WishartLab/heatmapper2/wiki/assets/Spatial.png" alt="Image"; style="max-width:500px;">
				 
			<br><br><h3>Format</h3>
			<i>Input data can be uploaded in three different formats:</i><br>
			<b>1 - AnnData</b><br>
			A preprocessed file generated by AnnData with the file extension <b>.h5ad</b> can be visualized on its own without any other files. This is the format the Heatmapper will output if you download a spatial heatmap table. <br><br>
			<b>2 - Space Ranger / Visium</b><br>
			You will need to upload a folder containing the following files:
				<ul>
				<li>.h5 counts file</li>
				<li>Two tissue_{hires/lowres}_image.png images (from the <i>spatial</i> folder)</li>
				<li>scalefactors_json.json (from the <i>spatial</i> folder)</li>
				<li>tissue_positions.csv (from the <i>spatial</i> folder)</li>
				</ul>
			See <a href="https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/output-overview">here</a> for an explanation of Space Ranger output.
			<br><br>
			<b>3 - NanoString</b><br>
			You will need to upload a folder containing the following files:
			<ul>
				<li>A count file that ends in _exprMat_file.csv</li>
				<li>A meta file that ends in _metadata_file.csv</li>
				<li>An optional FOV file that contains fov in the file name.</li>
				<li>The contents of the CellLabels folder (.tif)</li>
				<li>The contents of the CellComposite folder (.png or .jpg)</li>
				</ul>
			<br>
			<br><h3>Interface</h3>
			Click on the '?' icon beside sidebar options to read more about them.
		""")


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
		# https://squidpy.readthedocs.io/en/stable/api/squidpy.pl.spatial_segment.html
		pl.spatial_segment(
			adata,
			color=count,
			library_key="fov",  # 
			seg_cell_id="cell_ID",
			library_id=id,  #
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
			dpi=dpi,
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
			dpi=dpi,
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
			
			# add placeholder message if no data uploaded
			if adata is None:
				return CreateErrorImg("No data to display!\n\nPlease upload your data or select an example data set in the sidebar.", "#027bc2", file)
			
			# generate heat maps
			try:
				if input.SourceFile() == "Example" or input.UploadType() == "Visium":
					return GenerateVisium(adata, file, p)
				elif input.UploadType() == "NanoString":
					return GenerateNanoString(adata, file, p)
			except:
				return CreateErrorImg("Spatial heat map could not be rendered.", "#027bc2", file)


	@output
	@render.image(delete_file=True)
	def Heatmap(): return GenerateHeatmap()


	@output
	@render.image(delete_file=True)
	@reactive.event(input.Update)
	def HeatmapReactive(): return GenerateHeatmap()


	@output
	@render.plot
	def Centrality():
		with ui.Progress() as p:

			p.inc(message="Loading input...")
			adata = Data()
			score = config.Score()

			# add placeholder message if no data uploaded
			if adata is None: 
				# create image with error text
				fig, ax = subplots()
				ax.text(0, 50, "No data to display! \nPlease upload your data or select an example data set in the sidebar.", color="#027bc2", fontsize=8)
				ax.set_xlim(0, 200)
				ax.set_ylim(0, 100)
				# make axes transparent
				[ax.spines[side].set_alpha(0.0) for side in ["top", "bottom", "left", "right"]]
				ax.tick_params(axis='both', which='both', reset=False, color=[0,0,0,0], labelcolor=[0,0,0,0])
				return fig

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
			# add placeholder message if no data uploaded
			if adata is None: 
				# create image with error text
				fig, ax = subplots()
				ax.text(0, 50, "No data to display! \nPlease upload your data or select an example data set in the sidebar.", color="#027bc2", fontsize=8)
				ax.set_xlim(0, 200)
				ax.set_ylim(0, 100)
				# make axes transparent
				[ax.spines[side].set_alpha(0.0) for side in ["top", "bottom", "left", "right"]]
				ax.tick_params(axis='both', which='both', reset=False, color=[0,0,0,0], labelcolor=[0,0,0,0])
				return fig

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
			# add placeholder message if no data uploaded
			if adata is None: 
				# create image with error text
				fig, ax = subplots()
				ax.text(0, 50, "No data to display! \nPlease upload your data or select an example data set in the sidebar.", color="#027bc2", fontsize=8)
				ax.set_xlim(0, 200)
				ax.set_ylim(0, 100)
				# make axes transparent
				[ax.spines[side].set_alpha(0.0) for side in ["top", "bottom", "left", "right"]]
				ax.tick_params(axis='both', which='both', reset=False, color=[0,0,0,0], labelcolor=[0,0,0,0])
				return fig

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
		if adata is None: 
			Error("The downloaded .h5ad file is empty! Please upload your data or select an example data set in the sidebar.")
			return
		temp = NamedTemporaryFile()
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
	"""),

	ui.panel_title(title=None, window_title="Spatial"),
	NavBar(),

	ui.layout_sidebar(
		ui.sidebar(
			FileSelection(
				examples={
					"visium_hne_adata.h5ad": "Example 1",
				},
				types=[".h5", ".png", ".csv", ".json", ".h5ad", ".jpg", ".tif"],
				multiple=True,
				default="Upload",
				project="Spatial"
			),

			ui.input_select(id="UploadType", label=None, choices=["Visium", "NanoString"]),

			ui.panel_conditional(
				"input.MainTab === 'TableTab'",
				config.State.UI(ui.input_select, id="State", label="Table", choices={"obs": "Observations", "var": "Variable"}),
				config.Type.UI(ui.input_radio_buttons, make_inline=False, id="Type", label="Datatype", choices=["Integer", "Float", "String"], inline=True),
				ui.input_action_button(id="Reset", label="Reset Values"),
				ui.download_button(id="DownloadTable", label="Download h5ad File"),
			),

			ui.panel_conditional(
				"input.MainTab != 'TableTab'",

				Update(),

				ui.tooltip(ui.HTML("<b>Minimum Count Filtering</b>"), "Values below the minimum count will not be displayed"),
				ui.div(
					Inlineify(ui.input_numeric, id="GeneCount", label="Min Gene Count", min=0, value=400),
					ui.popover(
						ui.span(TooltipIcon,),
						"Only display genes that are expressed in at least this many cells. Values below the minimum gene count will not be displayed. A higher minimum gene count will likely exclude more genes, but speed up rendering.",
						placement="right",
						id="GeneCount_tooltip",
					),
					style="display: inline-flex; gap: 5px;",
				),
				ui.div(
					Inlineify(ui.input_numeric, id="CellCount", label="Min Cell Count", min=0, value=100),
					ui.popover(
						ui.span(TooltipIcon,),
						"Only display cells that have at least this many RNA counts. Values below the minimum cell count will not be displayed. A higher minimum cell count will likely exclude more values, but speed up rendering.",
						placement="right",
						id="CellCount_tooltip",
					),
					style="display: inline-flex; gap: 5px;",
				),

				ui.HTML("<b>Keys</b>"),
				config.Keys.UI(ui.input_select, id="Keys", label="Annotation Keys", choices=[], selectize=True, multiple=True, tooltip=ui.HTML('Select annotation keys to plot. More than one key can be specified, with each key plotted separately and displayed next to each other. Start typing in the name of a key to search for it. <br>Read more <a href="https://squidpy.readthedocs.io/en/stable/api/squidpy.pl.spatial_scatter.html" target="_blank">here</a>.')),
				config.Count.UI(ui.input_select, id="Count", label="Nanostring Count", choices=[], tooltip="NanoString files only - Select count values to plot."),
			),

			ui.panel_conditional(
				"input.MainTab === 'HeatmapTab'",
				ui.HTML("<b>Heatmap</b>"),
				config.Statistic.UI(ui.input_select, id="Statistic", label="Visium Statistic", choices={"moran": "Moran's I", "sepal": "Sepal", "geary": "Geary's C"}, tooltip=ui.HTML('''Visium files only - Select a statistic to plot. <br>Read more: <br><a href="https://en.wikipedia.org/wiki/Moran%27s_I" target="_blank">Moran\'s I</a> <br><a href="https://academic.oup.com/bioinformatics/article/37/17/2644/6168120?login=true" target="_blank">Sepal</a> <br><a href="https://en.wikipedia.org/wiki/Geary%27s_C" target="_blank">Geary's C</a>''')),
				config.ColorMap.UI(ui.input_select, id="ColorMap", label="Color Map", choices=ColorMaps + ["Spring", "Summer", "Autumn", "Winter"], tooltip="Select a color scheme."),
				config.Shape.UI(ui.input_select, id="Shape", label="Data Shape", choices=["Circle", "Square", "Hex"], tooltip="Change the shape of each data point"),
				config.Columns.UI(ui.input_slider, id="Columns", label="# of Columns", min=1, max=10, step=1, tooltip="Specify how many plots to display side by side per row."),
				config.Spacing.UI(ui.input_slider, id="Spacing", label="Column Spacing", min=0.0, max=1.0, step=0.1, tooltip="Specify the spacing between plots."),

				ui.HTML("<b>Opacity</b>"),
				config.ImgOpacity.UI(ui.input_slider, id="ImgOpacity", label="Image Opacity", min=0.0, max=1.0, step=0.1, tooltip="Change the opacity of the background image. 1.0 indicates full opacity, while lower values make the background image more transparent."),
				config.Opacity.UI(ui.input_slider, id="Opacity", label="Data Opacity", min=0.0, max=1.0, step=0.1, tooltip="Change the opacity of the data points. 1.0 indicates full opacity, while lower values make the background image more visible."),

				ui.HTML("<b>Image Settings</b>"),
				config.Size.UI(ui.input_numeric, id="Size", label="Heatmap Size", min=1, tooltip="Change the width (in pixels) of the heatmap on your screen."),
				config.DPI.UI(ui.input_numeric, id="DPI", label="Resolution (DPI)", min=1, tooltip="Specify the resolution of the image in pixels per inch. Higher DPI values result in higher quality images, but larger file sizes. This setting affects the heatmap on screen as well as the downloaded plot."),

				ui.HTML("<b>Features</b>"),
				config.Features.UI(
					ui.input_checkbox_group, 
					make_inline=False, 
					id="Features", 
					label=None, 	
					choices=["Image", "Legend", "Frame"],
					tooltip=ui.HTML("Image toggles the visibility of the background image. <br>Legend toggles the visibility of the sidebar color legend. <br>Frame toggles the visibility of a frame around the heatmap with x and y axis titles."),
				),

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
					tooltip="Specify which centrality score to plot",
				),
			),

			ui.panel_conditional(
				"input.MainTab === 'Ripley'",
				ui.HTML("<b>Ripley</b>"),
				config.Function.UI(ui.input_select, id="Function", label="Function", choices=["L", "F", "G"]),
				config.Distance.UI(ui.input_select, id="Distance", label="Distance", choices=DistanceMethods, tooltip="Select a distance metric used to compute the function"),
			),


			ui.panel_conditional(
				"input.MainTab === 'Occurrence'",
				ui.HTML("<b>Co-Occurrence</b>"),
				config.CoCluster.UI(ui.input_select, id="CoCluster", label="Group", choices=[]),
				config.OccurrenceGraph.UI(ui.input_select, id="OccurrenceGraph", label="Graph", choices=["Scatter", "Line"], tooltip="Choose to visualize co-occurrence as a scatter or line plot"),
				config.Interval.UI(ui.input_slider, id="Interval", label="Interval", min=1, max=100, step=1, tooltip="Define the interval at which co-occurrence is computed"),
				config.Splits.UI(ui.input_slider, id="Splits", label="Splits", min=0, max=10, step=0, tooltip="Define the number of splits in which to divide spatial coordinates (if 0, Heatmapper selects a value automatically)"),
			),
			padding="10px",
			gap="20px",
			width="300px",
		),

		# Add the main interface tabs.
		MainTab(
			ui.nav_panel("Centrality Scores", ui.output_plot("Centrality", height="86vh"), value="Centrality"),
			ui.nav_panel("Ripley's Function", ui.output_plot("Ripley", height="86vh"), value="Ripley"),
			ui.nav_panel("Co-occurrence", ui.output_plot("Occurrence", height="86vh"), value="Occurrence")
		),
		height="86vh",
	)
)

app = App(app_ui, server)
