#
# Heatmapper
# Spectral
#
# This file contains the ShinyLive application for Spectral Heatmapper.
# It can be run with the following command within this directory:
#		shinylive export . [site]
# Where [site] is the destination of the site folder.
#
# If you would rather deploy the application as a PyShiny application,
# run the following command within this directory:
#		shiny run
#
#

# ShinyLive needs this imported
import regex

from shiny import App, reactive, render, ui
from matplotlib.pyplot import subplots, colorbar, style, get_cmap
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from tempfile import NamedTemporaryFile
from io import BytesIO
from pymzml.run import Reader
from pandas import DataFrame
from scipy.spatial.distance import squareform
from scipy.interpolate import griddata
from numpy import zeros, unique, array, concatenate, asarray, hstack, column_stack, newaxis, full_like, zeros_like, meshgrid, full, where, linspace, log1p

from shared import Cache, MainTab, NavBar, FileSelection, Filter, ColumnType, TableOptions, InitializeConfig, ColorMaps, Update, Msg, File, InterpolationMethods, Error

try:
	from user import config
except ImportError:
	from config import config


def server(input, output, session):

	# Information regarding example files.
	Info = {
		"1min.mzml": "An Example mzML from https://github.com/HUPO-PSI/mzML"
	}


	def Hash():
		"""
		@brief Compute a Hash String for the Main Heatmap
		@return A hash of all the inputs used for the heatmap.
		"""
		tab = input.MainTab()
		if tab == "HeatmapTab":
			return [
				File(input),
				config.ColorMap(),
				config.Features(),
				config.TextSize(),
				config.ID(),
				config.Peaks(),
				config.DPI(),
				config.Elevation(),
				config.Zoom(),
				config.Rotation(),
				config.Dimension(),
				input.mode(),
			]
		elif tab == "SimilarityTab":
			return [
				File(input),
				config.ColorMap(),
				config.Features(),
				config.TextSize(),
				config.ID(),
				config.DPI(),
				config.Interpolation(),
				input.mode(),
			]


	def HandleData(path, p=None):
		"""
		@brief A custom Data Handler for the Cache.
		@param path: The path to the file
		@returns A data object from the cache.
		@info This Data Handler supports .mzML files
		"""
		if path.suffix.lower() == ".mzml": return Reader(path.resolve(), build_index_from_scratch=True)
		else: return None
	DataCache = Cache("spectral", DataHandler=HandleData)
	Data = reactive.value(None)
	Valid = reactive.value(False)

	InitializeConfig(config, input)


	@reactive.effect
	@reactive.event(input.SourceFile, input.File, input.Example, input.Reset)
	async def UpdateData():
		Data.set((await DataCache.Load(input, p=ui.Progress(), default=None)))
		Valid.set(False)
		DataCache.Invalidate(File(input))

		reader = Data()
		if reader is None: return

		ids = set()
		first = None
		for spectra in reader:
			ids.add(spectra.ID)
			if first is None: first = spectra.ID
		ui.update_select(id="ID", selected=[first], choices=list(ids))


	def GetData(): return Table.data_view() if Valid() else Data()


	@output
	@render.data_frame
	def Table(): Valid.set(True); return render.DataGrid(Data(), editable=True)


	@Table.set_patch_fn
	def UpdateTable(*, patch: render.CellPatch) -> render.CellValue:
		if config.Type() == "Integer": value = int(patch["value"])
		elif config.Type() == "Float": value = float(patch["value"])
		else: value = patch["value"]
		DataCache.Invalidate(File(input))
		return value


	def GenerateSimilarity():
		"""
		@brief Generate the Similarity matrix.
		"""
		
		# Sometimes the MainTab doesn't update immediately, despite calling the function.
		# We just return nothing if we aren't actually on the tab to avoid redundant calculation.
		if input.MainTab() != "SimilarityTab": return

		inputs = Hash()
		if not DataCache.In(inputs):
			with ui.Progress() as p:
				p.inc(message="Loading input...")

				reader = GetData()
				if reader is None: return

				# Get all the spectra the user wants.
				distances = {}
				indices = [int(i) for i in config.ID()]
				spectra = []
				for s in reader:
					if s.ID in indices: spectra.append(s)

				for s in spectra:
					distances[s.ID] = {}
					for s2 in spectra:
						p.inc(message=f"Computing similarity of spectra {s.ID} and {s2.ID}")
						
						# If they're the same spectra, set it to 1.
						if s.ID == s2.ID:
							distances[s.ID][s2.ID] = 1.0
						
						# If the second spectra exists, used that pre-calculated result.
						elif s2.ID in distances:
							distances[s.ID][s2.ID] = distances[s2.ID][s.ID]
							
						# Otherwise call the similarity function.
						else:
							distances[s.ID][s2.ID] = s.similarity_to(s2)

				p.inc(message="Plotting")
				df = DataFrame(distances, columns=indices, index=indices)
				fig, ax = subplots()
				interpolation = config.Interpolation().lower()
				plot = ax.imshow(df, cmap=config.ColorMap().lower(), interpolation=interpolation, aspect="equal")

				# Visibility of features
				try:
					if "legend" in input.Features():
						cbar = colorbar(plot, ax=ax, label="Value", pad=0.1)
						cbar.ax.tick_params(labelsize=config.TextSize())
				except Exception: pass

				if "y" in config.Features():
					ax.set_yticklabels(df.columns)
					ax.set_yticks(range(len(df.columns)))
					ax.tick_params(axis="y", labelsize=config.TextSize())
				else: ax.set_yticklabels([])

				if "x" in config.Features():
					ax.set_xticklabels(df.columns)
					ax.set_xticks(range(len(df.columns)))
					ax.tick_params(axis="x", labelsize=config.TextSize())
				else: ax.set_xticklabels([])

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


	def GenerateHeatmap():
		"""
		@brief Generate the main heatmap
		"""
		
		
		if input.MainTab() != "HeatmapTab": return

		inputs = Hash()
		if not DataCache.In(inputs):
			with ui.Progress() as p:
				p.inc(message="Loading input...")
				reader = GetData()
				if reader is None: return

				fig, ax = subplots(subplot_kw={"projection": "3d"})
				cmap = get_cmap(config.ColorMap().lower())

				# We additionally cache interpolation.
				interpolation_cache = [File(input), config.Peaks(), config.Dimension(), "Interpolation"]
				if not DataCache.In(interpolation_cache):

					peaks = config.Peaks().lower()
					
					x_min, x_max, y_min, y_max, z_min, z_max = 0, 0, 0, 0, 0, 0
					values = []
					intensities = []
					rts = []

					# For each spectra, get its values intensities, and rt.
					for spectrum in reader:
						p.inc(message=f"Reading Spectra {spectrum.ID}")
						rt = spectrum.scan_time[0]
						for mz, intensity in spectrum.peaks(peaks):
							values.append(mz)
							intensities.append(intensity)
							rts.append(rt)
					if not values:
						Error("No Spectra in File!")
						return

					# Get the min and max of each list.
					vm, vM, rm, rM = min(values), max(values), min(rts), max(rts)

					# Create a grid for mz and rt
					p.inc(message="Interpolating")
					dimension = config.Dimension()
					mz_grid, rt_grid = meshgrid(
						linspace(vm, vM, dimension),
						linspace(rm, rM, dimension)
					)

					# Interpolate.
					intensity_grid = griddata((values, rts), intensities, (mz_grid, rt_grid), method='cubic')
					intensity_grid[intensity_grid < 0] = 0

					DataCache.Store([mz_grid, rt_grid, intensity_grid, vm, vM, rm, rM], interpolation_cache)

				else:
					mz_grid, rt_grid, intensity_grid, vm, vM, rm, rM = DataCache.Get(interpolation_cache)

				p.inc(message="Plotting")
				plot = ax.plot_surface(mz_grid, rt_grid, intensity_grid, cmap=cmap)

				ax.set_xlim(vm, vM)
				ax.set_ylim(rm, rM)
				
				# This causes issues.
				#ax.set_zlim(0, iM)

				ax.view_init(elev=config.Elevation(), azim=config.Rotation())
				ax.set_box_aspect(None, zoom=config.Zoom())

				if plot is None: return

				# Visibility of features
				try:
					if "legend" in input.Features():
						cbar = colorbar(plot, ax=ax, label="Value", pad=0.2)
						cbar.ax.tick_params(labelsize=config.TextSize())
				except Exception: pass

				if "y" in config.Features():
					ax.tick_params(axis="y", labelsize=config.TextSize())
					ax.set_ylabel("Retention Time")
				else: ax.set_yticklabels([])

				if "x" in config.Features():
					ax.tick_params(axis="x", labelsize=config.TextSize())
					ax.set_xlabel("MZ")
				else: ax.set_xticklabels([])

				if "z" in config.Features():
					ax.tick_params(axis="z", labelsize=config.TextSize())
					ax.set_zlabel("Intensity")
				else:
					ax.set_zticklabels([])

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
	def Similarity(): return GenerateSimilarity()


	@output
	@render.image(delete_file=True)
	@reactive.event(input.Update)
	def HeatmapReactive(): return GenerateHeatmap()


	@reactive.effect
	@reactive.event(input.ExampleInfo)
	def ExampleInfo():
		Msg(ui.HTML(Info[input.Example()]))


	@render.download(filename="table.csv")
	def DownloadTable(): yield GetData().to_string()


	@render.download(filename="heatmap.png")
	def DownloadHeatmap(): yield DataCache.Get(Hash())


app_ui = ui.page_fluid(

	NavBar(),

	ui.layout_sidebar(
		ui.sidebar(

			FileSelection(examples={"1min.mzml": "Example 1"}, types=[".mzml", ".mzML"], project="Spectral"),

			TableOptions(config),

			ui.panel_conditional(
				"input.MainTab != 'TableTab'",

				Update(),


				ui.HTML("<b>Heatmap</b>"),

				config.ID.UI(ui.input_select, id="ID", label="ID", selectize=True, multiple=True, choices=[0], conditional="input.MainTab === 'SimilarityTab'"),

				config.TextSize.UI(ui.input_numeric, id="TextSize", label="Text", min=1, max=50, step=1),
				config.ColorMap.UI(ui.input_select, id="ColorMap", label="Map", choices=ColorMaps),

				config.Peaks.UI(ui.input_select, id="Peaks", label="Peak Type", choices=["Raw", "Centroided", "Reprofiled"], conditional="input.MainTab === 'HeatmapTab'"),

				config.Interpolation.UI(ui.input_select, id="Interpolation", label="Inter", choices=InterpolationMethods, conditional="input.MainTab === 'SimilarityTab'"),

				config.Dimension.UI(ui.input_numeric, id="Dimension", label="Contour Size", conditional="input.MainTab === 'HeatmapTab'", min=1),


				ui.HTML("<b>3D</b>"),
				config.Elevation.UI(ui.input_numeric, id="Elevation",	label="Elevation"),
				config.Rotation.UI(ui.input_numeric, id="Rotation",	label="Rotation", step=1, min=1),
				config.Zoom.UI(ui.input_numeric, id="Zoom",	label="Zoom", step=1, min=1),


				ui.HTML("<b>Image Settings</b>"),
				config.Size.UI(ui.input_numeric, id="Size", label="Size", min=1),
				config.DPI.UI(ui.input_numeric, id="DPI", label="DPI", min=1),

				# Customize what aspects of the heatmap are visible
				ui.HTML("<b>Features</b>"),
				config.Features.UI(ui.input_checkbox_group, make_inline=False, id="Features", label=None,
						choices={"x": "X Labels", "y": "Y Labels", "z": "Z Labels", "legend": "Legend"}
				),

				ui.download_button(id="DownloadHeatmap", label="Download"),
			),
			padding="10px",
			gap="20px",
			width="250px",
		),

		# Add the main interface tabs.
		MainTab(
			ui.nav_panel("Similarity", ui.output_plot("Similarity", height="90vh"), value="SimilarityTab"),
		),
	)
)

app = App(app_ui, server)
