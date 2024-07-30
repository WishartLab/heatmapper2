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
from numpy import zeros, unique, array, concatenate, asarray, hstack, column_stack, newaxis, full_like

from shared import Cache, MainTab, NavBar, FileSelection, Filter, ColumnType, TableOptions, InitializeConfig, ColorMaps, Update, Msg, File, InterpolationMethods

try:
	from user import config
except ImportError:
	from config import config


def server(input, output, session):

	# Information regarding example files.
	Info = {
		"1min.mzml": "An Example mzML from https://github.com/HUPO-PSI/mzML"
	}

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
		ui.update_select(id="Index", selected=[0], choices=list(range(reader.get_spectrum_count())))


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


	def ColoredLine(x, y, c, ax, z=None, **lc_kwargs):
		"""
		@brief Generate a MatPlotLib X/Y axis, colored with list c.
		@param x The X values
		@param y The Y values
		@pararm c The color values.
		@param ax The axis
		@param lc_kwargs Keyword arguments to be passed to the LineCollection
		@info Borrowed from https://matplotlib.org/stable/gallery/lines_bars_and_markers/multicolored_line.html
		"""

		# Default the capstyle to butt so that the line segments smoothly line up
		default_kwargs = {"capstyle": "butt"}
		default_kwargs.update(lc_kwargs)

		cmap = get_cmap(config.ColorMap().lower())
		norm = Normalize(min(c), max(c))
		colors = cmap(norm(c))

		# Compute the midpoints of the line segments. Include the first and last points
		# twice so we don't need any special syntax later to handle them.
		x, y = asarray(x), asarray(y)

		x_midpts = hstack((x[0], 0.5 * (x[1:] + x[:-1]), x[-1]))
		y_midpts = hstack((y[0], 0.5 * (y[1:] + y[:-1]), y[-1]))

		if z is None:
			coord_start = column_stack((x_midpts[:-1], y_midpts[:-1]))[:, newaxis, :]
			coord_mid = column_stack((x, y))[:, newaxis, :]
			coord_end = column_stack((x_midpts[1:], y_midpts[1:]))[:, newaxis, :]
			segments = concatenate((coord_start, coord_mid, coord_end), axis=1)
			lc = LineCollection(segments, linewidths=[config.Width()] * len(x), colors=colors, **default_kwargs)

		else:
			z = asarray(z)
			z_midpts = hstack((z[0], 0.5 * (z[1:] + z[:-1]), z[-1]))
			coord_start = column_stack((x_midpts[:-1], y_midpts[:-1], z_midpts[:-1]))[:, newaxis, :]
			coord_mid = column_stack((x, y, z))[:, newaxis, :]
			coord_end = column_stack((x_midpts[1:], y_midpts[1:], z_midpts[1:]))[:, newaxis, :]
			segments = concatenate((coord_start, coord_mid, coord_end), axis=1)
			lc = Line3DCollection(segments, linewidths=[config.Width()] * len(x), colors=colors, **default_kwargs)


		lc.set_array(c)	# set the colors of each segment
		ax.add_collection(lc)
		return lc


	def SpectraHeatmap(reader, p):

		peaks = config.Peaks().lower()
		indicies = [int(i) for i in config.Index()]

		if len(indicies) == 1:
			spectra = indicies[0]

			# Collect spectrum attributes
			mz_values = []
			intensities = []


			# This is the "Pythonic" way of doing things, because
			# Heaven forbid we just index the values we want.
			for i, spectrum in enumerate(reader):
				print(spectrum.scan_time)
				if i == spectra:
					for mz, intensity in spectrum.peaks(peaks):
						mz_values.append(mz)
						intensities.append(intensity)
					break

			p.inc(message="Plotting...")
			color = input.mode()
			with style.context('dark_background' if color == "dark" else "default"):
				fig, ax = subplots()
				plot = ColoredLine(mz_values, intensities, intensities, ax)
				ax.set_xlim(min(mz_values), max(mz_values))
				ax.set_ylim(min(intensities), max(intensities))

				return fig, ax, plot

		else:
			fig, ax = subplots(subplot_kw={"projection": "3d"})

			x_min, x_max, y_min, y_max, z_min, z_max = 0, 0, 0, len(indicies), 0, 0

			for i, spectrum in enumerate(reader):
				if i in indicies:
					mz_values = []
					intensities = []
					for mz, intensity in spectrum.peaks(peaks):
						mz_values.append(mz)
						intensities.append(intensity)

					x_min = min(x_min, min(mz_values))
					x_max = max(x_max, max(mz_values))
					z_min = min(z_min, min(intensities))
					z_max = max(z_max, max(intensities))

					position = [float(indicies.index(i))] * len(mz_values)

					cmap = get_cmap(config.ColorMap().lower())

					# Plot the line for each spectrum
					plot = ColoredLine(mz_values, position, intensities, ax, cmap=config.ColorMap().lower(), z=intensities)

			ax.set_xlim(x_min, x_max)
			ax.set_ylim(y_min, y_max)
			ax.set_zlim(z_min, z_max)

			if "z" in config.Features():
				ax.tick_params(axis="z", labelsize=config.TextSize())
			else:
				ax.set_zticklabels([])
			return fig, ax, plot


	def GenerateSimilarity():
		inputs = [
			File(input),
			config.ColorMap(),
			config.Features(),
			config.TextSize(),
			config.Index(),
			config.DPI(),
			config.Interpolation(),
			input.mode(),
		]

		if not DataCache.In(inputs):
			with ui.Progress() as p:
				p.inc(message="Loading input...")

				reader = GetData()
				if reader is None: return

				distances = {}

				indices = [int(i) for i in config.Index()]

				for i, s in enumerate(reader):
					if i in indices:
						distances[i] = {}
						for x, s2 in enumerate(reader):
							if x in indices:
								if x == i:
									distances[i][x] = 1.0
								elif x in distances:
									distances[i][x] = distances[x][i]
								else:
									distances[i][x] = s.similarity_to(s2)

				df = DataFrame(distances)

				fig, ax = subplots()
				interpolation = config.Interpolation().lower()
				plot = ax.imshow(df, cmap=config.ColorMap().lower(), interpolation=interpolation, aspect="equal")

				# Visibility of features
				try:
					if "legend" in input.Features():
						cbar = colorbar(plot, ax=ax, label="Value", pad=0.1)
						cbar.ax.tick_params(labelsize=config.TextSize())
				except Exception: pass

				if "y" in config.Features(): ax.tick_params(axis="y", labelsize=config.TextSize())
				else: ax.set_yticklabels([])

				if "x" in config.Features(): ax.tick_params(axis="x", labelsize=config.TextSize())
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
		inputs = [
			File(input),
			config.ColorMap(),
			config.Opacity(),
			config.Features(),
			config.TextSize(),
			config.Index(),
			config.Peaks(),
			config.DPI(),
			config.Width(),
			input.mode(),
		]

		if not DataCache.In(inputs):
			with ui.Progress() as p:
				p.inc(message="Loading input...")
				reader = GetData()
				if reader is None: return

				fig, ax, plot = SpectraHeatmap(reader, p)
				if plot is None: return

				# Visibility of features
				try:
					if "legend" in input.Features():
						cbar = colorbar(plot, ax=ax, label="Value", pad=0.1)
						cbar.ax.tick_params(labelsize=config.TextSize())
				except Exception: pass

				if "y" in config.Features(): ax.tick_params(axis="y", labelsize=config.TextSize())
				else: ax.set_yticklabels([])

				if "x" in config.Features(): ax.tick_params(axis="x", labelsize=config.TextSize())
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
	def DownloadHeatmap():
		yield DataCache.Get([
			File(input),
			config.ColorMap(),
			config.Opacity(),
			config.Features(),
			config.TextSize(),
			config.Index(),
			config.Peaks(),
			config.DPI(),
			config.Width(),
			input.mode(),
			input.MainTab()
		])


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

				config.Index.UI(ui.input_select, id="Index", label="Index", selectize=True, multiple=True, choices=[0]),


				config.TextSize.UI(ui.input_numeric, id="TextSize", label="Text", min=1, max=50, step=1),
				config.ColorMap.UI(ui.input_select, id="ColorMap", label="Map", choices=ColorMaps),
				config.Opacity.UI(ui.input_numeric, id="Opacity", label="Opacity", min=0.0, max=1.0, step=0.1),

				config.Peaks.UI(ui.input_select, id="Peaks", label="Peak Type", choices=["Raw", "Centroided", "Reprofiled"], conditional="input.MainTab === 'HeatmapTab'"),
				config.Width.UI(ui.input_numeric, id="Width", label="Width", min=0.1, max=10.0, step=0.1, conditional="input.MainTab === 'HeatmapTab'"),

				config.Interpolation.UI(ui.input_select, id="Interpolation", label="Inter", choices=InterpolationMethods, conditional="input.MainTab === 'SimilarityTab'"),


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
