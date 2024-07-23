#
# Heatmapper
# Geocoordinate
#
# This file contains the ShinyLive application for Geocoordinate Heatmapper.
# It can be run with the following command within this directory:
#		shinylive export . [site]
# Where [site] is the destination of the site folder.
#
# If you would rather deploy the application as a PyShiny application,
# run the following command within this directory:
#		shiny run
#
#


from shiny import App, reactive, render, ui
from folium import Map as FoliumMap, Circle
from folium.plugins import HeatMap, HeatMapWithTime
from folium.raster_layers import ImageOverlay
from tempfile import NamedTemporaryFile
from matplotlib.pyplot import subplots
from scipy.stats import gaussian_kde
from scipy.io import netcdf_file
from numpy import vstack, convolve, ones, float32
from branca.colormap import LinearColormap
from pandas import DataFrame
from io import StringIO

from shared import Cache, NavBar, MainTab, FileSelection, Filter, ColumnType, TableOptions, InitializeConfig, Error, Update, Msg, File

try:
	from user import config
except ImportError:
	from config import config

# Fine, Shiny
import branca, certifi, xyzservices, requests


def server(input, output, session):

	Info = {
		"example1.txt": "This example dataset shows deaths from a cholera outbreak in 1854. John Snow used this data in conjunction with local pump locations as evidence that cholera is spread by contaminated water. A digitised version of the data is available online, courtesy of Robin Wilson (robin@rtwilson.com).",
		"example2.txt": "This example data set shows bike thefts in Vancouver in 2011. The data was obtained from a 2013 Vancouver Sun blog post by Chad Skelton.",
		"example3.txt": "This example data set shows the location of traffic signals in Toronto. The data was obtained from Toronto Open Data. The idea to use this data set comes from this R-bloggers post by Myles Harrison.",
		"example1.csv": "Random data",
		"example21.csv": "A parsed version of the Northeast and North Central Pacific hurricane database (HURDAT2) 2000-2022, available at https://www.nhc.noaa.gov/data/",
		"example3.csv": "Recorded mean temperature (F) in the USA in 2023 as measured by the EPA, available at https://aqs.epa.gov/aqsweb/airdata/FileFormats.html#_daily_summary_files",
		"test.txt": "NASA Temperature Anaomolies from 1980-2024: https://data.giss.nasa.gov/tmp/gistemp/NMAPS/tmp_GHCNv4_ERSSTv5_1200km_Anom_6_2024_2024_1951_1980_100_180_90_0_2_/amaps.txt"
	}


	def NCDataFrame(nc_path, p=None):
		"""
		@brief Takes a netcdf_file path and returns an appropriate DataFrame
		@param nc: The netcdf_file path
		@returns	A dataframe with Latitude, Longitude, and optional Time alongwith
							all variables bounded by those values.
		"""

		nc = netcdf_file(nc_path, mmap=False)
		if p: p.inc(message="Extracting keys...")
		dimensions = list(nc.dimensions.keys())
		latitude, longitude, time = Filter(dimensions, ColumnType.Latitude), Filter(dimensions, ColumnType.Longitude), Filter(dimensions, ColumnType.Time)

		if not latitude or not longitude:
			Error("NC file requires a latitude and longitude!")
			return None

		df = {
			latitude: [],
			longitude: [],
		}

		latitudes = nc.variables[latitude][:]
		longitudes = nc.variables[longitude][:]
		if time:
			df[time] = []
			times = nc.variables[time][:]
			time_index = config.TimeColumn()
			if not time_index: return None

			# You'll definitely regret this ;)
			if time_index == "All":
				time_indicies = range(len(times))
			else:
				time_indicies = [int(time_index)]
			expected_shape = (len(times), len(latitudes), len(longitudes))
		else:
			expected_shape = (len(latitudes), len(longitudes))
		variables = []
		for variable in list(nc.variables.keys()):
			if nc.variables[variable].shape == expected_shape:
				variables.append(variable)
				df[variable] = []

		for lat_index in range(len(latitudes)):
			if p: p.inc(message=f"Extracting Values at Latitude: {lat_index}...")
			for lon_index in range(len(longitudes)):
				for variable in variables:
					if time:
						for time_index in time_indicies:
							df[latitude].append(latitudes[lat_index])
							df[longitude].append(longitudes[lon_index])
							df[time].append(times[time_index])
							df[variable].append(nc.variables[variable][time_index, lat_index, lon_index])
					else:
						df[latitude].append(latitudes[lat_index])
						df[longitude].append(longitudes[lon_index])
						df[variable].append(nc.variables[variable][lat_index, lon_index])

		if p: p.inc(message="Converting into DataFrame...")
		nc.close()
		return DataFrame(df)


	def HandleData(path, p=None):
		"""
		@brief A custom Data Handler for the Cache.
		@param path: Path to the file
		@returns A data object from the cache.
		@info This Data Handler supports .nc files
		"""
		if path.suffix == ".nc":
			print(path)
			return path.resolve()
		else: return DataCache.DefaultHandler(path, p)
	DataCache = Cache("geocoordinate", DataHandler=HandleData)
	Data = reactive.value(None)
	Valid = reactive.value(False)

	InitializeConfig(config, input)


	@reactive.effect
	@reactive.event(input.SourceFile, input.File, input.Example, input.Reset, ignore_init=True)
	async def UpdateData():
		try:
			with ui.Progress() as p:
				Data.set((await DataCache.Load(input, p=p)));
				Valid.set(False)

				if File(input).endswith(".nc"):
					nc = netcdf_file(Data(), mmap=False)
					dimensions = list(nc.variables.keys())
					time = Filter(dimensions, ColumnType.Time)
					if time is not None:
						ui.update_select(id="TimeColumn", choices=list(range(0, len(nc.variables[time][:]))) + ["All"], selected=0)
					ui.update_select(id="ValueColumn", choices=dimensions + ["Uniform"])

				else:
					columns = Data().columns
					time = Filter(columns, ColumnType.Time, good=["None"], id="TimeColumn")
					value = Filter(columns, ColumnType.Value, good=["Uniform"], id="ValueColumn")
					if time == value:
						ui.update_select(id="ValueColumn", selected=columns[1] if len(columns) > 1 else None)
				DataCache.Invalidate(File(input))
		except Exception as e:
			Error(f"Failed to load file", e)


	def GetData(): return Table.data_view() if Valid() else Data()


	def GenerateMap(df, map, v_col, lon_col, lat_col):
		"""
		@brief Generates a standard heatmap
		@param df The DataFrame containing the data
		@param map The folium map to attach the heatmap to.
		"""

		opacity = config.Opacity()
		radius = config.Radius()
		blur = config.Blur()
		render = config.RenderMode()

		latitude = df[lat_col]
		longitude = df[lon_col]
		values = df[v_col]

		# Calculate kernel density estimation
		if "KDE" in config.Features():
			stack = vstack([longitude, latitude])
			kde = gaussian_kde(stack)
			density = kde(stack)
			df[v_col] = values + density * 0.1

		if render == "Raster":
			HeatMap(list(zip(df[lat_col], df[lon_col], df[v_col])),
			min_opacity=opacity,
			max_zoom=0,
			radius=radius,
			blur=blur).add_to(map)

		elif render == "Vector":
				# Define a linear colormap
				colormap = LinearColormap(
					colors=['#8000ff', '#00bfff', '#00ff80', '#ffff00', '#ff8000', '#ff0000'],
					vmin=df[v_col].min(),
					vmax=df[v_col].max()
				)

				# Add CircleMarkers to the map for each data point, applying colors based on values
				for index, row in df.iterrows():
					value = row[v_col]
					color = colormap(value)
					Circle(
						location=[row[lat_col], row[lon_col]],
						radius=radius,
						color=color,
						fill=True,
						opacity=opacity,
						fill_opacity=opacity,
						stroke=False
					).add_to(map)
		map.fit_bounds(map.get_bounds())


	def GenerateTemporalMap(df, map, t_col, v_col, lon_col, lat_col):
		"""
		@brief Generates a temporal heatmap
		@param df The DataFrame containing the data
		@param map The folium map to attach the heatmap to.
		"""

		# Ensure we have a valid time column
		if t_col not in df or v_col not in df: return

		# Sort by time so we can work linearly.
		df = df.sort_values(by=t_col)

		# Normalize
		if v_col != "Uniform":
			values = df[v_col]
			df[v_col] = (values - values.min()) / (values.max() - values.min())

		# Group data by time
		data = []
		for time, group_df in df.groupby(t_col):
			time_slice = []
			for _, row in group_df.iterrows():
				lat = row[lat_col]
				lon = row[lon_col]
				value = 1 if v_col == "Uniform" else row[v_col]
				time_slice.append([lat, lon, value])
			data.append(time_slice)

		radius = config.Radius() // 2
		opacity = config.Opacity()
		blur = config.Blur() / 30

		# Make the heamap
		HeatMapWithTime(
			data,
			index=df[t_col].drop_duplicates().to_list(),
			radius=radius,
			min_opacity=opacity,
			blur=blur,
			max_speed=60).add_to(map)


	@output
	@render.data_frame
	def Table():
		df = Data()
		if File(input).endswith(".nc"):
			df = NCDataFrame(df)

		try:
			Valid.set(True); return render.DataGrid(df, editable=True)
		except Exception as e:
			Error(f"Failed to render table", e)


	@Table.set_patch_fn
	def UpdateTable(*, patch: render.CellPatch) -> render.CellValue:
		if config.Type() == "Integer": value = int(patch["value"])
		elif config.Type() == "Float": value = float(patch["value"])
		else: value = patch["value"]
		DataCache.Invalidate(File(input))
		return value


	def GenerateHeatmap():
		with ui.Progress() as p:
			p.inc(message="Loading input...")
			df = GetData()
			if df is None: return

			if File(input).endswith(".nc") and not Valid():
				df = NCDataFrame(GetData())
				if df is None: return None
				nc = True
			else:
				df = df.copy(deep=True)
				nc = False

			p.inc(message="Formatting...")

			lon_col = Filter(df.columns, ColumnType.Longitude)
			lat_col = Filter(df.columns, ColumnType.Latitude)
			if lat_col is None or lon_col is None: return

			v_col = config.ValueColumn()
			if v_col == "Uniform":
				df["Heatmapper_Uniform_Values"] = [1] * len(df[lat_col])
				v_col = "Heatmapper_Uniform_Values"

			t_col = config.TimeColumn()

			map = FoliumMap((df[lat_col][0], df[lon_col][0]), tiles=config.MapType())

			p.inc(message="Dropping Invalid Values...")
			if config.ROI():
				to_drop = []
				l, u = config.Min(), config.Max()
				for index, value in zip(df.index, df[v_col]):
					if value < l or value > u:
						if config.ROI_Mode() == "Remove": to_drop.append(index)
						elif config.ROI_Mode() == "Round": df.at[index, v_col] = u if value > u else l
				df = df.drop(to_drop)
				if len(df) == 0:
					Error("No locations! Ensure Key Column and Key Properties are correct, and your ROI is properly set!")
					return

			if "Smoothing" in config.Features():
				p.inc(message="Smoothing...")
				df = df.sort_values(by=[lat_col, lon_col])
				df[v_col] = convolve(df[v_col], ones(5) / 5 , mode='same')

				"""
				# Reduce the amount of data
				df[lon_col] = df[lon_col].round(3)
				df[lat_col] = df[lat_col].round(3)

				df = df.sort_values(by=[lon_col, lat_col])
				grid = df.pivot(index=lat_col, columns=lon_col, values=v_col)

				conv_kernel = ones((3, 3)) / 9
				smoothed_grid = convolve2d(grid, conv_kernel, mode='same', boundary='fill', fillvalue=0)

				smoothed_df = DataFrame(smoothed_grid, index=grid.index, columns=grid.columns)
				smoothed_df = smoothed_df.reset_index()
				smoothed_df = melt(smoothed_df, id_vars=lat_col, var_name=lon_col, value_name=v_col)
				smoothed_df = smoothed_df.sort_values(by=[lon_col, lat_col])
				df = smoothed_df.reset_index(drop=True)
				"""

			# Generate the right heatmap.
			p.inc(message="Plotting...")
			if t_col != "None" and not nc: GenerateTemporalMap(df, map, t_col, v_col, lon_col, lat_col)
			else: GenerateMap(df, map, v_col, lon_col, lat_col)
			return map


	@output
	@render.ui
	def Heatmap():
		try:
			return GenerateHeatmap()
		except Exception as e:
			Error(f"Failed to generate heatmap", e)


	@output
	@render.ui
	@reactive.event(input.Update)
	def HeatmapReactive():
		try:
			return GenerateHeatmap()
		except Exception as e:
			Error(f"Failed to generate heatmap", e)
			raise e


	@reactive.effect
	@reactive.event(input.ExampleInfo)
	def ExampleInfo():
		Msg(ui.HTML(Info[input.Example()]))


	@render.download(filename="table.csv")
	def DownloadTable(): yield GetData().to_string()


	@render.download(filename="heatmap.html")
	def DownloadHeatmap(): m = LoadMap(); yield m.get_root().render()


app_ui = ui.page_fluid(

	NavBar(),

	ui.layout_sidebar(
		ui.sidebar(

			FileSelection(
				examples={
					"example1.txt": "Example 1",
					"example2.txt": "Example 2",
					"example3.txt": "Example 3",
					"example1.csv": "Example 4",
					"example21.csv": "Example 5",
					"example3.csv": "Example 6",
					"test.txt": "Example 7",
				},
				types=[".csv", ".txt", ".dat", ".tsv", ".tab", ".xlsx", ".xls", ".odf", ".nc"],
				project="Geocoordinate"
			),

			TableOptions(config),

			ui.panel_conditional(
				"input.MainTab === 'HeatmapTab'",

				Update(),

				ui.HTML("<b>Columns</b>"),
				config.TimeColumn.UI(ui.input_select, id="TimeColumn", label="Time", choices=[], multiple=False),
				config.ValueColumn.UI(ui.input_select, id="ValueColumn", label="Value", choices=[], multiple=False),

				ui.HTML("<b>Heatmap</b>"),
				config.RenderMode.UI(ui.input_select, id="RenderMode", label="Render", choices=["Raster", "Vector"]),
				config.MapType.UI(ui.input_select,id="MapType", label="Map", choices={"CartoDB Positron": "CartoDB", "OpenStreetMap": "OSM"}),

				config.Radius.UI(ui.input_numeric, id="Radius", label="Size", min=5),

				config.Opacity.UI(ui.input_numeric, id="Opacity", label="Opacity", min=0.0, max=1.0, step=0.01),
				config.Blur.UI(ui.input_numeric, id="Blur", label="Blurring", min=1, max=30, step=1),


				ui.HTML("<b>Range of Interest</b>"),
				config.ROI.UI(ui.input_checkbox, make_inline=False, id="ROI", label="Enable (Lower/Upper)"),
				config.ROI_Mode.UI(ui.input_radio_buttons, make_inline=False, id="ROI_Mode", label=None, choices=["Remove", "Round"], inline=True),
				ui.layout_columns(
					config.Min.UI(ui.input_numeric,make_inline=False, id="Min", label=None, min=0),
					config.Max.UI(ui.input_numeric, make_inline=False, id="Max", label=None, min=0),
				),

				ui.HTML("<b>Features</b>"),
				config.Features.UI(
					ui.input_checkbox_group, id="Features", make_inline=False, label=None,
					choices=["Smoothing", "KDE"],
				),

				# Add the download buttons.
				ui.download_button("DownloadHeatmap", "Heatmap")
			),
			padding="10px",
			gap="20px",
			width="250px",
		),

		MainTab(m_type=ui.output_ui),
		height="90vh",
	)
)

app = App(app_ui, server)
