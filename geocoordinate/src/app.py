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


from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from folium import Map as FoliumMap
from folium.plugins import HeatMap, HeatMapWithTime
from pandas import DataFrame

from shared import Table, Cache, NavBar, FileSelection, Filter, ColumnType, FillColumnSelection

# Fine, Shiny
import branca, certifi, xyzservices


def server(input: Inputs, output: Outputs, session: Session):

	Info = {
		"example1.txt": "This example dataset shows deaths from a cholera outbreak in 1854. John Snow used this data in conjunction with local pump locations as evidence that cholera is spread by contaminated water. A digitised version of the data is available online, courtesy of Robin Wilson (robin@rtwilson.com).",
		"example2.txt": "This example data set shows bike thefts in Vancouver in 2011. The data was obtained from a 2013 Vancouver Sun blog post by Chad Skelton.",
		"example3.txt": "This example data set shows the location of traffic signals in Toronto. The data was obtained from Toronto Open Data. The idea to use this data set comes from this R-bloggers post by Myles Harrison.",
		"example1.csv": "Random data",
		"example21.csv": "A parsed version of the Northeast and North Central Pacific hurricane database (HURDAT2) 2000-2022, available at https://www.nhc.noaa.gov/data/",
		"example3.csv": "Recorded mean temperature (F) in the USA in 2023 as measured by the EPA, available at https://aqs.epa.gov/aqsweb/airdata/FileFormats.html#_daily_summary_files",

	}

	DataCache = Cache("geocoordinate")


	def GenerateMap(df, map, v_col, lon_col, lat_col):
		"""
		@brief Generates a standard heatmap
		@param df The DataFrame containing the data
		@param map The folium map to attach the heatmap to.
		"""

		# Get each column
		longitudes = df[lon_col].tolist()
		latitudes = df[lat_col].tolist()
		values = [1] * len(latitudes) if input.Uniform() else df[v_col].tolist()

		HeatMap(list(zip(latitudes, longitudes, values)),
		min_opacity=input.Opacity(),
		radius=input.Radius(),
		blur=input.Blur()).add_to(map)
		map.fit_bounds(map.get_bounds())


	def GenerateTemporalMap(df, map, t_col, v_col, lon_col, lat_col):
		"""
		@brief Generates a temporal heatmap
		@param df The DataFrame containing the data
		@param map The folium map to attach the heatmap to.
		"""

		# Sort by time so we can work linearly.
		df = df.sort_values(by=t_col)

		# Normalize
		if not input.Uniform():
			values = df[v_col]
			df[v_col] = (values - values.min()) / (values.max() - values.min())

		# Group data by time
		data = []
		for time, group_df in df.groupby(t_col):
			time_slice = []
			for _, row in group_df.iterrows():
				lat = row[lat_col]
				lon = row[lon_col]
				value = 1 if input.Uniform() else row[v_col]
				time_slice.append([lat, lon, value])
			data.append(time_slice)

		# Make the heamap
		HeatMapWithTime(
			data,
			index=df[t_col].drop_duplicates().to_list(),
			radius=input.Radius(),
			min_opacity=input.Opacity(),
			blur=input.Blur(),
			max_speed=60).add_to(map)


	async def LoadMap():
		"""
		@brief Generates a map with the provided information
		@returns the Folium.Map
		"""

		df = await DataCache.Load(input)

		# Give a placeholder map if nothing is selected, which should never really be the case.
		if df.empty: return FoliumMap((53.5213, -113.5213), tiles=input.MapType(), zoom_start=15)

		if not input.Uniform():
			v_col = input.ValueColumn()
			if v_col not in df: return
		else:
			v_col = None

		lon_col = Filter(df.columns, ColumnType.Longitude, only_one=True)
		lat_col = Filter(df.columns, ColumnType.Latitude, only_one=True)

		map = FoliumMap((df[lat_col][0], df[lon_col][0]), tiles=input.MapType())

		if type(input.ROI()) is tuple and v_col is not None:
			m, M = df[v_col].min(), df[v_col].max()
			l, u = input.ROI()[0], input.ROI()[1]
			if l >= m and u <= M:
				oob = []
				for index, row in df.iterrows():
					v = row[v_col]
					if v < l or v > u: oob.append(index)
				df = df.drop(oob)

		# Generate the right heatmap.
		if input.Temporal():
			t_col = input.TimeColumn()
			if t_col not in df: return
			GenerateTemporalMap(df, map, t_col, v_col, lon_col, lat_col)
		else:
			GenerateMap(df, map, v_col, lon_col, lat_col)
		return map


	@output
	@render.data_frame
	@reactive.event(input.Update, input.Reset, input.Example, input.File, ignore_none=False, ignore_init=False)
	async def LoadedTable(): return await DataCache.Load(input)


	@output
	@render.ui
	@reactive.event(input.Update, input.Reset, input.Example, input.File, input.TimeColumn, input.ValueColumn, input.Temporal, input.MapType, input.Opacity, input.Radius, input.Blur, input.Uniform, input.ROI, ignore_none=False, ignore_init=False)
	async def Map(): return await LoadMap()


	@output
	@render.text
	def ExampleInfo(): return Info[input.Example()]


	@render.download(filename="table.csv")
	async def DownloadTable(): df = await DataCache.Load(input); yield df.to_string()


	@render.download(filename="heatmap.html")
	async def DownloadHeatmap(): m = await LoadMap(); yield m.get_root().render()


	@reactive.Effect
	@reactive.event(input.Update)
	async def Update(): await DataCache.Update(input)


	@reactive.Effect
	@reactive.event(input.Reset)
	async def Reset(): await DataCache.Purge(input)


	@reactive.Effect
	@reactive.event(input.TableRow, input.TableCol, input.Example, input.File, input.Reset, input.Update)
	async def UpdateTableValue():
		"""
		@brief Updates the label for the Value input to display the current value.
		"""
		df = await DataCache.Load(input)

		rows, columns = df.shape
		row, column = int(input.TableRow()), int(input.TableCol())

		if 0 <= row <= rows and 0 <= column <= columns:
			ui.update_text(id="TableVal", label="Value (" + str(df.iloc[row, column]) + ")"),


	@reactive.Effect
	@reactive.event(input.Example, input.File, input.Reset, input.Update, input.Temporal, input.Uniform)
	async def UpdateColumns():
		df = await DataCache.Load(input)
		if not input.Uniform():
			if not FillColumnSelection(df.columns, ColumnType.Value, "ValueColumn"):
				ui.update_checkbox(id="Uniform", value=True)
		if input.Temporal(): FillColumnSelection(df.columns, ColumnType.Time, "TimeColumn")


	@reactive.Effect
	@reactive.event(input.Temporal)
	def UpdateSliders():
		"""
		@brief Heatmapper and HeatmapperWithTime have different bounds for these rules, so we
		update on the fly.
		"""
		if input.Temporal():
			ui.update_slider(id="Opacity", value=0.6, min=0.0, max=1.0, step=0.1)
			ui.update_slider(id="Radius", value=15, min=1, max=50, step=1)
			ui.update_slider(id="Blur", value=0.8, min=0.0, max=1.0, step=0.1)
		else:
			ui.update_slider(id="Opacity", value=0.5, min=0.0, max=1.0, step=0.1),
			ui.update_slider(id="Radius", value=25, min=5, max=50, step=5),
			ui.update_slider(id="Blur", value=15, min=1, max=30, step=1),


	@reactive.Effect
	@reactive.event(input.Update, input.Reset, input.Example, input.File, input.ValueColumn,ignore_none=False, ignore_init=False)
	async def UpdateROI():
		df = await DataCache.Load(input)
		v_col = input.ValueColumn()

		if v_col not in df: ui.update_slider(id="ROI", value=0, min=0, max=0)
		else:
			m, M = int(df[v_col].min()), int(df[v_col].max())
			ui.update_slider(id="ROI", value=(m, M), min=m, max=M)


app_ui = ui.page_fluid(

	NavBar("Geocoordinate"),

	ui.layout_sidebar(
		ui.sidebar(

			FileSelection(
				examples={
					"example1.txt": "Example 1",
					"example2.txt": "Example 2",
					"example3.txt": "Example 3",
					"example1.csv": "Example 4",
					"example21.csv": "Example 5",
					"example3.csv": "Example 6"
				},
				types=[".csv", ".txt", ".xlsx"]
			),

			ui.input_checkbox(id="Temporal", label="Temporal Heatmap"),
			ui.input_checkbox(id="Uniform", label="Uniform Values"),

			# Only provide a temporal column if we're working with time
			ui.panel_conditional(
				"input.Temporal",
				ui.input_select(id="TimeColumn", label="Time Column", choices=[], multiple=False),
			),

			# Only give an option if we aren't working with a uniform value range.
			ui.panel_conditional(
				"!input.Uniform",
				ui.input_select(id="ValueColumn", label="Value Column", choices=[], multiple=False),
			),

			# Only OpenStreatMap and CartoDB Positron seem to work.
			ui.input_radio_buttons(id="MapType", label="Map Type", choices=["OpenStreetMap", "CartoDB Positron"], selected="CartoDB Positron"),

			ui.input_slider(id="Opacity", label="Heatmap Opacity", value=0.5, min=0.0, max=1.0, step=0.1),
			ui.input_slider(id="Radius", label="Size of Points", value=25, min=5, max=50, step=5),
			ui.input_slider(id="Blur", label="Blurring", value=15, min=1, max=30, step=1),

			ui.input_slider(id="ROI", label="Range of Interest", value=(0,0), min=0, max=100),

			# Add the download buttons.
			ui.download_button("DownloadHeatmap", "Heatmap"),
			ui.download_button("DownloadTable", "Table"),
		),

		# Add the main interface tabs.
		ui.navset_tab(
				ui.nav_panel("Interactive", ui.output_ui("Map")),
				Table,
		),
	)
)

app = App(app_ui, server)