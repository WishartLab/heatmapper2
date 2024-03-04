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

from shared import Table, Cache, NavBar, FileSelection

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


	def GenerateMap(df, map):
		"""
		@brief Generates a standard heatmap
		@param df The DataFrame containing the data
		@param map The folium map to attach the heatmap to.
		"""

		default_value = input.ValueColumn()
		if default_value is None or default_value not in df: return

		# Get the long and lat.
		longitudes = df["Longitude"].tolist()
		latitudes = df["Latitude"].tolist()

		values = df[default_value].tolist()

		HeatMap(list(zip(latitudes, longitudes, values)),
		min_opacity=input.Opacity(),
		radius=input.Radius(),
		blur=input.Blur()).add_to(map)
		map.fit_bounds(map.get_bounds())


	def GenerateTemporalMap(df, map):
		"""
		@brief Generates a temporal heatmap
		@param df The DataFrame containing the data
		@param map The folium map to attach the heatmap to.
		"""

		default_value = input.ValueColumn()
		if default_value is None or default_value not in df: return
		default_time = input.TimeColumn()
		if default_time is None or default_time not in df: return

		# Sort by time so we can work linearly.
		df = df.sort_values(by=default_time)

		# Normalize
		values = df[default_value]
		df[default_value] = (values - values.min()) / (values.max() - values.min())

		# Group data by time
		data = []
		for time, group_df in df.groupby(default_time):
			time_slice = []
			for _, row in group_df.iterrows():
				lat = row["Latitude"]
				lon = row["Longitude"]
				value = row[default_value]
				time_slice.append([lat, lon, value])
			data.append(time_slice)

		# Make the heamap
		HeatMapWithTime(
			data,
			index=df[default_time].drop_duplicates().to_list(),
			radius=input.Radius(),
			min_opacity=input.Opacity(),
			blur=input.Blur(),
			max_speed=60,).add_to(map)


	async def LoadMap():
		"""
		@brief Generates a map with the provided information
		@returns the Folium.Map
		"""

		df = await DataCache.Load(input)

		# Give a placeholder map if nothing is selected, which should never really be the case.
		if df.empty: return FoliumMap((53.5213, -113.5213), tiles=input.MapType(), zoom_start=15)


		map = FoliumMap((df["Latitude"][0], df["Longitude"][0]), tiles=input.MapType())

		# Generate the right heatmap.
		if input.Temporal(): GenerateTemporalMap(df, map)
		else: GenerateMap(df, map)
		return map


	@output
	@render.data_frame
	@reactive.event(input.Update, input.Reset, input.Example, input.File, ignore_none=False, ignore_init=False)
	async def LoadedTable(): return await DataCache.Load(input)


	@output
	@render.ui
	@reactive.event(input.Update, input.Reset, input.Example, input.File, input.TimeColumn, input.ValueColumn, input.Temporal, input.MapType, input.Opacity, input.Radius, input.Blur, ignore_none=False, ignore_init=False)
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
	@reactive.event(input.Example, input.File, input.Reset, input.Update)
	async def UpdateColumns():

		# Give options for the key and value columns
		df = await DataCache.Load(input)
		choices = df.columns.tolist()
		if choices:

			default_time = None
			for time in ["Time", "Date"]:
				if time in df: default_time = time; break
			if not default_time:
				columns = df.columns.tolist()
				for n in ["Longitude", "Latitude", "Weight", "Intensity", "Value"]:
					if n in columns:
						columns.remove(n)
				if columns:
					default_time = columns[0]
				else:
					default_time = df.columns[0]

			default_value = None
			for value in ["Weight", "Intensity", "Value"]:
				if value in df: default_value = value; break
			if not default_value:
				columns = df.columns.tolist()
				for n in ["Longitude", "Latitude", "Time", "Date"]:
					if n in columns:
						columns.remove(n)
				if columns:
					default_value = columns[0]
				else:
					default_value = df.columns[0]

			ui.update_select(id="TimeColumn", choices=choices, selected=default_time)
			ui.update_select(id="ValueColumn", choices=choices, selected=default_value)


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

			ui.input_checkbox(id="Temporal", label="Temporal Choropleth"),

			# Only provide a temporal column if we're working with time
			ui.panel_conditional(
				"input.Temporal",
				ui.input_select(id="TimeColumn", label="Time Column", choices=[], multiple=False),
			),

			ui.input_select(id="ValueColumn", label="Value Column", choices=[], multiple=False),

			# Only OpenStreatMap and CartoDB Positron seem to work.
			ui.input_radio_buttons(id="MapType", label="Map Type", choices=["OpenStreetMap", "CartoDB Positron"], selected="CartoDB Positron"),

			ui.input_slider(id="Opacity", label="Heatmap Opacity", value=0.5, min=0.0, max=1.0, step=0.1),
			ui.input_slider(id="Radius", label="Size of Points", value=25, min=5, max=50, step=5),
			ui.input_slider(id="Blur", label="Blurring", value=15, min=1, max=30, step=1),

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