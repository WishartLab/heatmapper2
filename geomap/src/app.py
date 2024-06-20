#
# Heatmapper
# Geomap
#
# This file contains the ShinyLive application for Geomap Heatmapper.
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
from folium import Map as FoliumMap, Choropleth
from folium.plugins import TimeSliderChoropleth
from pandas import DataFrame, to_datetime
from branca.colormap import linear
from pathlib import Path
from json import loads
from datetime import datetime
from time import mktime

from shared import Cache, NavBar, MainTab, FileSelection, Pyodide, Filter, ColumnType, TableOptions, Raw, InitializeConfig, ColorMaps, Error, Update, Msg
from geojson import Mappings

try:
	from user import config
except ImportError:
	from config import config

# Fine, Shiny
import branca, certifi, xyzservices

URL = f"{Raw}/geomap/data/" if Pyodide else "../data/"

def server(input, output, session):

	Info = {
		"example1.txt": "This example file is from the Open Data Portal. The data is from a carbon monoxide emissions study conducted by Environment Canada. The three columns represent results from 1990, 2000, and 2013.",
		"example2.txt": "This example file is from Statistics Canada. The data is adapted from New cases and age-standardized rate for primary cancer (based on the February 2014 CCR tabulation file), by cancer type and sex, Canada, provinces and territories. The columns represent new cancer cases (age-standardized rate per 100,000 population) from 2006 to 2010.",
		"example3.txt": "This example file is from the U.S. Centers for Disease Control and Prevention. The data is from Diagnosed Diabetes, Age Adjusted Rate (per 100), Adults - Total, 2013.",
		"example6.csv": "COVID 19 information reported by the Canadian Government, available at https://open.canada.ca/data/en/dataset/261c32ab-4cfd-4f81-9dea-7b64065690dc/resource/39434379-45a1-43d5-aea7-a7a50113c291",
		"https://media.githubusercontent.com/media/WishartLab/heatmapper2/main/geomap/example_input/owid-covid-data.csv": "Global COVID 19 Statistics from Our World in Data."
	}

	def HandleData(path):
		"""
		@brief A custom Data Handler for the Cache.
		@param path: Path to the file
		@returns A data object from the cache.
		@info This Data Handler supports geojson files as json
		"""
		if path.suffix == ".geojson": return loads(path.open().read())
		else: return DataCache.DefaultHandler(path)
	DataCache = Cache("geomap", DataHandler=HandleData)
	Data = reactive.value(None)
	Valid = reactive.value(False)
	JSON = reactive.value(None)

	InitializeConfig(config, input)


	@reactive.effect
	@reactive.event(input.SourceFile, input.File, input.Example, input.Reset)
	async def UpdateData():
		Data.set((await DataCache.Load(input, p=ui.Progress())));
		Valid.set(False)

		columns = Data().columns
		key = Filter(columns, ColumnType.Name, id="KeyColumn")
		val = Filter(columns, ColumnType.Value, id="ValueColumn", all=True)
		if val:
			choice = 0
			while choice < len(val) and val[choice] == key: choice += 1
			ui.update_select(id="ValueColumn", selected=val[choice])


	@reactive.effect
	@reactive.event(input.JSONUpload, input.JSONSelection, input.JSONFile)
	async def UpdateGeoJSON():
		JSON.set(await DataCache.Load(
			input,
			source_file=input.JSONUpload(),
			example_file=input.JSONSelection(),
			source=URL,
			input_switch=input.JSONFile(),
			example="Provided",
			default=None,
			p=ui.Progress(),
			p_name="GeoJSON"
		))

		geojson = JSON()
		if geojson is None: return
		properties = list(geojson['features'][0]['properties'].keys())
		Filter(properties, ColumnType.NameGeoJSON, id="KeyProperty")


	def GetData(): return Table.data_view() if Valid() else Data()


	def LoadChoropleth(df, map, geojson, k_col, v_col, k_prop, p):
		"""
		@brief Applies a Choropleth to a Folium Map
		@param df: The DataFrame that contains information to plot
		@param map: The Folium map
		@param geojson: The geojson information that contains territory information
		@param k_vol: The name of the column within df that contains names
		@param v_col: the name of the column within df that contains the values to plot.
		"""

		colormap = config.ColorMap().lower()
		opacity = config.Opacity()
		bins = config.Bins()

		Choropleth(
				geo_data=geojson,
				name="choropleth",
				data=df,
				columns=[k_col, v_col],
				key_on=f"feature.properties.{k_prop}",
				fill_color=colormap,
				fill_opacity=opacity,
				line_opacity=opacity,
				legend_name="Legend",
				bins=bins
		).add_to(map)


	def LoadTemporalChoropleth(df, map, geojson, k_col, v_col, k_prop, p):
		"""
		@brief Applies a TimeSliderChoropleth to a Folium map
		@param df: The DataFrame that contains data to plot
		@param map: The Folium map
		@param geojson: The geojson information that contains territory information
		@param k_col: The name of the column within df that contains names
		@param v_col: The name of the column within df that contains the values to plot.
		@info df can either contain a Time column, or all non key-columns will be handled as time columns.
		"""

		def Timestamp(time):
			year, month, day = 1970, 1, 1
			try:
				if len(time) >= 1: year = int(time[0])
				if len(time) >= 2: month = int(time[1])
				if len(time) >= 3: day = int(time[2])
			except ValueError: pass
			return str(round(mktime(datetime(year, month, day).timetuple())))


		# Check if we have a dedicated time column, or separate columns for each time slot.
		column = Filter(df.columns, ColumnType.Time)

		color = config.ColorMap()
		if color == "Inferno": cmap = linear.inferno.scale
		elif color == "Magma": cmap = linear.magma.scale
		elif color == "Plasma": cmap = linear.plasma.scale
		elif color == "Viridis": cmap = linear.viridis.scale
		elif color == "Cividis": cmap = linear.cividis.scale

		m, M = df[v_col].min(), df[v_col].max()

		colormap = cmap(m, M)
		style = {}

		if column:
			p.inc(message="Grouping by Date...")
			grouped = df.groupby(k_col)

			p.inc(message="Formatting...")
			for i, (name, group) in enumerate(grouped):
				style[i] = {}
				for time, value in zip(group[column], group[v_col]):
					timestamp = Timestamp(time.split("-"))
					if timestamp is None: return

					style[i][timestamp] = {'color': colormap(value), 'opacity': config.Opacity()}
				for feature in geojson["features"]:
							if feature["properties"][k_prop] == name:
									feature["id"] = i
		else:
			p.inc(message="Formatting...")
			time_columns = df.columns.drop(k_col)
			for i, row in df.iterrows():
					name = row[k_col]
					style[i] = {}
					for time in time_columns:
						value = row[time]

						time = time.split(" ")
						timestamp = Timestamp(time if len(time) == 1 else [time[0]])
						if timestamp is None: return

						style[i][timestamp] = {'color': colormap(value), 'opacity': config.Opacity()}
					for feature in geojson["features"]:
						if feature["properties"][k_prop] == name:
								feature["id"] = i

		p.inc(message="Creating Choropleth...")
		TimeSliderChoropleth(
				data=geojson,
				styledict=style,
		).add_to(map)
		colormap.add_to(map)


	@output
	@render.data_frame
	def Table(): Valid.set(True); return render.DataGrid(Data(), editable=True)


	@Table.set_patch_fn
	def UpdateTable(*, patch: render.CellPatch) -> render.CellValue:
		if config.Type() == "Integer": value = int(patch["value"])
		elif config.Type() == "Float": value = float(patch["value"])
		else: value = patch["value"]
		return value


	def GenerateHeatmap():
		with ui.Progress() as p:

			p.inc(message="Loading input...")
			df = GetData()
			if df is None: return

			p.inc(message="Loading GeoJSON...")
			try:
				geojson = JSON()
				properties = list(geojson['features'][0]['properties'].keys())
			except Exception:
				return

			p.inc(message="Formatting...")
			k_col, v_col, k_prop = config.KeyColumn(), config.ValueColumn(), config.KeyProperty()
			if k_col not in df or v_col not in df or k_prop not in properties: return

			map_type = config.MapType()

			# Give a placeholder map if nothing is selected, which should never really be the case.
			if df.empty or geojson is None: return FoliumMap((53.5213, -113.5213), tiles=map_type, zoom_start=15)

			# Create map
			map = FoliumMap(tiles=map_type)

			p.inc(message="Dropping Invalid Values...")
			names = []
			for feature in geojson["features"]:
				names.append(feature["properties"][k_prop])

			to_drop = []
			l, u = config.Min(), config.Max()
			for index, key, value in zip(df.index, df[k_col], df[v_col]):
				if key not in names: to_drop.append(index)
				elif config.ROI() and (value < l or value > u):
					if config.ROI_Mode() == "Remove": to_drop.append(index)
					elif config.ROI_Mode() == "Round": df.at[index, v_col] = u if value > u else l
			df = df.drop(to_drop)
			if len(df) == 0:
				Error("No locations! Ensure Key Column and Key Properties are correct, and your ROI is properly set!")
				return

			# Load the choropleth.
			p.inc(message="Plotting...")
			if input.Temporal(): LoadTemporalChoropleth(df, map, geojson, k_col, v_col, k_prop, p)
			else: LoadChoropleth(df, map, geojson, k_col, v_col, k_prop, p)

			map.fit_bounds(map.get_bounds())
			return map

	@output
	@render.ui
	def Heatmap(): return GenerateHeatmap()

	@output
	@render.ui
	@reactive.event(input.Update)
	def HeatmapReactive(): return GenerateHeatmap()


	@output
	@render.data_frame
	def GeoJSON():
		try:
			geojson = JSON()
			names = [feature['properties'][config.KeyProperty()] for feature in geojson['features']]
			return DataFrame({config.KeyProperty(): names})
		except Exception:
			Error("Could not render the GeoJSON table!")


	@reactive.effect
	@reactive.event(input.ExampleInfo)
	def ExampleInfo():
		Msg(ui.HTML(Info[input.Example()]))


	@render.download(filename="table.csv")
	def DownloadTable(): yield GetData().to_string()


	@render.download(filename="heatmap.html")
	def DownloadHeatmap(): yield LoadMap().get_root().render()


app_ui = ui.page_fluid(

	NavBar(),

	ui.layout_sidebar(
		ui.sidebar(

			FileSelection(
				examples={
				"example1.txt": "Example 1",
				"example2.txt": "Example 2",
				"example3.txt": "Example 3",
				"example6.csv": "Example 4",
				"https://heatmapper2.ca/geomap/example_input/owid-covid-data.csv": "Example 5"},
				types=[".csv", ".txt", ".dat", ".tsv", ".tab", ".xlsx", ".xls", ".odf"],
				project="Geomap"
			),

			ui.input_radio_buttons(id="JSONFile", label="Specify a GeoJSON File", choices=["Provided", "Upload"], selected="Provided", inline=True),
			ui.panel_conditional(
				"input.JSONFile === 'Upload'",
				ui.input_file("JSONUpload", None, accept=[".geojson"], multiple=False),
			),
			ui.panel_conditional(
				"input.JSONFile === 'Provided'",
				ui.input_select(id="JSONSelection", label=None, choices=Mappings, multiple=False, selected="canada.geojson"),
			),


			TableOptions(config),

			ui.panel_conditional(
				"input.MainTab === 'HeatmapTab'",

				Update(),

				ui.HTML("<b>Columns/Properties</b>"),
				config.KeyColumn.UI(ui.input_select, id="KeyColumn", label="Key", choices=[]),
				config.ValueColumn.UI(ui.input_select, id="ValueColumn", label="Value", choices=[]),
				config.KeyProperty.UI(ui.input_select, id="KeyProperty", label="GeoJSON", choices=[]),

				ui.HTML("<b>Heatmap</b>"),
				config.Temporal.UI(ui.input_checkbox, id="Temporal", label="Temporal"),
				config.MapType.UI(ui.input_select, id="MapType", label="Map", choices={"CartoDB Positron": "CartoDB", "OpenStreetMap": "OSM"}),
				config.Opacity.UI(ui.input_numeric, id="Opacity", label="Opacity", min=0.0, max=1.0, step=0.1),

				ui.HTML("<b>Colors</b>"),
				config.ColorMap.UI(ui.input_select, id="ColorMap", label="Map", choices=ColorMaps),
				config.Bins.UI(ui.input_numeric, id="Bins", label="Number", min=3, max=253, step=1),

				ui.HTML("<b>Range of Interest</b>"),
				config.ROI.UI(ui.input_checkbox, make_inline=False, id="ROI", label="Enable (Lower/Upper)"),
				config.ROI_Mode.UI(ui.input_radio_buttons, make_inline=False, id="ROI_Mode", label=None, choices=["Remove", "Round"], inline=True),
				ui.layout_columns(
					config.Min.UI(ui.input_numeric,make_inline=False, id="Min", label=None, min=0),
					config.Max.UI(ui.input_numeric, make_inline=False, id="Max", label=None, min=0),
				),

				ui.download_button(id="DownloadHeatmap", label="Download"),
			),
			padding="10px",
			gap="20px",
			width="250px",
		),

		MainTab(ui.nav_panel("GeoJSON", ui.output_data_frame("GeoJSON")), m_type=ui.output_ui),
		height="90vh",
	)
)

app = App(app_ui, server)
