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

from shared import Cache, NavBar, MainTab, FileSelection, Pyodide, Filter, ColumnType, TableOptions, Raw
from geojson import Mappings

# Fine, Shiny
import branca, certifi, xyzservices

URL = f"{Raw}/geomap/data/" if Pyodide else "../data/"

def server(input, output, session):

	Info = {
		"example1.txt": "This example file is from the Open Data Portal. The data is from a carbon monoxide emissions study conducted by Environment Canada. The three columns represent results from 1990, 2000, and 2013.",
		"example2.txt": "This example file is from Statistics Canada. The data is adapted from New cases and age-standardized rate for primary cancer (based on the February 2014 CCR tabulation file), by cancer type and sex, Canada, provinces and territories. The columns represent new cancer cases (age-standardized rate per 100,000 population) from 2006 to 2010.",
		"example3.txt": "This example file is from the U.S. Centers for Disease Control and Prevention. The data is from Diagnosed Diabetes, Age Adjusted Rate (per 100), Adults - Total, 2013.",
		"example6.csv": "COVID 19 information reported by the Canadian Government, available at https://open.canada.ca/data/en/dataset/261c32ab-4cfd-4f81-9dea-7b64065690dc/resource/39434379-45a1-43d5-aea7-a7a50113c291"
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

	@reactive.effect
	@reactive.event(input.SourceFile, input.File, input.Example, input.Reset)
	async def UpdateData(): Data.set((await DataCache.Load(input))); Valid.set(False)


	@reactive.effect
	@reactive.event(input.JSONUpload, input.JSONSelection, input.JSONFile)
	async def UpdateGeoJSON(): JSON.set(await DataCache.Load(
			input,
			source_file=input.JSONUpload(),
			example_file=input.JSONSelection(),
			source=URL,
			input_switch=input.JSONFile(),
			default=None
		))


	def GetData(): return Table.data_view() if Valid() else Data()


	def LoadChoropleth(df, map, geojson, k_col, v_col):
		"""
		@brief Applies a Choropleth to a Folium Map
		@param df: The DataFrame that contains information to plot
		@param map: The Folium map
		@param geojson: The geojson information that contains territory information
		@param k_vol: The name of the column within df that contains names
		@param v_col: the name of the column within df that contains the values to plot.
		"""
		Choropleth(
				geo_data=geojson,
				name="choropleth",
				data=df,
				columns=[k_col, v_col],
				key_on="feature.properties.name",
				fill_color=input.ColorMap().lower(),
				fill_opacity=input.Opacity(),
				line_opacity=input.Opacity(),
				legend_name="Legend",
				bins=input.Bins()
		).add_to(map)


	def LoadTemporalChoropleth(df, map, geojson, k_col, v_col):
		"""
		@brief Applies a TimeSliderChoropleth to a Folium map
		@param df: The DataFrame that contains data to plot
		@param map: The Folium map
		@param geojson: The geojson information that contains territory information
		@param k_col: The name of the column within df that contains names
		@param v_col: The name of the column within df that contains the values to plot.
		@info df can either contain a Time column, or all non key-columns will be handled as time columns.
		"""
		# Check if we have a dedicated time column, or separate columns for each time slot.
		column = Filter(df.columns, ColumnType.Time, bad = [k_col, v_col], only_one=True, reject_unknown=True)
		values = v_col if column else df.columns[1:]

		color = input.ColorMap()
		if color == "Inferno": cmap = linear.inferno.scale
		elif color == "Magma": cmap = linear.magma.scale
		elif color == "Plasma": cmap = linear.plasma.scale
		elif color == "Viridis": cmap = linear.viridis.scale

		m, M = df[values].values.min(), df[values].values.max()
		colormap = cmap(m, M)

		style = {}

		if column:
			grouped = df.groupby(k_col)

			for i, (name, group) in enumerate(grouped):
				style[i] = {}
				for a, row in group.iterrows():

					# If the year isn't parsable, pray that it just works without parsing :)
					try: year = str(to_datetime(row[column].split()[0]).timestamp()).split('.')[0]
					except Exception: year = row[column]

					style[i][year] = {'color': colormap(row[v_col]), 'opacity': input.Opacity()}
				for feature in geojson["features"]:
							if feature["properties"]["name"] == name:
									feature["id"] = i

		else:
			# Convert DataFrame to style format
			for i, row in df.iterrows():
					name = row[k_col]
					style[i] = {}
					for year in df.columns[1:]:
						value = row[year]

						# If the year isn't parsable, pray that it just works without parsing :)
						try: year = str(to_datetime(year.split()[0]).timestamp()).split('.')[0]
						except Exception: pass

						style[i][year] = {'color': colormap(value), 'opacity': input.Opacity()}
					for feature in geojson["features"]:
						if feature["properties"]["name"] == name:
								feature["id"] = i

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
		if input.Type() == "Integer": value = int(patch["value"])
		elif input.Type() == "Float": value = float(patch["value"])
		else: value = patch["value"]
		return value


	@output
	@render.ui
	def Heatmap():
		with ui.Progress() as p:

			p.inc(message="Loading input...")
			df = GetData()
			if df is None: return

			p.inc(message="Loading GeoJSON...")
			geojson = JSON()

			p.inc(message="Formatting...")
			k_col, v_col = input.KeyColumn(), input.ValueColumn()
			if k_col not in df or v_col not in df: return

			# Give a placeholder map if nothing is selected, which should never really be the case.
			if df.empty or geojson is None: return FoliumMap((53.5213, -113.5213), tiles=input.MapType(), zoom_start=15)

			# Create map
			map = FoliumMap(tiles=input.MapType())

			if type(input.ROI()) is tuple:
				m, M = df[v_col].min(), df[v_col].max()
				l, u = input.ROI()[0], input.ROI()[1]
				if l >= m and u <= M:
					oob = []
					for index, row in df.iterrows():
						v = row[v_col]
						if v < l or v > u: oob.append(index)
					df = df.drop(oob)

			# Load the choropleth.
			p.inc(message="Plotting...")
			if input.Temporal(): LoadTemporalChoropleth(df, map, geojson, k_col, v_col)
			else: LoadChoropleth(df, map, geojson, k_col, v_col)

			map.fit_bounds(map.get_bounds())
			return map

	@output
	@render.data_frame
	def GeoJSON():
		geojson = JSON()
		if geojson is None: return
		names = [feature['properties']['name'] for feature in geojson['features']]
		return DataFrame({'Name': names})


	@output
	@render.text
	@reactive.event(input.SourceFile, input.Example)
	def ExampleInfo(): return Info[input.Example()]


	@render.download(filename="table.csv")
	def DownloadTable(): yield GetData().to_string()


	@render.download(filename="heatmap.html")
	def DownloadHeatmap(): yield LoadMap().get_root().render()


	@reactive.Effect
	def UpdateColumns():
		df = GetData()
		key = Filter(df.columns, ColumnType.Name, only_one=True, ui_element="KeyColumn")
		Filter(df.columns, ColumnType.Value, bad=[key], ui_element="ValueColumn")


	@reactive.Effect
	def UpdateROI():
		df = GetData()
		v_col = input.ValueColumn()
		if v_col not in df: ui.update_slider(id="ROI", value=0, min=0, max=0)
		else:
			m, M = int(df[v_col].min()), int(df[v_col].max())
			ui.update_slider(id="ROI", value=(m, M), min=m, max=M)


app_ui = ui.page_fluid(

	NavBar(),

	ui.layout_sidebar(
		ui.sidebar(

			FileSelection(
				examples={"example1.txt": "Example 1", "example2.txt": "Example 2", "example3.txt": "Example 3", "example6.csv": "Example 4"},
				types=[".csv", ".txt", ".dat", ".tsv", ".tab", ".xlsx", ".xls", ".odf"],
				project="Geomap"
			),

			TableOptions(),

			ui.input_radio_buttons(id="JSONFile", label="Specify a GeoJSON File", choices=["Provided", "Upload"], selected="Provided", inline=True),
			ui.panel_conditional(
				"input.JSONFile === 'Upload'",
				ui.input_file("JSONUpload", "Choose a File", accept=[".geojson"], multiple=False),
			),
			ui.panel_conditional(
				"input.JSONFile === 'Provided'",
				ui.input_select(id="JSONSelection", label=None, choices=Mappings, multiple=False, selected="canada.geojson"),
			),

			ui.input_checkbox(id="Temporal", label="Temporal Choropleth"),

			ui.input_select(id="KeyColumn", label="Key", choices=[], multiple=False),
			ui.input_select(id="ValueColumn", label="Value", choices=[], multiple=False),

			# Only OpenStreatMap and CartoDB Positron seem to work.
			ui.input_radio_buttons(id="MapType", label="Map Type", choices=["OpenStreetMap", "CartoDB Positron"], selected="CartoDB Positron"),

			ui.input_select(id="ColorMap", label="Color Map", choices=["Inferno", "Magma", "Plasma", "Viridis"], selected="Viridis"),
			ui.input_slider(id="Opacity", label="Heatmap Opacity", value=0.5, min=0.0, max=1.0, step=0.1),
			ui.input_slider(id="Bins", label="Number of Colors", value=8, min=3, max=8, step=1),
			ui.input_slider(id="ROI", label="Range of Interest", value=(0,0), min=0, max=100),

			# Add the download buttons.
			"Download",
			ui.download_button("DownloadHeatmap", "Heatmap"),
			ui.download_button("DownloadTable", "Table"),

			id="SidebarPanel",
		),

		MainTab(ui.nav_panel("GeoJSON", ui.output_data_frame("GeoJSON")), m_type=ui.output_ui),
	)
)

app = App(app_ui, server)
