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

from shared import Cache, NavBar, MainTab, FileSelection, Pyodide, Filter, ColumnType, TableOptions, Raw, InitializeConfig, ColorMaps, Error, Update, Msg, File
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
		"example1.txt": '<u>Input type:</u> .txt Data <br><u>Contents:</u> Data from a carbon monoxide emissions study conducted by Environment Canada. The three columns show carbon monoxide emissions (in kilotonnes) from 1990, 2005, and 2022. The data represents emissions from human activities only (e.g. emissions from fires and vegetation are excluded). Data is excluded for Nunavut in 1990, as it was a part of the Northwest Territories until 1999. <br><u>Source:</u> <a href="https://open.canada.ca/data/en/dataset/8c3d9333-fd03-4f32-b991-580e4fa6c343/resource/d1c16e1c-abd4-47f7-abd1-218a9e4848d5"; target="_blank">Open Data Portal</a>',
		"example2.txt": '''<u>Input type:</u> .txt Data <br><u>Contents:</u> Data from Statistics Canada: 'New cases and age-standardized rates of primary cancer, by cancer type and sex'. The columns represent new cancer cases (age-standardized rate per 100,000 population) each year from 2013 to 2017 in Canadian provinces and territories.<br><u>Source:</u> <a href="https://doi.org/10.25318/1310074701-eng"; target="_blank">Statistics Canada</a>''',
		"example6.csv": '''<u>Input type:</u> .csv Data<br> <u>Contents:</u> COVID-19 information by province or territory, reported by the Canadian Government, from February 8, 2020 to February 17, 2024. The 'Value Column' dropdown indicates which column is visualized. <br><u>Source:</u> <a href="https://open.canada.ca/data/en/dataset/261c32ab-4cfd-4f81-9dea-7b64065690dc/resource/39434379-45a1-43d5-aea7-a7a50113c291"; target="_blank">Open Data Portal</a>''',
		#"https://media.githubusercontent.com/media/WishartLab/heatmapper2/main/geomap/example_input/owid-covid-data.csv": "File type: csv<br>Contents: Global COVID 19 Statistics.<br>Source: Our World in Data"
	}

	def HandleData(path, p=None):
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
		p = ui.Progress()
		try:
			Data.set((await DataCache.Load(input, p=p)))
			Valid.set(False)

			columns = Data().columns
			key = Filter(columns, ColumnType.Name, id="KeyColumn")
			val = Filter(columns, ColumnType.Value, id="ValueColumn", all=True)
			if val:
				choice = 0
				while choice < len(val) and val[choice] == key: choice += 1
				ui.update_select(id="ValueColumn", selected=val[choice])
			DataCache.Invalidate(File(input))
		except:
			p.close()
			Error("File could not be loaded!\nData can be uploaded as a .csv, .tsv, .txt, .xslx, .dat, .tab, or .odf file. GeoJSON files must follow GeoJSON standards.")
			return


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
		elif color == "Plasma": cmap = linear.plasma.scale
		elif color == "Viridis": cmap = linear.viridis.scale
		#elif color == "Cividis": cmap = linear.cividis.scale

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
	def Table(): 
		df = Data()
		print(f"Data: \n{df}")
		if df is None or df.empty:
			return DataFrame({"Note": ["No data to display! Please upload your data or select an example data set in the sidebar."]})
		try:
			grid = render.DataGrid(df, editable=True)
			Valid.set(True)
			return grid
		except:
			return DataFrame({'Note': [ui.HTML('Could not render data! <br><br>Please make sure the names in your dataset match the names in a corresponding GeoJSON. Your input data should also have a column of values, and optionally, time. <br><br>See more formatting info <a href="https://github.com/WishartLab/heatmapper2/wiki/Format#geomap:~:text=work%20with%20it.-,Geomap,-Geomap%20has%20two"; target="_blank"; rel=”noopener noreferrer”>here</a>.')]})


	@Table.set_patch_fn
	def UpdateTable(*, patch: render.CellPatch) -> render.CellValue:
		if config.Type() == "Integer": value = int(patch["value"])
		elif config.Type() == "Float": value = float(patch["value"])
		else: value = patch["value"]
		DataCache.Invalidate(File(input))
		return value
	
	# Info text in welcome tab
	@render.ui
	def Welcome():
		return ui.HTML("""
			<h1>Geomap</h1>
			Geomap displays values based on geographical boundaries, such as country, state, or province. Upload a data file and specify a GeoJSON in the sidebar to get started, or select 'Example' to check out a pre-loaded example. Navigate to the 'Heatmap' tab to see the heatmap, 'Table' to look at the input data, or 'GeoJSON' to see the geographical boundaries available in the currently selected GeoJSON file.
			
			<br><br>
			<img src="https://github.com/WishartLab/heatmapper2/wiki/assets/Geomap.png" alt="Geomap"; style="max-width:500px;">
				 
			<br><br><h3>Format</h3>
			Geomap requires a data file as well as a GeoJSON.
			<br>
			<i>Input data can be formatted as follows:</i>
				 <ul>
				 <li>A 'Name' column and 'Value' column(s), where names in the 'Name' column match available geographical boundaries in the currently selected GeoJSON. Select which value column to display using the 'Value' dropdown, if there is more than one.</li>
				 <li><u>Temporal format 1:</u> A 'Name' column, 'Value' column(s), and a 'Time' column. Rows will be grouped by time and plotted linearly. Names in the 'Name' column should match available geographical boundaries in the currently selected GeoJSON. Select which value column to display using the 'Value' dropdown, if there is more than one. (See Example 3)</li>
				 <li><u>Temporal format 2:</u> A 'Name' column, and multiple 'Time' columns, each containing the value of the associated name at that time (i.e. each row contains a name, and multiple values of that name at different time points). 'Time' column names are parsed such that all characters up to the first whitespace indicate the time (e.g. '1990 [emissions in kilotonnes]' becomes '1990'). See Example 2. </li>
				 </ul>
			<br>
			<i>Geomap heatmaps can be generated from the following file formats:</i>
			<table style="border-spacing: 100px";>
			<tr>
				<th>Table Files</th>
				<th>GeoJSON Files</th>
			</tr>
			<tr>
				<td style="padding-right:50px;">
					<li>.csv</li>
					<li>.dat</li>
					<li>.odf</li>
					<li>.tab</li>
					<li>.tsv</li>
					<li>.txt</li>
					<li>.xls</li>
					<li>.xlsx</li>
				</td>
				<td style="vertical-align:top;">
					<li>standard .geojson files, see <a href="https://geojson.org/"; target="_blank"; rel=”noopener noreferrer”>geojson.org</a></li>
				</td>
			</tr>
			</table>
				 
			<br><h3>Interface</h3>
			Remember to select or upload a GeoJSON with boundaries that match the names in your data.
		
			<br>Click on the '?' icon beside sidebar options to read more about them.
		""")


	def GenerateHeatmap():
		with ui.Progress() as p:

			p.inc(message="Loading input...")
			df = GetData()
			if df is None:
				print("df is None")
				return ui.HTML("No data to display! Please upload your data or select an example data set in the sidebar.")

			p.inc(message="Loading GeoJSON...")
			try:
				geojson = JSON()
				properties = list(geojson['features'][0]['properties'].keys())
			except Exception:
				return ui.HTML('Make sure a GeoJSON is selected in the sidebar, <br>or upload your own following the <a href="https://geojson.org/"; target="_blank"; rel=”noopener noreferrer”>GeoJSON format</a>.')

			p.inc(message="Formatting...")
			k_col, v_col, k_prop = config.KeyColumn(), config.ValueColumn(), config.KeyProperty()
			if k_col not in df or v_col not in df or k_prop not in properties: return ui.HTML("Data could not be displayed. <br>Please upload a Table file and a GeoJSON, or select an example data set in the sidebar. <br><br><i>Uploaded Table files should include: <br>a Key column (e.g. 'name', 'continent', 'country', 'location') <br>and a Value column (e.g. 'value', 'weight', 'intensity')</i>")

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
				Error("No locations found")
				return ui.HTML("No locations were found. <br>Please ensure your data table contains a name column, whose values match a property in the GeoJSON. <br><br>Still having trouble? Make sure your range of interest does not exclude all values!")

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
		geojson = JSON()
		if geojson is None:
			return DataFrame({'Note':[ui.HTML('Hmmm, we could not render the GeoJSON table.<br><br>Make sure a GeoJSON is selected in the sidebar, <br>or upload your own following the <a href="https://geojson.org/"; target="_blank"; rel=”noopener noreferrer”>GeoJSON format</a>.')]})
		try:	
			names = [feature['properties'][config.KeyProperty()] for feature in geojson['features']]
			return DataFrame({config.KeyProperty(): names})
		except:
			Error("Could not render the GeoJSON table!")
			return DataFrame({'Note':[ui.HTML('Hmmm, we could not render the GeoJSON table.<br><br>Make sure a GeoJSON is selected in the sidebar, <br>or upload your own following the <a href="https://geojson.org/"; target="_blank"; rel=”noopener noreferrer”>GeoJSON format</a>.')]})


	@reactive.effect
	@reactive.event(input.ExampleInfo)
	def ExampleInfo():
		Msg(ui.HTML(Info[input.Example()]))
		

	@render.download(filename=lambda: f"table{config.TableType()}")
	def DownloadTable(): 
		data = GetData()
		
		# return error if no data to download
		if data.empty:
			Error("The downloaded table is empty! Please upload your data or select an example data set in the sidebar.")
		
		file_contents = data.to_string()
		yield file_contents


	@render.download(filename="heatmap.html")
	def DownloadHeatmap(): yield GenerateHeatmap().get_root().render()


	@render.download(filename=f"heatmap{config.SettingType()}")
	def DownloadSettings(): 
		'''
		Save current settings to a table file.
		'''
		json = input.JSONSelection()
		if json is None:
			json = input.JSONUpload()

		yield f"Data Filename:\t{File(input)}\nGeoJSON File:\t{json}\nName Column:\t{config.KeyColumn()}\nValue Column:\t{config.ValueColumn()}\nGeoJSON Property:\t{config.KeyProperty()}\nMap Type:\t{config.MapType()}\nData Opacity:\t{config.Opacity()}\nColor Map:\t{config.ColorMap()}\nNumber of Color Bins:\t{config.Bins()}\nRange of Interest:\t{config.ROI()}\nRange of Interest Mode:\t{config.ROI_Mode()}\nRange of Interest Min:\t{config.Min()}\nRange of Interest Max:\t{config.Max()}"


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

	ui.panel_title(title=None, window_title="Geomap"),
	NavBar(),

	ui.layout_sidebar(
		ui.sidebar(

			FileSelection(
				examples={
				"example1.txt": "Ex1: CO Emissions",
				"example2.txt": "Ex2: Cancer Cases",
				"example6.csv": "Ex3: COVID-19"},
				#"https://heatmapper2.ca/geomap/example_input/owid-covid-data.csv": "Example 5"},
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
				config.KeyColumn.UI(ui.input_select, id="KeyColumn", label="Name Column", choices=[], tooltip="Specify a column in your data that contains location names. These location names must correspond to location names in the GeoJSON file. Click on the 'GeoJSON' tab in the main view area to see location names in the currently selected GeoJSON file."),
				config.ValueColumn.UI(ui.input_select, id="ValueColumn", label="Value Column", choices=[], tooltip="Specify a column containing the data to plot. If 'Temporal' is selected, this column is ignored if it does not have a corresponding column with time values."),
				config.KeyProperty.UI(ui.input_select, id="KeyProperty", label="GeoJSON", choices=[], tooltip="Select a property in the GeoJSON file that corresponds to the location names in your data. Click on the 'GeoJSON' tab in the main view area to see available properties in the currently selected GeoJSON file."),

				ui.HTML("<b>Heatmap</b>"),
				config.Temporal.UI(ui.input_checkbox, id="Temporal", label="Temporal", tooltip="Specify if the input data should be interpreted over time, which can be navigated with a time slider embedded into the map. Temporal data must have an explicit time column, or separate columns for each time period. "),
				config.MapType.UI(ui.input_select, id="MapType", label="Background Map", choices={"CartoDB Positron": "CartoDB", "OpenStreetMap": "OSM"}, tooltip="Specify the background map to plot your data on. CartoDB is a simpler map, while OSM is more highly annotated."),
				config.Opacity.UI(ui.input_slider, id="Opacity", label="Heatmap Opacity", min=0.0, max=1.0, step=0.1, tooltip="Specify the opacity of the heatmap. 1.0 indicates full opacity, while lower values make the background map more visible."),

				ui.HTML("<b>Colors</b>"),
				config.ColorMap.UI(ui.input_select, id="ColorMap", label="Color Map", choices=ColorMaps),
				config.Bins.UI(ui.input_numeric, id="Bins", label="# of Color Bins", min=3, max=253, step=1, tooltip="Specify the number of color bins to use. A higher number of color bins results in a smoother gradient between neighbouring values. Fewer bins results in more distinct colors. This feature does not apply to temporal heatmaps."),

				ui.HTML("<b>Range of Interest</b>"),
				config.ROI.UI(ui.input_checkbox, make_inline=False, id="ROI", label="Enable (Lower/Upper)", tooltip="Define a minimum and maximum bound (inclusive) for data. Select 'Remove' to ignore all values outside of the range. Select 'Round' to round values outside of the range to the maximum or minimum value."),
				config.ROI_Mode.UI(ui.input_radio_buttons, make_inline=False, id="ROI_Mode", label=None, choices=["Remove", "Round"], inline=True, tooltip="Remove data points outside the range of interest, or round them to the maximum or minimum value"),
				ui.layout_columns(
					config.Min.UI(ui.input_numeric,make_inline=False, id="Min", label=None, min=0, tooltip="Minimum displayed value in range of interest (inclusive)."),
					config.Max.UI(ui.input_numeric, make_inline=False, id="Max", label=None, min=0, tooltip="Maximum displayed value in range of interest (inclusive)."),
				),

				ui.download_button(id="DownloadHeatmap", label="Download HTML"),

				config.SettingType.UI(ui.input_radio_buttons, make_inline=False, 
				id="SettingType", label="Settings File Type", choices=[".txt", ".csv", ".tsv", ".xlsx"], inline=True),
				ui.download_button(id="DownloadSettings", label="Download Current Settings"),
			),
			padding="10px",
			gap="20px",
			width="300px",
		),
		MainTab(ui.nav_panel("GeoJSON", ui.output_data_frame("GeoJSON")), m_type=ui.output_ui),
		height="86vh",
	)
)

app = App(app_ui, server)
