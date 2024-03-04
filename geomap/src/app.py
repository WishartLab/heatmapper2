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

from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from folium import Map as FoliumMap, Choropleth
from folium.plugins import TimeSliderChoropleth
from pandas import DataFrame, to_datetime
from branca.colormap import linear
from json import loads
from shapely.geometry import shape

from shared import Table, Cache, NavBar, FileSelection, Pyodide, Filter, TemporalColumns

# Fine, Shiny
import branca, certifi, xyzservices


# Generated from dictionary.sh
Mappings = { "africa.geojson": "Africa", "akron.geojson": "Akron", "alameda.geojson": "Alameda", "albany.geojson": "Albany", "albuquerque.geojson": "Albuquerque", "amsterdam.geojson": "Amsterdam", "amusement-parks.geojson": "Amusement Parks", "anchorage.geojson": "Anchorage", "angers.geojson": "Angers", "angers-loire-metropole.geojson": "Angers Loire Metropole", "antwerp.geojson": "Antwerp", "apulia.geojson": "Apulia", "arlingtonva.geojson": "Arlingtonva", "asia.geojson": "Asia", "athens.geojson": "Athens", "atlanta.geojson": "Atlanta", "augsburg.geojson": "Augsburg", "austin.geojson": "Austin", "australia.geojson": "Australia", "austria-oberoesterreich.geojson": "Austria Oberoesterreich", "austria-states.geojson": "Austria States", "austria-steiermark.geojson": "Austria Steiermark", "bad-belzig.geojson": "Bad Belzig", "badenwuerttemberg-kreise.geojson": "Badenwuerttemberg Kreise", "baltimore.geojson": "Baltimore", "bari.geojson": "Bari", "basel.geojson": "Basel", "bayern.geojson": "Bayern", "belgium-arrondissements.geojson": "Belgium Arrondissements", "berlin.geojson": "Berlin", "bern-districts.geojson": "Bern Districts", "bern-quarters.geojson": "Bern Quarters", "birmingham.geojson": "Birmingham", "blacksburg.geojson": "Blacksburg", "blumenau.geojson": "Blumenau", "bogota.geojson": "Bogota", "boston.geojson": "Boston", "brandenburg.geojson": "Brandenburg", "brandenburg-municipalities.geojson": "Brandenburg Municipalities", "braunschweig.geojson": "Braunschweig", "brazil-states.geojson": "Brazil States", "bremen.geojson": "Bremen", "bronx.geojson": "Bronx", "brooklyn.geojson": "Brooklyn", "buenos-aires.geojson": "Buenos Aires", "calgary.geojson": "Calgary", "california-counties.geojson": "California Counties", "california-vista-points.geojson": "California Vista Points", "caltrain-stations.geojson": "Caltrain Stations", "canada.geojson": "Canada", "canberra.geojson": "Canberra", "caribbean-islands.geojson": "Caribbean Islands", "chapel-hill.geojson": "Chapel Hill", "charlotte.geojson": "Charlotte", "charlottesville.geojson": "Charlottesville", "chemnitz.geojson": "Chemnitz", "chesapeake.geojson": "Chesapeake", "chicago.geojson": "Chicago", "china.geojson": "China", "cincinnati.geojson": "Cincinnati", "cleveland.geojson": "Cleveland", "cologne.geojson": "Cologne", "colorado-counties.geojson": "Colorado Counties", "columbus.geojson": "Columbus", "copenhagen.geojson": "Copenhagen", "cuba.geojson": "Cuba", "dallas.geojson": "Dallas", "dane-county-municipalities.geojson": "Dane County Municipalities", "denmark-municipalities.geojson": "Denmark Municipalities", "denver.geojson": "Denver", "des-moines.geojson": "Des Moines", "detroit.geojson": "Detroit", "dresden.geojson": "Dresden", "dublin.geojson": "Dublin", "duesseldorf.geojson": "Duesseldorf", "durham.geojson": "Durham", "edmonton.geojson": "Edmonton", "eindhoven.geojson": "Eindhoven", "enschede.geojson": "Enschede", "esztergom.geojson": "Esztergom", "europe-1914.geojson": "Europe 1914", "europe-1938.geojson": "Europe 1938", "europe-capitals.geojson": "Europe Capitals", "europe.geojson": "Europe", "fairbanks.geojson": "Fairbanks", "fargo.geojson": "Fargo", "fort-lauderdale.geojson": "Fort Lauderdale", "france-departments.geojson": "France Departments", "france-regions.geojson": "France Regions", "frankfurt-main.geojson": "Frankfurt Main", "freiburg.geojson": "Freiburg", "geneva.geojson": "Geneva", "germany-capitals.geojson": "Germany Capitals", "germany.geojson": "Germany", "ghent.geojson": "Ghent", "gisborne.geojson": "Gisborne", "grand-rapids.geojson": "Grand Rapids", "greece-prefectures.geojson": "Greece Prefectures", "greece-regions.geojson": "Greece Regions", "hamburg.geojson": "Hamburg", "hampton.geojson": "Hampton", "hartford.geojson": "Hartford", "henderson.geojson": "Henderson", "honolulu.geojson": "Honolulu", "houston.geojson": "Houston", "hungary.geojson": "Hungary", "illinois-counties.geojson": "Illinois Counties", "india.geojson": "India", "indianapolis.geojson": "Indianapolis", "iran-provinces.geojson": "Iran Provinces", "ireland-counties.geojson": "Ireland Counties", "isle-of-man.geojson": "Isle Of Man", "italy-provinces.geojson": "Italy Provinces", "italy-regions.geojson": "Italy Regions", "james-city-county.geojson": "James City County", "japan.geojson": "Japan", "kaiserslautern.geojson": "Kaiserslautern", "kansas-city.geojson": "Kansas City", "korea.geojson": "Korea", "las-vegas.geojson": "Las Vegas", "leipzig.geojson": "Leipzig", "le-mans-cantons.geojson": "Le Mans Cantons", "lexington.geojson": "Lexington", "liberia-central.geojson": "Liberia Central", "liberia-east.geojson": "Liberia East", "liberia.geojson": "Liberia", "liberia-west.geojson": "Liberia West", "lombardy.geojson": "Lombardy", "london.geojson": "London", "london-underground.geojson": "London Underground", "long-beach.geojson": "Long Beach", "los-angeles-county.geojson": "Los Angeles County", "los-angeles.geojson": "Los Angeles", "louisville.geojson": "Louisville", "luxembourg-cantons.geojson": "Luxembourg Cantons", "luxembourg-communes.geojson": "Luxembourg Communes", "luzern.geojson": "Luzern", "macon.geojson": "Macon", "madrid-districts.geojson": "Madrid Districts", "madrid.geojson": "Madrid", "malaysia.geojson": "Malaysia", "manhattan-bridges.geojson": "Manhattan Bridges", "manhattan.geojson": "Manhattan", "melbourne.geojson": "Melbourne", "mexico.geojson": "Mexico", "miami.geojson": "Miami", "middle_east_countries.geojson": "Middle_east_countries", "milan.geojson": "Milan", "milwaukee.geojson": "Milwaukee", "minneapolis-cities.geojson": "Minneapolis Cities", "minneapolis.geojson": "Minneapolis", "mississauga.geojson": "Mississauga", "montreal.geojson": "Montreal", "moscow.geojson": "Moscow", "muenster.geojson": "Muenster", "new-haven.geojson": "New Haven", "new-orleans.geojson": "New Orleans", "new-york-areas-of-interest.geojson": "New York Areas Of Interest", "new-york-city-boroughs.geojson": "New York City Boroughs", "new-york-counties.geojson": "New York Counties", "nordrhein-westfalen.geojson": "Nordrhein Westfalen", "norfolk.geojson": "Norfolk", "north-america.geojson": "North America", "north-carolina-cities.geojson": "North Carolina Cities", "oakland.geojson": "Oakland", "oceania.geojson": "Oceania", "oklahoma-cities.geojson": "Oklahoma Cities", "oklahoma-counties.geojson": "Oklahoma Counties", "olympia.geojson": "Olympia", "oman.geojson": "Oman", "oman-provinces.geojson": "Oman Provinces", "orlando.geojson": "Orlando", "pakistan.geojson": "Pakistan", "paris.geojson": "Paris", "peaks.geojson": "Peaks", "philadelphia.geojson": "Philadelphia", "phoenix.geojson": "Phoenix", "pittsburgh.geojson": "Pittsburgh", "poland.geojson": "Poland", "poland-parks.geojson": "Poland Parks", "porirua.geojson": "Porirua", "portland.geojson": "Portland", "portugal.geojson": "Portugal", "potsdam.geojson": "Potsdam", "prague.geojson": "Prague", "providence.geojson": "Providence", "quebec.geojson": "Quebec", "queens.geojson": "Queens", "raleigh.geojson": "Raleigh", "red-deer.geojson": "Red Deer", "richmond.geojson": "Richmond", "riga.geojson": "Riga", "rio-de-janeiro.geojson": "Rio De Janeiro", "rochester.geojson": "Rochester", "rockville.geojson": "Rockville", "roller-coasters-fastest-steel.geojson": "Roller Coasters Fastest Steel", "romania.geojson": "Romania", "rome-rioni.geojson": "Rome Rioni", "rotterdam.geojson": "Rotterdam", "russia.geojson": "Russia", "sacramento.geojson": "Sacramento", "salt-lake-city.geojson": "Salt Lake City", "san-antonio.geojson": "San Antonio", "san-diego.geojson": "San Diego", "san-francisco.geojson": "San Francisco", "san-jose.geojson": "San Jose", "saskatoon.geojson": "Saskatoon", "savannah.geojson": "Savannah", "seattle.geojson": "Seattle", "seoul.geojson": "Seoul", "serbia.geojson": "Serbia", "silicon-valley.geojson": "Silicon Valley", "south-africa.geojson": "South Africa", "south-america.geojson": "South America", "southeast-asia.geojson": "Southeast Asia", "spain-communities.geojson": "Spain Communities", "spain-provinces.geojson": "Spain Provinces", "springfield.geojson": "Springfield", "stamford.geojson": "Stamford", "staten-island.geojson": "Staten Island", "st-louis.geojson": "St Louis", "st-petersburg.geojson": "St Petersburg", "surrey.geojson": "Surrey", "sweden-counties.geojson": "Sweden Counties", "switzerland.geojson": "Switzerland", "sydney.geojson": "Sydney", "szczecin.geojson": "Szczecin", "taiwan.geojson": "Taiwan", "tampa.geojson": "Tampa", "the-hague.geojson": "The Hague", "the-netherlands.geojson": "The Netherlands", "thessaloniki.geojson": "Thessaloniki", "toronto.geojson": "Toronto", "tucson.geojson": "Tucson", "turkey.geojson": "Turkey", "turku.geojson": "Turku", "ulm.geojson": "Ulm", "united-kingdom.geojson": "United Kingdom", "united-kingdom-regions.geojson": "United Kingdom Regions", "united-states-1810.geojson": "United States 1810", "united-states-big-cities.geojson": "United States Big Cities", "united-states.geojson": "United States", "united-states-international-airports.geojson": "United States International Airports", "united-states-mlb-stadiums.geojson": "United States Mlb Stadiums", "unna.geojson": "Unna", "utrecht.geojson": "Utrecht", "vancouver.geojson": "Vancouver", "venice.geojson": "Venice", "venlo.geojson": "Venlo", "vermont-counties.geojson": "Vermont Counties", "vienna.geojson": "Vienna", "villetta.geojson": "Villetta", "washington.geojson": "Washington", "wellington.geojson": "Wellington", "west-linn.geojson": "West Linn", "west-palm-beach.geojson": "West Palm Beach", "wiesenburg.geojson": "Wiesenburg", "williamsburg.geojson": "Williamsburg", "windsor.geojson": "Windsor", "winterthur.geojson": "Winterthur", "zurich-city.geojson": "Zurich City", "zurich.geojson": "Zurich", "world.geojson": "World"}

def server(input: Inputs, output: Outputs, session: Session):



	Info = {
		"example1.txt": "This example file is from the Open Data Portal. The data is from a carbon monoxide emissions study conducted by Environment Canada. The three columns represent results from 1990, 2000, and 2013.",
		"example2.txt": "This example file is from Statistics Canada. The data is adapted from New cases and age-standardized rate for primary cancer (based on the February 2014 CCR tabulation file), by cancer type and sex, Canada, provinces and territories. The columns represent new cancer cases (age-standardized rate per 100,000 population) from 2006 to 2010.",
		"example3.txt": "This example file is from the U.S. Centers for Disease Control and Prevention. The data is from Diagnosed Diabetes, Age Adjusted Rate (per 100), Adults - Total, 2013.",
		"example6.csv": "COVID 19 information reported by the Canadian Government, available at https://open.canada.ca/data/en/dataset/261c32ab-4cfd-4f81-9dea-7b64065690dc/resource/39434379-45a1-43d5-aea7-a7a50113c291"
	}

	DataCache = Cache("geomap")


	def LoadJSON():
		"""
		@brief Returns the GeoJSON depending on whether the user wants to use a provided one, or their own.
		@returns Either the path to the uploaded file, or the URL to the one provided by us (Folium supports both)
		"""

		URL = "https://raw.githubusercontent.com/kkernick/kkernick.github.io/main/geomap/data/" if Pyodide else "../data/"

		if input.JSONFile() == "Upload":
			file: list[FileInfo] | None = input.JSONUpload()
			if file is None:
				return URL + "canada.geojson"
			return file[0]["datapath"]
		else:
			return URL + input.JSONSelection()


	def LoadChoropleth(df, map):
		key, value = input.KeyColumn(), input.ValueColumn()

		# Add the heatmap and return.
		Choropleth(
				geo_data=LoadJSON(),
				name="choropleth",
				data=df,
				columns=[key, value],
				key_on="feature.properties.name",
				fill_color=input.ColorMap().lower(),
				fill_opacity=input.Opacity(),
				line_opacity=input.Opacity(),
				legend_name="Legend",
				bins=input.Bins()
		).add_to(map)


	async def LoadTemporalChoropleth(df, map):
		geojson = await DataCache.Download(LoadJSON())
		geojson = loads(geojson.decode('utf-8'))

		key, value = input.KeyColumn(), input.ValueColumn()

		# Check if we have a dedicated time column, or separate columns for each time slot.
		column = TemporalColumns.intersection(set(df.columns))

		values = value if column else df.columns[1:]
		match input.ColorMap():
			case "Inferno": cmap = linear.inferno.scale
			case "Magma": cmap = linear.magma.scale
			case "Plasma": cmap = linear.plasma.scale
			case "Viridis": cmap = linear.viridis.scale

		m, M = df[values].values.min(), df[values].values.max()
		colormap = cmap(m, M)

		style = {}

		if column:
			temporal_column = Filter(columns=df.columns, good_columns=TemporalColumns, bad_columns=[value])
			grouped = df.groupby(key)

			for i, (name, group) in enumerate(grouped):
				style[i] = {}
				for a, row in group.iterrows():

					# If the year isn't parsable, pray that it just works without parsing :)
					try: year = str(to_datetime(row[temporal_column].split()[0]).timestamp()).split('.')[0]
					except Exception: pass

					style[i][year] = {'color': colormap(row[value]), 'opacity': input.Opacity()}
				for feature in geojson["features"]:
							if feature["properties"]["name"] == name:
									feature["id"] = i

		else:
			# Convert DataFrame to style format
			for i, row in df.iterrows():
					name = row[key]
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


	async def LoadMap():
		"""
		@brief Generates a map with the provided information
		@returns the Folium.Map
		"""

		df = await DataCache.Load(input, copy=True)
		key, value = input.KeyColumn(), input.ValueColumn()

		# If the columns aren't defined, or aren't valid, don't do anything.
		if key not in df or value not in df: return

		# Give a placeholder map if nothing is selected, which should never really be the case.
		if df.empty: return FoliumMap((53.5213, -113.5213), tiles=input.MapType(), zoom_start=15)

		# Create map
		map = FoliumMap(tiles=input.MapType())

		# Load the choropleth.
		if input.Temporal(): await LoadTemporalChoropleth(df, map)
		else: LoadChoropleth(df, map)

		map.fit_bounds(map.get_bounds())
		return map


	@output
	@render.data_frame
	@reactive.event(input.Update, input.Reset, input.Example, input.File, ignore_none=False, ignore_init=False)
	async def LoadedTable(): return await DataCache.Load(input)


	@output
	@render.ui
	@reactive.event(input.Update, input.Reset, input.Example, input.File, input.KeyColumn, input.ValueColumn, input.JSONSelection, input.JSONUpload, input.Temporal, input.MapType, input.ColorMap, input.Opacity, input.Bins, ignore_none=False, ignore_init=False)
	async def Map(): return await LoadMap()


	@output
	@render.text
	def ExampleInfo(): return Info[input.Example()]


	@render.download(filename="table.csv")
	async def DownloadTable(): df = await DataCache.Load(input); yield df.to_string()


	@render.download(filename="heatmap.html")
	async def DownloadHeatmap(): m = await DataCache.Load(input); yield m.get_root().render()


	@reactive.Effect
	async def UpdateColumns():

		# Give options for the key and value columns
		df = await DataCache.Load(input)
		choices = df.columns.tolist()
		if choices:
			default_key = input.KeyColumn() if input.KeyColumn() is not None else df.columns[0]
			default_value = input.ValueColumn() if input.ValueColumn() is not None else df.columns[1]

			ui.update_select(id="KeyColumn", choices=choices, selected=default_key)
			ui.update_select(id="ValueColumn", choices=choices, selected=default_value)


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


app_ui = ui.page_fluid(

	NavBar("Geomap"),

	ui.layout_sidebar(
		ui.sidebar(

			FileSelection(
				examples={"example1.txt": "Example 1", "example2.txt": "Example 2", "example3.txt": "Example 3", "example6.csv": "Example 4"},
				types=[".csv", ".txt", ".xlsx"]
			),

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

			# Add the download buttons.
			"Download",
			ui.download_button("DownloadHeatmap", "Heatmap"),
			ui.download_button("DownloadTable", "Table"),

			id="SidebarPanel",
		),

		# Add the main interface tabs.
		ui.navset_tab(
				ui.nav_panel("Interactive", ui.output_ui("Map")),
				Table
		),
	)
)

app = App(app_ui, server)