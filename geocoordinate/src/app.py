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
from folium import Map as FoliumMap
from folium.plugins import HeatMap, HeatMapWithTime
from folium.raster_layers import ImageOverlay
from tempfile import NamedTemporaryFile
from matplotlib.pyplot import subplots
from scipy.stats import gaussian_kde
from numpy import vstack

from shared import Cache, NavBar, MainTab, FileSelection, Filter, ColumnType, TableOptions

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

	}

	DataCache = Cache("geocoordinate")
	Data = reactive.value(None)
	Valid = reactive.value(False)

	@reactive.effect
	@reactive.event(input.SourceFile, input.File, input.Example, input.Reset)
	async def UpdateData(): Data.set((await DataCache.Load(input))); Valid.set(False)


	def GetData(): return Table.data_view() if Valid() else Data()


	def GenerateMap(df, map, v_col, lon_col, lat_col):
		"""
		@brief Generates a standard heatmap
		@param df The DataFrame containing the data
		@param map The folium map to attach the heatmap to.
		"""

		opacity = input.Opacity()
		radius = input.Radius()
		blur = input.Blur()

		if input.Uniform():
			df["Value"] = [1] * len(df[lat_col])
			v_col = "Value"

		if input.Scale():
			HeatMap(list(zip(df[lat_col], df[lon_col], df[v_col])),
			min_opacity=opacity,
			radius=radius,
			blur=blur).add_to(map)

		else:
			latitude = df[lat_col]
			longitude = df[lon_col]
			values = df[v_col]

			# Threshold for density estimation
			threshold = 0.1
			stack = vstack([longitude, latitude])

			# Calculate kernel density estimation
			kde = gaussian_kde(stack)
			density = kde(stack)
			values = values + density * threshold

			# Generate the contour
			fig, ax = subplots()
			contour = ax.scatter(
				longitude, 
				latitude, 
				c=values * 100, 
				s=[radius] * len(values), 
				cmap="jet", 
				alpha=opacity, 
				linewidths=0
			)

			ax.axis('off')  # Turn off axes
			ax.get_xaxis().set_visible(False)  # Hide x-axis
			ax.get_yaxis().set_visible(False)  # Hide y-axis

			# Save the contour as an image, add to the map
			temp = NamedTemporaryFile(suffix='.png')
			fig.savefig(temp, format='png', bbox_inches='tight', pad_inches=0, transparent=True, dpi=500)
			image_overlay = ImageOverlay(
				image=temp.name,
				bounds=[[min(latitude), min(longitude)], [max(latitude), max(longitude)]],
				opacity=opacity,
				pixelated=False,
			)
			image_overlay.add_to(map)

		map.fit_bounds(map.get_bounds())


	def GenerateTemporalMap(df, map, t_col, v_col, lon_col, lat_col):
		"""
		@brief Generates a temporal heatmap
		@param df The DataFrame containing the data
		@param map The folium map to attach the heatmap to.
		"""

		# Ensure we have a valid time column
		if t_col not in df: return

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

		radius = input.Radius()
		opacity = input.Opacity()
		blur = input.Blur()

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

			# Set the Value Column Accordingly (Helper functions handle None)
			p.inc(message="Formatting...")
			if not input.Uniform():
				v_col = input.ValueColumn()
				if v_col not in df: return
			else:
				v_col = None

			# Get lat and lon, generate the map
			lon_col = Filter(df.columns, ColumnType.Longitude, only_one=True)
			lat_col = Filter(df.columns, ColumnType.Latitude, only_one=True)
			map = FoliumMap((df[lat_col][0], df[lon_col][0]), tiles=input.MapType())

			# IF ROI is defined, and we have a column of values to work with.
			if type(input.ROI()) is tuple and v_col is not None:

				# Get the min and max of the data, and the lower and upper bound
				m, M = df[v_col].min(), df[v_col].max()
				l, u = input.ROI()[0], input.ROI()[1]

				# Avoid race condition of new data being uploaded, but the ROI not being updated in tandem.
				if l >= m and u <= M:

					# Get values outside the upper and lower bound, drop them.
					oob = []
					for index, row in df.iterrows():
						v = row[v_col]
						if v < l or v > u: oob.append(index)
					df = df.drop(oob)

			# Generate the right heatmap.
			p.inc(message="Plotting...")
			if input.Temporal(): 
				t_col = input.TimeColumn()
				GenerateTemporalMap(df, map, t_col, v_col, lon_col, lat_col)
			else: GenerateMap(df.copy(deep=True), map, v_col, lon_col, lat_col)
			return map


	@output
	@render.text
	@reactive.event(input.SourceFile, input.Example)
	def ExampleInfo(): return Info[input.Example()]


	@render.download(filename="table.csv")
	def DownloadTable(): yield GetData().to_string()


	@render.download(filename="heatmap.html")
	def DownloadHeatmap(): m = LoadMap(); yield m.get_root().render()


	@reactive.Effect
	def UpdateColumns():
		df = GetData()
		if not input.Uniform():
			if not Filter(df.columns, ColumnType.Value, ui_element="ValueColumn"):
				ui.update_checkbox(id="Uniform", value=True)
		if input.Temporal(): Filter(df.columns, ColumnType.Time, ui_element="TimeColumn")


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
			ui.update_slider(id="Opacity", value=0.7, min=0.0, max=1.0, step=0.1)
			ui.update_slider(id="Radius", value=25, min=5, max=100, step=5)
			ui.update_slider(id="Blur", value=15, min=1, max=30, step=1)


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
				examples={
					"example1.txt": "Example 1",
					"example2.txt": "Example 2",
					"example3.txt": "Example 3",
					"example1.csv": "Example 4",
					"example21.csv": "Example 5",
					"example3.csv": "Example 6"
				},
				types=[".csv", ".txt", ".dat", ".tsv", ".tab", ".xlsx", ".xls", ".odf"],
				project="Geocoordinate"
			),

			TableOptions(),

			ui.input_checkbox(id="Temporal", label="Temporal Heatmap"),
			ui.input_checkbox(id="Uniform", label="Uniform Values"),
			ui.input_checkbox(id="Scale", label="Scale Intensity with Zoom", value=True),

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

			ui.input_slider(id="Opacity", label="Heatmap Opacity", value=0.7, min=0.0, max=1.0, step=0.1),
			ui.input_slider(id="Radius", label="Size of Points", value=25, min=5, max=100, step=5),
			ui.input_slider(id="Blur", label="Blurring", value=15, min=1, max=30, step=1),
			ui.input_slider(id="ROI", label="Range of Interest", value=(0,0), min=0, max=100),

			# Add the download buttons.
			ui.download_button("DownloadHeatmap", "Heatmap"),
			ui.download_button("DownloadTable", "Table"),
		),

		MainTab(m_type=ui.output_ui),
	)
)

app = App(app_ui, server)
