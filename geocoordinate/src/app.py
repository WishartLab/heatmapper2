#
# Heatmapper
# Geocoordinate
#
# This file contains the ShinyLive application for Geocoordinate Heatmapper.
# It can be run with the following command within this directory:
#        shinylive export . [site]
# Where [site] is the destination of the site folder.
#
# If you would rather deploy the application as a PyShiny application,
# run the following command within this directory:
#        shiny run
#
#


from shiny import App, reactive, render, ui
from folium import Map as FoliumMap, Circle, Rectangle
from folium.plugins import HeatMap, HeatMapWithTime
from folium.raster_layers import ImageOverlay
from tempfile import NamedTemporaryFile
from matplotlib.pyplot import subplots
from scipy.stats import gaussian_kde
from scipy.interpolate import griddata
from scipy.io import netcdf_file
from numpy import vstack, convolve, ones, float32, linspace, meshgrid, isnan
from branca.colormap import LinearColormap
from pandas import DataFrame
from io import StringIO
from math import sqrt

from shared import Cache, NavBar, MainTab, FileSelection, Filter, ColumnType, TableOptions, InitializeConfig, Error, Update, Msg, File, MapTiles

try:
    from user import config
except ImportError:
    from config import config

# Fine, Shiny
import branca, certifi, xyzservices, requests


def server(input, output, session):

    Info = {
        "example1.txt": '<u>Input type:</u> .txt Data<br><u>Contents:</u> Deaths from a cholera outbreak in 1854. John Snow used this data in conjunction with local pump locations as evidence that cholera is spread by contaminated water. <br><u>Source:</u> A digitised version of the data is available <a href="https://blog.rtwilson.com/john-snows-cholera-data-in-more-formats/"; target="_blank">courtesy of Robin Wilson</a>.',
        "example3.txt": '<u>Input type:</u> .txt Data<br><u>Contents:</u> The location of traffic signals in Toronto, Canada.<br><u>Source:</u> <a href="https://open.toronto.ca/dataset/traffic-signals-tabular/"; target="_blank">Toronto Open Data</a>',
        "example21.csv": '<u>Input type:</u> .csv Data <br><u>Contents:</u> A parsed version of the Northeast and North Central Pacific hurricane database (HURDAT2) 2000-2022. Data includes the date and time, location, and intensity of each hurricane. This example can be viewed as a temporal heat map.<br><u>Source:</u> <a href="https://www.nhc.noaa.gov/data/"; target="_blank">www.nhc.noaa.gov/data/</a>',
        "NYCtaxi-2013.csv": '''<u>Input type:</u> .csv Data <br><u>Contents:</u> This example contains a subset of 2013 taxi trip data from the NYC Taxi and Limousine Commission. The columns are 'latitude', 'longitude', 'value', and 'time'. The time column is the pickup time of the taxi service, and the number in the value column is the number of passengers in that taxi. This example can be viewed as a temporal heat map.<br><u>Source:</u> <a href="https://www.nyc.gov/site/tlc/about/tlc-trip-record-data.page"; target="_blank">www.nyc.gov</a>'''
    }


    def NCDataFrame(nc_path, p=None):
        """
        @brief Takes a netcdf_file path and returns an appropriate DataFrame
        @param nc: The netcdf_file path
        @returns    A dataframe with Latitude, Longitude, and optional Time alongwith
                            all variables bounded by those values.
        """

        nc = netcdf_file(nc_path, mmap=False)
        if p: p.inc(message="Extracting keys...")
        dimensions = list(nc.dimensions.keys())
        latitude, longitude, time = Filter(dimensions, ColumnType.Latitude), Filter(dimensions, ColumnType.Longitude), Filter(dimensions, ColumnType.Time)

        if not latitude or not longitude:
            Error("NC files require both a latitude and longitude field.")
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
                Data.set((await DataCache.Load(input, p=p)))
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
            Error(ui.HTML('File could not be loaded! Data can be uploaded as a .csv, .tsv, .txt, .xslx, .dat, .tab, or .odf file. Error'), e)


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
        if "Color by Density" in config.Features():
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

                    if config.RenderShape() == "Circle":
                        Circle(
                            location=[row[lat_col], row[lon_col]],
                            radius=radius,
                            color=color,
                            fill=True,
                            opacity=opacity,
                            fill_opacity=opacity,
                            stroke=False
                        ).add_to(map)
                    else:
                        lat, lon = row[lat_col], row[lon_col]
                        rect_radius = radius / 100000
                        Rectangle(
                            bounds=[[lat - rect_radius, lon - rect_radius], [lat + rect_radius, lon + rect_radius]],
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
        # display placeholder message if no table uploaded yet
        if df is None:
            return DataFrame({"Note": ["No data to display! Please upload your data or select an example data set in the sidebar."]})

        if File(input).endswith(".nc"):
            df = NCDataFrame(df)

        try:
            Valid.set(True)
            return render.DataGrid(df, editable=True)
        except Exception as e:
            return DataFrame({'Note': ['Table could not be rendered. Ensure your data is properly formatted, and your data file is a compatible type. <a href="https://github.com/WishartLab/heatmapper2/wiki/Format#geocoordinate:~:text=the%20Table%20tab)-,Geocoordinate,-Geocoordinate%20takes%20a"; target="_blank"; rel=”noopener noreferrer;>Read more</a>.']})


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
            <h1>Geocoordinate</h1>
            This heatmap maps values onto geospatial coordinates (latitude and longitude). Upload a data file in the sidebar to get started, or select 'Example' to check out a pre-loaded example. Navigate to the 'Heatmap' tab to see the heatmap, or 'Table' to look at the input data.
                 
            <br><br>
            <img src="https://github.com/WishartLab/heatmapper2/wiki/assets/Geocoordinate.png" alt="Image"; style="max-width:500px;">
                 
            <br><br><h3>Format</h3>
            <i>Input data can be formatted as follows:</i>
                <ul>
                <li><u>Static Heatmaps</u> require 'Latitude' and 'Longitude' columns. An optional 'Value' column can be used to color data points based on value - a uniform value is applied if this column does not exist (see examples '1: Cholera Deaths' and '2: Traffic Signals').</li>
                <li><u>Temporal Heatmaps</u> require 'Latitude' and 'Longitude' columns plus an additional 'Time' column. An optional 'Value' column can be used to color data points based on value - a uniform value is applied if this column does not exist (see examples '3: Hurricanes' and '4: NYC Taxi Trips').</li>
                </ul>
                 
            <i>Geocoordinate heatmaps can be generated from the following file formats:</i>
                <li>.csv</li>
                <li>.dat</li>
                <li>.odf</li>
                <li>.tab</li>
                <li>.tsv</li>
                <li>.txt</li>
                <li>.xls</li>
                <li>.xlsx</li>
            
            <br><h3>Interface</h3>
            Click on the '?' icon beside sidebar options to read more about them.
        """)


    def GenerateHeatmap():
        with ui.Progress() as p:
            p.inc(message="Loading input...")
            df = GetData()
            if df is None: 
                return ui.HTML('No data to display! <br>Please upload your data or select an example data set in the sidebar.')

            if File(input).endswith(".nc") and not Valid():
                df = NCDataFrame(GetData())
                if df is None: 
                    return ui.HTML('The heat map could not be rendered. <br><br>Please ensure your .nc input file is properly formatted. <br>Input data should contain latitude and longitude columns, with optional value and time columns. <br>More information on formatting is available in the <a href="https://github.com/WishartLab/heatmapper2/wiki/Format#geocoordinate:~:text=the%20Table%20tab)-,Geocoordinate,-Geocoordinate%20takes%20a"; target="_blank"; rel=”noopener noreferrer;>Wiki</a>.')
                nc = True
            else:
                df = df.copy(deep=True)
                nc = False

            p.inc(message="Formatting...")

            lon_col = Filter(df.columns, ColumnType.Longitude)
            lat_col = Filter(df.columns, ColumnType.Latitude)
            if lat_col is None or lon_col is None: 
                return ui.HTML('The heat map could not be rendered. <br><br>Please ensure your input data contains a latitude column (named "latitude" or "lat"), and a longitude column (named "longitude", "long", or "lon"). Column names are case-insensitive. <br>More information on formatting is available in the <a href="https://github.com/WishartLab/heatmapper2/wiki/Format#geocoordinate:~:text=the%20Table%20tab)-,Geocoordinate,-Geocoordinate%20takes%20a"; target="_blank"; rel=”noopener noreferrer;>Wiki</a>.')

            v_col = config.ValueColumn()
            if v_col == "Uniform":
                df["Heatmapper_Uniform_Values"] = [1] * len(df[lat_col])
                v_col = "Heatmapper_Uniform_Values"

            t_col = config.TimeColumn()

            tile_info = MapTiles[config.MapType()]
            map = FoliumMap((df[lat_col][0], df[lon_col][0]), tiles=tile_info["tiles"], attr=tile_info["attr"])

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
                    Error("No locations to display! Check your Range of Interest and ensure the Value Column is properly set.")
                    return ui.HTML("No locations to display! Check your Range of Interest and ensure the Value Column is properly set.")

            if config.Interpolation() != 1:
                p.inc(message="Interpolating...")
                resolution = config.Interpolation()

                points = df[[lat_col, lon_col]].values
                values = df[v_col].values

                # Create a grid based on the resolution
                lon_min, lon_max = df[lon_col].min(), df[lon_col].max()
                lat_min, lat_max = df[lat_col].min(), df[lat_col].max()

                lon_new = linspace(lon_min, lon_max, num=int(sqrt(len(values) * resolution)) + 1)
                lat_new = linspace(lat_min, lat_max, num=int(sqrt(len(values) * resolution)) + 1)

                lon_grid, lat_grid = meshgrid(lon_new, lat_new)

                # Interpolate the values
                grid_z = griddata(points, values, (lat_grid, lon_grid), method='cubic')

                # Flatten the grid for the new DataFrame
                lon_flat = lon_grid.flatten()
                lat_flat = lat_grid.flatten()
                value_flat = grid_z.flatten()

                # Remove NaN values that could be introduced by interpolation
                #valid_mask = ~isnan(value_flat)
                #lon_flat = lon_flat[valid_mask]
                #lat_flat = lat_flat[valid_mask]
                #value_flat = value_flat[valid_mask]

                # Create a new DataFrame with interpolated values
                df = DataFrame({
                        lon_col: lon_flat,
                        lat_col: lat_flat,
                        v_col: value_flat
                })
                print(df)

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
        except KeyError:
            return ui.HTML('The heat map could not be rendered due to a Key Error. <br><br>Please ensure your input data contains latitude and longitude columns, with optional value and time columns. <br>More information on formatting is available in the <a href="https://github.com/WishartLab/heatmapper2/wiki/Format#geocoordinate:~:text=the%20Table%20tab)-,Geocoordinate,-Geocoordinate%20takes%20a"; target="_blank"; rel=”noopener noreferrer;>Wiki</a>.')
        except Exception as e:
            #Error(f"Failed to generate heatmap", e)
            return ui.HTML('The heat map could not be rendered. <br><br>Please ensure your input data contains latitude and longitude columns, with optional value and time columns. <br>More information on formatting is available in the <a href="https://github.com/WishartLab/heatmapper2/wiki/Format#geocoordinate:~:text=the%20Table%20tab)-,Geocoordinate,-Geocoordinate%20takes%20a"; target="_blank"; rel=”noopener noreferrer;>Wiki</a>.')


    @output
    @render.ui
    @reactive.event(input.Update)
    def HeatmapReactive():
        try:
            return GenerateHeatmap()
        except Exception as e:
            return ui.HTML('The heat map could not be rendered. <br><br>Please ensure your input data contains latitude and longitude columns, with optional value and time columns. <br>More information on formatting is available in the <a href="https://github.com/WishartLab/heatmapper2/wiki/Format#geocoordinate:~:text=the%20Table%20tab)-,Geocoordinate,-Geocoordinate%20takes%20a"; target="_blank"; rel=”noopener noreferrer;>Wiki</a>.')


    @reactive.effect
    @reactive.event(input.ExampleInfo)
    def ExampleInfo():
        Msg(ui.HTML(Info[input.Example()]))


    @render.download(filename=lambda: f"table{config.TableType()}")
    def DownloadTable(): 
        data = GetData()
        
        # return error if no data to download
        if data is None:
            Error("The downloaded table is empty! Please upload your data or select an example data set in the sidebar.")
        
        file_contents = data.to_string()
        yield file_contents


    @render.download(filename="heatmap.html")
    def DownloadHeatmap(): m = GenerateHeatmap(); yield m.get_root().render()


    @render.download(filename=lambda: f"settings{config.SettingType()}")
    def DownloadSettings(): 
        """
        Download a table file containing current config settings
        """
        yield f"Time Column:\t{config.TimeColumn()}\nValue Column:\t{config.ValueColumn()}\nRender Mode:\t{config.RenderMode()}\nVector Shape:\t{config.RenderShape()}\nMap Type:\t{config.MapType()}\nData Point Size:\t{config.Radius()}\nOpacity:\t{config.Opacity()}\nBlurring:\t{config.Blur()}\nRange of Interest:\t{config.ROI()}\nRange of Interest Mode:\t{config.ROI_Mode()}\nRange of Interest Min:\t{config.Min()}\nRange of Interest Max:\t{config.Max()}\nColour by Density:\t{config.Features()}"


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

    ui.panel_title(title=None, window_title="Geocoordinate"),
    NavBar(),

    ui.layout_sidebar(
        ui.sidebar(

            FileSelection(
                examples={
                    "example1.txt": "1: Cholera Deaths",
                    "example3.txt": "2: Traffic Signals",
                    "example21.csv": "3: Hurricanes",
                    "NYCtaxi-2013.csv": "4: NYC Taxi Trips",
                },
                types=[".csv", ".txt", ".dat", ".tsv", ".tab", ".xlsx", ".xls", ".odf", ".nc"],
                project="Geocoordinate"
            ),

            TableOptions(config),

            ui.panel_conditional(
                "input.MainTab === 'HeatmapTab'",

                Update(),

                ui.HTML("<b>Columns</b>"),
                config.TimeColumn.UI(ui.input_select, id="TimeColumn", label="Time Column", choices=[], multiple=False, tooltip="Optional: Specify a time column to plot data over time. If an explicit time column is specified, data can be visualized temporally with a media-player-like interface (play, pause, rewind, and frame speed options). If 'None' is selected, the heatmap will be static."),
                config.ValueColumn.UI(ui.input_select, id="ValueColumn", label="Value Column", choices=[], multiple=False, tooltip="If a column from the input data is specified, values from that column will be associated with each latitude, longitude point, and the point will be colored based on its value. If 'Uniform' is selected, data points will be assigned a uniform value and colored uniformly on the map."),

                ui.HTML("<b>Heatmap</b>"),
                config.RenderMode.UI(ui.input_select, id="RenderMode", label="Render Mode", choices=["Raster", "Vector"], tooltip="Display data as discrete vector points, or a smooth raster shape (vector does not apply to temporal heatmaps). The intensity of raster points scales when the map is zoomed in or out. Vector points maintain a constant intensity regardless of zoom, but are more computationally expensive."),
                config.RenderShape.UI(ui.input_select, id="RenderShape", label="Vector Shape", choices=["Circle", "Rectangle"], tooltip="Specify the shape of vector points. Rectangular points are useful for contiguous data (like temperature or rainfall), while circular points are useful for discrete data (like disease cases or wildlife sightings)."),
                config.MapType.UI(ui.input_select,id="MapType", label="Map Type", choices={"CartoDB Voyager": "CartoDB Voyager", "OpenStreetMap": "OSM"}, tooltip="Specify the background map to plot your data on. CartoDB Voyager is a detailed, authenticated CartoDB style, while OSM is more highly annotated."),

                config.Radius.UI(ui.input_numeric, id="Radius", label="Data Point Size", min=5, tooltip="Specify how large each data point should be on the map."),

                config.Opacity.UI(ui.input_slider, id="Opacity", label="Opacity", min=0.0, max=1.0, step=0.1, tooltip="Specify the opacity of the heatmap. 1.0 indicates full opacity, while lower values make the background map more visible."),
                config.Blur.UI(ui.input_slider, id="Blur", label="Blurring", min=1, max=30, step=1, tooltip="Specify how much neighbouring points bleed into one another. Higher values make the heatmap appear more homogeneous, while lower values emphasize individual points. This applies to raster heatmaps only."),

                # TODO: FIX INTERPOLATION
                # config.Interpolation.UI(ui.input_numeric, id="Interpolation", label="Inter", min=1, max=10, step=0.1, tooltip="Calculate intermediate values between points. This can lead to artifacts if data is not contiguous. (METHOD? APPLIES TO VECTOR AND RASTER?)"),

                ui.HTML("<b>Range of Interest</b>"),
                config.ROI.UI(ui.input_checkbox, make_inline=False, id="ROI", label="Enable Range of Interest", tooltip="Define a minimum and maximum bound (inclusive) for data points. Select 'Remove' to ignore all data points outside of the range. Select 'Round' to round data points outside of the range to the maximum or minimum value. This setting is not applicable if 'Uniform' values are used."),
                config.ROI_Mode.UI(ui.input_radio_buttons, make_inline=False, id="ROI_Mode", label=None, choices=["Remove", "Round"], inline=True, tooltip="Remove data points outside the range of interest, or round them to the maximum or minimum value"),
                ui.layout_columns(
                    config.Min.UI(ui.input_numeric,make_inline=False, id="Min", label=None, min=0, tooltip="Minimum displayed value in range of interest (inclusive)."),
                    config.Max.UI(ui.input_numeric, make_inline=False, id="Max", label=None, min=0, tooltip="Maximum displayed value in range of interest (inclusive)."),
                ),

                ui.HTML("<b>Features</b>"),
                config.Features.UI(
                    ui.input_checkbox_group, id="Features", make_inline=False, label=None,
                    choices=["Color by Density"], selected=None, tooltip="Visualize the distribution of data points, rather than their assigned value. Red indicates higher density areas while indigo indicates lower density areas. Density is calculated using Gaussian Kernal Density Estimation (KDE)."
                ),

                # Button to download heat map
                ui.download_button("DownloadHeatmap", "Download HTML"),

                # Button to download current config settings
                config.SettingType.UI(ui.input_radio_buttons, make_inline=False, 
                id="SettingType", label="Settings File Type", choices=[".txt", ".csv", ".tsv", ".xlsx"], inline=True),
                ui.download_button(id="DownloadSettings", label="Download Current Settings"),
            ),
            padding="10px",
            gap="20px",
            width="300px",
        ),

        MainTab(m_type=ui.output_ui),
        height="86vh",
    )
)

app = App(app_ui, server)
