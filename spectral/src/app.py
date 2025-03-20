#
# Heatmapper
# Spectral
#
# This file contains the ShinyLive application for Spectral Heatmapper.
# It can be run with the following command within this directory:
#        shinylive export . [site]
# Where [site] is the destination of the site folder.
#
# If you would rather deploy the application as a PyShiny application,
# run the following command within this directory:
#        shiny run
#
#

# ShinyLive needs this imported
import regex

from io import BytesIO
from matplotlib.pyplot import subplots, colorbar, style, get_cmap, close as fig_close
#from matplotlib.collections import LineCollection
#from matplotlib.colors import Normalize
#from mpl_toolkits.mplot3d.art3d import Line3DCollection
from numpy import zeros, unique, array, concatenate, asarray, hstack, column_stack, newaxis, full_like, zeros_like, meshgrid, full, where, linspace, log1p
from pandas import DataFrame
from plotly.graph_objects import Surface, Figure, Scatter3d
from plotly.io import renderers as r
from pymzml.run import Reader
#from scipy.spatial.distance import squareform
from scipy.interpolate import griddata
from shiny import App, reactive, render, ui
from tempfile import NamedTemporaryFile

from shared import Cache, MainTab, NavBar, FileSelection, Filter, ColumnType, TableOptions, InitializeConfig, Update, Msg, File, InterpolationMethods, Error

try:
    from user import config
except ImportError:
    from config import config


# Suppress intermediate Plotly displays
r.default = None


def server(input, output, session):
    # Information regarding example files.
    Info = {
        #"1min.mzml": "An Example mzML from https://github.com/HUPO-PSI/mzML"
        "BSA1-subset.mzML": "A subset of 8 spectra from the OpenMS Bovine Serum Albumin sample<br>https://github.com/OpenMS/OpenMS/tree/develop/share/OpenMS/examples/BSA"
    }


    def Hash():
        """
        @brief Compute a Hash String for the Main Heatmap
        @return A hash of all the inputs used for the heatmap.
        """
        tab = input.MainTab()
        if tab == "HeatmapTab":
            return [
                File(input),
                config.ColorMap(),
                config.Features(),
                config.TextSize(),
                config.ID(),
                config.Peaks(),
                config.DPI(),
                config.DimensionRT(),
                config.DimensionMZ(),
                input.mode(),
            ]
        elif tab == "SimilarityTab":
            return [
                File(input),
                config.ColorMap(),
                config.Features(),
                config.TextSize(),
                config.ID(),
                config.DPI(),
                config.Interpolation(),
                input.mode(),
            ]


    def CreateErrorImg(text, color, type="text"):
        """
        @brief Generates an image of the provided text
        @param text: The text to display as an error
        @param color: Hex color code for the text
        @param inputs: A list of all the inputs for caching (from HashString())
        @returns 
        """
        if type == "img":
            # create image with error text
            file = NamedTemporaryFile(delete=False, suffix=".png")
            fig, ax = subplots()
            ax.text(0, 50, text, color=color, fontsize=12)
            ax.set_xlim(0, 200)
            ax.set_ylim(0, 100)
            # make axes transparent
            [ax.spines[side].set_alpha(0.0) for side in ["top", "bottom", "left", "right"]]
            ax.tick_params(axis='both', which='both', reset=False, color=[0,0,0,0], labelcolor=[0,0,0,0])
            # get image for display
            fig.savefig(file.name, format="png", dpi=100, bbox_inches="tight")
            fig_close(fig)
            img: types.ImgData = {"src": file.name, "width": "400px"}
            return img
        else:
            return ui.HTML(text)


    def HandleData(path, p=None):
        """
        @brief A custom Data Handler for the Cache.
        @param path: The path to the file
        @returns A data object from the cache.
        @info This Data Handler supports .mzML files
        """
        if path.suffix.lower() == ".mzml": 
            return Reader(path.resolve(), build_index_from_scratch=True)
        else: 
            return None
    DataCache = Cache("spectral", DataHandler=HandleData)
    Data = reactive.value(None)
    Valid = reactive.value(False)
    InitializeConfig(config, input)


    @reactive.effect
    @reactive.event(input.SourceFile, input.File, input.Example, input.Reset)
    async def UpdateData():
        p = ui.Progress()
        try:
            Data.set((await DataCache.Load(input, p=p, default=None)))
            Valid.set(False)
            DataCache.Invalidate(File(input))
        except:
            p.close()
            Error("File could not be loaded!\nPlease upload data as a .mzml file.")
            return

        reader = Data()
        if reader is None: return

        ids = set()
        first = None
        try:
            for spectra in reader:
                ids.add(spectra.ID)
                print(f"added {spectra.ID}")
                if first is None: first = spectra.ID
            ui.update_select(id="ID", selected=[first], choices=list(ids))
        except:
            p.close()
            Error("File could not be parsed. \nPlease check the formatting of your .mzml file.")
            return


    @output
    @render.data_frame
    def Table(): 
        Valid.set(True)
        table_data = []
        # Data() is a pymzml.run.Reader object
        data = Data()
        
        # display placeholder message if no file has been uploaded
        if data is None:
            return DataFrame({"Note": ["No data to display! Please upload your data or select an example data set in the sidebar."]})

        # goal is to turn the pymzml.run.Reader object into a dataframe
        for spectrum in data:
            if spectrum.ms_level:
                polarity = None
                if spectrum["negative scan"]:
                    polarity = "negative"
                elif spectrum["positive scan"]:
                    polarity = "positive"

                spectrum_dict = {
                    "ID": spectrum.ID,
                    "MS Level": spectrum.ms_level,
                    "RT (minutes)": spectrum.scan_time_in_minutes(),
                    "Polarity": polarity,
                    "Raw Peaks": len(spectrum.peaks("raw")),
                    "Centroided Peaks": len(spectrum.peaks("centroided")),
                    "Reprofiled Peaks": len(spectrum.peaks("reprofiled")),
                    "Highest Intensity": spectrum.highest_peaks(1)[0],
                    "Extreme Values (mz)": spectrum.extreme_values("mz"),
                    #"Total Ion Current": spectrum.tic,
                }
                table_data.append(spectrum_dict)
        df = DataFrame(table_data)
        return render.DataGrid(df, editable=False)


    # Info text in welcome tab
    @render.ui
    def Welcome():
        return ui.HTML("""
            <h1>Spectral Heatmaps</h1>
            Spectral heatmaps display Mass Spectrometry data in a 3D plot, or compute a similarity matrix between spectra. <br><br>
            Upload a .mzML file in the sidebar to get started, or select 'Example' to check out a pre-loaded example. <br><br>
            Navigate to the 'Heatmap' tab to see the heatmap, 'Similarity' to view the similarity matrix, or 'Table' to look at the input data.
                 
            <br><br>
            <img src="https://github.com/WishartLab/heatmapper2/wiki/assets/Spectral.png" alt="Image"; style="max-width:500px;">
                 
            <br><br><h3>Format</h3>
            Uploaded files should have the extension <b>.mzML</b>.<br>
                 
            <br><h3>Interface</h3>
            Click on the '?' icon beside sidebar options to read more about them.
            """)


    def GenerateSimilarity():
        """
        @brief Generate the Similarity matrix.
        """
        
        # Sometimes the MainTab doesn't update immediately, despite calling the function.
        # We just return nothing if we aren't actually on the tab to avoid redundant calculation.
        if input.MainTab() != "SimilarityTab": return

        inputs = Hash()
        if not DataCache.In(inputs):
            with ui.Progress() as p:
                p.inc(message="Loading input...")
                reader = Data()
                
                # display placeholder message if no data has been uploaded
                if reader is None: 
                    return CreateErrorImg("No data to display! \nPlease upload your data or select an example data set in the sidebar.", "#027bc2", "img")

                # Get all the spectra the user wants.
                distances = {}
                #indices = [int(i) for i in config.ID()]
                indices = [i for i in config.ID()]
                spectra = []
                for s in reader:
                    if str(s.ID) in indices: spectra.append(s)

                for s in spectra:
                    sid1 = str(s.ID)
                    distances[sid1] = {}
                    for s2 in spectra:
                        p.inc(message=f"Computing similarity of spectra {s.ID} and {s2.ID}")
                        sid2 = str(s2.ID)
                        
                        # If they're the same spectra, set it to 1.
                        if s.ID == s2.ID:
                            distances[sid1][sid2] = 1.0
                        
                        # If the second spectra exists, use that pre-calculated result.
                        elif sid2 in distances:
                            distances[sid1][sid2] = distances[sid2][sid1]
                            
                        # Otherwise call the similarity function.
                        else:
                            distances[sid1][sid2] = s.similarity_to(s2)

                p.inc(message="Plotting")
                df = DataFrame(data=distances, columns=indices, index=indices, dtype="float")
                fig, ax = subplots()
                interpolation = config.Interpolation().lower()
                plot = ax.imshow(df, cmap=config.ColorMap().lower(), interpolation=interpolation, aspect="equal")

                # Visibility of features
                try:
                    if "legend" in input.Features():
                        cbar = colorbar(plot, ax=ax, label="Value", pad=0.1)
                        cbar.ax.tick_params(labelsize=config.TextSize())
                except Exception: pass

                if "y" in config.Features():
                    ax.set_yticklabels(df.columns)
                    ax.set_yticks(range(len(df.columns)))
                    ax.tick_params(axis="y", labelsize=config.TextSize())
                else: ax.set_yticklabels([])

                if "x" in config.Features():
                    ax.set_xticklabels(df.columns)
                    ax.set_xticks(range(len(df.columns)))
                    ax.tick_params(axis="x", labelsize=config.TextSize())
                else: ax.set_xticklabels([])

                b = BytesIO()
                dpi = config.DPI()
                if dpi < 5: dpi = 5
                elif dpi > 1500: dpi = 1500
                fig.savefig(b, format="png", dpi=dpi)
                b.seek(0)
                DataCache.Store(b.read(), inputs)

        b = DataCache.Get(inputs)
        with NamedTemporaryFile(delete=False, suffix=".png") as temp:
            temp.write(b)
            temp.close()
            img: types.ImgData = {"src": temp.name, "height": f"{config.Size()}vh"}
            return img


    def GenerateHeatmap():
        """
        @brief Generate the main heatmap
        """
        
        if input.MainTab() != "HeatmapTab": return

        inputs = Hash()
        if not DataCache.In(inputs):
            with ui.Progress() as p:
                p.inc(message="Loading input...")
                reader = Data()

                # display placeholder message if no data is uploaded
                if reader is None: 
                    return CreateErrorImg("No data to display!<br>Please upload your data or select an example data set in the sidebar.", "#027bc2", "text")

                # We additionally cache interpolation.
                rt_dimension = config.DimensionRT()+1
                mz_dimension = config.DimensionMZ()
                interpolation_cache = [File(input), config.Peaks(), rt_dimension, mz_dimension, "Interpolation"]
                print(f"rt_dim, mz_dim:\t{rt_dimension}, {mz_dimension}")
                if not DataCache.In(interpolation_cache):

                    # try to find based on selected peak type, parse normally if that fails
                    #peaks = config.Peaks().lower()
                    
                    mz_values = []
                    rt_values = []
                    intensities = []
                    
                    # iterate through all spectra, 
                    # get m/z, intensities, and rt for MS1 spectra
                    for spectrum in reader:
                        p.inc(message=f"Reading Spectra {spectrum.ID}")
                        try:
                            if spectrum.ms_level == 1:
                                retention_time = spectrum.scan_time_in_minutes()
                                
                                # Access the peaks arrays directly
                                mzs = spectrum.mz
                                intens = spectrum.i

                                if mzs is not None and intens is not None:
                                    mz_values.extend(mzs)
                                    rt_values.extend([retention_time] * len(mzs))
                                    intensities.extend(intens)
                        except Exception as e:
                            p.inc(message=f"Warning: Could not process spectrum: {e}")
                            continue
            
                    # ??? does this work ???
                    if not mz_values:
                        Error("No Spectra in File!")
                        return CreateErrorImg("No spectra could be identified in the input file.", "#027bc2", "text")
                        
                    # get min and max
                    intensity_min = min(intensities)
                    intensity_max = max(intensities)
                    print(f"---INTENSITIES---\n{intensity_min}\t{intensity_max}")
                    vm, vM, rm, rM = min(mz_values), max(mz_values), min(rt_values), max(rt_values)
                    print(f"---MZ---\n{vm}\t{vM}")
                    print(f"---RT---\n{rm}\t{rM}")

                    # Create a grid for mz and rt
                    p.inc(message="Interpolating")
                    rt_dimension = config.DimensionRT()+1
                    mz_dimension = config.DimensionMZ()
                    # linspace returns 'dimension' evenly spaced samples, calculated over the interval min, max
                    mz_grid = linspace(vm, vM, mz_dimension)
                    rt_grid = linspace(rm, rM, rt_dimension)
                    mz_mesh, rt_mesh = meshgrid(mz_grid, rt_grid)

                    # Interpolate using SciPy
                    #points = np.column_stack((mz_values, rt_values))
                    intensity_grid = griddata(
                        points=(mz_values, rt_values), 
                        values=intensities, 
                        xi=(mz_mesh, rt_mesh), 
                        method='cubic',  # linear, nearest, cubic
                        )
                    intensity_grid[intensity_grid < 0] = 0

                    DataCache.Store([mz_mesh, rt_mesh, intensity_grid, vm, vM, rm, rM], interpolation_cache)

                else:
                    mz_mesh, rt_mesh, intensity_grid, vm, vM, rm, rM = DataCache.Get(interpolation_cache)

                p.inc(message="Plotting")
                fig = Figure(data=[Surface(
					x=mz_mesh,
                    y=rt_mesh,
                    z=intensity_grid,
                    colorscale=config.ColorMap().lower(),
                    colorbar=dict(title='Intensity'),
                    showscale='legend' in config.Features()
				)])

                fig.update_layout(
                    title="",
                    scene=dict(
                        xaxis=dict(
                            title=dict(
                                text='m/z' if "x" in config.Features() else "",
                                font=dict(size=config.TextSize())
                            ),
                            showticklabels="x" in config.Features()
                        ),
                        yaxis=dict(
                            title=dict(
                                text='Retention Time (min)' if "y" in config.Features() else "",
                                font=dict(size=config.TextSize())
                            ),
                            showticklabels="y" in config.Features()
                        ),
                        zaxis=dict(
                            title=dict(
                                text='Intensity' if "z" in config.Features() else "",
                                font=dict(size=config.TextSize())
                            ),
                            showticklabels="z" in config.Features()
                        ),
                        camera=dict(
                            up=dict(x=0, y=0, z=1),
                            center=dict(x=0, y=0, z=0),
                            eye=dict(x=1.5, y=1.5, z=1.0)
                        )
                    ),
                    width=1000,
                    height=600,
                    margin=dict(l=0, r=0, b=0, t=40),
                    template='plotly_white',
                )

                #fig.update_layout(coloraxis_showscale='legend' in config.Features())

                # set legend if "legend" in input.Features()
                # set axis label, legend text size with config.TextSize()

                if fig is None: 
                    return CreateErrorImg("The input file could not be plotted.", "#027bc2", "text")

                p.inc(message="Exporting...")
                DataCache.Store(fig.to_html(), inputs)

        return DataCache.Get(inputs)


    # @output
    # @render.image(delete_file=True)
    # def Heatmap(): return GenerateHeatmap()


    @output
    @render.image(delete_file=True)
    def Similarity(): return GenerateSimilarity()


    # @output
    # @render.image(delete_file=True)
    # @reactive.event(input.Update)
    # def HeatmapReactive(): return GenerateHeatmap()
    @output
    @render.ui
    def HeatmapReactive(): return ui.HTML(GenerateHeatmap())

    @output
    @render.ui
    def Heatmap(): return ui.HTML(GenerateHeatmap())


    @reactive.effect
    @reactive.event(input.ExampleInfo)
    def ExampleInfo():
        Msg(ui.HTML(Info[input.Example()]))


    @render.download(filename=lambda: f"table{config.TableType()}")
    def DownloadTable(): 
        data = Data()
        
        # return error if no data to download
        if data is None:
            Error("The downloaded table is empty! Please upload your data or select an example data set in the sidebar.")
            return
        
        output_data = []
        for spectrum in data:
            if spectrum.ms_level:
                polarity = None
                if spectrum["negative scan"]:
                    polarity = "negative"
                elif spectrum["positive scan"]:
                    polarity = "positive"

                spectrum_dict = {
                    "ID": spectrum.ID,
                    "MS Level": spectrum.ms_level,
                    "RT (minutes)": spectrum.scan_time_in_minutes(),
                    "Polarity": polarity,
                    "Raw Peaks": len(spectrum.peaks("raw")),
                    "Centroided Peaks": len(spectrum.peaks("centroided")),
                    "Reprofiled Peaks": len(spectrum.peaks("reprofiled")),
                    "Highest Intensity": spectrum.highest_peaks(1)[0],
                    "Extreme Values (mz)": spectrum.extreme_values("mz"),
                }
                output_data.append(spectrum_dict)
        yield str(output_data)


    @reactive.effect
    @reactive.event(input.MainTab)
    def _():
        tab = input.MainTab()
        choices = [".png", ".jpg", ".html"]
        label="Download File Type"

        if tab == "SimilarityTab":
            choices = [".png", ".jpg"]
            selected = [".png"]
        elif tab == "HeatmapTab":
            choices = [".html"]
            selected = [".html"]
        else:
            label=""
            choices = []
            selected = []

        ui.update_radio_buttons(
            id="HeatmapType",
            label=label,
            choices=choices,
            inline=True,
            selected=selected
        )


    @render.download(filename=lambda: f"heatmap{config.HeatmapType()}")
    def DownloadHeatmap():
        data = DataCache.Get(Hash())
        if data is None:
            Error("You are trying to download an empty file! \nPlease upload your data or select an example data set in the sidebar.")
            if config.HeatmapType() == ".html":
                data = "Uhoh, no data to display! Please upload your data or select an example data set in the Heatmapper2 application."
            else:
                return
        yield data


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

    ui.panel_title(title=None, window_title="Spectral"),
    NavBar(),

    ui.layout_sidebar(
        ui.sidebar(

            FileSelection(examples={"BSA1-subset.mzML": "Example 1"}, types=[".mzml", ".mzML"], project="Spectral"),

            TableOptions(config),

            ui.panel_conditional(
                "input.MainTab != 'TableTab'",

                Update(),


                ui.HTML("<b>Heatmap</b>"),

                config.ID.UI(ui.input_select, id="ID", label="ID", selectize=True, multiple=True, choices=[0], conditional="input.MainTab === 'SimilarityTab'", tooltip="Select the IDs of the spectra whose similarity you would like to plot."),

                config.TextSize.UI(ui.input_numeric, id="TextSize", label="Text Size", min=1, max=50, step=1, tooltip="Change the text size of all axis labels. Axis labels can be toggled on and off in the 'Features' section at the bottom of this sidebar."),
                config.ColorMap.UI(
                    ui.input_select, id="ColorMap", 
                    label="Color Map", 
                    choices={
                        "viridis": "Viridis",
                        "inferno": "Inferno",
                        "plasma": "Plasma",
                        "turbo": "Rainbow",
                        }, 
                    tooltip="Select a color scheme for the heatmap."),

                config.Peaks.UI(ui.input_select, id="Peaks", label="Peak Type", choices=["Raw", "Centroided", "Reprofiled"], conditional="input.MainTab === 'HeatmapTab'", tooltip=ui.HTML('Select a peak type from the input file to display. <br>Raw visualizes unprocessed data. <br>Centroided peaks have reduced noise. <br>Reprofiled peaks have been smoothed. <br><a href="https://academic.oup.com/bioinformatics/article/28/7/1052/209917" target="_blank">Read more here</a>.')),

                config.Interpolation.UI(ui.input_select, id="Interpolation", label="Inter", choices=InterpolationMethods, conditional="input.MainTab === 'SimilarityTab'", tooltip="Specify an interpolation algorithm to apply to the similarity heat map image. This can cause values to bleed together and appear smoother."),

                config.DimensionRT.UI(ui.input_numeric, id="DimensionRT", label="Interpolate RT", conditional="input.MainTab === 'HeatmapTab'", min=1, tooltip="Specify the level of interpolation to use for generating retention time contours. Lower values generate the curve using fewer data points, which improves computation time but decreases accuracy. Higher values increase computation time, but result in smoother and more accurate curves."),
                
                config.DimensionMZ.UI(ui.input_numeric, id="DimensionMZ", label="Interpolate m/z", conditional="input.MainTab === 'HeatmapTab'", min=1, tooltip="Specify the level of interpolation for the m/z axis. Lower values generate curves across fewer m/z values, which improves computation time but decreases accuracy. Higher values increase computation time, but result in a more accurate graph."),

                config.Size.UI(ui.input_numeric, id="Size", label="Heatmap Size", min=1, conditional="input.MainTab === 'SimilarityTab'", tooltip="Change the width (in pixels) of the heatmap on your screen."),
                config.DPI.UI(ui.input_numeric, id="DPI", label="Resolution (DPI)", min=5, conditional="input.MainTab === 'SimilarityTab'", tooltip="Specify the resolution of the image in pixels per inch. Higher DPI values result in higher quality images, but larger file sizes. This setting affects the heatmap on screen as well as the downloaded plot."),


                # Customize what aspects of the heatmap are visible
                ui.HTML("<b>Features</b>"),
                config.Features.UI(
                    ui.input_checkbox_group, 
                    make_inline=False, 
                    id="Features", 
                    label=None,
                    choices={"x": "X Labels", "y": "Y Labels", "z": "Z Labels", "legend": "Legend"},
                    tooltip=ui.HTML('X and Y labels toggle the data labels along their respective axes. <br><br>Z labels toggles the data labels along the Z axis if rendering as a 3D plot. <br><br>Legend displays a colorbar legend on the heatmap.'),
                ),

                config.HeatmapType.UI(ui.input_radio_buttons, make_inline=False, id="HeatmapType", label="Download File Type", choices=[".png", ".jpg", ".html"], inline=True),
                ui.download_button(id="DownloadHeatmap", label="Download Heatmap"),
            ),
            padding="10px",
            gap="20px",
            width="300px",
        ),

        # Add the main interface tabs.
        MainTab(
            ui.nav_panel("Similarity", ui.output_plot("Similarity", height="86vh"), value="SimilarityTab"),
            m_type=ui.output_ui,
        ),
    )
)

app = App(app_ui, server)
