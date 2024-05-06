#
# Heatmapper
# Shared
#
# This file contains shared functionality between Heatmapper applications. It is not a standalone application.
# Due to the way ShinyLive exports applications, this file is symlinked into each project to reduce redundancy.
#

from shiny import ui
from shiny.types import FileInfo
from pandas import DataFrame, read_csv, read_excel, read_table
from io import BytesIO
from tempfile import NamedTemporaryFile
from sys import modules
from pathlib import Path
from enum import Enum
from os.path import exists
from copy import deepcopy

# Used for fetching web resources in a variety of fashions.
URL = "https://wishartlab.github.io/heatmapper2"
Raw = "https://raw.githubusercontent.com/WishartLab/heatmapper2/main"

# Define the Server and Port of the Shiny instances (Port is incremented)
# Change these if Heatmapper is running on a server.
Server = "http://35.208.86.138"
Port = 8000

# Detect the running environment
if "pyodide" in modules:
	from pyodide.http import pyfetch
	Pyodide = True
else:
	Pyodide = False

# MatPlotLib Colors
Colors = ["Blue", "Orange", "Green", "Red", "Purple", "Brown", "Pink", "Gray", "Olive", "Cyan", "White", "Yellow"]

class ColumnType(Enum): Time = 0; Name = 1; Value = 2; Longitude = 3; Latitude = 4; X = 5; Y = 6; Z = 7; Cluster = 8; Free = 9; Spatial = 10;
Columns = {
	ColumnType.Time: {"time", "date", "year"},
	ColumnType.Name: {"name", "orf", "uniqid", "face", "triangle"},
	ColumnType.Value: {"value", "weight", "intensity", "in_tissue"},
	ColumnType.Longitude: {"longitude", "long"},
	ColumnType.Latitude: {"latitude", "lat"},
	ColumnType.X: {"x"},
	ColumnType.Y: {"y"},
	ColumnType.Z: {"z"},
	ColumnType.Cluster: {"cell type", "celltype_mapped_refined", "cluster", "cell_class", "cell_subclass", "cell_cluster"},
	ColumnType.Free: {None},
	ColumnType.Spatial: {"spatial"}
}


def Filter(columns, ctype: ColumnType, good: list = [], bad: list = [], only_one=False, ui_element=None, reject_unknown=False):
	"""
	@brief Filters available column names based on what input we want
	@param columns: The columns of the DataFrame (Usually just df.columns)
	@param ctype: The type of column we're looking for (Look at the ColumnType Enum)
	@param good: A list of column names on top of those defined by the type to be included
	@param bad: A list of column names on top of those defined by the type to be excluded from the result.
	@param only_one: Only return a single result, so the variable can be used immediately.
	@param ui_element: An optional Shiny selection input to update.
	@param reject_unknown: Only include columns explicitly defined
	@return: A list of column names to use.
	"""

	# Fold cases
	folded = [column.lower() for column in columns]

	# Add and remove what user asked for, filtering None
	options = set(folded)
	if bad: options -= set([b.lower() for b in bad if b])
	if good: options &= set([g.lower() for g in good if g])

	# Take an intersection of our columns and the type we want. If there is a match, return those
	# Otherwise, remove all columns we know it shouldn't be, and return that instead.
	intersection = options & Columns[ctype]
	if intersection or reject_unknown: options = intersection
	else:
		for key, value in Columns.items():
			if key != ctype: options -= value

	# Get the valid indices, and sort them in ascending order
	indices = [folded.index(value) for value in options]
	indices.sort()

	# Get the original column names, without case-folding, and return as a list.
	reassembled = [columns[index] for index in indices]
	if not reassembled: return None

	# Update a UI element, if one was provided
	if ui_element is not None: ui.update_select(id=ui_element, choices=reassembled, selected=reassembled[0])
	return reassembled[0] if only_one else reassembled


class Cache:
	"""
	@brief A class that encompasses fetching/storing web resources.
	"""

	@staticmethod
	def HandleDataFrame(path, function):
		"""
		@brief Handle DataFrame's
		@param i: The binary of the file
		@param function: The pandas function to use to read the file.
		@returns A DataFrame
		"""

		# Read the table once.
		df = function(path.resolve()).fillna(0)

		# If the first column value is a float, we assume it's data, and not column names.
		# Re-read the DataFrame with generic column names instead
		try:
			float(df.columns[0])
			df = function(path.resolve(), header=None, names=[f"Column {i}" for i in range(df.shape[1])])
		except ValueError: pass
		return df


	@staticmethod
	def DefaultHandler(path):
		"""
		@brief The default handler. It can handle csv, xlsx, and defaults all other files to read_table
		@param n: The path to the file
		@returns: An object, if the provided file is supported, None otherwise.
		"""

		suffix = path.suffix
		if suffix == ".csv": return Cache.HandleDataFrame(path, read_csv)
		elif suffix == ".xlsx" or suffix == ".xls" or suffix == ".odf": return Cache.HandleDataFrame(path, read_excel)
		else: return Cache.HandleDataFrame(path, read_table)


	async def _remote(self, url):
			r = await pyfetch(url);
			if not r.ok: return None
			return await r.bytes()


	async def _local(self, url):
		if not exists(url): return None
		return Path(url)


	def __init__(self, project, DataHandler = DefaultHandler):
		"""
		@brief Initialize an instance of the Cache object.
		@param project: The name of the project. This is used to fetch web resources.
		@param DataHandler:	The function that should be called to process files. It should
												take a name, and a binary stream, and return a DataFrame.
		"""

		# The primary is the unprocessed, fetched web resources
		self._primary = {}

		# The objects are anything that applications want to store
		self._objects = {}

		# The data handler for processing the binary files.
		self._handler = DataHandler

		# If we're in a Pyodide environment, we fetch resources from the web.
		if Pyodide:
			self._download = lambda url: self._remote(url)
			self._source = f"{Raw}/{project}/example_input/"

		# Otherwise, we fetch locally.
		else:
			self._download = lambda url: self._local(url)
			self._source = "../example_input/"


	async def Load(self, input, source_file=None, example_file=None, source=None, input_switch=None, default=DataFrame()):
		"""
		@brief Caches whatever the user has currently uploaded/selection, returning the identifier within the secondary cache.
		@param input: The Shiny input variable. Importantly, these must be defined:
			input.File: The uploaded file
			input.Example: The selected example file
			input.SourceFile: Whether the user wants "Upload" or "Example"
		@param source_file: The input ID that should be used to fetch the file (Defaults to input.File() if None)
		@param example_file: The input ID that should be used to fetch th example (Defaults to input.Example() if None)
		@param input_switch:	The input ID to check for Upload/Example/Other. The value is compared against "Upload" for user
													uploaded items, and defaults to fetching example_file otherwise. (Defaults to input.SourceFile())
		@param default:	The object that should be returned if files cannot be fetched. Ensures that Load will always return an
										object, avoiding the needing to check output. Defaults to a DataFrame. The object should be able to
										initialize without arguments.
		@param return_n: Return the filename for post-processing.
		"""

		if source_file is None: source_file = input.File()
		if example_file is None: example_file = input.Example()
		if source is None: source = self._source
		if input_switch is None: input_switch = input.SourceFile()

		# Grab an uploaded file, if its done, or grab an example (Using a cache to prevent redownload)
		if input_switch == "Upload":
			file: list[FileInfo] | None = source_file
			if file is None: return default

			# The datapath can be immediately used to load examples, but we explicitly need to use
			# Local as a user uploaded file will always be fetched on disk.
			n = str(file[0]["datapath"])
			path = Path(n)

		# Example files, conversely, can be on disk or on a server depending on whether we're in a WASM environment.
		else:
			n = str(source + example_file)
			raw = await self._download(n)

			# WASM needs a temporary file, but they are deleted out of their scope.
			if Pyodide:
				temp = NamedTemporaryFile(suffix=Path(n).suffix); 
				temp.write(BytesIO(raw).read()); temp.seek(0)
				if n not in self._primary: self._primary[n] = self._handler(Path(temp.name))
			elif n not in self._primary: self._primary[n] = self._handler(raw)

		# If the object cannot be copied, then we can just return it directly
		try:
			return deepcopy(self._primary[n])
		except AttributeError:
			return self._primary[n]


	def Store(self, object, inputs):
		"""
		@brief Store arbitrary data in the Cache.
		@param object: The object to store
		@param inputs: A list of values that compose a hash of the object.
		"""
		h = "".join(str(i) for i in inputs)
		self._objects[h] = object


	def Get(self, inputs):
		"""
		@brief Retrieve arbitrary data in the Cache.
		@param inputs: A list of values that compose a hash of the object.
		"""
		h = "".join(str(i) for i in inputs)
		if h in self._objects:
				return self._objects[h]
		else: return None

	def In(inputs):
		h = "".join(str(i) for i in inputs)
		return h in self._objects


def NavBar():
	"""
	@brief Returns a Navigation Bar for each project, with the current project selected.
	@returns A list, containing a ui.panel_title, and a ui.navset_bar.
	"""

	Sources = {
		"expression": f"{URL}/expression/site/index.html" if Pyodide else f"{Server}:{Port}",
		"pairwise": f"{URL}/pairwise/site/index.html" if Pyodide else f"{Server}:{Port + 1}",
		"image": f"{URL}/image/site/index.html" if Pyodide else f"{Server}:{Port + 2}",
		"geomap": f"{URL}/geomap/site/index.html" if Pyodide else f"{Server}:{Port + 3}",
		"geocoordinate": f"{URL}/geocoordinate/site/index.html" if Pyodide else f"{Server}:{Port + 4}",
		"3d": f"{Server}:{Port + 5}",
		"spatial": f"{Server}:{Port + 6}",
	}

	return (
		ui.panel_title(title=None, window_title="Heatmapper"),
		ui.navset_bar(
			ui.nav_control(ui.HTML(f'<a href="{Sources["expression"]}" target="_blank" rel="noopener noreferrer">Expression</a>')),
			ui.nav_control(ui.HTML(f'<a href="{Sources["pairwise"]}" target="_blank" rel="noopener noreferrer">Pairwise</a>')),
			ui.nav_control(ui.HTML(f'<a href="{Sources["image"]}" target="_blank" rel="noopener noreferrer">Image</a>')),
			ui.nav_control(ui.HTML(f'<a href="{Sources["geomap"]}" target="_blank" rel="noopener noreferrer">Geomap</a>')),
			ui.nav_control(ui.HTML(f'<a href="{Sources["geocoordinate"]}" target="_blank" rel="noopener noreferrer">Geocoordinate</a>')),
			ui.nav_control(ui.HTML(f'<a href="{Sources["3d"]}" target="_blank" rel="noopener noreferrer">3D</a>')),
			ui.nav_control(ui.HTML(f'<a href="{Sources["spatial"]}" target="_blank" rel="noopener noreferrer">Spatial</a>')),
			ui.nav_control(ui.HTML('<a href=https://github.com/WishartLab/heatmapper2/wiki target="_blank" rel="noopener noreferrer">About</a>')),
			ui.nav_spacer(),
			ui.nav_control(ui.input_dark_mode(id="mode")),
			title="Heatmapper",
		),
	)


def FileSelection(examples, types, upload_label="Choose a File", multiple=False, default="Example", project="Overview"):
	"""
	@brief Returns the file selection dialog for the user to upload/select an example
	@param examples: Either a list of example file names, or a dictionary mapping
	@param types: The valid file extensions for user uploaded files.
	@param upload_label: The label for the upload input. Useful to define specifically what kind of files are needed
	@param multiple: Whether to accept multiple files. 
	@param default: Whether to start on the example, or upload dialog
	@param project: The name of a project, to specify a specified header within the Interface documentation
	@returns A list, containing the necessary ui components for uploading/selecting
	@info The returns elements are named:
		input.SourceFile: The ui.input_radio_buttons for whether the user wants to choose an "Example" or "Upload"
		input.File: The ui.input_file for user uploaded files.
		input.Example: The ui.input_select for an example file selection
	@info multiple=True is not handled properly by the Cache. You will need to create a function that properly handles
		each file (See spatial for an implementation)
	@info If you're examples are large files, or require significant computation, you may want to switch it to Upload instead.
	"""

	# If the user needs help with the formatting.
	return [
	ui.layout_columns(
		ui.HTML("<a href=https://github.com/WishartLab/heatmapper2/wiki/Format target='_blank' rel='noopener noreferrer'>Format</a>"),
		ui.HTML(f"<a href='https://github.com/WishartLab/heatmapper2/wiki/Interface#{project}' target='_blank' rel='noopener noreferrer'>Help</a>"),
		col_widths=[6,6]
	),

	# Specify whether to use example files, or upload one.
	ui.input_radio_buttons(id="SourceFile", label="Specify a Source File", choices=["Example", "Upload"], selected=default, inline=True),

	# Only display an input dialog if the user is one Upload
	ui.panel_conditional(
		"input.SourceFile === 'Upload'",
		ui.input_file("File", upload_label, accept=types, multiple=multiple),
	),

	# Otherwise, add the example selection and an info button.
	ui.panel_conditional(
		"input.SourceFile === 'Example'",
		ui.layout_columns(
			ui.input_select(id="Example", label=None, choices=examples, multiple=False),
			ui.popover(ui.input_action_link(id="ExampleInfoButton", label="Info"), ui.output_text("ExampleInfo")),
			col_widths=[7,3],
		)
	),
	]


def TableOptions():
	"""
	@brief Return the options for Table Manipulation.
	@returns A conditional panel that provides a DataType, and a ResetButton.
	"""
	return  ui.panel_conditional(
		"input.MainTab === 'TableTab'",
		ui.input_radio_buttons(id="Type", label="Datatype", choices=["Integer", "Float", "String"], inline=True),
		ui.input_action_button(id="Reset", label="Reset Values"),
	),


def ColorMap():
	"""
	@brief Returns a ColorMap input Selection
	@returns	A list of UI elements. Firstly, a header that contained a toggle for custom ColorMaps.
						Then, two conditional panels based on the status of the toggle. If the user wants
						custom color maps, then provide a selectize.js selection box that allows for multiple
						selections. These options are used to define the gradient, low to high. If not, just a
						collection of predefined maps, where the keys must be split on spaces to generate a map
						that MatPlotLib can use.
	"""
	return [
		ui.layout_columns("Color Map", ui.input_checkbox(id="Custom", label="Custom")),
		ui.panel_conditional(
			"input.Custom",
				ui.input_select(
				id="CustomColors",
				label=None,
				choices=Colors,
				selected=["Blue", "White", "Yellow"],
				multiple=True,
				selectize=True,
			),
		),
		ui.panel_conditional(
			"!input.Custom",
			ui.input_select(id="ColorMap", label=None, choices={
					"Blue White Yellow": "Blue/Yellow",
					"Red Black Green": "Red/Green",
					"Pink White Green": "Pink/Green",
					"Blue Green Yellow": "Blue/Green/Yellow",
					"Black Gray White": "Grayscale",
					"Red Orange Yellow Green Blue Indigo Violet": "Rainbow",
				}
			),
		)]


def MainTab(*args, m_type=ui.output_plot):
	return ui.navset_tab(
		ui.nav_panel("Heatmap", m_type(id="Heatmap", height="90vh"), value="HeatmapTab"),
		ui.nav_panel("Table", ui.output_data_frame(id="Table"), value="TableTab"),
		*args,
		id="MainTab"
	)
