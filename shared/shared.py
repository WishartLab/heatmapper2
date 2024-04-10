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
		if url not in self._primary:
			r = await pyfetch(url);
			if not r.ok: return None
			else: self._primary[url] = await r.bytes()
		return self._primary[url]


	async def _local(self, url):
		if not exists(url): return None
		return Path(url)


	def _local_sync(self, url):
		if not exists(url): return None
		return Path(url)


	@staticmethod
	def GetHash(*args): return "".join([str(option) for option in args])


	def __init__(self, project, DataHandler = DefaultHandler):
		"""
		@brief Initialize an instance of the Cache object.
		@param project: The name of the project. This is used to fetch web resources.
		@param DataHandler:	The function that should be called to process files. It should
												take a name, and a binary stream, and return a DataFrame.
		"""

		# The primary cache now serves as a file agnostic cache, containing the raw bytes of files.
		self._primary = {}

		# The Secondary cache now serves as the transformed output through the handler. There is now
		# no need to specify mutability because the primary cache doesn't contain data that can be changed.
		# It serves solely as a cache for the Handler if the user throws out whatever is in the secondary.
		self._secondary = {}

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


	def SyncLoad(self, input, source_file=None, example_file=None, source=None, input_switch=None, default=DataFrame(), return_n=False):
		"""
		@brief A synchronous loading function that only supports local files
		@info See Load() for more information
		"""
		if source_file is None: source_file = input.File()
		if example_file is None: example_file = input.Example()
		if source is None: source = self._source
		if input_switch is None: input_switch = input.SourceFile()

		# Grab an uploaded file, if its done, or grab an example (Using a cache to prevent redownload)
		if input_switch == "Upload":
			file: list[FileInfo] | None = source_file
			if file is None: return (None, default) if return_n else default

			# The datapath can be immediately used to load examples, but we explicitly need to use
			# Local as a user uploaded file will always be fetched on disk.
			n = str(file[0]["datapath"])
			path = Path(n)

		# Example files, conversely, can be on disk or on a server depending on whether we're in a WASM environment.
		else:
			n = str(source + example_file)
			path = self._local_sync(n)
		
		# If the secondary cache hasn't been populated (Or was purge by the user), populate it.
		if n not in self._secondary:
			self._secondary[n] = self._handler(path)

		return (n, self._secondary[n]) if return_n else self._secondary[n]


	async def Load(self, input, source_file=None, example_file=None, source=None, input_switch=None, default=DataFrame(), return_n=False):
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
			if file is None: return (None, default) if return_n else default

			# The datapath can be immediately used to load examples, but we explicitly need to use
			# Local as a user uploaded file will always be fetched on disk.
			n = str(file[0]["datapath"])
			path = Path(n)

		# Example files, conversely, can be on disk or on a server depending on whether we're in a WASM environment.
		else:
			n = str(source + example_file)
			raw = await self._download(n)
			if type(raw) is bytes:
				temp = NamedTemporaryFile(suffix=Path(n).suffix); temp.write(BytesIO(raw).read()); temp.seek(0)
				path = Path(temp.name)
			else: path = raw

		# If the secondary cache hasn't been populated (Or was purge by the user), populate it.
		if n not in self._secondary:
			self._secondary[n] = self._handler(path)

		return (n, self._secondary[n]) if return_n else self._secondary[n]


	async def Update(self, input):
		"""
		@brief Updates information within the secondary cache based on user selection
		@param input: The Shiny input. Importantly, these must be defined:
			input.TableRow: The row to modify
			input.TableCol: The column to modify
			input.TableVal: What the user wants to set as the new value
		@info This function should be called on a reactive hook for a "Update" button.
		"""

		# Get the data
		df = await self.Load(input)
		if not type(df) is DataFrame: return

		row_count, column_count = df.shape
		row, column = input.TableRow(), input.TableCol()

		# So long as row and column are sane, update.
		if row < row_count and column < column_count:
			try:
				if input.Type() == "Integer": df.iloc[row, column] = int(input.TableVal())
				elif input.Type() == "Float": df.iloc[row, column] = float(input.TableVal())
				else: df.iloc[row, column] = input.TableVal()
			except ValueError: pass


	async def Purge(self, input, source_file=None, example_file=None, source=None):
		"""
		@brief Purges the secondary cache of whatever the user has uploaded/selected
		@param input: The Shiny input. See N() for required objects.
		@param source_file: The source ID, defaults to input.File()
		@param example_file: The example ID, defaults to input.Example()
		@param source: The path that should be appending to the path for fetching. 
		@info This function should be called on a reactive hook for a "Reset" button.
		"""

		if source_file is None: source_file = input.File()
		if example_file is None: example_file = input.Example()
		if source is None: source = self._source

		if input.SourceFile() == "Upload":
			file: list[FileInfo] | None = source_file
			if file is None: return None
			n = file[0]["datapath"]
		else: n = source + example_file
		del self._secondary[n]


	def Hash(self, hash, data):
		"""
		@brief Stores an object with a group of options
		@param hash: The hash created by Cache.GetHash()
		@param data: The object that should be stored in the cache.
		"""
		self._secondary[hash] = data


	def Hashed(self, hash):
		"""
		@brief Checks if a particular set of options are cached.
		@param hash: The hash created by Cache.GetHash()
		@returns The object, if it exists, None otherwise.
		"""
		if hash in self._secondary:
			return self._secondary[hash]
		else:
			return None



def TableValueUpdate(df, input):
	"""
	@brief Updates the value displayed in the TableVal based on the current selection
	@param df The DataFrame
	@param input The shiny input
	"""

	if not df.empty:
		rows, columns = df.shape
		try:
			row, column = int(input.TableRow()), int(input.TableCol())
			if 0 <= row <= rows and 0 <= column <= columns:
				ui.update_text(id="TableVal", label="Value (" + str(df.iloc[row, column]) + ")")
		except TypeError: pass


def NavBar(current):
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
				ui.nav_panel(ui.HTML(f'<a href="{Sources["expression"]}" target="_blank" rel="noopener noreferrer">Expression</a>'), value="Expression"),
				ui.nav_panel(ui.HTML(f'<a href="{Sources["pairwise"]}" target="_blank" rel="noopener noreferrer">Pairwise</a>'), value="Pairwise"),
				ui.nav_panel(ui.HTML(f'<a href="{Sources["image"]}" target="_blank" rel="noopener noreferrer">Image</a>'), value="Image"),
				ui.nav_panel(ui.HTML(f'<a href="{Sources["geomap"]}" target="_blank" rel="noopener noreferrer">Geomap</a>'), value="Geomap"),
				ui.nav_panel(ui.HTML(f'<a href="{Sources["geocoordinate"]}" target="_blank" rel="noopener noreferrer">Geocoordinate</a>'), value="Geocoordinate"),
				ui.nav_panel(ui.HTML(f'<a href="{Sources["3d"]}" target="_blank" rel="noopener noreferrer">3D</a>'), value="3D"),
				ui.nav_panel(ui.HTML(f'<a href="{Sources["spatial"]}" target="_blank" rel="noopener noreferrer">Spatial</a>'), value="Spatial"),
				ui.nav_panel(ui.HTML('<a href=https://github.com/WishartLab/heatmapper2/wiki target="_blank" rel="noopener noreferrer">About</a>'), value="About"),
				title="Heatmapper",
				selected=current,
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
	return [ui.HTML('<a href=https://github.com/WishartLab/heatmapper2/wiki/Format target="_blank" rel="noopener noreferrer">Data Format</a>'),

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
			col_widths=[8,2],
		)
	),

	ui.br(),
	ui.HTML(f"<a href='https://github.com/WishartLab/heatmapper2/wiki/Interface#{project}' target='_blank' rel='noopener noreferrer'>Info on Settings</a>"),
	]


def MainTab(*args, m_type=ui.output_plot):
	return ui.navset_tab(
		ui.nav_panel("Interactive", m_type("Heatmap", height="75vh"), value="Interactive"),
		ui.nav_panel("Table",
			ui.layout_columns(
				ui.input_numeric("TableRow", "Row", 0, min=0),
				ui.input_numeric("TableCol", "Column", 0, min=0),
				ui.input_text("TableVal", "Value", 0),
				ui.input_select(id="Type", label="Datatype", choices=["Integer", "Float", "String"]),
				col_widths=[2,2,6,2],
			),
			ui.layout_columns(
				ui.input_action_button("Update", "Update"),
				ui.input_action_button("Reset", "Reset Values"),
			),
			ui.output_data_frame("LoadedTable")
		),
		*args,
		id="MainTab"
	)
