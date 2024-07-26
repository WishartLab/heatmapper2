#
# Heatmapper
# Shared
#
# This file contains shared functionality between Heatmapper applications. It is not a standalone application.
# Due to the way ShinyLive exports applications, this file is symlinked into each project to reduce redundancy.
#

from shiny import ui, reactive
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
Server = "http://heatmapper2.ca"
Port = 8000

# Detect the running environment
if "pyodide" in modules:
	from pyodide.http import pyfetch
	Pyodide = True
	async def fetch(url):
		response = await pyfetch(url)
		if response.ok: return (await response.bytes())
		else:
			Error("Could not download file")
			return None
else:
	from urllib.request import urlopen
	Pyodide = False
	async def fetch(url):
		try: return urlopen(url).read()
		except Exception:
			Error("Could not download file!")
			return None

# Shared Values
Colors = ["Blue", "Orange", "Green", "Red", "Purple", "Brown", "Pink", "Gray", "Olive", "Cyan", "White", "Yellow"]
DistanceMethods = ["Braycurtis", "Canberra", "Chebyshev", "Cityblock", "Correlation", "Cosine", "Dice", "Euclidean", "Hamming", "Jaccard", "Jensenshannon", "Kulczynski1", "Matching", "Minkowski", "Rogerstanimoto", "Russellrao", "Seuclidean", "Sokalmichener", "Sokalsneath", "Sqeuclidean", "Yule"]
InterpolationMethods = ["None", "Antialiased", "Nearest", "Bilinear", "Bicubic", "Spline16", "Spline36", "Hanning", "Hamming", "Hermite", "Kaiser", "Quadric", "Catrom", "Gaussian", "Bessel", "Mitchell", "Sinc", "Lanczos", "Blackman"]
ClusteringMethods = ["Single", "Complete", "Average", "Weighted", "Centroid", "Median", "Ward"]
ColorMaps = ["Viridis", "Plasma", "Inferno", "Magma"]

class ColumnType(Enum): Time = 0; Name = 1; Value = 2; Longitude = 3; Latitude = 4; X = 5; Y = 6; Z = 7; Cluster = 8; Free = 9; Spatial = 10; NameGeoJSON = 11; FOV = 12; Count = 13
Columns = {
	ColumnType.Time: {"time", "date", "year"},
	ColumnType.Name: {"name", "orf", "uniqid", "face", "triangle", "iso_code", "continent", "country", "location"},
	ColumnType.Value: {"value", "weight", "intensity", "in_tissue"},
	ColumnType.Longitude: {"longitude", "long", "lon"},
	ColumnType.Latitude: {"latitude", "lat"},
	ColumnType.X: {"x"},
	ColumnType.Y: {"y"},
	ColumnType.Z: {"z"},
	ColumnType.Cluster: {"cluster"},
	ColumnType.Free: {None},
	ColumnType.Spatial: {"spatial"},
	ColumnType.NameGeoJSON: {"name", "admin", "iso_a3", "iso_a2", "iso"},
	ColumnType.FOV: {"fov"},
	ColumnType.Count: {'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 'n_counts'}
	}


def Filter(columns, ctype: ColumnType, good: list = [], id=None, all=False, remove_unknown=False):
	"""
	@brief Filters available column names based on what input we want
	@param columns: The columns of the DataFrame (Usually just df.columns)
	@param ctype: The type of column we're looking for (Look at the ColumnType Enum)
	@param good: A list of column names on top of those defined by the type to be included
	@param id: An element id to update with a new value.
	@param all: Return all matches columns
	@return: A list of column names to use.

	@info This purpose of this function is to try and remove irrelevant columns from user selection,
	but returning everything if by filtering so we remove all the columns. In essence, it folds the case
	of all columns, and performs a set intersection on the required column type. This set is then returned
	to the case of the original columns, and then good and bad are applied (Therefore, they are case-sensitive)
	Since both good and bad are applied after the intersection, they don't need to be valid names (So long as)
	the application can handle that exception. Look at Geocoordinate to see how it uses a "None" and "Uniform"
	value in the good list, despite these values both not a valid ValueColumn, and not existing in the data.

	The logic for the UI updating can be confusing, but in essence we don't just want to return the good
	list, because that means we removed all actual columns. If this happens, we return all the columns, and
	add the good list onto the START (So its the default), that way users can choose a column if Heatmapper
	doesn't like their column names.
	"""

	# Fold cases
	folded = []
	for column in columns:
		try:
			folded.append(column.lower())
		except Exception:
			folded.append(column)
	options = set(folded)
	if ctype != ColumnType.Free: options &= Columns[ctype]

	if remove_unknown:
		for type in Columns:
			if type != ctype: options -= Columns[type]

	indices = [folded.index(value) for value in options]; indices.sort()
	reassembled = [columns[index] for index in indices] + good

	if id:
		if reassembled == good:
			options = set(folded)
			for type in Columns:
				if type != ctype: options -= Columns[type]
			indices = [folded.index(value) for value in options]; indices.sort()
			reassembled = good + [columns[index] for index in indices]
		ui.update_select(id=id, choices=reassembled)
	if all: return reassembled
	return reassembled[0] if reassembled and len(reassembled) > 0 else None


class Cache:
	"""
	@brief A class that encompasses fetching/storing web resources.
	"""

	@staticmethod
	def HandleDataFrame(path, function, p=None):
		"""
		@brief Handle DataFrame's
		@param i: The binary of the file
		@param function: The pandas function to use to read the file.
		@returns A DataFrame
		"""

		# Read the table once.
		if p: p.inc(message="Reading Table...")
		df = function(path.resolve()).fillna(0)

		# If the first column value is a float, we assume it's data, and not column names.
		# Re-read the DataFrame with generic column names instead
		try:
			float(df.columns[0])
			if p: p.inc(message="Generating incidices...")
			df = function(path.resolve(), header=None, names=[f"Column {i}" for i in range(df.shape[1])])
		except ValueError: pass
		return df


	@staticmethod
	def DefaultHandler(path, p=None):
		"""
		@brief The default handler. It can handle CSVs, Excel files, Tables, and all other files will simply
		be stored as strings of the file content

		@param n: The path to the file
		@returns: An object.
		"""

		suffix = path.suffix
		if suffix == ".csv": return Cache.HandleDataFrame(path, read_csv, p)
		elif suffix in {".xlsx", ".xls", ".odf"}: return Cache.HandleDataFrame(path, read_excel, p)
		elif suffix in {".txt", ".dat", ".tsv", ".tab"}: return Cache.HandleDataFrame(path, read_table, p)
		else: return open(path.resolve(), "r").read()


	async def _local(self, url):
		if not exists(url): return None
		return open(url, "rb").read()


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
			self._download = fetch
			self._source = f"{Raw}/{project}/example_input/"

		# Otherwise, we fetch locally.
		else:
			self._download = self._local
			self._source = "../example_input/"


	async def Download(self, n, p=None):
		"""
		@brief Downloads any arbitrary URL and stores it in the cache
		@param n: The URL name
		@returns The handled data

		"""
		if n not in self._primary:
			raw = await (fetch(n) if n.startswith("https://") else self._download(n))
			if raw is None: return None

			path = Path(n)
			if path.is_file():
				self._primary[n] = self._handler(path, p)
			else:
				temp = NamedTemporaryFile(suffix=Path(n).suffix);
				temp.write(raw);
				temp.seek(0)
				self._primary[n] = self._handler(Path(temp.name), p)
		try:
			return deepcopy(self._primary[n])
		except Exception:
			return self._primary[n]


	async def Load(self,
		input,
		source_file=None,
		example_file=None,
		source=None,
		input_switch=None,
		upload="Upload",
		example="Example",
		default=DataFrame(),
		p=None,
		p_name="file",
		wasm=True,
		wasm_blacklist=tuple()
		):

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
		@param upload: The value of the input_switch such that we should fetch a source file from source_file
		@param example: The value of the input_switch such that we should fetch an example from example_file
		@param default:	The object that should be returned if files cannot be fetched. Ensures that Load will always return an
										object, avoiding the needing to check output. Defaults to a DataFrame. The object should be able to
										initialize without arguments.
		@param p: A progress bar to increment; optional.
		@param p_name: What we're fetching, to be displayed in the progress bar; optional
		@param wasm: Whether this fetch can run in WebAssembly
		@param wasm_blacklist: A tuple of file extensions that should not be fetched in WebAssembly.
		@returns The data if it exists; default if no file can be found; 0 if there's a WebAssembly violation
		"""

		if not wasm and Pyodide:
			if p: p.close()
			return 0

		if source_file is None: source_file = input.File()
		if example_file is None: example_file = input.Example()
		if source is None: source = self._source
		if input_switch is None: input_switch = input.SourceFile()

		# Grab an uploaded file, if its done, or grab an example (Using a cache to prevent redownload)
		if input_switch == upload:
			if p: p.inc(message=f"Loading Uploaded {p_name}...")
			file: list[FileInfo] | None = source_file
			if file is None:
				if p: p.close()
				return default

			if p: p.inc(message=f"Handling {p_name}...")
			# The datapath can be immediately used to load examples, but we explicitly need to use
			# Local as a user uploaded file will always be fetched on disk.
			n = str(file[0]["datapath"])
			if n.endswith(wasm_blacklist) and Pyodide:
				if p: p.close()
				return 0
			if n not in self._primary: self._primary[n] = self._handler(Path(n), p)

		# Example files, conversely, can be on disk or on a server depending on whether we're in a WASM environment.
		elif input_switch == example:
			if p: p.inc(message=f"Fetching {p_name}...")

			# If we explicitly provide a URL, use it, but only in Pyodide (We still assume the file exists on disk when running
			# in server-mode).
			if example_file.startswith("https://"):
				n = example_file if Pyodide else str(source + example_file.split("/")[-1])
			else:
				n = str(source + example_file)

			if n.endswith(wasm_blacklist) and Pyodide:
				if p: p.close()
				return 0

			if p: p.close()
			return await self.Download(n)

		# If the application has a unique method of input (IE 3D's ID, don't handle it.)
		else:
			if p: p.close()
			return None

		if p: p.close()
		# If the object cannot be copied, then we can just return it directly
		try:
			return deepcopy(self._primary[n])
		except Exception:
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

	def In(self, inputs):
		h = "".join(str(i) for i in inputs)
		return h in self._objects


	def Invalidate(self, input):
		invalid = []
		for key, value in self._objects.items():
			if input in key:
				invalid.append(key)
		for i in invalid: del self._objects[i]


def NavBar():
	"""
	@brief Returns a Navigation Bar for each project, with the current project selected.
	@returns A list, containing a ui.panel_title, and a ui.navset_bar.
	"""

	Sources = {
		"expression": f"{URL}/site/expression/index.html" if Pyodide else f"{Server}:{Port}",
		"pairwise": f"{URL}/site/pairwise/index.html" if Pyodide else f"{Server}:{Port + 1}",
		"image": f"{URL}/site/image/index.html" if Pyodide else f"{Server}:{Port + 2}",
		"geomap": f"{URL}/site/geomap/index.html" if Pyodide else f"{Server}:{Port + 3}",
		"geocoordinate": f"{URL}/site/geocoordinate/index.html" if Pyodide else f"{Server}:{Port + 4}",
		"3d": f"{URL}/site/3d/index.html" if Pyodide else f"{Server}:{Port + 5}",
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
			title=ui.HTML(
				f'<a href="{URL}" target="_blank" rel="noopener noreferrer"> \
					<img src="{Raw}/site/logo.png" alt="Heatmapper"> \
				</a>'),
		),
	)


def FileSelection(examples, types, upload_label=None, multiple=False, default="Upload", project="Overview", extras=[]):
	"""
	@brief Returns the file selection dialog for the user to upload/select an example
	@param examples: Either a list of example file names, or a dictionary mapping
	@param types: The valid file extensions for user uploaded files.
	@param upload_label: The label for the upload input. Useful to define specifically what kind of files are needed
	@param multiple: Whether to accept multiple files.
	@param default: Whether to start on the example, or upload dialog
	@param project: The name of a project, to specify a specified header within the Interface documentation
	@param extras: Extra options for giving the application information no render. You are responsible for handling it.
	@returns A list, containing the necessary ui components for uploading/selecting
	@info The returns elements are named:
		input.SourceFile: The ui.input_radio_buttons for whether the user wants to choose an "Example" or "Upload"
		input.File: The ui.input_file for user uploaded files.
		input.Example: The ui.input_select for an example file selection
	@info multiple=True is not handled properly by the Cache. You will need to create a function that properly handles
		each file (See spatial for an implementation)
	@info If you're examples are large files, or require significant computation, you may want to switch it to Upload instead.
	"""

	return [
		ui.HTML(f"<a href='https://github.com/WishartLab/heatmapper2/wiki/Interface#{project}' target='_blank' rel='noopener noreferrer'>Help</a>"),

		ui.input_radio_buttons(
			id="SourceFile",
			label=ui.HTML("Specify a File (<a href=https://github.com/WishartLab/heatmapper2/wiki/Format target='_blank' rel='noopener noreferrer'>Format</a>)"),
			choices=["Example", "Upload"] + extras,
			selected=default,
			inline=True
	),

		# Only display an input dialog if the user is one Upload
		ui.panel_conditional(
			"input.SourceFile === 'Upload'",
			ui.input_file(id="File", label=None, accept=types, multiple=multiple),
		),
		ui.panel_conditional(
			"input.SourceFile === 'Example'",
			Inlineify(ui.input_select, id="Example", label=ui.input_action_link(id="ExampleInfo", label="Example"), choices=examples),
		),
	]


def TableOptions(config):
	"""
	@brief Return the options for Table Manipulation.
	@returns A conditional panel that provides a DataType, and a ResetButton.
	"""
	return  ui.panel_conditional(
		"input.MainTab === 'TableTab'",
		config.Type.UI(ui.input_radio_buttons, make_inline=False, id="Type", label="Datatype", choices=["Integer", "Float", "String"], inline=True),
		ui.input_action_button(id="Reset", label="Reset Values"),
		ui.download_button(id="DownloadTable", label="Download Table"),
	),


def MainTab(*args, m_type=ui.output_plot):
	return ui.navset_tab(
		ui.nav_panel("Heatmap",
			ui.panel_conditional("input.UpdateToggle", m_type(id="Heatmap")),
			ui.panel_conditional("!input.UpdateToggle", m_type(id="HeatmapReactive")),
			value="HeatmapTab"
		),
		ui.nav_panel("Table", ui.output_data_frame(id="Table"), value="TableTab"),
		*args,
		id="MainTab"
	)


def Inlineify(ui_element, widths=[4,8], gap="20px", **kwargs):
	label = kwargs["label"]
	kwargs["label"] = None
	return ui.layout_columns(
		label,
		ui_element(**kwargs),
		col_widths=widths,
		gap=gap,
	)

class Config:
	"""
	@brief A configuration entry.
	"""

	def __init__(self, visible=True, **kwargs):
		"""
		@brief Create a configuration entry.
		@param default: The default value for an input.
		@param visible: Whether the input should be shown in the sidebar
		@param **kwargs: Arguments to be passed to the input.
		"""
		self.visible = visible
		self.kwargs = kwargs
		if "selected" in kwargs:
			self.default = kwargs["selected"]
		elif "value" in kwargs:
			self.default = kwargs["value"]
		else:
			self.default = None
		self.resolve = None

	def __call__(self):
		try:
			resolved = self.resolve()
			return self.default if resolved is None else resolved
		except Exception:
			return self.default


	def Resolve(self, input):
		self.resolve = input


	def UI(self, ui_element, make_inline=True, widths=[4,8], gap="20px", conditional=None, *args, **kwargs):
		"""
		@brief Displays the configured UI.
		@param ui The Shiny interface element to use.
		@parram **kwargs: Additional arguments to be passed to the input.
		@note	keyword arguments passed to the Config object during initialization will overrule
					arguments passed to this function. Duplicates are allowed.
		"""

		combined = self.kwargs

		for key in kwargs.keys():
			combined[key] = kwargs[key]

		if "selected" in combined: combined["selected"] = self()
		elif "value" in combined: combined["value"] = self()

		if self.visible:
			if make_inline and "label" in combined:
				element = Inlineify(ui_element, widths, gap, **combined)
			else: element = ui_element(*args, **combined)

			# There doesn't seem any good way to remove the conditional panel spacing.
			# Rather than having conditional configurations stick out due to inconsistent spacing
			# Just give all the panels the same spacing by putting them in true conditions.
			return ui.panel_conditional("1 === 1" if conditional is None else conditional, element, _add_ws=False)


class ConfigHandler(dict):
	"""
	@brief: A dictionary that can be accessed with dots, and can automatically resolve.
	"""

	__getattr__ = dict.get
	__setattr__ = dict.__setitem__
	__delattr__ = dict.__delitem__


	def Resolve(self, input):
		"""
		@brief Resolves all stored objects.
		@param input The input to use for resolving.
		"""
		for conf, var in self.items():
			var.Resolve(input[conf])


def InitializeConfig(config, input):
	"""
	@brief Initializes the configuration variable.
	@param config: The configuration variable
	@param input: The Shiny input

	This function will update each configuration's resolve member, so that
	if
	"""
	for conf, var in config.items(): var.Resolve(input[conf])


def Error(message, exception=None):
	if exception:
		message = f"{message} due to {type(exception).__name__}: {exception}"

	return ui.notification_show(ui=message, type="error", duration=3)

def Msg(message): return ui.notification_show(ui=message, type="default", duration=3)


def Update(): return ui.input_action_button(
		id="Update",
		label=ui.layout_columns(
			ui.panel_conditional("input.UpdateToggle", "Auto"),
			"Update",
			ui.input_switch(id="UpdateToggle", label=None, value=True),
			col_widths=[1,9,1],
			gap="1px",
			height="1px", 	# Make it as small as possible
		)
	)


def File(input):
	if input.SourceFile() == "Upload":
		return "None" if input.File() is None else input.File()[0]["name"]
	return input.Example()
