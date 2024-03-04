#
# Heatmapper
# Shared
#
# This file contains shared functionality between Heatmapper applications. It is not a standalone application.
# Due to the way ShinyLive exports applications, this file is symlinked into each project to reduce redundancy.
#

from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from shiny.types import FileInfo
from pandas import DataFrame, read_csv, read_excel, read_table
from io import BytesIO
from sys import modules
from copy import deepcopy
from pathlib import Path

# If pyodide is found, we're running WebAssembly.
if "pyodide" in modules:
	from pyodide.http import pyfetch
	Pyodide = True
# Otherwise,
else:
	from os.path import exists
	Pyodide = False


TemporalColumns = {"time", "date"}


def Filter(columns, good_columns, bad_columns):
	ret = None
	for column in columns:
		lowered = column.lower()
		if lowered in good_columns:
			return column
		elif lowered in bad_columns:
			break
		elif ret is not None:
			ret = column
	return ret


class Cache:
	"""
	@brief A class that encompasses fetching/storing web resources.
	"""

	@staticmethod
	def DefaultHandler(n, i):
		"""
		@brief The default handler. It can handle csv, xlsx, and defaults all other files to read_table
		@param n: The name of the file. We use this for pattern matching against the suffix.
		@param i: The binary of the file (Either via read() or BytesIO())
		@returns: A null-filled DataFrame.
		"""
		match Path(n).suffix:
			case ".csv": df = read_csv(i)
			case ".xlsx": df = read_excel(i)
			case _: df = read_table(i)
		return df.fillna(0)


	@staticmethod
	async def Remote(url): r = await pyfetch(url); return await r.bytes() if r.ok else None


	@staticmethod
	async def Local(url): return open(url, "rb").read() if exists(url) else None


	def __init__(self, project, DataHandler = DefaultHandler):
		"""
		@brief Initialize an instance of the Cache object.
		@param project: The name of the project. This is used to fetch web resources.
		@param DataHandler:	The function that should be called to process files. It should
												take a name, and a binary stream, and return a DataFrame.
		"""

		# The primary cache is immutable, and is used when the resource has not been fetched before.
		self._primary = {}

		# The secondary cache is mutable, and is populated by the primary cache. Purge deletes from here.
		self._secondary = {}

		# The data handler for processing the binary files.
		self._handler = DataHandler

		# If we're in a Pyodide environment, we fetch resources from the web.
		if Pyodide:
			self.Download = lambda url: Cache.Remote(url)
			self.Source = "https://raw.githubusercontent.com/kkernick/kkernick.github.io/main/{}/example_input/".format(project)

		# Otherwise, we fetch locally.
		else:
			self.Download = lambda url: Cache.Local(url)
			self.Source = "../example_input/"


	async def Load(self, input, copy=False):
		n = await self.N(input);
		df = DataFrame() if n is None else self._secondary[n]
		return deepcopy(df) if copy else df


	async def N(self, input):
		"""
		@brief Caches whatever the user has currently uploaded/selection, returning the identifier within the secondary cache.
		@param input: The Shiny input variable. Importantly, these must be defined:
			input.File: The uploaded file
			input.Example: The selected example file
			input.SourceFile: Whether the user wants "Upload" or "Example"
		@returns: The identifier. You should probably use Load() unless you need this.
		"""

		# Grab an uploaded file, if its done, or grab an example (Using a cache to prevent redownload)
		if input.SourceFile() == "Upload":
			file: list[FileInfo] | None = input.File()
			if file is None: return None
			n = file[0]["name"]

			# Populate the base cache, if we need to
			if n not in self._primary: self._primary[n] = self._handler(n, read(file[0]["datapath"], "wb"))

		else:
			n = input.Example()
			if n not in self._primary: self._primary[n] = self._handler(n, BytesIO(await self.Download(self.Source + n)))
		if n not in self._secondary: self._secondary[n] = deepcopy(self._primary[n])
		return n


	def Cache(self): return self._secondary


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
		row_count, column_count = df.shape
		row, column = input.TableRow(), input.TableCol()

		# So long as row and column are sane, update.
		if row < row_count and column < column_count:
			match input.Type():
				case "Integer": df.iloc[row, column] = int(input.TableVal())
				case "Float": df.iloc[row, column] = float(input.TableVal())
				case "String": df.iloc[row, column] = input.TableVal()


	async def Purge(self, input):
		"""
		@brief Purges the secondary cache of whatever the user has uploaded/selected
		@param input: The Shiny input. See N() for required objects.
		@info This function should be called on a reactive hook for a "Reset" button.
		"""
		if input.SourceFile() == "Upload":
			file: list[FileInfo] | None = input.File()
			if file is None: return None
			n = file[0]["name"]
		else:
			n = input.Example()
		del self._secondary[n]


def NavBar(current):
	"""
	@brief Returns a Navigation Bar for each project, with the current project selected.
	@returns A list, containing a ui.panel_title, and a ui.navset_bar.
	"""

	return [
			ui.panel_title(title=None, window_title="Heatmapper"),

		ui.navset_bar(
				ui.nav_panel(ui.HTML('<a href=https://kkernick.github.io/expression/site/index.html>Expression</a>'), value="Expression"),
				ui.nav_panel(ui.HTML('<a href=https://kkernick.github.io/pairwise/site/index.html>Pairwise</a>'), value="Pairwise"),
				ui.nav_panel(ui.HTML('<a href=https://kkernick.github.io/image/site/index.html>Image</a>'), value="Image"),
				ui.nav_panel(ui.HTML('<a href=https://kkernick.github.io/geomap/site/index.html>Geomap</a>'), value="Geomap"),
				ui.nav_panel(ui.HTML('<a href=https://kkernick.github.io/geocoordinate/site/index.html>Geocoordinate</a>'), value="Geocoordinate"),
				ui.nav_panel(ui.HTML('<a href=https://kkernick.github.io/about/site/index.html>About</a>'), value="About"),
				title="Heatmapper",
				selected=current,
		)
	]


def FileSelection(examples, types):
	"""
	@brief Returns the file selection dialog for the user to upload/select an example
	@param examples: Either a list of example file names, or a dictionary mapping
	@param types: The valid file extensions for user uploaded files.
	@returns A list, containing the necessary ui components for uploading/selecting
	@info The returns elements are named:
		input.SourceFile: The ui.input_radio_buttons for whether the user wants to choose an "Example" or "Upload"
		input.File: The ui.input_file for user uploaded files.
		input.Example: The ui.input_select for an example file selection
	"""

	# If the user needs help with the formatting.
	return [ui.HTML('<a href=https://kkernick.github.io/about/site/index.html>Data Format</a>'),

	# Specify whether to use example files, or upload one.
	ui.input_radio_buttons(id="SourceFile", label="Specify a Source File", choices=["Example", "Upload"], selected="Example", inline=True),

	# Only display an input dialog if the user is one Upload
	ui.panel_conditional(
		"input.SourceFile === 'Upload'",
		ui.input_file("File", "Choose a File", accept=types, multiple=False),
	),

	# Otherwise, add the example selection and an info button.
	ui.panel_conditional(
		"input.SourceFile === 'Example'",
		ui.layout_columns(
			ui.input_select(id="Example", label=None, choices=examples, multiple=False),
			ui.popover(ui.input_action_link(id="ExampleInfoButton", label="Info"), ui.output_text("ExampleInfo")),
			col_widths=[10,2],
		)
	)]


# The Table element
Table = ui.nav_panel("Table",
	ui.layout_columns(
		ui.input_numeric("TableRow", "Row", 0),
		ui.input_numeric("TableCol", "Column", 0),
		ui.input_text("TableVal", "Value", 0),
		ui.input_select(id="Type", label="Datatype", choices=["Integer", "Float", "String"]),
		col_widths=[2,2,6,2],
	),
	ui.layout_columns(
		ui.input_action_button("Update", "Update"),
		ui.input_action_button("Reset", "Reset Values"),
	),
	ui.output_data_frame("LoadedTable"),
)
