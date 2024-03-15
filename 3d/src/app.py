#
# Heatmapper
# 3D
#
# This file contains the ShinyLive application for 3D Heatmapper.
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
from numpy import linspace, zeros, float32, tan, cos, sin, radians, array, pi, ones_like, c_, ones, argsort
from scipy.interpolate import interp1d
from matplotlib.pyplot import subplots, get_cmap
from matplotlib.collections import PolyCollection
from pywavefront import Wavefront
from io import BytesIO
from pathlib import Path
from tempfile import NamedTemporaryFile

# Shared functions
from shared import Cache, MainTab, NavBar, FileSelection, Filter, ColumnType, TableValueUpdate


def server(input, output, session):

	# Information regarding example files.
	Info = {
		"example1.csv": {
			"Object": "bunny.obj",
			"Description": "A bunny, mapped with the distance from the starting Z position"
		}
	}

	def HandleData(n, i):
		"""
		@brief A custom Data Handler for the Cache.
		@param n: The name of the file
		@param i: The source of the file. It can be a path to a file (string) or a BytesIO object.
		@returns A data object from the cache.
		@info This Data Handler supports object files as Wavefront files
		"""
		match Path(n).suffix:
			case ".obj":
				# If the file is sourced from a BytesIO, we need to store it in a file temporarily.
				if type(i) is BytesIO:
					temp = NamedTemporaryFile(); temp.write(i.read())
					i = temp.name
				return Wavefront(i, create_materials=True, collect_faces=True)
			case _: return DataCache.DefaultHandler(n, i)
	DataCache = Cache("3d", DataHandler=HandleData)


	def Frustum(left, right, bottom, top, znear, zfar):
		"""
		@brief Returns a viewing volume giving clipping planes
		@param left: left clipping plane
		@param right: right clipping plane
		@param bottom: bottom clipping plane
		@param top: top clipping plane
		@param znear: near clipping plame
		@param zfar: far clipping plane
		@returns The viewing volume
		@info https://matplotlib.org/matplotblog/posts/custom-3d-engine/
		"""
		M = zeros((4, 4), dtype=float32)
		M[0, 0] = +2.0 * znear / (right - left)
		M[1, 1] = +2.0 * znear / (top - bottom)
		M[2, 2] = -(zfar + znear) / (zfar - znear)
		M[0, 2] = (right + left) / (right - left)
		M[2, 1] = (top + bottom) / (top - bottom)
		M[2, 3] = -2.0 * znear * zfar / (zfar - znear)
		M[3, 2] = -1.0
		return M


	def Perspective(fovy, aspect, znear, zfar):
		"""
		@brief Return a perspective projection
		@param fovy: field of View
		@param aspect: aspect
		@param znear: near clipping plane
		@param zfar: far clipping plane
		@return The viewing volume
		@info https://matplotlib.org/matplotblog/posts/custom-3d-engine/
		"""

		h = tan(0.5*radians(fovy)) * znear
		w = h * aspect
		return Frustum(-w, w, -h, h, znear, zfar)


	def Translate(x, y, z):
		"""
		@brief Apply a translation to move the model around"
		@param x: The x translation
		@param y: The y translation
		@param z: The z translation
		@return The translated numpy array.
		@info https://matplotlib.org/matplotblog/posts/custom-3d-engine/
		"""

		return array([[1, 0, 0, x],
											[0, 1, 0, y],
											[0, 0, 1, z],
											[0, 0, 0, 1]], dtype=float)


	def Y(theta):
		"""
		@param Rotate about the Y axis
		@param theta: The angle to rotate
		@info https://matplotlib.org/matplotblog/posts/custom-3d-engine/
		"""

		t = pi * theta / 180
		c, s = cos(t), sin(t)
		return array([[1, 0,  0, 0],
											[0, c, -s, 0],
											[0, s,  c, 0],
											[0, 0,  0, 1]], dtype=float)


	def X(theta):
		"""
		@param Rotate about the X axis
		@param theta: The angle to rotate
		@info https://matplotlib.org/matplotblog/posts/custom-3d-engine/
		"""

		t = pi * theta / 180
		c, s = cos(t), sin(t)
		return  array([[ c, 0, s, 0],
											[ 0, 1, 0, 0],
											[-s, 0, c, 0],
											[ 0, 0, 0, 1]], dtype=float)


	def DataFrameHeatmap(V, df):
		"""
		@brief Generate a colormap to apply to the model via a DataFrame
		@param V: the set of vertices associated with the model
		@param df: The DataFrame
		@returns The generated colormap if df is defined, a uniform map otherwise.
		"""

		# If the DataFrame is None, just generate a uniform mapping.
		if df.empty: return get_cmap(input.ColorMap().lower())([0.5 for _ in range(len(V))])

		# If there is a name column, make sure the triangles are in order.
		name_col = Filter(df.columns, ColumnType.Name, only_one=True)
		if name_col is not None: df.sort_values(name_col)

		# Normalize the data, make a colormap
		data = df[Filter(df.columns, ColumnType.Value, only_one=True)]
		m, M = data.min(), data.max()
		data = (data-m)/(M-m)

		# Ensure that the size of the data matches the triangles of the object.
		if len(data) != len(V):
			new_values = linspace(data.min(), data.max(), len(V))
			interp_func = interp1d(data, ones_like(data))
			data = interp_func(new_values)
		return get_cmap(input.ColorMap().lower())(data)


	async def GenerateHeatmap():
		"""
		@brief Generates a Heatmap.
		"""

		# Get the source and model. We don't need a source file, we do need a model.
		model = await DataCache.Load(
			input,
			source_file=input.Object(),
			example_file=Info[input.Example()]["Object"],
			default=None
		)
		source = await DataCache.Load(input)
		if model is None: return

		# Get the faces and vertices
		V, F = array(model.vertices), array(model.mesh_list[0].faces)

		# Magic math I don't understand ;)
		V = (V-(V.max(0)+V.min(0))/2) / max(V.max(0)-V.min(0))
		model = Y(input.Y()) @ X(input.X())
		view = Translate(0,0,-3.5)
		proj = Perspective(input.Zoom(), 1, 1, 100)
		MVP = proj @ view @ model
		V = c_[V, ones(len(V))]  @ MVP.T
		V /= V[:,3].reshape(-1,1)
		V = V[F]
		T =  V[:,:,:2]
		Z = -V[:,:,2].mean(axis=1)

		# Generate a heatmap
		C = DataFrameHeatmap(V, source)

		# Render back to front
		I = argsort(Z)
		T, C = T[I,:], C[I,:]

		# Setup the figure.
		fig, ax = subplots(figsize=(6,6))
		ax.set_xlim([-1,+1])
		ax.set_ylim([-1,+1])
		ax.axis('off')
		ax.set_aspect('equal')
		collection = PolyCollection(T, closed=True, linewidth=0.1, facecolor=C, edgecolor="black")
		ax.add_collection(collection)
		return ax


	@output
	@render.data_frame
	@reactive.event(input.SourceFile, input.File, input.Example, input.File, input.Update, input.Reset)
	async def LoadedTable(): return await DataCache.Load(input)


	@output
	@render.plot
	@reactive.event(input.SourceFile, input.File, input.Example, input.Object, input.Update, input.Reset, input.ColorMap, input.X, input.Y, input.Zoom)
	async def Heatmap(): return await GenerateHeatmap()


	@output
	@render.text
	@reactive.event(input.SourceFile, input.Example)
	def ExampleInfo(): return Info[input.Example()]["Description"]


	@render.download(filename="table.csv")
	async def DownloadTable():
		df = await DataCache.Load(input);
		if df is not None:
			yield df.to_string()


	@reactive.Effect
	@reactive.event(input.Update)
	async def Update(): await DataCache.Update(input)


	@reactive.Effect
	@reactive.event(input.Reset)
	async def Reset(): await DataCache.Purge(input)


	@reactive.Effect
	@reactive.event(input.SourceFile, input.File, input.Example, input.TableRow, input.TableCol, input.Update, input.Reset)
	async def UpdateTableValue():
		df = await DataCache.Load(input)
		if df is not None:
			TableValueUpdate(df, input)


app_ui = ui.page_fluid(

	NavBar("3D"),

	ui.layout_sidebar(
		ui.sidebar(

			FileSelection(examples={"example1.csv": "Example 1"}, types=[".csv", ".txt", ".xlsx"]),

			ui.panel_conditional("input.SourceFile === 'Upload'", ui.input_file("Object", "Choose an Object File", accept=[".obj"], multiple=False)),

			# Set the ColorMap used.
			ui.input_select(id="ColorMap", label="Color Map", choices=["Viridis", "Plasma", "Inferno", "Magma", "Cividis"], selected="Viridis"),

			ui.input_slider(id="X", label="X Rotation", value=45, min=0, max=360, step=5),
			ui.input_slider(id="Y", label="Y Rotation", value=20, min=0, max=360, step=5),
			ui.input_slider(id="Zoom", label="Zoom", value=25, min=0, max=50, step=5),

			# Add the download buttons.
			ui.download_button("DownloadTable", "Download Table"),
		),

		# Add the main interface tabs.
		MainTab(),
	)
)

app = App(app_ui, server)
