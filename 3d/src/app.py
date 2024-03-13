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

from shiny import App, Inputs, Outputs, Session, reactive, render, ui
import numpy as np
from matplotlib.pyplot import subplots
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from pandas import DataFrame
from pywavefront import Wavefront
from io import BytesIO, StringIO
from shared import Cache, MainTab, NavBar, FileSelection, Filter, ColumnType, FillColumnSelection, TableValueUpdate


def server(input: Inputs, output: Outputs, session: Session):

	# Information regarding example files.
	Info = {
		"example1.csv": {
			"Object": "bunny.obj",
			"Description": "A bunny, mapped with random data"
		}
	}

	DataCache = Cache("3d")


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
		M = np.zeros((4, 4), dtype=np.float32)
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
		h = np.tan(0.5*np.radians(fovy)) * znear
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
		return np.array([[1, 0, 0, x],
											[0, 1, 0, y],
											[0, 0, 1, z],
											[0, 0, 0, 1]], dtype=float)


	def Y(theta):
		"""
		@param Rotate about the Y axis
		@param theta: The angle to rotate
		@info https://matplotlib.org/matplotblog/posts/custom-3d-engine/
		"""
		t = np.pi * theta / 180
		c, s = np.cos(t), np.sin(t)
		return np.array([[1, 0,  0, 0],
											[0, c, -s, 0],
											[0, s,  c, 0],
											[0, 0,  0, 1]], dtype=float)


	def X(theta):
		"""
		@param Rotate about the X axis
		@param theta: The angle to rotate
		@info https://matplotlib.org/matplotblog/posts/custom-3d-engine/
		"""
		t = np.pi * theta / 180
		c, s = np.cos(t), np.sin(t)
		return  np.array([[ c, 0, s, 0],
											[ 0, 1, 0, 0],
											[-s, 0, c, 0],
											[ 0, 0, 0, 1]], dtype=float)


	async def LoadObject():
		"""
		@brief Loads the image to render behind the heatmap.
		@returns an Image object, if an image is specified, otherwise None.
		"""

		# Grab an uploaded file, if its done, or grab an example (Using a cache to prevent redownload)
		if input.SourceFile() == "Upload":
			file: list[FileInfo] | None = input.Object()
			return None if file is None else Wavefront(file[0]["datapath"], create_materials=True, collect_faces=True)
		else:
			n = Info[input.Example()]["Object"]
			cache = DataCache.Cache()
			if n not in cache:
				# Unfortunate, but we need to write to a temporary file.
				open("temp.obj", "wb").write(await DataCache.Download(DataCache.Source + n))
				cache[n] = Wavefront("temp.obj", create_materials=True, collect_faces=True)
			return cache[n]


	async def GenerateHeatmap():
		model = await LoadObject()
		df = await DataCache.Load(input)

		# If there is a name column, make sure the triangles are in order.
		name_col = Filter(df.columns, ColumnType.Name, only_one=True)
		if name_col is not None: df.sort_values(name_col)

		# Get the faces and vertices
		V, F = np.array(model.vertices), np.array(model.mesh_list[0].faces)

		# Magic math I don't understand ;)
		V = (V-(V.max(0)+V.min(0))/2) / max(V.max(0)-V.min(0))
		model = Y(input.Y()) @ X(input.X())
		view = Translate(0,0,-3.5)
		proj = Perspective(input.Zoom(), 1, 1, 100)
		MVP = proj @ view @ model
		V = np.c_[V, np.ones(len(V))]  @ MVP.T
		V /= V[:,3].reshape(-1,1)
		V = V[F]
		T = V[:,:,:2]

		# Get the value of each face, normalize, create a colormap
		Z = df[Filter(df.columns, ColumnType.Value, only_one=True)]
		zmin, zmax = Z.min(), Z.max()
		Z = (Z-zmin)/(zmax-zmin)
		C = plt.get_cmap(input.ColorMap().lower())(Z)

		# Render back to front so it looks pretty.
		I = np.argsort(Z)
		T, C = T[I,:], C[I,:]

		# Rendering
		fig, ax = plt.subplots(figsize=(6,6))
		# Set the axis limits and remove the frame
		ax.set_xlim([-1,+1])
		ax.set_ylim([-1,+1])
		ax.axis('off')

		# Set the aspect ratio to be equal
		ax.set_aspect('equal')
		collection = PolyCollection(T, closed=True, linewidth=0.1, facecolor=C, edgecolor="black")
		ax.add_collection(collection)
		return ax

	@output
	@render.data_frame
	@reactive.event(input.Update, input.Reset, input.Example, input.File, ignore_none=False, ignore_init=False)
	async def LoadedTable(): return await DataCache.Load(input)


	@output
	@render.plot
	@reactive.event(input.Update, input.Reset, input.Example, input.File, input.ColorMap, input.X, input.Y, input.Zoom, input.Object, ignore_none=False, ignore_init=False)
	async def Heatmap(): return await GenerateHeatmap()


	@output
	@render.text
	def ExampleInfo(): return Info[input.Example()]["Description"]

	@render.download(filename="table.csv")
	async def DownloadTable(): df = await DataCache.Load(input); yield df.to_string()


	@reactive.Effect
	@reactive.event(input.Update)
	async def Update(): await DataCache.Update(input)


	@reactive.Effect
	@reactive.event(input.Reset)
	async def Reset(): await DataCache.Purge(input)


	@reactive.Effect
	@reactive.event(input.TableRow, input.TableCol, input.Example, input.File, input.Reset, input.Update)
	async def UpdateTableValue(): TableValueUpdate(await DataCache.Load(input), input)


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