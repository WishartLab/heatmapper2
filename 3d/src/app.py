#
# Heatmapper
# 3D
#
# This file contains the Shiny application for 3D Heatmapper.
# It can be run with the following command within this directory:
#		shiny run
#
# Exporting via ShinyLive is not currently supported, as pyvista
# is not yet available in the Pyodide environment. Required libraries
# include: openmpi, verdict, glew, alongside python libraries in requirements.txt
# WebGL is required for this application.
#

from numpy.core.multiarray import MAXDIMS
from shiny import App, reactive, render, ui
from pandas import DataFrame, read_table
from Bio.PDB import PDBParser, PDBIO
from io import StringIO
from imgkit import from_file as convert
from numpy import mean
from numpy.linalg import norm
from pathlib import Path

# Shared functions
from shared import Cache, MainTab, NavBar, FileSelection, Filter, ColumnType, TableOptions, InitializeConfig, ColorMaps, Update, Pyodide, Error, Msg, File

if not Pyodide: from pyvista import Plotter, plotting, read_texture, read as VistaRead

from py3Dmol import view

try:
	from user import config
except ImportError:
	from config import config


def server(input, output, session):

	# Information regarding example files.
	Info = {
		"example1.csv": {
			"Object": "bunny.obj",
			"Description": '<u>Input type:</u> .csv data, .obj object<br><u>Contents:</u> A bunny object file, mapped with random data from a .csv file. <br><u>Source:</u> <a href="https://github.com/alecjacobson/common-3d-test-models/tree/master/data"; target="_blank">github.com</a>'
		},
		"texture.jpg": {
			"Object": "FinalBaseMesh.obj",
			"Description": '<u>Input type:</u> .jpg data, .obj object<br><u>Contents:</u> A human model with a sample heatmap texture applied. The heat map is a flat image wrapped around the object.<br><u>Source:</u> <a href="https://free3d.com/3d-model/male-base-mesh-6682.html"; target="_blank">free3d.com</a>'
		},
		"4K8X.pdb": {
			"Object": None,
			"Description": '<u>Input type:</u> PDB File<br><u>Contents:</u> Binary complex of 9N DNA polymerase in the replicative state, from organism Thermococcus sp. 9oN-7. The data was determined by X-ray diffraction, with a resolution of 2.28 angstrom. <br><u>Source:</u> <a href="https://www.rcsb.org/structure/4k8x"; target="_blank">www.rcsb.org</a>'
		}
	}

	# colour scheme options
	Schemes = ["Residue #", "Reverse Residue #", "B-factor", "RMSF", "RMSD", "2ndary Structure", "pLDDT"]


	def HandleData(path, p=None):
		"""
		@brief A custom Data Handler for the Cache.
		@param n: The Path object to the file.
		@returns A data object from the cache.
		@info This Data Handler supports object files, and images as textures.
		"""

		suffix = path.suffix
		if suffix == ".obj": 
			return VistaRead(path.resolve())
		elif suffix == ".png" or suffix == ".jpg": 
			return read_texture(path.resolve())
		else: return DataCache.DefaultHandler(path)
	
	DataCache = Cache("3d", DataHandler=HandleData)

	Data = reactive.value(None)
	Valid = reactive.value(False)
	Object = reactive.value(None)

	InitializeConfig(config, input)


	@reactive.effect
	@reactive.event(input.SourceFile, input.File, input.Example, input.Reset, input.ID)
	async def UpdateData():
		"""
		@brief Updates the Data variables when input changes.
		@info When any relevant reactive input changes, this function requests data from the Cache,
		and then invalidates the current data in the table.
		"""
		if input.SourceFile() == "PDB-ID":
			try:
				data = await DataCache.Download(f"https://files.rcsb.org/view/{input.ID()}.pdb")
				Data.set(data)
			except:
				Error("File could not be loaded!\nPlease ensure you are using a valid PDB ID, and are connected to the internet.")
				return
		else:
			p = ui.Progress()
			try:
				Data.set((await DataCache.Load(input, default=None, p=p, wasm_blacklist=(".csv", ".txt", ".dat", ".tsv", ".tab", ".xlsx", ".xls", ".odf", ".png", ".jpg"))))
			except:
				p.close()
				Error(ui.HTML('File could not be loaded!<br>Protein files must be in .pdb format. <br>Object data can be uploaded as a table file, image file, or .obj file. <a href="https://github.com/WishartLab/heatmapper2/wiki/Format#3d"; target="_blank">Read more</a>'))
				return
		Valid.set(False)
		DataCache.Invalidate(File(input))


	@reactive.effect
	@reactive.event(input.SourceFile, input.Object, input.Example)
	async def UpdateObject():
		"""
		@brief Updates the Object variable when input changes.
		@info When the Object's source changes, this function fetches the most up to date information from the Cache.
		@info This function will not do anything in WebAssembly.
		"""
		if type(GetData()) != str:
			Object.set(await DataCache.Load(input,
				source_file=input.Object(),
				example_file=Info[input.Example()]["Object"],
				default=None,
				p=ui.Progress(),
				p_name="object",
				wasm=False
			))


	@reactive.effect
	@reactive.event(input.OptFile)
	def MergeOptFile():
		"""
		@brief Moves RMSD, RMSF, B-Factor or pLDDT data from an input file to the B-factor column of the PDB data
		"""
		# ignore optional file if 3D model is not PDB
		data = GetData()
		if type(data) != str:
			return
		
		# load optional file as a dataframe
		opt_file = input.OptFile()
		if opt_file is None:
			Error("Optional data could not be merged with the PDB file - Please check your data formatting.")
			return
		
		n = str(opt_file[0]["datapath"])
		# don't load files incompatible with WASM 
		# if n.endswith(blacklist) and Pyodide: return
		path = Path(n)
		try:
			opt_data = Cache.DefaultHandler(path)
		except:
			Error("File could not be loaded!\nPlease ensure you are using a properly formatted table file.")
			return
		
		# find residue name, number, chain, and value columns in data
		name_cols = ["name", "residue", "res_name"]
		num_cols = ["num", "number", "res_number"]
		chain_cols = ["chain", "letter"]
		val_cols = ["rmsd", "rmsf", "plddt", "bfactor"]

		cols = [col.lower() for col in opt_data.columns]
		name = next((col for col in name_cols if col in cols), None)
		num = next((col for col in num_cols if col in cols), None)
		chain = next((col for col in chain_cols if col in cols), None)
		val = next((col for col in val_cols if col in cols), None)
		
		if name is None or num is None or chain is None or val is None:
			Error("Additional data could not be merged! Please check your column names and formatting.")
			return
		
		# if input file uses 1 letter codes, convert to 3 letter
		aa_codes = {
			'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
			'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
			'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
			'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
		}

		# map values to atom number and chain ID
		map = {}
		for index, row in opt_data.iterrows():
			new_chain_id = row[chain].strip()
			new_res_num = row[num]
			# change 1 letter code to 3 letter if applicable
			aa = row[name].strip()
			if aa in aa_codes:
				new_name = aa_codes[aa]
			else:
				new_name = aa
			key = tuple(i for i in (new_res_num, new_chain_id, new_name) if i is not None)
			map[key] = float(row[val])
		
		# update pdb string
		updated_pdb = ""
		for line in data.splitlines():
			# handle multiple models
			if line.startswith("MODEL"):
				updated_pdb += line + "\n"
				continue
			elif line.startswith("ENDMDL"):
				updated_pdb += line + "\n"
				continue

			if line.startswith("ATOM"):
				number = int(line[22:26].strip())  # residue number
				#number = int(line[6:11].strip())  # atom number
				chain_id = line[21:22].strip()
				res_name = line[17:20].strip()
				
				# match with chain ID or residue number
				key_3 = (number, chain_id, res_name)
				key_2 = (number, chain_id)
				key_1 = (number,)
				if key_3 in map:
					new_val = map[key_3]
				elif key_2 in map:
					new_val = map[key_2]
				elif key_1 in map:
					new_val = map[key_1]
				else:  # no match found
					updated_pdb += line + "\n"
					continue

				# replace existing b-factor with new value
				# b-factor is in columns 60:66
				new_line = f"{line[:60]}{new_val:6.2f}{line[66:]}\n"
				updated_pdb += new_line

			else:
				updated_pdb += line + "\n"
		# overwrite old PDB string with new one
		Data.set(updated_pdb)
		Valid.set(False)
		DataCache.Invalidate(File(input))


	def GetData(): return Table.data_view() if Valid() else Data()


	@output
	@render.data_frame
	def Table():
		data = Data()
		
		# placeholder message if no data has been successfully uploaded
		if data is None:
			return DataFrame({"Note": ["No data to display! Please upload your data or select an example data set in the sidebar."]})
		
		# warning message if input is an image, not a table
		if isinstance(data, plotting.texture.Texture):
			df = DataFrame({"Note": ["This heatmap is mapping an image file (.png or .jpg) onto the 3D surface. There is no table data to display."]})
			return df
		
		# display residue numbers, B-factor data from PDB files
		elif isinstance(data, str):
			output_data = []
			selection = config.OptType()
			col_name = "B-factor"
			if selection == "rmsf":
				col_name = "RMSF"
			elif selection == "rmsd":
				col_name = "RMSD"
			elif selection == "plddt":
				col_name = "pLDDT"

			for line in data.splitlines():
				if line.startswith("ATOM"):
					residue_num = int(line[22:26].strip())
					residue_name = line[17:20].strip()
					atom_num = int(line[6:11].strip())
					b_factor = float(line[60:66].strip())
					output_data.append((residue_num, residue_name, atom_num, b_factor))
			df = DataFrame(output_data, columns=["Residue Number", "Residue Name", "Atom Number", col_name])
			return df
		
		# display table data
		else:
			try:
				grid = render.DataGrid(Data(), editable=True)
				Valid.set(True)  # use this table as data to generate heatmap
				return grid
			except TypeError:
				return DataFrame({"Note": ["Please ensure your uploaded file is properly formatted. The provided input format cannot be rendered."]})


	@Table.set_patch_fn
	def UpdateTable(*, patch: render.CellPatch) -> render.CellValue:
		if input.Type() == "Integer": value = int(patch["value"])
		elif input.Type() == "Float": value = float(patch["value"])
		else: value = patch["value"]
		DataCache.Invalidate(File(input))
		return value
	

	# Info text in welcome tab
	@render.ui
	def Welcome():
		return ui.HTML("""
			<h1>3D Heatmaps</h1>
			3D heatmaps render PDB models, or display data or textures on a 3D model. 
			<br><br>
			Upload a data file or specify a PDB ID in the sidebar to get started, or select 'Example' to check out a pre-loaded example.
			<br><br>	 
			Navigate to the 'Heatmap' tab to see the heatmap, or 'Table' to look at the input data.
				 
			<br><br>
			<img src="https://github.com/WishartLab/heatmapper2/wiki/assets/3D.png" alt="Image"; style="max-width:500px;">
				 
			<br><br><h3>Format</h3>
			<i>3D heatmaps can be created in two different ways:</i><br><br>
			<b>1 - PDB</b><br>
			Upload a <b>.pdb</b> file, or select 'ID' in the sidebar and enter a <b>PDB ID</b> (see example PDB 4K8X).<br>
			You can optionally upload an <b>additional table file</b> containing RMSF, RMSD, B-factor, or pDDLT values for the protein. These values will replace any values currently in the B-factor column of the PDB file.<br>
			<br><i>Optional additional files MUST include the following columns (case-insensitive):</i><br>
					<li><u>Name column:</u> "NAME", "RESIDUE", or "RES_NAME"</li>
					<li><u>Number column:</u> "NUM", "NUMBER", or "RES_NUMBER"</li>
					<li><u>Chain column:</u> "CHAIN", or "LETTER"</li>
					<li><u>Value column:</u> "RMSD", "RMSF", "PLDDT", or "BFACTOR"</li>
			<br><br>
			<b>2 - Object Files</b><br>
			Input an .obj file and either a table file or an image. If an image is used, it will be mapped onto the surface of the object (see Ex: 3D Human).<br>
			If a table file is used, values in a 'Value' column will be mapped to each face of the model. If a 'Name' column is provided, it should contain the numerical values of the faces, otherwise the values will be applied linearly (see Ex: 3D Bunny).
			<br><br>
			<table style="border-spacing: 100px";>
			<tr>
				<th>Table Files</th>
				<th>Image Files</th>
			</tr>
			<tr>
				<td style="padding-right:75px;">
					<li>.csv</li>
					<li>.dat</li>
					<li>.odf</li>
					<li>.tab</li>
					<li>.tsv</li>
					<li>.txt</li>
					<li>.xls</li>
					<li>.xlsx</li>
				</td>
				<td style="vertical-align:top;">
					<li>.jpg</li>
					<li>.png</li>
				</td>
			</tr>
			</table> 
			<br><br>
			<br><h3>Interface</h3>
			Click on the '?' icon beside sidebar options to read more about them.
			
		""")


	def PDBViewer(source, p):
		"""
		@brief Returns an HTML string containing the Py3DMol.js viewer of the source.
		@param source: The source containing PDB data as a string.
		@param p: The progress bar.
		@param image: bool indicating if output should be generated as a png for download only
		@returns An HTML string that should be wrapped with ui.HTML
		"""

		# Used for caching.
		global_inputs = [
			config.ModelType(),
			input.File() if input.SourceFile() == "Upload" else input.ID() if input.SourceFile() == "PDB-ID" else input.Example(),
			config.Size(),
			config.ColorScheme(),
			config.PFeatures(),
			config.Thickness(),
			config.Width(),
			config.Opacity(),
			config.Model(),
			config.OptType(),
		]

		if not DataCache.In(global_inputs):

			parser = PDBParser()
			structure = parser.get_structure("protein", StringIO(source))
			model = config.Model()

			def GenerateScheme(source, initial_scheme, function_declared=False, model=0):
				"""
				@brief Py3DMol has a color, colorscheme, and colorfunc attribute. This function puts the right one in without cluttering the interface with three different options.
				@param initial_scheme: The value of the scheme. Spectrum is a color, B-Color requires a function, all others
				use a colorscheme
				@param function_declared: Since there is a colorscheme for both the heatmap and the structure, we can accidentally
				redefine the same JavaScript function twice if they both use the same B-Color scheme. This avoids that.
				"""
				# dict mapping 
				scheme_dict = {
					"Residue #": "residue",
					"Reverse Residue #": "reverse",
					"B-factor": "b-factor",
					"2ndary Structure": "ssJmol",
					"RMSF": "rmsf",
					"RMSD": "rmsd",
					"pLDDT": "plddt",
				}
				
				prop = "color"
				scheme = scheme_dict[initial_scheme]

				####### colour using B-Factor data
				if scheme == "b-factor":

					# Initial weights
					darkblue, blue, lightblue, white, orange, red = 5, 10, 15, 20, 40, 50
					scheme = "Scheme"

					# Declare the function.
					if not function_declared:
						viewer.startjs += f"""\n
							let {scheme} = function(atom) {{
								if (atom.b < {darkblue}) return "darkblue"
								if (atom.b < {blue}) return "blue"
								else if (atom.b < {lightblue}) return "lightblue"
								else if (atom.b < {white}) return "white"
								else if (atom.b < {orange}) return "orange"
								else if (atom.b < {red}) return "red"
								else return "darkred"
							}}\n"""
					prop = "colorfunc"
				

				####### colour by pLDDT value
				elif scheme == "plddt":
					if config.OptType() != "plddt":
						# remind user that pLDDT requires additional input
						Msg("Is your model all grey or all red? \nUpload an additional data file with pLDDT values to visualize plDDT.")

					# Initial weights
					red, orange, yellow, lightblue, blue = 10, 50, 70, 90, 95
					scheme = "Scheme"

					# Declare the function.
					if not function_declared:
						viewer.startjs += f"""\n
							let {scheme} = function(atom) {{
								if (atom.b == 0) return "grey"
								else if (atom.b < {red}) return "red"
								else if (atom.b < {orange}) return "orange"
								else if (atom.b < {yellow}) return "yellow"
								else if (atom.b < {lightblue}) return "lightblue"
								else if (atom.b < {blue}) return "blue"
								else return "darkblue"
							}}\n"""
					prop = "colorfunc"
				

				####### colour by RMSD value
				elif scheme == "rmsd":
					
					if config.OptType() != "rmsd":
						# remind user that RMSD requires additional input
						Msg("Is your model all grey or all red? \nUpload an additional data file with RMSD values to visualize RMSD.")
					
					# Initial weights
					darkblue, blue, lightblue, white, yellow, orange, red = 0.5, 1.0, 1.5, 2, 3, 4, 5
					scheme = "rmsd"

					if not function_declared:
						viewer.startjs += f"""\n
							let {scheme} = function(atom) {{
								if (atom.b == 0) return "grey"
								else if (atom.b < {darkblue}) return "darkblue"
								else if (atom.b < {blue}) return "blue"
								else if (atom.b < {lightblue}) return "#73c9ff"
								else if (atom.b < {white}) return "white"
								else if (atom.b < {yellow}) return "#fff27d"
								else if (atom.b < {orange}) return "#ff6200"
								else if (atom.b < {red}) return "red"
								else return "darkred"
							}}\n"""
					prop = "colorfunc"


				####### colour using RMSF data
				elif scheme == "rmsf":
					
					# Initial weights
					darkblue, blue, lightblue, white, yellow, orange, red = 0.5, 1.0, 1.5, 2, 3, 4, 5
					scheme = "rmsf"
					prop = "colorfunc"

					# if additional RMSF data file added, just return scheme
					if config.OptType() == "rmsf":
						if not function_declared:
							viewer.startjs += f"""\n
								let {scheme} = function(atom) {{
									if (atom.b == 0) return "grey"
									else if (atom.b < {darkblue}) return "darkblue"
									else if (atom.b < {blue}) return "blue"
									else if (atom.b < {lightblue}) return "#73c9ff"
									else if (atom.b < {white}) return "white"
									else if (atom.b < {yellow}) return "#fff27d"
									else if (atom.b < {orange}) return "#ff6200"
									else if (atom.b < {red}) return "red"
									else return "darkred"
								}}\n"""
					
					# if no RMSF data, calculate
					else:
						if len(structure) == 1:
							Error("RMSF requires a PDB with more than one model to compute difference! Or, upload an additional data file with RMSF values.")
							if not function_declared:
								viewer.startjs += f"""\n
									let {scheme} = function(atom) {{
										if (atom.b == 0) return "grey"
										else if (atom.b < {darkblue}) return "darkblue"
										else if (atom.b < {blue}) return "blue"
										else if (atom.b < {lightblue}) return "#73c9ff"
										else if (atom.b < {white}) return "white"
										else if (atom.b < {yellow}) return "#fff27d"
										else if (atom.b < {orange}) return "#ff6200"
										else if (atom.b < {red}) return "red"
										else return "darkred"
									}}\n"""
							return source, prop, scheme
					
						Msg("Calculating RMSF...")
						# List of atom names of interest
						atom_names_of_interest = ["C", "CA", "N"]

						entry = [input.File() if input.SourceFile() == "Upload" else input.ID() if input.SourceFile() == "PDB-ID" else input.Example(), model]
						if not DataCache.In(entry):
							main_model = structure[model]
							for chain in main_model:
								for residue in chain:
									for atom in residue:
										if atom.get_id() in atom_names_of_interest:
											distances = []
											for model in structure:
												if model != main_model:
													try:
														corresponding_atom = model[chain.id][residue.id][atom.get_id()]
														distance = norm(atom.coord - corresponding_atom.coord)
														distances.append(distance)
													except KeyError: continue
											# Calculate mean distance
											if distances: atom.set_bfactor(mean(distances))
											output = StringIO()

							output = StringIO()
							io = PDBIO()
							io.set_structure(main_model)
							io.save(output)
							DataCache.Store(output.getvalue(), entry)
							output.close()

						source = DataCache.Get(entry)

						if not function_declared:
							viewer.startjs += f"""\n
								let {scheme} = function(atom) {{
									if (atom.b == 0) return "grey"
									else if (atom.b < {darkblue}) return "darkblue"
									else if (atom.b < {blue}) return "blue"
									else if (atom.b < {lightblue}) return "#73c9ff"
									else if (atom.b < {white}) return "white"
									else if (atom.b < {yellow}) return "yellow"
									else if (atom.b < {orange}) return "orange"
									else if (atom.b < {red}) return "red"
									else return "darkred"
								}}\n"""
						prop = "colorfunc"

				
				####### colour by residue number
				elif scheme == "residue":
					# get number of residues
					uniq_res = {int(line[22:26].strip()) for line in source.splitlines() if line.startswith("ATOM")}
					max = len(uniq_res)
					prop = "colorscheme"
					scheme = {
						"prop": "resi",  # index, b, resi
			   			"gradient": "linear", 
						"colors": ["red", "orange", "yellow", "green", "blue", "purple"],
						"min": 0, 
						"max": max
					}
				
				####### colour by reverse residue number
				elif scheme == "reverse":
					# get number of residues
					uniq_res = {int(line[22:26].strip()) for line in source.splitlines() if line.startswith("ATOM")}
					max = len(uniq_res)
					prop = "colorscheme"
					scheme = {
						"prop": "resi",  # index, b, resi
			   			"gradient": "linear", 
						"colors": ["purple", "blue", "green", "yellow", "orange", "red"],
						"min": 0, 
						"max": max
					}

				####### colour by secondary structure
				elif scheme == "ssJmol": 
					prop = "colorscheme"
					
				return source, prop, scheme

			viewer = view(width=f"{input.Size()}vw", height=f"{input.Size()}vh")
			source, heatmap_property, heatname_name = GenerateScheme(source, config.ColorScheme(), model=model)

			viewer.addModelsAsFrames(source)
			viewer.zoomTo()


			p.inc(message="Styling...")
			viewer.setStyle({"cartoon": {
				heatmap_property: heatname_name,
				"style": "trace" if "Simplified View" in config.PFeatures() else "rectangle",
				"thickness": config.Thickness(),
				"tubes": "Helices as Tubes" in config.PFeatures(),
				"width": config.Width(),
				"opacity": config.Opacity(),
				"scale": 1,
			}})

			if heatmap_property == "colorfunc": viewer.startjs = viewer.startjs.replace(f'"{heatname_name}"', f'{heatname_name}')

			p.inc(message="Exporting...")
			DataCache.Store(viewer.write_html(), global_inputs)

		return DataCache.Get(global_inputs)


	def ModelViewer(source, p):
		"""
		@brief Generates an HTML string of the PyVista Model viewer.
		@param source: The data to be applied to the object.
		@param p: The progress bar.
		@param image: bool indicating if output should be generated as a png for download only
		@returns An HTML string that should be wrapped with ui.HTML
		@info Object will also need to be defined.
		"""
		if Pyodide:
			# Error message pop-up
			Error(f"The WebAssembly version of Heatmapper2 does not support object rendering! Please use the Server version (server.heatmapper2.ca/3d) for this functionality.")
			# Error message in Heatmap tab
			return "The WebAssembly version of Heatmapper2 does not support object rendering! Please use the Server version (server.heatmapper2.ca/3d) for this functionality."

		# For Caching.
		inputs = [
			File(input),
			config.Style(),
			config.Opacity(),
			config.Features(),
			config.ColorMap(),
			config.Colors(),
		]

		if not DataCache.In(inputs):
			model = Object()
			if model is None: return

			p.inc(message="Plotting...")
			pl = Plotter()

			style = config.Style().lower()
			opacity = config.Opacity()
			features = config.Features()
			cmap = config.ColorMap().lower()
			colors = config.Colors()

			# If there's no source, just render the model
			if source is None:
				pl.add_mesh(
					model,
					style=style,
					opacity=opacity,
					show_edges="Edges" in features,
					lighting="Lighting" in features,
					smooth_shading="Smooth Shading" in features,
				)

			# If are data source is a table, render it as a heatmap.
			elif type(source) is DataFrame:
				values = source[Filter(source.columns, ColumnType.Name)]
				try:
					pl.add_mesh(
						model,
						scalars=values,
						style=style,
						cmap=cmap,
						opacity=opacity,
						n_colors=colors,
						show_edges="Edges" in features,
						lighting="Lighting" in features,
						smooth_shading="Smooth Shading" in features,
					)
				except NotImplementedError:
					return "Make sure you have uploaded a Table or Image file, as well as an Object file!"
				except ValueError:
					return "The number of rows in the Table file must match either the number of cells in the Object file, or the number of points in the Object file!"

			# If we have a texture, map it.
			elif type(source) is plotting.texture.Texture:
				try:
					mesh = model.texture_map_to_plane()
					pl.add_mesh(mesh, texture=source)
				except:
					return "Make sure you have uploaded a Table or Image file, as well as an Object file!"

			# Exporting as None returns the HTML as a file handle, which we read.
			p.inc(message="Exporting...")
			DataCache.Store(pl.export_html(filename=None).read(), inputs)
		return DataCache.Get(inputs)


	def GenerateHeatmap():
		"""
		@brief Generates the Heatmap based on input, returns HTML
		@returns HTML of the heatmap.
		"""
		with ui.Progress() as p:

			# Get the model and data.
			p.inc(message="Loading input...")

			source = GetData()
			if source is None: return "No data to display! <br>Please upload your data or select an example data set in the sidebar."

			if type(source) == str:
				try:
					viewer = PDBViewer(source, p)
					return viewer
				except:
					return "3D model could not be generated, please check the format of your PDB file."
				
			try:
				model_viewer = ModelViewer(source, p)
				return model_viewer
			except:
				return "3D model could not be generated, please check file formatting."


	@output
	@render.ui
	def HeatmapReactive(): return ui.HTML(GenerateHeatmap())


	@reactive.effect
	@reactive.event(input.ExampleInfo)
	def ExampleInfo():
		Msg(ui.HTML(Info[input.Example()]["Description"]))


	@render.download(filename="table.csv")
	def DownloadTable():
		df = GetData()
		# don't download if there is no data
		if df is None:
			Error("The downloaded table is empty! Please upload your data or select an example data set in the sidebar.")
		# if data is already a string
		elif isinstance(df, str):
			yield df
		else:
			yield df.to_string()


	@render.download(filename=lambda: f"heatmap{config.HeatmapType()}")
	def DownloadHeatmap():
		if config.HeatmapType() == ".png":
			# yield convert(GenerateHeatmap(), "heatmap.png")
			pass
		else: 
			html = GenerateHeatmap()
			if html is not None:
				yield html


	@output
	@render.ui
	def ConditionalElements():
		elements = []
		data = GetData()

		if data is None: return

		if type(data) == str or input.SourceFile() == "PDB-ID":
			elements += [
				# add upload box for additional data file
				ui.panel_conditional("input.SourceFile === 'Upload' || input.SourceFile === 'PDB-ID'", ui.input_file("OptFile", "Add Optional B-factor, RMSD, or RMSF Data", accept=[".csv", ".txt", ".dat", ".tsv", ".tab", ".xlsx", ".xls", ".odf"], multiple=False, placeholder='Optional Data')),
				config.OptType.UI(ui.input_radio_buttons, make_inline=True, id="OptType", label="Optional Data Type:", choices={"bfactor": "B-Factor", "rmsf": "RMSF", "rmsd": "RMSD", "plddt": "pLDDT"}, tooltip="Specify the type of data that has been uploaded as an optional additional file. Optional files should be a table file with 4 columns: atom number, residue name, chain letter, and value (B-factor, RMSF, RMSD, or pLDDT data)."),
				
				ui.HTML("<b>Customization</b>"),
				config.ColorScheme.UI(ui.input_select, id="ColorScheme", label="Color Scheme", choices=Schemes, tooltip=ui.HTML('Define the coloring of the model. <br><b>Residue #</b> - Apply a rainbow gradient based on residue number. <br><b>Reverse Residue #</b> - Apply a reversed rainbow gradient based on residue number. <br><b>B-Factor</b> - Color by B-factor. Low values are blue, and high values are red. <br><b>RMSF</b> - Color by Root Mean Square Fluctuation. Low values are blue, and high values are red. RMSF requires a PDB with more than one model to compute difference. <br><b>RMSD</b> - Color by Root Mean Square Deviation. Low values are blue, and high values are red. RMSD values must be provided in an additional table file. <br><b>2ndary Structure</b> - Color by secondary structure using the ssJmol coloring scheme. <br><b>pLDDT</b> - Color by predicted Local Distance Difference Test confidence. Low confidence values are red, and high confidence values are blue. ')),			
				
				config.Model.UI(ui.input_numeric, id="Model", label="PDB Model", min=0, tooltip="Select which model to use from the PDB file, if the PDB contains multiple models."),
				
				config.Opacity.UI(ui.input_slider, id="Opacity", label="Model Opacity", min=0.0, max=1.0, step=0.1, tooltip=ui.HTML('Specify the opacity of the <i>model</i>. 1.0 indicates full opacity, while lower values make the model more transparent.')),	
				config.Thickness.UI(ui.input_slider, id="Thickness", label="Ribbon Thickness", min=0, max=10, step=0.1, tooltip="Specify the thickness of the visualized components. Lower values make components thinner, while higher values (to a maximum of 10) make components thicker."),
				config.Width.UI(ui.input_slider, id="Width", label="Ribbon Width", min=0, max=10, step=0.1, tooltip=ui.HTML("Specify the width of the visualized components. Lower values make components narrower, while higher values (to a maximum of 10) make components wider. <br>If 'Simplified View' is selected below, Ribbon Width is ignored.")),
				config.PFeatures.UI(ui.input_checkbox_group, make_inline=False, id="PFeatures", label=None, choices=["Helices as Tubes", "Simplified View"], tooltip=ui.HTML('''Helices as Tubes - display alpha helices as simple cylinders. <br><br>Simplified View - draw the model as a simple outline. This overrides the 'Helices as Tubes' feature.''')),

		
				config.Size.UI(ui.input_numeric, id="Size", label="View Size", min=1, max=100, step=1, tooltip="Change the size of the viewer in your browser."),
			]

		else:
			elements.append(ui.panel_conditional("input.SourceFile === 'Upload'", ui.input_file("Object", "Choose an Object File", accept=[".obj"], multiple=False)))
			if type(data) == DataFrame:
				elements += [
					ui.HTML("<b>Heatmap</b>"),
					config.Opacity.UI(ui.input_slider, id="Opacity", label="Opacity", min=0.0, max=1.0, step=0.1, tooltip="Specify the opacity of the heatmap. 1.0 indicates full opacity, while lower values make the heatmap more transparent."),
					config.Style.UI(ui.input_select, id="Style", label="Style", choices=["Surface", "Wireframe", "Points"], tooltip="Specify how to render the model. Surface visualizes data as triangles on a surface. Wireframe displays a wireframe of the outer geometry. Points displays values as a collection of points."),
					ui.HTML("<b>Colors</b>"),
					config.Colors.UI(ui.input_numeric, id="Colors", label="Number", value=256, min=1, step=1, tooltip="Specify the number of colors to use to visualize data."),
					config.ColorMap.UI(ui.input_select, id="ColorMap", label="Map", choices=ColorMaps, tooltip="Specify a color scheme to use."),
					ui.HTML("<b>Features</b>"),
					config.Features.UI(ui.input_checkbox_group, make_inline=False, id="Features", label=None, choices=["Edges", "Lighting", "Smooth Shading",], tooltip=ui.HTML("Edges - display the wireframe edges on top of the surface visualization. Does not apply to Wireframe or Points. <br>Lighting - visualize the model with an external light source. Lighting may affect color accuracy. Lighting must be enabled for smooth shading to be applied. <br>Smooth Shading - smooth shadows on the surface of the model.")),
			]

		return elements

	@output
	@render.ui
	def GetInputTypes():
		if config.ModelType() =="Object":
			return FileSelection(
				examples={
					"example1.csv": "Ex: 3D Bunny", 
					"texture.jpg": "Ex: 3D Human"},
				types=[".csv", ".txt", ".dat", ".tsv", ".tab", ".xlsx", ".xls", ".odf", ".png", ".jpg"],
				project="3D",
			)
		else:
			return FileSelection(
				examples={"4K8X.pdb": "PDB 4K8X"},
				types=[".pdb"],
				project="3D",
				extras=["PDB-ID"],
			)


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

	ui.panel_title(title=None, window_title="3D"),
	NavBar(),

	ui.layout_sidebar(
		ui.sidebar(
			ui.HTML("<br>"),
			#Update(),

			config.ModelType.UI(ui.input_select, id="ModelType", label="Choose 3D Model Format", choices=["Protein", "Object"], tooltip=ui.HTML('Create a Protein heat map using a .pdb file and optional Table file, or an Object heat map using an .obj file and an Image or Table file.')),
			
			ui.output_ui(id="GetInputTypes"),

			ui.panel_conditional(
				"input.SourceFile === 'PDB-ID'",
				ui.input_text(id="ID", value="1upp", label="PDB ID"),
			),

			TableOptions(config),

			ui.panel_conditional(
				"input.MainTab === 'HeatmapTab'",

				ui.output_ui(id="ConditionalElements"),

				config.HeatmapType.UI(ui.input_radio_buttons, make_inline=False, id="HeatmapType", label="Download File Type", choices=[".html"], inline=True),
				ui.download_button(id="DownloadHeatmap", label="Download HTML"),
			),
			padding="10px",
			gap="20px",
			width="300px",
		),

		# Add the main interface tabs.
		MainTab(m_type=ui.output_ui),
		height="86vh",
	)
)

app = App(app_ui, server)
