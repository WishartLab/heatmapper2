#
# Heatmapper 
# Pairwise Configuration
#
# This file contains configuration for Pairwise. 


from shared import Config, ConfigHandler

config = ConfigHandler({

	# "Distance", "Correlation"
	"Opacity": Config(value=1.0),

	# Any numerical value between 1 and 256.
	"Colors": Config(value=256),

	# "Viridis", "Plasma", "Inferno", "Magma", "Cividis"
	"ColorMap": Config(selected="Viridis"),

	# "Integer" "Float" "String"
	"Type": Config(selected="Integer"),

	# "Surface", "Wireframe", "Points"
	"Style": Config(selected="Surface"),

	# Combination of "Edges", "Lighting", "Interpolation", "Smooth Shading"
	"Features": Config(selected=["Lighting", "Interpolation", "Smooth Shading"]),

	# No value, just toggle visibility.
	"DownloadTable": Config(),
})