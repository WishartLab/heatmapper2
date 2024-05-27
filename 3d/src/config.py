#
# Heatmapper 
# 3D Configuration
#
# This file contains configuration for 3D. 


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

	# "Cartoon", "Stick", "Sphere", "Line"
	"PStyle": Config(selected="Cartoon"),

	# See 3D's ConditionalElements function for details
	"ColorScheme": Config(selected="spectrum"),

	# Can vary depending on PStyle. See ConditionalElements for details.
	"PFeatures": Config(selected=[]),

	# "Trace", "Oval", "Rectangle", "Parabola", "Edged"
	"PCStyle": Config(selected="Rectangle"),

	# Any value from 0-10
	"Thickness": Config(value=0.4),

	# Any value from 0-10
	"Width": Config(value=3),

	# Any value from 0-10
	"Radius": Config(value=0.25),

	# Any value from 0-10
	"Scale": Config(value=1),

	# Any numerical value
	"Size": Config(value=600),

	"HeatmapType": Config(selected="Builtin")

})