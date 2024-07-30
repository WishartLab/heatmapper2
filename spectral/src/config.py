#
# Heatmapper
# Image Configuration
#
# This file contains configuration for Image.


from shared import Config, ConfigHandler

config = ConfigHandler({

	# Dependent on input
	"Index": Config(selected=[0]),
	"Peaks": Config(selected="Raw"),
	"Width": Config(value=1.0),
	"Interpolation": Config(selected="Nearest"),

	# Any value in between 1-50
	"TextSize": Config(value=8),

	# Any value between 0.0-1.0
	"Opacity": Config(value=0.5),

	# See shared.py for ColorMaps
	"ColorMap": Config(selected="Viridis"),

	# Any combination of "x", "y", "legend"
	"Features": Config(selected=["legend"]),

	# "Integer" "Float" "String"
	"Type": Config(selected="Integer"),

	# Any value greater than 1.
	"Size": Config(value=700),

	# Any value greater than 1.
	"DPI": Config(value=300),
})
