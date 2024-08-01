#
# Heatmapper
# Image Configuration
#
# This file contains configuration for Image.


from shared import Config, ConfigHandler

config = ConfigHandler({

	# Dependent on input
	"ID": Config(selected=[]),
	"Peaks": Config(selected="Raw"),
	"Interpolation": Config(selected="Nearest"),
	"Dimension": Config(value=100),

	"Elevation": Config(value=45),
	"Rotation": Config(value=45),
	"Zoom": Config(value=1),

	# Any value in between 1-50
	"TextSize": Config(value=8),

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
