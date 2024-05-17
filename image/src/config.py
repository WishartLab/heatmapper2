#
# Heatmapper 
# Image Configuration
#
# This file contains configuration for Image. 


from shared import Config, ConfigHandler

config = ConfigHandler({

	# Any value in between 1-50
	"TextSize": Config(value=8),

	# Any value between 0.0-1.0
	"Opacity": Config(value=0.5),

	# See shared.py for ColorMaps
	"ColorMap": Config(selected="Viridis"),

	# One of "MPL2005", "MPL2014", "Serial", "Threaded"
	"Algorithm": Config(selected="MPL2014"),

	# Any value between 1-100
	"Levels": Config(value=20),

	# Any combination of "x", "y", "legend"
	"Features": Config(selected=["legend"]),

	# "Integer" "Float" "String"
	"Type": Config(selected="Integer"),

	# Any value greater than 1.
	"Size": Config(value=700),

	# Any value greater than 1.
	"DPI": Config(value=300),
})