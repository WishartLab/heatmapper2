#
# Heatmapper 
# Geomap Configuration
#
# This file contains configuration for Geomap. 


from shared import Config, ConfigHandler

config = ConfigHandler({

	# Temporal Heatmaps. True or False
	"Temporal": Config(value=False),

	# See expression's config for explanation of dynamic inputs.
	"KeyColumn": Config(),
	"ValueColumn": Config(),
	"KeyProperty": Config(),

	# "CartoDB Positron", "OpenStreetMap"
	"MapType": Config(selected="CartoDB Positron"),

	# See shared.py for ColorMaps
	"ColorMap": Config(selected="Viridis"),

	# Any floating value between 0.0 and 1.0
	"Opacity": Config(value=0.5),

	# Any number from 3-8
	"Bins": Config(value=5),

	# Allow toggling Range of Interest
	"ROI": Config(value=False),

	# Any number; the minimum/maximum bound for ROI
	"Min": Config(value=0),
	"Max": Config(value=0),

	# "Integer" "Float" "String"
	"Type": Config(selected="Integer"),

	# No value, just toggle visibility.
	"DownloadTable": Config(),
})