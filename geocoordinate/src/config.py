#
# Heatmapper
# Geocoordinate Configuration
#
# This file contains configuration for Geocoordinate.


from shared import Config, ConfigHandler

config = ConfigHandler({

	# See expression's config for explanation of dynamic inputs.
	"TimeColumn": Config(),
	"ValueColumn": Config(),

	"RenderMode": Config(selected="Raster"),
	"RenderShape": Config(selected="Circle"),

	# "CartoDB Positron", "OpenStreetMap"
	"MapType": Config(selected="CartoDB Positron"),

	# Any floating value between 0.0 and 1.0
	"Opacity": Config(value=0.7),

	# Any numerical value from 5-100
	"Radius": Config(value=25),

	# Any numerical value from 1-30
	"Blur": Config(value=15),

	# Allow toggling Range of Interest
	"ROI": Config(value=False),

	# "Remove" or "Round"
	"ROI_Mode": Config(selected="Remove"),

	# Any number; the minimum/maximum bound for ROI
	"Min": Config(value=0),
	"Max": Config(value=0),

	# "Temporal", "Uniform", "Scaled"
	"Features": Config(selected=["KDE"]),

	# "Integer" "Float" "String"
	"Type": Config(selected="Integer"),

	# No value, just toggle visibility.
	"DownloadTable": Config(),
})
