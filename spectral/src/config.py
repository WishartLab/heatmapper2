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
	"DimensionRT": Config(value=100),
    "DimensionMZ": Config(value=495),

	# Any value in between 1-50
	"TextSize": Config(value=8),

	# See shared.py for ColorMaps
	# Plotly built-in continuous color scales:
    # https://plotly.com/python/builtin-colorscales/
	"ColorMap": Config(selected="Viridis"),

	# Any combination of "x", "y", "legend"
	"Features": Config(selected=["legend", "x", "y"]),

	# "Integer" "Float" "String"
	"Type": Config(selected="Integer"),

	# Any value greater than 1.
	"Size": Config(value=700),

	# Any value greater than 1.
	"DPI": Config(value=300),
    
	# ".txt" ".csv" ".tsv" ".xlsx"
    "TableType": Config(selected=".txt"),
    
	# ".html", ".png", ".jpg"
	"HeatmapType": Config(selected=".png"),
    
	# ".txt" ".csv" ".tsv" ".xlsx"
	"SettingType": Config(selected=".txt"),
})
