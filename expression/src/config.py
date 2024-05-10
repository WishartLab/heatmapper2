#
# Heatmapper 
# Pairwise Configuration
#
# This file contains configuration for Pairwise. 


from shared import Config, ConfigHandler

config = ConfigHandler({

	# Dynamic input's like columns change based on input.
	# If the input is visible and the selected value exists, it will be used unconditionally.
	# If the selected value doesn't exist, a value will be determined based on the input.
	# If the dynamic input is not set to be visible, the selected value is used irregardless
	# Of whether its present in the data.
	# See shared.py for ColumnType.Name for possible values.
	"NameColumn": Config(selected="UNIQID"),

	# Number between 1-50
	"ClusterMethod": Config(selected="Average"),

	# See shared.py for DistanceMethods
	"DistanceMethod": Config(selected="Euclidean"),

	# Any numerical value from 1-50
	"TextSize": Config(value=8),

	# "Row", "Column", "None"
	"ScaleType": Config(selected="Row"),

	# See shared.py for InterpolationMethods
	"Interpolation": Config(selected="None"),

	# A list of any of "row", "col", "x", "y", "legend"
	"Features": Config(selected=["row", "col", "x", "y", "legend"]),

	# True or False to Use custom color maps or pre-defined ones.
	"Custom": Config(value=False),

	# A list of colors. See shared.py for Colors
	"CustomColors": Config(selected=["Blue", "White", "Yellow"]),
	
	# See shared.py ColorMap for options. 
	"ColorMap": Config(selected="Blue White Yellow"),

	# Any number between 3-100 to define amount of bins for color mapping
	"Bins": Config(value=50),

	# "Top", "Bottom", "Left", "Right"
	"Orientation": Config(selected="Left"),

	# "Integer" "Float" "String"
	"Type": Config(selected="Float"),
})