#
# Heatmapper
# Pairwise Configuration
#
# This file contains configuration for Pairwise.


from shared import Config, ConfigHandler

config = ConfigHandler({


	# "Distance", "Correlation"
	"MatrixType": Config(selected="Distance"),
	"HeightMatrix": Config(selected="Distance"),

	"Elevation": Config(value=90),
	"Rotation": Config(value=0),
	"Zoom": Config(value=1),
	"InterpolationLevels": Config(value=1),
	"MinScale": Config(value=False),
	"Opacity": Config(value=1.0),

	# Number between 1-50
	"TextSize": Config(value=8),

	# Check shared.py for DistanceMethods.
	"DistanceMethod": Config(selected="Euclidean"),

	# "Pearson", "Kendall", "Spearman"
	"CorrelationMethod": Config(selected="Pearson"),

	# Check shared.py for InterpolationMethods
	"Interpolation": Config(selected="Nearest"),

	# Any combination of "x", "y", "label", "legend" in a list
	"Features": Config(selected=["legend"]),

	# Any alphanumeric string that corresponds to a PDB chain.
	"Chain": Config(value="A"),

	# Any K-Mer length from 3-5
	"K": Config(value=3),

	# True or False to Use custom color maps or pre-defined ones.
	"Custom": Config(value=False),

	# See shared.py ColorMap for options.
	"ColorMap": Config(selected="Blue White Yellow"),

	# Any number between 3-100 to define amount of bins for color mapping
	"Bins": Config(value=50),

	# Any value greater than 1.
	"Size": Config(value=600),

	# Any value greater than 1.
	"DPI": Config(value=150),

	# "Integer" "Float" "String"
	"Type": Config(selected="Integer"),
})
