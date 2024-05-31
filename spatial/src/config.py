#
# Heatmapper
# Spatial Configuration
#
# This file contains configuration for Spatial.


from shared import Config, ConfigHandler

config = ConfigHandler({


	# "obs", "var"
	"TableType": Config(selected="obs"),

	# "moran", "sepal", "geary"
	"Statistic": Config(selected="moran"),

	# See expression's config for description of dynamic inputs
	"Keys": Config(),

	# See shared.py for ColorMaps
	"ColorMap": Config(selected="Viridis"),

	# "Circle", "Square", "Hex"
	"Shape": Config(selected="Hex"),

	# A value from 0.0-1.0
	"ImgOpacity": Config(value=1.0),

	# A value from 0.0-1.0
	"Opacity": Config(value=1.0),

	# Any value from 1-10
	"Columns": Config(value=2),

	# Any value from 0.0-1.0
	"Spacing": Config(value=0.1),

	# Any value greater than 1.
	"Size": Config(value=600),

	# Any value greater than 1.
	"DPI": Config(value=200),

	# Any of "Image", "Legend", "Frame"
	"Features": Config(selected=["Image", "Legend"]),

	# One of "closeness_centrality", "average_clustering", "degree_centrality"
	"Score": Config(selected="closeness_centrality"),

	# "L", "F", "G"
	"Function": Config(selected="L"),

	# See shared.py DistanceMethods
	"Distance": Config(selected="Euclidean"),

	# "Scatter", "Line"
	"OccurrenceGraph": Config(selected="Scatter"),

	# See expression's configuration for dynamic inputs
	"Cluster": Config(),

	# Any value from 1-50
	"Interval": Config(value=50),

	# Any value from 0-10
	"Splits": Config(value=0),

	# No value, just toggle visibility.
	"DownloadTable": Config(),
})
