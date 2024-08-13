# Heatmapper 2

Welcome to Heatmapper! This documentation refers to [Heatmapper 2](https://heatmapper2.ca). If you’re using the older, R-Based Heatmapper, it’s located here: [http://heatmapper.ca](http://heatmapper.ca).

If you’ve used the legacy Heatmapper, you may want to look at [Running](https://github.com/WishartLab/heatmapper2/wiki/Running) to see what’s new between versions, or if you want to better understand the core technologies used by the project. In summary, there are three ways to run Heatmapper:
1. As a Web-Assembly, Client-Side application run in your browser, hosted at https://heatmapper2.ca
2. As a PyShiny, Server-Side application run in your browser, hosted at http://35.208.86.138
3. As a Python application from your computer. See [Running](https://github.com/WishartLab/heatmapper2/wiki/Running) for details.

If your browser and system supports it, Heatmapper will run as a client-side, WASM application within your browser by default. Some applications, such as the Spatial Heatmapper, can only be run on the server, but Heatmapper will automatically switch the source based on the selected application and capabilities of the system.

If you’re looking for documentation, it’s [here](https://github.com/WishartLab/heatmapper2/wiki). For each heatmap, you can find general information about the heatmap is for, the input files that it supports, and a description of the settings and features that can be modified within the interface. For an outline of the interface, and an explanation of the various settings for each application, see [Interface](https://github.com/WishartLab/heatmapper2/wiki/Interface). For an outline of the support file formats each application supports, see [Format](https://github.com/WishartLab/heatmapper2/wiki/Format).

## Heatmap Overview

| Heatmap                                                               | Operating Mode         |                                                                                                                                                                             | Screenshot                                                                                |
| --------------------------------------------------------------------- | ---------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------- |
| [Expression](https://heatmapper2.ca/expression/site/index.html)       | All Modes              | Displays unclustered (or previously clustered) data from transcriptomic (microarray or RNAseq), proteomic or metabolomic experiments.                                       | ![Expression](https://github.com/WishartLab/heatmapper2/wiki/assets/Expression.png)       |
| [Pairwise](https://heatmapper2.ca/pairwise/site/index.html)           | All Modes              | Displays two sub-classes of heatmaps: display all pairwise distances between the points in a data set or display correlations between all pairs of variables in a data set. | ![Pairwise](https://github.com/WishartLab/heatmapper2/wiki/assets/Pairwise.png)           |
| [Image](https://heatmapper2.ca/image/site/index.html)                 | All Modes              | Display a heatmap over any image.                                                                                                                                           | ![Image](https://github.com/WishartLab/heatmapper2/wiki/assets/Image.png)                 |
| [Geomap](https://heatmapper2.ca/geomap/site/index.html)               | All Modes              | Display a heatmap based on country, state, province etc. political boundaries.                                                                                              | ![Geomap](https://github.com/WishartLab/heatmapper2/wiki/assets/Geomap.png)               |
| [Geocoordinate](https://heatmapper2.ca/geocoordinate/site/index.html) | All Modes              | Display data on geospatial coordinates (latitude and longitude).                                                                                                            | ![Geocoordinate](https://github.com/WishartLab/heatmapper2/wiki/assets/Geocoordinate.png) |
| [3D](https://heatmapper2.ca/3D/site/index.html)                       | Partial                | Displays data or textures on a 3D, interactive model.                                                                                                                       | ![3D](https://github.com/WishartLab/heatmapper2/wiki/assets/3D.png)                       |
| [Spatial](http://35.208.86.138:8006)                                  | No WebAssembly Support | Display spatial molecular data, and visualize various related metrics.                                                                                                      | ![Spatial](https://github.com/WishartLab/heatmapper2/wiki/assets/Spatial.png)             |
| [Spectral](https://heatmapper2.ca/spectral/site/index.html)           | All Modes              | Display Mass Spectrometry data.                                                                                                                                             | ![Spectral](https://github.com/WishartLab/heatmapper2/wiki/assets/Spectral.png)           |

## References

Thanks to the following libraries for making Heatmapper possible:
* [Shiny](https://shiny.posit.co/py/)
* [Matplotlib](https://matplotlib.org/)
* [SciPy](https://scipy.org/)
* [Pandas](https://pandas.pydata.org/)
* [Biopython](https://biopython.org/)
* [Folium](https://python-visualization.github.io/folium/latest/)
* [PyVista](https://pyvista.org/)
* [Trame](https://kitware.github.io/trame/)
* [Squidpy](https://squidpy.readthedocs.io/en/stable/#)
* [AnnData](https://anndata.readthedocs.io/en/latest/)
* [Scanpy](https://scanpy.readthedocs.io/en/stable/)
* [Scikit-Misc](https://has2k1.github.io/scikit-misc/stable/)
