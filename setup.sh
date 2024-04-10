#!/bin/bash

echo "Creating Virtual Environment..."
python -m venv venv

echo "Installing Dependencies"
source venv/bin/activate
pip install shiny shinylive matplotlib scipy pandas biopython folium branca certifi xyzservices pyvista trame vtk trame-vtk meshio pywebview trame-vuetify nest_asyncio squidpy anndata scanpy scikit-misc

echo "Setting up Repository"
git clone https://github.com/WishartLab/heatmapper2.git
cd heatmapper2
git lfs install
git lfs checkout
git lfs fetch

echo "Done!"