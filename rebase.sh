#!/bin/bash

# Remove cache dirs
for dir in $(find -type d -name "__pycache__"); do
	rm -r $dir
done

# Define the list of subfolder names
subfolders=("expression" "geocoordinate" "geomap" "image" "pairwise" "3d" "spatial")
stylesheet='    <link rel="stylesheet" href="./../shinylive/heatmapper.css" />'

# Iterate through each subfolder
for folder in "${subfolders[@]}"; do
    # Check if the subfolder exists
    if [ -d "$folder" ]; then
    		cp shared/shared.py $folder/src/shared.py
        shinylive export --subdir $folder $folder/src site
        sed -i "/<link rel=\"stylesheet\" href=\"\.\/\.\.\/shinylive\/shinylive.css\" \/>/a $stylesheet" "site/$folder/index.html"
    else
        # Print a message if the subfolder doesn't exist
        echo "Subfolder '$folder' not found."
    fi
done
