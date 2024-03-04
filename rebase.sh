#!/bin/bash

# Define the list of subfolder names
subfolders=("expression" "geocoordinate" "geomap" "image" "pairwise")

# Iterate through each subfolder
for folder in "${subfolders[@]}"; do
    # Check if the subfolder exists
    if [ -d "$folder" ]; then
        # Change directory to the subfolder
        cd "$folder" || exit
        # Run the command in the subfolder
        shinylive export src site
        # Change back to the parent directory
        cd ..
    else
        # Print a message if the subfolder doesn't exist
        echo "Subfolder '$folder' not found."
    fi
done
