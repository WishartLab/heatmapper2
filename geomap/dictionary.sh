#!/bin/bash

# Function to titleize a string
titleize() {
    echo "$1" | awk '{for(i=1;i<=NF;i++){$i=toupper(substr($i,1,1)) tolower(substr($i,2))} print}'
}

# Initialize an empty Python dictionary string
python_dict="{"

# Iterate over all files in the current directory
for file in data/*; do
    # Check if the file is a regular file
    if [ -f "$file" ]; then
        # Extract filename without extension
        filename=$(basename -- "$file")
        filename_noext="${filename%.*}"

        # Replace dashes with spaces and titleize
        modified=$(titleize "${filename_noext//-/ }")

        # Append to Python dictionary string
        python_dict="$python_dict \"$filename\": \"$modified\","
    fi
done

# Remove the trailing comma and close the dictionary
python_dict="${python_dict%,}"
python_dict="$python_dict }"

# Print the Python dictionary
echo "$python_dict"
