 #!/bin/bash

# Rather than modify the system, we can just install pip things user-side, and modify the path
PATH="$PATH:$HOME/.local/bin"

# Applications that should be started
subfolders=("expression" "pairwise" "image" "geomap" "geocoordinate" "3d" "spatial")

# The initial port. It's incremented for each application.
port=8000

# Iterate through each subfolder
for folder in "${subfolders[@]}"; do

	# Kill any process that is using the port
	fuser -k $port/tcp

	cd $folder/src
	echo $(pwd)

	# Split into a separate process, at the specified port, and reload the application
	# If the source files change, so we don't need to manually teardown.
	shiny run -r -p $port &

	# Increment port, return to root
	port=$((port+1))
	cd ../..
done
