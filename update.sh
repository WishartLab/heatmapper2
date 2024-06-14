#!/bin/sh

while true; do
	echo "Pulling for Updates..."
	git pull
	echo "Done!"
	sleep 600
done
