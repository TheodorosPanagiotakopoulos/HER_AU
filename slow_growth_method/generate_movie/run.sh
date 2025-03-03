#! /bin/bash

a=`ls -d */`
for i in $a
do
	cd $i
	if ls frame*.png 1> /dev/null 2>&1; then
    	echo "Files found! Performing action..."
		python add_time.py
		python movie.py
	else
    	echo "No matching files found."
	fi
	cd ../
done
