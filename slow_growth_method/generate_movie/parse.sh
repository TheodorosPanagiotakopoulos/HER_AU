#!/bin/bash

for i in */; do
    cd "$i"

    if [[ -f "OUTCAR.gz" && ! -f "OUTCAR" ]]; then
        path=$(pwd)
		rm -f DATA.js MOVIE.traj
        echo "Working dir path: $path"
        echo "OUTCAR.gz found!"
        echo -n "Unzipping OUTCAR.gz..."
        gzip -d OUTCAR.gz && echo " Done."
        echo "Removing binary characters from OUTCAR if they exist..."
        src_dir="/home/theodoros/PROJ_ElectroCat/theodoros/HER/Au/HER_Au/slow_grow_method/CH3NH3/3_CH3NH3/CH3NH3_splitting/3_CH3NH3_40_H2O_v2"
        if [[ ! -f "fix.py" ]]; then
            cp "$src_dir/fix.py" .
            echo "Copied fix.py"
        fi
        if [[ ! -f "parse.py" ]]; then
            cp "$src_dir/parse.py" .
            echo "Copied parse.py"
        fi
        python3 fix.py
        echo -n "Zipping OUTCAR..."
        gzip OUTCAR && echo " Done."
        if [[ -f "CONTCAR" ]]; then
            line=$(sed -n '7p' CONTCAR)
            python3 parse.py "$line"
        fi
    elif [[ ! -f "OUTCAR.gz" || (-f "OUTCAR.gz" && -f "OUTCAR") ]]; then
        echo -e "\e[5mSkipping directory (OUTCAR.gz missing or both OUTCAR and OUTCAR.gz exist): $PWD\e[0m"
    fi
    cd ../
done
