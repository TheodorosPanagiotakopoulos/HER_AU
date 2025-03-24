# Slow Growth Method Simulations

This repository contains scripts to parse and generate images from slow growth method simulations.

## Directory Structure

Your simulations should be organized in the following format:

root_directory/  
│── slow_growth_method/  
│   ├── parse.py  
│   ├── parse_CH3NH3.py  (if the molecule is CH3NH3)  
│   ├── fix.py            (if the molecule is CH3NH3)  
│   ├── generate_image.py  
│   ├── get_data.py  
│   ├── get_mols.py  
│── RUN1/  
│── RUN2/  
│── ...  
│── RUNN/

## Instructions

### 1. Run `parse.py`

Navigate to the `slow_growth_method` directory and run:

python parse.py

This script processes the simulation data from the RUNX folders.

Note: If the cation is CH3NH3, use parse_CH3NH3.py instead of parse.py.

In that case:

    First, copy fix.py into every run folder (RUN1, RUN2, ..., RUNN).

    Then, navigate to the main directory (where all the RUN folders and parse.sh are located), and run:

This will execute parse_CH3NH3.py inside each RUN folder automatically.
###2. Generate Images

In the slow_growth_method directory, run:

python generate_image.py

Make sure that get_data.py and get_mols.py are also present in the same directory, as they are required for generating images.
Requirements

Ensure you have the necessary dependencies installed. You can install them using:

pip install -r requirements.txt

This setup ensures a smooth workflow from parsing simulation data to generating the final images or movie. Happy simulating! 🚀
