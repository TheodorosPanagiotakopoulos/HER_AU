Slow Growth Method Simulations

This repository contains scripts to parse, process, and generate movies from slow growth method simulations.
Directory Structure

Your simulations should be organized in the following format:

root_directory/  
│── slow_growth_method/  
│   ├── parse.py  
│   ├── parse_CH3NH3.py  
│   ├── generate_image.py  
│   ├── get_data.py  
│   ├── get_mols.py  
│   ├── add_time.py  
│   ├── movie.py  
│   ├── fix.py  
│── RUN1/  
│── RUN2/  
│── ...  
│── RUNN/  
│── automate_movie.py  
│── run.sh  
│── parse.sh

Instructions
1. Run parse.py

Navigate to the slow_growth_method directory and execute:

    python parse.py

This script processes the simulation data from the RUNX folders.

Note: If the cation is CH3NH3, you must use parse_CH3NH3.py instead of parse.py.

In that case:

- First, copy fix.py into every run folder (RUN1, RUN2, ..., RUNN).

- Then, go to the main directory (where you see all the RUN folders and parse.sh) and run:

bash parse.sh

This will execute parse_CH3NH3.py inside each RUN folder automatically.
2. Generate Images

In the slow_growth_method directory, run:

    python generate_image.py

Make sure get_data.py and get_mols.py are present in the same directory, as they are required for generating images.
3. Add Time and Process Movie

Ensure that add_time.py and movie.py are present in the slow_growth_method directory. These scripts are required for the next steps in movie creation.
4. Run automate_movie.py

Move one folder back, to the main directory where all the slow growth method simulations (RUN1, RUN2, ...) are stored. Then, execute:

    python automate_movie.py

This script automates the movie creation process. It requires run.sh in the same directory to function properly.
5. Modify the Output Path

Before running automate_movie.py, ensure that the path where the generated video should be saved is correctly set in the script.
Requirements

Ensure you have the necessary dependencies installed. If needed, install them using:

    pip install -r requirements.txt

