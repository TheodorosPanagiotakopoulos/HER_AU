# Slow Growth Method Simulations

This repository contains scripts to parse, process, and generate movies from slow growth method simulations.

## Directory Structure

Your simulations should be organized in the following format:

```
root_directory/
│── slow_growth_method/
│   ├── parse.py
│   ├── generate_image.py
│   ├── get_data.py
│   ├── get_mols.py
│   ├── add_time.py
│   ├── movie.py
│── RUN1/
│── RUN2/
│── ...
│── RUNN/
│── automate_movie.py
│── run.sh
```

## Instructions

### 1. Run `parse.py`  
Navigate to the `slow_growth_method` directory and execute:

```bash
python parse.py
```

This script processes the simulation data from the RUNX folders.

### 2. Generate Images

In the `slow_growth_method` directory, run:

```bash
python generate_image.py
```

Make sure `get_data.py` and `get_mols.py` are present in the same directory, as they are required for generating images.

### 3. Add Time and Process Movie

Ensure that `add_time.py` and `movie.py` are present in the `slow_growth_method` directory. These scripts are needed for the next steps.

### 4. Run `automate_movie.py`

Move one folder back, where all the slow growth method simulations (`RUN1`, `RUN2`, ...) are stored. Then, execute:

```bash
python automate_movie.py
```

This script automates the movie creation process. It requires `run.sh` in the same directory to function properly.

### 5. Modify the Output Path

Before running `automate_movie.py`, ensure that the path where the generated video should be saved is correctly set in the script.

## Requirements

Ensure you have the necessary dependencies installed. If needed, install them using:

```bash
pip install -r requirements.txt
```

(If `requirements.txt` is missing, manually install any required libraries.)

This setup ensures a smooth workflow from parsing simulation data to generating a final movie. Happy simulating! 🚀

