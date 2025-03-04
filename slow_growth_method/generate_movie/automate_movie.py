import os
import re
import sys
import glob
import subprocess


# Checks if any PNG files matching the pattern "frame*.png" exist in the given simulation directory.
# path_to_SG_simulation: The path to the directory where the simulation frames are stored.
# Returns: True if at least one matching PNG file is found, otherwise False.
def any_frame_png_exists( path_to_SG_simulation ):
	return bool( glob.glob( os.path.join( path_to_SG_simulation , "frame*.png" ) ) )

# Retrieves a list of directories in the current working directory that start with a given prefix.
# prefix: A string that the directory names should start with. Default is an empty string, which matches all directories.
# Returns: A list of directory names matching the specified prefix.
def get_dirs( prefix = "" ):
    path = os.getcwd()
    return [ i for i in os.listdir( path ) if i.startswith( prefix ) ]

# Modifies a specific pattern in a given file by replacing its assigned value.
# file: The name or path of the file to be modified. Default is "movie.py".
# pattern: The pattern to search for in the file. Default is 'output = "'.
# change_to: The new value to replace the existing one within the pattern. Default is an empty string.
# Returns: None. The function directly modifies the file in place.
def change_pattern( file = "movie.py", pattern='output = "', change_to="" ):
	file_path = file
	with open( file_path, "r" ) as file:
		lines = file.readlines()

	with open( file_path, "w" ) as file:
		for line in lines:
			if line.strip().startswith( pattern ):
				line = re.sub( r'output = ".*?"', 'output = "' + change_to + '"', line )
			file.write( line )

# Updates movie output filenames and manages simulation directories.
# path_to_cations_dir: The path to the cations directory containing simulation data.
# The function:
# - Extracts a prefix from the directory path.
# - Identifies relevant subdirectories matching the prefix.
# - Removes outdated output files from each directory.
# - Checks for existing frame PNG files and updates the output filename pattern accordingly.
# - Executes the simulation script.
# - Copies the generated MP4 files to a designated location.
# Returns: None. The function modifies files and directories in place.
def change_movie_output( path_to_cations_dir ):
	prefix = path_to_cations_dir.split( "/" )[ -2 ]
	print( "prefix = ", prefix )
	dirs = get_dirs( prefix )
	number_of_cations =  path_to_cations_dir.split( "/" )[ -2 ][ 0 ] 
	mapping = path_to_cations_dir.split( "/" )[ -1 ]
	h_map = { "H2O_splitting_from_NH4_hydration_shell" : "_NH4_hyd_cell_", "H2O_splitting_NOT_from_NH4_hydration_shell" : "_NH4_NO_hyd_cell_", "NH4_splitting" : "_splitting", "shuttling" : "_shuttling" }
	for i in dirs:
		new_path = os.path.join( path_to_cations_dir, i )
		os.chdir( new_path )
		os.system( "rm -rf outdir *.mp4" )
		ending = new_path.split( "/" )[ -1 ].split( "_" )[ -1 ] 
		frames = any_frame_png_exists( new_path )
		if frames:
			print( new_path )
			new_pattern = "MD_" + str( number_of_cations ) + h_map[ mapping ] + ending
			change_pattern( change_to = new_pattern )
			print( "frames found" )
		else:
			print( new_path )
			print( "No frames found, moving to the next directory" )
	os.chdir( path_to_cations_dir )
	os.system( "./run.sh" )
	for i in dirs:
		new_path = os.path.join( path_to_cations_dir, i )
		os.chdir( new_path )
		os.system( "cp *.mp4 /home/theodoros/PROJ_ElectroCat/theodoros/HER/Au/HER_Au/slow_grow_method/NH4/movies_NH4" )

# Generates movies by updating and processing simulation output files.
# path_to_cations_dir: The path to the directory containing cation simulation data.
# The function:
# - Calls change_movie_output to update filenames and manage simulation outputs.
# - Ensures the necessary modifications are applied before generating the movie files.
# Returns: None. The function modifies files and directories in place.
def generate_movie( path_to_cations_dir ):
	change_movie_output( path_to_cations_dir )


if __name__ == "__main__" :
	path = os.getcwd()
	generate_movie( path_to_cations_dir = path )
