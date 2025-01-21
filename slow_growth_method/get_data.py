import os
import gzip
import glob
import sys
import pandas as pd
import numpy as np
import json
import getpass
from ase.io import read

# Replaces '~' in the given path with the user's home directory (/home/<username>)
# path: The file path that may start with '~'
def convert( path ):
	username=getpass.getuser()
	return path.replace( '~', '/home/' + username )

# Splits a file path into its components using "/" as the delimiter.
# path_to_SG_calculation: The full file path to be split into individual directory or file components.
def split_path( path_to_SG_calculation ):
        return path_to_SG_calculation.split( "/" )

# Sorts a dictionary by its values in ascending order and returns a new sorted dictionary.
# dictionary: The dictionary to be sorted, with keys and their corresponding values.
def sort_dict( dictionary ):
	keys = list( dictionary.keys() )
	values = list( dictionary.values() )
	sorted_value_index = np.argsort( values )
	sorted_dict = { keys[ i ]: values[ i ] for i in sorted_value_index }
	return sorted_dict

# Flattens a 2D matrix (list of lists) into a single list containing all elements in row-wise order.
# matrix: A 2D list (matrix) to be flattened into a 1D list.
def flatten_matrix( matrix ):
	flat_list = list()
	for i in matrix:
		flat_list.extend( i )
	return flat_list

# Retrieves a list of all directories or files in the current directory that match the pattern "RUN*".
# Returns: A list of matching items or an empty list if no matches are found.
def get_RUNs():
		runs = glob.glob( "RUN*" )
		if runs:
			return runs
		else:
			return []

# Extracts the "H" index from the first line of an ICONST file.
# iconst: The path to the ICONST file to process.
# to_print: A flag (string) to control whether the "H index" is printed ("True" for printing, default is "False").
# Returns: The "H" index as an integer (adjusted to zero-based indexing).
def get_H_from_ICONST( iconst, to_print = "False"):
	last_line = list()
	with open( iconst ) as file:
		lines = [ line.rstrip() for line in file ]
	last_line.append( lines[ 0 ] )
	H_idx = last_line[ 0 ].split( " " )[ -2 ]
	if to_print == "True":
		print( "H index: ", H_idx )
	return int( H_idx ) - 1

# Calculates the minimum initial distance between the "H" atom and any "Au" atom in the system.
# path_to_SG_calculation: The path to the SG calculation directory.
# to_print: A flag (string) to control whether the minimum distance is printed ("True" for printing, default is "False").
# Returns: The minimum distance between the "H" atom and any "Au" atom, rounded to 2 decimal places.
def get_initial_H_Au_distance( path_to_SG_calculation, to_print = "False" ):
	distance_H_to_Au = list()
	os.chdir( path_to_SG_calculation )
	runs = get_RUNs()
	H_idx = get_H_from_ICONST( "ICONST" )
	if not runs:
		system = read( "POSCAR" )
	else:
		system = read( "RUN1/POSCAR" )
	au_indices = [ i for i, j in enumerate( system ) if j.symbol == "Au" ]
	for au_idx in au_indices:
		distance_H_to_Au.append( np.linalg.norm( system.positions[ au_idx ] - system.positions[ H_idx ] ) )
	min_H_Au_dist = round( min( distance_H_to_Au ), 2 )
	if to_print == "True":
		print( "min_H_Au_dist: ", min_H_Au_dist )	
	return min_H_Au_dist

################################# SLOW GROWTH METHOD #################################

# Processes a "REPORT" file or "REPORT.gz" file to extract and parse lines containing "cc" and "b_m".
# Returns: Two lists extracted from the file - one containing specific columns from lines with "cc" and another with "b_m".
def get_cc_bm():
	files = glob.glob( "REPORT*" )
	if "REPORT" in files:
		print( "REPORT found" )
		with open( files[ 0 ], "rb" ) as file:
			lines = file.readlines()
		cc = [ line.decode( "utf-8", errors = "ignore" ).strip() for line in lines if b"cc" in line ]	
		b_m = [ line.decode( "utf-8", errors = "ignore" ).strip() for line in lines if b"b_m" in line ]
	elif "REPORT.gz" in files:
		with gzip.open( files[ 0 ], "rb" ) as file:
			lines = file.readlines()
		cc = [ line.decode( "utf-8", errors = "ignore" ).strip() for line in lines if b"cc" in line ]	
		b_m = [ line.decode( "utf-8", errors = "ignore" ).strip() for line in lines if b"b_m" in line ]
	else:
		print( "REPORT NOT found" )
	df_cc = pd.DataFrame( [ row.split() for row in cc ] )
	df_bm = pd.DataFrame( [ row.split() for row in b_m ] )
	#print( df_bm.to_string() )
	return list( df_cc[ 3 ] ), list( df_bm[ 1 ] )

# Collects and processes "cc" and "b_m" data from a specified directory containing SG calculations.
# path_to_SG_calculation: The path to the directory containing SG calculation data.
# Returns: Two lists of floats - one for "cc" and another for "b_m". Returns (None, None) if no valid data is found.
def collect_cc_and_bm( path_to_SG_calculation ):
	if os.path.exists( path_to_SG_calculation ):
		os.chdir( path_to_SG_calculation )
		runs = get_RUNs()
		if not runs and os.path.isfile( path_to_SG_calculation + "/OUTCAR" ) == False:
			return None, None
		if not runs:
			print( "No RUN directories" )
			CC, B_M = get_cc_bm()
			return [ float( i ) for i in  CC ], [ float( i ) for i in  B_M ]
		else:
			print( runs )
			CC = list()
			B_M = list()
			for i in range( 0, len( runs ) ):
				new_path = path_to_SG_calculation + "/RUN" + str( i + 1 )
				os.chdir( new_path )
				cc, b_m = get_cc_bm()
				CC.append( cc )
				B_M.append( b_m )
			os.chdir( path_to_SG_calculation )
			cc, b_m = get_cc_bm()
			CC.append( cc )
			B_M.append( b_m )
			return [ float( i ) for i in flatten_matrix( CC ) ], [ float( i ) for i in flatten_matrix( B_M ) ]
	else:
		print( "Path not found" )

# Computes the free energy (tg) based on the "cc" and "b_m" data from a specified SG calculation directory.
# path_to_SG_calculation: The path to the directory containing SG calculation data.
# Returns: Two lists - one for "cc" (cumulative charge) and another for "tg" (free energy). Returns (None, None) if no valid data is found.
def get_free_energy( path_to_SG_calculation ):
	tg = [ 0.0 ]
	cc, b_m = collect_cc_and_bm( path_to_SG_calculation )
	if ( not cc ) or ( not b_m ):
		return None, None
	for i in range( 1, len( cc ) ):
		gg = 0.5 * ( cc[ i ]  -  cc[ i - 1 ] ) * ( b_m[ i ]  +  b_m[ i - 1 ] )
		tg.append( tg[ -1 ] + gg )
	return cc, tg

# Calculates the energy barrier from the free energy data of an SG calculation directory.
# path_to_SG_calculation: The path to the directory containing SG calculation data.
# Returns: The energy barrier (rounded to 2 decimal places) or None if free energy data is unavailable.
def get_barrier( path_to_SG_calculation ):
	#do not allow printing from get free_energy() 
    original_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    try:
        x, y = get_free_energy( path_to_SG_calculation )
        if x == None or y == None:
            return None
        else:
            barrier = round( max( y ) - min( y ), 2 )
            return barrier
    finally:
            sys.stdout.close()
            sys.stdout = original_stdout

# Updates a dictionary with energy barrier values for a specific SG calculation path.
# path_to_SG_calculation: The path to the SG calculation directory.
# barriers_dict: The dictionary to store barrier values, keyed by the directory path. It is initally empty
# Returns: The updated dictionary with the barrier value for the specified path.
def get_barriers_to_dictionary( path_to_SG_calculation, barriers_dict ):
	os.chdir( path_to_SG_calculation )
	path_to_list = split_path( path_to_SG_calculation )
	barrier = get_barriers( path_to_SG_calculation )
	key = path_to_list[ -3 ] + "/" + path_to_list[ -2 ] +  "/" + path_to_list[ -1 ]
	if barrier is None:
		barriers_dict[ key ] = "Not started yet"
	else:
		barriers_dict[ key ] = barrier
	#print( barriers_dict )
	return barriers_dict

################################# FROM DATABASE #################################

# Loads a JSON database from a specified file.
# database_file: The path to the JSON file containing the database.
# Returns: The loaded database as a Python dictionary.
def load_database( database_file ):
	with open( database_file, "r" ) as file:
		database = json.load( file )
	return database

# Extracts and filters barriers from a database based on a specific key, only considering entries marked as "Good".
# database: The database containing barrier information.
# val: The key in the database to filter and process.
# Returns: A sorted dataframe of barriers with keys as "barrier_<name>" and values as the barrier values.
def get_barrier_from_db( database, val, to_print = False ):
	data = pd.DataFrame()
	filtered_data = {}
	names = list()
	barriers = list()
	distances = list()
	for key, value in database[ val ].items():
		#aux_key =""
		if value[ "note" ] == "Good":
			#print( value[ "path" ] )
			if value[ "path" ].split("/")[ -2 ] == "H2O_from_hydration_shell_splitting" or value[ "path" ].split("/")[ -2 ] == "H2O_splitting_from_CH3NH3_hydration_shell":
				aux_key = "hyd_shell"
			elif value[ "path" ].split("/")[ -2 ] == "free_H2O_splitting" or value[ "path" ].split("/")[ -2 ] == "H2O_splitting_NOT_from_CH3NH3_hydration_shell":
				aux_key = "NO_hyd_shell"
			elif value[ "path" ].split("/")[ -2 ] == "CH3NH3_splitting":
				aux_key = value[ "path" ].split("/")[ -2 ]
			elif value[ "path" ].split("/")[ -2 ] == "shuttling":
				aux_key = value[ "path" ].split("/")[ -2 ] + "_shuttling"
			filtered_data[  "bar_" + aux_key + "_" + value[ "path" ].split("/")[ -1 ] ] = get_barrier( convert( value[ "path" ] ) )
			filtered_data[  "Dist(H-Au)_" + value[ "path" ].split("/")[ -1 ] ] = get_initial_H_Au_distance( convert( value[ "path" ] ) )
	sorted_dict = sort_dict( filtered_data )
	for key in sorted_dict:
		if key.startswith( "bar" ):
			names.append( key )
			barriers.append( sorted_dict[ key ] )
		else:
			distances.append( sorted_dict[ key ] )
	data[ "CONF" ] = names
	data[ "barrier" ] = barriers
	data[ "H-Au_distance" ] = distances
	if to_print == True:
		print( data.to_string() )
	return data 

if __name__ == "__main__":
	path = "/home/theodoros/PROJ_ElectroCat/theodoros/HER/Au/HER_Au/database/"
	data = load_database( path + "database_for_theo.js" )

	#Na_5_hyd = get_barrier_from_db( data, "5_Na_H2O_dissociation_from_hydration_shell", to_print = True )

	#Na_5_No_hyd = get_barrier_from_db( data, "5_Na_H2O_dissociation_NOT_from_hydration_shell", to_print = True )
	
	#Na_3_hyd = get_barrier_from_db( data, "3_Na_H2O_dissociation_from_hydration_shell", to_print = True )

	#Na_3_No_hyd = get_barrier_from_db( data, "3_Na_H2O_dissociation_NOT_from_hydration_shell", to_print = True )
	
	#Na_1_hyd = get_barrier_from_db( data, "1_Na_H2O_dissociation_from_hydration_shell", to_print = True )

	#Na_1_No_hyd = get_barrier_from_db( data, "1_Na_H2O_dissociation_NOT_from_hydration_shell", to_print = True )

	CH3NH3_5_hyd = get_barrier_from_db( data, "5_CH3NH3_H2O_dissociation_from_hydration_shell", to_print = True )
	
	CH3NH3_5_NO_hyd = get_barrier_from_db( data, "5_CH3NH3_H2O_dissociation_NOT_from_hydration_shell", to_print = True )
	
	CH3NH3_5_splitting = get_barrier_from_db( data, "5_CH3NH3_spliting", to_print = True )

	CH3NH3_5_shuttling = get_barrier_from_db( data, "5_CH3NH3_shuttling", to_print = True )
