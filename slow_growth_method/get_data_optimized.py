import os
import gzip
import glob
import sys
import pandas as pd
import numpy as np
import json
import getpass
from ase.io import read
import get_mols

def get_element_indices( system, element, verbose = False ):
	return  [ i for i, j in enumerate( system ) if j.symbol == element ]

# Replaces '~' in the given path with the user's home directory (/home/<username>)
# path: The file path that may start with '~'
def convert( path ):
	username = getpass.getuser()
	return path.replace( '~', '/home/' + username )

# Splits a file path into its components using "/" as the delimiter.
# path_to_SG_calculation: The full file path to be split into individual directory or file components.
def split_path( path_to_SG_calculation ):
        return path_to_SG_calculation.split( "/" )

# Sorts a dictionary by its values in ascending order and returns a new sorted dictionary.
# dictionary: The dictionary to be sorted, with keys and their corresponding values.
def sort_dict( dictionary ):
	return dict( sorted( dictionary.items(), key = lambda item: item[ 1 ] ) )

# Flattens a 2D matrix (list of lists) into a single list containing all elements in row-wise order.
# matrix: A 2D list (matrix) to be flattened into a 1D list.
def flatten_matrix( matrix ):
	return[ element for row in matrix for element in row ]

# Retrieves a list of all directories or files in the current directory that match the pattern "RUN*".
# Returns: A list of matching items or an empty list if no matches are found.
def get_RUNs( path_to_simulation ):
	path_to_simulation = os.path.normpath( path_to_simulation )
	runs = glob.glob( os.path.join( path_to_simulation, "RUN*" ) )
	runs.sort( key = lambda x: int( os.path.basename( x ).replace( "RUN", "" ) ) )
	return runs

# Extracts the "H" index from the first line of an ICONST file.
# iconst: The path to the ICONST file to process.
# to_print: A flag (string) to control whether the "H index" is printed ("True" for printing, default is "False").
# Returns: The "H" index as an integer (adjusted to zero-based indexing).
def get_H_from_ICONST( iconst, verbose = False ):
	first_line = list()
	with open( iconst ) as file:
		lines = [ line.rstrip() for line in file ]
	first_line.append( lines[ 0 ] )
	H_idx = first_line[ 0 ].split( " " )[ -2 ]
	if verbose:
		print( "H index: ", H_idx )
	return int( H_idx ) - 1

# Parses an ICONST file to extract atomic indices based on the number of lines in the file.
# iconst: Path to the ICONST file.
# to_print: If True, prints detailed index information. Defaults to False.
# Returns: 
# - If the ICONST file has 2 lines, returns a tuple (H_idx, O_idx) with zero-based indices.
# - If the ICONST file has 3 lines, returns a tuple (O_H2O_idx, H_H2O_idx, H_cation_idx) with zero-based indices.
# Raises: ValueError if the ICONST file contains less than 2 lines or more than 3 lines.
def get_data_ICONST( iconst, verbose = False ):
	first_line = list()
	second_line = list()
	with open( iconst ) as file:
		lines = [ line.rstrip() for line in file ]
	if len( lines ) < 2:
		raise ValueError( "ICONST file must contain at least two lines." )
	elif len( lines ) == 2:
		with open( iconst ) as file:
			lines = [ line.rstrip() for line in file ]
		first_line =  lines[ 0 ].split( " " ) 
		O_idx = first_line[ 1 ]
		H_idx = first_line[ 2 ]
		if verbose: 
			print( "ICONST H index: ", H_idx )
			print( "Actual H index: ", int( H_idx ) - 1 ) 
			print( "ICONST O index: ", O_idx )
			print( "Actual O index: ", int( O_idx ) - 1 )
		return int( H_idx ) - 1, int( O_idx ) - 1
	elif len( lines ) == 3:
		with open( iconst ) as file:
			lines = [ line.rstrip() for line in file ]
		first_line = lines[ 0 ].split( " " )
		second_line = lines[ 1 ].split( " " )
		O_H2O_idx = first_line[ 1 ]
		H_H2O_idx = first_line[ 2 ]
		H_cation_idx = second_line[ 2 ]
		if verbose:
			print( "ICONST O(H2O) index: ", O_H2O_idx )
			print( "ICONST H(H2O) index: ", H_H2O_idx )
			print( "ICONST H(cation) index: ", H_cation_idx )
			print( "Actual O(H2O) index: ", int( O_H2O_idx ) -1 )
			print( "Actual H(H2O) index: ", int( H_H2O_idx ) -1 )
			print( "Actual H(cation) index: ", int( H_cation_idx ) -1 )
		return int( O_H2O_idx ) -1, int( H_H2O_idx ) -1, int( H_cation_idx ) -1
	else:
		raise ValueError( "ICONST file must NOT contain more than three lines." )

def get_standarized_ICONST_data( iconst_data ):
	if len( iconst_data ) == 2:
		H_idx, O_idx = iconst_data
		return H_idx, O_idx, np.nan
	elif len( iconst_data ) == 3:
		O_idx, H_idx, H_cation_idx = data_ICONST
		return H_idx, O_idx, H_cation_idx
	else:
		raise ValueError( "ICONST file must have 2 or 3 lines" )

# Calculates the minimum initial distance between the "H" atom and any "Au" atom in the system.
# path_to_SG_calculation: The path to the SG calculation directory.
# to_print: A flag (string) to control whether the minimum distance is printed ("True" for printing, default is "False").
# Returns: The minimum distance between the "H" atom and any "Au" atom, rounded to 2 decimal places.
def get_initial_H_Au_distance( path_to_SG_calculation, verbose = False ):
	distance_H_to_Au = list()
	os.chdir( path_to_SG_calculation )
	runs = get_RUNs( path_to_SG_calculation )
	data_ICONST = get_data_ICONST( "ICONST" )
	H_idx, O_idx, H_cation_idx = get_standarized_ICONST_data( data_ICONST )	
	if not runs:
		system = read( "POSCAR" )
	else:
		system = read( "RUN1/POSCAR" )
	au_indices = get_element_indices( system, "Au" )
	for au_idx in au_indices:
		distance_H_to_Au.append( np.linalg.norm( system.positions[ au_idx ] - system.positions[ H_idx ] ) )
	min_H_Au_dist = round( min( distance_H_to_Au ), 2 )
	if verbose:
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
		runs = get_RUNs( path_to_SG_calculation )
		if not runs and not (os.path.isfile(path_to_SG_calculation + "/OUTCAR") or os.path.isfile(path_to_SG_calculation + "/OUTCAR.gz")):
			return None, None
		#if not runs and os.path.isfile( path_to_SG_calculation + "/OUTCAR" ) == False:
		#	return None, None
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

################################# FROM ICONST #################################

# Retrieves the initial atomic structure system for a given simulation path.
# path_to_simulation: The path to the simulation directory.
# Returns: 
# - The initial atomic structure as an ASE atoms object.
# - If no `RUN*` directories are present, reads the system from `POSCAR` in the base path.
# - If `RUN*` directories are present, reads the system from `POSCAR` in the `RUN1` directory.
def get_initial_system( path_to_simulation ):
	runs = get_RUNs( path_to_simulation )
	path = convert( path_to_simulation )
	if not runs:
		initial_system =  path +  "/POSCAR" 
	else:
		initial_system =  path + "/RUN1/POSCAR" 
	return initial_system 

def get_O_cation_min_distance(cation, path, O_idx, cation_list):
	runs = get_RUNs(path)
	if not runs:
		system = read(path + "/POSCAR")
	else:
		first_run = runs[0].split("/")[-1]
		system = read(path + "/" + first_run + "/POSCAR")

	O_position = system.positions[O_idx]
	distances = list()

	if cation == "Na":
		for cation_idx in cation_list:
			cation_position = system.positions[cation_idx]
			distances.append(np.linalg.norm(O_position - cation_position))
		return round(min(distances), 3), np.nan, np.nan

	elif cation in ["N-NH4", "N-CH3NH3"]:
		cation_distances = list()
		H_distances = list()
		H_bond_info = ""

		closest_H_distance = float("inf")  # Track the global closest H
		closest_H_idx = None  # Track the index of the closest H

		for cation_group in cation_list:
			N_idx = cation_group[0]
			N_position = system.positions[N_idx]
			cation_distances.append(np.linalg.norm(O_position - N_position))

			H_indices = cation_group[1:]

			for H_idx in H_indices:
				H_position = system.positions[H_idx]
				H_distance = np.linalg.norm(O_position - H_position)
				H_distances.append(H_distance)

				# Update global closest H
				if H_distance < closest_H_distance:
					closest_H_distance = H_distance
					closest_H_idx = H_idx

		# Only set H_bond_info if we actually found a valid H index
		if closest_H_idx is not None:
			H_bond_info = f"{O_idx}-{closest_H_idx}"

		return round(min(cation_distances), 3), round(closest_H_distance, 3), H_bond_info

	else:
		raise ValueError("Unsupported cation type. Supported types are: 'Na', 'N-NH4', 'N-CH3NH3'.")



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
# Returns: A sorted dictionary of barriers with keys as "barrier_<name>" and values as the barrier values.
def get_barrier_from_db( database, val, fixed_length = 45, to_print = False ):
	data = pd.DataFrame()
	filtered_data = {}
	names = list()
	barriers = list()
	distances_H_Au = list()
	distances_O_cation = list()
	distances_O_H = list()
	status = list()
	bond_data = list()
	ICONST_indices = list()
	aux_key = None
	suggetion = list()

	for key, value in database[ val ].items():
		if value[ "note" ] == "Good" or value[ "note" ] == "Bad":
			stat = value[ "note" ] #new
			status.append( stat )
			path_key = value[ "path" ].split( "/" )[ -2 ]
			if path_key in [
				"H2O_from_hydration_shell_splitting",
				"H2O_splitting_from_NH4_hydration_shell",
				"H2O_splitting_from_CH3NH3_hydration_shell"
			]:
				aux_key = "hyd_shell"
			elif path_key in [
				"free_H2O_splitting",
				"H2O_splitting_NOT_from_NH4_hydration_shell",
				"H2O_splitting_NOT_from_CH3NH3_hydration_shell"
			]:
				aux_key = "NO_hyd_shell"
			elif path_key in ["NH4_splitting", "CH3NH3_splitting"]:
				aux_key = path_key
			elif path_key == "shuttling":
				aux_key = f"{path_key}_shuttling"
			else:
				aux_key = ""

			if not aux_key:
				print("Nothing to report here")
				return None
			
			initial_system = get_initial_system( convert( value[ "path" ] ) )
			if "Na" in f"bar_{aux_key}_{value['path'].split('/')[-1]}":
				cation = "Na"
				cations = get_mols.get_Na_mols(initial_system)
			elif "NH4" in f"bar_{aux_key}_{value['path'].split('/')[-1]}":
				cation = "N-NH4"
				cations = get_mols.get_NH4_mols(initial_system)
			elif "CH3NH3" in f"bar_{aux_key}_{value['path'].split('/')[-1]}":
				cation = "N-CH3NH3"
				cations = get_mols.get_CH3NH3_mols(initial_system)

			data_ICONST = get_data_ICONST( convert( value[ "path" ] ) + "/ICONST" )
			
			H_cation_idx = None
			if len( data_ICONST ) == 2:
				H_idx, O_idx = data_ICONST
				ICONST_idx = f"{O_idx}-{H_idx}"
			elif len( data_ICONST ) == 3:
				O_idx, H_H2O_idx, H_cation_idx = data_ICONST
				ICONST_idx = f"{O_idx}-{H_H2O_idx}-{H_cation_idx}"
			
			ICONST_indices.append( ICONST_idx )

			filtered_data[ f"bar_{aux_key}_{value['path'].split('/')[-1]}" ] = get_barrier(
				convert( value[ "path" ] )
			)
			filtered_data[ f"Dist(H-Au)_{value['path'].split('/')[-1]}" ] = get_initial_H_Au_distance(
				convert( value[ "path" ] )
			)
			filtered_data[ f"Dist(O-cation)_{aux_key}_{value['path'].split( '/' )[-1] }" ] = get_O_cation_min_distance(
				cation, convert( value[ "path" ] ), O_idx, cations
			)[ 0 ]
			filtered_data[ f"Dist(O-H)_{aux_key}_{value['path'].split( '/' )[-1] }" ] = get_O_cation_min_distance(
                cation, convert( value[ "path" ] ), O_idx, cations
            )[ 1 ]
			filtered_data[ f"Bond_data_{aux_key}_{value['path'].split( '/' )[-1] }" ] = get_O_cation_min_distance(
                cation, convert( value[ "path" ] ), O_idx, cations
            )[ 2 ]

	if aux_key is None:
		print( "No valid aux_key found in database for: ", val, "\n" )
		return None
	
	sorted_dict = sort_dict( filtered_data )
	

	for key in filtered_data:
		if key.startswith( "bar" ):
			names.append( key )
			barriers.append( sorted_dict[ key ] )
		elif key.startswith( "Dist(H-Au)_" ):
			distances_H_Au.append( sorted_dict[ key ] )
		elif key.startswith( "Dist(O-cation)_" ):
			distances_O_cation.append( sorted_dict[ key ] )
		elif key.startswith( "Dist(O-H)" ):
			distances_O_H.append( sorted_dict[ key ] )
		elif key.startswith( "Bond_data_" ):
			bond_data.append( sorted_dict[ key ] )

	data[ "CONF" ] = names
	data[ "barrier" ] = barriers
	data[ "H-Au_distance" ] = distances_H_Au
	if ( aux_key and "splitting" not in aux_key ):
		#data[ "O-cation_distance" ] = distances_O_cation
		data[ "O-H_distance" ] = distances_O_H
	
	data[ "O-H-index" ] = bond_data
	data[ "ICONST_idx" ] = ICONST_indices
	data[ "status" ] = status

	sorted_data = data.sort_values( by = "barrier" )

	sorted_data[ "CONF" ] = sorted_data[ "CONF" ].apply( lambda x: x[ : fixed_length ].ljust( fixed_length ) )
	sorted_data[ "propose" ] = np.nan	

	pd.set_option( "display.colheader_justify", "left" ) 
	sorted_data.style.set_properties( **{ "text-align": "center" } ) 

	for i, j in zip( sorted_data[ "O-H-index" ], sorted_data[ "ICONST_idx" ] ):
		if i.split( "-" )[ -1 ] != j.split( "-" )[ -1]:
			suggested_value =  str( i.split( "-" )[ 0 ] ) + " " + str( j.split( "-" )[ 1 ] ) + " " + str( i.split( "-" )[ -1] )
		else:
			suggested_value = np.nan
		suggetion.append( suggested_value )
		
	sorted_data[ "propose" ] = suggetion
	

	if to_print:
		print( sorted_data.to_string(), "\n" )

	return sorted_data


if __name__ == "__main__":
	path = "/home/theodoros/PROJ_ElectroCat/theodoros/HER/Au/HER_Au/database/"
	data = load_database( path + "database_for_theo.js" )

	#Na_1_hyd = get_barrier_from_db( data, "1_Na_H2O_dissociation_from_hydration_shell", to_print = True )

	#Na_1_No_hyd = get_barrier_from_db( data, "1_Na_H2O_dissociation_NOT_from_hydration_shell", to_print = True )	

	#Na_3_hyd = get_barrier_from_db( data, "3_Na_H2O_dissociation_from_hydration_shell", to_print = True )

	#Na_3_No_hyd = get_barrier_from_db( data, "3_Na_H2O_dissociation_NOT_from_hydration_shell", to_print = True )
	
	#Na_5_hyd = get_barrier_from_db( data, "5_Na_H2O_dissociation_from_hydration_shell", to_print = True )

	#Na_5_No_hyd = get_barrier_from_db( data, "5_Na_H2O_dissociation_NOT_from_hydration_shell", to_print = True )

	
	#NH4_1_hyd = get_barrier_from_db( data, "1_NH4_H2O_dissociation_from_hydration_shell", to_print = True )
	
	#NH4_1_NO_hyd = get_barrier_from_db( data, "1_NH4_H2O_dissociation_NOT_from_hydration_shell", to_print = True )
	
	#NH4_1_splitting = get_barrier_from_db( data, "1_NH4_spliting", to_print = True )

	NH4_1_shuttling = get_barrier_from_db( data, "1_NH4_shuttling", to_print = True )


	#NH4_3_hyd = get_barrier_from_db( data, "3_NH4_H2O_dissociation_from_hydration_shell", to_print = True )

	#NH4_3_NO_hyd = get_barrier_from_db( data, "3_NH4_H2O_dissociation_NOT_from_hydration_shell", to_print = True )

	#NH4_3_splitting = get_barrier_from_db( data, "3_NH4_spliting", to_print = True )

	#NH4_3_shuttling = get_barrier_from_db( data, "3_NH4_shuttling", to_print = True )


	#CH3NH3_1_hyd = get_barrier_from_db( data, "1_CH3NH3_H2O_dissociation_from_hydration_shell", to_print = True )

	#CH3NH3_1_NO_hyd = get_barrier_from_db( data, "1_CH3NH3_H2O_dissociation_NOT_from_hydration_shell", to_print = True )

	#CH3NH3_1_splitting = get_barrier_from_db( data, "1_CH3NH3_spliting", to_print = True )

	#CH3NH3_1_shuttling = get_barrier_from_db( data, "1_CH3NH3_shuttling", to_print = True )


	#CH3NH3_3_hyd = get_barrier_from_db( data, "3_CH3NH3_H2O_dissociation_from_hydration_shell", to_print = True )

	#CH3NH3_3_NO_hyd = get_barrier_from_db( data, "3_CH3NH3_H2O_dissociation_NOT_from_hydration_shell", to_print = True )

	#CH3NH3_3_splitting = get_barrier_from_db( data, "3_CH3NH3_spliting", to_print = True )

	#CH3NH3_3_shuttling = get_barrier_from_db( data, "3_CH3NH3_shuttling", to_print = True )


	#CH3NH3_5_hyd = get_barrier_from_db( data, "5_CH3NH3_H2O_dissociation_from_hydration_shell", to_print = True )
	
	#CH3NH3_5_NO_hyd = get_barrier_from_db( data, "5_CH3NH3_H2O_dissociation_NOT_from_hydration_shell", to_print = True )
	
	#CH3NH3_5_splitting = get_barrier_from_db( data, "5_CH3NH3_spliting", to_print = True )

	#CH3NH3_5_shuttling = get_barrier_from_db( data, "5_CH3NH3_shuttling", to_print = True )
