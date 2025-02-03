kimport os
import sys
import gzip
import glob
import json
import getpass
import numpy as np
import pandas as pd
from ase.io import read
import get_mols_updated as get_mols

# Retrieves the indices of all atoms of a specified element in the given atomic system.
# system: The atomic system (a list of atoms or an ASE Atoms object).
# element: The chemical symbol of the element whose indices are to be found.
# verbose: A flag (boolean) to control whether additional information is printed (default is False).
# Returns: A list of indices corresponding to the positions of the specified element in the system.
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
# verbose: A flag (string) to control whether the "H index" is printed ("True" for printing, default is "False").
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


################################# FROM ICONST #################################

# Parses an ICONST file to extract atomic indices based on the number of lines in the file.
# iconst: Path to the ICONST file.
# verbose: If True, prints detailed index information. Defaults to False.
# Returns: 
# - If the ICONST file has 2 lines, returns a tuple (H_idx, O_idx) with zero-based indices.
# - If the ICONST file has 3 lines, returns a tuple (O_H2O_idx, H_H2O_idx, H_cation_idx) with zero-based indices.
# Raises: ValueError if the ICONST file contains less than 2 lines or more than 3 lines.
def get_data_ICONST( path_to_SG_simulation, verbose = False ):
	iconst_path = os.path.join( path_to_SG_simulation, "ICONST" )
	with open( iconst_path ) as file:
		lines = [ line.rstrip() for line in file ]

	if len( lines ) < 2:
		raise ValueError( "ICONST file must contain at least two lines." )
	elif len(lines) > 3:
		raise ValueError( "ICONST file must NOT contain more than three lines." )

	first_line = lines[ 0 ].split()
	O_idx = int( first_line[ 1 ] ) - 1  
	H_idx = int( first_line[ 2 ] ) - 1  

	if verbose:
		print( "ICONST H index:", H_idx + 1 )
		print( "Actual H index:", H_idx )
		print( "ICONST O index:", O_idx + 1 )
		print( "Actual O index:", O_idx )

	if len( lines ) == 2:
		return int( H_idx ), int( O_idx )

	second_line = lines[ 1 ].split()
	H_cation_idx = int( second_line[ 2 ] ) - 1

	if verbose:
		print( "ICONST O(H2O) index:", O_idx + 1)
		print( "ICONST H(H2O) index:", H_idx + 1)
		print( "ICONST H(cation) index:", H_cation_idx + 1)
		print( "Actual O(H2O) index:", O_idx)
		print( "Actual H(H2O) index:", H_idx)
		print( "Actual H(cation) index:", H_cation_idx)
	return int( O_idx ), int( H_idx ), int( H_cation_idx )

# Standardizes the ICONST data by ensuring it follows a consistent format.
# iconst_data: A list containing either 2 or 3 indices representing atomic relationships.
# Returns: A tuple containing (H_idx, O_idx, H_cation_idx). If only two indices are provided, H_cation_idx is set to NaN.
# Raises: ValueError if the input does not contain exactly 2 or 3 elements.
def get_standarized_ICONST_data( iconst_data ):
	if len( iconst_data ) == 2:
		H_idx, O_idx = iconst_data
		return H_idx, O_idx, np.nan
	elif len( iconst_data ) == 3:
		O_idx, H_idx, H_cation_idx = iconst_data
		return H_idx, O_idx, H_cation_idx
	else:
		raise ValueError( "ICONST file must have 2 or 3 lines" )

# Calculates the minimum initial distance between the "H" atom and any "Au" atom in the system.
# path_to_SG_calculation: The path to the SG calculation directory.
# verbose: A flag (string) to control whether the minimum distance is printed ("True" for printing, default is "False").
# Returns: The minimum distance between the "H" atom and any "Au" atom, rounded to 2 decimal places.
def get_initial_H_Au_distance( path_to_SG_simulation, verbose = False ):
	distance_H_to_Au = list()
	os.chdir( path_to_SG_simulation )
	runs = get_RUNs( path_to_SG_simulation )
	data_ICONST = get_data_ICONST( path_to_SG_simulation )
	H_idx, O_idx, H_cation_idx = get_standarized_ICONST_data( data_ICONST )

	if not runs:
		system = read( "POSCAR" )
	else:
		system = read( "RUN1/POSCAR" )
	au_indices = get_element_indices( system, "Au" )

	for au_idx in au_indices:
		distance_H_to_Au.append( ( np.linalg.norm( system.positions[ au_idx ] - system.positions[ H_idx ] ), au_idx ) )

	min_H_Au_dist, Au_idx = min( distance_H_to_Au, key=lambda x: x[ 0 ] )
	min_H_Au_dist = round( min_H_Au_dist, 2 )

	if verbose:
		print("min_H_Au_dist:", min_H_Au_dist, "| Closest Au_idx:", Au_idx )

	return min_H_Au_dist,  str( H_idx ) + "-" + str( Au_idx )

# Retrieves the initial system configuration from a SG simulation.
# path_to_SG_simulation: The path to the SG simulation directory.
# Returns: The initial atomic structure read from the POSCAR file, either from the first run in the directory or directly from the main path if no runs are found.
def get_initial_system( path_to_SG_simulation ):
	runs = get_RUNs( path_to_SG_simulation )
	path = convert( path_to_SG_simulation )
	if not runs:
		return read( os.path.join( path_to_SG_simulation, "POSCAR" ) )
	else:
		return read( os.path.join( path_to_SG_simulation, "RUN1/POSCAR" ) )

# Gets the minimum Euclidean distance from a reference position.
# reference_position: NumPy array of the reference point.
# positions: List or NumPy array of positions.
# Returns: Minimum distance, rounded to 2 decimals.
def get_min_distance( reference_position, positions ):
	distances = [ np.linalg.norm( reference_position - pos ) for pos in positions ]
	return round( min( distances ), 2 )

# Gets the minimum distances between oxygen and cation/hydrogen.
# path_to_SG_simulation: Path to SG simulation.
# cation_type: "Na", "N-NH4", or "N-CH3NH3".
# verbose: Prints debug info if True.
# Returns: (O–cation distance, O–H distance, H-bond info).
def get_distances( path_to_SG_simulation, cation_type, verbose = False ):
	system = get_initial_system(path_to_SG_simulation )

	iconst_data = get_data_ICONST(path_to_SG_simulation, verbose)
	if len(iconst_data) == 2:
		H_idx, O_idx = iconst_data
		H_cation_idx = np.nan
	else:
		O_idx, H_idx, H_cation_idx = iconst_data

	O_position = system.positions[O_idx]
	closest_H_distance = float("inf")
	closest_H_idx = None

	runs = get_RUNs( path_to_SG_simulation )
	if cation_type == "Na":
		if not runs:
			cation_list = get_mols.get_Na_mols(path_to_SG_simulation)
		else:
			cation_list = get_mols.get_Na_mols(path_to_SG_simulation + "/RUN1")
	elif cation_type == "N-NH4":
		if not runs:
			cation_list = get_mols.get_NH4_mols(path_to_SG_simulation)
		else:
			cation_list = get_mols.get_NH4_mols(path_to_SG_simulation + "/RUN1")
	elif cation_type == "N-CH3NH3":
		if not runs:
			cation_list = get_mols.get_CH3NH3_mols(path_to_SG_simulation)
		else:
			cation_list = get_mols.get_CH3NH3_mols(path_to_SG_simulation + "/RUN1")
	else:
		raise ValueError( "Unsupported cation type. Supported types are: 'Na', 'N-NH4', 'N-CH3NH3'." )

	central_atom_positions = [ system.positions[ cation_group[ 0 ] ] for cation_group in cation_list ]
	min_cation_distance = get_min_distance( O_position, central_atom_positions )
	if cation_type in [ "N-NH4", "N-CH3NH3" ]:
		H_positions = list()
		H_indices_list = list()
		for cation_group in cation_list:
			for H_idx in cation_group[ 1: ]:
				H_positions.append( system.positions[ H_idx ] )
				H_indices_list.append( H_idx )

		if H_positions:
			distances = [ np.linalg.norm( O_position - H_pos ) for H_pos in H_positions ]
			closest_idx_in_list = int( np.argmin( distances ) )
			closest_H_distance = distances[ closest_idx_in_list ]
			closest_H_idx = H_indices_list[ closest_idx_in_list ]
		else:
			closest_H_distance = np.nan
			closest_H_idx = np.nan

	if cation_type == "Na":
		return round( min_cation_distance, 2 ), np.nan, np.nan
	else:
		H_bond_info = f"{O_idx}-{closest_H_idx}" if closest_H_idx is not None else np.nan
		return round( min_cation_distance, 2 ), round( closest_H_distance, 2 ), H_bond_info

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
	return list( df_cc[ 3 ] ), list( df_bm[ 1 ] )

# Collects and processes "cc" and "b_m" data from a specified directory containing SG calculations.
# path_to_SG_calculation: The path to the directory containing SG calculation data.
# Returns: Two lists of floats - one for "cc" and another for "b_m". Returns (None, None) if no valid data is found.
def collect_cc_and_bm( path_to_SG_calculation ):
	if os.path.exists( path_to_SG_calculation ):
		os.chdir( path_to_SG_calculation )
		runs = get_RUNs( path_to_SG_calculation )
		if not runs and not ( os.path.isfile( path_to_SG_calculation + "/OUTCAR" ) or os.path.isfile( path_to_SG_calculation + "/OUTCAR.gz" ) ):
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
    sys.stdout = open( os.devnull, "w" )
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
	return barriers_dict

################################# FROM DATABASE #################################

# Loads a JSON database from a specified file.
# database_file: The path to the JSON file containing the database.
# Returns: The loaded database as a Python dictionary.
def load_database( database_file ):
	with open( database_file, "r" ) as file:
		database = json.load( file )
	return database

# Determines the auxiliary key based on the provided path key.
# path_key: A string representing the specific path key.
# Returns: A corresponding auxiliary key as a string, or None if no match is found.
def get_aux_key( path_key ):
	if path_key in [ "H2O_from_hydration_shell_splitting", "H2O_splitting_from_NH4_hydration_shell", "H2O_splitting_from_CH3NH3_hydration_shell" ]:
		return "hyd_shell"
	elif path_key in [ "free_H2O_splitting", "H2O_splitting_NOT_from_NH4_hydration_shell", "H2O_splitting_NOT_from_CH3NH3_hydration_shell" ]:
		return "NO_hyd_shell"
	elif path_key in [ "NH4_splitting", "CH3NH3_splitting" ]:
		return path_key
	elif path_key == "shuttling":
		return f"{path_key}_shuttling"
	else:
		return None

# Processes a single database entry and extracts relevant data.
# value: A dictionary containing details about a specific database entry.
# filtered_data: Unused parameter in the function, can be removed if unnecessary.
# Returns: A dictionary of processed data and the ICONST index.
def process_database_entry( value, filtered_data ):
	path_key = value[ "path" ].split( "/" )[ -2 ]
	aux_key = get_aux_key( path_key )
	if not aux_key:
		print( "Nothing to report here" )
		return {}, None

	path_to_SG_simulation = convert( value[ "path" ] )
	initial_system = get_initial_system( path_to_SG_simulation )
	data_ICONST = get_data_ICONST( path_to_SG_simulation, verbose = False )

	if len( data_ICONST ) == 2:
		H_idx, O_idx = data_ICONST
		ICONST_idx = f"{O_idx}-{H_idx}"
	elif len( data_ICONST ) == 3:
		O_idx, H_H2O_idx, H_cation_idx = data_ICONST
		ICONST_idx = f"{O_idx}-{H_H2O_idx}-{H_cation_idx}"
	else:
		raise ValueError( "Invalid ICONST data format." )

	if "Na" in f"bar_{aux_key}_{value[ 'path' ].split( '/' )[ -1 ] }":
		cation = "Na"
	elif "NH4" in f"bar_{aux_key}_{value[ 'path' ].split( '/' )[ -1 ] }":
		cation = "N-NH4"
	elif "CH3NH3" in f"bar_{aux_key}_{value[ 'path' ].split( '/' )[ -1 ] }":
		cation = "N-CH3NH3"
	else:
		raise ValueError( "Unsupported cation type." )

	min_cation_distance, closest_H_distance, H_bond_info = get_distances( path_to_SG_simulation, cation )
	min_H_Au_dist, H_Au_idx = get_initial_H_Au_distance( path_to_SG_simulation )

	updated_data = {
		f"bar_{aux_key}_{value['path'].split( '/' )[ -1 ] }": get_barrier(path_to_SG_simulation),
		f"Dist(H-Au)_{value['path'].split( '/' )[ -1 ] }": min_H_Au_dist,
		f"Index(H-Au)_{value[ 'path' ].split( '/' )[ -1 ] }": H_Au_idx,
		f"Dist(O-cation)_{aux_key}_{value[ 'path' ].split( '/' )[ -1 ] }": min_cation_distance,
		f"Dist(O-H)_{aux_key}_{value[ 'path' ].split( '/' )[ -1 ] }": closest_H_distance,
		f"Bond_data_{aux_key}_{value[ 'path' ].split( '/' )[ -1 ] }": H_bond_info,
	}

	return updated_data, ICONST_idx

# Creates a pandas DataFrame from the filtered data.
# filtered_data: A dictionary containing processed data.
# status: A status label to be added to the DataFrame.
# aux_key: A key to determine if specific data should be included.
# Returns: A sorted pandas DataFrame based on barrier values.
def create_dataframe( filtered_data, status, aux_key ):
	data = pd.DataFrame()
	names, barriers, distances_H_Au, index_H_Au, distances_O_cation, distances_O_H, bond_data = [], [], [], [], [], [], []

	for key in filtered_data:
		if key.startswith( "bar" ):
			names.append( key )
			barriers.append( filtered_data[ key ] )
		elif key.startswith( "Dist(H-Au)_" ):
			distances_H_Au.append( filtered_data[ key ] )
		elif key.startswith( "Index(H-Au)_" ):
			index_H_Au.append( filtered_data[ key ] )
		elif key.startswith( "Dist(O-cation)_" ):
			distances_O_cation.append( filtered_data[ key ] )
		elif key.startswith( "Dist(O-H)" ):
			distances_O_H.append( filtered_data[ key ] ) 
		elif key.startswith( "Bond_data_" ):
			bond_data.append( filtered_data[ key ] )

	data[ "CONF" ] = names
	data[ "barrier" ] = barriers
	data[ "D(H-Au)" ] = distances_H_Au
	data[ "I(H-Au)" ] = index_H_Au
	data[ "D(O-Na)" ] = distances_O_cation
	if aux_key and "splitting" not in aux_key:
		data[ "D(O-H)" ] = distances_O_H
	data[ "I(O-H)" ] = bond_data
	data[ "status" ] = status
	

	pd.set_option( "display.colheader_justify", "left" )
	data.style.set_properties( **{ "text-align": "center" } )
	#data = data.sort_values( by = "barrier" ).reset_index( drop = True )

	return data

#The function will add a suggestion only if the hydrogen (H) of the nitrogen (N) group in the cation, 
#which is closer to the oxygen (O) of the H2O molecule that will dissociate, is different from the hydrogen (H) of the nitrogen (N) group 
#in the cation listed in the ICONST file. In other words, it will suggest performing a simulation using the hydrogen (H) of the nitrogen (N)
#group in the cation that is closer to the H₂O molecule that will dissociate.
# Adds suggestions based on a comparison between "I(O-H)" and "ICONST".
# sorted_data: A pandas DataFrame containing sorted data with "I(O-H)" and "ICONST".
# Returns: The updated DataFrame with a new "suggest" column containing suggestions or NaN values.
def add_suggestions( sorted_data ):
	suggestions = list()
	for i, j in zip( sorted_data[ "I(O-H)" ], sorted_data[ "ICONST" ] ):
		if isinstance( i, float ) and np.isnan( i ):
			suggestions.append( np.nan )
		elif isinstance( i, str ):
			if i.split( "-" )[ - 1 ] != j.split( "-")[ -1 ] :
				suggestions.append( f"{ i.split( '-' )[ 0 ] } { j.split( '-' )[ 1 ] } { i.split( '-' )[ -1 ] }" )
			else:
				suggestions.append( np.nan )
		else:
			suggestions.append( np.nan )
	sorted_data[ "suggest" ] = suggestions
	return sorted_data

# Cleans the given DataFrame based on conditions in the "CONF" column.
# - Drops specific columns based on whether "CONF" contains certain substrings.
# - Removes columns that contain only NaN values.
# dataframe: A pandas DataFrame containing sorted data with a "CONF" column.
# Returns: The cleaned DataFrame with unnecessary columns removed.
def get_sorted_data_cleaned( dataframe ):
	conf_series = dataframe[ "CONF" ].astype( str )
	drop_columns = list()

	if not conf_series.str.contains( "Na", na = False ).any():
		drop_columns.append( "D(O-Na)" )

	if conf_series.str.contains( "splitting", na = False).any():
		drop_columns.extend( [ "suggest", "I(O-H)" ] )

	if conf_series.str.contains( "hyd_shell", na = False ).any():
		drop_columns.append( "suggest" )

	dataframe = dataframe.drop( columns = [ col for col in drop_columns if col in dataframe.columns ] )
	dataframe = dataframe.dropna( axis = 1, how = "all" )

	return dataframe


# Retrieves barrier data from the database and processes it into a structured DataFrame.
# database: A dictionary containing the database with relevant entries.
# val: The key for the specific dataset to be retrieved from the database.
# fixed_length: The length to which the "CONF" column entries should be padded/truncated. Default is 45.
# verbose: If True, prints the DataFrame. Default is False.
# Returns: A sorted pandas DataFrame containing processed data, or None if no valid data is found.
def get_barrier_from_db( database, val, fixed_length = 43, verbose = False):
	filtered_data = {}
	status = list()
	ICONST_indices = list()

	for key, value in database[ val ].items():
		if value[ "note" ] in [ "Bad" ]: #[ "Good", "Bad" ]:
			status.append( value[ "note" ] )
			updated_data, ICONST_idx = process_database_entry( value, filtered_data )
			filtered_data.update( updated_data )
			ICONST_indices.append( ICONST_idx )
	path_key = value[ "path" ].split( "/" )[ -2 ]

	if not filtered_data:
		print( "No valid data found in database for:", val, "\n" )
		return None

	aux_key = get_aux_key( path_key )
	sorted_data = create_dataframe( filtered_data, status, aux_key )
	sorted_data[ "CONF" ] = sorted_data[ "CONF" ].apply( lambda x: x[ :fixed_length ].ljust( fixed_length ) )
	sorted_data[ "ICONST" ] = ICONST_indices
	
	sorted_data = add_suggestions( sorted_data )
	sorted_data = sorted_data.sort_values( by = "barrier" ).reset_index( drop = True )
	
	cols = [ col for col in sorted_data.columns if col != "status" ] + [ "status" ]
	sorted_data = sorted_data.reindex( columns = cols )
	
	sorted_data = get_sorted_data_cleaned( sorted_data )

	if verbose:
		print( sorted_data.to_string(), "\n")

	return sorted_data

if __name__ == "__main__":
	path = "/home/theodoros/PROJ_ElectroCat/theodoros/HER/Au/HER_Au/database/"
	data = load_database( path + "database_for_theo.js" )	

	#Na_1_hyd = get_barrier_from_db( data, "1_Na_H2O_dissociation_from_hydration_shell", verbose = True )

	#Na_1_No_hyd = get_barrier_from_db( data, "1_Na_H2O_dissociation_NOT_from_hydration_shell", verbose = True )	


	#Na_3_hyd = get_barrier_from_db( data, "3_Na_H2O_dissociation_from_hydration_shell", verbose = True )

	#Na_3_No_hyd = get_barrier_from_db( data, "3_Na_H2O_dissociation_NOT_from_hydration_shell", verbose = True )
	

	#Na_5_hyd = get_barrier_from_db( data, "5_Na_H2O_dissociation_from_hydration_shell", verbose = True )
	
	#Na_5_No_hyd = get_barrier_from_db( data, "5_Na_H2O_dissociation_NOT_from_hydration_shell", verbose = True )
	

	#NH4_1_hyd = get_barrier_from_db( data, "1_NH4_H2O_dissociation_from_hydration_shell", verbose = True )
	
	#NH4_1_NO_hyd = get_barrier_from_db( data, "1_NH4_H2O_dissociation_NOT_from_hydration_shell", verbose = True )
	
	#NH4_1_splitting = get_barrier_from_db( data, "1_NH4_spliting", verbose = True )

	#NH4_1_shuttling = get_barrier_from_db( data, "1_NH4_shuttling", verbose = True )

	
	#NH4_3_hyd = get_barrier_from_db( data, "3_NH4_H2O_dissociation_from_hydration_shell", verbose = True )

	#NH4_3_NO_hyd = get_barrier_from_db( data, "3_NH4_H2O_dissociation_NOT_from_hydration_shell", verbose = True )

	#NH4_3_splitting = get_barrier_from_db( data, "3_NH4_spliting", verbose = True )

	#NH4_3_shuttling = get_barrier_from_db( data, "3_NH4_shuttling", verbose = True )


	#CH3NH3_1_hyd = get_barrier_from_db( data, "1_CH3NH3_H2O_dissociation_from_hydration_shell", verbose = True )

	#CH3NH3_1_NO_hyd = get_barrier_from_db( data, "1_CH3NH3_H2O_dissociation_NOT_from_hydration_shell", verbose = True )

	#CH3NH3_1_splitting = get_barrier_from_db( data, "1_CH3NH3_spliting", verbose = True )

	#CH3NH3_1_shuttling = get_barrier_from_db( data, "1_CH3NH3_shuttling", verbose = True )


	CH3NH3_3_hyd = get_barrier_from_db( data, "3_CH3NH3_H2O_dissociation_from_hydration_shell", verbose = True )

	CH3NH3_3_NO_hyd = get_barrier_from_db( data, "3_CH3NH3_H2O_dissociation_NOT_from_hydration_shell", verbose = True )

	CH3NH3_3_splitting = get_barrier_from_db( data, "3_CH3NH3_spliting", verbose = True )

	CH3NH3_3_shuttling = get_barrier_from_db( data, "3_CH3NH3_shuttling", verbose = True )


	#CH3NH3_5_hyd = get_barrier_from_db( data, "5_CH3NH3_H2O_dissociation_from_hydration_shell", verbose = True )
	
	#CH3NH3_5_NO_hyd = get_barrier_from_db( data, "5_CH3NH3_H2O_dissociation_NOT_from_hydration_shell", verbose = True )
	
	#CH3NH3_5_splitting = get_barrier_from_db( data, "5_CH3NH3_spliting", verbose = True )

	#CH3NH3_5_shuttling = get_barrier_from_db( data, "5_CH3NH3_shuttling", verbose = True )

