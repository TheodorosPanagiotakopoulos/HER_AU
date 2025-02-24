import os
import numpy as np

def get_dirs():
	path = os.getcwd()
	return sorted( [name for name in os.listdir(path) if os.path.isdir(os.path.join(path, name)) and name.startswith("RUN")], key=lambda x: int(x[3:]) if x[3:].isdigit() else float('inf') )


def get_data_ICONST( path_to_SG_simulation, verbose = False ):
	iconst_path = os.path.join( path_to_SG_simulation, "RUN1/ICONST" )
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
def get_standarized_ICONST_data( path_to_SG_simulation ):
	iconst_data = get_data_ICONST( path_to_SG_simulation )
	print( path_to_SG_simulation )
	cations = [ "NH4", "CH3NH3" ]
	if len( iconst_data ) == 2 and any( cation + "_splitting" not in path_to_SG_simulation for cation in cations ):
		H_idx, O_idx = iconst_data
		return H_idx, O_idx, np.nan
	elif len( iconst_data ) == 2 and any(cation + "_splitting" in path_to_SG_simulation for cation in cations):
		N_idx, H_idx = iconst_data
		return H_idx, N_idx, np.nan
	elif len( iconst_data ) == 3:
		O_idx, H_idx, H_cation_idx = iconst_data
		return H_idx, O_idx, H_cation_idx
	else:
		raise ValueError( "ICONST file must have 2 or 3 lines" )

get_standarized_ICONST_data( os.getcwd() )
