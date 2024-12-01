import numpy as np
from ase.io import read
from scipy.spatial.distance import cdist

def get_H2O_mols( poscar, threshold = 1.2 ):
	system = read( poscar )	
	oxigen_indices = [ i for i, j in enumerate( system ) if j.symbol == "O" ]
	hydrogen_indices = [ i for i, j in enumerate( system ) if j.symbol == "H" ]
	
	H2O_mols = list()
	for i in oxigen_indices:
		distances = system.get_distances( i, hydrogen_indices, mic = True )
		close_hydrogens = [ hydrogen_indices[ i ] for i, j in enumerate( distances ) if j < threshold ]
		if len( close_hydrogens ) == 2:
			H2O_mols.append( [ close_hydrogens[ 0 ], i, close_hydrogens[ 1 ] ] )
	return ( H2O_mols )

def get_NH4_mols( poscar, threshold = 1.2 ):
	system = read( poscar )
	nitrogen_indices = [i for i, j in enumerate( system ) if j.symbol == "N" ]
	hydrogen_indices = [i for i, j in enumerate( system ) if j.symbol == "H" ]

	NH4_mols = list()
	for i in nitrogen_indices:
		distances = system.get_distances( i , hydrogen_indices, mic = True )
		close_hydrogens = [ hydrogen_indices[ i ] for i, j in enumerate( distances ) if j < threshold ]
		if len( close_hydrogens ) == 4:
			NH4_mols.append( [ i ] + close_hydrogens )
	return NH4_mols


def get_CH3NH3_mols( poscar, threshold = 1.2 ):
	system = read( poscar )
	nitrogen_indices = [ i for i, j in enumerate(system) if j.symbol == "N" ]
	hydrogen_indices = [ i for i, j in enumerate(system) if j.symbol == "H" ]
	carbon_indices = [ i for i, j in enumerate(system) if j.symbol == "C" ]
	CH3NH3_mols = list()
	for i in nitrogen_indices:
		distances = system.get_distances( i, hydrogen_indices, mic=True )
		close_hydrogens = [ hydrogen_indices[ i ] for i, j in enumerate( distances ) if j < threshold ]
		if len( close_hydrogens ) == 3:
			distances_to_carbon = system.get_distances( i, carbon_indices, mic = True )
			close_carbon = [ carbon_indices[ i ] for i, j in enumerate( distances_to_carbon ) if j < 1.55 ]
			if len( close_carbon ) == 1:
				distances_to_hydrogen_from_carbon = system.get_distances(close_carbon[0], hydrogen_indices, mic=True)
				close_hydrogens_from_carbon = [ hydrogen_indices[ i ] for i, j in enumerate( distances_to_hydrogen_from_carbon ) if j < threshold ]
				if len( close_hydrogens_from_carbon ) == 3:
					CH3NH3_mols.append( [ i ] + close_hydrogens + close_carbon + close_hydrogens_from_carbon )
	#print( "CH3NH3_mols = ", CH3NH3_mols )
	return CH3NH3_mols


def get_closest_H2O_to_surface( poscar, H2O_mols, distance_threshold = 2.6 ):
	system = read( "POSCAR" )
	au_indices = [ i for i, j in enumerate( system  ) if j.symbol == "Au" ]
	au_positions = system.positions[ au_indices ]
	H2O_close_to_electrode = list()
	for i in H2O_mols:
		h1_idx, o_idx, h2_idx = i
		H2O_positions = system.positions[ [ h1_idx, o_idx, h2_idx ] ]
		distances = cdist( H2O_positions, au_positions )
		min_distance = np.min( distances )
		if min_distance < distance_threshold:
			H2O_close_to_electrode.append( i )	
	#print( H2O_close_to_electrode )
	return H2O_close_to_electrode

def get_H2O_close_to_NH4( system_path, H2O_close_to_electrode, NH4_mols, threshold=2.5 ):
	system = read( system_path )
	results = []

	for h2o in H2O_close_to_electrode:
		h1_idx, o_idx, h2_idx = h2o
		H2O_positions = system.positions[[h1_idx, o_idx, h2_idx]]

		for nh4 in NH4_mols:
			nh4_positions = system.positions[nh4]
			distances = cdist(H2O_positions, nh4_positions)
			min_distance_idx = np.unravel_index(np.argmin(distances), distances.shape)
			min_distance = distances[min_distance_idx]

			if min_distance < threshold:
				h2o_atoms = [system.symbols[idx] for idx in [h1_idx, o_idx, h2_idx]]
				nh4_atoms = [system.symbols[idx] for idx in nh4]

				h2o_atom_idx = min_distance_idx[0]
				nh4_atom_idx = min_distance_idx[1]
				h2o_symbol = h2o_atoms[h2o_atom_idx]
				nh4_symbol = nh4_atoms[nh4_atom_idx]
				h2o_idx = [h1_idx, o_idx, h2_idx][h2o_atom_idx]
				nh4_idx = nh4[nh4_atom_idx]

				description = f"symbol {h2o_symbol} = {h2o_idx} - symbol {nh4_symbol} = {nh4_idx}"

				results.append( ( [ h1_idx, o_idx, h2_idx ], nh4, description, round( min_distance, 3 ) ) )
				break
	results.sort( key = lambda x: x[ 3 ] )
	for result in results:
		print( result )
	return results


def get_H2O_close_to_CH3NH3( system, H2O_close_to_electrode, CH3NH3_mols, threshold = 2.5 ):
	system = read( system )
	results = list()
	for h2o in H2O_close_to_electrode:
		h1_idx, o_idx, h2_idx = h2o
		H2O_positions = system.positions[[h1_idx, o_idx, h2_idx]]
		for ch3nh3 in CH3NH3_mols:
			ch3nh3_positions = system.positions[ ch3nh3 ]
			distances = cdist( H2O_positions, ch3nh3_positions )
			min_distance_idx = np.unravel_index( np.argmin( distances ), distances.shape )
			min_distance = distances[ min_distance_idx ]
			if min_distance < threshold:
				h2o_atoms = [system.symbols[idx] for idx in [h1_idx, o_idx, h2_idx]]
				ch3nh3_atoms = [system.symbols[idx] for idx in ch3nh3]
				h2o_atom_idx = min_distance_idx[ 0 ]
				ch3nh3_atom_idx = min_distance_idx[ 1 ]
				h2o_symbol = h2o_atoms[ h2o_atom_idx ]
				ch3nh3_symbol = ch3nh3_atoms[ ch3nh3_atom_idx ]
				h2o_idx = [ h1_idx, o_idx, h2_idx ][ h2o_atom_idx ]
				ch3nh3_idx = ch3nh3[ ch3nh3_atom_idx ]
				description = f"symbol {h2o_symbol} = {h2o_idx} - symbol {ch3nh3_symbol} = {ch3nh3_idx}"
				results.append( ( [ h1_idx, o_idx, h2_idx ], ch3nh3, description, round( min_distance, 3 ) ) )
				break
	results.sort( key = lambda x: x [ 3 ] )
	for result in results:
		print( result )
	return results



if __name__ == "__main__":
	H2O_mols = get_H2O_mols( "POSCAR" )
	CH3NH3_mols = get_CH3NH3_mols( "POSCAR" )
	H2O_close = get_closest_H2O_to_surface( "POSCAR", H2O_mols )
	CH3NH3_close_to_H2O = get_H2O_close_to_CH3NH3( "POSCAR", H2O_close, CH3NH3_mols )
	#NH4_mols = get_NH4_mols( "POSCAR" )
	#H2O_close = get_closest_H2O_to_surface( "POSCAR", H2O_mols )
	#get_H2O_close_to_NH4_with_distances( "POSCAR", H2O_close, NH4_mols )
	#print( "---------" )
	#get_H2O_close_to_NH4_with_distances_v2( "POSCAR", H2O_close, NH4_mols )

