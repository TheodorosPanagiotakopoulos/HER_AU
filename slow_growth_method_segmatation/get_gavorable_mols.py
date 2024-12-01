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

def get_H2O_close_to_NH4_with_distances(system, H2O_close_to_electrode, NH4_mols, threshold=2.5):
	system = read( "POSCAR" )
	results = list()
	for h2o in H2O_close_to_electrode:
		h1_idx, o_idx, h2_idx = h2o
		H2O_positions = system.positions[ [ h1_idx, o_idx, h2_idx ] ]
		for nh4 in NH4_mols:
			nh4_positions = system.positions[ nh4 ]
			distances = cdist(H2O_positions, nh4_positions)
			min_distance = np.min( distances )
			if min_distance < threshold:
				results.append( ( h2o, nh4, min_distance ) )
				break
	for i in results:
		print( i )
	return results



if __name__ == "__main__":
	H2O_mols = get_H2O_mols( "POSCAR" )
	NH4_mols = get_NH4_mols( "POSCAR" )
	H2O_close = get_closest_H2O_to_surface( "POSCAR", H2O_mols )
	get_H2O_close_to_NH4_with_distances( "POSCAR", H2O_close, NH4_mols )
