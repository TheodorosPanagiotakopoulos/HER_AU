import pandas as pd
import numpy as np
from ase.io import read
from scipy.spatial.distance import cdist

################################## H2O molecules ##################################

#H2O molecules in the system
def get_H2O_mols( poscar, threshold = 1.2, to_print = "False" ):
	system = read( poscar )	
	oxigen_indices = [ i for i, j in enumerate( system ) if j.symbol == "O" ]
	hydrogen_indices = [ i for i, j in enumerate( system ) if j.symbol == "H" ]
	
	H2O_mols = list()
	for i in oxigen_indices:
		distances = system.get_distances( i, hydrogen_indices, mic = True )
		close_hydrogens = [ hydrogen_indices[ i ] for i, j in enumerate( distances ) if j < threshold ]
		if len( close_hydrogens ) == 2:
			H2O_mols.append( [ close_hydrogens[ 0 ], i, close_hydrogens[ 1 ] ] )
	if to_print == "True":
		print( H2O_mols )
	return ( H2O_mols )

#H2O molecules within the surface threshold
#For H2O dissociation (best)
#H2O -> OH + H*
def get_H2O_within_surface_threshold(poscar, H2O_mols, distance_threshold=2.6):
    system = read(poscar)
    au_indices = [i for i, atom in enumerate(system) if atom.symbol == "Au"]
    au_positions = system.positions[au_indices]

    results = []

    for h2o in H2O_mols:
        h1_idx, o_idx, h2_idx = h2o
        H2O_positions = system.positions[[h1_idx, o_idx, h2_idx]]

        H_positions = [system.positions[h1_idx], system.positions[h2_idx]]
        distances_to_Au = cdist(H_positions, au_positions)
        min_dist_idx = np.argmin(distances_to_Au)
        min_distance = distances_to_Au.flatten()[min_dist_idx]

        if min_distance < distance_threshold:
            closest_H_idx = h1_idx if min_dist_idx // len(au_indices) == 0 else h2_idx
            closest_Au_idx = au_indices[min_dist_idx % len(au_indices)]

            results.append({
                "H2O": f"[{h2o[0]}, {h2o[1]}, {h2o[2]}]",
                "Closest H": closest_H_idx,
                "Closest Au": closest_Au_idx,
                "Distance": round(min_distance, 3)
            })

    # Sort results by Distance
    results = sorted(results, key=lambda x: x["Distance"])

    # Convert to DataFrame for clean display
    df = pd.DataFrame(results)

    # Print the DataFrame
    print(df)

    return [list(map(int, h2o.strip("[]").split(", "))) for h2o in df["H2O"]]



#Closest H2O molecule to the electrode based on the hydrogen atom of the H2O molecule nearest to the electrode
#The function returns only one H2O molecule
#For H2O dissociation
#H2O -> OH + H*
def get_closest_H2O_to_electrode( poscar, H2O_close_to_electrode, to_print = "False" ):
	system = read( poscar )
	au_indices = [ i for i, atom in enumerate(system) if atom.symbol == "Au" ]
	na_indices = [ i for i, atom in enumerate(system) if atom.symbol == "Na" ]

	au_positions = system.positions[ au_indices ]
	na_positions = system.positions[ na_indices ]
	na_bonding_threshold = 2.65

	closest_H2O = None
	min_distance = float( "inf" )

	for h2o in H2O_close_to_electrode:
		h1_idx, o_idx, h2_idx = h2o
		H2O_positions = system.positions[ [ h1_idx, o_idx, h2_idx ] ]

		distances_to_na = cdist( H2O_positions, na_positions )
		if np.any( distances_to_na < na_bonding_threshold ):
			continue 

		distances = cdist( H2O_positions, au_positions )
		min_distance_idx = np.unravel_index( np.argmin( distances ), distances.shape )
		distance = distances[ min_distance_idx ]

		if distance < min_distance:
			min_distance = distance
			h2_idx_closest = [ h1_idx, o_idx, h2_idx ][ min_distance_idx[ 0 ] ]
			au_idx_closest = au_indices[ min_distance_idx[ 1 ] ]
			closest_H2O = ( [ h1_idx, o_idx, h2_idx ], h2_idx_closest, au_idx_closest, round( min_distance, 3 ) )
	if to_print == "True":
		print( closest_H2O )
	return closest_H2O


################################## Na molecules ##################################

#H2O molecules that are close to electrode and belong to Na hydration shell
#For H2O dissociation that belongs to Na hydration shell
#H2O -> OH + H*
def get_H2O_near_electrode_from_Na_hydration_shell( poscar, H2O_close_to_electrode, distance_threshold = 2.6, to_print = "False" ):
	system = read( poscar )
	au_indices = [ i for i, atom in enumerate(system) if atom.symbol == "Au" ]
	na_indices = [ i for i, atom in enumerate(system) if atom.symbol == "Na" ]
	au_positions = system.positions[ au_indices ]

	results = list()

	for h2o in H2O_close_to_electrode:
		h1_idx, o_idx, h2_idx = h2o
		H2O_positions = system.positions[ [ h1_idx, o_idx, h2_idx ] ]

		is_attached_to_Na = any( np.linalg.norm(system.positions[ o_idx ] - system.positions[ na_idx ] ) < distance_threshold for na_idx in na_indices )
		if not is_attached_to_Na:
			continue
		distances_h1 = cdist( [ system.positions[ h1_idx ] ], au_positions ).flatten()
		distances_h2 = cdist( [ system.positions[ h2_idx ] ], au_positions ).flatten()

		min_h1_distance = np.min( distances_h1 )
		min_h2_distance = np.min( distances_h2 )

		if min_h1_distance < min_h2_distance:
			closest_H_idx = h1_idx
			min_distance = min_h1_distance
			closest_Au_idx = au_indices[ np.argmin( distances_h1 ) ]
		else:
			closest_H_idx = h2_idx
			min_distance = min_h2_distance
			closest_Au_idx = au_indices[ np.argmin(distances_h2 ) ]

		if min_distance < distance_threshold:
			results.append( [ [o_idx, h1_idx, h2_idx], f"O: {o_idx}", f"H: {closest_H_idx}", f"Au_closest_to_H_idx: {closest_Au_idx}", f"distance: {round(min_distance, 3)}" ] )
	results.sort(key = lambda x: float( x[ -1 ].split( ': ')[ 1 ] ) )
	if to_print == "True":
		for i in results:
			print( i )
	return results

#H2O molecules that are close to electrode and do NOT belong to Na hydration shell 
#For H2O dissociation that does NOT belong to Na hydration shell
#H2O -> OH + H*
def get_H2O_near_electrode_NOT_from_Na_hydration_shell( poscar, H2O_close_to_electrode, distance_threshold = 3.6, to_print = "False" ):
	system = read( poscar )
	au_indices = [ i for i, atom in enumerate( system ) if atom.symbol == "Au" ]
	na_indices = [ i for i, atom in enumerate( system ) if atom.symbol == "Na" ]
	au_positions = system.positions[ au_indices ]

	results = list()

	for h2o in H2O_close_to_electrode:
		h1_idx, o_idx, h2_idx = h2o
		H2O_positions = system.positions[ [ h1_idx, o_idx, h2_idx ] ]

		is_attached_to_Na = any( np.linalg.norm( system.positions[ o_idx ] - system.positions[ na_idx ] ) < distance_threshold for na_idx in na_indices )
		if is_attached_to_Na:
			continue
		distances_h1 = cdist( [ system.positions[ h1_idx ] ], au_positions ).flatten()
		distances_h2 = cdist( [ system.positions[ h2_idx ] ], au_positions ).flatten()

		min_h1_distance = np.min( distances_h1 )
		min_h2_distance = np.min( distances_h2 )

		if min_h1_distance < min_h2_distance:
			closest_H_idx = h1_idx
			min_distance = min_h1_distance
			closest_Au_idx = au_indices[ np.argmin( distances_h1 ) ]
		else:
			closest_H_idx = h2_idx
			min_distance = min_h2_distance
			closest_Au_idx = au_indices[ np.argmin( distances_h2 ) ]
		if min_distance < distance_threshold:
			results.append( [ [o_idx, h1_idx, h2_idx], f"O: {o_idx}", f"H: {closest_H_idx}", f"Au_closest_to_H_idx: {closest_Au_idx}", f"distance: {round(min_distance, 3)}" ] )
	results.sort( key = lambda x: float( x[ -1 ].split( ': ' )[ 1 ] ) )
	if to_print == "True":
		for i in results:
			print( i )
	return results


################################## NH4 molecules ##################################

#NH4 molecules in the system
def get_NH4_mols( poscar, threshold = 1.2, to_print = "False" ):
	system = read( poscar )
	nitrogen_indices = [i for i, j in enumerate( system ) if j.symbol == "N" ]
	hydrogen_indices = [i for i, j in enumerate( system ) if j.symbol == "H" ]

	NH4_mols = list()
	for i in nitrogen_indices:
		distances = system.get_distances( i , hydrogen_indices, mic = True )
		close_hydrogens = [ hydrogen_indices[ i ] for i, j in enumerate( distances ) if j < threshold ]
		if len( close_hydrogens ) == 4:
			NH4_mols.append( [ i ] + close_hydrogens )
	if to_print == "True":
		print( NH4_mols )
	return NH4_mols

#NH4 molecules within the surface threshold
#For NH4 split
#NH4 -> NH3 + H* (best)
def get_NH4_within_surface_threshold( poscar, NH4_mols, distance_threshold = 5.6, to_print = "False" ):
	system = read( poscar )
	au_indices = [ i for i, atom in enumerate(system) if atom.symbol == "Au" ]
	au_positions = system.positions[ au_indices ]

	NH4_close_to_electrode = list()
	results = list()

	for mol in NH4_mols:
		N_idx, H1_N, H2_N, H3_N, H4_N = mol

		NH4_H_indices = [ H1_N, H2_N, H3_N, H4_N ]
		NH4_H_positions = system.positions[ NH4_H_indices ]

		distances_to_Au = cdist( NH4_H_positions, au_positions )
		min_dist_idx = np.argmin( distances_to_Au )
		min_distance = distances_to_Au.flatten()[ min_dist_idx ]

		if min_distance < distance_threshold:
			closest_H_idx = NH4_H_indices[ min_dist_idx // len( au_indices ) ]
			closest_Au_idx = au_indices[ min_dist_idx % len( au_indices ) ]

			results.append((
				min_distance,
				[ N_idx, H1_N, H2_N, H3_N, H4_N ],
				closest_H_idx,
				closest_Au_idx
			))
			NH4_close_to_electrode.append( mol )

	results.sort( key = lambda x: x[0] )
	if to_print == "True":
		for distance, mol, closest_H_idx, closest_Au_idx in results:
			print( f"[{mol[0]}, {mol[1]}, {mol[2]}, {mol[3]}, {mol[4]}]\t" f"H(N): {closest_H_idx}\t" f"Au: {closest_Au_idx}\t" f"Dist: {round(distance, 3)}" )
	return NH4_close_to_electrode

#Closest NH4 molecule to the electrode based on the hydrogen atom of the NH4 molecule nearest to the electrode
#The function returns only one NH4 molecule
#For NH4 split
#NH4 -> NH3 + H*
def get_NH4_closest_to_electrode( poscar, NH4_mols, to_print = "False" ):
	system = read( poscar )
	au_indices = [ i for i, atom in enumerate( system ) if atom.symbol == "Au" ]

	au_positions = system.positions[ au_indices ]
	closest_NH4 = None
	min_distance = float( "inf" )

	for nh4 in NH4_mols:
		n_idx, h1_idx, h2_idx, h3_idx, h4_idx = nh4
		NH4_positions = system.positions[ [n_idx, h1_idx, h2_idx, h3_idx, h4_idx] ]

		H_positions = NH4_positions[ 1: ]
		distances = cdist( H_positions, au_positions )
		min_distance_idx = np.unravel_index( np.argmin( distances ), distances.shape )
		distance = distances[ min_distance_idx ]

		if distance < min_distance:
			min_distance = distance
			closest_h_idx = [ h1_idx, h2_idx, h3_idx, h4_idx ][ min_distance_idx[ 0 ] ]
			au_idx_closest = au_indices[ min_distance_idx[ 1 ] ]
			closest_NH4 = ( [ n_idx, h1_idx, h2_idx, h3_idx, h4_idx ], closest_h_idx, au_idx_closest, round(min_distance, 3) )
	if to_print == "True":
		print( closest_NH4 )
	return closest_NH4

#H2O molecules close to surface and close to NH4
#For shuttling (Best)
def get_H2O_close_to_surface_and_NH4( poscar, H2O_close_to_electrode, NH4_mols, max_distance_to_Au = 4.5, max_distance_O_to_H_of_NH4 = 3.8, to_print = "False" ):
	system = read( poscar )
	Au_indices = [ i for i, atom in enumerate( system ) if atom.symbol == "Au" ]
	results = list()

	for h2o in H2O_close_to_electrode:
		h1_idx, o_idx, h2_idx = h2o
		O_position = system.positions[ o_idx ]
		H2O_H_positions = [ system.positions[ h1_idx ], system.positions[ h2_idx ] ]

		distances_to_Au = cdist( H2O_H_positions, system.positions[ Au_indices ] )
		min_H2O_Au_dist_idx = np.unravel_index( np.argmin( distances_to_Au ), distances_to_Au.shape )
		min_H2O_Au_dist = distances_to_Au[ min_H2O_Au_dist_idx ]
		H_idx_of_H2O_closest_to_th_electrode = h1_idx if min_H2O_Au_dist_idx[ 0 ] == 0 else h2_idx
		closest_Au_idx = Au_indices[ min_H2O_Au_dist_idx[ 1 ] ]

		for nh4 in NH4_mols:
			N_idx = next( i for i in nh4 if system[ i ].symbol == "N" )
			NH4_H_indices = [ i for i in nh4 if system[ i ].symbol == "H" ]
			distance_O_to_N = np.linalg.norm( O_position - system.positions[ N_idx ] )

			if distance_O_to_N < max_distance_O_to_H_of_NH4:
				NH4_H_positions = system.positions[ NH4_H_indices ]

				distances_to_O = cdist( [ O_position ], NH4_H_positions )
				min_distance_idx = np.unravel_index( np.argmin( distances_to_O ), distances_to_O.shape )
				min_distance = distances_to_O[ min_distance_idx ]
				H_of_NH4_idx_closest_to_O_of_H2O = NH4_H_indices[ min_distance_idx[ 1 ] ]

				if min_distance < max_distance_to_Au:
					results.append({
						"H2O": [o_idx, h1_idx, h2_idx],
						"NH4": [N_idx] + NH4_H_indices,
						"O": o_idx,
						"H2O_closest": {
							"H_index": H_idx_of_H2O_closest_to_th_electrode,
							"distance_to_Au": round(min_H2O_Au_dist, 3),
							"Au_idx": closest_Au_idx
						},
						"NH4_closest": {
							"H_of_NH4_idx_closest_to_O_of_H2O": H_of_NH4_idx_closest_to_O_of_H2O,
							"distance_to_O": round(min_distance, 3)
						}
					})

	results.sort( key = lambda x: x[ "NH4_closest" ][ "distance_to_O" ] )
	if to_print == "True":
		for result in results:
			print( result )
			print( "\n" )
	return results

#H2O molecules close to NH4, up to a threshold
#For shuttling
def get_H2O_close_to_NH4( psocar, H2O_close_to_electrode, NH4_mols, threshold = 2.5, to_print = "False" ):
	system = read( poscar )
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
	if to_print == "True":
		print( results )
	return results


################################## CH3NH3 molecules ##################################

#CH3NH3 molecules in the system
def get_CH3NH3_mols( poscar, threshold = 1.2, to_print = "False" ):
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
	if to_print == "True":
		print( "CH3NH3_mols = ", CH3NH3_mols )
	return CH3NH3_mols

#CH3NH3 molecules within the surface threshold
##For CH3NH3 dissociation (best)
#CH3NH3 -> CH3NH2 + H*
def get_CH3NH3_within_surface_threshold( poscar, CH3NH3_mols, distance_threshold = 3.6, to_print = "False" ):
	system = read( poscar )
	au_indices = [ i for i, atom in enumerate(system) if atom.symbol == "Au" ]
	au_positions = system.positions[ au_indices ]

	CH3NH3_close_to_electrode = list()
	results = list()

	for mol in CH3NH3_mols:
		N_idx, H1_N, H2_N, H3_N, C_idx, H1_C, H2_C, H3_C = mol
		
		NH3_H_indices = [ H1_N, H2_N, H3_N ]
		NH3_H_positions = system.positions[ NH3_H_indices ]

		distances_to_Au = cdist( NH3_H_positions, au_positions )
		min_dist_idx = np.argmin( distances_to_Au )
		min_distance = distances_to_Au.flatten()[ min_dist_idx ]

		if min_distance < distance_threshold:
			closest_H_idx = NH3_H_indices[ min_dist_idx // len( au_indices ) ]
			closest_Au_idx = au_indices[ min_dist_idx % len( au_indices ) ]

			results.append((
				min_distance,
				[N_idx, H1_N, H2_N, H3_N, C_idx, H1_C, H2_C, H3_C],
				closest_H_idx,
				closest_Au_idx
			))
			CH3NH3_close_to_electrode.append( mol )

	results.sort(key = lambda x: x[ 0 ] )
	if to_print == "True":
		for distance, mol, closest_H_idx, closest_Au_idx in results:
			print( f"[{mol[0]}, {mol[1]}, {mol[2]}, {mol[3]}, {mol[4]}, {mol[5]}, {mol[6]}, {mol[7]}] \t" f"H(N): {closest_H_idx} \t" f"Au: {closest_Au_idx} \t" f"Dist: {round( distance, 3 ) }" )
	return CH3NH3_close_to_electrode

#Closest CH3NH3 molecule to the electrode based on the hydrogen atom of the NH4 molecule nearest to the electrode
#The function returns only one NH4 molecule
#For CH3NH3 split
#CH3NH3 -> CH3NH2 + H*
def get_CH3NH3_closest_to_electrode( poscar, CH3NH3_mols, to_print = "False" ):
	system = read( poscar )
	au_indices = [ i for i, atom in enumerate( system ) if atom.symbol == "Au" ]

	au_positions = system.positions[ au_indices ]
	closest_CH3NH3 = None
	min_distance = float( "inf" )

	for ch3nh3 in CH3NH3_mols:
		n_idx, h1_idx, h2_idx, h3_idx, c_idx, h4_idx, h5_idx, h6_idx = ch3nh3
		CH3NH3_positions = system.positions[ [n_idx, h1_idx, h2_idx, h3_idx, c_idx, h4_idx, h5_idx, h6_idx] ]

		H_positions = CH3NH3_positions[ 1: ]
		distances = cdist( H_positions, au_positions )
		min_distance_idx = np.unravel_index( np.argmin( distances ), distances.shape )
		distance = distances[ min_distance_idx ]

		if distance < min_distance:
			min_distance = distance
			closest_h_idx_of_NH3_group = [ h1_idx, h2_idx, h3_idx ][ min_distance_idx[ 0 ] ]
			au_idx_closest = au_indices[ min_distance_idx[ 1 ] ]
			closest_CH3NH3 = ( [ n_idx, h1_idx, h2_idx, h3_idx, c_idx, h4_idx, h5_idx, h6_idx ], closest_h_idx_of_NH3_group, au_idx_closest, round(min_distance, 3 ) )
	if to_print == "True":
		print( closest_CH3NH3 )
	return closest_CH3NH3

#H2O molecules close to surface and close to CH3NH3
#For shuttling (Best)
def get_H2O_close_to_surface_and_CH3NH3( poscar, H2O_close_to_electrode, CH3NH3_mols, max_distance_to_Au = 5.5, max_distance_of_H_of_NH3_group_of_CH3NH3_to_O_of_H2O = 4, to_print = "False" ):
	system = read( poscar )
	Au_indices = [ i for i, atom in enumerate(system) if atom.symbol == "Au" ]
	results = list()

	for h2o in H2O_close_to_electrode:
		h1_idx, o_idx, h2_idx = h2o
		O_position = system.positions[ o_idx ]
		H2O_H_positions = [ system.positions[ h1_idx ], system.positions[ h2_idx ] ]

		distances_to_Au = cdist( H2O_H_positions, system.positions[ Au_indices ] )
		min_H2O_Au_dist_idx = np.unravel_index(np.argmin( distances_to_Au ), distances_to_Au.shape )
		min_H2O_Au_dist = distances_to_Au[ min_H2O_Au_dist_idx ]
		H_idx_of_H2O_closest_to_th_electrode = h1_idx if min_H2O_Au_dist_idx[ 0 ] == 0 else h2_idx
		closest_Au_idx = Au_indices[ min_H2O_Au_dist_idx[ 1 ] ]

		for ch3nh3 in CH3NH3_mols:
			N_idx = next( i for i in ch3nh3 if system[i].symbol == "N" )
			NH3_H_indices = [ i for i in ch3nh3 if system[i].symbol == "H" and np.linalg.norm( system.positions[ i ] - system.positions[ N_idx ] ) < 1.2 ]
			distances_to_O = cdist( [ O_position ], system.positions[ NH3_H_indices ] )
			min_distance_idx = np.unravel_index( np.argmin( distances_to_O ), distances_to_O.shape )
			min_distance = distances_to_O[ min_distance_idx ]
			H_of_NH3_idx_closest_to_O_of_H2O = NH3_H_indices[ min_distance_idx[ 1 ] ]

			if min_distance < max_distance_of_H_of_NH3_group_of_CH3NH3_to_O_of_H2O:
				results.append({
					"H2O": [o_idx, h1_idx, h2_idx],
					"CH3NH3": [N_idx] + NH3_H_indices,
					"O": o_idx,
					"H2O_closest": {
						"H_index": H_idx_of_H2O_closest_to_th_electrode,
						"distance_to_Au": round(min_H2O_Au_dist, 3),
						"Au_idx": closest_Au_idx
					},
					"CH3NH3_closest": {
						"H_of_NH3_idx_closest_to_O_of_H2O": H_of_NH3_idx_closest_to_O_of_H2O,
						"distance_to_O": round(min_distance, 3)
					}
				})

	results.sort( key = lambda x: x[ "CH3NH3_closest" ][ "distance_to_O" ] )
	if to_print == "True":
		for result in results:
			print( result )
			print( "\n" )
	return results

#H2O molecules close to CH3NH3, up to a threshold
#For shuttling
def get_H2O_close_to_CH3NH3( system, H2O_close_to_electrode, CH3NH3_mols, threshold = 2.5, to_print = "False" ):
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
				h2o_atoms = [system.symbols[ idx ] for idx in [ h1_idx, o_idx, h2_idx ] ]
				ch3nh3_atoms = [system.symbols[ idx ] for idx in ch3nh3 ]
				h2o_atom_idx = min_distance_idx[ 0 ]
				ch3nh3_atom_idx = min_distance_idx[ 1 ]
				h2o_symbol = h2o_atoms[ h2o_atom_idx ]
				ch3nh3_symbol = ch3nh3_atoms[ ch3nh3_atom_idx ]
				h2o_idx = [ h1_idx, o_idx, h2_idx ][ h2o_atom_idx ]
				ch3nh3_idx = ch3nh3[ ch3nh3_atom_idx ]
				description = f"symbol { h2o_symbol } = { h2o_idx } - symbol { ch3nh3_symbol } = { ch3nh3_idx }"
				results.append( ( [ h1_idx, o_idx, h2_idx ], ch3nh3, description, round( min_distance, 3 ) ) )
				break
	results.sort( key = lambda x: x [ 3 ] )
	if to_print == "True":
		print( results )
	return results


if __name__ == "__main__":
	CH3NH3_mols = get_CH3NH3_mols( "POSCAR" )
	get_CH3NH3_within_surface_threshold( "POSCAR", CH3NH3_mols, to_print = "True" )
