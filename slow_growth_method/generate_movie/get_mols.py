import os
import pandas as pd
import numpy as np
from ase.io import read
from scipy.spatial.distance import cdist
import glob

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
def get_H2O_within_surface_threshold( poscar, H2O_mols, distance_threshold = 2.6 ):
	system = read( poscar )
	au_indices = [ i for i, atom in enumerate( system ) if atom.symbol == "Au" ]
	au_positions = system.positions[au_indices]

	results = list()

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
			results.append({ "H2O": f"[{h2o[0]}, {h2o[1]}, {h2o[2]}]", "Closest H": closest_H_idx, "Closest Au": closest_Au_idx, "Distance": round(min_distance, 3 ) } )
	results = sorted( results, key=lambda x: x[ "Distance" ] )
	df = pd.DataFrame( results )
	print( df )

	return [ list( map( int, h2o.strip( "[]" ).split(", ") ) ) for h2o in df[ "H2O" ] ]


################################## Na molecules ##################################

#Na molecules in the system
def get_Na_mols( poscar, to_print = "False" ):
	system = read( poscar )
	na_indices = [ i for i, atom in enumerate(system) if atom.symbol == "Na" ]
	if to_print == "True":
		print( na_indices )
	return na_indices

#Get H2O mols in the Na hydration shell 
#H2O -> OH + H*
def get_Na_hydration_shell(poscar, H2O_mols, Na_atoms, distance_threshold=2.6, to_print="False"):
	system = read(poscar)
	au_indices = [i for i, atom in enumerate(system) if atom.symbol == "Au"]
	au_positions = system.positions[au_indices]

	final_results = list()

	for na_idx in Na_atoms:
		molecule_results = list()

		for h2o in H2O_mols:
			h1_idx, o_idx, h2_idx = h2o

			distance_to_oxygen = np.linalg.norm(system.positions[o_idx] - system.positions[na_idx])
			if distance_to_oxygen > distance_threshold:
				continue

			distances_h1_to_au = cdist([system.positions[h1_idx]], au_positions).flatten()
			distances_h2_to_au = cdist([system.positions[h2_idx]], au_positions).flatten()

			min_h1_distance_to_au = np.min(distances_h1_to_au)
			min_h2_distance_to_au = np.min(distances_h2_to_au)

			closest_au_h1_idx = au_indices[np.argmin(distances_h1_to_au)]
			closest_au_h2_idx = au_indices[np.argmin(distances_h2_to_au)]

			molecule_results.append({ "[Na]": [na_idx], "[H1, O, H2]": [h1_idx, o_idx, h2_idx], "[Au]": [closest_au_h1_idx, closest_au_h2_idx], "[H1 - Au]": [f"{h1_idx} - {closest_au_h1_idx} = {round(min_h1_distance_to_au, 3)}"], "[H2 - Au]": [f"{h2_idx} - {closest_au_h2_idx} = {round(min_h2_distance_to_au, 3)}"], "[Na - O H2O]": round(distance_to_oxygen, 3), } )

		molecule_results.sort(key=lambda x: float(x["[H1 - Au]"][0].split(" = ")[1]))
		final_results.extend( molecule_results )

		if to_print == "True":
			for result in molecule_results:
				print(result)
			print("\n")

	return final_results

#Get H2O mols NOT in the Na hydration shell 
#H2O -> OH + H*
def get_non_Na_hydration_shell( poscar, H2O_mols, Na_atoms, distance_threshold = 2.7, to_print = "False" ):
    system = read( poscar )
    au_indices = [ i for i, atom in enumerate( system ) if atom.symbol == "Au" ]
    au_positions = system.positions[ au_indices ]
    non_hydration_H2O = list()

    for h2o in H2O_mols:
        h1_idx, o_idx, h2_idx = h2o
        is_in_hydration_shell = any(
            np.linalg.norm( system.positions[ o_idx ] - system.positions[ na_idx ] ) <= distance_threshold
            for na_idx in Na_atoms
        )
        if not is_in_hydration_shell:
            distances_h1_to_au = cdist( [ system.positions[ h1_idx ] ], au_positions ).flatten()
            distances_h2_to_au = cdist( [ system.positions[ h2_idx ] ], au_positions ).flatten()
            closest_au_h1_idx = au_indices[ np.argmin( distances_h1_to_au ) ]
            closest_au_h2_idx = au_indices[ np.argmin( distances_h2_to_au ) ]
            min_h1_distance_to_au = np.min( distances_h1_to_au )
            min_h2_distance_to_au = np.min( distances_h2_to_au )
            non_hydration_H2O.append({
                "H2O": [ h1_idx, o_idx, h2_idx ],
                f"{h1_idx}-{closest_au_h1_idx}": round( min_h1_distance_to_au, 3 ),
                f"{h2_idx}-{closest_au_h2_idx}": round( min_h2_distance_to_au, 3 )
            })
    sorted_list = sorted( non_hydration_H2O, key=lambda d: list(d.values() )[ 1 ] )
    if to_print == "True":
        print( "H2O molecules not in Na hydration shell:" )
        for h2o in sorted_list:
            print( h2o )
    return non_hydration_H2O


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
#NH4 -> NH3 + H*
def get_NH4_within_surface_threshold(poscar, NH4_mols, distance_threshold = 5.6, to_print = "False" ):
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
		closest_info = list()
		for h_idx, h_distances in enumerate( distances_to_Au ):
			min_dist_idx = np.argmin( h_distances)
			min_distance = h_distances[ min_dist_idx ]
			closest_Au_idx = au_indices[ min_dist_idx ]
			closest_info.append( ( NH4_H_indices[ h_idx ], closest_Au_idx, min_distance ) )
		closest_info.sort( key=lambda x: x[ 2 ] )
		if any( info[ 2 ] < distance_threshold for info in closest_info ):
			results.append( (
				[N_idx, H1_N, H2_N, H3_N, H4_N],
				closest_info
			) )
			NH4_close_to_electrode.append( mol )
	results.sort( key=lambda x: min( info[ 2 ] for info in x[ 1 ] ) )
	if to_print == "True":
		for mol, closest_info in results:
			print( f"Molecule: {mol}" )
			for h_idx, au_idx, dist in closest_info:
				print( f"\tH(N): {h_idx}\tAu: {au_idx}\tDist: {round(dist, 3)}" )
	return NH4_close_to_electrode


#H2O molecules that belong to NH4 hydration shell
#For H2O dissociation ( H2O -> OH + H* ) if H of H2O is close to electrode
def get_NH4_hydration_shell( poscar, H2O_mols, NH4_molecules, distance_threshold = 3.1, to_print = "False" ):
	system = read( poscar )
	au_indices = [ i for i, atom in enumerate( system ) if atom.symbol == "Au" ]
	au_positions = system.positions[ au_indices ]

	final_results = list()

	for nh4 in NH4_molecules:
		n_idx, h1_nh4_idx, h2_nh4_idx, h3_nh4_idx, h4_nh4_idx = nh4

		molecule_results = list()

		for h2o in H2O_mols:
			h1_idx, o_idx, h2_idx = h2o

			distance_n_to_o = np.linalg.norm( system.positions[ n_idx ] - system.positions[ o_idx ] )

			if distance_n_to_o > distance_threshold:
				continue

			distances_h1_to_au = cdist( [ system.positions[ h1_idx ] ], au_positions ).flatten()
			distances_h2_to_au = cdist( [ system.positions[ h2_idx ] ], au_positions ).flatten()

			min_h1_distance_to_au = np.min( distances_h1_to_au )
			min_h2_distance_to_au = np.min( distances_h2_to_au )

			closest_au1_idx = au_indices[ np.argmin(distances_h1_to_au ) ]
			closest_au2_idx = au_indices[ np.argmin(distances_h2_to_au ) ]

			n_attached_h_distances = {}
			for h_idx in [h1_nh4_idx, h2_nh4_idx, h3_nh4_idx, h4_nh4_idx]:
				distances_to_au = cdist( [ system.positions[ h_idx ] ], au_positions ).flatten()
				min_distance_to_au = np.min( distances_to_au )
				closest_au_idx = au_indices[ np.argmin(distances_to_au ) ]
				n_attached_h_distances[ f"H-Au" ] = f"{h_idx} - {closest_au_idx} = {round(min_distance_to_au, 3)}"

			molecule_results.append({
				"[N, H1, H2, H3, H4]": [ n_idx, h1_nh4_idx, h2_nh4_idx, h3_nh4_idx, h4_nh4_idx ],
				"[H1, O, H2]": [ h1_idx, o_idx, h2_idx ],
				"H1-Au1": f"{h1_idx} - {closest_au1_idx} = {round(min_h1_distance_to_au, 3)}",
				"H2-Au2": f"{h2_idx} - {closest_au2_idx} = {round(min_h2_distance_to_au, 3)}",
				"N-O Distance": round( distance_n_to_o, 3 )
			})

		molecule_results.sort( key=lambda x: x[ "N-O Distance" ] )
		final_results.extend( molecule_results )

		if to_print == "True":
			for result in molecule_results:
				print( result )
			print( "\n" )

	return final_results

#H2O molecules that not belong to NH4 hydration shel
#H2O -> OH + H*
def get_non_NH4_hydration_shell(poscar, H2O_mols, NH4_molecules, distance_threshold=3.2, to_print="False"):
	system = read(poscar)
	au_indices = [i for i, atom in enumerate(system) if atom.symbol == "Au"]
	au_positions = system.positions[au_indices]

	non_hydration_H2O = list()

	for h2o in H2O_mols:
		h1_idx, o_idx, h2_idx = h2o
		is_in_hydration_shell = False

		for nh4 in NH4_molecules:
			n_idx, h1_nh4, h2_nh4, h3_nh4, h4_nh4 = nh4

			distance_n_to_o = np.linalg.norm(system.positions[n_idx] - system.positions[o_idx])

			if distance_n_to_o <= distance_threshold:
				is_in_hydration_shell = True
				break

		if not is_in_hydration_shell:
			distances_h1_to_au = cdist([system.positions[h1_idx]], au_positions).flatten()
			distances_h2_to_au = cdist([system.positions[h2_idx]], au_positions).flatten()

			closest_au_h1_idx = au_indices[np.argmin(distances_h1_to_au)]
			closest_au_h2_idx = au_indices[np.argmin(distances_h2_to_au)]

			min_h1_distance_to_au = np.min(distances_h1_to_au)
			min_h2_distance_to_au = np.min(distances_h2_to_au)

			non_hydration_H2O.append({
				"H2O": [h1_idx, o_idx, h2_idx],
				f"{h1_idx}-{closest_au_h1_idx}": round(min_h1_distance_to_au, 3),
				f"{h2_idx}-{closest_au_h2_idx}": round(min_h2_distance_to_au, 3)
			})

	sorted_list = sorted(non_hydration_H2O, key=lambda d: list(d.values())[1])

	if to_print == "True":
		print("H2O molecules not in NH4 hydration shell:")
		for h2o in sorted_list:
			print(h2o)

	return non_hydration_H2O

#H2O molecules that belong to NH4 hydration shell
#For shuttling
def get_NH4_hydration_shell_shuttling( poscar, H2O_mols, NH4_molecules, distance_threshold = 2.6, to_print="False" ):
	system = read( poscar )
	au_indices = [ i for i, atom in enumerate( system ) if atom.symbol == "Au" ]
	au_positions = system.positions[ au_indices ]

	final_results = list()

	for nh4 in NH4_molecules:
		n_idx, h1_nh4, h2_nh4, h3_nh4, h4_nh4 = nh4
		NH4_H_indices = [ h1_nh4, h2_nh4, h3_nh4, h4_nh4 ]

		molecule_results = list()

		for h2o in H2O_mols:
			h1_idx, o_idx, h2_idx = h2o

			distance_n_to_o = np.linalg.norm( system.positions[ n_idx ] - system.positions[ o_idx ] )

			# Check if the N-O distance is within 4 A
			if distance_n_to_o > 4:
				continue

			closest_nh4_h_to_h2o = None
			min_distance_to_nh4_h = float( 'inf' )
			all_distances_to_o = {}

			for h_nh4_idx in NH4_H_indices:
				distance = np.linalg.norm( system.positions[ o_idx ] - system.positions[ h_nh4_idx ] )
				all_distances_to_o[ h_nh4_idx ] = f"{h_nh4_idx}\t{round(distance, 3)}"
				if distance < distance_threshold and distance < min_distance_to_nh4_h:
					closest_nh4_h_to_h2o = h_nh4_idx
					min_distance_to_nh4_h = distance

			if closest_nh4_h_to_h2o is None and distance_n_to_o > 4:
				continue

			distances_h1_to_au = cdist( [ system.positions[ h1_idx ] ], au_positions ).flatten()
			distances_h2_to_au = cdist( [ system.positions[ h2_idx ] ], au_positions ).flatten()

			min_h1_distance_to_au = np.min( distances_h1_to_au )
			min_h2_distance_to_au = np.min( distances_h2_to_au )

			closest_au1_idx = au_indices[ np.argmin( distances_h1_to_au ) ]
			closest_au2_idx = au_indices[ np.argmin( distances_h2_to_au ) ]

			distances_of_nh4_h_to_o = {}
			for h_idx in NH4_H_indices:
				if h_idx in all_distances_to_o:
					split_distance = all_distances_to_o[ h_idx ].split( '\t' )[ 1 ]
					distances_of_nh4_h_to_o[ h_idx ] = f"H:{h_idx}-O:{o_idx}= {split_distance}"

			molecule_results.append({
				"[N, H1, H2, H3, H4]": [ n_idx, h1_nh4, h2_nh4, h3_nh4, h4_nh4],
				"[H1, O, H2]": [h1_idx, o_idx, h2_idx ],
				"H1-Au1": f"{h1_idx} - {closest_au1_idx} = {round(min_h1_distance_to_au, 3)}",
				"H2-Au2": f"{h2_idx} - {closest_au2_idx} = {round(min_h2_distance_to_au, 3)}",
				"N-O Distance": round( distance_n_to_o, 3 ),
				"Distances of NH4-H to O": distances_of_nh4_h_to_o
			})

		molecule_results.sort( key=lambda x: float(x[ "H1-Au1" ].split( '= ' )[ 1 ] ) )
		final_results.extend( molecule_results )

		if to_print == "True":
			for result in molecule_results:
				print( result )
			print( "\n" )

	return final_results


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
def get_CH3NH3_within_surface_threshold(poscar, CH3NH3_mols, to_print="False"):
	system = read( poscar )
	au_indices = [ i for i, atom in enumerate( system) if atom.symbol == "Au"]
	au_positions = system.positions[ au_indices ]

	results = list()

	for mol in CH3NH3_mols:
		N_idx, H1_N, H2_N, H3_N, C_idx, H1_C, H2_C, H3_C = mol
		NH3_H_indices = [ H1_N, H2_N, H3_N ]
		NH3_H_positions = system.positions[ NH3_H_indices ]

		distances_to_Au = cdist( NH3_H_positions, au_positions )
		h_results = list()

		for idx, H_idx in enumerate( NH3_H_indices ):
			h_distances_to_Au = distances_to_Au[ idx ]
			closest_au_idx = np.argmin( h_distances_to_Au )
			closest_distance = h_distances_to_Au[ closest_au_idx ]

			h_results.append( ( closest_distance, H_idx, au_indices[ closest_au_idx ] ) )

		h_results.sort( key=lambda x: x[0] )  # Sort H results by distance
		results.append( (mol, h_results) )

	if to_print == "True":
		for mol, h_results in results:
			for dist, H_idx, Au_idx in h_results:
				print(
					f"[{mol[0]}, {mol[1]}, {mol[2]}, {mol[3]}, {mol[4]}, {mol[5]}, {mol[6]}, {mol[7]}] \t"
					f"H(N): {H_idx} \t"
					f"Au: {Au_idx} \t"
					f"Dist: {round(dist, 3)}"
				)
			print() 
	return results


#Get H2O mols in the CH3NH3 hydration shell 
#For shuttling
#For H2O -> OH + H* if H of H2O is close to electrode
def get_CH3NH3_hydration_shell(poscar, H2O_mols, CH3NH3_molecules, distance_threshold = 3.2, to_print="False"):
	system = read( poscar )
	au_indices = [ i for i, atom in enumerate( system ) if atom.symbol == "Au" ]
	au_positions = system.positions[ au_indices ]

	final_results = list()

	for ch3nh3 in CH3NH3_molecules:
		n_idx, h1_nh3_idx, h2_nh3_idx, h3_nh3_idx, c_idx, h4_nh3_idx, h5_nh3_idx, h6_nh3_idx = ch3nh3

		molecule_results = list()

		for h2o in H2O_mols:
			h1_idx, o_idx, h2_idx = h2o

			distance_n_to_o = np.linalg.norm( system.positions[ n_idx ] - system.positions[ o_idx ] )

			if distance_n_to_o > distance_threshold:
				continue

			distances_h1_to_au = cdist( [ system.positions[ h1_idx ] ], au_positions ).flatten()
			distances_h2_to_au = cdist( [ system.positions[ h2_idx ] ], au_positions ).flatten()

			min_h1_distance_to_au = np.min( distances_h1_to_au )
			min_h2_distance_to_au = np.min( distances_h2_to_au )

			closest_au1_idx = au_indices[ np.argmin(distances_h1_to_au ) ]
			closest_au2_idx = au_indices[ np.argmin(distances_h2_to_au ) ]

			n_attached_h_distances = {}
			for h_idx in [h4_nh3_idx, h5_nh3_idx, h6_nh3_idx]:
				distances_to_au = cdist( [ system.positions[ h_idx ] ], au_positions ).flatten()
				min_distance_to_au = np.min( distances_to_au )
				closest_au_idx = au_indices[ np.argmin(distances_to_au ) ]
				n_attached_h_distances[f"H-Au"] = f"{h_idx} - {closest_au_idx} = {round(min_distance_to_au, 3)}"

			molecule_results.append({
				"[C, H1, H2, H3, N, H4, H5, H6]": [c_idx, h1_nh3_idx, h2_nh3_idx, h3_nh3_idx, n_idx, h4_nh3_idx, h5_nh3_idx, h6_nh3_idx],
				"[H1, O, H2]": [h1_idx, o_idx, h2_idx],
				"H1-Au1": f"{h1_idx} - {closest_au1_idx} = {round(min_h1_distance_to_au, 3)}",
				"H2-Au2": f"{h2_idx} - {closest_au2_idx} = {round(min_h2_distance_to_au, 3)}",
				"N-O Distance": round(distance_n_to_o, 3)
			})

		molecule_results.sort( key=lambda x: x[ "N-O Distance" ] )
		final_results.extend( molecule_results )

		if to_print == "True":
			for result in molecule_results:
				print(result)
			print("\n")

	return final_results

#H2O molecules that not belong to CH3NH3 hydration shel
#H2O -> OH + H*
def get_non_CH3NH3_hydration_shell(poscar, H2O_mols, CH3NH3_molecules, distance_threshold = 3.2, to_print="False"):
	system = read( poscar )
	au_indices = [ i for i, atom in enumerate( system ) if atom.symbol == "Au" ]
	au_positions = system.positions[ au_indices ]

	non_hydration_H2O = list()

	for h2o in H2O_mols:
		h1_idx, o_idx, h2_idx = h2o
		is_in_hydration_shell = False

		for ch3nh3 in CH3NH3_molecules:
			n_idx, h1_ch3nh3, h2_ch3nh3, h3_ch3nh3, c_idx, h4_ch3nh3, h5_ch3nh3, h6_ch3nh3 = ch3nh3

			distance_n_to_o = np.linalg.norm( system.positions[ n_idx ] - system.positions[ o_idx ] )

			if distance_n_to_o <= distance_threshold:
				is_in_hydration_shell = True
				break

		if not is_in_hydration_shell:
			distances_h1_to_au = cdist( [ system.positions[ h1_idx ] ], au_positions ).flatten()
			distances_h2_to_au = cdist( [ system.positions[ h2_idx ] ], au_positions ).flatten()

			closest_au_h1_idx = au_indices[ np.argmin( distances_h1_to_au ) ]
			closest_au_h2_idx = au_indices[ np.argmin( distances_h2_to_au ) ]

			min_h1_distance_to_au = np.min( distances_h1_to_au )
			min_h2_distance_to_au = np.min( distances_h2_to_au )

			non_hydration_H2O.append({
				"H2O": [h1_idx, o_idx, h2_idx],
				f"{h1_idx}-{closest_au_h1_idx}": round(min_h1_distance_to_au, 3),
				f"{h2_idx}-{closest_au_h2_idx}": round(min_h2_distance_to_au, 3)
			})

	sorted_list = sorted( non_hydration_H2O, key=lambda d: list( d.values() )[ 1 ] )
	
	if to_print == "True":
		print( "H2O molecules not in CH3NH3 hydration shell:" )
		for h2o in sorted_list:
			print( h2o )

	return non_hydration_H2O


#H2O molecules that belong to NH4 hydration shell
#For shuttling
def get_CH3NH3_hydration_shell_shuttling(poscar, H2O_mols, CH3NH3_molecules, distance_threshold = 2.6, to_print="False"):
	system = read( poscar )
	au_indices = [ i for i, atom in enumerate( system ) if atom.symbol == "Au" ]
	au_positions = system.positions[ au_indices ]

	final_results = list()

	for ch3nh3 in CH3NH3_molecules:
		n_idx, h1_ch3nh3, h2_ch3nh3, h3_ch3nh3, c_idx, h4_ch3nh3, h5_ch3nh3, h6_ch3nh3 = ch3nh3
		NH3_H_indices = [ h1_ch3nh3, h2_ch3nh3, h3_ch3nh3 ]
		CH3_H_indices = [ h4_ch3nh3, h5_ch3nh3, h6_ch3nh3 ]

		molecule_results = list()

		for h2o in H2O_mols:
			h1_idx, o_idx, h2_idx = h2o

			distance_n_to_o = np.linalg.norm( system.positions[ n_idx ] - system.positions[ o_idx ] )

			# Check if the N-O distance is within 4 A
			if distance_n_to_o > 4:
				continue

			closest_ch3nh3_h_to_h2o = None
			min_distance_to_ch3nh3_h = float( 'inf' )
			all_distances_to_o = {}

			for h_ch3nh3_idx in NH3_H_indices:
				distance = np.linalg.norm( system.positions[ o_idx ] - system.positions[ h_ch3nh3_idx ] )
				all_distances_to_o[h_ch3nh3_idx] = f"{h_ch3nh3_idx}\t{round(distance, 3)}"
				if distance < distance_threshold and distance < min_distance_to_ch3nh3_h:
					closest_ch3nh3_h_to_h2o = h_ch3nh3_idx
					min_distance_to_ch3nh3_h = distance

			if closest_ch3nh3_h_to_h2o is None and distance_n_to_o > 4:
				continue

			distances_h1_to_au = cdist( [ system.positions[ h1_idx ] ], au_positions ).flatten()
			distances_h2_to_au = cdist( [ system.positions[ h2_idx] ], au_positions ).flatten()

			min_h1_distance_to_au = np.min( distances_h1_to_au )
			min_h2_distance_to_au = np.min( distances_h2_to_au )

			closest_au1_idx = au_indices[ np.argmin(distances_h1_to_au ) ]
			closest_au2_idx = au_indices[ np.argmin(distances_h2_to_au ) ]

			distances_of_ch3nh3_h_to_o = {}
			for h_idx in NH3_H_indices:
				if h_idx in all_distances_to_o:
					split_distance = all_distances_to_o[ h_idx ].split('\t')[ 1 ]
					distances_of_ch3nh3_h_to_o[ h_idx ] = f"H:{h_idx}-O:{o_idx}= {split_distance}"

			molecule_results.append({
				"[N, H1, H2, H3, C, H4, H5, H6]": [n_idx, h1_ch3nh3, h2_ch3nh3, h3_ch3nh3, c_idx, h4_ch3nh3, h5_ch3nh3, h6_ch3nh3],
				"[H1, O, H2]": [h1_idx, o_idx, h2_idx],
				"H1-Au1": f"{h1_idx} - {closest_au1_idx} = {round(min_h1_distance_to_au, 3)}",
				"H2-Au2": f"{h2_idx} - {closest_au2_idx} = {round(min_h2_distance_to_au, 3)}",
				"N-O Distance": round(distance_n_to_o, 3),
				"Distances of CH3NH3-H to O": distances_of_ch3nh3_h_to_o
			})

		molecule_results.sort( key=lambda x: float( x[ "H1-Au1" ].split( '= ' )[ 1 ] ) )
		final_results.extend( molecule_results )

		if to_print == "True":
			for result in molecule_results:
				print( result )
			print("\n")

	return final_results

################################## Get status ##################################
# Retrieves a list of all directories or files in the current directory that match the pattern "RUN*".
# Returns: A list of matching items or an empty list if no matches are found.
def get_RUNs( path_to_simulation ):
	path_to_simulation = os.path.normpath(path_to_simulation )
	runs = glob.glob( os.path.join( path_to_simulation, "RUN*" ) )
	return runs

#Finds the index of the Hydrogen (H) atom in an ICONST file that is moving toward a specified surface.
def get_H_from_ICONST( iconst, to_print = False ):
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
		H_idx = first_line[ 2 ]
		if to_print: 
			print( "ICONST H index: ", H_idx )
			print( "Actual H index: ", int( H_idx ) - 1 ) 
		return int( H_idx ) - 1
	elif len( lines ) == 3:
		with open( iconst ) as file:
			lines = [ line.rstrip() for line in file ]
		first_line = lines[ 0 ].split( " " )
		second_line = lines[ 1 ].split( " " )
		O_H2O_idx = first_line[ 1 ]
		H_H2O_idx = first_line[ 2 ]
		H_cation_idx = second_line[ 2 ]
		if to_print:
			print( "ICONST O(H2O) index: ", O_H2O_idx )
			print( "ICONST H(H2O) index: ", H_H2O_idx )
			print( "ICONST H(cation) index: ", H_cation_idx )
			print( "Actual O(H2O) index: ", int( O_H2O_idx ) -1 )
			print( "Actual H(H2O) index: ", int( H_H2O_idx ) -1 )
			print( "Actual H(cation) index: ", int( H_cation_idx ) -1 )
		return int( O_H2O_idx ) -1, int( H_H2O_idx ) -1, int( H_cation_idx ) -1
	else:
		raise ValueError( "ICONST file must NOT contain more than three lines." )

#check if the H atom from ICONST indeed moved to Au
def get_status( path_to_simulation, threshold_distance = 2.2 ):
	result = False
	if not os.path.isfile( path_to_simulation + "/CONTCAR" ):
		return result
	system = read( path_to_simulation + "/CONTCAR" )
	runs = get_RUNs( path_to_simulation )
	if not runs:
		initial_system = read( path_to_simulation + "/POSCAR")
	else:
		initial_system = read( path_to_simulation + "/RUN1/POSCAR")
	au_indices = [ i for i, j in enumerate( system ) if j.symbol == "Au" ]
	H_info = get_H_from_ICONST( path_to_simulation + "/ICONST" )
	if isinstance( H_info, int ):
		H_idx = H_info
		distance_H_to_Au = [ np.linalg.norm(system.positions[ au_idx ] - system.positions[ H_idx ] ) for au_idx in au_indices ]
		min_H_Au_dist = round( min( distance_H_to_Au ), 3 )
		if min_H_Au_dist < threshold_distance:
			result = True
	elif isinstance( H_info, tuple ):
		O_idx, H_H2O_idx, H_cation_idx = H_info
		distance_H_to_Au = [ np.linalg.norm( system.positions[ au_idx ] - system.positions[ H_H2O_idx ] ) for au_idx in au_indices ]
		min_H_Au_dist = round( min( distance_H_to_Au ), 3 )
		distance_O_to_H_cation = np.linalg.norm( system.positions[ O_idx ] - system.positions[ H_cation_idx ] )
		initial_distance_O_to_H_cation = np.linalg.norm( initial_system.positions[ O_idx ] - initial_system.positions[ H_cation_idx ] )
		#print( "initial dist: ", initial_distance_O_to_H_cation )
		#print( "min H Au dist: ", min_H_Au_dist )
		if min_H_Au_dist < threshold_distance and distance_O_to_H_cation < threshold_distance and initial_distance_O_to_H_cation < threshold_distance:
			result = True
	return result


if __name__ == "__main__":
	loc = os.getcwd()
	status = get_status( loc )
	print( status )
