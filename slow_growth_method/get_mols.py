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
def get_Na_hydration_shell( poscar, H2O_mols, Na_atoms, distance_threshold = 2.6, to_print = "False" ):
	system = read( poscar )
	au_indices = [ i for i, atom in enumerate( system ) if atom.symbol == "Au" ]
	au_positions = system.positions[ au_indices ]

	final_results = list()

	for na_idx in Na_atoms:  
		molecule_results = list()

		for h2o in H2O_mols:
			h1_idx, o_idx, h2_idx = h2o

			distance_to_oxygen = np.linalg.norm(system.positions[o_idx] - system.positions[na_idx])
			if distance_to_oxygen > distance_threshold:
				continue

			distances_h1_to_au = cdist( [ system.positions[ h1_idx ] ], au_positions ).flatten()
			distances_h2_to_au = cdist( [ system.positions[ h2_idx ] ], au_positions ).flatten()

			min_h1_distance_to_au = np.min( distances_h1_to_au )
			min_h2_distance_to_au = np.min( distances_h2_to_au )

			if min_h1_distance_to_au < min_h2_distance_to_au:
				closest_h2o_h_idx = h1_idx
				min_distance_to_au = min_h1_distance_to_au
				closest_au_idx = au_indices[ np.argmin( distances_h1_to_au ) ]
			else:
				closest_h2o_h_idx = h2_idx
				min_distance_to_au = min_h2_distance_to_au
				closest_au_idx = au_indices[ np.argmin( distances_h2_to_au ) ]

			molecule_results.append({
				"[Na]": [na_idx],
				"[H1, O, H2]": [h1_idx, o_idx, h2_idx],
				"[Au]": [closest_au_idx],
				"[H - Au]": [f"{closest_h2o_h_idx} - {closest_au_idx} = {round(min_distance_to_au, 3)}"],
				"[Na - O H2O]": round(distance_to_oxygen, 3),
				"O-H": [o_idx, closest_h2o_h_idx ]
			})

		molecule_results.sort( key=lambda x: float( x["[H - Au]" ][ 0 ].split( '= ' )[ 1 ] ) )
		final_results.extend( molecule_results )

		if to_print == "True":
			for result in molecule_results:
				print( result )
			print( "\n" )

	return final_results

#Get H2O mols NOT in the Na hydration shell 
#H2O -> OH + H*
def get_non_Na_hydration_shell(poscar, H2O_mols, Na_atoms, distance_threshold=2.6, to_print="False"):
	system = read( poscar )
	au_indices = [ i for i, atom in enumerate(system) if atom.symbol == "Au" ]
	au_positions = system.positions[ au_indices ]
	non_hydration_H2O = list()

	for h2o in H2O_mols:
		h1_idx, o_idx, h2_idx = h2o

		is_in_hydration_shell = False
		for na_idx in Na_atoms:
			distance_to_oxygen = np.linalg.norm(system.positions[o_idx] - system.positions[na_idx])
			if distance_to_oxygen <= distance_threshold:
				is_in_hydration_shell = True
				break  

		if not is_in_hydration_shell:
			distances_h1_to_au = cdist( [ system.positions[ h1_idx ] ], au_positions ).flatten()
			distances_h2_to_au = cdist( [ system.positions[ h2_idx ] ], au_positions ).flatten()

			closest_au_h1_idx = au_indices[ np.argmin(distances_h1_to_au ) ]
			closest_au_h2_idx = au_indices[ np.argmin(distances_h2_to_au ) ]

			min_h1_distance_to_au = np.min( distances_h1_to_au )
			min_h2_distance_to_au = np.min( distances_h2_to_au )

			non_hydration_H2O.append({
				"H2O": [ h1_idx, o_idx, h2_idx ],
				f"{h1_idx}-{closest_au_h1_idx}": round(min_h1_distance_to_au, 3),
				f"{h2_idx}-{closest_au_h2_idx}": round(min_h2_distance_to_au, 3)
			})
	sorted_list = sorted( non_hydration_H2O, key=lambda d: list( d.values() )[ 1 ] )
	if to_print == "True":
		print("H2O molecules not in Na hydration shell:")
		for h2o in sorted_list:
			print(h2o)
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

#H2O molecules that belong to NH4 hydration shell
#For shuttling always
#For H2O -> OH + H* if H of H2O is close to electrode
def get_NH4_hydration_shell( poscar, H2O_mols, NH4_molecules, distance_threshold=2.6, to_print = "False" ):
	system = read( poscar )
	au_indices = [ i for i, atom in enumerate( system ) if atom.symbol == "Au" ]
	au_positions = system.positions[ au_indices ]

	final_results = list()

	for nh4 in NH4_molecules:
		n_idx, h1_nh4_idx, h2_nh4_idx, h3_nh4_idx, h4_nh4_idx = nh4
		NH4_H_indices = [ h1_nh4_idx, h2_nh4_idx, h3_nh4_idx, h4_nh4_idx ]

		molecule_results = list()

		for h2o in H2O_mols:
			h1_idx, o_idx, h2_idx = h2o

			closest_nh4_h_to_h2o = None
			min_distance_to_nh4_h = float( 'inf' )
			for h_nh4_idx in NH4_H_indices:
				distance = np.linalg.norm( system.positions[ o_idx ] - system.positions[ h_nh4_idx ] )
				if distance < distance_threshold and distance < min_distance_to_nh4_h:
					closest_nh4_h_to_h2o = h_nh4_idx
					min_distance_to_nh4_h = distance

			if closest_nh4_h_to_h2o is None:
				continue

			distances_h1_to_au = cdist( [ system.positions[ h1_idx ] ], au_positions ).flatten()
			distances_h2_to_au = cdist( [ system.positions[ h2_idx ] ], au_positions ).flatten()

			min_h1_distance_to_au = np.min( distances_h1_to_au )
			min_h2_distance_to_au = np.min( distances_h2_to_au )

			if min_h1_distance_to_au < min_h2_distance_to_au:
				closest_h2o_h_idx = h1_idx
				min_distance_to_au = min_h1_distance_to_au
				closest_au_idx = au_indices[ np.argmin( distances_h1_to_au ) ]
			else:
				closest_h2o_h_idx = h2_idx
				min_distance_to_au = min_h2_distance_to_au
				closest_au_idx = au_indices[ np.argmin( distances_h2_to_au ) ]
			molecule_results.append({
				"[N, H1, H2, H3, H4]": [n_idx, h1_nh4_idx, h2_nh4_idx, h3_nh4_idx, h4_nh4_idx],
				"[H1, O, H2]": [h1_idx, o_idx, h2_idx],
				"[Au]": [closest_au_idx],
				"[H - Au]": [f"{closest_h2o_h_idx} - {closest_au_idx} = {round(min_distance_to_au, 3)}"],
				"[H - O NH4-H2O]": round(min_distance_to_nh4_h, 3),
				"O-H-H_NH4": [o_idx, closest_h2o_h_idx, closest_nh4_h_to_h2o]
			})
		molecule_results.sort(key=lambda x: float( x[ "[H - Au]" ][ 0 ].split( '= ' )[ 1 ] ) )
		final_results.extend( molecule_results )
		if to_print == "True":
			for result in molecule_results:
				print( result )
			print( "\n" )
	return final_results

#H2O molecules that not belong to NH4 hydration shel
#H2O -> OH + H*
def get_non_NH4_hydration(poscar, H2O_mols, NH4_molecules, distance_threshold=2.6, to_print="False"):
	system = read(poscar)
	au_indices = [i for i, atom in enumerate(system) if atom.symbol == "Au"]
	au_positions = system.positions[au_indices]
	non_hydration_H2O = list()

	for h2o in H2O_mols:
		h1_idx, o_idx, h2_idx = h2o
		is_in_hydration_shell = False
		for nh4 in NH4_molecules:
			n_idx, h1_nh4, h2_nh4, h3_nh4, h4_nh4 = nh4
			NH4_H_indices = [h1_nh4, h2_nh4, h3_nh4, h4_nh4]

			for h_nh4_idx in NH4_H_indices:
				distance = np.linalg.norm(system.positions[o_idx] - system.positions[h_nh4_idx])
				if distance <= distance_threshold:
					is_in_hydration_shell = True
					break 

			if is_in_hydration_shell:
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
	print( type( non_hydration_H2O ) )
	sorted_list = sorted( non_hydration_H2O, key=lambda d: list( d.values() )[ 1 ] )
	if to_print == "True":
		print("H2O molecules not in NH4 hydration shell:")
		for h2o in sorted_list:
			print(h2o)

	return non_hydration_H2O

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

#Get H2O mols in the CH3NH3 hydration shell 
#For shuttling
#For H2O -> OH + H* if H of H2O is close to electrode
def get_CH3NH3_hydration_shell(poscar, H2O_mols, CH3NH3_molecules, distance_threshold=2.6, to_print="False"):
	system = read(poscar)
	au_indices = [i for i, atom in enumerate(system) if atom.symbol == "Au"]
	au_positions = system.positions[au_indices]

	final_results = []

	for ch3nh3 in CH3NH3_molecules:
		n_idx, h1_ch3nh3_idx, h2_ch3nh3_idx, h3_ch3nh3_idx, c_idx, h4_ch3nh3_idx, h5_ch3nh3_idx, h6_ch3nh3_idx = ch3nh3
		CH3NH3_H_indices = [h1_ch3nh3_idx, h2_ch3nh3_idx, h3_ch3nh3_idx, h4_ch3nh3_idx, h5_ch3nh3_idx, h6_ch3nh3_idx]

		molecule_results = list()

		for h2o in H2O_mols:
			h1_idx, o_idx, h2_idx = h2o

			closest_ch3nh3_h_to_h2o = None
			min_distance_to_ch3nh3_h = float('inf')
			for h_ch3nh3_idx in CH3NH3_H_indices:
				distance = np.linalg.norm(system.positions[o_idx] - system.positions[h_ch3nh3_idx])
				if distance < distance_threshold and distance < min_distance_to_ch3nh3_h:
					closest_ch3nh3_h_to_h2o = h_ch3nh3_idx
					min_distance_to_ch3nh3_h = distance

			if closest_ch3nh3_h_to_h2o is None:
				continue

			distances_h1_to_au = cdist([system.positions[h1_idx]], au_positions).flatten()
			distances_h2_to_au = cdist([system.positions[h2_idx]], au_positions).flatten()

			min_h1_distance_to_au = np.min(distances_h1_to_au)
			min_h2_distance_to_au = np.min(distances_h2_to_au)

			if min_h1_distance_to_au < min_h2_distance_to_au:
				closest_h2o_h_idx = h1_idx
				min_distance_to_au = min_h1_distance_to_au
				closest_au_idx = au_indices[np.argmin(distances_h1_to_au)]
			else:
				closest_h2o_h_idx = h2_idx
				min_distance_to_au = min_h2_distance_to_au
				closest_au_idx = au_indices[np.argmin(distances_h2_to_au)]

			molecule_results.append({
				"[N, H1, H2, H3, C, H4, H5, H6]": [n_idx, h1_ch3nh3_idx, h2_ch3nh3_idx, h3_ch3nh3_idx, c_idx, h4_ch3nh3_idx, h5_ch3nh3_idx, h6_ch3nh3_idx],
				"[H1, O, H2]": [h1_idx, o_idx, h2_idx],
				"[Au]": [closest_au_idx],
				"[H - Au]": [f"{closest_h2o_h_idx} - {closest_au_idx} = {round(min_distance_to_au, 3)}"],
				"[H - O CH3NH3-H2O]": round(min_distance_to_ch3nh3_h, 3),
				"O-H-H_CH3NH3": [o_idx, closest_h2o_h_idx, closest_ch3nh3_h_to_h2o]
			})

		molecule_results.sort(key=lambda x: float(x["[H - Au]"][0].split('= ')[1]))
		final_results.extend(molecule_results)

		if to_print == "True":
			for result in molecule_results:
				print(result)
			print("\n") 

	return final_results

#H2O molecules that not belong to CH3NH3 hydration shel
#H2O -> OH + H*
def get_non_CH3NH3_hydration_shell( poscar, H2O_mols, CH3NH3_molecules, distance_threshold = 2.6, to_print = "False" ):
	system = read( poscar )
	au_indices = [ i for i, atom in enumerate( system ) if atom.symbol == "Au" ]
	au_positions = system.positions[ au_indices ]

	non_hydration_H2O = list()

	for h2o in H2O_mols:
		h1_idx, o_idx, h2_idx = h2o
		is_in_hydration_shell = False
		for ch3nh3 in CH3NH3_molecules:
			n_idx, h1_ch3nh3, h2_ch3nh3, h3_ch3nh3, c_idx, h4_ch3nh3, h5_ch3nh3, h6_ch3nh3 = ch3nh3
			CH3NH3_H_indices = [ h1_ch3nh3, h2_ch3nh3, h3_ch3nh3, h4_ch3nh3, h5_ch3nh3, h6_ch3nh3 ]
			for h_ch3nh3_idx in CH3NH3_H_indices:
				distance = np.linalg.norm( system.positions[ o_idx ] - system.positions[ h_ch3nh3_idx ] )
				if distance <= distance_threshold:
					is_in_hydration_shell = True
					break  

			if is_in_hydration_shell:
				break  
		if not is_in_hydration_shell:
			distances_h1_to_au = cdist( [ system.positions[ h1_idx ] ], au_positions ).flatten()
			distances_h2_to_au = cdist( [ system.positions[ h2_idx ] ], au_positions ).flatten()

			closest_au_h1_idx = au_indices[ np.argmin( distances_h1_to_au ) ]
			closest_au_h2_idx = au_indices[ np.argmin(distances_h2_to_au ) ]

			min_h1_distance_to_au = np.min( distances_h1_to_au )
			min_h2_distance_to_au = np.min( distances_h2_to_au )
			non_hydration_H2O.append({
				"H2O": [ h1_idx, o_idx, h2_idx ],
				f"{h1_idx}-{closest_au_h1_idx}": round( min_h1_distance_to_au, 3 ),
				f"{h2_idx}-{closest_au_h2_idx}": round( min_h2_distance_to_au, 3 )
			})
	sorted_list = sort( non_hydration_H2O, key = lambda d: list( d.values()[ 1 ] ) )
	if to_print == "True":
		print("H2O molecules not in CH3NH3 hydration shell:")
		for h2o in sorted_list:
			print( h2o )

	return non_hydration_H2O


if __name__ == "__main__":
	H2O_mols = get_H2O_mols( "POSCAR" )
	Na_mols = get_Na_mols( "POSCAR" )
	get_non_Na_hydration_shell( "POSCAR", H2O_mols, Na_mols, to_print = "True" )
