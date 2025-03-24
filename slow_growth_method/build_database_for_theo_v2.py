import json as js
import os
import getpass
import get_mols


username=getpass.getuser()

keywords = [ 'path', 'note' ]
def db_add_key( database, key ):
	if key in database.keys():
		print ( "key %s already exit" %(key) )
	else:
		database[ key ] = {}

def db_add( database, key, name, input ):
	if name in database[ key ].keys():
		print ( "ERROR: %s already exit in %s" %(str( name ),  key) )
		exit()
	if 'path' not in input.keys():
		print ( "ERROR: path is required %s" %(input) )
		exit()
    
	database[ key ][ name ] = input 

def db_update( database, key, name, input ):
	if 'path' not in input.keys():
		print ( "ERROR: path is required %s" %(input) )
		exit()
	if name not in database[ key ].keys():
		print ( "ERROR: %s does not exit in %s" %(str( name ),  key) )
		exit()

	if key not in database.keys():
		print ( "ERROR: %s does not exit" %( key) )
		exit()

	database[ key ][ name ] = input

def db_save( database, fout ):
	with open( fout, 'w' ) as f:
		js.dump( database, f, indent = 4 )

def db_load( database, fin ):
	with open( fin, 'r' ) as f:
		database  = js.load( f )

def convert( path ):
	return path.replace( '~', '/home/' + username )

def add_to_database( loc, key, folders, note = "" ):
	db_add_key( database, key )
	for name in [ x for x in os.listdir( convert( loc ) ) if folders in x ]:
		path = loc + name
		loc_ = convert( loc )  +  name 
		runs = get_mols.get_RUNs( loc_ )
		loc_status = loc_ + "/RUN" + str( len( runs ) )
		loc_contcar = loc_ + "/RUN" + str( len( runs ) ) + "/CONTCAR"
		if not runs:
			loc_contcar = loc_ + "/CONTCAR" 
			loc_status = loc_ 
		note = "Bad"
		if ( os.path.isfile( loc_contcar ) ) == True:
			status = get_mols.get_status( loc_status )
			if status == True:
				note = "Good"
			input = { 'path': loc_, 'note': note }
			db_add( database, key, name, input )
		else:
			input = { 'path': loc_, 'note': note }
			db_add( database, key, name, input )
	db_save( database, 'database_for_theo.js' )

#for relaxed systems and MD with and without Voltage for reaching stability
def add_to_database_v2( loc, key ):
	db_add_key( database, key )
	for name in [ x for x in os.listdir( convert( loc ) ) if os.path.isdir( convert( loc ) +  "/" + x )  ]:
		path = loc + name
		input = { 'path': path, 'note': "" }
		db_add( database, key, name, input )
	db_save( database, 'database_for_theo.js' )

if __name__ == "__main__":
	database = {}        
	
	K_relax = "/home/theodoros/PROJ_ElectroCat/theodoros/HER/Au/HER_Au/relaxation/K"
	add_to_database_v2( K_relax, "K_relax" )	

	Na_relax = "/home/theodoros/PROJ_ElectroCat/theodoros/HER/Au/HER_Au/relaxation/Na"
	add_to_database_v2( Na_relax, "Na_relax" )

	NH4_relax = "/home/theodoros/PROJ_ElectroCat/theodoros/HER/Au/HER_Au/relaxation/NH4"
	add_to_database_v2( NH4_relax, "NH4_relax" )

	CH3NH3_relax = "/home/theodoros/PROJ_ElectroCat/theodoros/HER/Au/HER_Au/relaxation/CH3NH3"
	add_to_database_v2( CH3NH3_relax, "CH3NH3_relax" )

	Na_MD = "/home/theodoros/PROJ_ElectroCat/theodoros/HER/Au/HER_Au/MD/Na"
	add_to_database_v2( Na_MD, "Na_MD" )

	NH4_MD = "/home/theodoros/PROJ_ElectroCat/theodoros/HER/Au/HER_Au/MD/NH4"
	add_to_database_v2( NH4_MD, "NH4_MD" )

	CH3NH3_MD = "/home/theodoros/PROJ_ElectroCat/theodoros/HER/Au/HER_Au/MD/CH3NH3"
	add_to_database_v2( CH3NH3_MD, "CH3NH3_MD" )

	Na_MD_voltage = "/home/theodoros/PROJ_ElectroCat/theodoros/HER/Au/HER_Au/MD_voltage/Na"
	add_to_database_v2( Na_MD_voltage, "Na_MD_voltage" )

	NH4_MD_voltage = "/home/theodoros/PROJ_ElectroCat/theodoros/HER/Au/HER_Au/MD_voltage/NH4"
	add_to_database_v2( NH4_MD_voltage, "NH4_MD_voltage" )

	CH3NH3_MD_voltage = "/home/theodoros/PROJ_ElectroCat/theodoros/HER/Au/HER_Au/MD_voltage/CH3NH3"
	add_to_database_v2( CH3NH3_MD_voltage, "CH3NH3_MD_voltage" )


	base_Na = "~/PROJ_ElectroCat/theodoros/HER/Au/HER_Au/slow_grow_method/Na/"
	Na_NashMap = { "A" : "/free_H2O_splitting/",  "B" :  "/H2O_from_hydration_shell_splitting/" }
	Na_key_HashMap = { "A" : "_Na_H2O_dissociation_NOT_from_hydration_shell", "B" : "_Na_H2O_dissociation_from_hydration_shell" }
	Na_name = "_Na_40_H2O_v"
	for i in [ 1, 3, 5 ]:
		add_to_database( base_Na + str( i ) + "_Na" + Na_NashMap[ "A" ], str( i ) + Na_key_HashMap[ "A" ],  str( i ) + Na_name )
		add_to_database( base_Na + str( i ) + "_Na" + Na_NashMap[ "B" ], str( i ) + Na_key_HashMap[ "B" ],  str( i ) + Na_name )

	base_NH4 = "~/PROJ_ElectroCat/theodoros/HER/Au/HER_Au/slow_grow_method/NH4/"
	NH4_HashMap = { "A" : "/H2O_splitting_from_NH4_hydration_shell/",  "B" :  "/H2O_splitting_NOT_from_NH4_hydration_shell/", "C" : "/NH4_splitting/", "D" : "/shuttling/" }
	NH4_key_HashMap = { "A" : "_NH4_H2O_dissociation_from_hydration_shell", "B" : "_NH4_H2O_dissociation_NOT_from_hydration_shell", "C" :"_NH4_spliting", "D" : "_NH4_shuttling" }
	NH4_name = "_NH4_40_H2O_v"
	for i in [ 1, 3, 5 ]:
		add_to_database( base_NH4 + str( i ) + "_NH4" + NH4_HashMap[ "A" ], str( i ) + NH4_key_HashMap[ "A" ],  str( i ) + NH4_name )
		add_to_database( base_NH4 + str( i ) + "_NH4" + NH4_HashMap[ "B" ], str( i ) + NH4_key_HashMap[ "B" ],  str( i ) + NH4_name )
		add_to_database( base_NH4 + str( i ) + "_NH4" + NH4_HashMap[ "C" ], str( i ) + NH4_key_HashMap[ "C" ],  str( i ) + NH4_name )
		add_to_database( base_NH4 + str( i ) + "_NH4" + NH4_HashMap[ "D" ], str( i ) + NH4_key_HashMap[ "D" ],  str( i ) + NH4_name )

	base_CH3NH3 = "~/PROJ_ElectroCat/theodoros/HER/Au/HER_Au/slow_grow_method/CH3NH3/"
	CH3NH3_HashMap = { "A" : "/H2O_splitting_from_CH3NH3_hydration_shell/",  "B" :  "/H2O_splitting_NOT_from_CH3NH3_hydration_shell/", "C" : "/CH3NH3_splitting/", "D" : "/shuttling/" }
	CH3NH3_key_HashMap = { "A" : "_CH3NH3_H2O_dissociation_from_hydration_shell", "B" : "_CH3NH3_H2O_dissociation_NOT_from_hydration_shell", "C" :"_CH3NH3_spliting", "D" : "_CH3NH3_shuttling" }
	CH3NH3_name = "_CH3NH3_40_H2O_v"
	for i in [ 1, 3, 5 ]:
		add_to_database( base_CH3NH3 + str( i ) + "_CH3NH3" + CH3NH3_HashMap[ "A" ], str( i ) + CH3NH3_key_HashMap[ "A" ],  str( i ) + CH3NH3_name )
		add_to_database( base_CH3NH3 + str( i ) + "_CH3NH3" + CH3NH3_HashMap[ "B" ], str( i ) + CH3NH3_key_HashMap[ "B" ],  str( i ) + CH3NH3_name )
		add_to_database( base_CH3NH3 + str( i ) + "_CH3NH3" + CH3NH3_HashMap[ "C" ], str( i ) + CH3NH3_key_HashMap[ "C" ],  str( i ) + CH3NH3_name )
		add_to_database( base_CH3NH3 + str( i ) + "_CH3NH3" + CH3NH3_HashMap[ "D" ], str( i ) + CH3NH3_key_HashMap[ "D" ],  str( i ) + CH3NH3_name )
