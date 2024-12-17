import os
import gzip
import glob
import sys
import pandas as pd
import numpy as np

def flatten_matrix( matrix ):
	flat_list = list()
	for i in matrix:
		flat_list.extend( i )
	return flat_list

def get_RUNs():
		runs = glob.glob( "RUN*" )
		if runs:
			return runs
		else:
			return []

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

def get_free_energy( path_to_SG_calculation ):
	tg = [ 0.0 ]
	cc, b_m = collect_cc_and_bm( path_to_SG_calculation )
	if ( not cc ) or ( not b_m ):
		return None, None
	for i in range( 1, len( cc ) ):
		gg = 0.5 * ( cc[ i ]  -  cc[ i - 1 ] ) * ( b_m[ i ]  +  b_m[ i - 1 ] )
		tg.append( tg[ -1 ] + gg )
	return cc, tg

def get_barriers( path ):
	#do not allow printing from get free_energy() 
    original_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    try:
        x, y = get_free_energy( path )
        if x == None or y == None:
            return None
        else:
            barrier = round( max( y ) - min( y ), 2 )
            return barrier
    finally:
            sys.stdout.close()
            sys.stdout = original_stdout

if __name__ == "__main__":
	path = "/home/theodoros/PROJ_ElectroCat/theodoros/HER/Au/HER_Au/slow_grow_method/Na/1_Na/free_H2O_splitting/1_Na_40_H2O_v3"
	#df1, df2 = get_cc_bm()
	x, y = get_free_energy( path )
	diff = round( max( y ) - min( y ), 2 )
	print( diff )
