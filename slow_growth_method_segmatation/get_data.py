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

def get_cc_bm_in_df():
	reports = [f"REPORT.{i}" for i in range(1, 51)]
	df = pd.DataFrame()
	for i in reports:
		print( "Parsing" + i )
		#print(f"Parsing REPORT.{i}")
		with open( i, "rb" ) as file:
			lines = file.readlines()
		cc = [ line.decode( "utf-8", errors = "ignore" ).strip() for line in lines if b"cc" in line ]	
		b_m = [ line.decode( "utf-8", errors = "ignore" ).strip() for line in lines if b"b_m" in line ]
		df_cc = pd.DataFrame( [ row.split() for row in cc ] )[ 3 ]
		df_bm = pd.DataFrame( [ row.split() for row in b_m ] )[ 1 ]
		merged_df = pd.concat( [ df_cc, df_bm ], axis = 1 )
		df = pd.concat([df, merged_df], ignore_index = True )
	df.rename(columns={ 3 : "cc", 1 : 'b_m'}, inplace=True)
	print("Final DataFrame:")
	#print( df.to_string() )
	for i in [ "cc", "b_m" ]:
		df[ i ] = df[ i ].astype(float)
	print( type( df["cc"][ 1 ] ))
	return df

def get_free_energy( path_to_SG_calculation ):
	os.chdir( path_to_SG_calculation )
	tg = [ 0.0 ]
	df = get_cc_bm_in_df()
	cc = df[ "cc" ]
	b_m = df[ "b_m" ]
	for i in range( 1, len( cc ) ):
		gg = 0.5 * ( cc[ i ]  -  cc[ i - 1 ] ) * ( b_m[ i ]  +  b_m[ i - 1 ] )
		tg.append( tg[ -1 ] + gg )
	print( len( tg ) )
	print( len( cc ) )
	return cc, tg

if __name__ == "__main__":
	path = "/home/theodoros/PROJ_ElectroCat/theodoros/HER/Au/HER_Au/slow_grow_method/tutorial/"
	cc, tg = get_free_energy( path )
