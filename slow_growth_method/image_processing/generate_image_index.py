import os
import matplotlib.pyplot as plt
from ase.io.trajectory import Trajectory
from ase.io import read
import json as js
import numpy as np
from figure import get_figure_2
import data

os.system( "rm -rf frame_*" )
images = range( 20000 ) #[  0, 1050, 2000 ]
loc = os.getcwd() + "/"
dirs = data.get_dirs()
dirs = [ loc + x for x in dirs ]
count = -1
step  = 25 #00
start = -1

print( dirs )
### Indices of atoms to bi highlighted should be in this section

H_shuttle = [ ]  #Index of H that shuttle. Will be in cyan
H_highlight = [ 184 ] # Index of H to be highlighted
O_highlight = [ ]     # Index of O to be highlighted
N_highlight = [ 188 ]    # 188 Index of N to be highlighted
center_atom = 91   #x, y of this atom will be shift to 0,0, it will be at the center of image

###

system = read( loc + '/RUN1/POSCAR' )
for dir in dirs:
	with open( dir + '/DATA.js', 'r' ) as f:
		data = js.load( f )
	npoints = len( data.keys() )
	traj = Trajectory( dir + '/MOVIE.traj', 'r' )
	for i in range( len( traj ) ):
		count += 1
		if count > start and count % step == 0:
			print( count )
			img = traj[ i ]
			ref = img[ center_atom ].position * 1.
			ref[ 1 ] += 1.
			for at in O_highlight:
				img[ at ].symbol = 'S'
			for at in H_shuttle:
				img[ at ].symbol = 'Ne'
			for at in H_highlight:
				img[ at ].symbol = 'He'
			for at in N_highlight:
				img[ at ].symbol = 'P' #P

			center = (img.cell[ 0, : ] + img.cell[ 1, : ] ) /2.
			center[ 2 ] = 0.
			img.positions[ :, 0 ] -= ref[ 0 ]
			img.positions[ :, 1 ] -= ref[ 1 ]
			img.positions[ :, 0 ] += center[ 0 ]
			img.positions[ :, 1 ] += center[ 1 ]
			img.wrap()
			if count in images:
				prefix = str( count )
				fout = 'frame_' + prefix
				if not os.path.isfile( fout + '.png' ):
					get_figure_2( img, fout, rot = "180z,-80x,10y", w = 11, h = 11 )
			os.system( 'rm -f *.ini *.pov' )
