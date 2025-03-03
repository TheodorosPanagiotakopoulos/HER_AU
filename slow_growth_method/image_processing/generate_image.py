import data
import os
import matplotlib.pyplot as plt
from ase.io.trajectory import Trajectory
from ase.io import read
import json as js
import numpy as np
from figure import get_figure_2

os.system( "rm -rf frame_*" )
loc = os.getcwd() + "/"
dirs = data.get_dirs()
dirs = [ loc + x for x in dirs ]
count = -1
step  = 25
start = -1

idx1, idx2, idx3 = data.get_standarized_ICONST_data( loc )
print( idx1, idx2, idx3 )

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
            ref = img[ 64 ].position * 1.
            #img[ 103 ].symbol = 'S'
            #img[ 181 ].symbol = 'He'
            #img[ 178 ].symbol = 'He'
            img[ 189 ].symbol = 'He'
            img[ 190 ].symbol = 'He'
            img[ 191 ].symbol = 'He'
            img[ 192 ].symbol = 'He'
            
            img[ 194 ].symbol = 'Ge'
            img[ 195 ].symbol = 'Ge'
            img[ 196 ].symbol = 'Ge'
            img[ 197 ].symbol = 'Ge'
            
            img[ 184 ].symbol = 'Li'
            img[ 185 ].symbol = 'Li'
            img[ 186 ].symbol = 'Li'
            img[ 187 ].symbol = 'Li'
            
            img[ idx1 ].symbol = 'He'
            #img[ 120 ].symbol = 'He'
            #img[ 74 ].symbol = 'S'
            if system[ idx2 ].symbol != "N": 
                img[ idx2 ].symbol = 'S'

            center = (img.cell[ 0, : ] + img.cell[ 1, : ] ) /2.
            center[ 2 ] = 0.
            img.positions[ :, 0 ] -= ref[ 0 ]
            img.positions[ :, 1 ] -= ref[ 1 ]
            img.positions[ :, 0 ] += center[ 0 ]
            img.positions[ :, 1 ] += center[ 1 ]
            img.wrap()
            prefix = str( count )
            fout = 'frame_' + prefix
            if not os.path.isfile( fout + '.png' ):
                get_figure_2( img, fout, rot = "180z,-80x", w = 13, h = 12 )
            os.system( 'rm -f *.ini *.pov' )
            '''
            img[79].symbol = 'S'
            img[80].symbol = 'S'
            ref = img[ 98 ].position * 1.
            center = (img.cell[ 0, : ] + img.cell[ 1, : ] ) /2.
            center[ 2 ] = 0.
            img.positions[ :, 0 ] -= ref[ 0 ]
            img.positions[ :, 1 ] -= ref[ 1 ]
            img.positions[ :, 0 ] += center[ 0 ]
            img.positions[ :, 1 ] += center[ 1 ]
            #img.wrap()
            prefix = str( count )
            fout = 'frame_' + prefix
            if not os.path.isfile( fout + '.png' ):
                get_figure_2( img, fout, rot = "0z,-90x", w = 20, h = 20 )
            os.system( 'rm -f *.ini *.pov' )
            '''
