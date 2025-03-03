from ase.io          import read, write
from dlePy.supercell import supercell, remove_atoms_outside
import numpy as np
import sys
from ase.io.pov import get_bondpairs, get_hydrogenbonds, set_high_bondorder_pairs
from ase.visualize import view
import matplotlib.cm as cm
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

def get_atom_setting( ):
    with open( '/shared/apps/VESTA-x86_64/elements.ini', 'r' ) as f:
        lines = f.readlines( )
    elements = { }
    for line in lines:
        tmp = line.split( )
        color = ( float( tmp[ 5 ] ), float( tmp[ 6 ] ), float( tmp[ 7 ] ) )
        rcov   = float( tmp[ 2 ] )
        rvdw   = float( tmp[ 3 ] )
        r   = float( tmp[ 4 ] )
        element = tmp[ 1 ]
        elements[ element ] = { }
        elements[ element ][ 'r'    ]  = r
        elements[ element ][ 'rvdw' ]  = rvdw
        elements[ element ][ 'rcov' ]  = rcov
        elements[ element ][ 'color' ] = color
    return elements

def get_figure( sys0, fout, rot = "-30x", w = 12, h = 11, alpha_top = 0.5, alpha_bot = 0.5 ):
    system  = supercell( sys0, 1, 3, 1 )
    center = 0.5 * system.cell[ 0 ] + 0.5* system.cell[ 1 ] + 0.5* system.cell[ 2 ]
    for i in range( 3 ):
        system.positions[ :, i ] -= center[ i ]
    system_cp = system.copy()
    dels = []
    for at in system:
        if at.position[2] >  0.5:
            dels.append( at.index )

    for at in sorted( dels, reverse = True ):
        del ( system[ at ] )

    get_figure_2( system, fout + '_bot', rot = rot, w = w, h = h )

    system = system_cp.copy()
    dels = []
    for at in system:
        if at.position[2] >  0.4 or  at.position[2] < -2.:
            dels.append( at.index )

    for at in sorted( dels, reverse = True ):
        del ( system[ at ] )

    get_figure_2( system, fout + '_center', rot = rot, w = w, h = h )

    system = system_cp.copy()

    dels = []
    for at in system:
        if at.position[2] >  2.5 or  at.position[2] < 0.5:
            dels.append( at.index )

    for at in sorted( dels, reverse = True ):
        del ( system[ at ] )

    get_figure_2( system, fout + '_top', rot = rot, w = w, h = h )

    fig = plt.figure( figsize = ( 3.5, 3.5 ) )
    ax = fig.add_subplot( 1, 1, 1 )
    img = mpimg.imread( fout + '_bot.png' )
    ax.imshow( img, zorder = 0 )
    print ( np.max( img ) )
    img[ :, :, : ] = 0. 
    ax.imshow( img[:,:,0], vmin = 0, vmax = 1, cmap = cm.gray, alpha = alpha_bot, zorder = 1 )
    img = mpimg.imread( fout + '_center.png' )
    ax.imshow( img, zorder = 2 )
    img = mpimg.imread( fout + '_top.png' )
    ax.imshow( img, alpha = alpha_top, zorder = 3 )
    ax.axis( 'off' )
    plt.subplots_adjust(left=0.0, right=1., top=1., bottom=0.0, wspace=0.00, hspace= 0.00)
    plt.savefig( fout + '.png', dpi = 600 )
    

def get_figure_2( sys0, fout, rot = "-30x", w = 15, h = 15 ):
    system  = supercell( sys0, 3, 1, 1 )
    center = 1.8*sys0.cell[ 0 ] + 0.5*sys0.cell[ 1 ] + 0.5* sys0.cell[ 2 ]	
    #center = 1.5*sys0.cell[ 0 ] + 0.5*sys0.cell[ 1 ] + 0.5* sys0.cell[ 2 ]
    #center = 0.5*sys0.cell[ 0 ] + 0.5*sys0.cell[ 1 ] + 0.5* sys0.cell[ 2 ]
    #center[0] += sys0.cell[ 0, 0 ]
    #center[2] -= 3.
    #center[2] += 3
    #center[1] += 10
    for i in range( 3 ):
        system.positions[ :, i ] -= center[ i ]
    
    elements = get_atom_setting( )
    elements[ 'S' ][ 'color' ] = ( 1, 0, 0 )
    elements[ 'O' ][ 'color' ] = ( 1, 0.7, 0.7 )
    elements[ 'H' ][ 'color' ] = ( 1, 1, 1 )
    elements[ 'He' ][ 'color' ] = ( 1, 0, 1 )
    elements[ 'S' ][ 'rcov' ] = elements[ 'O' ][ 'rcov' ]
    elements[ 'He' ][ 'rcov' ] = elements[ 'H' ][ 'rcov' ]
    elements[ 'C' ][ 'rcov' ] = elements[ 'O' ][ 'rcov' ]
    elements[ 'Na' ][ 'rcov' ] = elements[ 'Na' ][ 'rcov' ] * 0.6

    radii = [ ]
    colors = [ ]
    for i, at in enumerate( system ):
        radius = elements[ at.symbol ][ 'rcov' ]
        if at.symbol in [ 'O', 'H', 'S', 'He', 'N' ]:
            radius = radius/3.
        radii.append( radius )
        color = elements[ at.symbol ][ 'color' ]
        alpha = 0
        if at.symbol in [ 'C', 'S', 'K' ]:
            alpha =1 
        factor = 0
        colors.append( ( factor + color[ 0 ], factor + color[ 1 ], factor + color[ 2 ], alpha ) )


    bond_pairs = get_bondpairs( system, radius = 0.6 )
    hydrogenbond1 = get_hydrogenbonds( system, atype1 = 'Na', atype2= [ 'O', 'S' ], radius = 4, rhbondrange = (2.3, 4.0) )

    hydrogenbond2 = get_hydrogenbonds( system, atype1 = 'H', atype2= [ 'O', 'S', 'N' ], radius = 5, rhbondrange = (1.3, 2.6) )
    hydrogenbond3 = get_hydrogenbonds( system, atype1 = 'He', atype2= [ 'O', 'S', 'N' ], radius = 5, rhbondrange = (1.3, 2.6) )
    #hydrogenbond2 = get_hydrogenbonds( system, atype1 = 'K', atype2= [ 'O', 'S' ], radius = 2, rhbondrange = (1.3, 3.0) )
    hydrogenbond  = {}
    for key in hydrogenbond1.keys():
        hydrogenbond[ key ] = hydrogenbond1[ key ]
    for key in hydrogenbond2.keys():
        hydrogenbond[ key ] = hydrogenbond2[ key ]
    for key in hydrogenbond3.keys():
        hydrogenbond[ key ] = hydrogenbond3[ key ]

    bond_pairs = set_high_bondorder_pairs(bond_pairs, high_bondorder_pairs=hydrogenbond)

    width  = w
    height = h
    bbox = (
           -width/2,
           -height/2-2.5,
           width/2,
           height/2-2.5,
           )
    #for rotx in [ 0, 15, 30, 45, 60, 75, 90 ]:
    #    rot = '90z,' + str( rotx ) + 'x' 
    write( fout+ '.pov', system, format = 'pov', run_povray = True,
           canvas_width = 1000,    # Set width, in pixel
           radii = radii,              # Set radius 
           bondatoms = bond_pairs, # Display bonds
           bbox  = bbox,
           colors = colors,        # Set colors
           celllinewidth = 0.0,
           rotation = rot,
           hydrogenbond = { 'ndots':5, 'color' : [0., 1, 0.], 'rdot':0.05 }
           )
    

if __name__ == "__main__":
    #system = read( 'TS_POSCAR' )
    #for rotx in [ 0, 15, 30, 45, 60, 75, 90 ]:
    #    rot = '90z,-' + str( rotx ) + 'x' 
    #    get_figure( system, 'TS_' + str( rotx ) , rot = rot, h = 15, w = 15, alpha_top = 0.80, alpha_bot = 0.5 ) 
    for poscar in [ 'IS' , 'TS', 'FS' ]:
        system = read( poscar + '_POSCAR' )
        rot = '90z,0x'
        get_figure( system, poscar, rot = rot, h = 15, w = 15, alpha_top = 0.80, alpha_bot = 0.5 )
