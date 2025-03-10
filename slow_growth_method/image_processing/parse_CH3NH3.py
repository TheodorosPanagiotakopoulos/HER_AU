import os
import sys

images = []
data = {}
paths =  [x for x in os.listdir( '.' ) if 'RUN' in x ] + ['./']
#paths =  [x for x in os.listdir( '.' ) if 'RUN' in x ]
ionspertype = ''
if len( sys.argv ) > 1:
    ionspertype = sys.argv[1]

print( "ions per type = " + ionspertype )
for path in paths:
    if not os.path.isfile( path + '/MOVIE.traj' ):
        images = []
        data = {}
        outcar = path + '/OUTCAR.gz'
        images, data = read_OUTCAR( outcar , maxiter = 20000, images = images, data = data, ionspertype = ionspertype  )
        traj = Trajectory( path + '/MOVIE.traj', 'w')
        add_to_traj( traj=traj, images=images, first = 0, last = 20000)
        with open( path + '/DATA.js', 'w' ) as f:
            js.dump( data, f )
