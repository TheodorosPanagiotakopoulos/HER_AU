import os
import tppy
from ase.io import read, write
from ase.visualize import view

atoms =  [ 85, 144, 191 ]
atoms_aux = [ [ atoms[ 0 ] , atoms[ 1 ] ], [ atoms[ 0 ] , atoms[ 2 ] ] ]

path_to_MD_CONTCAR = "//home/theodoros/PROJ_ElectroCat/theodoros/HER/Au/HER_Au/MD_voltage/NH4/5_NH4_40_H2O/RUN11/"
struc = read( path_to_MD_CONTCAR + "CONTCAR" )

#struc.wrap()
tppy.get_ICONST( atoms_aux, " "  )
struc = tppy.shift_center( struc, atoms[ 0 ] )

view( struc )
write ( 'POSCAR', struc , direct=True, format='vasp', sort = False )
struc = tppy.add_velocity( path_to_MD_CONTCAR )
os.system( "echo PBE | ~/.potcar.j" )

#tppy.change_atomic_number( struc, atoms )
