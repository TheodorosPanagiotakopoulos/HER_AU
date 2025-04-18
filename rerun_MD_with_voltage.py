import os
import re
import glob
import shutil

def get_path_and_dirs():
    path = os.getcwd()
    dirs = [ i for i in os.listdir( "." ) if os.path.isdir( i ) == True  ]
    return path, dirs

def make_dirs_and_copy():
    dirs = [ i for i in os.listdir( "." ) if os.path.isdir( i ) == True  ]
    if dirs:
        new_dir_name = f"RUN{len(dirs) + 1 }"
    else:
        new_dir_name = "RUN1"
    os.mkdir( new_dir_name )
    new_dirs = [ i for i in os.listdir( "." ) if os.path.isdir( i ) == True  ]
    new_dirs.sort( key = os.path.getctime )
    target_dir =  new_dirs[ - 1 ]
    files = [ i for i in glob.glob( "*" ) if os.path.isfile( i ) ]
    for i in files:
        if i in [ "OUTCAR", "REPORT", "vasprun.xml", "XDATCAR" ]:
            os.system( f"gzip {i}" )
    files = [ i for i in glob.glob( "*" ) if os.path.isfile( i ) ]
    for i in files:
        shutil.copy( i, target_dir )
    for i in [ "OUTCAR.gz", "REPORT.gz", "vasprun.xml.gz", "XDATCAR.gz", "POSCAR" ]:
        os.remove( i )
    os.rename( "CONTCAR", "POSCAR" )

def get_OUTCAR_NELECT():
    files = [ i for i in glob.glob( "*" ) if os.path.isfile( i ) ]
    if "OUTCAR" in files:
        print( "OUTCAR found" )
        with open( "OUTCAR", "rb" ) as file:
            lines = file.readlines()
        nelect_lines = [line.decode('utf-8', errors='ignore').strip() for line in lines if b"NELECTCURRENT" in line]
        last_line = nelect_lines[ - 1 ].strip()
        if nelect_lines:
            match =  re.search( r'NELECTCURRENT\s+([\d.-]+)', last_line )
            #print( "NELECTCURRENT = ", match.group( 1 )  )
            return  float( match.group( 1 ) )
        else:
            print( "NELECT not found yet" )
            return
    else:
        print( "OUTCAR not found" )
        return

def update_INCAR( val ):
    with open( "INCAR", 'r') as file:
        lines = file.readlines()

    NELECT_present = False
    for i, line in enumerate( lines ):
        if line.startswith( 'NELECT' ):
            lines[ i ] = f'NELECT = {val}\n'
            NELECT_present = True
            break

    if not NELECT_present:
        lines.append(f'NELECT = {val}\n')

    with open( "INCAR" , 'w' ) as file:
        file.writelines( lines )

def parse_INCAR():
    path, dirs = get_path_and_dirs()
    for i in dirs:
        os.chdir( path + "/" + i  )
        path_new = os.getcwd()
        print( "-" * len( path_new ) )
        print( os.getcwd()  )
        NELECTCURRENT = get_OUTCAR_NELECT()
        print( "nelect =", NELECTCURRENT )
        make_dirs_and_copy()
        update_INCAR( NELECTCURRENT )

def submit_job():
    cpath = os.getcwd()
    parse_INCAR()
    os.chdir( cpath )
    os.system( "sbatch job" )


#run in the directory where the job has stopped
if __name__ == "__main__":
    nelect = get_OUTCAR_NELECT()
    print( nelect ) 
    make_dirs_and_copy()
    update_INCAR( nelect )
