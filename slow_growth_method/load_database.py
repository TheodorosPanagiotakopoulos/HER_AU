import json

def load_database_to_dict( database_file ):
    with open( database_file, "r" ) as file:
        database = json.load( file )
    return database

def get_values( val ):
    for key, value in  database[ val ].items():
        print( key, value )

if __name__ == "__main__":
    file = "database_for_theo.js"
    database = load_database_to_dict( file )
    for key in database.keys():
        print( key )

    #for value in database[ "1_CH3NH3_H2O_dissociation_NOT_from_hydration_shell" ]:
    #   print(  value )

    get_values( "3_CH3NH3_H2O_dissociation_from_hydration_shell" )
