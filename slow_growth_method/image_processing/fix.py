import string

def is_printable(line):
    if not line.strip():
        return True
    threshold = 0.8
    printable_chars = set(string.printable)
    num_printable = sum(1 for char in line if char in printable_chars)
    return (num_printable / max(len(line), 1)) >= threshold

def clean_outcar(filename="OUTCAR"):
    try:
        with open(filename, "r", encoding="latin-1", errors="replace") as file:
            lines = file.readlines()

        keywords = {"SIGMA_RC_K", "CORE_C"}
        removed_anything = False

        for i in range( len( lines ) - 1 ):
            if any( keyword in lines[ i ] for keyword in keywords ):
                next_line = lines[ i + 1 ].strip()

                if is_printable(next_line):
                    print( "No changes to be done" )
                else:
                    print( f"Removing: {repr(next_line)}" )
                    lines[ i + 1 ] = ""
                    removed_anything = True

        if removed_anything:
            with open( filename, "w", encoding = "latin-1" ) as cleaned_file:
                cleaned_file.writelines( lines )
            print("\n Cleaned file saved in 'OUTCAR' (original file overwritten).")
        else:
            print("\n Nothing to remove. 'OUTCAR' remains unchanged." )
    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found." )
    except Exception as e:
        print(f"An error occurred: {e}" )

clean_outcar()
