from database.sequence_database import SequenceDatabase
from database.load_database import create_database

def normalize_database(db):

    """
    Allows for flexibility in the input provided as a Database, such as:

    FASTA File sequences: 
    
    -Detects that a FILE was passed, processes it using FASTA-specific parsing and extracts each sequence within the FASTA file into a new SequenceDatabase object;
    Again, this maintains consistency of utilizing one Data type.

    Sequences stored in a list:

    - Detects that the input DB is in the form of a list data structure, then iterates through it, adding each sequence into a new SequenceDatabase object before returning this Object to continue with the pipeline's processes.

    SequenceDatabase:

    - If input is already in the form of a SequenceDatabase object, it detects this and just returns the object itself.

    If none of the above:

    - Raises a TypeError to indicate that either the pipeline has not yet accommodated for the passed data structure, or that the user has passed an invalid Database format. 

    """
    if isinstance(db, SequenceDatabase):
        return db
    
    elif isinstance(db, dict):
        return SequenceDatabase(db)
    
    elif isinstance(db, list):
        db_dict = []
        
        id = "id"
        for num, seq in enumerate(db, start = 1):
            db_dict.append({
                "id": id + str(num),
                "sequence": seq
            })

        return SequenceDatabase(db_dict)
    
    elif isinstance(db, str):
        return create_database(db)
    else:
        raise TypeError
