from database.sequence_database import SequenceDatabase
from database.load_database import create_database

def normalize_database(DB):

    if isinstance(DB, SequenceDatabase):
        return DB
    
    elif isinstance(DB, list):
        return SequenceDatabase(DB)
    elif isinstance(DB, str):
        return create_database(DB)
    else:
        raise TypeError
