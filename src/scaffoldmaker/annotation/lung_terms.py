"""
Common resource for lung annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
lungs_terms = [
    ("Trachea", "FMA:00", "UBERON:0000000"),
    ("RMB", "FMA:01", "UBERON:0000001"),
    ("LMB", "FMA:02", "UBERON:0000002"),
    ("Left lung", "FMA:03", "UBERON:0000003")
]

def get_lungs_term(name : str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in lungs_terms:
        if name in term:
            return ( term[0], term[1] )
    raise NameError("Lungs annotation term '" + name + "' not found.")
