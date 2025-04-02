def reversecompl(seq:str):
    """Renvoie le brin complémentaire d'une séquence.
    entrée seq : sequence de nucléotides (brin sens)
    sortie     : sequence de nucléotides (brin complementaire)
    >>> reversecompl('AACGTGGCA')
    'TGCCACGTT'
    """
    compl = {'A': 'T', 'C': 'G', 'G': 'C', 'T':'A'}
    
    rev = list(seq)
    rev.reverse()
    
    res = ""
    for x in rev : res += compl[x]
    return res