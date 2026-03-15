from itertools import product
from collections import Counter

def generate_kmers(length) :
    nt = ["A","C","G","T"]
    aux = product(nt, repeat=length)
    return ["".join(x) for x in aux]

def to_kmers(seq:str, length) :
    return [seq[i:i+length] for i in range(len(seq)- length + 1)]

def kmers_freq(seq, length) :
    kmers = to_kmers(seq, length)
    return dict(Counter(kmers))

def remove_homo(kmers, limit) :
    toRem = set()

    for k in kmers :
        prev = ""; cpt = 0
        for x in k :
            if x == prev : cpt += 1
            if cpt >= limit : 
                toRem.add(k)
                break
            prev = x
            
    return set(kmers)-toRem

def remove_hetero(kmers, limit) :
    toRem = set()

    nt = ["A","C","G","T"]    
    dinuc = {a + b for a,b in set(product(nt, repeat=2)) if a != b}

    for k in kmers : 
        aux = dict()
        for d in dinuc :
            aux[d] = int(k.count(d))

        for _,freq in aux.items() :
            if freq >= limit : toRem.add(k)
    return set(kmers)-toRem

# to eliminate motifs of which A and T make up more than a certain proportion
def remove_ta(motifs, limit) : # limit is a proportion (0 < limit < 1)
    result = set()
    
    for motif in motifs :
        cpt = 0
        a,t = False, False
        for x in motif :
            if x == 'A' :
                cpt += 1
                a = True
            elif x == 'T' :
                cpt += 1
                t = True
        if not (cpt/len(motif) > limit) and (a and t) : result.add(motif)
    return result

def candidate_kmers(k, homo, hetero, ta) :
    return list(remove_ta(remove_hetero(remove_homo(generate_kmers(k), homo), hetero), ta))