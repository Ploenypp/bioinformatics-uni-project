from math import floor
from textwrap import wrap
from collections import Counter
from suffix_trees import STree

import re

nuc = ('A', 'C', 'G', 'T')

def readFasta(fastaFileName):
    """
    Lit un fichier fasta
    entrée fastaFileName: nom du fichier fasta
    sortie séquences: liste contenant toutes les séquences du fichier
    """
    
    sequence = ""
    sequences_list = []
    prev_header = ""
    header = ""
        
    for line in open(fastaFileName):
        string = line.strip()
        if string[0] != ">":
            if prev_header != header:
                prev_header = header
            sequence = sequence + string
        else:
            header = string
            if sequence != "":
                sequences_list.append(sequence)
                sequence = ""

    sequences_list.append(sequence)
    return sequences_list

# PREPARING CANDIDATE KMERS
def searchMotifs(k:int, sequences:list):
    """
    Cherche les motifs de taille k dans un ensemble de séquences
    entrée k : taille du motif
    entrée séquences : liste de séquences
    sortie motifs: dictionnaire de motifs, clé = motif, valeur = fréquence d'observation
    >>>searchMotifs(3, ['TAAGTAA', 'TATAA', 'CTATC'])
    {'TAA': 3, 'AAG': 1, 'AGT': 1, 'GTA': 1, 'TAT': 2, 'ATA': 1, 'CTA': 1, 'ATC': 1}
    """
    
    motifs  = {}

    for seq in sequences :
        for i in range(len(seq)-k+1) :
            aux = seq[i:i+k]
            if aux in motifs : motifs[aux] += 1
            else : motifs[aux] = 1
    return motifs

# eliminating uninformative kmers
def removeLowComplexeHomo(motifs:list, m:int):
    """
    Enlève les motifs peu complexe ayant m fois le même nucléotide
    entrée motifs: liste de motifs
    entrée m: taille de repetition de nucléotide
    sortie motifsClean: liste de motifs sans les motifs peu complexe
    >>>removeLowComplexeHomo(['TAA', 'AAG', 'AGT', 'GTA', 'TAT', 'ATA', 'CTA', 'ATC'], 2)
    ['AGT', 'GTA', 'CTA', 'ATC']
    """

    motifsNotClean = []

    for motif in motifs :
        freq = Counter(motif)
        for x,f in freq.items() :
            if f >= m : motifsNotClean.append(motif)
    return list(set(motifs)-set(motifsNotClean))

def removeLowComplexeHetero(motifs:list, n:int, variation = "yes"):
    """
    Enlève les motifs peu complexe ayant n fois un dinucléotide
    entrée motifs: liste de motifs, clé = motif, valeur = fréquence d'observation
    entrée n: taille de répétition de dinucléotide
    entrée variation : string, si "yes" permettre variation d'un nucléotide
    sortie motifsClean: liste de motifs sans les motifs peu complexe
    >>>removeLowComplexeHetero(['GGTTTGG', 'TGAGTTA', 'TGCCGTG', 'AGAGAGA', 'TCACCGA', 'TTGGTAT', 'AGGGTGG', 'TGGCTTA', 'AGAGTAG', 'GCCCCTC'], 3, variation = "yes")
    ['GGTTTGG', 'TGAGTTA', 'TGCCGTG', 'AGAGAGA', 'TCACCGA', 'TTGGTAT', 'TGGCTTA', 'AGAGTAG']
    """
    
    motifsClean = []
    trouve = []

    # Créer une liste de dinucléotides sans mono dinucléotides
    all_di_nuc = [nuc_1 + nuc_2 for nuc_1 in nuc for nuc_2 in nuc if nuc_1 != nuc_2]

    # Pour chaque motif de la liste
    for motif in motifs:
        if variation == "yes":
            # Permettre une varation
            trouve = [nuc_di for nuc_di in all_di_nuc if len(re.findall(re.compile("(?=("+ nuc_di +"))"), motif)) >= n]

        else:
            # Pas permettre de varations
            trouve = [nuc_di for nuc_di in all_di_nuc if len(re.findall(re.compile("(?=("+ nuc_di +")"+ nuc_di[0] +")"), motif)) >= n]

        # Après avoir evalué tous les dinucléotides, évaluer si on a trouvé de dinucleotides
        if len(trouve) == 0:
            # Ajouter juste les motifs que n'ont pas eu des dinucléotides trouves
            motifsClean.append(motif)
        # Vider la liste
        trouve = []

    return motifsClean

def removeTARich(motifs:list, p:int):
    """
    Enlève les motifs contenant > p de T et A en combination
    entrée motifs: liste de motifs
    entrée p: proportion de nucleótides T et A
    sortie motifsClean: liste de motifs sans les motifs TA rich
    """

    motifsClean = []

    for m in motifs :
        cpt = 0
        a,t = False,False
        for x in m :
            if x == 'A' or x == 'a' :
                cpt += 1
                a = True
            if x == 'T' or x == 't' :
                cpt += 1
                t = True
        if not (cpt/len(m) > p) and (a and t) : motifsClean.append(m)

    return motifsClean

# SORTING & RANKING CANDIDATES
def getTopMotifs(motifs:dict, top, decroissant = True):
    """
    Renvoie les top motifs le plus fréquent
    entrée motifs: dictionnaire de motifs, clé = motif, valeur = fréquence d'observation
    entrée top : les top plus fréquent, si top == "all", ne filtrez pas et renvoyer le dictionnaire complet trié
    entrée decroissant : bool, si True le dictionnaire est trié par ordre décroissant de valeur
    sortie motifsfreq: dictionnaire contenant les top motifs les plus fréquents, clé = motif, valeur = fréquence d'observation
    >>>getTopMotifs({'TAA': 3, 'AAG': 1, 'AGT': 1, 'GTA': 1, 'TAT': 2, 'ATA': 1, 'CTA': 1, 'ATC': 1}, 2)
    {'TAA': 3, 'TAT': 2}
    """

    lst = [(key,freq) for key,freq in motifs.items()]
    lst.sort(key=lambda x : x[1], reverse=decroissant)

    if top == "all" : return dict(lst)
    return dict(lst[:top])

# HASHTABLE
def hashTable(sequences:list, kmersV:list):
    """
    Cherche les motifs de taille k dans un ensemble de séquences avec la méthode de Hash Table
    entrée séquences : liste de séquences
    entrée kmersV: liste de Kmers valides à chercher
    sortie motifs_freq_dict: dictionnaire de motifs, clé = motif, valeur = fréquence d'observation
    >>>hashTable(['TAAGTAA', 'TATAA', 'CTATC'], ['TAA', 'AAG', 'AGT', 'GTA', 'TAT', 'ATA', 'CTA', 'ATC'])
    {'TAA': 3, 'AAG': 1, 'AGT': 1, 'GTA': 1, 'TAT': 2, 'ATA': 1, 'CTA': 1, 'ATC': 1}
    """
    
    motifs_freq_dict  = {}
    k = len(kmersV[0])
    for seq in sequences :
        for i in range(len(seq)-k+1) :
            aux = seq[i:i+k]
            if aux in motifs_freq_dict : motifs_freq_dict[aux] += 1
            else : motifs_freq_dict[aux] = 1
    
    return {key:freq for key,freq in motifs_freq_dict.items() if key in kmersV}

# MEDIAN STRING
def hamDistance(str1:str, str2:str):
    """
    Calcul la distance de Hamming entre deux chaînes de caractères
    entrée str1: chaîne de caractères
    entrée str2: chaîne de caractères
    sortie distance: distance de Hamming
    >>>hamDistance("TTGGTAT", "TTGCTAA")
    2
    """
    
    distance = 0
    for i in range(len(str1)) :
        if str1[i] != str2 [i] : distance += 1
        
    return distance

def totalDistance(motif:str, sequences, k):
    """
    Calcul la totalDistance
    entrée motif: motif à comparer, chaîne de caractères
    entrée sequences: liste de séquences
    entrée k: taille du motif
    sortie total_distance: somme de distance de hamming minimal
    """
    
    total_distance = 0
    for seq in sequences :
        dist = []
        for i in range(len(seq)-k+1) :
            aux = seq[i:i+k].upper()
            dist.append(hamDistance(motif,aux))
        total_distance += min(dist)

    return total_distance

def medianStringSearch(sequences, kmersV):
    """
    Implement l'algorithme MedianStringSearch
    entrée séquences : liste de séquences
    entrée kmersV: Liste de Kmers à chercher
    sortie motif_dist_dict: un dictionnaire contenant les motifs et leurs distances
    """
    
    motif_dist_dict = {k:totalDistance(k,sequences,len(kmersV[0])) for k in kmersV}

    return getTopMotifs(motif_dist_dict,len(motif_dist_dict),False)

# SUFFIX TREES
def constructTree(sequences):
    """
    construis un abre de suffixes
    entrée sequences : liste de séquences
    sortie suffix_tree : arbre de suffixes
    """    

    suffix_tree = STree.STree(sequences)

    return suffix_tree

def getSeeds(kmersV, v):
    """
    entrée kmersV : liste de motifs à chercher
    entrée v : nombre de variations dans les motifs
    sortie k : taille de k-mer
    sortie seed_size : taille de seed calculé
    sortie seed_nb : nombre de seeds
    """
    
    k = len(kmersV[0])
    
    seed_nb = v+1
    
    if (k%seed_nb != 0) :
        seed_size = floor(k/seed_nb) + 1
    else :
        seed_size = floor(k/seed_nb)
    
    return k, seed_size, seed_nb

def findSeedStarPos(kmer, seed_size, seed_nb, stree):
    """
    Renvoie tous les positions de départ des seeds de k-mer
    entrée kmer : chaîne de caractères avec le motif à analyser
    entrée seed_size : taille de seed calculé
    entrée seed_nb : nombre de seeds contenus dans le k-mer
    entrée stree : sufix tree des séquénces à analyser
    sortie all_candidate_starts : liste de positions de départ qui donnent un match partiel
    """
    
    if len(kmer) < seed_size : return []

    seeds = {s for s in wrap(kmer,seed_size) if len(s) >= seed_size}
    all_candidates = {x for s in seeds for x in stree.find_all(s) if x >= 0}
    candidate_starts = sorted(list({x for x in all_candidates if x >= 0}))
        
    return candidate_starts

def findKmerCandidates(all_candidate_starts, k, kmer, sequences, v):
    """
    entrée all_candidate_starts : liste de positions de départ qui donnent un match partiel
    entrée k : taille de k-mer
    enréee kmer : chaîne de caractères avec le motif à analyser
    entrée sequences : liste de séquences
    entrée v : nombre de variations dans les motifs
    sortie : liste de séquence qui donnent un match avec un nombre inferieur ou égal de variations (v)
    """
    
    all_candidate_kmer = set()
    for c in all_candidate_starts :
        if hamDistance(kmer,sequences[c:c+k]) <= v : 
            all_candidate_kmer.add(sequences[c:c+k])
    return all_candidate_kmer

def inexactMatch(kmersV, sequences, stree, v):
    """
    cherche de motifs variables dans un suffix tree
    entrée kmersV : liste de motifs à chercher
    entrée sequences : liste de séquences
    entrée stree : suffix tree
    entrée v : nombre de variations dans les motifs
    sortie motif_occur : dictionnaire clés = motif ; value = nombre d'occurrences
    sortie motif_seq : dictionnaire clés = motif ; value = liste de motifs variables
    """

    motif_occur = dict()
    motifs_seq = dict()

    motifs_seq = {}
    for motif in kmersV :
        k, seed_size, seed_nb = getSeeds([motif],v)
        pos = findSeedStarPos(motif,seed_size,seed_nb,stree)
        candidates = findKmerCandidates(pos,k,motif,sequences,v)

        for aux in candidates :
            motif_occur[aux] = len(stree.find_all(aux))
        
        motifs_seq[motif] = candidates
    
    return motif_occur, motifs_seq