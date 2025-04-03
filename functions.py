from math import floor
from textwrap import wrap
from collections import Counter
from suffix_trees import STree

import numpy as np
import re

nuc = ('A', 'C', 'G', 'T')
#  READING FILES
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

def readJaspar(jasparFileName) :
    file = open(jasparFileName)
    lines = [re.sub(r'\s+', ' ', f.strip()) for f in file]

    mat = []

    i = 1
    while i < len(lines) :
        line = lines[i]

        lst = []
        aux = ""

        for x in line :
            try :
                n = int(x)
            except ValueError :
                if len(aux) > 0 : lst.append(int(aux))
                aux = ""
            else : 
                aux += x

        mat.append(lst)
        i += 1
    return mat

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
    
    cpt = 0
    limit = len(str1)
    if len(str1) >= len(str2) : limit = len(str2)
    for i in range(limit) :
        if str1[i] != str2[i] : cpt += 1
    return cpt

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

# LOCATING THE FIXATION SITE 

# Indexing motifs of a fixed length
def indexTable(m, sequence):
    """
    Indexer les positions d'occurrences de tous les mots de taille k dans une sequence
    entrée m : taille du mot à chercher dans le motif m <= k
    entrée sequence : chaine de caractère représentant une sequence d'ADN
    sortie indexes : dictionaire où les clés sont les mots et les valeurs les positions dans la sequence
    """
    indexes  = {}
    for i in range(len(sequence)) :
        if i + m > len(sequence) : break
        
        aux = sequence[i:i+m] 
        if aux in indexes : indexes[aux].append(i)
        else : indexes[aux] = [i]

    return indexes

def chercherWithIndexTable(m, table, sequence, motif, maxVar):
    """
    chercher les positions d'un motif dans une séquence en admettant au maximum maxVar variations
    entrée m : taille du mot à chercher dans le motif m <= k, le meme utilise pour indexer les sequences
    entrée table : dictionaire où les clés sont les mots et les valeurs les positions dans la sequence
    entrée sequence : chaine de caractère représentant une sequence d'ADN
    entrée motif : chaine de caractère représentant le motif à chercher
    entrée maxVar : le maximum variations entre le motif et un mot de taille k dans la sequence
    sortie motifPos : dictionnaire où les clés sont les motifs trouvé et les valeurs leurs positions dans la sequence.
    """
    k = len(motif)
    motif = motif.lower()
    motifPos = dict()

    for cand,pos in table.items() :
        if cand.lower() in motif :
            aux_dict = dict()
            for p in pos :
                tmp = sequence[p:p+k].lower()
                if hamDistance(motif,tmp) < maxVar : 
                    if tmp in aux_dict : aux_dict[tmp].append(p)
                    else : aux_dict[tmp] = [p]
            if len(aux_dict) > 0 :
                for x,lst in aux_dict.items() :
                    motifPos[x] = lst

    return motifPos

def findMotifData(sequences, motif, m, maxVar):
    """
    chercher les positions d'un motif dans un ensemble de séquence d'ADN en admettant un maximum de variations
    entrée m : taille du mot à chercher dans le motif m <= k, le meme utilise pour indexer les sequences
    entrée sequences : list contenant les sequence d'ADN
    entrée motif : chaine de caractère représentant le motif à chercher
    entrée maxVar : le maximum variations entre le motif et un mot de taille k dans la sequence
    sortie posList : list contenant les positions dans les sequences ou se trouve le motif.
    """
    posDict_seq = dict()
    i = 0
    while i < len(sequences) :
        idx = indexTable(m, sequences[i])
        posDict = chercherWithIndexTable(m, idx, sequences[i], motif, maxVar)
        if len(posDict) > 0 : posDict_seq[i] = posDict

        i += 1

    return posDict_seq

# Frequency Matrixes
def computing_pwm(M, cols):
    """
    Calcul la matrice de poids position à partir de la matrice de frequence
    entrée M : matrice de frequence
    sortie PWM : matrice de probabilites ou poids position
    """
    PWM = M+1
    sum_col = np.sum(PWM,axis=0)

    return PWM/sum_col

def f0_calcule(PWM, L):
    """
    Calcul les valeurs de probabilites d'un modele independant de positions (modele Null)
    entrée PWM : matrice de probabilites ou poids positions
    sortie  f_0 : vecteur contenant un modele independant de positions (modele Null)
    """
    return np.sum(PWM,axis=1)/L

def loglikehood(seq, PWM, f_0, L):
	"""
	Calcul le rapport de vraissemblance entre une sequence et une matrice de poids position
	entrée PWM : matrice de probabilites ou poids positions
	entrée f_0 : vecteur contenant le modele independant de positions (modele Null)
	entrée seq : une sequence d'ADN de taille k, ou k est le nombre de colonnes de PWM
	sortie ll : rapport de vraissemblance
	"""
	L_0 = len(seq)
	l_win = 0
	
	for i in range(L) :
		x = 3
		aux = seq[i].lower()
		if aux == 'a' : x = 0
		elif aux == 'c' : x = 1
		elif aux == 'g' : x = 2
	
		l_win += np.log2(PWM[x,i]/f_0[x])

	return l_win

def searchPWMOptmiseMotifs(sequences, k, PWM, f_0):
    """
    Cherche les positions dans un ensemble de séquence qui maxime le rapport de vraisemblance et elimine les motifs chevauchante
    entrée sequences : ensemble de séquence d'ADN
    entrée k : nombre de colonnes d'PWM
    entrée PWM : matrice de probabilités ou poids positions
    entrée f_0 : vecteur contenant le modèle indépendant de positions (modèle Null)
    sortie posList: liste contenant pour chaque séquence la/les positions ayant un rapport de vraisemblance positive

    """

    posList = []

    all_pos = dict()


    for i in range(len(sequences)) :
        seq = sequences[i]
        limit = len(seq)-k+1

        pos = []
        prev_pos, prev_log = -1,0

        for j in range(limit) :
            curr_log = loglikehood(seq[j:j+k],PWM,f_0,k)

            if curr_log > 0 :
                if j <= prev_pos+k and curr_log > prev_log : 
                    if prev_pos in pos : pos.remove(prev_pos)
                else : 
                    pos.append(j)
                prev_pos = j
                prev_log = curr_log

        if len(pos) > 0 :
            all_pos[i] = pos
            posList += pos
    
    return all_pos