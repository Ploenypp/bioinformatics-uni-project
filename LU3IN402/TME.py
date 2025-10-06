from collections import Counter 
import re
nuc = ('A', 'C', 'G', 'T')


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

def searchGivenMotif(motifsTrouve, motifSpecifique, decroissant = True):
    """
    Cherche un motif specifique dans un dictionnaire de motifs trouvés
    entrée motifsTrouve : dictionnaire de motifs, clé = motif, valeur = fréquence d'observation
    entrée motifSpecifique: un motif specifique à chercher
    entrée decroissant : bool, si True, le dictionnaire est trié par ordre décroissant de valeur
    sortie fréquence : la fréquence du motif
    sortie ranking : dans quelle position le motif a été trouvé
    >>>searchGivenMotif(test_motifs, "TAT")
    (2, 2)
    """
    
    ranking = 1
    frequence = 0
    
    ordered = sorted(motifsTrouve.items(),key=lambda item : item[1])
    if decroissant : ordered = reversed(ordered)
    ordered = dict(ordered)

    for (key,freq) in ordered.items() :
        if key == motifSpecifique :
            frequence = freq
            break
        ranking += 1
        
    return ranking, frequence

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

def gatherRevCompMotifs(motifs_dict):
    """
    Trouver les motifs reverse complementaire et rassembler leur résultats
    entrée motifs_dict : dictionaire de résultats de recherche de motifs
    sortie motifs_dict_gather : dictionaire de motifs, clé = motif, valeur = valeur motif fw, valeur motif rv
    """
    
    motifs_dict_gather = {}
    for motif in motifs_dict :
        rev = reversecompl(motif)
        if rev not in motifs_dict_gather.keys() :
            rVal = 0
            if rev in motifs_dict :
                rVal = motifs_dict[rev]
            if motifs_dict[motif] > rVal :
                motifs_dict_gather[motif] = (motifs_dict[motif],rVal)
            else : motifs_dict_gather[motif] = (rVal,motifs_dict[motif])
    
    return motifs_dict_gather

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