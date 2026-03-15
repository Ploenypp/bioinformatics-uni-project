import pandas as pd

def get_seqs_from_fasta(file_name) :
    f = open(file_name, "r")
    file = f.readlines()
    f.close()

    join_all = "".join(file)
    by_sample = join_all.split(">")[1:]

    seqs_dico = dict()
    for s in by_sample :
        info_dico = dict()
        info = s.split("\n")
        
        for label, val in [x.split(":") for x in info[0].split("\n")[0].split(" ")[1:]] :
            info_dico[label] = val
        info_dico["seq"] = "".join(info[1:-1]).upper()
        
        seqs_dico[info[0].split(" ")[0]] = info_dico

    return seqs_dico

def seqs_info_to_pd(seqs) :
    seq_id = list(seqs.keys())
    chromosomes = []
    starts = []
    stops = []
    sequences = []

    for s in seqs.values() :
        chromosomes.append(s["Chromosome"])
        starts.append(s["Start"])
        stops.append(s["Stop"])
        sequences.append(s["seq"])
    
    data = {"sequence id":seq_id,"chromosome":chromosomes,"start":starts,"stop":stops, "sequence":sequences}
    return pd.DataFrame(data)

def id_seq(seqs) :
    return {i:info["seq"] for i, info in seqs.items()}

def seq_only(seqs) :
    return [info["seq"] for _, info in seqs.items()]