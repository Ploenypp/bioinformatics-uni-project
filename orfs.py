def to_codon_list(seq:str) :
    return [seq[i:i+3] for i in range(0, len(seq)-2, 3)]

def complementary(seq:str) :
    compl = {"A":"T", "C":"G", "G":"C","T":"A"}
    return "".join([compl[x] for x in seq])

def reverse_complementary(seq:str) :
    return complementary(seq)[::-1]

def start_stop_pos(codons) :
    starts = [idx for idx, codon in enumerate(codons) if codon == "ATG"]
    stops = [idx for idx, codon in enumerate(codons) if codon in ["TAG","TAA","TGA"]]
    return starts, stops

def frames_strands(seqs_with_id) :
    res_dico = dict()
    for idx, seq in seqs_with_id.items() :
        aux_dico = dict()
        for i in range(1,4) :
            aux_dico[f"frame {i}"] = {"sense":seq,"antisense":reverse_complementary(seq)}
        res_dico[idx] = aux_dico
    return res_dico

def get_orfs_coor(seq) :
    starts, stops = start_stop_pos(to_codon_list(seq))
    coor = []; prev_start = -1; prev_stop = -1
    
    for i in range(len(starts)) :
        if starts[i] > prev_start and starts[i] > prev_stop :
            prev_start = starts[i]
            for j in range(len(stops)) :
                if stops[j] > starts[i] :
                    coor.append((starts[i], stops[j]))
                    prev_stop = stops[j]
                    break
    return coor

def get_all_orfs(all_frames_strands) :
    res_dico = dict()
    for idx, frames in all_frames_strands.items() :
        aux_dico = dict()
        for f, strand in frames.items() :
            aux_dico[f] = {"sense":get_orfs_coor(strand["sense"]), "antisense":get_orfs_coor(strand["antisense"])}
        res_dico[idx] = aux_dico
    return res_dico

def remove_null_orfs(all_orfs) :
    res_dico = dict()
    for idx, frames in all_orfs.items() :
        aux_dico = dict()
        for f, orfs in frames.items() :
            if len(orfs["sense"]) > 0 or len(orfs["sense"]) : aux_dico[f] = orfs
        if len(aux_dico) > 0 : res_dico[idx] = aux_dico
    return res_dico

def get_seq_from_coor(seq, start, stop) :
    return seq[start*3:stop*3+3]