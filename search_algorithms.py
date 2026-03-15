from seq_process import to_kmers

def hashtable(seqs:list, kmers:list) :
    n = len(kmers[0])
    kmer_freqs = dict()

    for s in seqs :
        seq_kmers = to_kmers(s,n)
        for k in kmers :
            if k not in kmer_freqs : 
                kmer_freqs[k] = seq_kmers.count(k)
            else : kmer_freqs[k] = kmer_freqs[k] + seq_kmers.count(k)

    return dict(sorted(kmer_freqs.items(),key=lambda item: item[1], reverse=True))

def hammingDistance(a, b) :
    return sum([1 for i in range(len(a)) if a[i] != b[i]])

def totalDistance(kmer, sequences) :
    k = len(kmer)
    total_dist = 0

    for seq in sequences :
        dist = [hammingDistance(kmer,seq[i:i+k]) for i in range(len(seq)-k+1)]
        total_dist += min(dist)

    return total_dist

def medianString(seqs:list, kmers:list) :
    return dict(sorted({k:totalDistance(k,seqs) for k in kmers}.items(), key=lambda item: item[1], reverse=False))