from collections import Counter

# COUNTING DNA NUCLEOTIDES
"""
Count Occurrences of DNA Nucleotide
:param fname: file containing DNA string s (length of at most 1000 nt)
:returns: four integers corresponding times nucleotides "A", "C", "G", "T" appear in s
"""
def count_dna(fname):
    with open(fname, "r") as s: 
        content = s.read()
        dna = list(content)     # turn string s into list of single nucleotides
        
        # METHOD #1
        # adenine_num = 0
        # cytosine_num = 0
        # guanine_num = 0
        # thymine_num = 0
        # for nt in dna:
        #     if nt == "A": 
        #         adenine_num += 1
        #     elif nt == "C":
        #         cytosine_num += 1
        #     elif nt == "G": 
        #         guanine_num += 1
        #     elif nt == "T": 
        #         thymine_num += 1
        # print(f"{adenine_num} {cytosine_num} {guanine_num} {thymine_num}")
        
        # METHOD #2
        counts = Counter(dna)
        return(f"{counts["A"]} {counts["C"]} {counts["G"]} {counts["T"]}")
                

# TRANSCRIBING DNA INTO RNA
"""
Transcribe DNA Sequence into Corresponding mRNA
:param fname: file containing a DNA string t with length of at most 1000 nt
:returns: transcribed mRNA string of t
"""
def transcribe_to_mrna(fname): 
    # convert all thymine nt to uracil nt 
    with open(fname, "r") as file:
        dna = file.read()
        mRNA = dna.replace("T", "U")
    return mRNA


# THE SECONDARY AND TETIARY STRUCTURE OF DNA
# secondary structure of dna refers to the double helix in which complementary strands
# run antiparallel of each other with base pariings A-T and G-C (storage of genetic info)
# tetiary structure of dna refers to the 3D arrangement of the double helix that allows
# you to understand folding/coiling patterns (access to genetic info)

# reverse complement counterpart represents the sequence of DNA strand (only one direction
# because of how dna synthesis occurs)
# this is essential for identifying potential binding sites, matching sequences, etc. 
"""
Reverse Complement Counterpart of DNA Strand
:param fname: file containing DNA string s of length at most 1000 bp (5'-3' strand)
:results: reverse complement strand to s (5'-3' strand)
"""
def reverse_complement_dna(fname): 
    with open(fname, "r") as file: 
        rev_dna = file.read()[::-1]
        
        # METHOD 1
        # rev_comp_dna = ""
        # for nt in rev_dna:
        #     if nt == "A": 
        #         rev_comp_dna += "T"
        #     elif nt == "C": 
        #         rev_comp_dna += "G"
        #     elif nt == "G": 
        #         rev_comp_dna += "C"
        #     elif nt == "T": 
        #         rev_comp_dna += "A"
        # return rev_comp_dna
        
        # METHOD 2
        bp_translation = str.maketrans("ACGT", "TGCA")
        rev_compl_dna = rev_dna.translate(bp_translation)
        return rev_compl_dna


# RABBIT AND RECURRENCE RELATIONS
# recurrance relation: refer to justin's recurrence relation notes
"""
Calculate Number of Rabbit Pairs Using Recurrence Relations
Using dynamic programming, calculate the fibonacci numbers as specified by its
recurrance relation (Fn=F_[n-1]+F_[n-2] with F_1 = F_2 = 1)
:param fname: file containing two positive integers
    :n: <=40 represents the number of months/generations of wascally rabbits
    :k: <=5 represents the litter size (pair) of each preproducing rabbit pair
:returns: total number of rabbit pairs after n months
"""
# assuming that rabbits live forever -> base model to analyze growth
def wascally_wabbits(fname):
    with open(fname, "r") as file: 
        data = file.read().split()
        
    months = int(data[0])       # n
    litter_size = int(data[1])  # k
    return wabbit_reccurence(months, litter_size)

# it takes 2 months for the wabbits to mature before they reproduce
def wabbit_reccurence(n, k):
    # F_1 = F_2 = 1
    if n == 1 or n == 2: 
        return 1
    elif n <= 0: 
        return 0
    else: 
        # Fn=F_[n-1]+F_[n-2] 
        # given a month, # of mature wabbit pairs = total # of wabbit pairs the month prior
        # so Fn=F_[n-1]+k(F_[n-2]) is the revised equation to account for diff litter sizes 
        return wabbit_reccurence(n-1, k) + k*wabbit_reccurence(n-2, k)
    

# COMPUTING GC CONTENT
"""
GC Content of DNA
Determine the Max GC Conent DNA Sequence 
GC content can be used to determine the identity of the organism the DNA is from and
gives insight to related specimens
:param fname: file containing the dna in FASTA format
:returns: maximum GC DNA sequence -> formatted as
    - DNA sequence id without ">"
    - GC content percentage (0.0001 absolute error)
"""
def gc_content(fname): 
    dna_dict = {}
    
    # incorrect bc it assumes that dna in fasta format
    # with open(fname, "r") as file:
    #     for id, dna_str in zip(file, file): 
    #         id = id.strip().replace(">", "")
    #         dna = dna_str.strip().upper()
    #         # calculate gc content
    #         gc_count = dna.count("C") + dna.count("G")
    #         gc_percent = (gc_count/len(dna)) * 100
    #         dna_dict[id] = gc_percent
    
    # corrected: takes all lines following the dna label till the next label
    with open(fname, "r") as file: 
        fasta_id = None
        dna = ""
        
        for line in file:
            if line.startswith(">"): 
                # if there is a dna sequence already being read for
                if fasta_id is not None:
                    gc_count = dna.count("C") + dna.count("G")
                    gc_percent = (gc_count/len(dna)) * 100
                    dna_dict[fasta_id] = gc_percent
                # get the id without ">"
                fasta_id = line.strip()[1:]
                dna = ""
            else:
                dna = dna + line.strip().upper()
        
        # save the last dna sequence data
        gc_count = dna.count("C") + dna.count("G")
        gc_percent = (gc_count/len(dna)) * 100
        dna_dict[fasta_id] = gc_percent
    
    # determine the key value pair of the max gc contnet  
    max_id = max(dna_dict, key=dna_dict.get)
    max_gc = dna_dict[max_id]
    print(f"{max_id}\n{max_gc:.6f}")

# COUNTING POINT MUTATIONS
"""
Count Point Mutations
Point mutations is a consequence of evolution (mitosis, meiosis, or generally any form
of DNA replication that may be caused by external disturbances) in which one a base is
replaced with another (transversion and transitions only in this case?? based on rosalind)
:param fname: contains two DNA sequences 
:returns: number of point mutations that occurred
"""
def pt_mutations(fname): 
    with open(fname, "r") as file: 
        dna_seqs = file.read().splitlines()
    # assume that the two sequences are of the same length
    return sum(1 for a, b in zip(dna_seqs[0], dna_seqs[1]) if a != b)

if __name__ == "__main__":
    print(count_dna("rosalind_dna.txt"))
    print(transcribe_to_mrna("rosalind_rna.txt"))
    print(reverse_complement_dna("rosalind_revc.txt"))
    print(wascally_wabbits("rosalind_fib.txt"))
    gc_content("rosalind_gc.txt")
    print(pt_mutations("rosalind_hamm.txt"))