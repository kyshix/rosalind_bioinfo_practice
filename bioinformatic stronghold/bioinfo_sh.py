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
Reverse
"""

if __name__ == "__main__":
    print(count_dna("rosalind_dna.txt"))
    print(transcribe_to_mrna("rosalind_rna.txt"))
    print(reverse_complement_dna("rosalind_revc.txt"))