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
        print(f"{counts["A"]} {counts["C"]} {counts["G"]} {counts["T"]}")
                
        
if __name__ == "__main__":
    count_dna("rosalind.txt")