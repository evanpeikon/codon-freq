# define dna sequence as string of nucleotides (i.e., A’s, T’s, G’s, and C’s). 
DNA_sequence = 'ATCATCATGATGAGGCCCATGCATTAC'

# create dictionary named codon_counts to store counts of ea/ unique codon in DNA_sequence
codon_counts = {}

# create for loop to iterate over DNA_sequence in steps of 3 
for i in range(0, len(DNA_sequence), 3):
    codon = DNA_sequence[i:(i + 3)] # extract ea/ 3 nucleotide long codon from DNA_sequence
    if codon in codon_counts: # If the extracted codon is already in the codon_counts dictionary....
        codon_counts[codon] += 1 # then, the code increments the count by 1. 
    else: # If the extracted codon isn't already in the codon_counts dictionary...
        codon_counts[codon] = 1 # initialize the count to 1. 

total_codons = sum(codon_counts.values()) # sum all values in the codon_counts dictionary and assign value to total_codons
codon_frequencies = {} # create empty dictionary to store frequencies of ea/ unique codon in dna sequence 

for codon, codon_count in codon_counts.items(): # for loop iterates over the codons and their counts stored in the codon_counts dictionary
    codon_frequencies[codon] = codon_count / total_codons # for ea/ codon, calculate the frequence by dividing its count by total_codons, then store free in codon_frequencies dict. 

print('Codon frequencies: %s' % codon_frequencies) # print frequency of each codon in codon_frequencies dict
