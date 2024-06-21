DNA_sequence = 'ATCATCATGATGAGGCCCATGCATTAC' # define dna seq as string of nucleotides (i.e., A’s, T’s, G’s, and C’s)

def translate_codon(codon): # create function that takes codon as input and returns corresponding amino acid 
  codon_table = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'}

def calculate_codon_usage_bias(DNA_sequence): # define function that takes DNA_sequence as input and returns dictionary containing codon useage bias
    codon_counts = {} # create empty dictionary to store counts of ea/ codon for ea/ amino acid 
  
    for i in range(0, len(DNA_sequence), 3): # iterate through dna sequence in steps of 3 nucleotides. 
        codon = DNA_sequence[i:(i + 3)] # extract 3 nucleotide codon from DNA_sequence
        amino_acid = translate_codon(codon) # translate codon function finds corresponding amino acid for ea/ codon. 
        if amino_acid not in codon_counts: # if amino acid NOT already in codon_counts dictionary (ie, not a real codon)...
            codon_counts[amino_acid] = {} # doesn't do anything...
        if codon in codon_counts[amino_acid]: # if amino acid is in codon_counts dict already...
            codon_counts[amino_acid][codon] += 1 # increment the count by 1 
        else: # if amino acid is real, but not in codon_counts dict already... 
            codon_counts[amino_acid][codon] = 1 # create new entry w/ count of 1

# after iterating through DNA_sequence, the codon_counts dict contains count of ea/ codon for ea/ amino acid. 
  
    codon_frequencies = {} # create dict to store freq of ea/ codon for ea/ amino acid 
  
    for amino_acid, counts in codon_counts.items(): # for loop calculates codon freq for ea/ amino acid
        total_counts = sum(counts.values()) # foe ea/ amino acid, calculate total counts of all codons 
        codon_frequencies[amino_acid] = {codon: count / total_counts for  codon, count in counts.items()} # compute freq of ea/ codon by dividing count by the total counts 
    return codon_frequencies # return dictionary 


codon_usage_bias = calculate_codon_usage_bias(DNA_sequence) # call function w/ DNA_sequence as argument 

for amino_acid, frequencies in codon_usage_bias.items(): # outer for loop iterates through items in codon_usage_bias dict using the variables amino_acid and frequencies to represent each amino acid and its associated codon frequency.
    print(f'Amino Acid: {amino_acid}') # print amino acid name
    for codon, frequency in frequencies.items(): # inner for loop iterates through items in freq dict containing codon-freq pair for current amino acid. 
        print(f'  Codon: {codon}, Frequency: {frequency:.2f}') # print codon name and associated freq
