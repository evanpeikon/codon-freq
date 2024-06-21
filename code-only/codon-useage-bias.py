DNA_sequence = 'ATCATCATGATGAGGCCCATGCATTAC' # input any dna seq as string 

def translate_codon(codon):    
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

def calculate_codon_usage_bias(DNA_sequence):
    codon_counts = {}
    for i in range(0, len(DNA_sequence), 3):
        codon = DNA_sequence[i:(i + 3)]
        amino_acid = translate_codon(codon)  
        if amino_acid not in codon_counts:
            codon_counts[amino_acid] = {}
        if codon in codon_counts[amino_acid]:
            codon_counts[amino_acid][codon] += 1
        else:
            codon_counts[amino_acid][codon] = 1
    codon_frequencies = {}
    for amino_acid, counts in codon_counts.items():
        total_counts = sum(counts.values())
        codon_frequencies[amino_acid] = {codon: count / total_counts for  codon, count in counts.items()}
    return codon_frequencies

codon_usage_bias = calculate_codon_usage_bias(DNA_sequence)
for amino_acid, frequencies in codon_usage_bias.items():
    print(f'Amino Acid: {amino_acid}')
    for codon, frequency in frequencies.items():
        print(f'  Codon: {codon}, Frequency: {frequency:.2f}')
