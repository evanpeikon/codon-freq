# 🧬 Finding Codon Frequency:

Codon frequency refers to the number of times a given codon, or nucleotide triple (i.e., ATC), occurs within a given DNA or RNA sequence, and it measures how often each codon appears relative to the total number of codons in the sequence.

A high codon frequency for a specific codon implies that it is used more often for coding a particular amino acid, while a low codon frequency indicates the opposite. In the code code block below, I'll show you how to calculate codon counts and codon frequency within a given DNA sequence:

```
DNA_sequence = 'ATCATCATGATGAGGCCCATGCATTAC' # input dna seq as a string of ATCG's
codon_counts = {} 

for i in range(0, len(DNA_sequence), 3): 
    codon = DNA_sequence[i:(i + 3)]
    if codon in codon_counts: 
        codon_counts[codon] += 1
    else: 
        codon_counts[codon] = 1 

total_codons = sum(codon_counts.values())
codon_frequencies = {}

for codon, codon_count in codon_counts.items(): 
    codon_frequencies[codon] = codon_count / total_codons 
print('Codon frequencies: %s' % codon_frequencies)
```

# 🧬 Codon Usage Bias

In many genomes you’ll observe codon usage bias, meaning that some synonymous codons, that code for the same amino acid, are used more frequently than others. For example, let’s take the amino acid leucine, which is encoded by the DNA codons TTA and TTG. If we have ten leucines in a protein sequence, we might expect to see each of the previously mentioned codons appear five times. However, some codons appear more or less frequently than expected, given the frequency of the amino acids they code for. This bias can result from natural selection and evolutionary forces. Additionally, some codons may be preferred due to translation efficiency, tRNA availability, or other factors.

In the code block below, I'll show you how to observe codon usage bias within a given DNA sequence:

```
DNA_sequence = 'ATCATCATGATGAGGCCCATGCATTAC'

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
```
