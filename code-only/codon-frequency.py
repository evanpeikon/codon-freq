DNA_sequence = 'ATCATCATGATGAGGCCCATGCATTAC' #input any dna sequence as string
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
