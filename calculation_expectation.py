import numpy as np

#    DNA storage channel codec
#    DNA sequence constraints
#        Homopolymer = 3


# Probability calculation
def prob_calculater(present_dna, last_dna):
    '''
    Select the last three characters in the DNA sequence and perform an inductive calculation
        DNA sequence combination of sequence length 3
        AAA meanings: The last 3 symbols are the DNA sequence of AAA
        Array map:
            0    'AAA': 0, 'AAC': 1, 'AAG': 2, 'AAT': 3, 'ACA': 4, 'ACC': 5, 'ACG': 6, 'ACT': 7,
            1    'AGA': 0, 'AGC': 1, 'AGG': 2, 'AGT': 3, 'ATA': 4, 'ATC': 5, 'ATG': 6, 'ATT': 7,
            2    'CAA': 0, 'CAC': 1, 'CAG': 2, 'CAT': 3, 'CCA': 4, 'CCC': 5, 'CCG': 6, 'CCT': 7,
            3    'CGA': 0, 'CGC': 1, 'CGG': 2, 'CGT': 3, 'CTA': 4, 'CTC': 5, 'CTG': 6, 'CTT': 7,
            4    'GAA': 0, 'GAC': 1, 'GAG': 2, 'GAT': 3, 'GCA': 4, 'GCC': 5, 'GCG': 6, 'GCT': 7,
            5    'GGA': 0, 'GGC': 1, 'GGG': 2, 'GGT': 3, 'GTA': 4, 'GTC': 5, 'GTG': 6, 'GTT': 7,
            6    'TAA': 0, 'TAC': 1, 'TAG': 2, 'TAT': 3, 'TCA': 4, 'TCC': 5, 'TCG': 6, 'TCT': 7,
            7    'TGA': 0, 'TGC': 1, 'TGG': 2, 'TGT': 3, 'TTA': 4, 'TTC': 5, 'TTG': 6, 'TTT': 7}
    '''

    # Get DNA sequence length
    m = int(present_dna.shape[0] / 2)

    # For the convenience of calculation, add a layer of zero to last_dna
    zero_ = np.zeros([1, present_dna.shape[1], present_dna.shape[2]])
    last_dna = np.insert(last_dna, last_dna.shape[0], values=zero_, axis=0)

    # Probability of getting DNA symbol based on the index
    for i in range(m, present_dna.shape[0]):
        present_dna[i, 0, 0] = 0.25 * (last_dna[i - 2, 2, 0] + last_dna[i - 2, 4, 0] + last_dna[i - 2, 6, 0])
        present_dna[i, 0, 1] = 0.25 * (last_dna[i - 2, 2, 0] + last_dna[i - 2, 4, 0] + last_dna[i - 2, 6, 0])\
                               + 0.5 * last_dna[i - 1, 0, 0]
        present_dna[i, 0, 2] = 0.25 * (last_dna[i - 2, 2, 0] + last_dna[i - 2, 4, 0] + last_dna[i - 2, 6, 0])\
                               + 0.5 * last_dna[i - 1, 0, 0]
        present_dna[i, 0, 3] = 0.25 * (last_dna[i - 2, 2, 0] + last_dna[i - 2, 4, 0] + last_dna[i - 2, 6, 0])
        present_dna[i, 0, 4] = 0.25 * (last_dna[i - 2, 0, 1] + last_dna[i - 2, 2, 1] + last_dna[i - 2, 4, 1]
                                       + last_dna[i - 2, 6, 1])
        present_dna[i, 0, 5] = 0.25 * (last_dna[i - 2, 0, 1] + last_dna[i - 2, 2, 1] + last_dna[i - 2, 4, 1]
                                       + last_dna[i - 2, 6, 1])
        present_dna[i, 0, 6] = 0.25 * (last_dna[i - 2, 0, 1] + last_dna[i - 2, 2, 1] + last_dna[i - 2, 4, 1]
                                       + last_dna[i - 2, 6, 1])
        present_dna[i, 0, 7] = 0.25 * (last_dna[i - 2, 0, 1] + last_dna[i - 2, 2, 1] + last_dna[i - 2, 4, 1]
                                       + last_dna[i - 2, 6, 1])
        present_dna[i, 1, 0] = 0.25 * (last_dna[i - 2, 0, 2] + last_dna[i - 2, 2, 2] + last_dna[i - 2, 4, 2]
                                       + last_dna[i - 2, 6, 2])
        present_dna[i, 1, 1] = 0.25 * (last_dna[i - 2, 0, 2] + last_dna[i - 2, 2, 2] + last_dna[i - 2, 4, 2]
                                       + last_dna[i - 2, 6, 2])
        present_dna[i, 1, 2] = 0.25 * (last_dna[i - 2, 0, 2] + last_dna[i - 2, 2, 2] + last_dna[i - 2, 4, 2]
                                       + last_dna[i - 2, 6, 2])
        present_dna[i, 1, 3] = 0.25 * (last_dna[i - 2, 0, 2] + last_dna[i - 2, 2, 2] + last_dna[i - 2, 4, 2]
                                       + last_dna[i - 2, 6, 2])
        present_dna[i, 1, 4] = 0.25 * (last_dna[i - 2, 0, 3] + last_dna[i - 2, 2, 3] + last_dna[i - 2, 4, 3]
                                       + last_dna[i - 2, 6, 3])
        present_dna[i, 1, 5] = 0.25 * (last_dna[i - 2, 0, 3] + last_dna[i - 2, 2, 3] + last_dna[i - 2, 4, 3]
                                       + last_dna[i - 2, 6, 3])
        present_dna[i, 1, 6] = 0.25 * (last_dna[i - 2, 0, 3] + last_dna[i - 2, 2, 3] + last_dna[i - 2, 4, 3]
                                       + last_dna[i - 2, 6, 3])
        present_dna[i, 1, 7] = 0.25 * (last_dna[i - 2, 0, 3] + last_dna[i - 2, 2, 3] + last_dna[i - 2, 4, 3]
                                       + last_dna[i - 2, 6, 3])
        present_dna[i, 2, 0] = 0.25 * (last_dna[i - 2, 0, 4] + last_dna[i - 2, 2, 4] + last_dna[i - 2, 4, 4]
                                       + last_dna[i - 2, 6, 4])
        present_dna[i, 2, 1] = 0.25 * (last_dna[i - 2, 0, 4] + last_dna[i - 2, 2, 4] + last_dna[i - 2, 4, 4]
                                       + last_dna[i - 2, 6, 4])
        present_dna[i, 2, 2] = 0.25 * (last_dna[i - 2, 0, 4] + last_dna[i - 2, 2, 4] + last_dna[i - 2, 4, 4]
                                       + last_dna[i - 2, 6, 4])
        present_dna[i, 2, 3] = 0.25 * (last_dna[i - 2, 0, 4] + last_dna[i - 2, 2, 4] + last_dna[i - 2, 4, 4]
                                       + last_dna[i - 2, 6, 4])
        present_dna[i, 2, 4] = 0.25 * (last_dna[i - 2, 0, 5] + last_dna[i - 2, 4, 5] + last_dna[i - 2, 6, 5])\
                               + 0.5 * (last_dna[i - 1, 2, 5])
        present_dna[i, 2, 5] = 0.25 * (last_dna[i - 2, 0, 5] + last_dna[i - 2, 4, 5] + last_dna[i - 2, 6, 5])
        present_dna[i, 2, 6] = 0.25 * (last_dna[i - 2, 0, 5] + last_dna[i - 2, 4, 5] + last_dna[i - 2, 6, 5])
        present_dna[i, 2, 7] = 0.25 * (last_dna[i - 2, 0, 5] + last_dna[i - 2, 4, 5] + last_dna[i - 2, 6, 5])\
                               + 0.5 * (last_dna[i - 1, 2, 5])
        present_dna[i, 3, 0] = 0.25 * (last_dna[i - 2, 0, 6] + last_dna[i - 2, 2, 6] + last_dna[i - 2, 4, 6]
                                       + last_dna[i - 2, 6, 6])
        present_dna[i, 3, 1] = 0.25 * (last_dna[i - 2, 0, 6] + last_dna[i - 2, 2, 6] + last_dna[i - 2, 4, 6]
                                       + last_dna[i - 2, 6, 6])
        present_dna[i, 3, 2] = 0.25 * (last_dna[i - 2, 0, 6] + last_dna[i - 2, 2, 6] + last_dna[i - 2, 4, 6]
                                       + last_dna[i - 2, 6, 6])
        present_dna[i, 3, 3] = 0.25 * (last_dna[i - 2, 0, 6] + last_dna[i - 2, 2, 6] + last_dna[i - 2, 4, 6]
                                       + last_dna[i - 2, 6, 6])
        present_dna[i, 3, 4] = 0.25 * (last_dna[i - 2, 0, 7] + last_dna[i - 2, 2, 7] + last_dna[i - 2, 4, 7]
                                       + last_dna[i - 2, 6, 7])
        present_dna[i, 3, 5] = 0.25 * (last_dna[i - 2, 0, 7] + last_dna[i - 2, 2, 7] + last_dna[i - 2, 4, 7]
                                       + last_dna[i - 2, 6, 7])
        present_dna[i, 3, 6] = 0.25 * (last_dna[i - 2, 0, 7] + last_dna[i - 2, 2, 7] + last_dna[i - 2, 4, 7]
                                       + last_dna[i - 2, 6, 7])
        present_dna[i, 3, 7] = 0.25 * (last_dna[i - 2, 0, 7] + last_dna[i - 2, 2, 7] + last_dna[i - 2, 4, 7]
                                       + last_dna[i - 2, 6, 7])
        present_dna[i, 4, 0] = 0.25 * (last_dna[i - 2, 1, 0] + last_dna[i - 2, 3, 0] + last_dna[i - 2, 5, 0]
                                       + last_dna[i - 2, 7, 0])
        present_dna[i, 4, 1] = 0.25 * (last_dna[i - 2, 1, 0] + last_dna[i - 2, 3, 0] + last_dna[i - 2, 5, 0]
                                       + last_dna[i - 2, 7, 0])
        present_dna[i, 4, 2] = 0.25 * (last_dna[i - 2, 1, 0] + last_dna[i - 2, 3, 0] + last_dna[i - 2, 5, 0]
                                       + last_dna[i - 2, 7, 0])
        present_dna[i, 4, 3] = 0.25 * (last_dna[i - 2, 1, 0] + last_dna[i - 2, 3, 0] + last_dna[i - 2, 5, 0]
                                       + last_dna[i - 2, 7, 0])
        present_dna[i, 4, 4] = 0.25 * (last_dna[i - 2, 1, 1] + last_dna[i - 2, 3, 1] + last_dna[i - 2, 5, 1]
                                       + last_dna[i - 2, 7, 1])
        present_dna[i, 4, 4] = 0.25 * (last_dna[i - 2, 1, 1] + last_dna[i - 2, 3, 1] + last_dna[i - 2, 5, 1]
                                       + last_dna[i - 2, 7, 1])
        present_dna[i, 4, 5] = 0.25 * (last_dna[i - 2, 1, 1] + last_dna[i - 2, 3, 1] + last_dna[i - 2, 5, 1]
                                       + last_dna[i - 2, 7, 1])
        present_dna[i, 4, 6] = 0.25 * (last_dna[i - 2, 1, 1] + last_dna[i - 2, 3, 1] + last_dna[i - 2, 5, 1]
                                       + last_dna[i - 2, 7, 1])
        present_dna[i, 4, 7] = 0.25 * (last_dna[i - 2, 1, 1] + last_dna[i - 2, 3, 1] + last_dna[i - 2, 5, 1]
                                       + last_dna[i - 2, 7, 1])
        present_dna[i, 5, 0] = 0.25 * (last_dna[i - 2, 1, 2] + last_dna[i - 2, 3, 2] + last_dna[i - 2, 7, 2])\
                               + 0.5 * last_dna[i - 1, 5, 2]
        present_dna[i, 5, 1] = 0.25 * (last_dna[i - 2, 1, 2] + last_dna[i - 2, 3, 2] + last_dna[i - 2, 7, 2])
        present_dna[i, 5, 2] = 0.25 * (last_dna[i - 2, 1, 2] + last_dna[i - 2, 3, 2] + last_dna[i - 2, 7, 2])
        present_dna[i, 5, 3] = 0.25 * (last_dna[i - 2, 1, 2] + last_dna[i - 2, 3, 2] + last_dna[i - 2, 7, 2])\
                               + 0.5 * last_dna[i - 1, 5, 2]
        present_dna[i, 5, 4] = 0.25 * (last_dna[i - 2, 1, 3] + last_dna[i - 2, 3, 3] + last_dna[i - 2, 5, 3]
                                       + last_dna[i - 2, 7, 3])
        present_dna[i, 5, 5] = 0.25 * (last_dna[i - 2, 1, 3] + last_dna[i - 2, 3, 3] + last_dna[i - 2, 5, 3]
                                       + last_dna[i - 2, 7, 3])
        present_dna[i, 5, 6] = 0.25 * (last_dna[i - 2, 1, 3] + last_dna[i - 2, 3, 3] + last_dna[i - 2, 5, 3]
                                       + last_dna[i - 2, 7, 3])
        present_dna[i, 5, 7] = 0.25 * (last_dna[i - 2, 1, 3] + last_dna[i - 2, 3, 3] + last_dna[i - 2, 5, 3]
                                       + last_dna[i - 2, 7, 3])
        present_dna[i, 6, 0] = 0.25 * (last_dna[i - 2, 1, 4] + last_dna[i - 2, 3, 4] + last_dna[i - 2, 5, 4]
                                       + last_dna[i - 2, 7, 4])
        present_dna[i, 6, 1] = 0.25 * (last_dna[i - 2, 1, 4] + last_dna[i - 2, 3, 4] + last_dna[i - 2, 5, 4]
                                       + last_dna[i - 2, 7, 4])
        present_dna[i, 6, 2] = 0.25 * (last_dna[i - 2, 1, 4] + last_dna[i - 2, 3, 4] + last_dna[i - 2, 5, 4]
                                       + last_dna[i - 2, 7, 4])
        present_dna[i, 6, 3] = 0.25 * (last_dna[i - 2, 1, 4] + last_dna[i - 2, 3, 4] + last_dna[i - 2, 5, 4]
                                       + last_dna[i - 2, 7, 4])
        present_dna[i, 6, 4] = 0.25 * (last_dna[i - 2, 1, 5] + last_dna[i - 2, 3, 5] + last_dna[i - 2, 5, 5]
                                       + last_dna[i - 2, 7, 5])
        present_dna[i, 6, 5] = 0.25 * (last_dna[i - 2, 1, 5] + last_dna[i - 2, 3, 5] + last_dna[i - 2, 5, 5]
                                       + last_dna[i - 2, 7, 5])
        present_dna[i, 6, 6] = 0.25 * (last_dna[i - 2, 1, 5] + last_dna[i - 2, 3, 5] + last_dna[i - 2, 5, 5]
                                       + last_dna[i - 2, 7, 5])
        present_dna[i, 6, 7] = 0.25 * (last_dna[i - 2, 1, 5] + last_dna[i - 2, 3, 5] + last_dna[i - 2, 5, 5]
                                       + last_dna[i - 2, 7, 5])
        present_dna[i, 7, 0] = 0.25 * (last_dna[i - 2, 1, 6] + last_dna[i - 2, 3, 6] + last_dna[i - 2, 5, 6]
                                       + last_dna[i - 2, 7, 6])
        present_dna[i, 7, 1] = 0.25 * (last_dna[i - 2, 1, 6] + last_dna[i - 2, 3, 6] + last_dna[i - 2, 5, 6]
                                       + last_dna[i - 2, 7, 6])
        present_dna[i, 7, 2] = 0.25 * (last_dna[i - 2, 1, 6] + last_dna[i - 2, 3, 6] + last_dna[i - 2, 5, 6]
                                       + last_dna[i - 2, 7, 6])
        present_dna[i, 7, 3] = 0.25 * (last_dna[i - 2, 1, 6] + last_dna[i - 2, 3, 6] + last_dna[i - 2, 5, 6]
                                       + last_dna[i - 2, 7, 6])
        present_dna[i, 7, 4] = 0.25 * (last_dna[i - 2, 1, 7] + last_dna[i - 2, 3, 7] + last_dna[i - 2, 5, 7])
        present_dna[i, 7, 5] = 0.25 * (last_dna[i - 2, 1, 7] + last_dna[i - 2, 3, 7] + last_dna[i - 2, 5, 7])\
                               + 0.5 * last_dna[i - 1, 7, 7]
        present_dna[i, 7, 6] = 0.25 * (last_dna[i - 2, 1, 7] + last_dna[i - 2, 3, 7] + last_dna[i - 2, 5, 7])\
                               + 0.5 * last_dna[i - 1, 7, 7]
        present_dna[i, 7, 7] = 0.25 * (last_dna[i - 2, 1, 7] + last_dna[i - 2, 3, 7] + last_dna[i - 2, 5, 7])

    return present_dna
