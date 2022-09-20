import numpy as np

#    DNA storage channel codec
#    DNA sequence constraints
#        Homopolymer = k
#    DNA_Map = ['A', 'C', 'G', 'T']


# Probability calculation
def prob_calculater(present_dna, last_dna, homopolymer=3):
    '''
    Select the last three characters in the DNA sequence and perform an inductive calculation
        DNA sequence combination of sequence length 3
        AAA meanings: The last 3 symbols are the DNA sequence of AAA
        Array map:
            0    'AAA': 0, 'AAC': 1, 'AAG': 2, 'AAT': 3,
            1    'ACA': 0, 'ACC': 1, 'ACG': 2, 'ACT': 3,
            2    'AGA': 0, 'AGC': 1, 'AGG': 2, 'AGT': 3,
            3    'ATA': 0, 'ATC': 1, 'ATG': 2, 'ATT': 3,
            4    'CAA': 0, 'CAC': 1, 'CAG': 2, 'CAT': 3,
            5    'CCA': 0, 'CCC': 1, 'CCG': 2, 'CCT': 3,
            6    'CGA': 0, 'CGC': 1, 'CGG': 2, 'CGT': 3,
            7    'CTA': 0, 'CTC': 1, 'CTG': 2, 'CTT': 3,
            8    'GAA': 0, 'GAC': 1, 'GAG': 2, 'GAT': 3,
            9    'GCA': 0, 'GCC': 1, 'GCG': 2, 'GCT': 3,
            10   'GGA': 0, 'GGC': 1, 'GGG': 2, 'GGT': 3,
            11   'GTA': 0, 'GTC': 1, 'GTG': 2, 'GTT': 3,
            12   'TAA': 0, 'TAC': 1, 'TAG': 2, 'TAT': 3,
            13   'TCA': 0, 'TCC': 1, 'TCG': 2, 'TCT': 3,
            14   'TGA': 0, 'TGC': 1, 'TGG': 2, 'TGT': 3,
            15   'TTA': 0, 'TTC': 1, 'TTG': 2, 'TTT': 3}

    Notice: if we have k homopolymer constraint (k > 3), select the last k characters in DNA
            this should be 4^{k-1} by 4 matrix
    '''

    # Get DNA sequence length
    m = int(present_dna.shape[0] / 2)

    # For the convenience of calculation, add an array of zero to last_dna (size:[1, 16, 4])
    zero_ = np.zeros([1, present_dna.shape[1], 4])
    last_dna = np.insert(last_dna, last_dna.shape[0], values=zero_, axis=0)

    # Get special case row number (AAAA, CCCC, GGGG, TTTT)
    gap = sum((4**i for i in range(homopolymer - 1)))
    acgt_row = np.arange(0, 4 ** (homopolymer - 1), gap)

    # Probability of getting DNA symbol based on the index
    for i in range(m, present_dna.shape[0]):
        for pre in range(4 ** (homopolymer - 1)):
            for last in range(4):
                tgt_last = pre % 4
                pre_ = pre // 4

                # basic case calculation
                #   for example: present_dna[i][ACA] = last_dna[i - 2][AAC + CAC + GAC + TAC]
                sum_ = 0
                for num in range(4):
                    sum_ += last_dna[i - 2, pre_ + num * 4 ** (homopolymer - 2), tgt_last]

                # Special Case based on DNA map
                if pre in acgt_row:
                    case_idx = np.where(acgt_row == pre)
                    case = case_idx[0][0]
                    sum_ = sum_ - last_dna[i - 2, pre_ + case * 4 ** (homopolymer-2), tgt_last]
                    # Case 1: AAA AAT CCC CCG GGC GGG TTA TTT
                    #   for example: present_dna[i][AAA] = 0.25 * last_dna[i - 2][CAA + GAA + TAA]
                    if last == case or last == 3 - case:
                        present_dna[i, pre, last] = 0.25 * sum_

                    # Case 2: AAC AAG CCA CCT GGA GGT TTC TTG
                    #   for example: present_dna[i][AAC] = 0.25 * last_dna[i - 2][CAA + GAA + TAA]
                    #                                      + 0.5 * last_dna[i - 1][AAA]
                    else:
                        value = last_dna[i - 1, pre_ + case * 4 ** (homopolymer - 2), tgt_last]
                        present_dna[i, pre, last] = 0.25 * sum_ + 0.5 * value

                # Case 3:
                #   for example: present_dna[i][ACA] = 0.25 * last_dna[i - 2][AAC + CAC + GAC + TAC]
                else:
                    present_dna[i, pre, last] = 0.25 * sum_

    return present_dna
