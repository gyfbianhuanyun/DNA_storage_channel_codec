import numpy as np

#    DNA storage channel codec
#    DNA sequence constraints
#        Homopolymer = k
DNA_MAP = ['A', 'C', 'G', 'T']


# Probability calculation
def prob_calculater(present_dna, last_dna, homopolymer=3, dna_map=DNA_MAP):
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

    # Get ACGT index in DNA map
    list_acgt_ = [None] * 4
    for idx in range(len(dna_map)):
        if dna_map[idx] == 'A':
            list_acgt_[0] = idx
        elif dna_map[idx] == 'C':
            list_acgt_[1] = idx
        elif dna_map[idx] == 'G':
            list_acgt_[2] = idx
        else:
            list_acgt_[3] = idx

    # Get special case row number (AAAA, CCCC, GGGG, TTTT)
    gap = sum((4**i for i in range(homopolymer-1)))
    special_case_row = np.arange(0, 4 ** (homopolymer - 1), gap)
    acgt_row = [special_case_row[list_acgt_[i]] for i in range(len(list_acgt_))]

    # Probability of getting DNA symbol based on the index
    for i in range(m, present_dna.shape[0]):
        for pre in range(4 ** (homopolymer - 1)):
            for last in range(4):
                tgt_last = pre % 4
                pre_ = pre // 4
                # Special Case based on DNA map
                # Case A
                if pre == acgt_row[0]:
                    # case 1: AAA or AAT
                    if last == list_acgt_[0] or last == list_acgt_[3]:
                        present_dna[i, pre, last] =\
                            0.25 * (last_dna[i - 2, pre_ + list_acgt_[1] * 4 ** (homopolymer-2), tgt_last]
                                    + last_dna[i - 2, pre_ + list_acgt_[2] * 4 ** (homopolymer-2), tgt_last]
                                    + last_dna[i - 2, pre_ + list_acgt_[3] * 4 ** (homopolymer-2), tgt_last])

                    # case 2: AAC or AAG
                    else:
                        present_dna[i, pre, last] = \
                            0.25 * (last_dna[i - 2, pre_ + list_acgt_[1] * 4 ** (homopolymer-2), tgt_last]
                                    + last_dna[i - 2, pre_ + list_acgt_[2] * 4 ** (homopolymer-2), tgt_last]
                                    + last_dna[i - 2, pre_ + list_acgt_[3] * 4 ** (homopolymer-2), tgt_last])\
                            + 0.5 * last_dna[i - 1, pre_ + list_acgt_[0] * 4 ** (homopolymer-2), tgt_last]

                # Case C
                elif pre == acgt_row[1]:
                    # case 1: CCC or CCG
                    if last == list_acgt_[1] or last == list_acgt_[2]:
                        present_dna[i, pre, last] = \
                            0.25 * (last_dna[i - 2, pre_ + list_acgt_[0] * 4 ** (homopolymer - 2), tgt_last]
                                    + last_dna[i - 2, pre_ + list_acgt_[2] * 4 ** (homopolymer - 2), tgt_last]
                                    + last_dna[i - 2, pre_ + list_acgt_[3] * 4 ** (homopolymer - 2), tgt_last])

                    # case 2: CCA or CCT
                    else:
                        present_dna[i, pre, last] = \
                            0.25 * (last_dna[i - 2, pre_ + list_acgt_[0] * 4 ** (homopolymer - 2), tgt_last]
                                    + last_dna[i - 2, pre_ + list_acgt_[2] * 4 ** (homopolymer - 2), tgt_last]
                                    + last_dna[i - 2, pre_ + list_acgt_[3] * 4 ** (homopolymer - 2), tgt_last])\
                            + 0.5 * last_dna[i-1, pre_ + list_acgt_[1] * 4 ** (homopolymer - 2), tgt_last]

                # Case G:
                elif pre == acgt_row[2]:
                    # case 1: GGG or GGC
                    if last == list_acgt_[1] or last == list_acgt_[2]:
                        present_dna[i, pre, last] = \
                            0.25 * (last_dna[i - 2, pre_ + list_acgt_[0] * 4 ** (homopolymer - 2), tgt_last]
                                    + last_dna[i - 2, pre_ + list_acgt_[1] * 4 ** (homopolymer - 2), tgt_last]
                                    + last_dna[i - 2, pre_ + list_acgt_[3] * 4 ** (homopolymer - 2), tgt_last])

                    # case 2: CCA or CCT
                    else:
                        present_dna[i, pre, last] = \
                            0.25 * (last_dna[i - 2, pre_ + list_acgt_[0] * 4 ** (homopolymer - 2), tgt_last]
                                    + last_dna[i - 2, pre_ + list_acgt_[1] * 4 ** (homopolymer - 2), tgt_last]
                                    + last_dna[i - 2, pre_ + list_acgt_[3] * 4 ** (homopolymer - 2), tgt_last])\
                            + 0.5 * last_dna[i-1, pre_ + list_acgt_[2] * 4 ** (homopolymer - 2), tgt_last]

                # Case T:
                elif pre == acgt_row[3]:
                    # case 1: TTT or TTA
                    if last == list_acgt_[0] or last == list_acgt_[3]:
                        present_dna[i, pre, last] = \
                            0.25 * (last_dna[i - 2, pre_ + list_acgt_[0] * 4 ** (homopolymer - 2), tgt_last]
                                    + last_dna[i - 2, pre_ + list_acgt_[1] * 4 ** (homopolymer - 2), tgt_last]
                                    + last_dna[i - 2, pre_ + list_acgt_[2] * 4 ** (homopolymer - 2), tgt_last])

                    # case 2: TTC or TTG
                    else:
                        present_dna[i, pre, last] = \
                            0.25 * (last_dna[i - 2, pre_ + list_acgt_[0] * 4 ** (homopolymer - 2), tgt_last]
                                    + last_dna[i - 2, pre_ + list_acgt_[1] * 4 ** (homopolymer - 2), tgt_last]
                                    + last_dna[i - 2, pre_ + list_acgt_[2] * 4 ** (homopolymer - 2), tgt_last]) \
                            + 0.5 * last_dna[i - 1, pre_ + list_acgt_[3] * 4 ** (homopolymer - 2), tgt_last]

                # Case 3:
                else:
                    present_dna[i, pre, last] = \
                        0.25 * (last_dna[i - 2, pre_, tgt_last]
                                + last_dna[i - 2, pre_ + 4 ** (homopolymer - 2), tgt_last]
                                + last_dna[i - 2, pre_ + 2 * 4 ** (homopolymer - 2), tgt_last]
                                + last_dna[i - 2, pre_ + 3 * 4 ** (homopolymer - 2), tgt_last])

    return present_dna
