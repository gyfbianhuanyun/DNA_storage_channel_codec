import numpy as np
from calculation_MP.calculation_expectation import prob_calculater
#    DNA storage channel codec
#    DNA sequence constraints

homopolymer = 5

# Baseline
dna_base = np.zeros([homopolymer * 2, 4 ** (homopolymer - 1), 4])
dna_base[-1, :] = 1 / (4 ** homopolymer)

target__dna_length__ = 100
last_dna = dna_base

for length in range(homopolymer+1, target__dna_length__+1):
    print(f'DNA length is {length}')
    _dna_ = np.zeros([length * 2, 4 ** (homopolymer-1), 4])
    _dna_ = prob_calculater(_dna_, last_dna, homopolymer)
    last_dna = _dna_
    print(f'\tSum of probability: {np.sum(_dna_)}')

    # Calculate Expectation
    expectation_ = 0
    for binary_length in range(_dna_.shape[0]):
        expectation_ += (binary_length + 1) * np.sum(_dna_[binary_length])
    print(f'\tExpectation: {expectation_}')

    # Calculate Mapping potential(bits/nt)
    print(f'\tMapping potential: {expectation_/length}(bits/nt)')
