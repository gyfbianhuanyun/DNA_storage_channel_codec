import numpy as np
from calculation_expectation import prob_calculater
#    DNA storage channel codec
#    DNA sequence constraints

Homopolymer = 3

# Baseline
DNA_base = np.zeros([Homopolymer * 2, 4 ** (Homopolymer - 1), 4])
DNA_base[-1, :] = 1 / (4 ** Homopolymer)

target__DNA_length__ = 100
last_DNA = DNA_base
for length in range(Homopolymer+1, target__DNA_length__+1):
    print(f'DNA length is {length}')
    _DNA_ = np.zeros([length * 2, 4 ** (Homopolymer-1), 4])
    _DNA_ = prob_calculater(_DNA_, last_DNA, Homopolymer)
    last_DNA = _DNA_
    print(f'\tSum of probability: {np.sum(_DNA_)}')

    # Calculate Expectation
    expectation_ = 0
    for binary_length in range(_DNA_.shape[0]):
        expectation_ += (binary_length + 1) * np.sum(_DNA_[binary_length])
    print(f'\tExpectation: {expectation_}')

    # Calculate Mapping potential(bits/nt)
    print(f'\tMapping potential: {expectation_/length}(bits/nt)')
