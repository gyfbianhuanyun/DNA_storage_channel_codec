import numpy as np
from calculation_expectation import prob_calculater

#    DNA storage channel codec
#    DNA sequence constraints
#        Homopolymer = 3

# Baseline: DNA length is 3
DNA_3 = np.zeros([6, 8, 8])
DNA_3[5, :] = 1 / 64

target__DNA_length__ = 100
last_DNA = DNA_3
for length in range(4, target__DNA_length__+1):
    print(f'DNA length is {length}')
    _DNA_ = np.zeros([length * 2, 8, 8])
    _DNA_ = prob_calculater(_DNA_, last_DNA)
    last_DNA = _DNA_
    print(f'\tSum of probability: {np.sum(_DNA_)}')

    # Calculate Expectation
    expectation_ = 0
    for binary_length in range(_DNA_.shape[0]):
        expectation_ += (binary_length + 1) * np.sum(_DNA_[binary_length])
    print(f'\tExpectation: {expectation_}')

    # Calculate Mapping potential(bits/nt)
    print(f'\tMapping potential: {expectation_/length}(bits/nt)')
