from codec_map import Encode_Map_b2d
import math
import numpy as np


# Encoding (satisfies homopolymer constraints)
# Get DNA sequence based on binary sequence
def encoder_b2d_homo(binary_data, homopolymer=3, codec_map=Encode_Map_b2d, dna_length=100):
    dna_data_all, flag, dna_seq_last = \
        homo_encoding(homopolymer, binary_data, dna_length, codec_map)

    # When the length of the DNA sequence does not meet the fixed length,
    # delete the binary sequence corresponding to the DNA sequence
    if dna_seq_last:
        binary_data = binary_data[:flag]

    return binary_data, dna_data_all


# Encoding (satisfies GC content constraints)
# Calculate the number of bases to add to meet GC content constraints
def encoder_b2d_gc(dna_data, gc_upper=0.5, gc_lower=0.5, dna_length=100):
    gc_content_list = []
    gc_count_list = []
    dna_data_array = np.array(dna_data)
    for i in range(len(dna_data)):
        # Count the number of bases 'C' and 'G'
        added_num_symbols, add_symbol, last_symbol, gc_count = \
            calculate_added_symbols(dna_data_array[i], gc_upper, gc_lower, dna_length)

        # Add bases to meet GC content constraint
        round_ = added_num_symbols // 2
        reminder_ = added_num_symbols % 2
        gc_content_list.append(added_num_symbols)
        gc_count_list.append(gc_count[0])

        if reminder_ == 0:
            add_bases = add_symbol * round_
            dna_data[i].append(add_bases)
        else:
            add_bases = add_symbol * round_ + last_symbol
            dna_data[i].append(add_bases)

    return dna_data, gc_content_list, gc_count_list


# Homopolymer encoding
def homo_encoding(homopolymer_constraint, binary_data, dna_length, codecmap=Encode_Map_b2d):
    dna_data = []
    dna_seq = []

    # Initial state
    homopolymer_count = 1
    check_base = None
    binary_symbol = ''
    flag = 0
    # Read data from binary list
    for i in range(0, len(binary_data)):
        # Determine whether the homopolymer constraints are met
        if homopolymer_count < homopolymer_constraint:
            # Convert binary to DNA
            binary_symbol += binary_data[i]
            if binary_symbol in codecmap.keys():
                present_base = codecmap[binary_symbol]
                dna_seq.append(present_base)

                # Restoring the initial state after conversion
                binary_symbol = ''

            # Not enough binary symbols (2 bits)
            else:
                continue

            # Check homopolymer constraints
            if check_base == present_base:
                homopolymer_count += 1
            else:
                # Restore the initial state and change the checking base
                homopolymer_count = 1
                check_base = present_base

        # Homopolymer coding
        else:
            # If homopolymer constraints are met, read 1 bit of binary data
            if check_base == 'A' or check_base == 'T':
                binary_symbol = 'A' + binary_data[i]
            else:
                binary_symbol = 'C' + binary_data[i]

            # Convert by codec map
            present_base = codecmap[binary_symbol]
            # Restore the initial state and change the checking base
            check_base = present_base
            homopolymer_count = 1
            binary_symbol = ''
            dna_seq.append(present_base)

        # Fixed DNA sequence length
        if len(dna_seq) == dna_length:
            dna_data.append(dna_seq)
            homopolymer_count = 1
            check_base = None
            dna_seq = []
            flag = i + 1

    return dna_data, flag, dna_seq


# Calculate the number of added symbols to meet GC content constraint
def calculate_added_symbols(dna_data, gc_upper, gc_lower, dna_length):
    # Count the number of bases 'C' and 'G'
    bases_, count_ = np.unique(dna_data, return_counts=True)

    if count_[np.where(bases_ == 'C')] or count_[np.where(bases_ == 'G')]:
        if count_[np.where(bases_ == 'C')]:
            _gc_count_c = count_[np.where(bases_ == 'C')]
        else:
            _gc_count_c = 0

        if count_[np.where(bases_ == 'G')]:
            _gc_count_g = count_[np.where(bases_ == 'G')]
        else:
            _gc_count_g = 0

        _gc_count = _gc_count_c + _gc_count_g

    elif count_[np.where(bases_ == 'A')] or count_[np.where(bases_ == 'T')]:
        if count_[np.where(bases_ == 'A')]:
            _gc_count_a = count_[np.where(bases_ == 'A')]
        else:
            _gc_count_a = 0

        if count_[np.where(bases_ == 'T')]:
            _gc_count_t = count_[np.where(bases_ == 'T')]
        else:
            _gc_count_t = 0

        _gc_count = dna_length - (_gc_count_a + _gc_count_t)

    else:
        raise ValueError("Encoded DNA sequence error, please check binary data")

    # Too few "C" and "G" bases
    if (_gc_count / dna_length) < gc_lower:
        added_symbol = math.ceil(abs((dna_length * gc_lower - _gc_count) / (1 - gc_lower)))
        if dna_data[-1] == 'C':
            add_symbol = 'GC'
            last_symbol = 'C'
        else:
            add_symbol = 'CG'
            last_symbol = 'G'

    # Too many "C" and "G" bases
    elif (_gc_count / dna_length) > gc_upper:
        added_symbol = int(abs((dna_length * gc_upper - _gc_count) / (1 - gc_upper)))
        if dna_data[-1] == 'A':
            add_symbol = 'TA'
            last_symbol = 'A'
        else:
            add_symbol = 'AT'
            last_symbol = 'T'

    # "C" and "G" bases satisfy constraints
    else:
        added_symbol = 0
        add_symbol = ''
        last_symbol = ''

    return added_symbol, add_symbol, last_symbol, _gc_count
