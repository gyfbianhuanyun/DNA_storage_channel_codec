from codec_map import Encode_Map_b2d, first_base_list
from utils import gen_binary_seq, running_progress
import math
import time
import numpy as np


# Encoding (satisfies homopolymer constraints)
# Get DNA sequence based on binary sequence
def encoder_b2d_homo(binary_data, homopolymer=3, codec_map=Encode_Map_b2d, dna_length=100):
    dna_data_all, flag, dna_seq_last = \
        homo_encoding(homopolymer, binary_data, dna_length, codec_map)

    # When the length of the DNA sequence does not meet the fixed length, pad the DNA base
    if dna_seq_last:
        dna_seq_last = padding_dna_sequence(dna_seq_last, dna_length)
        dna_data_all.append(dna_seq_last)

    return binary_data, dna_data_all


# Encoding (satisfies GC content constraints)
# Calculate the number of bases to add to meet GC content constraints
def encoder_b2d_gc(dna_data, gc_upper=0.4, gc_lower=0.6, dna_length=100):
    gc_content_list = []
    gc_count_list = []
    # dna_data_array = np.array(dna_data)
    for i in range(len(dna_data)):
        # Count the number of bases 'C' and 'G'
        added_num_symbols, add_symbol, last_symbol, gc_count = \
            calculate_added_symbols(dna_data[i], gc_upper, gc_lower, dna_length)

        # Add bases to meet GC content constraint
        round_ = added_num_symbols // 2
        reminder_ = added_num_symbols % 2
        gc_content_list.append(added_num_symbols)
        gc_count_list.append(gc_count)

        if reminder_ == 0:
            add_bases = add_symbol * round_
            dna_data[i].append(add_bases)
        else:
            add_bases = add_symbol * round_ + last_symbol
            dna_data[i].append(add_bases)

    return dna_data, gc_content_list, gc_count_list


# Encoding (Add random binary sequence to avoid excessive GC imbalance)
def encoder_b2d_random_base(binary_data, homopolymer=3, codec_map=Encode_Map_b2d, dna_length=100,
                            gc_upper=0.4, gc_lower=0.6, random_seed=555):
    # Generate random binary sequence based on fixed seed
    binary_base_list = gen_binary_seq(dna_length, seed=random_seed, times=2)

    gc_content_list = []
    gc_count_num_list = []
    dna_data = []
    binary_original_data = []

    # Set flag for reading data
    flag = 0

    # Record time
    start_time = time.time()
    all_data = len(binary_data)
    while flag < all_data:
        # Record the results
        compare_added_symbol_list = []
        add_symbol_list = []
        last_symbol_list = []
        gc_count_list = []
        flag_list = []
        dna_seq_list = []

        # Compare the 4 cases
        for i in range(len(binary_base_list)):
            binary_base = binary_base_list[i]

            if len(binary_base) <= all_data - flag:
                ori_data = binary_data[flag:len(binary_base) + flag]
            else:
                binary_base = binary_base[:all_data - flag]

            binary_data_addition = [str(int(x) ^ int(y)) for x, y in zip(binary_base, ori_data)]

            # Add the random binary sequence to the original sequence
            first_base = first_base_list[i]
            # Get a fixed-length DNA sequences and end-of-binary sequence encoding flag
            dna_data_one_seq, _flag_, remain_seq = \
                homo_encoding(homopolymer, binary_data_addition, dna_length, codec_map,
                              random_base_seq=True, check_base=first_base)

            # Add the first base
            if dna_data_one_seq:
                dna_data_one_seq[0].insert(0, first_base)
                dna_data_one_seq_array = np.array(dna_data_one_seq)
            else:
                break

            # Calculate the number of symbols that will be added to satisfy GC content constraints
            added_num_symbols, add_symbol, last_symbol, gc_count = \
                calculate_added_symbols(dna_data_one_seq_array[0], gc_upper, gc_lower, dna_length + 1)

            # Record the results
            compare_added_symbol_list.append(added_num_symbols)
            add_symbol_list.append(add_symbol)
            last_symbol_list.append(last_symbol)
            gc_count_list.append(gc_count)
            flag_list.append(_flag_)
            dna_seq_list.append(dna_data_one_seq)

        # Compare the results, select the cases with the least bases added
        if not compare_added_symbol_list:
            break
        index_min_ = compare_added_symbol_list.index(min(compare_added_symbol_list))
        added_num_symbols = compare_added_symbol_list[index_min_]
        add_symbol = add_symbol_list[index_min_]
        last_symbol = last_symbol_list[index_min_]
        gc_count = gc_count_list[index_min_]
        dna_data_one_seq = dna_seq_list[index_min_]
        flag_check = flag_list[index_min_]

        # Update and record binary data
        binary_encoded = binary_data[flag:flag+flag_check]
        flag += flag_check

        # Add bases to meet GC content constraint
        round_ = added_num_symbols // 2
        reminder_ = added_num_symbols % 2
        gc_content_list.append(added_num_symbols)
        gc_count_num_list.append(gc_count)

        if reminder_ == 0:
            add_bases = add_symbol * round_
            dna_data_one_seq[0].append(add_bases)
        else:
            add_bases = add_symbol * round_ + last_symbol
            dna_data_one_seq[0].append(add_bases)

        dna_data.append(dna_data_one_seq[0])
        binary_original_data.extend(binary_encoded)

        # Print running progress
        encoded_data = all_data - len(binary_data)
        running_progress(all_data, encoded_data, start_time)

    return binary_original_data, dna_data, gc_content_list, gc_count_num_list


# Homopolymer encoding
def homo_encoding(homopolymer_constraint, binary_data, dna_length, codecmap=Encode_Map_b2d,
                  random_base_seq=False, check_base=None):
    dna_data = []
    dna_seq = []

    # Initial state
    homopolymer_count = 1
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
            if random_base_seq:
                break

    return dna_data, flag, dna_seq


# Calculate the number of added symbols to meet GC content constraint
def calculate_added_symbols(dna_data, gc_upper, gc_lower, dna_length):
    # Count the number of bases 'C' and 'G'
    _gc_count = sum(1 for x in dna_data if x == 'C' or x == 'G')

    # Too few "C" and "G" bases
    if (_gc_count / dna_length) < gc_lower:
        added_symbol = math.ceil((dna_length * gc_lower - _gc_count) / (1 - gc_lower))
        if dna_data[-1] == 'C':
            add_symbol = 'GC'
            last_symbol = 'G'
        else:
            add_symbol = 'CG'
            last_symbol = 'C'

    # Too many "C" and "G" bases
    elif (_gc_count / dna_length) > gc_upper:
        added_symbol = math.ceil((_gc_count - dna_length * gc_upper) / gc_upper)
        if dna_data[-1] == 'A':
            add_symbol = 'TA'
            last_symbol = 'T'
        else:
            add_symbol = 'AT'
            last_symbol = 'A'

    # "C" and "G" bases satisfy constraints
    else:
        added_symbol = 0
        add_symbol = ''
        last_symbol = ''

    return added_symbol, add_symbol, last_symbol, _gc_count


# Padding the last DNA sequence
def padding_dna_sequence(sequence, target_length=100):
    # Count the number of padding bases
    padding_count = target_length - len(sequence)

    # Padding bases
    padding_sequence = ['C' if i % 2 == 0 else 'T' for i in range(padding_count)]

    # Add padding marker base
    encoded_sequence = sequence + ['A'] + padding_sequence

    return encoded_sequence
