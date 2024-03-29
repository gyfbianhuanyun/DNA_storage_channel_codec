from codec_map import Decode_Map_d2b, first_base_list
from utils import gen_binary_seq


# Get binary sequence based on DNA sequence
def decoder_d2b(dna_data, homopolymer=3, codec_map=Decode_Map_d2b, dna_length=100,
                random_base_seq=False, random_seed=555):
    binary_seq = []
    dna_num = len(dna_data)

    # Read data from DNA list
    for i in range(0, len(dna_data)):
        # Initial state
        homopolymer_count = 1
        if random_base_seq:
            start_idx = 1
            check_base = dna_data[i][0]
        else:
            start_idx = 0
            check_base = None

        if dna_num == 1:
            last_seq = True
            end_idx = dna_length
        else:
            # Padding sequence decoding
            if i == len(dna_data) - 1:
                last_seq = True
                dna_last_seq = dna_data[i][:dna_length + 1 + start_idx]
                # find padding index
                end_idx = len(dna_last_seq) - dna_last_seq[::-1].index('A') - 1
            else:
                last_seq = False
                end_idx = dna_length+start_idx

        _binary_seq_ = []
        for j in range(start_idx, end_idx):
            if j == len(dna_data[i]):
                break
            # Determine whether the homopolymer constraints are met
            if homopolymer_count < homopolymer:
                # Check the last 2bases to
                # determine whether the last binary number is encoded in 1bit or 2bits
                if last_seq and j == end_idx - 2:
                    last_symbol = dna_data[i][j + 1]
                    if last_symbol == 'C':
                        _binary_seq_.append('0')
                    elif last_symbol == 'T':
                        _binary_seq_.append('1')
                    elif last_symbol == 'G':
                        _binary_seq_.append('00')
                    else:
                        present_bits = codec_map[dna_data[i][j]]
                        _binary_seq_.append(present_bits)
                    break

                # Convert DNA to binary
                dna_symbol = dna_data[i][j]
                present_bits = codec_map[dna_symbol]
                _binary_seq_.append(present_bits)

                # Check homopolymer count
                if check_base == dna_symbol:
                    homopolymer_count += 1
                else:
                    # Restore the initial state and change the checking base
                    homopolymer_count = 1
                    check_base = dna_symbol

            # Homopolymer decoding
            else:
                # If homopolymer constraints are met, read next base of DNA data
                if check_base == 'A' or check_base == 'T':
                    dna_symbol = 'A' + dna_data[i][j]
                else:
                    dna_symbol = 'C' + dna_data[i][j]

                # Convert by codec map
                present_bits = codec_map[dna_symbol]
                # Restore the initial state and change the checking base
                check_base = dna_data[i][j]
                homopolymer_count = 1
                _binary_seq_.append(present_bits)

        binary_decoder = ''.join(_binary_seq_)
        _binary_seq_ = []
        _binary_seq_[:0] = binary_decoder

        if random_base_seq:
            index_ = first_base_list.index(dna_data[i][0])
            # XXX binary_length is at most 2 * dna_length

            if dna_length == -1:
                binary_length = 0
            else:
                if last_seq:
                    binary_length = (len(dna_data[i]) - 2) * 2
                else:
                    binary_length = len(dna_data[i]) * 2

            binary_base_list = gen_binary_seq(binary_length, seed=random_seed)
            binary_base = binary_base_list[index_]

            for l in range(len(_binary_seq_)):
                if binary_base[l] == _binary_seq_[l]:
                    _binary_seq_[l] = '0'
                else:
                    _binary_seq_[l] = '1'

        binary_seq.extend(_binary_seq_)

    return binary_seq
