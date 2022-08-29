from codec_map import Decode_Map_d2b


# Get binary sequence based on DNA sequence
def decoder_d2b(dna_data, homopolymer=3, codec_map=Decode_Map_d2b, dna_length=100):
    binary_seq = []

    # Initial state
    homopolymer_count = 1
    check_base = None

    # Read data from DNA list
    for i in range(0, len(dna_data)):
        for j in range(0, dna_length):
            if j == len(dna_data[i]):
                break
            # Determine whether the homopolymer constraints are met
            if homopolymer_count < homopolymer:
                # Convert DNA to binary
                dna_symbol = dna_data[i][j]
                present_bits = codec_map[dna_symbol]
                binary_seq.append(present_bits)

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
                binary_seq.append(present_bits)

        homopolymer_count = 1
        check_base = None

    return binary_seq
