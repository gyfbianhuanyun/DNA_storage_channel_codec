from Codec_map import Encode_Map_b2d


# Get DNA sequence based on binary sequence
def encoder_b2d(binary_data, homopolymer=3, codec_map=Encode_Map_b2d, dna_length=100):
    DNA_data_ = []
    dna_seq = []

    # Initial state
    homopolymer_count = 1
    check_base = None
    binary_symbol = ''

    # Read data from binary list
    for i in range(0, len(binary_data)):
        # Determine whether the homopolymer constraints are met
        if homopolymer_count < homopolymer:
            # Convert binary to DNA
            binary_symbol += binary_data[i]
            if binary_symbol in codec_map.keys():
                present_base = codec_map[binary_symbol]
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
            present_base = codec_map[binary_symbol]
            # Restore the initial state and change the checking base
            check_base = present_base
            homopolymer_count = 1
            binary_symbol = ''
            dna_seq.append(present_base)

        # Fixed DNA sequence length
        if len(dna_seq) == dna_length:
            DNA_data_.append(dna_seq)
            homopolymer_count = 1
            check_base = None
            dna_seq = []
    # When the DNA sequence length does not meet the fixed length
    if dna_seq:
        DNA_data_.append(dna_seq)

    return DNA_data_
