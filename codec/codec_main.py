from encoder import encoder_b2d_homo, encoder_b2d_gc, encoder_b2d_random_base
from decoder import decoder_d2b
from utils import write_data2file, read_dna, plot_gc_hist


def codec_processing(binary_data, opt):
    # Encoder
    print("Encode start")
    print("\tRead binary data")

    # Encoding (DNA data: list)
    binary_data_encoded, dna_data = encoder_b2d_homo(binary_data, homopolymer=opt.homopolymer_cons,
                                                     dna_length=opt.dna_length_fixed)

    dna_bases_num = len(dna_data) * opt.dna_length_fixed
    binary_bits = opt.binary_data_bits
    print("\tHomopolymer:")

    if opt.compression:
        print(f"\t\tOriginal file: {binary_bits} | compression file: {len(binary_data_encoded)}")
    else:
        print(f"\t\tOriginal file: {binary_bits} (without gzip compression)")

    print(f"\t\tMapping potential: {len(binary_data_encoded)/dna_bases_num}(bits/nt)")

    print(f"\tGC content:")
    if opt.random_base_seq:
        binary_data_encoded, dna_data, gc_content, gc_count = \
            encoder_b2d_random_base(binary_data, homopolymer=opt.homopolymer_cons, dna_length=opt.dna_length_fixed,
                                    gc_upper=opt.gc_cons_upper, gc_lower=opt.gc_cons_lower,
                                    random_seed=opt.random_seed)

    else:
        dna_data, gc_content, gc_count = \
            encoder_b2d_gc(dna_data, gc_upper=opt.gc_cons_upper, gc_lower=opt.gc_cons_lower,
                           dna_length=opt.dna_length_fixed)

    # Calculate the expected number of bases added when GC constraints are met
    sum_ = sum(gc_content)
    expected_gc = sum_ / len(gc_content)

    # Draw a histogram of the number of DNA bases added
    if opt.gc_hist:
        plot_gc_hist(gc_content)

    print(f"\t\tNum zero : {len(gc_content)} / {gc_content.count(0)}")
    print(f"\t\tMax added: {max(gc_content)}")
    print(f"\t\tTotal number of bases added: {sum_}")
    print(f"\t\tExpected value: {expected_gc}")

    dna_bases_num += sum_
    # Random binary sequence case: add the first base in each DNA sequence
    if opt.random_base_seq:
        dna_bases_num += len(gc_content)

    print(f"\tEncoding results:")
    print(f"\t\tMapping potential: {len(binary_data_encoded)/dna_bases_num}(bits/nt)")

    # Write DNA data to file
    write_data2file(dna_data, opt.encode_dna_filename)
    print("\tEncode finished")

    # Decoder
    print("Decode start")

    # Read data
    dna_data_ = read_dna(opt.encode_dna_filename)
    print("\tRead DNA data")

    # Decoding
    binary_decoder = decoder_d2b(dna_data_, homopolymer=opt.homopolymer_cons, dna_length=opt.dna_length_fixed,
                                 random_base_seq=opt.random_base_seq, random_seed=opt.random_seed)
    binary_decoder = ''.join(binary_decoder)
    binary_decoder_list = []
    binary_decoder_list[:0] = binary_decoder

    # Write binary data to file
    write_data2file(binary_decoder, opt.decode_binary_filename, mode='decode')

    print("\tDecode finished")

    if binary_decoder_list == binary_data_encoded:
        print("Codec success")
    else:
        print("Codec failed")

    print(f'Data information:\n'
          f'\tOriginal binary sequence (bits): {len(binary_data)}\n'
          f'\tEncoded  binary sequence (bits): {len(binary_data_encoded)}\n'
          f'\tDecoded  binary sequence (bits): {len(binary_decoder_list)}')
