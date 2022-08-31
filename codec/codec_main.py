from encoder import encoder_b2d_homo, encoder_b2d_gc
from decoder import decoder_d2b
from utils import write_data2file, read_dna


def codec_processing(binary_data, opt):
    # Encoder
    print("Encode start")
    print("\tRead binary data")

    # Encoding (DNA data: list)
    binary_data, dna_data = encoder_b2d_homo(binary_data, homopolymer=opt.homopolymer_cons,
                                             dna_length=opt.dna_length_fixed)

    dna_bases_num = len(dna_data) * opt.dna_length_fixed
    binary_bits = len(binary_data)
    print("\tHomopolymer:")
    print(f"\t\tMapping potential: {binary_bits/dna_bases_num}(bits/nt)")

    dna_data, gc_content, gc_count = encoder_b2d_gc(dna_data, gc_upper=opt.gc_cons_upper, gc_lower=opt.gc_cons_lower,
                                                    dna_length=opt.dna_length_fixed)

    # Calculate the expected number of bases added when GC constraints are met
    sum_ = sum(gc_content)
    expected_gc = sum_ / len(gc_content)
    print(f"\tGC content:")
    print(f"\t\tGC count : {gc_count}")
    print(f"\t\tAdded    : {gc_content}")
    print(f"\t\tNum zero : {len(gc_content)} / {gc_content.count(0)}")
    print(f"\t\tMax added: {max(gc_content)}")
    print(f"\t\tTotal number of bases added: {sum_}")
    print(f"\t\tExpected value: {expected_gc}")

    dna_bases_num += sum_
    print(f"\tEncoding results:")
    print(f"\t\tMapping potential: {binary_bits/dna_bases_num}(bits/nt)")

    # Write DNA data to file
    write_data2file(dna_data, opt.encode_dna_filename)
    # print("\tEncode finished")

    # Decoder
    print("Decode start")

    # Read data
    dna_data_ = read_dna(opt.encode_dna_filename)
    print("\tRead DNA data")

    # Decoding
    binary_decoder = decoder_d2b(dna_data_, homopolymer=opt.homopolymer_cons, dna_length=opt.dna_length_fixed)
    binary_decoder = ''.join(binary_decoder)
    binary_decoder_list = []
    binary_decoder_list[:0] = binary_decoder

    # Write binary data to file
    write_data2file(binary_decoder, opt.decode_binary_filename, mode='decode')

    print("\tDecode finished")

    if binary_decoder_list == binary_data:
        print("Codec success")
    else:
        print("Codec failed")
