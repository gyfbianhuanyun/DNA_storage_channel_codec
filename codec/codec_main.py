from encoder import encoder_b2d_homo, encoder_b2d_gc
from decoder import decoder_d2b
from utils import write_data2file, read_dna


def codec_processing(binary_data, homopolymer_, dna_length_, gc_upper, gc_lower, dna_filename, binary_filename):
    # Encoder
    print("Encode start")
    print("\tRead binary data")

    # Encoding (DNA data: list)
    dna_data = encoder_b2d_homo(binary_data, homopolymer=homopolymer_, dna_length=dna_length_)

    dna_bases_num = (len(dna_data) - 1) * dna_length_ + len(dna_data[-1])
    binary_bits = len(binary_data)
    print("\tHomopolymer:")
    print(f"\t\tMapping potential: {binary_bits/dna_bases_num}(bits/nt)")

    dna_data, gc_content, gc_count = encoder_b2d_gc(dna_data, gc_upper=gc_upper, gc_lower=gc_lower,
                                                    dna_length=dna_length_)

    # Calculate the expected number of bases added when GC constraints are met
    sum_ = sum(gc_content[:len(gc_content)-1])
    expected_gc = sum_ / (len(gc_content) - 1)
    print(f"\tGC content:")
    print(f"\t\tGC count : {gc_count}")
    print(f"\t\tAdded    : {gc_content}")
    print(f"\t\tNum zero : {len(gc_content)} / {gc_content.count(0)}")
    print(f"\t\tMax added: {max(gc_content)}")
    print(f"\t\tTotal number of bases added: {sum_}")
    print(f"\t\tExpected value: {expected_gc}")

    dna_bases_num = 100 * (len(gc_content) - 1) + sum(gc_content[:-1])
    print(f"\tEncoding results:")
    print(f"\t\tMapping potential: {binary_bits/dna_bases_num}(bits/nt)")

    # Write DNA data to file
    write_data2file(dna_data, dna_filename)
    print("\tEncode finished")

    # Decoder
    print("Decode start")

    # Read data
    dna_data_ = read_dna(dna_filename)
    print("\tRead DNA data")

    dna_bases_num = 0
    for i in range(len(dna_data_)):
        dna_bases_num += len(dna_data_[i])
    print(f"\t\tMapping potential: {binary_bits/dna_bases_num}(bits/nt)")

    # Decoding
    binary_decoder = decoder_d2b(dna_data_, homopolymer=homopolymer_, dna_length=dna_length_)
    binary_decoder = ''.join(binary_decoder)
    binary_decoder_list = []
    binary_decoder_list[:0] = binary_decoder

    # Write binary data to file
    write_data2file(binary_decoder, binary_filename, mode='decode')

    print("\tDecode finished")

    if binary_decoder_list == binary_data:
        print("Codec success")
    else:
        print("Codec failed")
