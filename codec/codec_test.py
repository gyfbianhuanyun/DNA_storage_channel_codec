import numpy as np

from encoder import encoder_b2d_homo, encoder_b2d_gc
from decoder import decoder_d2b
from utils import write_data2file, read_dna


# Information of DNA storage channel codec
homopolymer_ = 3
gc_lower = 0.4
gc_upper = 0.6
dna_length_ = 100

# Encoder
print("Encode start")

# Generate binary data: string
binary_symbol = ['0', '1']
binary_bits = 10000000
binary_list = np.random.choice(binary_symbol, binary_bits, p=[0.5, 0.5]).tolist()

print("\tRead binary data")

# Encoding (DNA data: list)
DNA_data = encoder_b2d_homo(binary_list, homopolymer=homopolymer_, dna_length=dna_length_)

DNA_bases_num = (len(DNA_data) - 1) * dna_length_ + len(DNA_data[-1])
print("\tHomopolymer:")
print(f"\t\tMapping potential: {binary_bits/DNA_bases_num}(bits/nt)")

DNA_data, gc_content, gc_count = encoder_b2d_gc(DNA_data, gc_upper=gc_upper, gc_lower=gc_lower,
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

# Write DNA data to file
DNA_filename = 'encoder_DNA_sequence'
write_data2file(DNA_data, DNA_filename)
print("\tEncode finished")

# Decoder
print("Decode start")

# Read data
dna_data_ = read_dna(DNA_filename)
print("\tRead DNA data")

DNA_bases_num = 0
for i in range(len(dna_data_)):
    DNA_bases_num += len(dna_data_[i])
print(f"\t\tMapping potential: {binary_bits/DNA_bases_num}(bits/nt)")

# Decoding
binary_decoder = decoder_d2b(dna_data_, homopolymer=homopolymer_, dna_length=dna_length_)
binary_decoder = ''.join(binary_decoder)
binary_decoder_list = []
binary_decoder_list[:0] = binary_decoder

# Write binary data to file
binary_filename = 'decoder_binary_sequence'
write_data2file(binary_decoder, binary_filename, mode='decode')

print("\tDecode finished")

if binary_decoder_list == binary_list:
    print("Codec success")
else:
    print("Codec failed")
