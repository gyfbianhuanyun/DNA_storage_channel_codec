from encoder import encoder_b2d
from decoder import decoder_d2b


# Encoder
print("Encode start")

# Binary data: string
binary_data_original = '0' * 1000
binary_list = []
binary_list[:0] = binary_data_original
print("\tRead binary data")

# DNA data: list
DNA_data = encoder_b2d(binary_list)
print("\tEncode finished")

# Decoder
print("Decode start")

# DNA data: list
dna_data_ = []
for i in range(len(DNA_data)):
    dna_data_.append(''.join(DNA_data[i]))
print("\tRead DNA data")

binary_decoder = ''.join(decoder_d2b(dna_data_))
print("\tDecode finished")

if binary_decoder == binary_data_original:
    print("Codec success")
else:
    print("Codec failed")
