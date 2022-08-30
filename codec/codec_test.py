import numpy as np
from codec_main import codec_processing


# Information of DNA storage channel codec
homopolymer_cons = 3
gc_cons_lower = 0.4
gc_cons_upper = 0.6
dna_length_fixed = 100
encode_dna_filename = 'encoder_DNA_sequence'
decode_binary_filename = 'decoder_binary_sequence'

# Generate binary data: string
binary_symbol = ['0', '1']
binary_bits = 10000000
binary_list = np.random.choice(binary_symbol, binary_bits, p=[0.5, 0.5]).tolist()

codec_processing(binary_list, homopolymer_cons, dna_length_fixed, gc_cons_upper, gc_cons_lower,
                 encode_dna_filename, decode_binary_filename)

