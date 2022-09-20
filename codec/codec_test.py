import argparse
import numpy as np
from codec_main import codec_processing


def main(opt):
    for i in range(opt.rounds):
        print(f'The {i + 1} rounds:')
        # Generate binary data: string
        binary_symbol = ['0', '1']
        binary_bits = opt.binary_data_bits
        binary_list = np.random.choice(binary_symbol, binary_bits, p=[0.5, 0.5]).tolist()

        # Run codec
        codec_processing(binary_list, opt)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # Information of DNA storage channel codec
    parser.add_argument("--binary_data_bits", type=int, default=10000000,
                        help="Amount of randomly generated binary data")
    parser.add_argument("--homopolymer_cons", type=int, default=3,
                        help="The DNA storage channel homopolymer constraints")
    parser.add_argument("--gc_cons_lower", type=float, default=0.4,
                        help="The lower bound probability of DNA storage channel GC content constraints")
    parser.add_argument("--gc_cons_upper", type=float, default=0.6,
                        help="The upper bound probability of DNA storage channel GC content constraints")
    parser.add_argument("--dna_length_fixed", type=int, default=100,
                        help="The fixed DNA sequence length of the DNA storage channel")
    parser.add_argument("--encode_dna_filename", type=str, default='encoder_DNA_sequence',
                        help="Encoded DNA sequence filename")
    parser.add_argument("--decode_binary_filename", type=str, default='decoder_binary_sequence',
                        help="Decoded binary sequence filename")
    parser.add_argument("--rounds", type=int, default=1,
                        help="Binary to DNA sequence encoding and decoding repeated rounds")
    parser.add_argument("--gc_hist", type=bool, default=False,
                        help="Draw a histogram of the number of DNA bases added under the constraint of GC content")

    dna_codec_opt = parser.parse_args()
    main(dna_codec_opt)
