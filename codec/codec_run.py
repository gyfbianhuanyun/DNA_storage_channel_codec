import argparse
import gzip
import numpy as np
import shutil

from codec_main import codec_processing
from utils import read_data


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1', 'True'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0', 'False'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def main(opt):
    for i in range(opt.rounds):
        print(f'The {i + 1} rounds:')
        # Generate binary data: string
        if not opt.data_filename:
            binary_symbol = ['0', '1']
            binary_bits = opt.binary_data_bits
            binary_list = np.random.choice(binary_symbol, binary_bits, p=[0.5, 0.5]).tolist()

        elif opt.data_filename:
            binary_list = read_data(opt.data_filename, opt.data_read_rb)
            opt.binary_data_bits = len(binary_list)

            # Compress data
            if opt.compression:
                with open(f'{opt.data_filename}', 'rb') as f_in:
                    with gzip.open(f'{opt.data_filename}.gz', 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                opt.data_filename += '.gz'
                binary_list = read_data(opt.data_filename, opt.data_read_rb)

        else:
            raise ValueError("Please check binary data")

        # Run codec
        codec_processing(binary_list, opt)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # Information of DNA storage channel codec
    # Data information
    parser.add_argument("--binary_data_bits", type=int, default=10000000,
                        help="Amount of binary data")
    parser.add_argument("--data_filename", type=str, default='test_file/plane.ppm',
                        help="Data file name")
    parser.add_argument("--data_read_rb", type=str2bool, default=True,
                        help="Whether to use 'rb' mode to read data")

    # DNA storage channel constraints
    parser.add_argument("--homopolymer_cons", type=int, default=3,
                        help="The DNA storage channel homopolymer constraints")
    parser.add_argument("--gc_cons_lower", type=float, default=0.4,
                        help="The lower bound probability of DNA storage channel GC content constraints")
    parser.add_argument("--gc_cons_upper", type=float, default=0.6,
                        help="The upper bound probability of DNA storage channel GC content constraints")
    parser.add_argument("--dna_length_fixed", type=int, default=100,
                        help="The fixed DNA sequence length of the DNA storage channel"
                             "If dna_length_fixed = -1, encode whole binary_data into single DNA strand")

    # Generated file in codec processing
    parser.add_argument("--encode_dna_filename", type=str, default='encoder_DNA_sequence',
                        help="Encoded DNA sequence filename")
    parser.add_argument("--decode_binary_filename", type=str, default='decoder_binary_sequence',
                        help="Decoded binary sequence filename")

    # Other codec processing options
    parser.add_argument("--rounds", type=int, default=1,
                        help="Binary to DNA sequence encoding and decoding repeated rounds")
    parser.add_argument("--gc_hist", type=str2bool, default=False,
                        help="Draw a histogram of the number of DNA bases added under the constraint of GC content")
    parser.add_argument("--compression", type=str2bool, default=False,
                        help="Whether to compress binary data")
    parser.add_argument("--random_base_seq", type=str2bool, default=False,
                        help="Add random binary sequence to avoid excessive GC imbalance issues")
    parser.add_argument("--random_seed", type=int, default=222,
                        help="The random seed to use when generating random binary sequences")

    dna_codec_opt = parser.parse_args()

    main(dna_codec_opt)
