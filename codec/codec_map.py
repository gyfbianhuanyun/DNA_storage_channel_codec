# Encoding map of binary sequence data as DNA sequence data
Encode_Map_b2d = {'00': 'A', '01': 'C', '10': 'G', '11': 'T',
                  'A0': 'C', 'A1': 'G', 'C0': 'A', 'C1': 'T'}

# Decoding map of DNA sequence data as binary sequence data
Decode_Map_d2b = {'A': '00', 'C': '01', 'G': '10', 'T': '11',
                  'AC': '0', 'AG': '1', 'CA': '0', 'CT': '1'}

# When the random_base_seq option is True,
# the base list used when referring to the sequence number of the random sequence
first_base_list = ['A', 'C', 'G', 'T']
