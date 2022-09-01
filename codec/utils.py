# read data from file
def read_dna(file_path):
    dna_data = []
    with open(f'{file_path}.txt', 'r') as f:
        for line in f.readlines():
            dna_data.append(line.strip('\n'))
    return dna_data


# Write to file
def write_data2file(data, file_path, mode='encode'):
    with open(f'{file_path}.txt', 'w') as f:
        if mode == 'encode':
            for i in range(len(data)):
                seq = ''.join(data[i])
                f.write(seq + '\n')
        else:
            f.write(data)
