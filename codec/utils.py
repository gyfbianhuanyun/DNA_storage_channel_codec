import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


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


# Draw a histogram of the number of DNA bases added under the constraint of GC content
def plot_gc_hist(_data_):
    data_array = np.array(_data_)
    data_count_list = []
    for i in range(data_array.max() + 1):
        data_count_list.append(np.sum(data_array == i))

    df = pd.DataFrame(data_count_list, columns=['Added symbol'])
    df.plot(kind='bar')
    plt.xlabel('The number of added DNA bases')
    plt.ylabel('The number of DNA sequence')
    plt.savefig('GC_histogram.pdf')
