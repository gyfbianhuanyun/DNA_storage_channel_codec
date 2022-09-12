import matplotlib.pyplot as plt
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
    df = pd.DataFrame(_data_, columns=['Add symbol'])

    df.plot(kind='bar')
    plt.xlabel('DNA sequence index')
    plt.ylabel('The number of added DNA bases')
    plt.xticks([])
    plt.savefig('GC_histogram.jpg')
