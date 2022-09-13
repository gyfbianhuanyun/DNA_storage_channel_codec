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

    ax = df.plot(kind='bar')
    plt.title('Distribution map of added DNA bases under the GC content constraint')
    plt.xlabel('The number of added DNA bases')
    plt.xlim(left=-1)
    plt.ylabel('The number of DNA sequence')

    # Add column value (The number of DNA sequence)
    rects = ax.patches
    for rect, label in zip(rects, data_count_list):
        ax.text(rect.get_x() + rect.get_width() / 2, rect.get_height() + max(data_count_list) / 50,
                label, ha="center", fontsize=8)

    plt.savefig('output_files/GC_histogram.pdf')
