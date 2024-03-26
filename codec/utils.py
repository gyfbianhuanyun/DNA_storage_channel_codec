import matplotlib.pyplot as plt
import pandas as pd
import time

# read data from file
from PIL import Image
import numpy as np


# read dna data from file
def read_dna(file_path):
    dna_data = []
    with open(f'{file_path}.txt', 'r') as f:
        for line in f.readlines():
            dna_data.append(line.strip('\n'))
    return dna_data


# read binary data from picture
def read_pic(file_path):
    # Read Image
    img = Image.open(file_path)
    # Convert Image to Numpy as array
    img = np.array(img)
    img_ = np.reshape(img, img.shape[0]*img.shape[1])
    img_ = "".join([format(binary, '#010b')[2:] for binary in img_])
    img_list = []
    img_list[:0] = img_
    return img_list


# read binary data from file or gzip file
def read_gz_rb(file_path):
    with open(file_path, 'rb') as f:
        content = f.read()

    bianry_data = "".join([format(binary, '#010b')[2:] for binary in content])
    binary_list = []
    binary_list[:0] = bianry_data

    return binary_list


# read binary data from file (txt)
def read_string(file_path):
    with open(file_path, 'r') as f:
        poem_ = f.read().replace('\n', '')

    poem_ = tobinary(poem_)
    poem_list = []
    poem_list[:0] = poem_

    return poem_list


# read data from file according to the path
def read_data(file_path, read_mode_rb):
    if read_mode_rb:
        data_list = read_gz_rb(file_path)

    else:
        if file_path.endswith(('.jpg', 'ppm', 'png', 'jpeg', 'gif', 'bmp')):
            data_list = read_pic(file_path)
        elif file_path.endswith('.txt'):
            data_list = read_string(file_path)
        elif file_path.endswith('.gz'):
            data_list = read_gz_rb(file_path)
        else:
            raise ValueError("Data load error, please check data file path or data file type"
                             "(jpg/ppm/txt/gz/gif/bmp)")

    return data_list


# write data to file
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


# String to binary
def tobinary(string):
    return "".join([format(ord(char), '#010b')[2:] for char in string])


# Binary to string
def tostring(binarystring):
    return "".join([chr(int(binarystring[i:i+8], 2)) for i in range(0, len(binarystring), 8)])


# Generate random binary sequence for base sequence
def gen_binary_seq(dna_length, seed, times=2):
    binary_symbol = ['0', '1']
    binary_bits = int(dna_length * times)
    np.random.seed(seed)
    binary_base_list = np.random.choice(binary_symbol, (4, binary_bits), p=[0.5, 0.5])
    return binary_base_list


# Print running progress
def running_progress(all_data, encoded_data, start_time):
    check_time = time.time()
    percent = format(encoded_data / all_data * 100, '.2f')
    cost_time = format(check_time - start_time, '.2f')
    if not float(percent):
        percent = 0.01
    estimated_time = format(100.0 / float(percent) * float(cost_time), '.2f')
    print(f'\t\t{percent}%, time:{cost_time}s, estimated time:{estimated_time}s', end='\r')
