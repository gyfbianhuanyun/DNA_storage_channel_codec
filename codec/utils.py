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


# read binary data from gzip file
def read_gz(file_path):
    with open(file_path, 'rb') as f:
        byte_content = f.read()

    zip_ = "".join([format(binary, '#010b')[2:] for binary in byte_content])
    zip_list = []
    zip_list[:0] = zip_

    return zip_list


# read binary data from file (txt)
def read_string(file_path):
    with open(file_path, 'r') as f:
        poem_ = f.read().replace('\n', '')

    poem_ = tobinary(poem_)
    poem_list = []
    poem_list[:0] = poem_

    return poem_list


# read data from file according to the path
def read_data(file_path):
    if file_path.endswith(('.jpg', 'ppm', 'png', 'jpeg')):
        data_list = read_pic(file_path)
    elif file_path.endswith('.txt'):
        data_list = read_string(file_path)
    elif file_path.endswith('.gz'):
        data_list = read_gz(file_path)
    else:
        raise ValueError("Data load error, please check data file path or data file type"
                         "(jpg/ppm/txt/gz)")

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


# String to binary
def tobinary(string):
    return "".join([format(ord(char), '#010b')[2:] for char in string])


# Binary to string
def tostring(binarystring):
    return "".join([chr(int(binarystring[i:i+8], 2)) for i in range(0, len(binarystring), 8)])
