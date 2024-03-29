# DNA_storage_channel_codec
We propose a simple and efficient way of DNA storage channel codec.

## Download
Clone code from GitHub though ssh
```
git@github.com:gyfbianhuanyun/DNA_storage_channel_codec.git
```

## Codec Processing
Binary to DNA sequence codec

1. Read data from file or random generate binary data
    Compress data or not
2. Encoding

    A. Homopolymer constraint encoding (Homopolymer encoder)

    B. GC content constraint encoding (GC content encoder)

    C. XOR encoding (XOR encoder)

    D. Calculate mapping potential

3. Decoding
4. Check result
    Whether the decoded binary data is equal to the original binary data


### 1. Read data
Our method can process multiple types of images(jpg, ppm, png, gif, bmp) or text(txt) data.

* For image data:
    
    We get pixel value and convert to binary data
 
* For text data:

    We get the ASCII value of the letter and convert it to binary data

* Or get the binary data directly from the file (optional):
    
    We use 'rb' mode to read the binary data of the file directly
 
Then we compress the data (lossless) using gzip (optional).

In the end, we get a list containing only binary data

### 2. Encoder

#### A. Homopolymer constraint encoder

* Step 1: Encode binary data into DNA base symbols according to ACGT diagrams {A: 00, C: 01, G: 10, T: 11}

* Step 2: When _h_ identical bases occur, read the next bit (h: homopolymer constraint)
        
    a. if the last base is A or T, then the next bit
            0: C, 1: G
    
    b. otherwise the last base is C or G, then the next 1 digit
            0: A, 1: T

* Step 3: Read the next 2 bits and repeat step 1

For example:
```
If the homopolymer constraint is 3
    it means that no more than 3 identical bases can occur

Binary data:   0101010101010
Step 1: DNA:    C C C
Step 2: The next bit is 0, so the next base is A
        DNA:    C C CA
Step 3: Repeat Step 1
        DNA:    C C CA G G G
Step 4: Encoding results
        DNA:    CCCAGGG  
```
#### B. GC content constraint encoder

Calculate the proportion of GC bases and add corresponding
bases to meet the constraint

For example:
```
If GC content constraint is 40% ~ 60%.
DNA    data:    CCCAGGG
Calculation:    0.4 <= 6 / (7 + x) <= 0.6
                  3 <= x <= 8
Added DNA  :    CCCAGGGATA
```

#### C. XOR encoder

To reduce the number of GC additions,
i.e. reduce the possibility of potentially long sequences.

* Step 1: Randomly generate 4 binary sequences using 4 bases to refer to them.
* Step 2: XORed with the original binary sequence respectively.
* Step 3: Encoding via homopolymer encoder and GC content encoder.
* Step 4: Select the case with the smallest number of additional GC bases
and place its corresponding base at the front of the sequence.

For example:
```
If GC content constraint is 40% ~ 60%.
   Homopolymer constraint is 3

Binary data     B:   01010101010100
Step 1: Generate 4 binary sequences
                A:   10011001101000
                C:   11100111010011
                G:   00100100111101
                T:   01110001100010 
Step 2: XOR with the original data
        B  XOR  A:   11001100111100
        B  XOR  C:   10110010000111
        B  XOR  G:   01110001101001
        B  XOR  T:   00100100110110
Step 3: Homopolymer encoder and GC content encoder
                A:   T A T A T T A  CGCGC
                C:   G T A G A C T  
                G:   C T A C G G C  AT
                T:   A G C A T C G  
Step 4: Select the smallest number of additional GC bases
              DNA:   G T A G A C T
        Place the 'C' at the front of the sequence
              DNA:   C G T A G A C T
```

#### D. Calculate mapping potential
$$ Mapping\  potential = {Binary\  bits \over DNA\  bases} $$

### 3. Decoder
A. If XOR encoder is used, check the first base to choice the generated binary sequence 

B. Select fixed-length DNA bases to remove added bases 

C. Decoding DNA sequence to binary sequence by ACGT diagrams

    If h identical bases appear, then the h+1th base is decoded through 2.Encoder.A.Homopolymer encoder.Step2 to obtain the binary sequence.

D. If the XOR encoder is used, the binary sequence is XORed
with the randomly generated sequence to obtain the final decoded sequence.


### 4. Check results
Check that the final decoded result is the same as the input binary sequence before encoding.

If the same, output `Codec success`, otherwise output `Codec failed`.

## Codec code
### Structure
#### 1. codec
DNA channel codec
```
Main python file
1. codec_run.py
    Codec runs python file containing arguments 
2. codec_main.py
    Python file containing codec processing
3. encoder.py
    Encoder part
4. decoder.py
    Decoder part
```

Run codec
```
Python codec_run.py --options information
```

#### 2. calculation_MP
Calculate the optimal mapping potential that our method can achieve
under given conditions (fixed DNA sequence length) and homopolymer constraints

### Options settings
```
# Data information
--binary_data_bits: Amount of binary data
--data_filename: Data file name
--data_read_rb: Use 'rb' mode to read binary data from file directly

# DNA storage channel constraints
--homopolymer_cons: The DNA storage channel homopolymer constraints
--gc_cons_lower: The lower bound probability of DNA storage channel GC content constraints
--gc_cons_upper: The upper bound probability of DNA storage channel GC content constraints
--dna_length_fixed: The fixed DNA sequence length of the DNA storage channel
                    If dna_length_fixed = -1, encode whole binary_data into single DNA strand

# Generated file in codec processing
--encode_dna_filename: Encoded DNA sequence filename
--decode_binary_filename: Decoded binary sequence filename

# Other codec processing options
--rounds: Binary to DNA sequence encoding and decoding repeated rounds
--compression: Whether to compress binary data
--gc_hist: Draw a histogram of the number of DNA bases added under the constraint of GC content
--random_base_seq: Add random binary sequence to avoid excessive GC imbalance issues (XOR encoder)
--random_seed: The seed of randomly generating binary sequence (XOR encoder)
```
