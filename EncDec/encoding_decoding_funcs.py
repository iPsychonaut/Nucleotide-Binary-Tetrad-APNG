# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 19:26:28 2023

@author: ian.michael.bollinger@gmail.com
"""
### ENCODING/DECODING FUNCTIONS

# Define the tetrabin encoding scheme for nucleotides
encoding_schemes = {'degenerate' : {'U': '0001', 'T': '0001', 'A': '0010', 'C': '0100', 'G': '1000', 
                                  'R': '1010', 'Y': '0101', 'K': '1001', 'M': '0110',
                                  'S': '1100', 'W': '0011', 'B': '1101', 'D': '1011',
                                  'H': '0111', 'V': '1110', 'N': '1111', ' ': '0000'},
                    'confidence' : {'U': '0001', 'u': '1110',
                                    'T': '0001', 't': '1110',
                                    'A': '0010', 'a': '1101',
                                    'C': '0100', 'c': '1011',
                                    'G': '1000', 'g': '0111',
                                    ' ': '0000',
                                    '<open1>': '0011', '<open2>': '1100',
                                    '<open3>': '1001', '<open4>': '0110',
                                    '<open5>': '1010', '<open6>': '0101',}}

def reverse_dict(input_dict: dict) -> dict:
    # Reverse the Keys and Values for a given Dictionary
    reversed_dict = {v: k for k, v in input_dict.items()}
    
    return reversed_dict


def ascii_bin_decode(input_string) -> str:
    # Split the binary string into groups of 8 bits
    byte_list = [input_string[i:i+8] for i in range(0, len(input_string), 8)]
    
    # Convert each group of 8 bits back to its decimal value
    decimal_list = [int(byte, 2) for byte in byte_list]
    
    # Convert each decimal value to its corresponding ASCII character
    char_list = [chr(decimal) for decimal in decimal_list]
    
    # Concatenate all the ASCII characters to form the original input string
    output_string = ''.join(char_list)
    
    return output_string


def ascii_bin_encode(input_string: str):
    # Convert each character in the input string to its corresponding 8-bit ASCII binary code
    encoded_ascii_bin = ''.join(format(ord(char), '08b') for char in input_string)
    
    return encoded_ascii_bin


def tetra_bin_encode(input_sequence: str) -> (str, str, str): 
    # Determine encoding scheme based on contents
    encoding_key, nucleotide_type = fasta_encoding_check(input_sequence)
    encoding_scheme = encoding_schemes[encoding_key]
    
    # Remove any new line characters
    input_sequence = input_sequence.replace('\n','')
       
    # Encode the nucleotide sequence using tetrabin encoding
    try:
        encoded_sequence = ''.join(encoding_scheme[n] for n in input_sequence)
    except KeyError:
        print(f'Invalid degenerate nucleotide sequence: {input_sequence}\nTrying confidence encoding')
        try:
            encoding_scheme = encoding_schemes['confidence']
            encoded_sequence = ''.join(encoding_scheme[n] for n in input_sequence)
        except KeyError:            
            print(f'Invalid nucleotide sequence: {input_sequence}\nTrying confidence encoding')
            return None
        
    return(encoded_sequence, nucleotide_type, encoding_key)

def tetra_bin_decode(final_encoded_string: str, encoding_scheme: dict) -> str:
    # Tetrabin decoding scheme
    tetrabin_decoding_scheme = reverse_dict(encoding_schemes[encoding_scheme])

    # Decode the description and sequence
    decoded_sequence = ''.join(tetrabin_decoding_scheme[final_encoded_string[i:i+4]] for i in range(0, len(final_encoded_string), 4))
    return decoded_sequence

def fasta_encoding_check(input_sequence: str) -> (str, str):   
    # Determine if case-based (confidence) nucleotide data
    if any(char in input_sequence for char in ['u', 't', 'a', 'c', 'g']):
        encoding_key = 'confidence'

    # Determine if letter-based (degenerate) nucleotide data
    elif any(char in input_sequence for char in ['R', 'Y', 'K', 'M', 'S', 'W', 'B', 'D', 'H', 'V', 'N']):
        encoding_key = 'degenerate'
        
    # Determine the type of nucleotide (DNA or RNA) in the sequence
    if any(char in input_sequence for char in ['U', 'u']):
        nucleotide_type = 'RNA'
    else:
        nucleotide_type = 'DNA'
    
    return (encoding_key, nucleotide_type)


if __name__ == '__main__':
    # dna_seq_degen = 'TACGBDHV'
    # rna_seq_degen = 'UACGBDHV'
    
    # dna_seq_conf = 'taCGTAcg'
    rna_seq_conf = 'uaCGUAcg'
    
    test_id_desc = '>200001.1 Example Organism'
    md5_test = 'asd3f578e53asdf5'
    
    input_sequence = rna_seq_conf
    
    encoded_binary, nucleotide_type, encoded_key = tetra_bin_encode(input_sequence)
    
    md5_desc_type_enc = f'{test_id_desc}<{md5_test}<{nucleotide_type}<{encoded_key}<'
    
    encoded_ascii_bin = ascii_bin_encode(md5_desc_type_enc)
    
    final_encoded_string = encoded_ascii_bin + encoded_binary
    
    encoded_string_split = final_encoded_string.split(ascii_bin_encode('<'), maxsplit = 4)
    
    decoded_id_desc = ascii_bin_decode(encoded_string_split[0])
    decoded_md5 = ascii_bin_decode(encoded_string_split[1])
    decoded_nucleotide_type = ascii_bin_decode(encoded_string_split[2])
    decoded_encoding_key = ascii_bin_decode(encoded_string_split[3])
    
    decoded_sequence = tetra_bin_decode(encoded_string_split[4], decoded_encoding_key)
    if decoded_nucleotide_type == 'RNA':
        if decoded_encoding_key == 'degenerate':
            decoded_sequence = decoded_sequence.replace('T','U')
        else:
            decoded_sequence = decoded_sequence.replace('t','u').replace('T','U')
            
    if test_id_desc == decoded_id_desc:
        if md5_test == decoded_md5:
            if input_sequence == decoded_sequence:
                print('ALL MODULE QC CHECKS PASSED')