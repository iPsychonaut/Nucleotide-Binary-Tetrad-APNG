# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 19:49:06 2023

@author: ian.michael.bollinger@gmail.com
"""
### QC FUNCTIONS
import hashlib
import pandas as pd
import sys
import os
from PIL import Image
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor

# SET DIRECTORIES & EXAMPLE DATA
working_directory = os.getcwd() 

# Add the absolute folder path to sys.path
if working_directory not in sys.path:
    sys.path.append(working_directory)

from EncDec.encoding_decoding_funcs import (tetra_bin_encode,
                                            ascii_bin_encode,
                                            ascii_bin_decode,
                                            tetra_bin_decode)
from NucImg.nucleotide_image_funcs import (get_largest_image_size,
                                           resize_image,
                                           process_tetrad_image,
                                           png_dir_apng_gen,
                                           split_apng,
                                           get_rgba_values)
from CustFasta.custom_fasta_funcs import (fasta_to_dataframe,
                                          reconstruct_fna_from_df)

md5_checksum_split = ascii_bin_encode('<')

def find_file_types(directory: str, file_type: str) -> list:
    # List all files in the directory
    all_files = os.listdir(directory)
    
    # Filter the files with 'chrom' in their name and the file_type extension
    chrom_file_types = [f for f in all_files if f.endswith(file_type)]
    
    return chrom_file_types

def md5_checksum(file_path: str) -> str:
    with open(file_path, 'rb') as f:
        file_data = f.read()
        md5 = hashlib.md5(file_data).hexdigest()
    return md5

def first_qc_check(input_fasta_file: str, first_check_index: int, output_encoded_image_path: str):
    generated_md5_checksum = md5_checksum(input_fasta_file)
    original_rgba_values, width = get_rgba_values(original_image_copy)
    encoded_rgba_values, _ = get_rgba_values(output_encoded_image_path)

    first_check_binary_data = ''.join(
        ''.join('0' if orig == enc else '1' for orig, enc in zip(orig_rgba, enc_rgba))
        for orig_rgba, enc_rgba in zip(original_rgba_values, encoded_rgba_values))

    if first_check_index == 0:
        pass
    else:
        encoded_string_split = first_check_binary_data.split(ascii_bin_encode('<'), maxsplit = 4)
        
        decoded_id_desc = ascii_bin_decode(encoded_string_split[0])
        decoded_md5_checksum = ascii_bin_decode(encoded_string_split[1])
        decoded_nucleotide_type = ascii_bin_decode(encoded_string_split[2])
        decoded_encoding_key = ascii_bin_decode(encoded_string_split[3])
        
        decoded_sequence = tetra_bin_decode(encoded_string_split[4], decoded_encoding_key)
        if decoded_nucleotide_type == 'RNA':
            if decoded_encoding_key == 'degenerate':
                decoded_sequence = decoded_sequence.replace('T','U')
            else:
                decoded_sequence = decoded_sequence.replace('t','u').replace('T','U')
    
    if decoded_md5_checksum == generated_md5_checksum:
        print(f'{output_encoded_image_path}\nPASSES FIRST QC: MD5 IMAGE ENCODING/DECODING')
    else:
        print(f'MD5 CHECKSUMS FAILED FIRST QC\nCHECK FILE INTEGRITY FOR {output_encoded_image_path}')

def second_qc_check(first_check_index: int, output_encoded_image_path: str, apng_image_list: list, original_rgba_values: list, input_fasta_file: str):
    generated_md5_checksum = md5_checksum(input_fasta_file)
    encoded_rgba_values, _ = get_rgba_values(output_encoded_image_path)    
    second_check_binary_data = ''
    difference_count = 0
    for i, rgba in enumerate(original_rgba_values):
        temp_tetrad = ''
        for n, item in enumerate(rgba):
            if rgba[n] == encoded_rgba_values[i][n]:        
                temp_tetrad += '0'
            else:
                temp_tetrad += '1'
                difference_count += 1
        second_check_binary_data += temp_tetrad
    if first_check_index == 0:
        pass
    else:    
        second_check_binary_data = remove_zero_chunks(second_check_binary_data)
        encoded_string_split = second_check_binary_data.split(ascii_bin_encode('<'), maxsplit = 4)
        try:
            decoded_ID_description = ascii_bin_decode(encoded_string_split[0])
            decoded_ID = decoded_ID_description.split(' ', 1)[0]
            decoded_description = decoded_ID_description.split(' ', 1)[1]
        except IndexError:
            decoded_ID_description = ascii_bin_decode(encoded_string_split[0][4:])
            decoded_ID = decoded_ID_description.split(' ', 1)[0]
            decoded_description = decoded_ID_description.split(' ', 1)[1]
        decoded_md5_checksum = ascii_bin_decode(encoded_string_split[1])
        decoded_nucleotide_type = ascii_bin_decode(encoded_string_split[2])
        decoded_encoding_key = ascii_bin_decode(encoded_string_split[3])
        decoded_sequence = tetra_bin_decode(encoded_string_split[4], decoded_encoding_key)
        if decoded_nucleotide_type == 'RNA':
            if decoded_encoding_key == 'degenerate':
                decoded_sequence = decoded_sequence.replace('T','U')
            elif decoded_encoding_key == 'confidence':
                decoded_sequence = decoded_sequence.replace('t','u').replace('T','U')
        
        if decoded_md5_checksum == generated_md5_checksum:
            print(f'{output_encoded_image_path}\nPASSES SECOND QC: MD5 APNG ENCODING/DECODING')
            return decoded_ID, decoded_description, decoded_sequence
        else:
            print(f'MD5 CHECKSUMS FAILED SECOND QC;\nCHECK FILE INTEGRITY FOR {output_encoded_image_path}')
            return None

def final_qc_check(extracted_results: list, output_fasta_file: str):
   
    # Filter out None values from results
    extracted_results = [extracted_result for extracted_result in extracted_results if extracted_result is not None]
    
    # Create DataFrame from results
    extracted_df = pd.DataFrame(extracted_results, columns=['ID', 'Description', 'Sequence']).set_index('ID')
    
    # Reconstruct FNA file from Results Dataframe
    reconstruct_fna_from_df(extracted_df, output_fasta_file)
    
    # Genearte an md5 checksum for the reconstructed file
    generated_md5_checksum = md5_checksum(output_fasta_file)
    
    if verify_file(output_fasta_file, generated_md5_checksum):
        print('PASSES FINAL QC: MD5 FNA DECODING CHECKSUM')
    else:
        # print('CORRECTING FOR ENDFILE-NEWLINE ERROR')
        # Read the file contents into a variable
        with open(output_fasta_file, 'r') as file:
            contents = file.readlines()
        
        # Check if the last character is a newline and remove it if it is
        if contents and contents[-1].endswith('\n'):
            contents[-1] = contents[-1].rstrip('\n')
        
        # Write the modified contents back to the file
        with open(output_fasta_file, 'w') as file:
            file.writelines(contents)
        if verify_file(output_fasta_file, generated_md5_checksum):
            print('PASSES FINAL QC: MD5 FNA DECODING CHECKSUM')
        else:
            print(f'MD5 CHECKSUMS FAILED FINAL QC;\nCHECK FILE INTEGRITY FOR {output_fasta_file}')
 
def verify_file(input_file_path: str, expected_md5_checksum: str) -> bool:
    # Calculate the input file's md5 checksum
    calculated_md5_checksum = md5_checksum(input_file_path)
    
    # Return TRUE/FALSE if it matched provided md5
    return calculated_md5_checksum == expected_md5_checksum

def remove_zero_chunks(input_string: str) -> str:
    # Remove buffer space from original image encoded as chunks of 0's
    chunks = [input_string[i:i + 12] for i in range(0, len(input_string), 12)]
    filtered_chunks = [chunk for chunk in chunks if chunk != '000000000000']
    
    return ''.join(filtered_chunks)

if __name__ == '__main__':
    # Set/Get Input Files
    if 'SPYDER_ARGS' in os.environ:  # Running in Spyder IDE REPLACE WITH YOUR OWN BEFORE TRYING TO RUN IN IDE
        working_directory = 'C:/Users/theda/OneDrive/Documents/Python/example_genome'
        
        # # EXAMPLE DATA 5 Chromosome Organism
        input_fasta_file = f'{working_directory}/GCA_001592805.2_ASM159280_genomic.fna'
        input_image_file = f'{working_directory}/Peltaster_fructicola.png'
        
        # EXAMPLE DATA SHORT FNA FILE
        # input_fasta_file = f'{working_directory}/small_ex2.fna'
        # input_image_file = f'{working_directory}/small_ex.png'
        
    else:  # Running in console
        parser = argparse.ArgumentParser(description="Encode an image with genomic data from a FASTA file.")
        parser.add_argument("arg1", help="Input FASTA File")
        parser.add_argument("arg2", help="Input Image File")
    
        args = parser.parse_args()
        
        input_fasta_file = args.arg1
        input_image_file = args.arg2

    # Generate Name Prefix
    output_name_prefix = os.path.splitext(os.path.basename(input_image_file))[0]
    
    output_directory = f'{working_directory}/output'
    examination_directory = f'{working_directory}/examination'
    
    # Check if output_directory exists, if not, create it
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
        print(f"Created directory: {output_directory}")
    else:
        print(f"Directory {output_directory} already exists")
    
    # Check if examination_directory exists, if not, create it
    if not os.path.exists(examination_directory):
        os.makedirs(examination_directory)
        print(f"Created directory: {examination_directory}")
    else:
        print(f"Directory {examination_directory} already exists")
    
    output_fasta_file = input_fasta_file.replace(working_directory, examination_directory)
    original_image_copy = input_image_file.replace(working_directory, output_directory)
    output_apng_file =  original_image_copy.replace('.png', '.apng')
    
    # Generate md5 Checksum based on input file
    generated_md5_checksum = md5_checksum(input_fasta_file)
    
    # Generate Dataframe from the Fasta File
    fasta_df = fasta_to_dataframe(input_fasta_file)
    for column in fasta_df.columns:
        fasta_df[column] = fasta_df[column] #.str.upper() ### POINT OF DEGENERACY ###
    
    # Open the original image and convert it to a palette-based format with 256 colors
    image = Image.open(input_image_file)
    palette_image = image.convert("P", palette=Image.ADAPTIVE, colors=256)
    
    # Save the palette image to the output directory
    palette_image.save(original_image_copy)
    
    # Determine the largest image size needed to fit the largest chromosome data
    binary_data_list = []
    for idx, row in tqdm(fasta_df.iterrows(), total=fasta_df.shape[0], desc='Processing Sequences', ncols=100):
        description = row['Description']
        sequence = row['Sequence']
        sequence_binary, nucleotide_type, encoding_key = tetra_bin_encode(sequence)
        md5_desc_type_enc = f'{description}<{generated_md5_checksum}<{nucleotide_type}<{encoding_key}<'
        md5_desc_type_enc_binary = ascii_bin_encode(md5_desc_type_enc)
        data_binary = md5_desc_type_enc_binary + sequence_binary
        binary_data_list.append(data_binary)
    
    # Determine the largest image needed for encoding
    max_width, max_height = get_largest_image_size(binary_data_list)
    
    # Resize the original image copy to the largest image size
    img_resized = resize_image(original_image_copy, max_width, max_height)
    img_resized.save(original_image_copy)
    
    # Encode all subsequent chromosome data using the resized image
    for idx, data_binary in tqdm(enumerate(binary_data_list), total=len(binary_data_list), desc="Encoding chromosomes", ncols=100):
        output_filename = f'{output_directory}/{output_name_prefix}_chrom_{idx + 1}.png'
        process_tetrad_image(original_image_copy, data_binary, output_filename)
    
    # Generate the APNG
    png_dir_apng_gen(output_directory, output_apng_file)

    # QUALITY CONTROL CHECKS
    split_apng(output_apng_file, examination_directory)

    # FIRST QC Check
    print('\nSTARTING FIRST QC CHECK: IMAGE ENCODING/DECODING')
    encoded_image_list = find_file_types(output_directory, '.png')
    encoded_image_list = [f'{output_directory}/{image_path}' for image_path in encoded_image_list]
    
    with ThreadPoolExecutor() as executor:
        # Use executor.map() to call first_qc_check with these arguments
        executor.map(first_qc_check, [input_fasta_file]*len(encoded_image_list), range(len(encoded_image_list)), encoded_image_list)

    # SECOND QC Check
    print('\nSTARTING SECOND QC CHECK: ANIMATED PNG ENCODING/DECODING')
    apng_image_list = find_file_types(examination_directory, '.png')
    apng_image_list = [f'{examination_directory}/{image_path}' for image_path in apng_image_list]
    original_rgba_values, width = get_rgba_values(apng_image_list[0])
    
    with ThreadPoolExecutor() as executor:
        extracted_results = list(executor.map(second_qc_check, range(len(encoded_image_list)), encoded_image_list, [apng_image_list] * len(encoded_image_list), [original_rgba_values] * len(encoded_image_list), [input_fasta_file] * len(encoded_image_list)))
    
    # FINAL QC Check
    print('\nSTARTING FINAL QC CHECK: FNA DECODING')
    final_qc_check(extracted_results, output_fasta_file)