# -*- coding: utf-8 -*-
"""
Created on Sat Apr  1 18:32:39 2023

@author: ian.michael.bollinger@gmail.com
"""
### MAIN FUNCTIONS
import os
import sys
import argparse
from PIL import Image
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor

# Get Working Directory and 
working_directory = os.getcwd()

# Add the absolute folder path to sys.path
if working_directory not in sys.path:
    sys.path.append(working_directory)

from EncDec.encoding_decoding_funcs import (tetra_bin_encode,
                                            ascii_bin_encode)
from NucImg.nucleotide_image_funcs import (get_largest_image_size,
                                           resize_image,
                                           process_tetrad_image,
                                           png_dir_apng_gen,
                                           split_apng,
                                           get_rgba_values)
from NucQC.nucleotide_qc_funcs import (md5_checksum,
                                       first_qc_check,
                                       second_qc_check,
                                       final_qc_check)
from CustFasta.custom_fasta_funcs import fasta_to_dataframe

def find_file_types(directory, file_type):
    # List all files in the directory
    all_files = os.listdir(directory)
    
    # Filter the files with 'chrom' in their name and the file_type extension
    chrom_file_types = [f for f in all_files if f.endswith(file_type)]
    
    return chrom_file_types

def main(working_directory: str):
    """
    Main function to encode an image with genomic data from a FASTA file.
    """
    # Set/Get Input Files
    working_directory = f'{working_directory}/tests'
    if 'SPYDER_ARGS' in os.environ:  # Running in Spyder IDE REPLACE WITH YOUR OWN BEFORE TRYING TO RUN IN IDE        
        # # EXAMPLE DATA 5 Chromosome Organism
        # input_fasta_file = f'{working_directory}/Peltaster_fructicola.fna'
        # input_image_file = f'{working_directory}/Peltaster_fructicola.png'
        
        # EXAMPLE DATA SHORT FNA FILE
        input_fasta_file = f'{working_directory}/small_ex2.fna'
        input_image_file = f'{working_directory}/small_ex.png'
        
    # Running in console
    else:  
        parser = argparse.ArgumentParser(description="Encode an image with genomic data from a FASTA file.")
        parser.add_argument("arg1", help="Input FASTA File")
        parser.add_argument("arg2", help="Input Image File")
    
        args = parser.parse_args()
        
        input_fasta_file = args.arg1
        input_image_file = args.arg2

    # Generate Name Prefix
    output_name_prefix = os.path.splitext(os.path.basename(input_image_file))[0]
    
    output_directory = f'{working_directory}/{output_name_prefix}/output'
    examination_directory = f'{working_directory}/{output_name_prefix}/examination'
    
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
    
    output_fasta_file = f'{examination_directory}/{input_fasta_file.split("/")[-1]}'
    original_image_copy = f'{output_directory}/{input_image_file.split("/")[-1]}'
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

if __name__ == '__main__':
    main(working_directory)