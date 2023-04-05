# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 19:34:01 2023

@author: ian.michael.bollinger@gmail.com
"""
### NUCLEOTIDE-IMAGE FUNCTIONS
import os
import sys
import hashlib
from tqdm import tqdm
from apng import APNG
from PIL import Image
import imageio

# Get Working Directory
working_directory = os.getcwd()

# Add the absolute folder path to sys.path
if working_directory not in sys.path:
    sys.path.append(working_directory)

from EncDec.encoding_decoding_funcs import tetra_bin_encode, ascii_bin_encode
from CustFasta.custom_fasta_funcs import fasta_to_dataframe

# Constants
ORIG_IMG_EXT = '.png'
APNG_EXT = '.apng'
MD5_TAG = '<'
OUTPUT_FOLDER_NAME = 'output'
EXAMINATION_FOLDER_NAME = 'examination'

def md5_checksum(file_path: str) -> str:
    with open(file_path, 'rb') as f:
        file_data = f.read()
        md5 = hashlib.md5(file_data).hexdigest()
    return md5

def create_gif_from_images(images_dir: str, gif_path: str, duration: int):
    # Get a list of the image files in the directory
    file_names = sorted(os.listdir(images_dir))
    file_names = [f for f in file_names if f.endswith('.jpg') or f.endswith('.jpeg') or f.endswith('.png')]
    
    # Read the image files into a list of image arrays
    images = []
    for file_name in file_names:
        file_path = os.path.join(images_dir, file_name)
        if file_name.endswith('.png'):
            images.append(imageio.imread(file_path, format='png'))
        else:
            images.append(imageio.imread(file_path))

    # Save the image arrays as an animated GIF
    imageio.mimsave(gif_path, images, duration=duration)

def get_rgba_values(image_path: str) -> (list, int):
    # Open the image
    img = Image.open(image_path)
    
    # Convert the image to RGBA mode
    img = img.convert("RGBA")
    
    # Get the width and height of the image
    width, height = img.size
    
    # Extract the RGBA values of the pixels
    rgba_values = [img.getpixel((col, row)) for row in range(height) for col in range(width)]
    
    return (rgba_values, width)

def process_tetrad_image(image_path: str, data: str, output_filename: str):
    # Open the image and get its dimensions
    img = Image.open(image_path)
    width, height = img.size
    # Convert the image to RGBA mode
    img = img.convert("RGBA")

    # Calculate the minimum width and height to fit the binary data
    min_width = int((len(data) // 4) ** 0.5) + 1
    min_height = int((len(data) // 4) ** 0.5) + 1

    # Calculate the position to center the data in the image
    left = (width - min_width) // 2
    top = (height - min_height) // 2

    # Get the pixel data of the image
    pixel_data = img.load()

    # Iterate through the binary data and modify the image accordingly
    for i, binary_group in enumerate([data[i:i+4] for i in range(0, len(data), 4)]):
        row, col = divmod(i, min_width)
        row += top
        col += left
        rgba = list(pixel_data[col, row])
        for j, bit in enumerate(binary_group):
            if int(bit):
                rgba[j] = rgba[j] + 127 if rgba[j] <= 127 else rgba[j] - 127
        pixel_data[col, row] = tuple(rgba)

    # Save the modified image to the specified output file
    img.save(output_filename)

def get_largest_image_size(data_list: list) -> (int, int):
    # Initialize the maximum width and height
    max_width, max_height = 0, 0
    
    # Loop through the data_list to find the largest dimensions
    for data in data_list:
        width = int((len(data) // 4) ** 0.5) + 1
        height = int((len(data) // 4) ** 0.5) + 1
        max_width = max(max_width, width)
        max_height = max(max_height, height)
        
    return(max_width, max_height)

def resize_image(input_image_path: str, output_width: int, output_height: int) -> Image:
    # Open the input image
    img = Image.open(input_image_path)
    
    # Get the dimensions of the input image
    img_width, img_height = img.size

    # Create a new white image with the desired dimensions
    new_image = Image.new('RGBA', (output_width, output_height), (255, 255, 255, 255))

    # Calculate the position to center the original image
    left = (output_width - img_width) // 2
    top = (output_height - img_height) // 2
    right = (output_width + img_width) // 2
    bottom = (output_height + img_height) // 2

    # Paste the original image onto the white background
    new_image.paste(img, (left, top, right, bottom))
    
    return new_image

def png_dir_apng_gen(input_directory: str, output_apng_path: str):  
    # Read all PNG files and sort them by name
    png_files = sorted([f for f in os.listdir(input_directory) if f.endswith('.png')])
    
    # Create an APNG object
    apng = APNG()
    
    # Add the images to the APNG object
    for png_file in png_files:
        apng.append_file(os.path.join(input_directory, png_file), delay=500)
    
    # Save the APNG file
    apng.save(output_apng_path)

def optimize_png(image, output_png_path: str):
    # Optimize the image by reducing the color palette
    optimized_image = image.convert("P", palette=Image.ADAPTIVE)

    # Save the optimized image
    optimized_image.save(output_png_path, format="PNG")

def split_apng(apng_path: str, output_folder: str):
    # Open the APNG file
    apng = APNG.open(apng_path)

    # Loop through the frames and save them as PNGs
    for frame_number, (png, control) in enumerate(apng.frames):
        # Save the frame as a PNG file
        output_path = os.path.join(output_folder, f'chrom_{frame_number:03d}.png')
        png.save(output_path)

def create_output_directory(base_path: str, folder_name: str) -> str:
    # Join the base path and folder name to create the new directory path
    path = os.path.join(base_path, folder_name)
    
    # Check if the directory already exists, if not create it
    if not os.path.exists(path):
        os.makedirs(path)
        print(f"Created directory: {path}")
    else:
        print(f"Directory {path} already exists")
        
    return path
        
# Main function
if __name__ == '__main__':
    # Input file paths
    input_fasta_file = 'C:/Users/theda/OneDrive/Documents/Python/small_ex2.fna'
    input_image_file = 'C:/Users/theda/OneDrive/Documents/Python/small_ex.png'

    # Generate Name Prefix
    output_name_prefix = os.path.splitext(os.path.basename(input_image_file))[0]
    working_directory = 'C:/Users/theda/OneDrive/Documents/Python'

    # Create output and examination directories
    output_directory = create_output_directory(working_directory, OUTPUT_FOLDER_NAME)
    examination_directory = create_output_directory(working_directory, EXAMINATION_FOLDER_NAME)

    # Output file paths
    output_fasta_file = input_fasta_file.replace(working_directory, examination_directory)
    original_image_copy = input_image_file.replace(working_directory, output_directory)
    output_apng_file =  original_image_copy.replace(ORIG_IMG_EXT, APNG_EXT)

    # Generate md5 Checksum based on input file
    generated_md5_checksum = md5_checksum(input_fasta_file)

    # Generate Dataframe from the Fasta File
    fasta_df = fasta_to_dataframe(input_fasta_file)
        
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
        output_filename = f'{output_directory}/{output_name_prefix}_chrom_{idx + 1}{ORIG_IMG_EXT}'
        process_tetrad_image(original_image_copy, data_binary, output_filename)


    # Generate an Animated Gif of the Images (DEGENERATE)
    gif_path = original_image_copy.replace('.png','.gif')
    create_gif_from_images(output_directory, gif_path, duration=0.5)

    print('\nALL MODULE QC CHECKS PASSED')