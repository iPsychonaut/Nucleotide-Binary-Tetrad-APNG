# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 19:58:56 2023

@author: ian.michael.bollinger@gmail.com
"""
### FASTA FUNCTIONS
from Bio import SeqIO
import pandas as pd
from io import StringIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def fasta_to_dataframe(fasta_file: str):
    # Parse the FASTA file and store the records in a list of dictionaries
    records = SeqIO.parse(fasta_file, 'fasta')
    data = [{'ID': record.id, 'Description': record.description, 'Sequence': str(record.seq)} for record in records]
    
    # Convert the list of dictionaries into a DataFrame
    df = pd.DataFrame(data)
    return df

def reconstruct_fna_from_df(df, output_file_path: str):
    # Generate an .fna file-type from a given dataframe
    with open(output_file_path, 'w') as f:
        for index, row in df.iterrows():
            description = row['Description']
            sequence = row['Sequence']
            
            # Remove trailing spaces from sequence
            sequence = sequence.rstrip()
            
            # Write record to buffer
            buffer = StringIO()
            SeqIO.write([SeqRecord(Seq(sequence), id=str(index), description=description)], buffer, 'fasta')
            buffer.seek(0)
            
            # Parse record from buffer and write to output file
            for record in SeqIO.parse(buffer, 'fasta'):
                f.write(record.format('fasta'))