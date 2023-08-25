#!/usr/bin/env python
# coding: utf-8

# In[1]:


from Bio import Entrez
Entrez.email = 'A.N.Other@example.com'

import Bio
import mygene
import csv
import re
import datetime
import os
import requests

from tqdm import tqdm
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation


# In[2]:


# TargetSeqLibrary class contains methods which construct a library of valid genes and their target sequences in the form of a CSV file.
# Genes which are not found in the NCBI database or inputted PAM is not found in select exons of a gene are added to a library of invalid genes in the form of a CSV file.
class TargetSeqLibrary:

    # Constructor for TargetSeqLibrary
    def __init__(self, input_CSV, input_species, input_PAM_seq, input_repeat_seq, input_length_of_target_seq, output_valid_CSV_file_path, output_invalid_CSV_file_path):
        self.CSV_of_genes = open(input_CSV, 'r')
        self.species = input_species
        self.PAM = input_PAM_seq
        self.repeat_seq = input_repeat_seq
        self.length_target_seq = input_length_of_target_seq
        self.valid_target_seq_CSV_path = output_valid_CSV_file_path
        self.invalid_gene_CSV_path = output_invalid_CSV_file_path
        
        self.valid_gene_dict = {}
        self.invalid_gene_dict = {'gene_ids': []}
    
    # Query gene accession numbers and exon ranges from NCBI for each gene in inputted CSV file
    def query_accession_numbers(self):
        mg = mygene.MyGeneInfo()
        exon_list = mg.querymany(self.CSV_of_genes, scopes = 'symbol', fields = 'exons', species = self.species, returnall = True)
        self.add_not_found_in_database_gene_to_target_seqs(exon_list['missing'])
        return exon_list
    
    # Add gene with valid target sequence into valid_gene_dict
    def add_gene_to_found_seqs(self, first_variant_of_found_gene):
        self.valid_gene_dict.update({first_variant_of_found_gene.gene_id : {'gRNA':first_variant_of_found_gene.return_gRNA_seq()}})
        
    # Add gene where suitable target sequence not found in select exons into invalid_gene_dict
    def add_not_found_target_seq_gene_to_target_seqs(self, first_variant_of_not_found_gene):
        if first_variant_of_not_found_gene.gene_id not in self.invalid_gene_dict['gene_ids']:
            self.invalid_gene_dict['gene_ids'].append(first_variant_of_not_found_gene.gene_id)
        
    # Add gene which was not found in NCBI database into invalid_gene_dict
    def add_not_found_in_database_gene_to_target_seqs(self, list_of_not_found_genes):
        for gene_id in list_of_not_found_genes:
            if gene_id not in self.invalid_gene_dict['gene_ids']:
                self.invalid_gene_dict['gene_ids'].append(gene_id)
    
    # Iterate through exon_list attained from NCBI, check if each gene contains a suitable target sequence, and add gene sequence to respective dictionary
    def check_every_gene_for_target_seq(self):
        # Wrap the loop with tqdm to display the progress bar
        for gene in tqdm(self.query_accession_numbers()['out'], total=len(self.query_accession_numbers()['out'])):
            if 'exons' in gene:
                all_exon_ranges = gene['exons']
                gene_id = gene['query']
                accession_number_of_first_variant = (all_exon_ranges[0])['transcript']
                first_variant_strand_directionality = (all_exon_ranges[0])['strand']
                first_variant = SingleGeneSeq(gene_id, accession_number_of_first_variant, first_variant_strand_directionality, self.PAM, self.repeat_seq, self.length_target_seq)

                target_seq = first_variant.search_for_PAM()
                
                # Check if a valid target sequence was found and if the length of found gRNA is the same as the inputted length of the target sequence plus the length of the inputted repeat sequence.
                # The second condition may occur if the PAM is found too close to the start or end of an exon, meaning that the target sequence will be shorter than the inputted length of the target sequence.
                if isinstance(target_seq, Bio.SeqRecord.SeqRecord) and first_variant.return_length_of_gRNA_seq() == (self.length_target_seq + len(self.repeat_seq)) :
                    self.add_gene_to_found_seqs(first_variant)
                    
                else: 
                    if first_variant.gene_id not in self.invalid_gene_dict['gene_ids']:
                        self.add_not_found_target_seq_gene_to_target_seqs(first_variant)

        self.save_CSV_of_valid_genes(self.valid_target_seq_CSV_path)
        self.save_CSV_of_invalid_genes(self.invalid_gene_CSV_path)
        
    # Write valid_gene_dict into a CSV 
    def save_CSV_of_valid_genes(self, file_path):
        try:
            with open(file_path, 'w', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(['gene id', 'gRNA'])
                for gene_id, values in self.valid_gene_dict.items():
                    writer.writerow([gene_id, values['gRNA']])
        except IOError:
            print("An error occurred while writing to a CSV file.")
    
    # Write invalid_gene_dict into a CSV      
    def save_CSV_of_invalid_genes(self, file_path):
        try:
            with open(file_path, 'w', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(['gene_id'])
                for gene_id in self.invalid_gene_dict['gene_ids']:
                    writer.writerow([gene_id])
        except IOError:
            print("An error occurred while writing to a CSV file.")


# In[3]:


# SingleGeneSeq class contains methods which construct a gene and find a valid target sequence for the gene's first variant.
class SingleGeneSeq:
  
    # Constructor for SingleGeneSeq
    def __init__(self, gene_id, accession_number_of_first_variant, first_variant_strand_directionality, input_PAM_seq, input_repeat_seq, input_length_of_target_seq):
        self.gene_id = gene_id
        self.accession_number_of_first_variant = accession_number_of_first_variant
        self.first_variant_strand_directionality = first_variant_strand_directionality
        self.input_PAM_seq = input_PAM_seq
        self.input_repeat_seq = input_repeat_seq
        self.input_length_of_target_seq = input_length_of_target_seq
        
        self.record = None
        self.start_pos = None
        self.prev_pos_of_PAM = 0
        self.seq_to_search = None
        self.range_of_first_exon_with_start_of_CDS = None
        self.target_seq = None
        self.target_seq_length = 0
        self.start_direction_of_PAM_search = None
        self.start_index = 0
        self.end_of_CDS_range_in_variant = 0
        
    # Retrieve details from NCBI based on a gene's first variant. These details include the exon and coding sequence range specific to the first variant.
    def fetch_details(self, accession_number_of_first_variant):
        handle = Entrez.efetch(db = "nucleotide", id = accession_number_of_first_variant, rettype = "gb", retmode = "text")
        details = SeqIO.read(handle, "genbank")
        handle.close()
        return details
 
    # Find first exon range which includes the beginning of the coding sequence
    def find_range_of_seq_to_search_with_CDS(self):
        start_of_CDS_range_of_variant = None
        exon_range_of_variant = None
        
        self.record = self.fetch_details(self.accession_number_of_first_variant)
        
        # Iterating through features of variant to find range of coding strand
        for feat in self.record.features:
            if feat.type == 'CDS':
                start_of_CDS_range_of_variant = feat.location.start + 1
                self.end_of_CDS_range_in_variant = feat.location.end
                break
            
        if start_of_CDS_range_of_variant is None or self.end_of_CDS_range_in_variant is None:
            return None

        # Iterating through features of variant to find range of exon which contains start of coding strand
        for feat in self.record.features:
            if feat.type == 'exon' and feat.location.end > start_of_CDS_range_of_variant:
                exon_range_of_variant = feat.location
                break
            
        if exon_range_of_variant != None:
            # Check that end of exon range is within the range of the CDS
            if exon_range_of_variant.end > self.end_of_CDS_range_in_variant:
                self.range_of_first_exon_with_start_of_CDS = FeatureLocation(start_of_CDS_range_of_variant, self.end_of_CDS_range_in_variant)
                return self.range_of_first_exon_with_start_of_CDS
            
            self.range_of_first_exon_with_start_of_CDS = FeatureLocation(start_of_CDS_range_of_variant, exon_range_of_variant.end)
            return self.range_of_first_exon_with_start_of_CDS
        
        else:
            return None
    
    
    def find_range_of_seq_to_search_past_CDS(self):
        next_exon_range_of_variant = None
        self.record = self.fetch_details(self.accession_number_of_first_variant)
        
        for feat in self.record.features:
            # Due to the transition from 1-based indexing of the NCBI database to 0-based indexing of Python for a 
            # FeatureLocation object, add 1 to the starting position of an exon to account for this when comparing two FeatureLocations.
            if feat.type == 'exon' and (feat.location.start + 1) > self.range_of_first_exon_with_start_of_CDS.end:
                next_exon_range_of_variant = feat.location
                break

        if next_exon_range_of_variant != None and next_exon_range_of_variant.start < self.end_of_CDS_range_in_variant:
            if next_exon_range_of_variant.end > self.end_of_CDS_range_in_variant:
                next_exon_location = FeatureLocation(next_exon_range_of_variant.start, self.end_of_CDS_range_in_variant)
                return next_exon_location
            
            next_exon_location = FeatureLocation(next_exon_range_of_variant.start, next_exon_range_of_variant.end)
            return next_exon_location
        
        else:
            return None
    
    # Extract the sequence of a specific range to look for PAM.
    def retrieve_seq_to_search_for_PAM(self):
        feat_location = self.find_range_of_seq_to_search_with_CDS()
        self.range_of_first_exon_with_start_of_CDS = feat_location

        if feat_location is None:
            return None
        
        # Extracting the selected exon sequence
        exon_seq_to_extract = SeqFeature(feat_location)
        self.seq_to_search = exon_seq_to_extract.extract(self.fetch_details(self.accession_number_of_first_variant)) 
        
        if self.seq_to_search is None:
            return None 
            
        return self.seq_to_search

    def search_for_PAM(self):
        PAM = self.input_PAM_seq
        
        # Searching for PAM from left of exon sequence
        CDS_exon_seq = self.retrieve_seq_to_search_for_PAM()
        
        if isinstance(self.seq_to_search, Bio.SeqRecord.SeqRecord):
            self.seq_to_search = CDS_exon_seq
            self.start_index = 0
            
            while True:
                self.start_pos = str(self.seq_to_search.seq[self.start_index:]).find(PAM)
                self.start_direction_of_PAM_search = "left_to_right"
                self.target_seq = self.find_target_seq()
                
                if self.start_pos == -1:
                    break
              
                if len(self.target_seq) == self.input_length_of_target_seq:
                    self.start_pos = self.start_pos + self.start_index
                    self.target_seq = self.find_target_seq()
                    break
                    
                self.start_index = self.start_index + self.start_pos + len(PAM)            
                
            # If PAM not found, search for reverse complement of PAM from left of exon sequence
            if self.start_pos == -1 or len(self.find_target_seq()) != self.input_length_of_target_seq:
                self.start_index = 0
                complement_PAM = Seq(PAM).complement()
                reversed_complement_PAM = str(complement_PAM[::-1])
                self.prev_pos_of_PAM = self.start_index

                while True:
                    self.start_pos = str(self.seq_to_search.seq[self.start_index:]).find(reversed_complement_PAM)
                    self.start_direction_of_PAM_search = "right_to_left"
                    self.target_seq = self.find_target_seq()
                    
                    if self.start_pos == -1:
                        break
                        
                    if len(self.target_seq) == self.input_length_of_target_seq:
                        self.start_pos = self.start_pos + self.start_index
                        self.target_seq = self.find_target_seq()
                        break
                        
                    else:
                        self.start_index = self.start_index + self.start_pos + len(reversed_complement_PAM)

                # Search next exon (if exists) if reverse complement of PAM not found
                if self.start_pos == -1 or len(self.find_target_seq()) != self.input_length_of_target_seq:
                    self.start_index = 0
                    self.prev_pos_of_PAM = self.start_pos
                    next_exon_location = self.find_range_of_seq_to_search_past_CDS()

                    # Extracting the next exon sequence
                    if next_exon_location != None:
                        next_exon_seq_to_extract = SeqFeature(next_exon_location)
                        self.seq_to_search = next_exon_seq_to_extract.extract(self.record)
                        
                    else:
                        return None

                    # Searching for PAM from left of next exon sequence
                    while True:
                        self.start_pos = str(self.seq_to_search.seq[self.start_index:]).find(PAM)
                        self.start_direction_of_PAM_search = "left_to_right"
                        self.target_seq = self.find_target_seq()
                        
                        if self.start_pos == -1:
                            break
                            
                        if len(self.target_seq) == self.input_length_of_target_seq:
                            self.start_pos = self.start_pos + self.start_index
                            self.target_seq = self.find_target_seq()
                            break
                                
                        self.start_index = self.start_index + self.start_pos + len(PAM)

                    # Search this same exon for reverse complement of PAM from left of exon sequence
                    if self.start_pos == -1 or len(self.find_target_seq()) != self.input_length_of_target_seq:
                        self.start_index = 0
                        
                        while True:
                            self.start_index = self.start_index + self.start_pos + len(reversed_complement_PAM)
                            self.prev_pos_of_PAM = self.start_index
                            self.start_pos = str(self.seq_to_search.seq[self.start_index:]).find(reversed_complement_PAM)
                            self.start_direction_of_PAM_search = "right_to_left"
    
                            self.target_seq = self.find_target_seq()
                                
                            if self.start_pos == -1:
                                break
                            
                            if len(self.target_seq) == self.input_length_of_target_seq:
                                self.start_pos = self.start_pos + self.start_index
                                break  
                      
                        if self.start_pos == -1 or len(self.target_seq) != self.input_length_of_target_seq:
                            return None
            
            return self.seq_to_search
        
        return None          

    # Cleave DNA from end/beginning of PAM to a given number of nucleotides past/before the PAM
    def find_target_seq(self):
        if self.start_pos != -1 and isinstance(self.seq_to_search, Bio.SeqRecord.SeqRecord):
            
            if self.start_direction_of_PAM_search == "left_to_right":
                target_seq_DNA = self.seq_to_search.seq[(self.start_pos + len(self.input_PAM_seq)) : (self.start_pos + len(self.input_PAM_seq) + self.input_length_of_target_seq)]
                return target_seq_DNA
            
            elif self.start_direction_of_PAM_search == "right_to_left":
                start_of_target_seq = (self.prev_pos_of_PAM + self.start_pos) - self.input_length_of_target_seq
                target_seq_DNA = (self.seq_to_search.seq[(start_of_target_seq) : (self.prev_pos_of_PAM + self.start_pos)])[::-1]
                return target_seq_DNA 
        
        return None
    
    # Transcribe the target DNA sequence to RNA
    def transcribe_to_RNA(self):
        if isinstance(self.target_seq, Bio.Seq.Seq):
            target_seq_RNA = ((self.target_seq).transcribe())
            return target_seq_RNA
        
        return None
    
    # Prepend the repeat sequence to the target RNA sequence, thereby creating gRNA 
    def add_repeat_seq(self):
        target_RNA_seq = self.transcribe_to_RNA()
        if target_RNA_seq != None:
            gRNA = self.input_repeat_seq + target_RNA_seq
            return gRNA
        
        return None
    
    # Return gRNA sequence
    def return_gRNA_seq(self):
        if self.add_repeat_seq() != None:
            return self.add_repeat_seq()
        
        return None
    
    # Return length of gRNA sequence
    def return_length_of_gRNA_seq(self):
        if self.return_gRNA_seq() != None:
            return len(self.return_gRNA_seq())
        
        return None


# In[4]:


def main():
    current_date_time = datetime.datetime.now()
    formatted_curr_date_time = current_date_time.strftime("%Y-%m-%d_%Hh-%Mm-%Ss")

    while True:
        input_CSV = input("Enter filepath of CSV file with genes to search: ")

        if os.path.exists(input_CSV):
            break
        else:
            print("Invalid filepath. Please try again.")
            
    output_valid_CSV_file_path = input("Enter desired filepath to save valid genes and gRNA sequences CSV: ") + formatted_curr_date_time + ".csv"
   
    output_invalid_CSV_file_path = input("Enter desired filepath to save invalid gene CSV: ") + formatted_curr_date_time + ".csv"

    while True:
        input_PAM_seq = input("Enter PAM to search for: ")

        allowed_letters = set("ATGC")
        uppercase_input = input_PAM_seq.upper()

        if set(uppercase_input) <= allowed_letters:
            break
        else:
            print("The input contains characters other than A, T, G, and C. Please try again.")
    
    while True:
        input_repeat_seq = input("Enter repeat sequence: ")

        allowed_letters = set("AUGC")
        uppercase_input = input_repeat_seq.upper()

        if set(uppercase_input) <= allowed_letters:
            break
        else:
            print("The input contains characters other than A, U, G, and C. Please try again.")
    
    while True:
        input_length_of_target_seq = input("Enter length of target sequence: ")

        try:
            input_length_of_target_seq = int(input_length_of_target_seq)
            break
        except ValueError:
            print("Please enter an integer.")
      
    while True:
        input_species = input("Enter species of collection of genes: ")
        target_gene_library = TargetSeqLibrary(input_CSV, input_species, input_PAM_seq, input_repeat_seq, input_length_of_target_seq, output_valid_CSV_file_path, output_invalid_CSV_file_path)

        try:
            target_gene_library.check_every_gene_for_target_seq()
            break
            
        except requests.exceptions.HTTPError as http_err:
            print("Invalid species. Please check all valid species in NCBI taxonomy database.")  
            
    target_gene_library = TargetSeqLibrary(input_CSV, input_species, input_PAM_seq, input_repeat_seq, input_length_of_target_seq, output_valid_CSV_file_path, output_invalid_CSV_file_path)
    target_gene_library.check_every_gene_for_target_seq()
 
if __name__ == "__main__":
    main()

    

