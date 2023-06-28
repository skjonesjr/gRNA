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
        for gene in self.query_accession_numbers()['out']:
            if 'exons' in gene:
                all_exon_ranges = gene['exons']
                gene_id = gene['query']
                accession_number_of_first_variant = (all_exon_ranges[0])['transcript']
                first_variant = SingleGeneSeq(gene_id, accession_number_of_first_variant, self.PAM, self.repeat_seq, self.length_target_seq)

                # Check if a valid target sequence was found and if the length of gRNA is the same as the inputted length of the target sequence plus the length of the inputted repeat sequence.
                # The second condition may occur if the PAM is found too close to the end of an exon, meaning that the target sequence will be shorter than the inputted length of the target sequence.
                if first_variant.search_for_PAM().seq != Seq('X') and first_variant.return_length_of_gRNA_seq() == (self.length_target_seq + len(self.repeat_seq)):
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
    def __init__(self, gene_id, accession_number_of_first_variant, input_PAM_seq, input_repeat_seq, input_length_of_target_seq):
        self.gene_id = gene_id
        self.accession_number_of_first_variant = accession_number_of_first_variant
        self.input_PAM_seq = input_PAM_seq
        self.input_repeat_seq = input_repeat_seq
        self.input_length_of_target_seq = input_length_of_target_seq
        
        self.record = None
        self.start_pos = None
        self.seq_to_search = None
        self.range_of_first_exon_with_start_of_CDS = None
        self.target_seq_length = 0
        
    # Retrieve details from NCBI based on a gene's first variant. These details include the exon and coding sequence range specific to the first variant.
    def fetch_details(self, accession_number_of_first_variant):
        handle = Entrez.efetch(db = "nucleotide", id = accession_number_of_first_variant, rettype = "gb", retmode = "text")
        details = SeqIO.read(handle, "genbank")
        handle.close()
        return details
 
    # Find first exon range which includes the beginning of the coding strand
    def find_range_of_seq_to_search_with_CDS(self):
        self.record = self.fetch_details(self.accession_number_of_first_variant)
        
        # Iterating through features of variant to find range of coding strand
        for feat in self.record.features:
            if feat.type == 'CDS':
                start_of_CDS_range_of_variant = feat.location.start
                break
                
        # Iterating through features of variant to find range of exon which contains beginning of coding strand
        for feat in self.record.features:
            if feat.type == 'exon' and feat.location.end > start_of_CDS_range_of_variant:
                exon_range_of_variant = feat.location
                break

        feature_location = FeatureLocation(start_of_CDS_range_of_variant, exon_range_of_variant.end)
    
        return feature_location
    
    def find_range_of_seq_to_search_past_CDS(self):
        self.record = self.fetch_details(self.accession_number_of_first_variant)
        
        for feat in self.record.features:
            next_exon_range_of_variant = None
            
            # Due to the transition from 1-based indexing of the NCBI database to 0-based indexing of Python for a FeatureLocation object, add 1 to the starting position of an exon to account for this when comparing two FeatureLocations.
            if feat.type == 'exon' and (feat.location.start + 1) > self.range_of_first_exon_with_start_of_CDS.end:
                next_exon_range_of_variant = feat.location
                break

        if next_exon_range_of_variant != None:
            next_exon_location = FeatureLocation(next_exon_range_of_variant.start, next_exon_range_of_variant.end)
            return next_exon_location
        else:
            return None
    
    # Extract the sequence of a specific range to look for PAM.
    def retrieve_seq_to_search_for_PAM(self):
        feat_location = self.find_range_of_seq_to_search_with_CDS()
        self.range_of_first_exon_with_start_of_CDS = feat_location

        # Extracting the selected exon sequence
        exon_seq_to_extract = SeqFeature(feat_location)
        self.seq_to_search = exon_seq_to_extract.extract(self.fetch_details(self.accession_number_of_first_variant))
    
        return self.seq_to_search
    
    # Look for PAM in extracted sequence. First, search for PAM from beginning of CDS to end of select exon.
    # If PAM not found, then, search for complement of PAM along the sequence.
    # If PAM still not found, search for PAM in the next exon.
    def search_for_PAM(self):
        PAM = self.input_PAM_seq
        
        # Searching for PAM from left of exon sequence
        CDS_exon_seq = self.retrieve_seq_to_search_for_PAM()
        self.start_pos = (CDS_exon_seq.seq).find(PAM)
        self.seq_to_search = CDS_exon_seq
        
        # Search for complement of PAM from right of exon sequence, if PAM not found
        if self.start_pos == -1:
            complement_PAM = Seq(PAM).complement()
            self.start_pos = (CDS_exon_seq.seq).rfind(complement_PAM)
            
        # Search next exon (if exists) if complement of PAM not found
            if self.start_pos == -1:
                next_exon_location = self.find_range_of_seq_to_search_past_CDS()
        
                if next_exon_location != None:
                    # Extracting the next exon sequence
                    next_exon_seq_to_extract = SeqFeature(next_exon_location)
                    self.seq_to_search = next_exon_seq_to_extract.extract(self.record)
                    # Searching for PAM from left of next exon sequence
                    self.start_pos = (self.seq_to_search.seq).find(PAM)

                if self.start_pos == -1:
                    self.seq_to_search.seq = Seq("X")
      
        return self.seq_to_search
    
    # Cleave DNA from end of PAM to a given number of nucleotides past the PAM
    def find_target_seq(self):
        if self.start_pos != None:
            target_seq_DNA = self.seq_to_search[(self.start_pos + len(self.input_PAM_seq)):(self.start_pos + len(self.input_PAM_seq) + self.input_length_of_target_seq)]

            return target_seq_DNA
    
    # Transcribe the target DNA sequence to RNA
    def transcribe_to_RNA(self):
        target_seq = self.find_target_seq()
        target_seq_RNA = (target_seq.seq).transcribe()
        
        return target_seq_RNA
    
    # Prepend the repeat sequence to the target RNA sequence, thereby creating gRNA 
    def add_repeat_seq(self):
        gRNA = self.input_repeat_seq + self.transcribe_to_RNA()
        
        return gRNA
    
    # Return gRNA sequence
    def return_gRNA_seq(self):
        return self.add_repeat_seq()
    
    # Return length of gRNA sequence
    def return_length_of_gRNA_seq(self):
        return len(self.return_gRNA_seq())


# In[4]:


def main():
    current_date_time = datetime.datetime.now()
    formatted_curr_date_time = current_date_time.strftime("%Y-%m-%d_%Hh-%Mm-%Ss")
    
    #input_CSV = r"C:\Users\ievut\OneDrive\Desktop\Oncogene Library Project\Some 30 Genes of Oncogene Library.csv" 
    input_CSV = input("Enter filepath of CSV file with genes to search: ")
    #output_valid_CSV_file_path = r"C:\Users\ievut\OneDrive\Desktop\Oncogene Library Project\valid_genes_" + formatted_curr_date_time + ".csv"
    output_valid_CSV_file_path = input("Enter desired filepath to save valid target sequence CSV: ") + formatted_curr_date_time + ".csv"
    #output_invalid_CSV_file_path = r"C:\Users\ievut\OneDrive\Desktop\Oncogene Library Project\invalid_genes_" + formatted_curr_date_time + ".csv"
    output_invalid_CSV_file_path = input("Enter desired filepath to save invalid gene CSV: ") + formatted_curr_date_time + ".csv"
    #input_PAM_seq = "TTTC"
    input_PAM_seq = input("Enter PAM to search for: ")
    #input_repeat_seq = 'UAAUUUCUACUAAGUGUAGAU'
    input_repeat_seq = input("Enter repeat sequence: ")
    #input_species = "human"
    input_species = input("Enter species of collection of genes: ")
    #input_length_of_target_seq = 20
    input_length_of_target_seq = int(input("Enter length of target sequence: "))
    
    target_gene_library = TargetSeqLibrary(input_CSV, input_species, input_PAM_seq, input_repeat_seq, input_length_of_target_seq, output_valid_CSV_file_path, output_invalid_CSV_file_path)
    target_gene_library.check_every_gene_for_target_seq()
    
if __name__ == "__main__":
    main()
    

