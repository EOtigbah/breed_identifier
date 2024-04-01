import math
import os
import datetime
from datetime import timezone
import re

import matplotlib.pyplot as plt
import pandas as pd


from Bio import  Phylo, SeqIO
from Bio.Align import MultipleSeqAlignment, PairwiseAligner, substitution_matrices
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor


from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas



"""
Prompt user to input the breed file and the query sequence file paths, 

Need to hashtage this out to get my pytests to run 

"""
def define_breed_file ():
    breed_file = input("Please enter fasta file path of the the breed database: ")
    try:
        if not breed_file.endswith('.fa') and not breed_file.endswith('.fasta'): # Checking if the input file is a FASTA file 
            raise ValueError("Input file must be a FASTA file, please upload a FASTA file (.fa or .fasta)") # Raise an error if input is not a Fasta file 
        return breed_file
    except ValueError as error:
        print(error)


def define_query_file ():
    query_file = input("Please enter fasta file path of the mystery breed: ")
    try:
        if not query_file.endswith('.fa') and not query_file.endswith('.fasta'): # Checking if the input file is a FASTA file 
            raise ValueError("Input file must be a FASTA file, please upload a FASTA file (.fa or .fasta)") # Raise an error if input is not a Fasta file 
        return query_file
    except ValueError as error:
        print(error)
breed_file = define_breed_file()
query_file = define_query_file()    

""" 
Use this to run pytests (instead of code above)
Please set your file paths manually
"""
# breed_file = '/Users/chalupa/Documents/GitHub/breed_identifier/dog_breeds.fa' <- replace with your own file path for pytest
# query_file = '/Users/chalupa/Documents/GitHub/breed_identifier/mystery.fa' <- replace with your own file path for pytest


def database_processing(breed_file):
    """ 
    This function processes the dog_breeds fasta file (the database used to to identify the breed 
    Here the gene_ID, breed and sequence will be extracted
    """


    # Expression using the re function to find a set of strings that match it: this tells python to match this "[breed="
    # The bracets () tells python what we want to retrive, '.' matches any character, '*' mathces any sequence of characters after 
    # '?' tells python to lazy match ie match as little characters as possible 
    breed_search = re.compile(r'\[breed=(.*?)\]') # ie pull"[breed =" and "]""

    # Initialize list to store breeds later
    breed_database = []

    with open(breed_file, 'r') as file:
        for record in SeqIO.parse(file, 'fasta'): #Open fasta file
            breed_retrive = breed_search.search(record.description)  #Search the description in each sequence and get the breed
            if breed_retrive: # If statment to check whether it can find the pattern 
                breed = breed_retrive.group(1) # If it does then we capture the second part of out search (works as a index) i
            record.annotations['breed'] = breed # This is telling python for the extracted breeds we what to put them as an annotation in our record 
            breed_database.append(record) # Appends the data, SeIO.parse automatically pulls thre record.ID and seq to append, (by setting breed to record annotations we also append the breed to the list )

    
        if not breed_database:  # Check if the breed database is empty and if it is raise an error 
            raise ValueError("There is no data in the FASTA file")

        return breed_database





def query_processing (query_file):
    """ 
    This function processes the query sequence ready for alignment 
    """

    # Open mystery seqeunce file -> I have set this to my query seqeunce 
    query_sequence = SeqIO.read(query_file, 'fasta')

    # To determine the breed of mystery sequence I will set the breed to of the mystery sequence to 'Mystery' for now
    query_sequence.annotations['breed'] = 'mystery'

    # I want to combine the query_sequence with the breed_database so I can perform a MSA need to make a list containing the query_sequence 
    query_sequence_data = [query_sequence]

    return query_sequence_data


def combining_data(query_sequence_data, breed_database):
    """ 
    Function to combine the breed database and the query seqeunce together 
    """

    # Combine the lists together using '+'
    combined_data = breed_database + query_sequence_data

    return combined_data

def MSA(combined_data):
    """ 
    Perform a MSA on the combined data 
    """
    
    # Perform MSA using the MultipleSeqAlignment function from biopython
    alignment = MultipleSeqAlignment(combined_data)

    return alignment



def calculate_alignment_length (alignment):
    """ 
    Function to get length of the alignment
    """
    # Get the length of alignment  
    alignment_length = alignment.get_alignment_length() 
    return(alignment_length)


def calculate_identity_percentage (alignment, query_sequence, alignment_length): # Set our inputs 
    """ 
    Function to calculate the identity score
    query seqeunce is used as our comparative seqeunce
    identity score formula = (match base count / alignment length) * 100  
    """
    # Initialize list to store alignment score 
    identity_scores = []
    
    # For loop --> to loop through the alignment 
    for record in alignment:
        # Retriving matching nucleotides (count the number of matching nucleotides in the alignment) 
        # for i in range means iterate over each postion in our alignment sequence
        match_count = sum(1 for i in range(alignment_length) if record[i] == query_sequence[i]) # if it matches the postion in the query_seqeunce then count it as 1 

        # Identity score calculation 
        identity_percentage = (match_count / alignment_length) * 100 

        
        breed = record.annotations.get('breed', 'Unknown') # Set breed to unknown if not present

        # Append the record.id, record breed and identity percentage to a list)
        identity_scores.append((record.id, [breed], identity_percentage))


    return identity_scores


def convert_data_df_output (alignment, identity_scores):
    """ 
    This funciton converts our lists to a pandas dataframe
    """

    # Initialize list, assemble our list contining the record_ID, breed, identity percentage and sequence  
    data = []

    # Set data catagories we want to pull and where we are pullign them from 
    for record, (record_id, breed, identity) in zip(alignment, identity_scores):
        record_data = {
            "record_id": str(record_id), 
            "breed": breed,
            "identity_(%)": identity,
            "sequence": str(record.seq) 
        }
        data.append(record_data)
    
    # Converts data to pandas dataframe 
    breed_match_dataframe = pd.DataFrame(data)

    return breed_match_dataframe


def convert_data_df_csv (alignment, identity_scores):
    """ 
    This funciton converts our lists to a pandas dataframe
    """

    # Initialize list, assemble our list contining the record_ID, breed, identity percentage and sequence 
    data = []

    # Set data catagories we want to pull and where we are pullign them from 
    for record, (record_id, breed, identity) in zip(alignment, identity_scores):
        record_data = {
        "record_id": str(record_id), 
        "breed": breed,
        "identity_(%)": identity,
        }
        data.append(record_data)
    
    # Converts data to pandas dataframe 
    breed_match_dataframe_csv = pd.DataFrame(data)

    return breed_match_dataframe_csv



def retrieve_match (breed_match_dataframe): 
    """ 
    This fucntion retrieves the breed match 
    Using the pandas df, we can filter the df agaisnt identity in descending order to get the Max identiry score
    Need to be careful since we loaded our query seqeunce into the alignment file as a control it will have a 100% match to itself 
    ---> to get the breed match the second row needs to be pulled 
    """

    # Sort values in descending order using the indentity percentage column 
    sorted_data = breed_match_dataframe.sort_values(by='identity_(%)', ascending=False)
    #print(sorted_data)

    #Use iloc to slice df and pull the second row
    breed_match = sorted_data.iloc[1]

    return breed_match


def get_matched_breed_sequence_str(mystery_breed_solved, breed_database):
    """" 
    This function pulls the breed match seqeunce
    This will be used later for e-value and p-value calculations
    Here the breed match seqeunce is acquired by its record_ID as this is the seqeucnes unique identifier
    """

    record_id = mystery_breed_solved.get('record_id')

    if record_id is None:
        print("No record ID found in mystery_breed_solved")
        return None

    # Find the matching record in the breed database
    matched_breed_record = None
    for record in breed_database:
        if record.id == record_id:
            matched_breed_record = record
            break

    # If the matching record is found, get its sequence
    if matched_breed_record:
        matched_breed_sequence = str(matched_breed_record.seq)
        return matched_breed_sequence
    # If not print error message 
    else:
        print("Error: matched breed record not found in the database.")
        return None


def get_query_sequence_str (query_sequence_data):
    query_sequence_str = ""
    for record in query_sequence_data:
        query_sequence_str += str(record.seq)
        return query_sequence_str
    






def get_differences (query_sequence_str, matched_breed_sequence_str):
    """ 
    This function gets the difference between the squery sequence and the breed matched seqeunce 
    """
    # Perform pairwise alignment
    aligner = PairwiseAligner() # Set aligner 
    aligner.mode = 'global' # Set mode

    match_alignment = aligner.align(query_sequence_str, matched_breed_sequence_str)

    # Set which each to their resprective index in the pairwaise alignment 
    aligned_query_sequence = match_alignment[0][0] 
    aligned_matched_sequence = match_alignment[0][1]

    # Initialize a list to store the differences
    differences = []

    # Compare the sequences nucleotide by nucleotide
    for i in range(len(aligned_query_sequence)): # Use a for loop --> loop over each postion of the alignment query seqeunce (use the length as our range)

        # If statement to get differences using "!=" not equal expression 
        if aligned_query_sequence[i] != aligned_matched_sequence[i]:
            # If the nucleotide in the query seqeunce does not equal the nucleotide in the matched seqeunce ---> append the: position(i), the nucleotide at postion i in the query and matched sequqnce  
            differences.append((i, aligned_query_sequence[i], aligned_matched_sequence[i])) 

            
    differences_df = pd.DataFrame(differences, columns=["position", "query_sequence", "matched_sequence"])
  
    return differences_df

def calculate_alignment_score(sequence1, sequence2):
    """ 
    Constructed a nucleotide substitution matrix
    This will be used to score the alignment match 
    The logic of this is 

    2 points if the nucleotides in the query and matched sequence match 
    -3 points if the mutation is a transtion mutation (this is a more likely mutation in terms of evolution)
    -4 points if the mutation is a tranversion mutation (this is less likely than the transation mutation) 
    -2 points if the posion is a space in one seqeunce with a nucleotide in the other seqeuence 
    0 points if the match is just a gap in both seqeunces 
    """
    

    # Set mutation scoring values
    nucleotide_match = 2
    nucleotide_mismatch_transition = -3
    nucleotide_mismatch_transversions = -4
    nucleotide_mismatch_gap = -5
    gap_gap = 0 

    # Make a dictionary with mutation scoring values to match up to mutation key 
    nucleotide_scoring_dict = {
        ('A', 'A'): nucleotide_match,
        ('C', 'C'): nucleotide_match,
        ('G', 'G'): nucleotide_match,
        ('T', 'T'): nucleotide_match,
        ('-', '-'): gap_gap,
        ('A', 'G'): nucleotide_mismatch_transition,
        ('C', 'T'): nucleotide_mismatch_transition,
        ('G', 'A'): nucleotide_mismatch_transition,
        ('T', 'C'): nucleotide_mismatch_transition,
        ('A', 'C'): nucleotide_mismatch_transversions,
        ('A', 'T'): nucleotide_mismatch_transversions,
        ('C', 'G'): nucleotide_mismatch_transversions,
        ('G', 'T'): nucleotide_mismatch_transversions,
        ('C', 'A'): nucleotide_mismatch_transversions,
        ('T', 'A'): nucleotide_mismatch_transversions,
        ('G', 'C'): nucleotide_mismatch_transversions,
        ('T', 'G'): nucleotide_mismatch_transversions,
        ('-', 'A'): nucleotide_mismatch_gap,
        ('-', 'C'): nucleotide_mismatch_gap,
        ('-', 'G'): nucleotide_mismatch_gap,
        ('-', 'T'): nucleotide_mismatch_gap,
        ('A', '-'): nucleotide_mismatch_gap,
        ('C', '-'): nucleotide_mismatch_gap,
        ('G', '-'): nucleotide_mismatch_gap,
        ('T', '-'): nucleotide_mismatch_gap,
    }

    # Make an array of the subsitution matrix using the substitution_matrices of biopython 
    nucleotide_scoring_matrix = substitution_matrices.Array(None, data=nucleotide_scoring_dict)

    # Calculate alignment score using the substitution matrix 
    alignment_score = 0 # Set counter to zero 

    # for loop to iterate over each nucleotide in each seqeunce 
    for nucleotide_1, nucleotide_2 in zip(sequence1, sequence2):

        # If statement to check status of each query seqeunce nucleotide with its corresponding nucleotide in the matched sequence 
        if (nucleotide_1, nucleotide_2) in nucleotide_scoring_matrix: # If the match, mismatch is in the matrix 
            alignment_score += nucleotide_scoring_matrix[(nucleotide_1, nucleotide_2)] # Add its respective score to our counter 

    
    return alignment_score

def calculate_bit_score (alignment_score):
    """
    Function to calculate bit score 
    bit score formula: (lambda_value * alignment_score - math.log(K_value)) / math.log(2) 
    bit formula, E-value and P-value formulas obtained from https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html#head3

    https://www.ncbi.nlm.nihgov/BLAST/tutorial/Altschul-3.html (for lambda and k values)
    lambda_value = 0.252
    K_value = 0.032

    """
    lambda_value = 0.252
    K_value = 0.035
    bit_score = (lambda_value * alignment_score - math.log(K_value)) / math.log(2)
    return bit_score

def calculate_e_value (query_sequence_str, bit_score, alignment):
    """
    Function to calculate e-value
    P-value calculation formula = 1 - math.exp(-e_value)
    """
    m = (len(alignment)-1)
    n = (len(query_sequence_str))
    e_value = m * n * (2 **(- bit_score))
    
    return e_value

def calculate_p_value (e_value):
    """
    Function to calculate p-value  
    P-value calculation formula = 1 - math.exp(-e_value)
    """
    p_value = 1 - math.exp(-e_value)
    return p_value



breed_database = database_processing(breed_file)

query_sequence_data = query_processing(query_file)

combined_data = combining_data(query_sequence_data, breed_database)

alignment = MSA(combined_data)


alignment_length = calculate_alignment_length(alignment)

identity_scores = calculate_identity_percentage(alignment, query_sequence_data[0], alignment_length)

breed_match_dataframe = convert_data_df_output(alignment, identity_scores)
breed_match_csv = convert_data_df_csv(alignment, identity_scores)

mystery_breed_solved = retrieve_match(breed_match_dataframe)
print(mystery_breed_solved)

matched_breed_sequence_str = get_matched_breed_sequence_str(mystery_breed_solved, breed_database)

query_sequence_str = get_query_sequence_str(query_sequence_data)



differences_df = get_differences(query_sequence_str, matched_breed_sequence_str)


alignment_score = calculate_alignment_score(query_sequence_str, matched_breed_sequence_str)
bit_score = calculate_bit_score(alignment_score)
e_value = calculate_e_value(query_sequence_str, bit_score, alignment)
p_value = calculate_p_value (e_value)
print(alignment_score,bit_score,e_value,p_value)



def phylo_tree_builder (alignment):
    """ Function to build the phylogenetic_tree builder 
    This was hard to figure out ):
    """
    calculator = DistanceCalculator('identity')
    dist_matrix = calculator.get_distance(alignment)

   
    constructor = DistanceTreeConstructor()
    UPGMATree = constructor.upgma(dist_matrix)
   
    scaling_factor = 2.0
    for clade in UPGMATree.find_clades():
        if clade.branch_length:
            clade.branch_length *= scaling_factor


    # Need to make a dictionary so I can match the recordID to the breed
    record_id_to_breed_dictionary = {}

    # Iterate over the breed_database and get the record_id and match it to the breed --> store it in the dictionary
    for record in alignment:
        record_id_to_breed_dictionary[record.id] = record.annotations['breed']


    fig, ax = plt.subplots(figsize=(60, 400))  # Used to adjust the figure size since there is a lot and wont fit into the standard size
    Phylo.draw(UPGMATree, axes=ax, do_show=False, label_func=lambda x: record_id_to_breed_dictionary.get(x.name, '')) # Name the clade branches according to data in the record_id_bred_dict

    return plt


"""
Finally we save all the data
1) Make a results folder (if it doesn't already exist)
2) Save each result as the time and the date that is was ran
3) Save the fasta file, a PDF for the results, a CSV file for the raw results and a Phylogenetic tree
"""

# Define the path where I am going to save my results to 
results_folder = 'results'

# I have decided that I will name each job as its own date and time 
job_date_time = job_date_time = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
folder_name = os.path.join(results_folder, f'Job_{job_date_time}')
os.makedirs(folder_name, exist_ok=True)  #I used exist_ok = True in the case that the folder exists this will avoid an error being raised

# Save fasta file to folder 
fasta_output_file= os.path.join(folder_name, 'Multiple_seqeunce_alignment.fasta')
SeqIO.write(alignment, fasta_output_file, "fasta")  

print(f'Multiple sequence alignment saved as a fasta file to: {fasta_output_file}')

# Saving to a PDF file 
pdf_output_file = os.path.join(folder_name, 'analysis_summary.pdf')
with open(pdf_output_file, 'wb') as f:

    # Define canvas size
    c = canvas.Canvas(f, pagesize=A4)

    # Set the font to bold
    c.setFont("Courier-Bold", 26)  

    c.drawString(100, 800, "Breed Match Results")
     
    c.setFont("Courier-Bold", 12)  

     
    c.setFont("Courier-Bold", 12)  
    c.drawString(100, 760, "Breed Match:")

     # Set it back 
    c.setFont("Courier", 12)
    c.drawString(100, 740, f'Record ID: {mystery_breed_solved['record_id']}')

    breed_print = mystery_breed_solved['breed'][0].strip("[]' ")
    c.drawString(100, 720, f'Breed: {breed_print}')
    c.drawString(100, 700, f'Identity Percentage (%): {mystery_breed_solved['identity_(%)']}')



    c.setFont("Courier-Bold", 12)  
    c.drawString(100, 660, "Statistics:")

   
    c.setFont("Courier", 12)


    c.drawString(100, 640, f'Alignment Score: {alignment_score}')
    c.drawString(100, 620, f'Bit Score: {bit_score}')
    c.drawString(100, 600, f'E-value: {e_value}')
    c.drawString(100, 580, f'P-value: {p_value}')

    # Set the font to bold
    c.setFont("Courier-Bold", 12)  
    c.drawString(100, 540, "Mismatches:")

    # Set it back 
    c.setFont("Courier", 12)

    # Convert differences_df to a string
    differences_df_str = differences_df.to_string(index=False)

    # Construct table for mismatches 
    x= 100
    y = 520

    # Split each line into individual lines
    lines = differences_df_str.split('\n')

    # Using a for loop -> iterate over each line and add it to the table 
    for line in lines:
        c.drawString(x, y, line)
        y -= 10  # This adjusts spacing between each line 

    c.save()
# Save results to a PDF file 
print(f'Summary of results saved to a PDF in: {pdf_output_file}')

# Save the dataframe of record_IDs breeds and identity percentages to a csv file 
csv_output_file = os.path.join(folder_name, 'analysis_csv.csv')
breed_match_csv.to_csv(csv_output_file, sep=',', index=False, encoding='utf-8')

print(f'Database containing record_ID, breeds and identity percentages saved to a CSV file: {csv_output_file}')


# Save phylogenetic tree output as a jpg
graph_output_file = os.path.join(folder_name, 'phylogenetic_tree.png')
phylo_tree_builder(alignment).savefig(graph_output_file)

print(f'Phylogenetic tree saved as png to: {graph_output_file}')
