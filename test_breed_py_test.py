
import os


# Define the target directory
target_directory = '/Users/chalupa/Documents/GitHub/breed_identifier/'

# Change the current directory to the target directory
os.chdir(target_directory)
current_directory = os.getcwd()
print("Current Directory:", current_directory)

# breed_file = '/Users/chalupa/Documents/GitHub/breed_identifier/dog_breeds.fa'

# # 
# query_file = '/Users/chalupa/Documents/GitHub/breed_identifier/mystery.fa'

import pytest
from breed_identifier import database_processing    
def test_database_processing():
    #Call the function in the test file 
    breed_file = './test_rottweiler.fa'
    breed_database = database_processing(breed_file)
    
    #Check that the breed was extracted 
    assert breed_database[0].annotations['breed'], 'Rottweiler cross'
test_database_processing()

# import unittest
# from breed_identifier_ import database_processing

# # class TestDatabaseProcessing(unittest.TestCase):
# #     def test_database_processing(self):
# #         # Call the function in the test file
# #         breed_file = './test_rottweiler.fa'
# #         breed_database = database_processing(breed_file)
        
# #         # Check that the breed was extracted
# #         self.assertEqual(breed_database[0].annotations['breed'], 'Rottweiler cross')

# # if __name__ == '__main__':
# #     unittest.main()

# from breed_identifier import query_processing
# def test_query_processing():
#     #Call the function in the test file 
#     query_file = './test_rottweiler.fa'
#     breed_database = query_processing (query_file)
    
#     #Check that the breed was extracted 
#     assert breed_database[0].annotations['breed'], 'mystery'
# test_query_processing()

# from breed_identifier import calculate_alignment_length 

# def test_calculate_alignment_length ():

#     from Bio import AlignIO

#     # Load up a test alignment of a known length (it only contains 1 sequence)
#     alignment = AlignIO.read('/Users/chalupa/Documents/GitHub/breed_identifier/dog_breeds_test.fasta', 'fasta')

#     # Run function to test
#     alignment_length = calculate_alignment_length(alignment)

#     # I know the length because I wrote the seqeunces
#     expected_length = 8

#     assert alignment_length == expected_length

# test_calculate_alignment_length()

# from breed_identifier import calculate_identity_percentage

# def test_calculate_identity_percentage():

#     from Bio import AlignIO

#     # Test alignment 
#     alignment = AlignIO.read('/Users/chalupa/Documents/GitHub/breed_identifier/dog_breeds_test.fasta', 'fasta')

#     # Alignment length set because I know it 
#     alignment_length = 8

#     # Set query_seqeunce to same seqeunce in alignment (to get 100% score)
#     query_sequence = 'GTTAATGT'
#     identity_score = calculate_identity_percentage(alignment, query_sequence, alignment_length)

#     expected_score = [('gb|CM023446.1|', ['Unknown'], 100.0)]

#     assert identity_score == expected_score

# test_calculate_identity_percentage()


# from breed_identifier import get_differences

# def test_get_differences ():

#     from Bio.Align import PairwiseAligner
#     import pandas as pd

#     query_sequence_str = 'ATGATG'
#     matched_breed_sequence_str = 'ATGATA'

#     differences_df = get_differences (query_sequence_str, matched_breed_sequence_str)
    
#     # Make a df of known results 
#     expected_differences = pd.DataFrame({"position": [5, 6],
#                                          "query_sequence": ["G", "-"],
#                                          "matched_sequence": ["-", "A"]
#                                          })

#     assert differences_df.equals(expected_differences) 


# test_get_differences ()

# from breed_identifier import calculate_alignment_score

# def test_calculate_alignment_score():
#     sequence1 = 'ATGATG'
#     sequence2 = 'ATGATA'

#     alignment_score = calculate_alignment_score (sequence1, sequence2)

#     expected_score = 7

#     assert alignment_score == expected_score 








