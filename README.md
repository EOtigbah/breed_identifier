# breed-identifier
A program that allows the the user to determine the breed of their dog using sequencing data 


To run the program you will be prompted to input the breed database (dog_breeds.fa) and the mystery breed (mystery.fa) please input the file paths of these files when prompted 

Please makesure you are in a known directory as a resuls folder will be created with the analysis of your query 

In order to run this program you will need to have certain modules on python
Below is a comprehesive list of modules:

### Modules Required
**Biopyhton:**
Phylo
SeqIO


**Bio.Align:**
MultipleSeqAlignment
PairwiseAligner
Substitution_matrices


**Bio.Phylo.TreeConstruction:**
DistanceCalculator
DistanceTreeConstructor 

**Bio.SeqRecord:**
SeqRecord

**Bio.Seq:**
Seq

**pandas**

**re**

**matplotlib**
**math**
**os**

**datetime**
timezone

**reportlab.lib.pagesizes**
A4

**reportlab.pdfgen**
canvas


Or you can copy this code into your console:

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







