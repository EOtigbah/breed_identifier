# breed-identifier
A program that allows the the user to determine the breed of their dog using sequencing data 

In order to run this program you will need to have certain modules on python
Below is a comprehesive list of modules:

### Modules Required
**Biopyhton:**
- SeqIO
- AlignIO
- Phylo
- SeqIO

**Bio.Align:**
MultipleSeqAlignment

**Bio.Phylo.TreeConstruction:**
DistanceCalculator
DistanceTreeConstructor 

**Bio.SeqRecord:**
SeqRecord

**Bio.Seq:**
Seq

**scipy.stats:**
ttest_1samp

**pandas**

**re**

**matplotlib**


**networkx**

**pydot**

Or you can copy this code into your console:

from Bio import pairwise2, Phylo, SeqIO
from Bio.Align import MultipleSeqAlignment, substitution_matrices
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import SeqIO

import math
import os
from datetime import datetime
import re

import matplotlib.pyplot as plt
import pandas as pd




from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas





