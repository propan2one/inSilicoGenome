from insilicogenome import __version__
from insilicogenome import insilicogenome
from insilicogenome import insilicodata
from Bio.Seq import Seq
import numpy as np

def test_version():
    assert __version__ == '0.1.0'

# def test_random_dnasequence():
#     size = 10
#     assert ((insilicogenome.random_dnasequence(size)) == 'TGCTTGATGG')
# 
# def test_replace_start_codons():
#     assert ((insilicogenome.replace_start_codons("ATGATGAAATTCCTTGCATTT")) == 'GCTTCTAAATTCCGAAGATTT')
# 
# 
# # Unit test
# FASTA="miniasm.fasta"
# range_start=3
# range_end=36
# >>> random_snv('miniasm.fasta', 3, 36)
# ['utg000001c', 5, '.', 'C', 'A', 45, 'PASS']