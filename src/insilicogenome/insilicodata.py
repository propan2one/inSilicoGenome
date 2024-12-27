import os
import csv
import getopt
import sys
import argparse
import numpy as np
from random import seed
from random import randint
from random import choice
from itertools import groupby
#from src.insilicogenome import insilicogenome
from insilicogenome import insilicogenome

help_message = """
usage: insilicogenome -o|--output <string> [-s|--size] <integer>
A Python utility to create artificial genomes, this tool is designed to be used in bioinformatics benchmarking programs.
Arguments :
    -h, --help          show this help message and exit.
    -i  --input         Fasta file [string]
    -b  --range_start   Begining of genomic feature [int]
    -e  --range_end     End of genomic feature [int]
    -c, --csv           Path to a csv file containing modification. [string]
    -s, --size          Size of the insertion and deletion. [int]
 Examples
    --------
    python ~/JD/ressources_bioinformatiques/insilicodata.py -i haplo01.fa -b 602 -e 2309 #Creat TSV file with variation
    python ~/JD/ressources_bioinformatiques/insilicodata.py -i haplo01.fa -c pMRT10279.tsv
"""

try:
    options, args = getopt.getopt(sys.argv[1:], "ho:i:s:b:e:c:s", ["help", "input=", "range_start=", "range_end=","csv=","size="])
except getopt.GetoptError as err:
    print(str(err))
    print(help_message)
    sys.exit(2)

#Default parameters
range_start=5
range_end=500

for opt, arg in options:
    if opt in ("-h", "--help"):
        print(help_message)
        sys.exit()
    elif opt in ("-i", "--input"):
        input_fasta = arg # attention c'est pas du input only
    elif opt in ("-b", "--range_start"):
        range_start = int(arg)
    elif opt in ("-e", "--range_end"):
        range_end = int(arg)
    elif opt in ("-c", "--csv"):
        csv_file = arg
    elif opt in ("-s", "--size"):
        size = int(arg)
        assert size<5000 ,'Error please provide a size value lower than 2 kbp'
    else:
        print(help_message)
        sys.exit(2)

# def write_fasta_genome(output, sequence):

def fasta_iter(FASTA):
    """
    Fasta file parser
    Credit from modified from Brent Pedersen and
         Ram & niuyw (https://www.biostars.org/p/710/)

    Parameters
    ----------
    FASTA : The name of the header and fasta file
      Fasta file

    Returns
    -------
    A generator object containing header name, sequence and header

    Examples
    --------
    >>> from insilicogenome import insilicodata
    >>> fasta = insilicogenome.fasta_iter(/path/genome.fasta)
    >>> for ff in fasta:
    ...     name, sequence, long_name = ff
    ...     print(sequence)
    
    """
    fin = open(FASTA, 'rb')
    faiter = (x[1] for x in groupby(fin, lambda line: str(line, 'utf-8')[0] == ">"))
    for header in faiter:
        headerStr = str(header.__next__(), 'utf-8')
        long_name = headerStr.strip().replace('>', '')
        name = long_name.split()[0]
        sequence = "".join(str(s, 'utf-8').strip() for s in faiter.__next__())
        yield (name, sequence, long_name)

def random_snv(FASTA, range_start=range_start,
            range_end=range_end, QUAL=45, FILTER='PASS', INFO=''):
    """
    Creat a list which contain information on SNP/SNP
    Optimize to parse a GFF3 file and create a VCF file in silico
    Parameters
    ----------
    FASTA : The name of the header and fasta file
      Fasta file, with the path (optionnal)
    
    Returns
    -------
    A list with: #CHROM POS ID REF ALT QUAL FILTER
    Examples
    --------
    >>> from insilicogenome import insilicodata
    >>> insilicodata.random_snv(/path/genome.fasta, range_start, range_end)
    >>> insilicodata.random_snv('miniasm.fasta', 3, 36)
    """
    if ((range_start>= 0) and (range_end>0) and (range_end>=range_start)):
        pass
    else:
        raise ValueError(
        "Problems with the range where the SNV occurs. "
        f"The value of range_start is '{range_start}'"
        f"The value of range_end is '{range_end}'")
    # add controle pour s'assurer que les valeur de range sont pas or range de seq
    if not os.path.isfile(FASTA):
        raise FileExistsError
    else:
        pass
    #fasta = insilicogenome.fasta_iter(FASTA)
    nucleotides=['A', 'T', 'C', 'G']
    #np.random.seed(1)
    fasta = fasta_iter(FASTA)
    for ff in fasta:
        name, sequence, long_name = ff
    row = []
    row.append(name) # CHROM
    if range_start == range_end:
        # This change allow multiple SNV
        POS=range_start
    else:
        POS=np.random.choice(range(range_start, range_end))
    row.append(POS) # POS
    row.append('.') # ID
    REF=sequence[POS-1] ### changing POS
    row.append(REF) # REF # Add strand option from GGF3
    row.append(''.join(np.random.choice([x for x in nucleotides if x not in REF], 1))) # ALT
    row.append(QUAL) # QUAL
    row.append(FILTER) # FILTER
    row.append(INFO) # INFO
    return row

def random_mnv(FASTA, range_start=range_start,
            range_end=range_end, QUAL=45, FILTER='PASS', INFO='', size=5):
    """
    Creat a list which contain information on SNP/SNP
    Optimize to parse a GFF3 file and create a VCF file in silico
    Parameters
    ----------
    FASTA : The name of the header and fasta file
      Fasta file, with the path (optionnal)
    
    Returns
    -------
    A list with: #CHROM POS ID REF ALT QUAL FILTER
    Examples
    --------
    >>> from insilicogenome import insilicodata
    >>> insilicodata.random_mnv(/path/genome.fasta, range_start, range_end)
    >>> insilicodata.random_mnv('miniasm.fasta', 3, 36)
    """
    if ((range_start>= 0) and (range_end>0) and (range_end>=range_start)):
        pass
    else:
        raise ValueError(
        "Problems with the range where the SNV occurs. "
        f"The value of range_start is '{range_start}'"
        f"The value of range_end is '{range_end}'")
    # add controle pour s'assurer que les valeur de range sont pas or range de seq
    if not os.path.isfile(FASTA):
        raise FileExistsError
    else:
        pass
    #fasta = insilicogenome.fasta_iter(FASTA)
    nucleotides=['A', 'T', 'C', 'G']
    #np.random.seed(1)
    fasta = fasta_iter(FASTA)
    for ff in fasta:
        name, sequence, long_name = ff
    row = []
    row.append(name) # CHROM
    if range_start == range_end:
        # This change allow multiple SNV
        POS=range_start
    else:
        POS=np.random.choice(range(range_start, range_end))    
    row.append(POS) # POS
    row.append('.') # ID
    REF='' # REF # Add strand option from GGF3
    for i in range(0, size):
        REF= REF + sequence[POS-1+i]
    row.append(REF) 
    ALT=''  # ALT
    for i in range(0, size):
        ALT= ALT + ''.join(np.random.choice([x for x in nucleotides if x not in REF[i]], 1))
    row.append(ALT)
    row.append(QUAL) # QUAL
    row.append(FILTER) # FILTER
    row.append(INFO) # INFO
    return row

def random_ins(FASTA, range_start=range_start,
            range_end=range_end, QUAL=45, FILTER='PASS', INFO='', size=5):
    """
    Creat a list which contain information on SNP/SNP
    Optimize to parse a GFF3 file and create a VCF file in silico
    Parameters
    ----------
    FASTA : The name of the header and fasta file
      Fasta file, with the path (optionnal)
    range_start : position where the range of insertion start
        int
    range_end : position where the range of insertion end
        int
    size : the size of the insertion
        int
    
    Returns
    -------
    A list with: #CHROM POS ID REF ALT QUAL FILTER
    Examples
    --------
    >>> from insilicogenome import insilicodata
    >>> insilicodata.random_ins(/path/genome.fasta, range_start, range_end, size)
    >>> insilicodata.random_ins('miniasm.fasta', 3, 36, 5)
    """
    if ((range_start>= 0) and (range_end>0) and (range_end>=range_start)):
        pass
    else:
        raise ValueError(
        "Problems with the range where the SNV occurs. "
        f"The value of range_start is '{range_start}'"
        f"The value of range_end is '{range_end}'")
    # add controle pour s'assurer que les valeur de range sont pas or range de seq
    if not os.path.isfile(FASTA):
        raise FileExistsError
    else:
        pass
    #fasta = insilicogenome.fasta_iter(FASTA)
    #np.random.seed(2)
    fasta = fasta_iter(FASTA)
    for ff in fasta:
        name, sequence, long_name = ff
    row = []
    row.append(name) # CHROM
    if range_start == range_end:
        # This change allow multiple SNV
        POS=range_start
    else:
        POS=np.random.choice(range(range_start, range_end))
    row.append(POS) # POS
    row.append('.') # ID
    REF=sequence[POS-1]
    row.append(REF) # REF # Add strand option from GGF3
    row.append(REF+''.join(np.random.choice(['A', 'T', 'C', 'G'], size,
                            p=[0.25, 0.25, 0.25, 0.25]))) # ALT
    row.append(QUAL) # QUAL
    row.append(FILTER) # FILTER
    row.append(INFO) # INFO
    return row

def random_del(FASTA, range_start=range_start,
            range_end=range_end, QUAL=45, FILTER='PASS', INFO='', size=5):
    """
    Creat a list which contain information on SNP/SNP
    Optimize to parse a GFF3 file and create a VCF file in silico
    Parameters
    ----------
    FASTA : The name of the header and fasta file
      Fasta file, with the path (optionnal)
    range_start : position where the range of insertion start
        int
    range_end : position where the range of insertion end
        int
    size : the size of the insertion
        int
    
    Returns
    -------
    A list with: #CHROM POS ID REF ALT QUAL FILTER
    Examples
    --------
    >>> from insilicogenome import insilicodata
    >>> insilicodata.random_del(/path/genome.fasta, range_start, range_end, size)
    >>> insilicodata.random_del('miniasm.fasta', 3, 36, 5)
    """
    if ((range_start>= 0) and (range_end>0) and (range_end>=range_start)):
        pass
    else:
        raise ValueError(
        "Problems with the range where the SNV occurs. "
        f"The value of range_start is '{range_start}'"
        f"The value of range_end is '{range_end}'")
    # add controle pour s'assurer que les valeur de range sont pas or range de seq
    if not os.path.isfile(FASTA):
        raise FileExistsError
    else:
        pass
    #fasta = insilicogenome.fasta_iter(FASTA)
    #np.random.seed(3)
    fasta = fasta_iter(FASTA)
    for ff in fasta:
        name, sequence, long_name = ff
    row = []
    row.append(name) # CHROM
    if range_start == range_end:
        # This change allow multiple SNV
        POS=range_start
    else:
        POS=np.random.choice(range(range_start, range_end))
    row.append(POS) # POS
    row.append('.') # ID
    REF=sequence[POS-1:POS+size]
    row.append(REF) # REF # Add strand option from GGF3
    row.append(sequence[POS-1]) # ALT
    row.append(QUAL) # QUAL
    row.append(FILTER) # FILTER
    row.append(INFO) # INFO
    return row

def generate_table_small_variation(FASTA, range_start=range_start,
            range_end=range_end, QUAL=45, FILTER='PASS', INFO='', size=5):
    """
    Creat table containing in silico variation close to VCF format
    
    Parameters
    ----------
    FASTA : The name of the header and fasta file
      Fasta file, with the path (optionnal)
    range_start : position where the range of insertion start
        int
    range_end : position where the range of insertion end
        int
    size : the size of the insertion
        int
    
    Returns
    -------
    A CSV file containing all mutated position

    Examples
    --------
    >>> from insilicogenome import insilicodata
    >>> insilicodata.generate_table_small_variation(/path/genome.fasta, range_start, range_end, size)
    >>> insilicodata.generate_table_small_variation('miniasm.fasta', 3, 36, 5)
    """
    if ((range_start>= 0) and (range_end>0) and (range_end>=range_start)):
        pass
    else:
        raise ValueError(
        "Problems with the range where the SNV occurs. "
        f"The value of range_start is '{range_start}'"
        f"The value of range_end is '{range_end}'")
    if not os.path.isfile(FASTA):
        raise FileExistsError
    else:
        pass
    #np.random.seed(4)
    output = os.path.join((os.path.dirname(FASTA)), # head of path
                           os.path.splitext(os.path.basename(FASTA))[0] + 
                           "_" + str(range_start) + "-" + str(range_end) + ".vcf")
    #dtype = {'names':['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'],
    #             'formats':[str] + [int] + 2*[str] + [int] + 2*[str]}
    rows = []
    rows.append(random_snv(FASTA, range_start, range_end, INFO="1SNP"))
    POS=np.random.choice(range(range_start, range_end))
    for x in range(0, size):
        rows.append(random_snv(FASTA, POS+x, POS+x, INFO=str(size)+"MNP")) # pb ? POS+x+1 ->POS+x+1
    rows.append(random_ins(FASTA, range_start, range_end, INFO="small_INS", size=1))
    rows.append(random_del(FASTA, range_start, range_end, INFO="small_DEL", size=1)) ## error in DEL
    np.random.seed(5)
    rows.append(random_ins(FASTA, range_start, range_end, INFO="INS", size=5))
    rows.append(random_del(FASTA, range_start, range_end, INFO="DEL", size=5))
    vcf = np.array(rows, dtype=object)
    np.savetxt(output,vcf[vcf[:, 1].argsort()], delimiter ="\t",  fmt ='% s')
    print(f"INFO:insilicogenome.app: Table {output} completed")
    return vcf

def create_variants(FASTA, vcf, **kwargs):
    """
    Create a fasta file with variation from a reference fasta file
    
    Parameters
    ----------
    FASTA : File containing header and genomic sequence
      Fasta file, with the path (optionnal)
    vcf : Table containing all variation in VCF format
      Numpy array, dtype='<U21'
    size : the size of the insertion
        int
    
    Returns
    -------
    A CSV file containing all mutated position
    Examples
    --------
    >>> from insilicogenome import insilicodata
    >>> vcf = np.genfromtxt(fname=tsv, delimiter="\t", dtype=("|U21"), comments="#") # tansform tsv to np.array
    >>> insilicodata.generate_table_small_variation(/path/genome.fasta, vcf)
    >>> insilicodata.generate_table_small_variation('miniasm.fasta', 3, 36, 5)
    """
    if "range_start" in kwargs:
        range_pos = "_" + str(range_start) + "_" + str(range_end)
    else:
        range_pos=""
    vcf = vcf[vcf[:, 1].argsort()[::-1]]
    vcflist = vcf.tolist()
    fasta = fasta_iter(FASTA)
    for ff in fasta:
        name, sequence, long_name = ff
    i=0
    while i < len(vcflist):
        if vcflist[i][3] == sequence[(int(vcflist[i][1])-1):(int(vcflist[i][1])+len(vcflist[i][3])-1)] :
            # Error SNP it's true
            pass
        else:
            print("Warning: insilicogenome.app: conflict detected\n"
            f"POS: {vcflist[i][1]} -> FASTA file: {sequence[(int(vcflist[i][1])-1):(int(vcflist[i][1])+len(vcflist[i][3])-1)]}\n"
            f"POS: {vcflist[i][1]} -> CSV file: {vcflist[i][3]}")
        sequence =sequence[0:(int(vcflist[i][1])-1)] + vcflist[i][4] + sequence[(int(vcflist[i][1])+len(vcflist[i][3])-1)::]
        i += 1
    output_fasta = name + range_pos + "_variant.fasta"
    insilicogenome.write_fasta_genome(output_fasta, sequence)
    print(f"INFO:insilicogenome.app: variant of {FASTA} have been create"
                f"Filename : {output_fasta}")

def multiple_FASTA_one_variation(FASTA, range_start, range_end, size):
    rows = []
    rows.append(random_snv(FASTA, range_start, range_end, INFO="1SNP")) # altered from orginial
    rows.append(random_mnv(FASTA, range_start, range_end, INFO=str(size)+"MNP", size=2))
    rows.append(random_ins(FASTA, range_start, range_end, INFO="1ins", size=1))
    rows.append(random_ins(FASTA, range_start, range_end, INFO="2ins", size=size))
    rows.append(random_del(FASTA, range_start, range_end, INFO="1del", size=1))
    rows.append(random_del(FASTA, range_start, range_end, INFO="2del", size=size))
    vcf = np.array(rows, dtype=object)
    output = os.path.join((os.path.dirname(FASTA)), # head of path
                               os.path.splitext(os.path.basename(FASTA))[0] + 
                               "_" + str(range_start) + "-" + str(range_end) + ".csv")
    np.savetxt(output,vcf[vcf[:, 1].argsort()], delimiter ="\t",  fmt ='% s')
    fasta = fasta_iter(FASTA)
    for ff in fasta:
        name, sequence, long_name = ff
    vcflist = vcf.tolist()
    for i in vcflist:
        sequence_tmp = sequence[0:(int(i[1])-1)] + i[4] + sequence[(int(i[1])+len(i[3])-1)::]
        output_fasta = name + "_" + str(i[7]) + "_POS"+ str(i[1]) + "_var" + ".fasta"
        insilicogenome.write_fasta_genome(output_fasta, sequence_tmp, description='variant')

def tsv_to_structured_array(tsv_file):
    """
    Take TSV file into numpy structured array

    Parameters
    ----------
    TSV : Tabular file corresponding to a VCF-like foramt
      TSV file, with the path (optionnal)
    
    Returns
    -------
    A numpy array

    Examples
    --------
    >>> from insilicogenome import insilicodata
    >>> vcf = insilicodata.tsv_to_structured_array(/path/file.tsv)
    """
    vcf=[]
    with open(tsv_file) as tsvfile:
        inputItems = list(csv.reader(tsvfile, delimiter="\t"))
        for rows in inputItems:
            vcf.append(rows)
    vcf = np.array(vcf, dtype=object)
    vcf[:,1]=vcf[:,1].astype('int')
    return vcf

def main():
    if opt not in ("-c", "--csv"):
        vcf = generate_table_small_variation(input_fasta, range_start=range_start, range_end=range_end)
        create_variants(input_fasta, vcf, range_start=range_start, range_end=range_end)
    else:
        print(f"INFO:insilicogenome.app: Table {csv_file} provided;"
              "Skipping generating variation file\n"
              f"Start transform {os.path.basename(input_fasta)}")
        vcf = tsv_to_structured_array(csv_file)
        create_variants(input_fasta, vcf)
    

#dtype = {'names':['col%i'%i for i in range(len(rows))],
#                 'formats':2*[str] + 2*[int] + 2*[np.int] + 2*[np.bool] + 3*[np.int]}

if __name__ == "__main__":
    main()
#vcf = vcf[vcf[:, 1].argsort()[::-1]]
#vcflist = vcf.tolist()
#fasta = fasta_iter(FASTA)
#for ff in fasta:
#    name, sequence, long_name = ff
#sequence[:start] + ''.join(codon_start) + sequence[(start+3)]
#if 
###Faut prendre en compte si c'est une deletion où une insertion peut être un systeme de taille genre si alt > ref etc..
## Check si la ref est bien la ref en amon et levé un warning c'est pas le cas
#
#
#
#int(vcflist[i][3])]
## position vcflist[i][1]
#
## on est en parsing antisens
#arr[..., 0]
#pos
#sequence[:int(vcflist[i][1])] + ''.join(codon_start) +
