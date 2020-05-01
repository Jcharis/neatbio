#!/usr/bin/env python
from collections import Counter

# Translate
CodonTable = {
    # 'M' - START, '*' - STOP
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TGT": "C", "TGC": "C",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TTT": "F", "TTC": "F",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATA": "I", "ATT": "I", "ATC": "I",
    "AAA": "K", "AAG": "K",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATG": "M",
    "AAT": "N", "AAC": "N",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "TAA": "*", "TAG": "*", "TGA": "*"
}


full_amino_acid_name ={'Alanine': 'Ala',
 'Cysteine': 'Cys',
 'Aspartic acid': 'Asp',
 'Glutamic acid': 'Glu',
 'Phenylalanine': 'Phe',
 'Glycine': 'Gly',
 'Histidine': 'His',
 'Isoleucine': 'Ile',
 'Lysine': 'Lys',
 'Leucine': 'Leu',
 'Methionine': 'Met',
 'Asparagine': 'Asn',
 'Proline': 'Pro',
 'Glutamine': 'Gln',
 'Arginine': 'Arg',
 'Serine': 'Ser',
 'Threonine': 'Thr',
 'Valine': 'Val',
 'Tryptophan': 'Trp',
 'Tyrosine': 'Tyr'}

aa3_to1_dict  = {'Ala': 'A',
 'Cys': 'C',
 'Asp': 'D',
 'Glu': 'E',
 'Phe': 'F',
 'Gly': 'G',
 'His': 'H',
 'Ile': 'I',
 'Lys': 'K',
 'Leu': 'L',
 'Met': 'M',
 'Asn': 'N',
 'Pro': 'P',
 'Gln': 'Q',
 'Arg': 'R',
 'Ser': 'S',
 'Thr': 'T',
 'Val': 'V',
 'Trp': 'W',
 'Tyr': 'Y'}


 
# GC Content
def gc_content(seq):
	result = float(str(seq).count('G') + str(seq).count('C'))/len(seq) *100
	return result

# AT Content
def at_content(seq):
	result = float(str(seq).count('A') + str(seq).count('T'))/len(seq) *100
	return result


	# Get Keys
def _get_key(val,my_dict):
	for key,value in my_dict.items():
		if val == value:
			return key

def _get_value(val,my_dict):
	for key,value in my_dict.items():
		if val == key:
			return value


# Method 1: Converting 1 to 3 letter
def convert_1to3(seq):
    term_list = []
    for i in seq:
        res = _get_key(i,aa3_to1_dict)
        term_list.append(res)
    return "".join(term_list)

def _kmers(seq,k=2):
    pair_list = []
    for i in range(0,len(seq),k):
        pair_list.append(seq[i:i+k])
    return pair_list

def convert_3to1(seq):
    term_list = []
    for i in _kmers(seq,k=3):
        res = _get_value(i,aa3_to1_dict)
        term_list.append(res)
    return ''.join(term_list)

def get_kmers(seq,k=2):
    pair_list = []
    for i in range(0,len(seq),k):
        pair_list.append(str(seq)[i:i+k])
    return pair_list


def hamming_distance(lhs,rhs):
    return len([(x,y) for x,y in zip(lhs,rhs) if x !=y])

def occurence(main_seq,sub_seq):
    start= 0
    indices =[]
    while True:
        start = main_seq.find(sub_seq,start)
        if start > 0:
            indices.append(start)
        else:
            break
        start +=1
    return indices


def get_acid_name(seq):
    term_list = []
    for i in _kmers(seq,k=3):
        res = _get_key(i,full_amino_acid_name)
        term_list.append(res)
    return ''.join(term_list)


def codon_frequency(seq, aminoacid):
    """Provides the frequency of each codon encoding a given aminoacid in a DNA sequence"""
    tmpList = []
    for i in range(0, len(seq) - 2, 3):
        if CodonTable[seq[i:i + 3]] == aminoacid:
            tmpList.append(seq[i:i + 3])

    freqDict = dict(Counter(tmpList))
    totalScore = sum(freqDict.values())
    for seq in freqDict:
        freqDict[seq] = round(freqDict[seq] / totalScore, 2)
    return freqDict



# For Dot plot
def _delta(x,y):
    return 0 if x == y else 1


def _M(seq1,seq2,i,j,k):
	return sum(_delta(x,y) for x,y in zip(seq1[i:i+k],seq2[j:j+k]))


def _makeMatrix(seq1,seq2,k):
    n = len(seq1)
    m = len(seq2)
    return [[_M(seq1,seq2,i,j,k) for j in range(m-k+1)] for i in range(n-k+1)]


def _plotMatrix(M,t, seq1, seq2, nonblank = chr(0x25A0), blank = ' '):
    print(' |' + seq2)
    print('-'*(2 + len(seq2)))
    for label,row in zip(seq1,M):
        line = ''.join(nonblank if s < t else blank for s in row)
        print(label + '|' + line)

def dotplot(seq1,seq2,k = 1,t = 1):
	"""Create a Simple Dotplot(Black and white) of 2 sequences """
	M = _makeMatrix(str(seq1),str(seq2),k)
	_plotMatrix(M, t, str(seq1),str(seq2)) #experiment with character choice