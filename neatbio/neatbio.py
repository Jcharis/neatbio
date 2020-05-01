#!/usr/bin/env python
import sys
import random

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

class Sequence(object):
	"""docstring for Sequence"""


	def __init__(self, seq=None):
		super(Sequence, self).__init__()
		self.seq = seq
		# Enforce string storage
		if not isinstance(self._validate_seq(seq), str):
			raise TypeError("The sequence data given to a Sequence object should "
				"be a string (not another Sequence object etc)","nor a non Nucleotide [A,C,T,G,U]")
	


	def __repr__(self):
		return 'Sequence(seq="{}")'.format(self.seq)

	def __str__(self):
		return self.seq 

	def __len__(self):
		return len(self.seq)

	def __contains__(self,sub_char):
		return sub_char in str(self)

	def __getitem__(self,index):
		if isinstance(index,int):
			return self.seq[index]
		else:
			return Sequence(self.seq[index])

	def __add__(self,other):
		if isinstance(other,Sequence):# works instead of a str
			return self.__class__(str(self.seq) + str(other.seq))
		else:
			# return str(self.seq) + str(other.seq)
			return 'Error: NotImplemented Error, Add a Valid Sequence'

	# Basic Fxn
	def count(self,subseq,start=0,end=sys.maxsize):
		return str(self).count(subseq,start,end)

	def find(self,subseq,start=0,end=sys.maxsize):
		return str(self).find(subseq,start,end)

	def rfind(self,subseq,start=0,end=sys.maxsize):
		return str(self).rfind(subseq,start,end)

	def index(self,subseq,start=0,end=sys.maxsize):
		return str(self).index(subseq,start,end)

	def rindex(self,subseq,start=0,end=sys.maxsize):
		return str(self).rindex(subseq,start,end)


	def _validate_seq(self,seq):
	    base_nucleotide = ["A","T","G","C","U"]
	    real_seq = seq.upper()
	    for base in real_seq:
	        if base not in base_nucleotide:
	            return False
	    return real_seq

	def _validate(self):
		base_nucleotide = ["A","T","G","C","U"]
		return set(base_nucleotide).issuperset(self.seq)
		

	def get_symbol_frequency(self):
	    base_dict = {"A":0,"T":0,"G":0,"C":0} # initial score
	    for base in self.seq:
	        if self._validate_seq(base) != False:
	            base_dict[base] +=1
	        else:
	            return "NucleotideError: {} not a nucleotide ['A,T,G,C']".format(base)
	    return base_dict

	def hamming_distance(self,other):
		if isinstance(other,Sequence):# works instead of a str
			return len([(x,y) for x,y in zip(self.seq,other) if x !=y])
		else:
			return f"Error: Not Implemented '{other}' Must be of a Sequence Object"

	def randomize(self):
		result = ''.join([random.choice(list(self.seq)) for x in range(len(self.seq)) ])
		return self.__init__(result) # Reinit as a valid seq

	def random_gen(self,length=9):
		result = ''.join([random.choice(list(self.seq)) for x in range(length) ])
		return result # Reinit as a valid seq
	
	# GC Content
	def gc(self):
	    result = float(str(self.seq).count('G') + str(self.seq).count('C'))/len(self.seq) *100
	    return result


	# GC Content
	def at(self):
	    result = float(str(self.seq).count('A') + str(self.seq).count('T'))/len(self.seq) *100
	    return result


	def complement(self):
	    base_pairs = {"A":"T","T":"A","G":"C","C":"G"}
	    comp_pairs = [base_pairs[a] for a in self.seq if a in base_pairs.keys()]
	    return "".join(comp_pairs)

	def reverse_complement(self):
	    base_pairs = {"A":"T","T":"A","G":"C","C":"G"}
	    comp_pairs = [base_pairs[a] for a in self.seq if a in base_pairs.keys()]
	    reverse_pairs = "".join(comp_pairs)[::-1]
	    return reverse_pairs

	def transcribe(self):
	    mrna_result = self.seq.replace("T","U")
	    return mrna_result

	def translate(self,start_pos=0):
	    amino_acids_list =[CodonTable[self.seq[pos:pos + 3]] for pos in range(start_pos,len(self.seq)-2,3)]
	    return "".join(amino_acids_list)



class ProteinSeq(object):
	"""docstring for ProteinSeq"""
	def __init__(self, seq):
		super(ProteinSeq, self).__init__()
		self.seq = seq
		
		# Enforce string storage
		if not isinstance(self._validate_protein_seq(seq), str):
			raise TypeError("The sequence data given to a Sequence object should "
				"be a string (not another Sequence object etc)","nor a non Nucleotide [A,C,T,G,U]")
	


	def __repr__(self):
		return 'ProteinSeq(seq="{}")'.format(self.seq)

	def __str__(self):
		return self.seq 

	def __len__(self):
		return len(self.seq)

	def __contains__(self,sub_char):
		return sub_char in str(self)

	def __getitem__(self,index):
		if isinstance(index,int):
			return self.seq[index]
		else:
			return Sequence(self.seq[index])

	def __add__(self,other):
		if isinstance(other,ProteinSeq):
			return self.__class__(str(self) + str(other))
		else:
			return 'Error: NotImplemented ,Use A Valid Protein Sequence'

	# Basic Fxn
	def count(self,subseq,start=0,end=sys.maxsize):
		return str(self).count(subseq,start,end)

	def find(self,subseq,start=0,end=sys.maxsize):
		return str(self).find(subseq,start,end)

	def rfind(self,subseq,start=0,end=sys.maxsize):
		return str(self).rfind(subseq,start,end)

	def index(self,subseq,start=0,end=sys.maxsize):
		return str(self).index(subseq,start,end)

	def rindex(self,subseq,start=0,end=sys.maxsize):
		return str(self).rindex(subseq,start,end)


	def _validate_protein_seq(self,seq):
	    codon_amino_acid = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
	    real_seq = seq.upper()
	    for base in real_seq:
	        if base not in codon_amino_acid:
	            return False
	    return real_seq
		

	def get_symbol_frequency(self):
	    codon_base_dict = {'A': 0,'C': 0,'D': 0,'E': 0,'F': 0,'G': 0,'H': 0,'I': 0,'K': 0,'L': 0,'M': 0,'N': 0,'P': 0,'Q': 0,'R': 0,'S': 0,'T': 0,'V': 0,'W': 0,'Y': 0} # initial score
	    for base in self.seq:
	        if self._validate_protein_seq(base) != False:
	            codon_base_dict[base] +=1
	        else:
	            return "AminoAcidError: {} not a valid Amino Acid Sequence".format(base)
	    return codon_base_dict


# # GC Content
# def gc_content(seq):
# 	result = float(str(seq).count('G') + str(seq).count('C'))/len(seq) *100
# 	return result

# # AT Content
# def at_content(seq):
# 	result = float(str(seq).count('A') + str(seq).count('T'))/len(seq) *100
# 	return result


# 	# Get Keys
# def get_key(val,my_dict):
# 	for key,value in my_dict.items():
# 		if val == value:
# 			return key

# def get_value(val,my_dict):
# 	for key,value in my_dict.items():
# 		if val == key:
# 			return value


# # Method 1: Converting 1 to 3 letter
# def convert_1to3(seq):
#     term_list = []
#     for i in seq:
#         res = get_key(i,aa3_to1_dict)
#         term_list.append(res)
#     return "".join(term_list)

# def _kmers(seq,k=2):
#     pair_list = []
#     for i in range(0,len(seq),k):
#         pair_list.append(seq[i:i+k])
#     return pair_list

# def convert_3to1(seq):
#     term_list = []
#     for i in _kmers(seq,k=3):
#         res = get_value(i,aa3_to1_dict)
#         term_list.append(res)
#     return ''.join(term_list)

# def get_kmers(seq,k=2):
#     pair_list = []
#     for i in range(0,len(seq),k):
#         pair_list.append(str(seq)[i:i+k])
#     return pair_list



# def _delta(x,y):
#     return 0 if x == y else 1


# def _M(seq1,seq2,i,j,k):
# 	return sum(_delta(x,y) for x,y in zip(seq1[i:i+k],seq2[j:j+k]))


# def _makeMatrix(seq1,seq2,k):
#     n = len(seq1)
#     m = len(seq2)
#     return [[_M(seq1,seq2,i,j,k) for j in range(m-k+1)] for i in range(n-k+1)]


# def _plotMatrix(M,t, seq1, seq2, nonblank = chr(0x25A0), blank = ' '):
#     print(' |' + seq2)
#     print('-'*(2 + len(seq2)))
#     for label,row in zip(seq1,M):
#         line = ''.join(nonblank if s < t else blank for s in row)
#         print(label + '|' + line)

# def dotplot(seq1,seq2,k = 1,t = 1):
# 	"""Create a Simple Dotplot(Black and white) of 2 sequences """
# 	M = _makeMatrix(str(seq1),str(seq2),k)
# 	_plotMatrix(M, t, str(seq1),str(seq2)) #experiment with character choice