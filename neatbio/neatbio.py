#!/usr/bin/env python
import random
import sys


# Translate
CodonTable = {
    # 'M' - START, '*' - STOP
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TGT": "C",
    "TGC": "C",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TTT": "F",
    "TTC": "F",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
    "CAT": "H",
    "CAC": "H",
    "ATA": "I",
    "ATT": "I",
    "ATC": "I",
    "AAA": "K",
    "AAG": "K",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATG": "M",
    "AAT": "N",
    "AAC": "N",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGA": "R",
    "AGG": "R",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "AGT": "S",
    "AGC": "S",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TGG": "W",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGA": "*",
}


class Sequence(object):
    """Create a Valid Sequence for DNA,RNA

    example: seq1 = Sequence('ATGC')

    """

    def __init__(self, seq=None,seqtype="DNA"):
        """Creat A Sequence Object
        
        Arguments:
        - seq - Sequence required string

        >>>import neatbio as nt
        >>>seq1 = nt.Sequence("ATGC")
        >>>seq1
        Sequence(seq="ATGC")

        """
        super(Sequence, self).__init__()
        self.seq = seq
        self.seqtype = seqtype
        self.is_dna = self.__check_sequence_type()
        assert self.is_dna, f'Sequence Type Not A {self.seqtype}, Please set seqtype as "RNA" '

        # To enforce a string storage
        if not isinstance(self.__validate_seq(seq), str):
            raise TypeError("The sequence data given to a Sequence object should "
                "be a string (not another Sequence object etc)","nor a non Nucleotide [A,C,T,G,U]")
        

    def __repr__(self):
        return "Sequence(seq='{}',seqtype='{}')".format(self.seq,self.seqtype)

    def __str__(self):
        return self.seq

    def __validate_seq(self, seq):
        base_nucleotide = ["A", "T", "G", "C", "U"]
        real_seq = seq.upper()
        for base in real_seq:
            if base not in base_nucleotide:
                return False
        return real_seq

    def __check_sequence_type(self):
        if "U" not in self.seq:
            return self.seqtype == "DNA"
        else:
            return self.seqtype == "RNA"
    

    def __len__(self):
        return len(self.seq)

    def __contains__(self, sub_char):
        return sub_char in str(self)

    def __getitem__(self, index):
        if isinstance(index, int):
            return self.seq[index]

    def __add__(self,other):
        if isinstance(other,Sequence):# works instead of a str
            return self.__class__(str(self.seq) + str(other.seq))
        else:
            # return str(self.seq) + str(other.seq)
            return 'Error: NotImplemented Error, Add a Valid Sequence'

    def __hash__(self):
        """Return a Hash of the sequences for string comparison"""
        return hash(str(self))



    ### Basic Fxn
    # This include basic string count,find,index

    def count(self, subseq, start=0, end=sys.maxsize):
        """Return the Count of the Number of Nucleotide in A Sequence"""
        return str(self).count(subseq, start, end)

    def find(self, subseq, start=0, end=sys.maxsize):
        """Find the Position of A  Nucleotide in A Sequence"""
        return str(self).find(subseq, start, end)

    def rfind(self, subseq, start=0, end=sys.maxsize):
        """Find the Position of A  Nucleotide in A Sequence From the Right"""
        return str(self).rfind(subseq, start, end)

    def index(self, subseq, start=0, end=sys.maxsize):
        """Find the Index/Position of A  Nucleotide in A Sequence"""
        return str(self).index(subseq, start, end)

    def rindex(self, subseq, start=0, end=sys.maxsize):
        """Find the Index/Position of A  Nucleotide in A Sequence From the Right"""
        return str(self).rindex(subseq, start, end)

    def encode(self, encoding="utf-8",errors="strict"):
        """Returns an encoded version of the sequence as a byte object
        >>>import neatbio as nt
        >>>nt.Sequence("ATGC").encode("ascii")
        b'ATGC'

        """
        return str(self).encode(encoding,errors)


    
    #### Main Functions For DNA/RNA Sequencing
    def get_symbol_frequency(self):
        """Get the Frequency of A Nucleotide in a Sequence """
        base_dict = {"A": 0, "T": 0, "G": 0, "C": 0}  # initial score
        if self.__validate_seq(str(self.seq)) != False:
            for base in self.seq:
                base_dict[base] += 1
        else:
            return "NucleotideError: {} not a nucleotide ['A,T,G,C']".format(base)
        return base_dict

    @property
    def gc(self):
        """GC Content of Sequence """
        result = (
            float(str(self.seq).count("G") + str(self.seq).count("C"))
            / len(self.seq)
            * 100
        )
        return result

    @property
    def at(self):
        """AT Content of Sequence """
        result = (
            float(str(self.seq).count("A") + str(self.seq).count("T"))
            / len(self.seq)
            * 100
        )
        return result

   
    def at_frequency(self):
        """AT Count/Frequence of Sequence"""
        result = str(self.seq).count("A") + str(self.seq).count("T")
        return result

 
    def gc_frequency(self):
        """GC Count/Frequence of Sequence """
        result = str(self.seq).count("G") + str(self.seq).count("C")
        return result
 


    def complement(self):
        """Return Complement of Sequence """
        base_pairs = {"A": "T", "T": "A", "G": "C", "C": "G"}
        comp_pairs = [base_pairs[a] for a in self.seq if a in base_pairs.keys()]
        result_complement = "".join(comp_pairs)
        return Sequence(result_complement)

    def reverse_complement(self):
        """Return A Reverse Complement of Sequence """
        base_pairs = {"A": "T", "T": "A", "G": "C", "C": "G"}
        comp_pairs = [base_pairs[a] for a in self.seq if a in base_pairs.keys()]
        reverse_pairs = "".join(comp_pairs)[::-1]
        return Sequence(reverse_pairs)

    def transcribe(self):
        """Transcribe Sequence into mRNA """
        return Sequence(self.seq.replace("T","U"),"RNA")

    def back_transcribe(self):
        """Transcribe mRNA Sequence back into DNA """
        dna_result = self.seq.replace("U", "T")
        return Sequence(dna_result)

    def translate(self, start_pos=0):
        """Translate Sequence into Protein/Amino Acids"""
        # Convert mrna to dna for codontable
        if 'U'in  self.seq:
            self.seq = self.seq.replace("U","T")
            amino_acids_list = [
            CodonTable[self.seq[pos : pos + 3]]
            for pos in range(start_pos, len(self.seq) - 2, 3)]
        else:
            amino_acids_list = [
            CodonTable[self.seq[pos : pos + 3]]
            for pos in range(start_pos, len(self.seq) - 2, 3)
            ]

        return "".join(amino_acids_list)

    def hamming_distance(self,other):
        """Returns the edit distance for similarity using Hamming Distance of 2 Equal Length Sequence"""
        if isinstance(other,Sequence):# works instead of a str
            return len([(x,y) for x,y in zip(self.seq,other) if x !=y])
        else:
            return f"Error: Not Implemented '{other}' Must be of a Sequence Object"

    def randomize(self):
        """Randomly Change the Location of Each Nucletide

        >>> import neatbio as nt
        >>> seq1 = nt.Sequence("ATGC")
        >>> seq1.randomize()
        >>> seq1
        Sequence(seq='CACT')

        """
        result = ''.join([random.choice(list(self.seq)) for x in range(len(self.seq)) ])
        return self.__init__(result) # Reinit as a valid seq

    def random_gen(self,length=9):
        """Randomly Generate A Sequence of N- Length"""
        result = ''.join([random.choice(list(self.seq)) for x in range(length) ])
        return result 


class ProteinSeq(object):
    def __init__(self, seq=None):
        """Creat A Protein Sequence Object
        
        Arguments:
        - seq - Sequence required string

        >>>import neatbio as nt
        >>>seq1 = nt.ProteinSeq("MIT")
        >>>seq1
        ProteinSeq(seq="MIT")

        """
        super(ProteinSeq, self).__init__()
        self.seq = seq

        # To enforce a string storage
        if not isinstance(self.__validate_protein_seq(seq), str):
            raise TypeError("The sequence data given to a Protein Sequence object should "
                "be a string (not another Protein Sequence object etc)","nor a non Amino Acid Code")
    

    def __repr__(self):
        return "ProteinSeq(seq='{}')".format(self.seq)

    def __str__(self):
        return self.seq

    def __validate_protein_seq(self, seq):
        base_amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y','*','-']
        real_seq = seq.upper()
        for base in real_seq:
            if base not in base_amino_acids:
                return False
        return real_seq

    def __len__(self):
        return len(self.seq)

    def __contains__(self, sub_char):
        return sub_char in str(self)

    def __getitem__(self, index):
        if isinstance(index, int):
            return self.seq[index]

    def __add__(self,other):
        if isinstance(other,ProteinSeq):# works instead of a str
            return self.__class__(str(self.seq) + str(other.seq))
        else:
            # return str(self.seq) + str(other.seq)
            return 'Error: NotImplemented Error, Add a Valid Protein Sequence'

    def __hash__(self):
        """Return a Hash of the sequences for string comparison"""
        return hash(str(self))

    # Get Keys
    def __get_codon_key(self,val):
        my_dict = CodonTable
        for key,value in my_dict.items():
            if val == value:
                return key


    def get_symbol_frequency(self):
        """Get the Frequency or Counts of Amino Acids"""
        codon_base_dict = {'A': 0,'C': 0,'D': 0,'E': 0,'F': 0,'G': 0,'H': 0,'I': 0,'K': 0,'L': 0,'M': 0,'N': 0,'P': 0,'Q': 0,'R': 0,'S': 0,'T': 0,'V': 0,'W': 0,'Y': 0,'*':0,'-':0} # initial score
        if self.__validate_protein_seq(base) != False:
            for base in self.seq:
                codon_base_dict[base] +=1
        else:
            return "AminoAcidError: {} not a valid Amino Acid Sequence".format(base)
        return codon_base_dict

    def hamming_distance(self,other):
        """Returns the edit distance for similarity using Hamming Distance of 2 Equal Length Sequence"""
        if isinstance(other,ProteinSeq):# works instead of a str
            return len([(x,y) for x,y in zip(self.seq,other) if x !=y])
        else:
            return f"Error: Not Implemented '{other}' Must be of a Protein Sequence Object"

    def back_translate(self):
        """Returns a Probable Nucleotide Sequence of A Protein

        example:
        >>> p2 = ProteinSeq('IKGLYPR')
        >>> p2.back_translate()
        'ATAAAAGGTTTATATCCTCGT'

        """

        base_nucleotide_list = []
        for i in self.seq:
            res = self.__get_codon_key(i)
            base_nucleotide_list.append(res)
        return Sequence(''.join(base_nucleotide_list))


    def get_amino_acid_percentage(self):
        """Returns a Dictionary of Percentage of Amino Acid for the Sequence"""
        percentages = {}

        aminoacid_freq = self.get_symbol_frequency()

        for acid in aminoacid_freq:
            percentages[acid] = aminoacid_freq[acid]/float(len(self.seq))

        return percentages

    def aromaticity(self):
        """Returns the Lobry Aromaticity using the relative frequency of 
        sum of Phe,Trp,Tyr
        """
        acid_percentage = self.get_amino_acid_percentage()
        aromaticity_score = sum(acid_percentage[acid] for acid in "YWF")
        return aromaticity_score

    def molar_coeff(self):
        """Returns the Molar Extinction Coefficient
        using Cys-Cys bond
        
        """
        aminoacid_freq = self.get_symbol_frequency()
        m_coef_reduced  = aminoacid_freq["W"] * 5500 + aminoacid_freq["Y"] * 1490
        m_coef_cystine_bond = m_coef_reduced + (aminoacid_freq["C"]// 2) * 125
        return {'Reduced':m_coef_reduced,"Cys-Cys Residue":m_coef_cystine_bond}


