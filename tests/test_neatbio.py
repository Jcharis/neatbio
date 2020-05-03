from neatbio import __version__
import neatbio as nt
from neatbio.sequtils import gc_content,at_content,get_acid_name 


def test_version():
    assert __version__ == '0.0.1'

def test_issequence():
	seq1 = nt.Sequence('ATGCTATGCTT')
	assert isinstance(seq1,type(seq1))

def test_complement():
	seq1 = nt.Sequence('ATGCTATGCTT')
	res_complement = seq1.complement()
	assert res_complement == 'TACGATACGAA'
	

def test_reverse_complement():
	seq1 = nt.Sequence('ATGCTATGCTT')
	res_complement = seq1.reverse_complement()
	assert res_complement == 'AAGCATAGCAT'

def test_translate():
	seq1 = nt.Sequence('ATGCTATGCTT')
	result = seq1.translate()
	assert result == 'MLC'
	
def test_transcription():
	seq1 = nt.Sequence('ATGCTATGCTT')
	result = seq1.transcribe()
	assert result == 'AUGCUAUGCUU'

def test_nucleotide_freq():
	seq1 = nt.Sequence('ATGCTATGCTT')
	result = seq1.get_symbol_frequency()
	assert result == {'A': 1, 'T': 0, 'G': 0, 'C': 0}

	

def test_gc_content():
	seq1 = nt.Sequence('ATGCTATGCTT')
	result = gc_content(seq1)
	assert result == 36.36363636363637

def test_isProteinsequence():
	seq1 = nt.ProteinSeq('MIT')
	assert isinstance(seq1,type(seq1))

def test_back_translate():
	protein1 = nt.ProteinSeq('IKGLYLPR')
	nucl_seq = protein1.back_translate()
	assert nucl_seq == 'ATAAAAGGTTTATATTTACCTCGT'

def test_get_acid_name():
	assert get_acid_name('Ala') == 'Alanine'


