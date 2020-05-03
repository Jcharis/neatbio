### neatbio
A Simple Yet Another Bioinformatics Library for DNA,RNA and Protein Sequencing

### Installation
```bash
pip install neatbio
```

### Benefits
+ Handling Sequences(DNA,RNA,Protein)
+ Protein Synthesis
+ Sequence Similarity
+ Kmers Generation and Kmer Distance
+ Probable Back Translation of Amino Acids
+ Reading FASTA files

### Why NeatBio?
+ NeatBio is yet another bioinformatics library along side powerful and popular bioinformatics libraries such as biopython,scikit-bio,biotite.
+ It is meant to complement these powerful library in a simple way.


### Usage
#### Handling Sequences
+ Neatbio offers the ability to analyse sequences for more insight

```python
>>>import neatbio as nt
>>> seq1 = nt.Sequence('ATGCATTGA')
>>> seq1.gc
33.33333333333333
>>> seq1.gc_frequency()
3
>>> seq1.at
66.66666666666666
>>> seq1.at_frequency()
6
>>> seq1.transcribe()
'AUGCAUUGA'
>>> mrna = seq1.transcribe()
>>> nt.Sequence(mrna).back_transcribe()
'ATGCATTGA'
>>> 

>>> seq1.translate()
'MH*'
>>> seq1.translate
```

#### Working with Proteins
```python
>>> protein1 = nt.ProteinSeq('MIT')
>>> protein1
ProteinSeq(seq='MIT')
>>> protein1.back_translate()
'ATGATAACT'

```
+ Note that the back_translate() function offers a probable sequence and not the exact
back-translation as multiple codons can represent the same amino acids.

#### Convert 3 Letter Amino Acid to 1 and vice versa
```python
>>> from neatbio.sequtils import convert_3to1,convert_1to3,get_acid_name
>>> convert_3to1('Ala')
'A'
>>> convert_1to3('L')
'Leu'

>>> get_acid_name('Ala')
'Alanine'
```

#### Generate DotPlot
```python
>>> import neatbio as nt 
>>> import neatbio.sequtils as utils
>>> seq1 = nt.Sequence('AGTCGTACT')
>>> seq2 = nt.Sequence('AGGCGCACT')
>>> 
>>> utils.dotplot(seq1,seq2)
 |AGGCGCACT
-----------
A|■     ■  
G| ■■ ■    
T|        ■
C|   ■ ■ ■ 
G| ■■ ■    
T|        ■
A|■     ■  
C|   ■ ■ ■ 
T|        ■
>>> 

```

#### Reading FASTA Files
```python
>>> import neatbio as nt 
>>> file1 = nt.read_fasta('sequence.fasta')
>>> file1['seqRecord']


>>> seq1 = nt.Sequence(file1['seqRecord'])
```


#### Documentation
+ Please read the [documentation](https://github.com/Jcharis/neatbio/wiki) for more information on what neatbio does and how to use is for your needs.

#### More Features To Add
+ sequence alignment
+ writing FASTA files
+ support for more file formats



#### Acknowledgements
   + Inspired by packages like BioPython,Scikit-Bio and Biotite

### NB
+ Contributions Are Welcomed
+ Notice a bug, please let us know.
+ Thanks A lot

### By
+ Jesse E.Agbe(JCharis)
+ Jesus Saves @JCharisTech
