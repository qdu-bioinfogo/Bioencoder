# Bioencoder

An amino acid sequence encoding toolbox for machine learning.

## Introduction

main features:

- Machine learning oriented
- Rich encoding varieties (NUM, BE, EAAC, AAINDEX, GACC, CKSAAP)
- Native support for big size `fasta` format
- Out-of-the-box

## Installation Tutorial

### Via PIP

```bash
python setup.py sdist
pip install ./dist/bioencoder-beta.tar.gz
```

## Usage

### Reading from a Fasta File

A standard `fasta` file like:

```fasta
>1|1
DGMRITLRDGCIVHLRASGNAPELRCYAEANLLNRAQDLVNTTLANIKKRC
>2|1
EGKLSMLQNTIKRLASLSTEEPVVICNDRHRFLVAEQLREIDKLANNIILE
```

To read and process the sequences to EEAC Embedding:

```python
from bioencoder import *
pos_data = "pos.fasta"
window_size = 12
pos_seqList,pos_labellist,pos_seqNamelist=get_data(pos_data,1,method="GAAC",window_size=window_size)
```

### Reading from a raw sequence

For example, A Str like`str='DGMRITLRDGCIVHLRASGNAPELRCYAEANLLNRAQDLVNTTLANIKKRC'`, Using bellow code to encoding the sequence to `EAAC` Embeding:

```python
from bioencoder.encoder import *
print(EAAC(seq,window=5))
```