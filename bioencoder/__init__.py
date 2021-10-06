from encoder import *
from utils import *
import numpy as np

def get_data(fname,label,method="None",window_size=15,encode_window_size=5):
    seqlist=[]
    labellist=[]
    seqNamelist=[]
    fa = FastaReader(fname)
    for seqName,seq_record in enumerate(fa.readAll()):
        seqNamelist.append(seqName)
        raw_seq = str(seq_record).replace("*","-")
        seq =  split_seq(raw_seq,window_size)
        result = seq
        if method == "BE":
            result = Binary(seq)
        elif method == "EAAC":
            result = EAAC(seq,encode_window_size)
        elif method == "GAAC":
            result = GAAC(seq)
        elif method == "AAINDEX":
            result = AAINDEX(seq)
        elif method == "NUM":
            result = NUM(seq)
        elif method == "CKSAAP":
            result = CKSAAP(seq,encode_window_size)
        seqlist.append(result)
        labellist.append(label)
    seqlist = np.array(seqlist)
    labellist = np.array(labellist)
    return seqlist,labellist,seqNamelist