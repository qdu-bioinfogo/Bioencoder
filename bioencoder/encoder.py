import re, sys, os, platform,warnings
import math
import numpy
from collections import Counter
from .utils import *
def NUM(sequence):
    AA = 'ARNDCQEGHILKMFPSTWYV-'
    code = []
    for aa in sequence.upper() :
        code.append(AA.find(aa))
    code = numpy.array(code).reshape(-1,1)
    return code
def Binary(sequence):
    #FINISHED
    AA = 'ARNDCQEGHILKMFPSTWYV'
    length = len(sequence)
    code=list()
    for aa in sequence:
        if aa == '-':
            code = code + [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            continue
        for aa1 in AA:
            tag = 1 if aa == aa1 else 0
            code.append(tag)
    code = numpy.array(code).reshape(length,20)
    return code
def EAAC(sequence,window=5):
    #FINISHED
    length = (len(sequence)-window+1)
    AA = 'ARNDCQEGHILKMFPSTWYV'
    code=list()
    for j in range(len(sequence)):
        if j < len(sequence) and j + window <= len(sequence):
            count = Counter(re.sub('-', '', sequence[j:j+window]))
            for key in count:
                count[key] = count[key] / len(re.sub('-', '', sequence[j:j+window]))
            for aa in AA:
                code.append(count[aa])
    code = numpy.array(code).reshape(length,20)
    return code
def GAAC(sequence):
    group = {
        'alphatic': 'GAVLMI',
        'aromatic': 'FYW',
        'postivecharge': 'KRH',
        'negativecharge': 'DE',
        'uncharge': 'STCPNQ'
    }
    code = []
    groupKey = group.keys()
    for i in range(len(sequence)):
        
        count = Counter(sequence)
        myDict = {}
        for key in groupKey:
            for aa in group[key]:
                myDict[key] = myDict.get(key, 0) + count[aa]
        for key in groupKey:
            code.append(myDict[key]/len(sequence))
    code = numpy.array(code).reshape(-1,len(groupKey))
    return code
def CKSAAP(sequence, gap=5):
    # ref  http://www.nohup.cc/article/110/
    AA = 'ARNDCQEGHILKMFPSTWYV'
    code=[]
    aaPairs = []
    for aa1 in AA:
        for aa2 in AA:
            aaPairs.append(aa1 + aa2)
    for g in range(gap+1):
            myDict = {}
            for pair in aaPairs:
                myDict[pair] = 0
            sum = 0
            for index1 in range(len(sequence)):
                index2 = index1 + g + 1
                if index1 < len(sequence) and index2 < len(sequence) and sequence[index1] in AA and sequence[index2] in AA:
                    myDict[sequence[index1] + sequence[index2]] = myDict[sequence[index1] + sequence[index2]] + 1
                    sum = sum + 1
            for pair in aaPairs:
                code.append(myDict[pair] / sum)
    code = numpy.array(code).reshape(-1,len(aaPairs))
    return code

def AAINDEX(sequence):
    #FINISHED
    AA = 'ARNDCQEGHILKMFPSTWYV'
    encodings = []
    fileAAindex = getdatadir()+'/data/AAindex.txt'
    with open(fileAAindex) as f:
        records = f.readlines()[1:]
    AAindex = []
    AAindexName = []
    code = []
    for i in records:
        AAindex.append(i.rstrip().split()[1:] if i.rstrip() != '' else None)
        AAindexName.append(i.rstrip().split()[0] if i.rstrip() != '' else None)
    single_length = len(AAindex)
    index = {}
    for i in range(len(AA)):
        index[AA[i]] = i
    for aa in sequence:
        if aa == '-':
            for j in AAindex:
                code.append(0)
            continue
        for j in AAindex:
            code.append(float(j[index[aa]]))
    code = numpy.array(code).reshape(len(sequence),len(AAindex))
    return code