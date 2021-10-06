import re, sys, os, platform,warnings
import math
from collections import Counter
from bioencoder.utils import getdatadir,Rvalue
def Binary(sequence):
    AA = 'ARNDCQEGHILKMFPSTWYV'
    code=list()
    for aa in sequence:
        if aa == '-':
            code = code + [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            continue
        for aa1 in AA:
            tag = 1 if aa == aa1 else 0
            code.append(tag)
    return code

def EAAC(sequence,window=5):
    if(len(sys.argv)==1):
        warnings.warn("Not set window size,default set to 5")
    AA = 'ARNDCQEGHILKMFPSTWYV'
    code=list()
    for j in range(len(sequence)):
        if j < len(sequence) and j + window <= len(sequence):
            count = Counter(re.sub('-', '', sequence[j:j+window]))
            for key in count:
                count[key] = count[key] / len(re.sub('-', '', sequence[j:j+window]))
            for aa in AA:
                code.append(count[aa])
    return code

def GAAC(sequence):
    group = {
        'alphatic': 'GAVLMI',
        'aromatic': 'FYW',
        'postivecharge': 'KRH',
        'negativecharge': 'DE',
        'uncharge': 'STCPNQ'
    }

    groupKey = group.keys()
    code = []
    count = Counter(sequence)
    myDict = {}
    for key in groupKey:
        for aa in group[key]:
            myDict[key] = myDict.get(key, 0) + count[aa]
        for key in groupKey:
            try:
                code.append(myDict[key]/len(sequence))
            except:
                code.append(0)
    return code

def CKSAAP(sequence, gap=5):
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
    return code

def AAINDEX(sequence):
    AA = 'ARNDCQEGHILKMFPSTWYV'
    encodings = []
    fileAAindex = getdatadir()+'/data/AAindex.txt'
    print(fileAAindex)
    with open(fileAAindex) as f:
        records = f.readlines()[1:]
    AAindex = []
    AAindexName = []
    code = []
    for i in records:
        AAindex.append(i.rstrip().split()[1:] if i.rstrip() != '' else None)
        AAindexName.append(i.rstrip().split()[0] if i.rstrip() != '' else None)

    index = {}
    for i in range(len(AA)):
        index[AA[i]] = i
    for aa in sequence:
        if aa == '-':
            for j in AAindex:
                code.append(0)
                continue
        for j in AAindex:
            code.append(j[index[aa]])
    return code

def PAAC(sequence,lambdaValue=30, w=0.05):
    dataFile = getdatadir()+'/data/PAAC.txt'
    with open(dataFile) as f:
        records = f.readlines()
    AA = ''.join(records[0].rstrip().split()[1:])
    AADict = {}
    for i in range(len(AA)):
        AADict[AA[i]] = i
    AAProperty = []
    AAPropertyNames = []
    for i in range(1, len(records)):
        array = records[i].rstrip().split() if records[i].rstrip() != '' else None
        AAProperty.append([float(j) for j in array[1:]])
        AAPropertyNames.append(array[0])
    AAProperty1 = []
    for i in AAProperty:
        meanI = sum(i) / 20
        fenmu = math.sqrt(sum([(j-meanI)**2 for j in i])/20)
        AAProperty1.append([(j-meanI)/fenmu for j in i])
    code = []
    theta = []
    for n in range(1, lambdaValue + 1):
        theta.append(sum([Rvalue(sequence[j], sequence[j + n], AADict, AAProperty1) for j in range(len(sequence) - n)]) / (len(sequence) - n))
        myDict = {}
        for aa in AA:
            myDict[aa] = sequence.count(aa)
        code = code + [myDict[aa] / (1 + w * sum(theta)) for aa in AA]
        code = code + [(w * j) / (1 + w * sum(theta)) for j in theta]
    return code