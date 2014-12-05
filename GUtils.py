#!/usr/bin/env python

import sys, re
from string import maketrans

def read_fasta(filename):
    read_list = []
    with open(filename, 'r') as f:
        lines = f.read().splitlines()
    curr_read = ""
    x = lines.pop(0)
    for line in lines:
        if line[0] == ">":
            read_list.append(curr_read)
            curr_read = ""
        else:
            curr_read += line
    read_list.append(curr_read)
    return read_list

def generate_proteins(read):
    rstring = string[::-1]
    intab = 'ACTG'
    outtab = 'TGAC'
    trantab = maketrans(intab,outtab)
    rstring = rstring.translate(trantab)

    string = re.sub('T','U', string)
    rstring = re.sub('T','U', rstring)
    
    codontab = {"UUU": "F", "UUC": "F", "UUA":"L", "UUG":"L", "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S", "UAU":"Y", "UAC":"Y",
                "UAA":-1, "UAG":-1, "UGU":"C", "UGC":"C", "UGA":-1, "UGG":"W", "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L", "CCU":"P",
                "CCC":"P", "CCA":"P", "CCG":"P", "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q", "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
                "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M", "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T", "AAU":"N", "AAC":"N", "AAA":"K",
                "AAG":"K", "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R", "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V", "GCU":"A", "GCC":"A",
                "GCA":"A", "GCG":"A", "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E", "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"}
    
    proteins = set()
    
    def extract_protein(rna_str):
        start_indecies = [m.start() for m in re.finditer('AUG',rna_str)]
        #print start_indecies
        for st in start_indecies:
            prot_str = "M"
            for i in range(st+3, len(rna_str),3):
                #print prot_str
                if i+3 >= len(rna_str):
                    break
                cond = codontab.get(rna_str[i:i+3])
                if cond == -1:
                    proteins.add(prot_str)
                    break
                else:
                    prot_str += cond
    extract_protein(rstring)
    extract_protein(string)
    return proteins
