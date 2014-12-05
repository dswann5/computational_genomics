#!/usr/bin/env python

import sys, argparse, fileinput, math, collections, re, GUtils
from collections import defaultdict
#from GUtils import GUtils

class FreqScore:
    """
    size - how many amino acids to consider for adjacency
    adjacency - location of "adjacent" amnio acids to look at
    alphabet - the alphabet to look at
    filename - the filename of the protein FASTA file """
    def __init__(self, size, adjacency, alphabet, filename, word_size):
        self.size = size
        self.adj = adjacency
        self.word_dict = defaultdict(list)
        possible_comb = []
        #create possible combinations of alphabet
        if size == 1:
            alphabet += "$#"
            for a in alphabet:
                possible_comb.append(a)
        else:
            for a in alphabet:
                self.__create_comp(size - 1, alphabet, a, possible_comb)
            self.__create_comp(size - 1, alphabet, "$", possible_comb)
        
        self.__create_comp(size, alphabet, "", possible_comb)
        self.counts = {}
        self.indiv_counts = {}
        self.__create_dict(possible_comb, size)

        #creates a list of proteins
        self.protein_list = GUtils.read_fasta(filename)
        #self.protein_list = self.__assemble_proteins(filename)

        #count the proteins probabilities and create word dictionary
        j = 0
        for protein in self.protein_list:
            self.__count_protein(size, adjacency, protein)
            self.__disect_words(protein, j, word_size)
            j += 1
        
        #smoothing protein counts
        for x in self.counts:
            self.counts[x] += 1
            self.indiv_counts[x[:size]] +=1
            
        #create probabilities of transitions
        self.__adapt_counts(size)
        #print protein_list
    
    def __create_comp(self, size, alphabet, poss_str, poss_comp):
        for a in alphabet:
            add_char = poss_str + a
            if size - 1 == 0:
                poss_comp.append(add_char)
            elif size - 1 == 1:
                self.__create_comp(size - 1, alphabet + "#", add_char, poss_comp)
            else:
                self.__create_comp(size - 1, alphabet, add_char, poss_comp)
                
    def __create_dict(self, poss_comp, size):
        counts = {}
        for x in poss_comp:
            if x[size-1] != "#":
                self.indiv_counts[x] = 0
                for y in poss_comp:
                    if y[0] != "$":
                        self.counts[self.__create_key(x,y)] = 0
            
    def __create_key(self, first, second):
        return first + "-" + second
    """
    def __assemble_proteins(self, filename):
        protein_list = []
        with open(filename, 'r') as f:
            lines = f.read().splitlines()
        curr_protein = ""
        x = lines.pop(0)
        for line in lines:
            if line[0] == ">":
                protein_list.append(curr_protein)
                curr_protein = ""
            else:
                curr_protein += line
        protein_list.append(curr_protein)
        return protein_list
    """
    def __count_protein(self, size, adj, protein):
        protein = "$" + protein + "#"
        for i in range(len(protein) - (2*size + adj) +1):
            start = protein[i:i+size]
            end = protein[i+size+adj:i+2*size+adj]
            self.counts[self.__create_key(start,end)] += 1
            self.indiv_counts[start] +=1
    
    def __disect_words(self, protein, index, word_size):
        for i in range(len(protein)-word_size):
            self.word_dict[protein[i:i+word_size]].append((index,i))
        
    
    def print_prot(self):
        print self.counts
        
    def print_index(self):
        print self.word_dict;
        
    def obtain_proteins(self, word):
        y = []
        for (index,start) in self.word_dict[word]:
            y.append((self.protein_list[index], start))
        return y
        
    def __adapt_counts(self,size):
        for x in self.counts:
            self.counts[x] = self.counts[x]/float(self.indiv_counts[x[:size]])
            
    def score_sequence(self,seq):
        score = 0
        for i in range(len(seq) - (2*size + adj) +1):
            start = protein[i:i+self.size]
            end = protein[i+self.size+adj:i+2*self.size+self.adj]
            score += math.log(self.counts[self.__create_key(start,end)])
        return score