#!/usr/bin/env python

import pickle, random, Queue

class ProteinIdentity:
    def __init__(self, protein, poss, curr_prot, curr_poss, score):
        self.protein = protein
        self.poss = poss
        self.curr_poss = curr_poss
        self.score = score
        self.curr_prot = curr_prot
        self.left_end = False
        self.right_end = False
    
    def __lt__(self, other):
        return self.score < other.score
    
    def __update_score(self):
        choice = []
        if self.left_end:
            choice.append[0]
        if self.right_end:
            choice.append[1]
        
        if not choice:
            #this is our boy
            return True
        
        x = random.choice(choice)
        if x == 0:
            score += edit(self.protein[self.curr_prot[0]],self.poss[self.curr_poss[0]])
            if self.curr_prot[0] > 0:
                self.curr_prot[0] -= 1
            if self.curr_poss[0] > 0:
                self.curr_poss[0] -= 1
        else:
            score += edit(self.protein[self.curr_prot[1]],self.poss[self.curr_poss[1]])
            if self.curr_prot[1] < len(self.protein) - 1:
                self.curr_prot[1] += 1
            if self.curr_poss[1] < len(self.poss) - 1:
                self.curr_poss[1] += 1
            
    def __at_left_end(self):
        if self.curr_prot[0] == 0 and self.curr_poss[0] == 0:
            self.left_end = True

    def __at_right_end(self):
        if self.curr_prot[1] == len(self.protein) - 1 and self.curr_poss[1] == len(self.poss) - 1:
            self.right_end = True

    def __substitution(self, first, second):
        

class PQueue:
 
    pqueue = Queue.PriorityQueue()

    def __init__(self, read_proteins, scheme):

        x = pickle.load(open(scheme,'r'))
        for prot in read_proteins:
            temp = "$" + prot + "#"
            prelim_score = x.score_sequence(temp)
        self.blossum = {}
        with open("blosum.txt", 'r') as f:
            lines = f.read().splitlines()
        columns = lines.pop(0).split()
        for line in lines:
            values = line.split()
            key = values.pop(0)
            for i in range(len(columns)):
                self.blossum[key+columns[i]] = int(values[i])
   
    def __substitution(self, first, second):
        return -self.blossum[first+second] #make it negative as lower values in PQueue are first
