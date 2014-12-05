#!/usr/bin/env python

import sys
        
def main(argv):
    blossum = {}
    with open("blosum.txt", 'r') as f:
            lines = f.read().splitlines()
    columns = lines.pop(0).split()
    for line in lines:
        values = line.split()
        key = values.pop(0)
        for i in range(len(columns)):
            blossum[key+columns[i]] = int(values[i])
        
    print blossum
    
if __name__ == '__main__':
    main(sys.argv)
