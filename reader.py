'''
Created on Sep 25, 2014

@author: xiong_yang
Description:
    Read haplotype and positions from MS simulation output
'''
import random

def read_hap(simfile, length):
    '''Read positions and haplotypes from the output of ms simulation'''
    positions = []
    haps = []
    alleles = []
    with open(simfile) as f:
        for i in range(5):
            f.readline()
        posline = f.readline()
        for p in posline.split()[1:]:
            positions.append(float(p) * length)
            allele = random.sample('ACGT', 2)
            alleles.append(allele)
        for line in f:
            hap = ''
            i = 0
            for char in line[:-1]:
                hap += alleles[i][int(char)]
                i += 1
            haps.append(hap)
    i = 0
    with open('{}.pos'.format(simfile), 'w') as out:
        for p in positions:
            out.write('{}\t{}\t{}\n'.format(p, alleles[i][0], alleles[i][1]))
            i += 1
    return positions, haps

