'''
Created on Sep 25, 2014

@author: xiong_yang

Description:
    Sequence and position generator
'''
import random

def gen_seq(size, pool='ACGT'):
    '''Generating a sequence with length of size and randomly choose character
     from character pool'''
    seq = ''
    for i in range(size):
        seq += random.sample(pool, 1)[0]
    return seq

def gen_pos(size, length):
    '''Randomly generating a series of positions with length of size within
    interval from [0, length], return a sorted list'''
    pos = []
    for i in range(size):
        pos.append(random.uniform(0, length))
    pos.sort()
    return pos