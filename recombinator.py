'''
Created on Sep 25, 2014

@author: xiong_yang
Description:
    Simulate the recombination process of two haplotypes in germline
    The recombination events occur follow Poisson process
'''

import random, math

def wait_time(lamb=1):
    '''Waiting for the first occurance follow distribution
    P(T<t) = 1 - exp(-lambda*t)
    '''
    prob = random.random()
    time = -math.log(1 - prob) / lamb
    return time

def break_points(length):
    '''Generating a set of break points along the length of time interval,
    following Poisson process'''
    bps = [0]
    while bps[-1] <= length:
        bps.append(bps[-1] + wait_time(length))
    return bps

def find_index(pos, positions):
    '''Find the index of position best match of given position in positions
    using binary search
    pos:    float value
    positions:  a list of float value'''
    length = len(positions) 
    if pos < positions[0]:
        return 0
    if pos > positions[-1]:
        return length
    left = 0
    right = length
    mid = (left + right + 1) / 2
    while right > left:
        if pos > positions[mid]:
            left = mid
        else:
            right = mid - 1
        mid = (left + right + 1) / 2
    if (left < length) and (abs(pos - positions[left]) > abs(positions[left + 1] - pos)):
        left = left + 1
    return left
                 
def recombine(haplo1, haplo2, positions):
    '''two haplotype recombinate in one generation. the positions are in Morgan.
    The recombinations occurred following Poisson process
    haplo1, haplo2: a character sequence
    positions:  a list of float value''' 
    length = positions[-1]
    bps = break_points(length)
    haps = ['', '']
    indexes = []
    for bp in bps:
        indexes.append(find_index(bp, positions))
    start = find_index(bps[0], positions) 
    copy1 = True
    for bp in bps[1:]:
        end = find_index(bp, positions)
        if copy1:
            haps[0] += haplo1[start:end]
            haps[1] += haplo2[start:end]
            copy1 = False
        else:
            haps[0] += haplo2[start:end]
            haps[1] += haplo1[start:end]
            copy1 = True
        start = end
    return haps
