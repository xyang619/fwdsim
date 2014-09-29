'''
Created on Sep 25, 2014

@author: xiong_yang
'''
import getopt, sys
from model import *
from reader import read_hap
    
def help():
    s = '''Description:
    Perform forward time simulation, two way admixture, in which
    implemented HI, CGF and GA model
Arguments:
    -h      print help message
    -f      simulation started haplotype in MS format <string>
    -s      generation of ancestral populations split   [int=2000]
    -p      split proportion of population 1    [float=0.5]
    -1      effective population size of ancestral population 1 [int=1000]
    -2      effective population size of ancestral population 2 [int=1000]
    -g      generation when admixture start [int=10]
    -m      proportion of admixture from ancestral population 1 [float=0.5]
    -d      admixture model, should be hi|ga|cgf    [string=hi]
    -3      effective population size of admixed population [int=1000]
    -n      number of haplotype to be sampled. npop1,npop2,nadm [int=50,50,50]
    -l      length of sequence in the unit of Morgan [float=1.0]
    -o      output file name [string=sim_${d}_${g}.txt]
    '''
    print(s)
    
def main():
    f = ''      #input file
    gs = 2000   #generation since two ancestral population split
    ps = 0.5    #proportion of samples split into ancestral population 1
    n1 = 1000   #effective population size of population 1
    n2 = 1000   #effective population size of population 2
    ga = 10     #generation since admixture event occured
    pa = 0.5    #proportion of ancestral contribution of population 1
    ma = 'hi'   #admixture model
    na = 1000   #effective population size of admixed population
    ns = []     #number of haplotypes sampled at the end of simulation
    l = 1.0     #the length of simulated sequence in Morgan
    of = ''     #output file
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hf:s:p:1:2:g:m:d:3:n:o:')
    except getopt.GetoptError as err:
        print(err)
        help()
        sys.exit(2)
    for o, a in opts:
        if o == '-h':
            help()
            sys.exit(1)
        elif o == '-f':
            f = a
        elif o == '-o':
            of = a
        elif o == '-s':
            gs = int(a)
        elif o == '-p':
            ps = float(a)
        elif o == '-1':
            n1 = int(a)
        elif o == '-2':
            n2 = int(a)
        elif o == '-g':
            ga = int(a)
        elif o == '-m':
            pa = float(a)
        elif o == '-d':
            if a not in ('hi', 'HI', 'ga', 'GA', 'cgf', 'CGF'):
                print('Warning! Model must be hi, ga or cgf')
                sys.exit(1)
            ma = a
        elif o == '-3':
            na = int(a)
        elif o == '-l':
            l = float(a)
        elif o == '-n':
            lst = a.split(',')
            if len(lst) < 3:
                print('Warning! Sample number set to default')
            else:
                for ele in lst:
                    ns.append(int(ele))
    if len(f) < 1:
        print('Error! Input file must be specified')
        help()
        sys.exit(1)
    if len(ns) < 3:
        print('Warning! Sample number set to default')
        ns = [50, 50, 50]
    if len(of) < 1:
        of = 'sim_{}_{}.txt'.format(ma, ga)
    positions, haps = read_hap(f, l) 
    pop1, pop2 = split(haps, ps)
    for i in range(len(pop2)):
        pop2[i] = pop2[i].lower()
    ip1 = evolve(pop1, positions, n1, gs - ga)
    ip2 = evolve(pop2, positions, n2, gs - ga)
    fp = None
    if ma in ('hi', 'HI'):
        fp = hi_model(ip1, n1, ip2, n2, positions, na, pa, ga)
    elif ma in ('ga', 'GA'):
        fp = ga_model(ip1, n1, ip2, n2, positions, na, pa, ga)
    elif ma in ('cgf', 'CGF'):
        fp = cgf_model(ip1, n1, ip2, n2, positions, na, pa, ga)
    sp = []  # numbers of samples
    i = 0
    for n in ns:
        sp.extend(sample(fp[i], n))
        i += 1
    with open(of, 'w') as out:
        out.write('split {} {} adm {} {} {} Ne {} {} {} sample {} {} {}\n'.format(gs, ps, ga, pa, ma, n1, n2, na, ns[0], ns[1], ns[2]))
        from time import gmtime, strftime
        out.write('{}\n'.format(strftime("%Y-%m-%d %H:%M:%S", gmtime())))
        out.write('\n')
        out.write('//\n')
        out.write('segsites: {}\n'.format(len(positions)))
        posline = 'positions:'
        for p in positions:
            posline += ' {:.8f}'.format(p)
        out.write(posline + '\n')
        for hap in sp:
            out.write(hap + '\n')

if __name__ == '__main__':
    main()
