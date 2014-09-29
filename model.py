'''
Created on Sep 25, 2014

@author: xiong_yang
Description:
    Implemented HI, GA and CGF admixture model    
'''
import random

from recombinator import recombine

def split(population, prop):
    '''The population split into two populations, with proportion prop split into
    population 1'''
    npop1 = int(len(population) * prop)
    npop2 = len(population) - npop1
    pop1 = []
    pop2 = []
    for i in range(npop1):
        pop1.extend(random.sample(population, 1))
    for i in range(npop2):
        pop2.extend(random.sample(population, 1))
    return pop1, pop2

def sample(population, nsample):
    '''Randomly sample nsample of haplotypes from population
    population: a list of haplotypes
    nsample:    an integer value'''
    samples = []
    for i in range(nsample):
        samples.extend(random.sample(population, 1))
    return samples
  
def forward_gen(population, positions, Ne):
    '''Randomly pair two haplotypes of the 2Ne haplotypes in population, and 
    recombine the pairedhaplotypes to generate population of next generation
    Ne: an integer value'''
    pops_next = []  # population of next generation
    n_pop = len(population)
    sample_pool = [i for i in range(n_pop)]
    for i in range(Ne):
        indexes = random.sample(sample_pool, 2)
        haps1 = population[indexes[0]]
        haps2 = population[indexes[1]]
        haps_next_gen = recombine(haps1, haps2, positions)
        pops_next.extend(haps_next_gen)
        sample_pool.remove(indexes[0])
        sample_pool.remove(indexes[1])
    return pops_next

def evolve(population, positions, Ne, gen=1):
    '''population envolve forward gen generations, with effective population 
    size remain the same
    gen:    an integer value'''
    pops_new = sample(population, 2 * Ne)
    for g in range(gen):
        pops_new = forward_gen(pops_new, positions, Ne)
    return pops_new

def hi_model(pop1, Ne1, pop2, Ne2, positions, Ne_adm, prop, gen):
    '''HI admixture model
    pop1, pop2: a list of haplotypes
    Ne1, Ne2, Ne_adm, gen: an integer value
    prop: a float value'''
    pops_adm = []
    npop1 = int(Ne_adm * prop)  # number of individual from population 1
    pops_adm.extend(sample(pop1, 2 * npop1))  # add the haplotype sampled from pop1
    pops_adm.extend(sample(pop2, 2 * (Ne_adm - npop1)))  # haplotype sampeld from pop2
    pops_adm = evolve(pops_adm, positions, Ne_adm, gen)  # evolove independently
    pop1 = evolve (pop1, positions, Ne1, gen)  # evolve of pop1
    pop2 = evolve (pop2, positions, Ne2, gen)  # evolve of pop2
    return pop1, pop2, pops_adm
    
def cgf_model(pop1, Ne1, pop2, Ne2, positions, Ne_adm, prop, gen):
    '''each generation receive alpha ancestry contribution from population 2,
    total contribution of pop1 is m, then alpha = 1- m**(1/gen)'''
    alpha = 1 - prop ** (1.0 / gen)
    pops_adm = hi_model(pop1, Ne1, pop2, Ne2, positions, Ne_adm, 1 - alpha, 1)[2]
    pop2 = evolve (pop2, positions, Ne2, gen)
    for g in range(gen - 1):
        pops_adm = hi_model(pops_adm, Ne_adm, pop2, Ne2, positions, Ne_adm, 1 - alpha , 1)[2]
        pop2 = evolve (pop2, positions, Ne2, gen)
    pop1 = evolve (pop1, positions, Ne1, gen)
    return pop1, pop2, pops_adm
    
def ga_model(pop1, Ne1, pop2, Ne2, positions, Ne_adm, prop, gen):
    '''each generation receive alpha and beta ancestry contributions from both
    parental populations, respectively.'''
    # First generation of admixed population
    pops_adm = []   
    npop1 = int(Ne_adm * prop)
    pops_adm.extend(sample(pop1, 2 * npop1))    
    pops_adm.extend(sample(pop2, 2 * (Ne_adm - npop1))) 
    pops_adm = evolve(pops_adm, positions, Ne_adm, 1)   
    pop1 = evolve (pop1, positions, Ne1, 1)
    pop2 = evolve (pop2, positions, Ne2, 1)
    # The following generations
    alpha = prop / gen
    beta = (1 - prop) / gen
    for g in range(gen - 1):
        tmp_haplo = []
        npop1 = int(Ne_adm * alpha)
        npop2 = int(Ne_adm * beta)
        nprev = Ne_adm - npop1 - npop2
        tmp_haplo.extend(sample(pop1, 2 * npop1))
        tmp_haplo.extend(sample(pop2, 2 * npop2))
        tmp_haplo.extend(sample(pops_adm, 2 * nprev))
        pops_adm = evolve(pops_adm, positions, Ne_adm, 1)
        pop1 = evolve (pop1, positions, Ne1, 1)
        pop2 = evolve (pop2, positions, Ne2, 1)
    return pop1, pop2, pops_adm