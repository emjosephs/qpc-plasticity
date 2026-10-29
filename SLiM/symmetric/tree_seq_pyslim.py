#!/usr/bin/env python
# coding: utf-8

import msprime, pyslim, tskit, sys, ast, os, re, time, numpy as np, argparse
os.chdir("/mnt/research/josephslab/Husain/cond_Qpc/symmetric/sim_output/trees/")
pyslim.warnings.simplefilter('ignore', msprime.TimeUnitsMismatchWarning)

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--seed", required = True)
parser.add_argument("-r", "--recom", required = True, type=float)
parser.add_argument("-N", "--NTot", required = True, type=int)
parser.add_argument("-mu", "--mutate", required = True, type=float)
parser.add_argument("-o", "--output", required = True, type=str)

args = vars(parser.parse_args())
seed = args['seed']
r = args['recom']
N = args['NTot']
mu = args['mutate']
output = args['output']

ts = tskit.load(seed + ".trees")

print(f"The tree sequence has {ts.num_trees} trees on a genome of length {ts.sequence_length},"
	f" {ts.num_individuals} individuals, {ts.num_samples} 'sample' genomes,"
	f" and {ts.num_mutations} mutations.")

metapopSide=10
demography = msprime.Demography.from_tree_sequence(ts)
for pop in demography.populations:
	# must set their effective population sizes
	pop.initial_size = 50

for x in range(1, metapopSide+1):
	for y in range(1, metapopSide+1):
		destID = x + (y-1) * metapopSide
		destPop = "".join(["p",str(destID)])
		if x>1:
			sourcePop="".join(["p",str(destID-1)])
			demography.set_migration_rate(rate=0.05,source=sourcePop,dest=destPop)
		if x<metapopSide:
			sourcePop="".join(["p",str(destID+1)])
			demography.set_migration_rate(rate=0.05,source=sourcePop,dest=destPop)
		if y>1:
			sourcePop="".join(["p",str(destID-metapopSide)])
			demography.set_migration_rate(rate=0.05,source=sourcePop,dest=destPop)
		if y<metapopSide:
			sourcePop="".join(["p",str(destID+metapopSide)])
			demography.set_migration_rate(rate=0.05,source=sourcePop,dest=destPop)


## BUILD RECOMBINATION MAP
C_num = 10
C_lengths = [100000] * C_num
rate_within = r
rate_between = 0.5

positions = [0]
rates = []
running = 0
for i in range(C_num):
    running += C_lengths[i]
    positions.append(running)
    rates.append(rate_within)
    if i < C_num - 1:
        positions.append(running + 1)
        rates.append(rate_between)

positions[-1] += 1  # end at last position + 1

rate_map = msprime.RateMap(position=positions, rate=rates)

tsRecap = pyslim.recapitate(
	ts,
	demography=demography,
	recombination_rate=rate_map,
	random_seed=seed
)

next_id = pyslim.next_slim_mutation_id(tsRecap)
mutated = msprime.sim_mutations(
	tsRecap,
	rate = mu,
	random_seed=seed,
	model=msprime.SLiMMutationModel(type=2, next_id=next_id),
	keep=True
)

print(f"The tree sequence now has {mutated.num_trees} trees,"
	f" and {mutated.num_mutations} mutations.")

# mutatedNuc = pyslim.generate_nucleotides(mutated)
# mutatedNuc = pyslim.convert_alleles(mutatedNuc)

with open("/mnt/research/josephslab/Husain/cond_Qpc/symmetric/sim_output/vcf/" + seed + "_mutated.vcf", "w") as vcf_file:
	mutated.write_vcf(vcf_file, allow_position_zero=True)
