import tskit
#pip install tskit
import msprime
#pip install msprime
from IPython.display import SVG
import numpy as np
import time
import matplotlib.pyplot as plt
from collections import Counter

#genealogy simuation
seq_length = 5_000_000 ##1MB
mu_microsat = 1e-5
mu_snp = 1e-8
mu_indel = 1e-9

#r = 10 ** np.random.uniform(-9, -7)
Ttime = np.arange(0, 6e3+700, 600)
transition_time = int(np.random.choice(Ttime))
#transition_time = int(3000.0)
print(f"Simulating with transition_time = {transition_time}")
alpha = 1.3
print(f"Simulating with alpha = {alpha}")
r = 1e-8
popsize = 5e4
print(f"Simulating with popsize = {popsize}")
#to simulate the sequence data, we need to define the proportion of sites of each marker type in the sequence
#need to make sure if the logic is correct
snpProp = 1/3
microsatProp = 1/3
indelProp = 1/3
#choose a mutation model for each mutation marker
EL2 = msprime.EL2(m=0.43, u=0.68, v=0.037, lo = 10, hi = 17, root_distribution=[0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125]) #microsatellites
JC = msprime.JC69() #snps

a = 0.2 #indel Ratio (insertion rate/ deletion rate)
mu = 0.1 #death rate

TKF = msprime.MatrixMutationModel(
    ["100", "101", "102", "103", "104"],
    root_distribution = [1/5, 1/5, 1/5, 1/5, 1/5],
    transition_matrix = [[1-a*mu, a*mu, 0.0, 0.0, 0.0],
                         [mu, 1-(2*a+1)*mu, 2*a*mu, 0.0, 0.0],
                         [0.0, 2*mu, 1-(3*a+2)*mu, 3*a*mu, 0.0],
                         [0.0, 0.0, 3*mu, 1-(4*a+3)*mu, 4*a*mu],
                         [0.0, 0.0, 0.0, 4*mu, 1-4*mu]]
                        )

ts = msprime.sim_ancestry(
    15,
    recombination_rate=r,
    sequence_length=seq_length,
    population_size=popsize,
    model=[
        msprime.BetaCoalescent(alpha=alpha, duration = transition_time),
        msprime.StandardCoalescent(),
    ]
)


#snps = []
#microsats = []
#indels = []
#snps_filtered = []
#microsats_filtered = []
#indels_filtered = []


#place mutation markers on simulated genealogies
ts_microsats = msprime.sim_mutations(ts, rate=mu_microsat, model=EL2)
#microsats.append(ts_microsats.num_mutations)
ts_snps = msprime.sim_mutations(ts, rate=mu_snp, model=JC)
#snps.append(ts_snps.num_mutations)
ts_indels = msprime.sim_mutations(ts, rate=mu_indel, model=TKF)
#indels.append(ts_indels.num_mutations)


#we designate the sites for each marker according to the proportion of sites of each marker
snp_positions = np.random.choice(seq_length, size=int(seq_length * snpProp), replace=False)
remain_positions = np.setdiff1d(np.arange(seq_length), snp_positions)
microsat_positions = np.random.choice(remain_positions, size=int(seq_length * microsatProp), replace=False)
indel_positions = np.setdiff1d(remain_positions, microsat_positions)


#deleting sites where microsats/snps shouldn't be placed (cause they're not the designated sites for each of them)
#microsats_to_remove = [site.id for site in ts_microsats.sites() if int(site.position) not in list(microsat_positions)]
microsat_positions_set = set(microsat_positions)
microsats_to_remove = [site.id for site in ts_microsats.sites() if int(site.position) not in microsat_positions_set]
ts_microsats_filtered = ts_microsats.delete_sites(microsats_to_remove)
#microsats_filtered.append(ts_microsats_filtered.num_mutations)

snp_positions_set = set(snp_positions)  # Convert to set
snps_to_remove = [site.id for site in ts_snps.sites() if int(site.position) not in snp_positions_set]
ts_snps_filtered = ts_snps.delete_sites(snps_to_remove)
#snps_filtered.append(ts_snps_filtered.num_mutations)

indel_positions_set = set(indel_positions)
indels_to_remove = [site.id for site in ts_indels.sites() if int(site.position) not in indel_positions_set]
ts_indels_filtered = ts_indels.delete_sites(indels_to_remove)
#indels_filtered.append(ts_indels_filtered.num_mutations)

#create the vcf files separately and merge them afterwards using bcftools
with open("JC69only_SIM.vcf", "w") as JC69only_vcf:
    ts_snps_filtered.write_vcf(JC69only_vcf, allow_position_zero=True)

with open("EL2_SIM.vcf", "w") as EL2_vcf:
    ts_microsats_filtered.write_vcf(EL2_vcf, allow_position_zero=True)

with open("JC69_SIM.vcf", "w") as JC69_vcf:
    ts_snps_filtered.write_vcf(JC69_vcf, allow_position_zero=True)

with open("TKF_SIM.vcf", "w") as TKF_vcf:
    ts_indels_filtered.write_vcf(TKF_vcf, allow_position_zero=True)


def write_map_file(ts_list, map_filename, scale_cM_per_bp=1e-6):
    sites = []
    for ts in ts_list:
        for site in ts.sites():
            sites.append(int(site.position))

    # Remove duplicates and sort positions
    sites = sorted(set(sites))

    with open(map_filename, "w") as f:
        for idx, pos in enumerate(sites):
            chrom = "1"  # must match VCF CHROM
            marker_id = f"marker{idx+1}"
            genetic_pos_cM = pos * scale_cM_per_bp
            f.write(f"{chrom}\t{marker_id}\t{genetic_pos_cM:.6f}\t{pos}\n")

#write_map_file([ts_snps_filtered, ts_microsats_filtered, ts_indels_filtered], "SIM.map")
write_map_file([ts_snps_filtered], "SNP.map")


#nohup python SIMgen.py > SIMgen_output.log 2>&1 &

#on terminal
#bgzip EL2_SIM.vcf
#bgzip JC69_SIM.vcf
#bgzip TKF_SIM.vcf
#bcftools index EL2_SIM.vcf.gz
#bcftools index JC69_SIM.vcf.gz
#bcftools index TKF_SIM.vcf.gz


#just concatenate since they're coming from same individuals
#bcftools concat EL2_SIM.vcf.gz JC69_SIM.vcf.gz TKF_SIM.vcf.gz -o merged_SIM.vcf
#bcftools sort merged_SIM.vcf -o merged_SIM.vcf #then sort them according to their positions
# SM
#bcftools concat EL2_SIM.vcf.gz JC69_SIM.vcf.gz -o SM_SIM.vcf
#bcftools sort SM_SIM.vcf -o SM_SIM.vcf

# SI
#bcftools concat JC69_SIM.vcf.gz TKF_SIM.vcf.gz -o SI_SIM.vcf
#bcftools sort SI_SIM.vcf -o SI_SIM.vcf
