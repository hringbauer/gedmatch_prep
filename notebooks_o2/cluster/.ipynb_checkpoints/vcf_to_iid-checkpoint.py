import numpy as np
import os  # For Saving to Folder
import pandas as pd
import os as os
import sys as sys


#####################################################
### Set Parameters
iid = "I14740"

#####################################################
### Run the Main Program

def create_vcf(iid="", ch=1, out_folder="/n/groups/reich/hringbauer/git/gedmatch_prep/output/vcf1000G_o2/"):
    """Create VCF for a single IID.
    Hard_coded Master VCF of Ali.
    Hard"""
    path_mastervcf = f"/n/data1/hms/genetics/reich/1000Genomes/ali/WholeGenomeImputation/imputed_r2/v54.1/chr{ch}.bcf"
    out_vcf = os.path.join(out_folder, f"{iid}.chr{ch}.bcf")
    
    c1 = f"bcftools view -s {iid} {path_mastervcf} > {out_vcf}"
    #c2 = f"bcftools index {out_vcf}"
    print(c1)
    #print(c2)
    os.system(c1)
    #os.system(c2)

    
if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise RuntimeError("Script needs argument (which chromosome to run)")
    
    ### Load Input Parameters
    ch = int(sys.argv[1])  # Get Parameter of python function
    print(f"Running Chromosome: {ch}...")
    
    create_vcf(iid=iid, ch=ch)