###########################
### Prepare 23andme txts from h5 files
### Harald Ringbauer, Feb 20th 2024

import numpy as np
import pandas as pd
import os as os
import sys as sys
import h5py
import itertools as it
from ancIBD.IO.prepare_h5 import vcf_to_1240K_hdf

###########################
###

def load_23andme_snps(path_23tsv ="/mnt/archgen/users/hringbauer/git/gedmatch_prep/data/23andme/snps.tsv", 
                      ch = 3, output=True):
    """Load 23andme SNPs for a specific chromosome and return dataframe"""

    df = pd.read_csv(path_23tsv, sep="\t", low_memory=False)

    df_ch = df[df["chr"]==str(ch)] # 23andme SNP file for Chr
    if output:
        print(f"Subset to {len(df_ch)} SNPs. ")
        print(f"Unique {len(set(df_ch.pos))} / {len(df_ch)} SNPs. ")
    return df_ch

def get_h5_df(h5_path="", output=True, iid_idx=0):
    """Convert hdf5 at h5_path to Pandas df with gt and gp info.
    Return df"""
    with h5py.File(h5_path, "r") as f: # Load for Sanity Check. See below!
        if output:
            print(list(f["variants"]))
            print(np.shape(f["calldata/GT"]))
        pos = f["variants/POS"][:]
        ref = f["variants/REF"][:]
        alt = f["variants/ALT"][:]
        #af_all = f["variants/RAF"][:]
        gp = f["calldata/GP"][:, :, :]
        gt = f["calldata/GT"][:, :, :]
        #m = f["variants/MAP"][:]
        samples = f["samples"][:]

    gp_max = np.max(gp[:,iid_idx,:], axis=1)
    gt2 = np.sum(gt[:,iid_idx,:], axis=1)

    frac_gp = np.mean(gp_max>0.99, axis=0)
    if output:
        print(f"Fraction maxGP>0.99: {frac_gp:.3f}")
    dft = pd.DataFrame({"pos":pos, "ref":ref.astype("str"), "alt":alt[:,0].astype("str"), "gt2":gt2, "gp":gp_max}) ### Prepare hdf5 dataframe
    return dft

def create_23gts(df, output=True):
    """Add 23gts column to dataframe with gt2 column.
    Return df"""
    gt23 = pd.Series(["--" for _ in range(len(df))]) # Initialize missing GTs
    idx0 = df["gt2"]==0
    
    gt23[idx0] =  df["ref"].str[0][idx0] * 2

    idx2 = df["gt2"]==2
    gt23[idx2] =  df["alt"].str[0][idx2] * 2

    idx1 = df["gt2"]==1
    gt23[idx1] = (df["ref"].str[0] + df["alt"].str[0])[idx1]
    gt23 = gt23.apply(sorted).apply("".join)
    df["gt23"]=gt23#
    
    if output:
        print(f"0 GTs: {np.sum(idx0)}")
        print(f"1 GTs: {np.sum(idx1)}")
        print(f"2 GTs: {np.sum(idx2)}")
    return df

def h5_to23me_df(chs = range(1,23), iid_idx = 0, output=True,
                 h5_path = f"/mnt/archgen/users/hringbauer/brienzi/data/h5.imputed.23andme/v0.2.imputed.",
                 path_23tsv='/mnt/archgen/users/hringbauer/git/gedmatch_prep/data/23andme/snps.tsv'):
    """Create 23andme formated dataframe of genotypes.
    iid_idx: Indivdiual index in h5 file"""
    dfts = []
    for ch in chs:
        h5_path_ch = f"{h5_path}ch{ch}.h5"
        
        df_ch = load_23andme_snps(ch = ch, path_23tsv=path_23tsv, output=output) # Load 23andme SNP set
        dft = get_h5_df(h5_path_ch, output=output, iid_idx=iid_idx)   # Load h5 data
        dfm = pd.merge(dft, df_ch, on="pos") ### Merge the files
        if output:
            print(f"Merged to {len(dfm)}/{len(df_ch)} SNPs")

        dfm = create_23gts(dfm, output=output) ### Create 23and SNP format column
        dfts.append(dfm)
    dft = pd.concat(dfts)
    return dft

def extract_23andme_header(path_23txt="", s="#"):
    """Extract 23andme Header from 23andme .txt file
    Return string of header
    path_23txt: Path to 23andme txt file"""

    with open(path_23txt,"r") as fi:
        h = []
        for ln in fi:
            if ln.startswith(s):
                h.append(ln[:])
            else:
                break
    h = "".join(h)
    return h

def save23_file(df, savepath ="",
                path_23 = "/mnt/archgen/users/hringbauer/git/gedmatch_prep/data/genome_Harald_Ringbauer_v4_Full_20181111230705.txt"):
    """Save 23andme file at savepath. 
    df: Input Genotype df.
    path_23: Where to copy the header from.
    savepath: Where to save to."""
    h = extract_23andme_header(path_23) # 1: Copy header
    dft = df[["rsid", "chr", "pos", "gt23"]] # Extract relevant columns

    ### Save all data
    with open(savepath, 'w') as f:
         f.write(h)
    dft.to_csv(savepath, header=False, mode="a", index=False, sep="\t")
    print(f"Successfully saved 23andmefile to: {savepath}")

def flag_failed_snps(df, maxgp=0.99, flag_to="--", output=True):
    """Flag failed SNPs genotypes.
    Modify df and return"""
    idx_fail = df["gp"]<maxgp
    df.loc[idx_fail, "gt23"] = flag_to
    
    if output:
        print(f"Flagged {np.sum(idx_fail)}/{len(idx_fail)} SNPs with maxGP < {maxgp} to: {flag_to}")
    return df

def generate_23_from_h5(chs =range(1,23), iid_idx=0, output=False, maxgp=0,
        h5_path = f"/mnt/archgen/users/hringbauer/brienzi/data/h5.imputed.23andme/v0.2.imputed.",
        path_23tsv='/mnt/archgen/users/hringbauer/git/gedmatch_prep/data/23andme/snps.tsv',
        savepath ="", flag_to="--"):
    """Generate 23andme text file from 23-imputed hdf5 file.
    Master Function that combines above detailled functions."""  
    ### 1) Create pandas dataframe from 23andme ref file and imputed h5
    dfm = h5_to23me_df(chs = chs, iid_idx=iid_idx, output=output,
             h5_path = h5_path, path_23tsv=path_23tsv)
    ### 2) Flag SNPs
    if maxgp>0:
        dfm = flag_failed_snps(dfm, maxgp=maxgp, flag_to=flag_to)
    ### 3) Save File
    save23_file(dfm, savepath)
    return dfm

##################################################
### Prepare H5s

def create_folder(path_file):
    """Create folder of path for file (if not existing)"""
    path_folder = os.path.dirname(path_file)
    if not os.path.exists(path_folder):
        os.makedirs(path_folder)
    
def get_path_imputed1000g_autorun(iid =""):
    """Get path of imputed VCF (all 1000G variants) of iid
    Return this path"""
    path_imp_vcf = f"/mnt/archgen/Autorun_eager/eager_outputs/TF/{iid[:3]}/{iid}/GTL_output/{iid}_imputed.vcf.gz" 
    assert(os.path.exists(path_imp_vcf))
    return path_imp_vcf

def prep_h5_autorun_23(iid=""):
    """Prepare HDF5 File for Individual IID in autorun.
    Hard Coded MPI EVA Folder Structure.
    Todo: Create VCF and HDF output path if not existing..."""
    path_vcf_imp = get_path_imputed1000g_autorun(iid =iid)
    
    for ch in range(1,23):
        print(f"Running Chromosome {ch}...")
        path_vcf = f"/mnt/archgen/users/hringbauer/git/gedmatch_prep/output/vcfs.23imp45/{iid}/ch{ch}.vcf.gz"
        path_h5 = f"/mnt/archgen/users/hringbauer/git/gedmatch_prep/output/h5.23imp45/{iid}/ch{ch}.h5"
        create_folder(path_vcf)
        create_folder(path_h5)
    
        vcf_to_1240K_hdf(in_vcf_path = path_vcf_imp,
                 path_vcf = path_vcf,
                 path_h5 = path_h5,
                 marker_path = f"/mnt/archgen/users/hringbauer/git/gedmatch_prep/data/23andme45/bcftools_ch{ch}.csv",
                 map_path = "", 
                 af_path = "",
                 col_sample_af = "AF_SAMPLE",
                 ch=ch, 
                 chunk_length=10000)