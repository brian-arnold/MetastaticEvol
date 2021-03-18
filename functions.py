import os
from collections import defaultdict
from snakemake import load_configfile

config = load_configfile('config.yaml')

# This file contains (1) functions used by the snakemake rules
# and (2) global variables accessed by these functions, e.g.
# the dictionaries that organize samples by patient

##################
# FUNCTIONS
##################

def get_samples_by_patient():
    tumors = defaultdict(dict)
    normals = defaultdict()
    with open(config["samples"], "r") as f:
        for l in f:
            if not l.startswith("PATIENT_ID"):
                info = l.strip().split() 
                patient = info[0]
                bam = info[2]
                if info[1] == "tumor":
                    sample = get_sample_name( bam )
                    tumors[patient][sample] = bam
                elif info[1] == "normal":
                    normals[patient] = bam
                else:
                    sys.exit("incorrect samples.txt file")
        return (tumors, normals)

def make_all_input_hatchet():
    # output is 1 vcf per tumor sample per patient, 
    # with samples grouped in directories according to patient names
    patients = set()
    paths = []
    with open(config["samples"], "r") as f:
        for l in f:
            if not l.startswith("PATIENT_ID"):
                info = l.strip().split() 
                if info[1] == "tumor":
                    patients.add( info[0] )

        for p in patients:
            paths.append( os.path.join( config["hatchet"]["xdir"], p, "rdr/tumor.1bed") )
        return paths

def make_all_input():
    # output is 1 vcf per tumor sample per patient, 
    # with samples grouped in directories according to patient names
    with open(config["samples"], "r") as f:
        outfiles = []
        for l in f:
            if not l.startswith("PATIENT_ID"):
                info = l.strip().split() 
                if info[1] == "tumor":
                    sample = get_sample_name( info[2] )
                    path = os.path.join( "varscan", info[0], f"{sample}.snp.vcf")
                    outfiles.append( path ) 
        return outfiles

def get_sample_name(x):
    sample = os.path.basename(x)
    sample = sample.replace(".bam", "") 
    return sample
        
def tumor_bam(wildcards):
    # tumors dict defined below
    tumor_in = tumors[wildcards.patient][wildcards.sample]
    return tumor_in

def tumor_bams(wildcards):
    # tumors dict defined below
    tumors_in = [ tumors[wildcards.patient][s] for s in tumors[wildcards.patient] ]
    return tumors_in

def normal_bam(wildcards):
    # normals dict defined below
    normal_in = normals[wildcards.patient]
    return normal_in

def get_tumors_hatchet(wildcards):
    # normals dict defined below
    tumors_in = [ tumors[wildcards.patient][s] for s in tumors[wildcards.patient] ]
    names = " ".join(tumors_in)
    return names 

def get_names_hatchet(wildcards):
    # normals dict defined below
    tumors_in = [ get_sample_name(t) for t in tumors[wildcards.patient] ]
    names = ["Normal"] + tumors_in
    return names 

##################
# GLOBALS
##################

# Data structures for grouping samples by patient:
# tumors[patient][sample] = bam
# normals[patient] = bam
tumors, normals = get_samples_by_patient() 




