#

from functions import *
configfile: "config.yaml"

rule all:
    input:
        # custom func instead of expand with wildcards 
        # bc patients may have variable/unequal samples
        make_all_input()

include: "rules/varscan.smk"

