#

from functions import *
import random
configfile: "config.yaml"

rule all:
    input:
        # custom func instead of expand with wildcards 
        # bc patients may have variable/unequal samples
        make_all_input_hatchet()

#include: "rules/varscan.smk"
include: "rules/hatchet.smk"

