#

rule hatchet_binBAM:
    input:
        normal = normal_bam,
        tumors = tumor_bams,
        ref = config['reference']
    output: 
        norm_bed = config['hatchet']['xdir'] + "/{patient}/rdr/normal.1bed",
        tumor_bed = config['hatchet']['xdir'] + "/{patient}/rdr/tumor.1bed",
        total = config['hatchet']['xdir'] + "/{patient}/rdr/total.tsv",
        bins_log = config['hatchet']['xdir'] + "/{patient}/rdr/bins.log"
        #norm_bed =  "{path}/{patient}/rdr/normal.1bed",
        #tumor_bed = "{path}/{patient}/rdr/tumor.1bed",
        #total = "{path}/{patient}/rdr/total.tsv",
        #bins_log = "{path}/{patient}/bins.log"
    params:
        tumors_string = get_tumors_hatchet, # space separated list of bam files
        all_names = get_names_hatchet, # space separated list of Normal and tumors
        bin = config['hatchet']['bin']
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 30000
    threads: 8
    conda:
        "../envs/hatchet.yml"
    shell:
        #"export GRB_LICENSE_FILE=\"/u/bjarnold//gurobi.lic\"\n"
        #"cd {input.xdir}\n"
        "python3 -m hatchet binBAM "
        "-N {input.normal} -T {params.tumors_string} -S {params.all_names} -b {params.bin} -g {input.ref} "
        "-j {threads} -O {output.norm_bed} -o {output.tumor_bed} -t {output.total} "
        "|& tee {output.bins_log}"

