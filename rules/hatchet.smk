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

rule hatchet_SNPCaller:
    # Note this uses bcftools to read a vcf url, which generates a .tbi file in the parent dir
    input:
        normal = normal_bam,
        ref = config['reference']
    output: 
        hatchet_SNPCaller_snps( "{patient}" ), # can't use functions as output? at least ones like input based on wilcards
        bafs_log = config['hatchet']['xdir'] + "/{patient}/baf/bafs.log"
        #bafs_log = config['hatchet']['xdir'] + "/{patient}/baf/bafs.log"
        #snps = directory(config['hatchet']['xdir'] + "/{patient}/snps"),
    params:
        minreads = config['hatchet']['minreads'],
        maxreads = config['hatchet']['maxreads'],
        list = get_hatchet_list(),
        snpdir = config['hatchet']['xdir'] + "/{patient}/snps",
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 30000
    threads: 8
    conda:
        "../envs/hatchet.yml"
    shell:
        #"export GRB_LICENSE_FILE=\"/u/bjarnold//gurobi.lic\"\n"
        "mkdir -p {params.snpdir}\n"
        "python3 -m hatchet SNPCaller "
        "-N {input.normal} -r {input.ref} -j {threads} "
        "-c {params.minreads} -C {params.maxreads} -R {params.list} "
        "-o {params.snpdir} |& tee {output.bafs_log}"

rule hatchet_deBAF:
    input:
        snps = hatchet_SNPCaller_snps( "{patient}" ), # not using wildcards here since function already exists for output above
        normal = normal_bam,
        tumors = tumor_bams,
        ref = config['reference']
    output: 
        norm_bed = config['hatchet']['xdir'] + "/{patient}/baf/normal.1bed",
        tumor_bed = config['hatchet']['xdir'] + "/{patient}/baf/tumor.1bed",
        bafs_log = config['hatchet']['xdir'] + "/{patient}/baf/bafs.log"
    params:
        minreads = config['hatchet']['minreads'],
        maxreads = config['hatchet']['maxreads'],
        tumors_string = get_tumors_hatchet, # space separated list of bam files
        all_names = get_names_hatchet, # space separated list of Normal and tumors
        snpdir = config['hatchet']['xdir'] + "/{patient}/snps/\*.vcf.gz",
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 30000
    threads: 8
    conda:
        "../envs/hatchet.yml"
    shell:
        #"export GRB_LICENSE_FILE=\"/u/bjarnold//gurobi.lic\"\n"
        #"cd {input.xdir}\n"
        "python3 -m hatchet deBAF "
        "-N {input.normal} -T {params.tumors_string} -S {params.all_names} -r {input.ref} "
        "-j {threads} -c {params.minreads} -C {params.maxreads} -L {input.snps} "
        "-O {output.norm_bed} -o {output.tumor_bed} |& tee {output.bafs_log}"
