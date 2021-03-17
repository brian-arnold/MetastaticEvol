#




#### REMEBER TO INCLUDE IN SNAKEFILE!!!!!!!!!!
# REMEMBER TO PUT GUROBI LICENSE SOMEWHERE????
#
#
#
#
rule hatchet_binBAM:
    input:
        tumor = tumor_bam,
        ref = config['reference']
    output: 
        mpile = "varscan/{patient}/{sample}.mpileup" 
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 3000
    conda:
        "../envs/hatchet.yml"
    shell:
        "samtools mpileup -f {input.ref} -q 1 {input.tumor} > {output}"

