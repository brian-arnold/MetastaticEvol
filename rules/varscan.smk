
rule varscan:
    input:
        tumor = get_tumor_in,
        normal = get_normal_in
    output: 
        vcf = "{patient}/{sample}.vcf" 
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 3000
    conda:
        "../envs/varscan.yml"
    shell:
        "touch {output} "

