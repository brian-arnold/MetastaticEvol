#

rule mpileup_tumors:
    input:
        tumor = tumor_bam,
        ref = config['reference']
    output: 
        mpile = "{patient}/{sample}.mpileup" 
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 3000
    conda:
        "../envs/varscan.yml"
    shell:
        "samtools mpileup -f {input.ref} -q 1 {input.tumor} > {output}"


rule mpileup_normals:
    input:
        normal = normal_bam,
        ref = config['reference']
    output: 
        mpile = "{patient}/normal.mpileup" 
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 3000
    conda:
        "../envs/varscan.yml"
    shell:
        "samtools mpileup -f {input.ref} -q 1 {input.normal} > {output}"

rule varscan:
    input:
        tumor_mpile = "{patient}/{sample}.mpileup",
        normal_mpile = "{patient}/normal.mpileup" 
    output: 
        vcf = "{patient}/{sample}.vcf" 
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 3000
    conda:
        "../envs/varscan.yml"
    shell:
        "varscan somatic {input.normal_mpile} {input.tumor_mpile} {output} --output-vcf 1"
