#

rule mpileup_tumors:
    input:
        tumor = tumor_bam,
        ref = config['reference']
    output: 
        mpile = "varscan/{patient}/{sample}.mpileup" 
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
        mpile = "varscan/{patient}/normal.mpileup" 
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 3000
    conda:
        "../envs/varscan.yml"
    shell:
        "samtools mpileup -f {input.ref} -q 1 {input.normal} > {output}"

rule varscan:
    input:
        tumor_mpile = "varscan/{patient}/{sample}.mpileup",
        normal_mpile = "varscan/{patient}/normal.mpileup"
    output: 
        vcf = "varscan/{patient}/{sample}.snp.vcf"
    params:
        name = "varscan/{patient}/{sample}"
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 3000
    conda:
        "../envs/varscan.yml"
    shell:
        "varscan somatic {input.normal_mpile} {input.tumor_mpile} {params.name} --output-vcf 1"
