#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.bam = '' // Path to the BAM file
params.bai = '' // Path to the BAI file
params.ref_tar = '' // Path to the reference tarball
params.pon = '' // Path to the panel of normals

workflow {
    
    // Main workflow
    def mutectOut = callMutations(
        params.bam,
        params.bai,
        params.ref_tar,
        params.pon
    )

    countVariants(mutectOut)
}

process callMutations {
    
    input:
    path bam
    path bai
    path ref_tar
    path pon
    
    output:
    path "mutect_out/*.vcf" into mutectVcf
    
    script:
    def bamBaseName = bam.baseName
    """
    tar -xvf ${ref_tar} -C ./ref
    REF_FASTA=$(find ./ref -name '*.fasta' -o -name '*.fa')
    mkdir -p mutect_out
    pbrun mutectcaller --ref $REF_FASTA --in-bam $bam --out-vcf mutect_out/${bamBaseName}_output.vcf --tumor-name TUMOR --normal-name NORMAL --pon $pon --num-gpus 1
    """
}

process countVariants {
    
    input:
    path vcfFile from mutectVcf
    
    output:
    path "variant_counts.txt"
    
    script:
    """
    bcftools view -v snps $vcfFile | grep -v '^#' | wc -l > snv_count.txt
    bcftools view -v indels $vcfFile | grep -v '^#' | wc -l > indel_count.txt
    echo -e "SNVs:\\t$(cat snv_count.txt)\\nIndels:\\t$(cat indel_count.txt)" > variant_counts.txt
    """
}

