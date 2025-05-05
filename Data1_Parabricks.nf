#!/usr/bin/env nextflow

// Define parameters
params.tumorFastqDir
params.outdir
params.reference
params.controlFastqDir = null
params.lowMemory = false

// Load paired FASTQ files for tumor samples
Channel
    .fromFilePairs("${params.tumorFastqDir}/*_R{1,2}.fastq.gz", checkIfExists: true)
    .set { tumorFastqPairs }

// Load paired FASTQ files for control samples if provided
if (params.controlFastqDir) {
    Channel
        .fromFilePairs("${params.controlFastqDir}/*_R{1,2}.fastq.gz", checkIfExists: true)
        .set { controlFastqPairs }
} else {
    controlFastqPairs = Channel.empty()
}

// Process: FQ2BAM for tumor samples
process fq2bamTumor {
    tag "${sample_id}"
    input:
        tuple val(sample_id), file(reads1), file(reads2) from tumorFastqPairs

    output:
        file "${sample_id}.bam" into tumorBamFiles

    script:
        """
        # Run fq2bam
        pbrun fq2bam --ref ${params.reference} --in-fq ${reads1} ${reads2} --out-bam ${sample_id}.bam
        """
        if (params.lowMemory) {
            """
            pbrun fq2bam --ref ${params.reference} --in-fq ${reads1} ${reads2} --out-bam ${sample_id}.bam --low-memory
            """
        }
}

// Process: FQ2BAM for control samples
process fq2bamControl {
    tag "${sample_id}"
    input:
        tuple val(sample_id), file(reads1), file(reads2) from controlFastqPairs

    output:
        file "${sample_id}.bam" into controlBamFiles

    script:
        """
        # Run fq2bam
        pbrun fq2bam --ref ${params.reference} --in-fq ${reads1} ${reads2} --out-bam ${sample_id}.bam
        """
        if (params.lowMemory) {
            """
            pbrun fq2bam --ref ${params.reference} --in-fq ${reads1} ${reads2} --out-bam ${sample_id}.bam --low-memory
            """
        }
}

// Process: Create Panel of Normals
process createPon {
    input:
        file(bamFile) from controlBamFiles

    output:
        file "${params.outdir}/pon.pon" into ponFile

    script:
        """
        # Create Panel of Normals
        pbrun prepon --in-pon-file ${params.outdir}/pon.vcf.gz --out-pon ${params.outdir}/pon.pon --in-bams ${controlBamFiles}
        """
}

// Process: Mutectcaller for tumor-only mode
process mutectcallerTumorOnly {
    tag "${sample_id}"
    input:
        file(bamFile) from tumorBamFiles

    output:
        file "${sample_id}.vcf" into tumorVcfFiles

    script:
        """
        # Extract sample name from BAM file
        sample_name=$(samtools view -H ${bamFile} | grep '@RG' | cut -d' ' -f10 | cut -d':' -f3)

        # Run mutectcaller
        pbrun mutectcaller --ref ${params.reference} --tumor-name ${sample_name} --in-tumor-bam ${bamFile} --out-vcf ${sample_id}.vcf
        """
}

// Process: Mutectcaller for tumor-normal mode with PON
process mutectcallerTumorNormalPon {
    tag "${sample_id}"
    input:
        tuple val(sample_id), file(tumorBam), file(normalBam) from tumorBamFiles.combine(controlBamFiles)
        file(ponFile) from ponFile

    output:
        file "${sample_id}.vcf" into tumorNormalPonVcfFiles

    script:
        """
        # Extract sample names from BAM files
        tumor_sample_name=$(samtools view -H ${tumorBam} | grep '@RG' | cut -d' ' -f10 | cut -d':' -f3)
        normal_sample_name=$(samtools view -H ${normalBam} | grep '@RG' | cut -d' ' -f10 | cut -d':' -f3)

        # Run mutectcaller with PON
        pbrun mutectcaller --ref ${params.reference} --tumor-name ${tumor_sample_name} --in-tumor-bam ${tumorBam} --in-normal-bam ${normalBam} --normal-name ${normal_sample_name} --pon ${ponFile} --out-vcf ${sample_id}.vcf
        """
}

// Workflow
workflow {
    fq2bamTumor(tumorFastqPairs)
    if (params.controlFastqDir) {
        fq2bamControl(controlFastqPairs)
        createPon(controlBamFiles)
        mutectcallerTumorNormalPon(tumorBamFiles.combine(controlBamFiles), ponFile)
    } else {
        mutectcallerTumorOnly(tumorBamFiles)
    }
}

// Publish results
tumorBamFiles.publishDir("${params.outdir}/tumor", mode: 'copy')
controlBamFiles.publishDir("${params.outdir}/control", mode: 'copy')
tumorVcfFiles.publishDir("${params.outdir}/tumor_vcfs", mode: 'copy')
tumorNormalPonVcfFiles.publishDir("${params.outdir}/tumor_normal_pon_vcfs", mode: 'copy')
