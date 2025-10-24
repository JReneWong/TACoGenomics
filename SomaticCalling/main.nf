nextflow.enable.dsl=2

// ---------- Params ----------
params.input_dir  = params.input_dir  ?: 'InputData'
params.outdir     = params.outdir     ?: 'results'
params.ref        = params.ref        ?: 'GenomeRef/GRCh38.d1.vd1.fa'
params.parabricks_image = params.parabricks_image ?: 'nvcr.io/nvidia/clara/clara-parabricks:4.4.0-1'
params.bwa_nstreams = (params.bwa_nstreams ?: 1) as int
params.low_memory   = (params.low_memory == null) ? true : (params.low_memory as boolean)
params.picard_image = params.picard_image ?: 'broadinstitute/picard:2.27.5'
params.dbsnp         = params.dbsnp         ?: 'GenomeRef/Homo_sapiens_assembly38.dbsnp138.vcf'
params.mills         = params.mills         ?: 'GenomeRef/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
params.known_indels  = params.known_indels  ?: 'GenomeRef/Homo_sapiens_assembly38.known_indels.vcf.gz'
params.outdir           = params.outdir           ?: 'results'
params.parabricks_image = params.parabricks_image ?: 'nvcr.io/nvidia/clara/clara-parabricks:4.4.0-1'
params.picard_image     = params.picard_image     ?: 'broadinstitute/picard:2.27.5'
params.gatk_image       = params.gatk_image       ?: 'broadinstitute/gatk:4.4.0.0'
params.dbsnp            = params.dbsnp            ?: 'GenomeRef/Homo_sapiens_assembly38.dbsnp138.vcf'
params.mills            = params.mills            ?: 'GenomeRef/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
params.known_indels     = params.known_indels     ?: 'GenomeRef/Homo_sapiens_assembly38.known_indels.vcf.gz'
params.bwa_nstreams     = params.bwa_nstreams     ?: 1
params.low_memory       = (params.low_memory in [true,false]) ? params.low_memory : true


// ---------- Process ----------
process FQ2BAM {
  // container comes from config; uncomment to hardcode:
  // container "${params.parabricks_image}"
  // containerOptions '--gpus all'

  publishDir "${params.outdir}/01_fq2bam", mode: 'copy'
  tag { "${sample}.${role}" }

  input:
    tuple val(sample), val(role), path(reads), path(ref_dir), val(ref_name)

  output:
    tuple val(sample), val(role), path("${sample}_${role}.raw.bam")

  script:
    def ref_fa = "${ref_dir}/${ref_name}"
    def r1 = reads[0]
    def r2 = reads[1]
    def lowMemFlag = params.low_memory ? '--low-memory' : ''
    """
    echo "== Debug: reference dir = ${ref_dir}"
    ls -l ${ref_dir} | head -n 60 || true
    pbrun --version || which pbrun
    nvidia-smi || true

    pbrun fq2bam \
      --ref ${ref_fa} \
      --in-fq ${r1} ${r2} \
      --out-bam ${sample}_${role}.raw.bam \
      ${lowMemFlag} \
      --bwa-nstreams ${params.bwa_nstreams}
    """
}

process ADD_RG {
  // Override the default Parabricks image *only for this step*
  container "${params.picard_image}"
  publishDir "${params.outdir}/02_readgroups", mode: 'copy'
  tag { "${sample}.${role}" }

  input:
    tuple val(sample), val(role), path(bam)

  output:
    tuple val(sample), val(role), path("${sample}_${role}.rg.bam")

  script:
    def rgid = (role == 'Tumor') ? "tm_${sample}" : "nm_${sample}"
    def rgpu = (role == 'Tumor') ? "tm_${sample}_unit1" : "nm_${sample}_unit1"
    """
    java -jar /usr/picard/picard.jar AddOrReplaceReadGroups \
      I=${bam} \
      O=${sample}_${role}.rg.bam \
      RGID=${rgid} \
      RGLB=lib1 \
      RGPL=ILLUMINA \
      RGPU=${rgpu} \
      RGSM=${role}
    """
}

process BQSR_TABLE {
  // Uses your default Parabricks container from docker_default.config
  publishDir "${params.outdir}/03_bqsr_tables", mode: 'copy'
  tag { "${sample}.${role}" }

  input:
    tuple val(sample), val(role), path(bam), path(ref_dir), val(ref_name), path(known_sites)

  output:
    // keep bam for the next step, and emit recal + ref info
    tuple val(sample), val(role), path(bam), path("${sample}_${role}.recal.txt"), path(ref_dir), val(ref_name)

  script:
    def ref_fa = "${ref_dir}/${ref_name}"
    // known_sites is a list-like staging of the three VCFs
    def ks = known_sites.collect { " --knownSites ${it}" }.join('')
    """
    pbrun bqsr \
      --ref ${ref_fa} \
      --in-bam ${bam} \
      --out-recal-file ${sample}_${role}.recal.txt \
      ${ks}
    """
}

process APPLY_BQSR {
  // Uses your default Parabricks container from docker_default.config
  publishDir "${params.outdir}/04_bqsr_bams", mode: 'copy'
  tag { "${sample}.${role}" }

  input:
    tuple val(sample), val(role), path(bam), path(recal), path(ref_dir), val(ref_name)

  output:
    // keep recal for downstream mutect; emit final bqsr BAM + recal
    tuple val(sample), val(role), path("${sample}_${role}.bqsr.bam"), path(recal)

  script:
    def ref_fa = "${ref_dir}/${ref_name}"
    """
    pbrun applybqsr \
      --ref ${ref_fa} \
      --in-bam ${bam} \
      --in-recal-file ${recal} \
      --out-bam ${sample}_${role}.bqsr.bam
    """
}

process MUTECTCALLER {
  publishDir "${params.outdir}/05_mutect", mode: 'copy'
  tag { sample }

  input:
    tuple val(sample),
          path(t_bam), path(t_recal),
          path(n_bam), path(n_recal),
          path(ref_dir), val(ref_name)

  output:
    // emit both VCF and its stats
    tuple val(sample), path("${sample}.mutect.vcf.gz"), path("${sample}.mutect.vcf.gz.stats")

  script:
    def ref_fa = "${ref_dir}/${ref_name}"
    """
    pbrun mutectcaller \
      --ref ${ref_fa} \
      --in-tumor-bam ${t_bam} \
      --in-normal-bam ${n_bam} \
      --in-tumor-recal-file ${t_recal} \
      --in-normal-recal-file ${n_recal} \
      --out-vcf ${sample}.mutect.vcf.gz \
      --tumor-name Tumor \
      --normal-name Normal
    """
}

process FILTER_MUTECT_CALLS {
  // Use GATK container for this step
  container "${params.gatk_image}"
  publishDir "${params.outdir}/06_filtered_vcf", mode: 'copy'
  tag { sample }

  input:
    tuple val(sample), path(vcf), path(stats), path(ref_dir), val(ref_name)

  output:
    path("${sample}.filtered.vcf")

  script:
    def ref_fa = "${ref_dir}/${ref_name}"
    """
    set -euo pipefail

    # 1) Index the bgzipped VCF so GATK can read it
    gatk IndexFeatureFile -I ${vcf}

    # 2) Filter Mutect calls (now the index is present)
    gatk FilterMutectCalls \
      -V ${vcf} \
      --stats ${stats} \
      -O ${sample}.filtered.vcf \
      -R ${ref_fa}
    """
}

// ---------- Workflow ----------
workflow {
  println "INPUT DIR: ${params.input_dir}"

  Channel.fromPath("${params.input_dir}/*.fastq.gz")
         .map { it.name }
         .view { "Found: ${it}" }

  // Pair Tumor/Normal R1/R2
  def paired = Channel
    .fromPath("${params.input_dir}/*.fastq.gz")
    .map { f ->
      def m = (f.name =~ /^(.+?)_(T|N)_.+\.R([12])\.fastq\.gz$/)
      assert m.matches() : "Filename doesn't match expected pattern: ${f.name}"
      tuple(m[0][1], (m[0][2]=='T' ? 'Tumor':'Normal'), (m[0][3] as int), f)
    }
    .groupTuple(by:[0,1])                         // -> sample, role, [mates], [files]
    .map { sample, role, mates, files ->
      def r = [1:null, 2:null]
      for (int i=0; i<mates.size(); i++) r[mates[i] as int] = files[i]
      assert r[1] && r[2] : "Missing R1/R2 for ${sample}/${role}"
      tuple(sample, role, [ r[1], r[2] ])         // -> (sample, role, [R1,R2])
    }

  // Compute ref dir + filename once, attach to each tuple
  def REF      = file(params.ref)
  def REF_DIR  = REF.parent
  def REF_NAME = REF.name

  def inputs = paired.map { s, role, reads -> tuple(s, role, reads, REF_DIR, REF_NAME) }

  // 1) fq2bam -> raw BAMs
  raw_bams = FQ2BAM(inputs)     // emits: (sample, role, <sample_role>.raw.bam)

  // 2) Add read groups with Picard -> RG-fixed BAMs
  rg_bams  = ADD_RG(raw_bams)   // emits: (sample, role, <sample_role>.rg.bam)

  // Build a staged list of known-sites files once
  def KNOWN_SITES = [ file(params.dbsnp), file(params.mills), file(params.known_indels) ]

  // Attach ref + known sites to each RG-fixed BAM
  def bqsr_inputs = rg_bams.map { s, role, bam -> tuple(s, role, bam, REF_DIR, REF_NAME, KNOWN_SITES) }

  // 3) Create BQSR recal tables (keeps bam in the tuple)
  def bqsr_out = BQSR_TABLE(bqsr_inputs)     // -> (sample, role, bam, recal.txt, ref_dir, ref_name)

  // 4) Apply BQSR
  def bqsr_bams = APPLY_BQSR(bqsr_out)       // -> (sample, role, sample_role.bqsr.bam, recal.txt)

  // Split Tumor & Normal and DO NOT carry ref info here
  def tumors  = bqsr_bams
    .filter { it[1] == 'Tumor'  }                 // (sample, role, bqsr.bam, recal)
    .map    { s, r, bam, recal -> tuple(s, bam, recal) }

  def normals = bqsr_bams
    .filter { it[1] == 'Normal' }
    .map    { s, r, bam, recal -> tuple(s, bam, recal) }

  // Join by sample -> (sample, t_bam, t_recal, n_bam, n_recal)
  def tn_joined = tumors.join(normals)

  // Attach the single ref_dir/ref_name once
  def mutect_inputs = tn_joined.map { s, t_bam, t_recal, n_bam, n_recal ->
    tuple(s, t_bam, t_recal, n_bam, n_recal, REF_DIR, REF_NAME)
  }

  // 5) Mutect caller
  def mutect_vcfs = MUTECTCALLER(mutect_inputs)   // -> (sample, .vcf.gz, .vcf.gz.stats)

  // 6) Filter Mutect calls (add ref info + pass stats)
  def filter_inputs = mutect_vcfs.map { s, vcf, stats ->
    tuple(s, vcf, stats, REF_DIR, REF_NAME)
  }

  FILTER_MUTECT_CALLS(filter_inputs)
}

