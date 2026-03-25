nextflow.enable.dsl = 2

import groovy.json.JsonOutput


def normalizeBool(value, boolean defaultValue = false) {
    if (value == null) {
        return defaultValue
    }
    if (value instanceof Boolean) {
        return value
    }
    value.toString().toLowerCase() in ['true', '1', 'yes', 'y']
}

def requireParam(String name) {
    if (!params[name]) {
        error "Parâmetro obrigatório ausente: --${name}"
    }
}

def requireFilePath(String label, def pathValue) {
    if (!pathValue) {
        error "Arquivo obrigatório ausente para ${label}"
    }
    file(pathValue, checkIfExists: true)
}

def resolveDelimitedSheet(String sheetPath) {
    def p = file(sheetPath, checkIfExists: true)
    def first = p.text.readLines().find { it?.trim() && !it.startsWith('#') }
    if (!first) {
        error "Planilha vazia: ${sheetPath}"
    }
    first.contains('\t') ? '\t' : ','
}

def loadRows(String sheetPath) {
    def sep = resolveDelimitedSheet(sheetPath)
    Channel
        .fromPath(sheetPath)
        .splitCsv(header: true, sep: sep)
        .filter { row -> row.any { k, v -> v?.toString()?.trim() } }
}

def toFastqTuple(Map row) {
    def sampleId = row.sample_id ?: row.sample ?: row.id
    def r1 = row.fastq_r1 ?: row.r1 ?: row.read1
    def r2 = row.fastq_r2 ?: row.r2 ?: row.read2
    if (!sampleId || !r1 || !r2) {
        error 'A sample_sheet FASTQ precisa das colunas sample_id, fastq_r1 e fastq_r2'
    }
    tuple(sampleId.toString(), file(r1.toString(), checkIfExists: true), file(r2.toString(), checkIfExists: true))
}

def toBamTuple(Map row) {
    def sampleId = row.sample_id ?: row.sample ?: row.id
    def bam = row.bam ?: row.cram
    if (!sampleId || !bam) {
        error 'A sample_sheet BAM precisa das colunas sample_id e bam/cram'
    }
    tuple(sampleId.toString(), file(bam.toString(), checkIfExists: true))
}

def toNormalBamTuple(Map row) {
    def sampleId = row.sample_id ?: row.sample ?: row.id
    def bam = row.bam ?: row.cram
    if (!sampleId || !bam) {
        error 'A normal_bamsheet precisa das colunas sample_id e bam/cram'
    }
    tuple(sampleId.toString(), file(bam.toString(), checkIfExists: true))
}

process FASTP {
    tag "$sample_id"
    label 'small'
    publishDir "${params.outdir}/qc/fastp", mode: params.publish_dir_mode

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}.trim.R1.fastq.gz"), path("${sample_id}.trim.R2.fastq.gz"), emit: reads
    path("${sample_id}.fastp.json"), emit: json
    path("${sample_id}.fastp.html"), emit: html

    script:
    """
    fastp \
        --in1 ${read1} \
        --in2 ${read2} \
        --out1 ${sample_id}.trim.R1.fastq.gz \
        --out2 ${sample_id}.trim.R2.fastq.gz \
        --thread ${task.cpus} \
        --json ${sample_id}.fastp.json \
        --html ${sample_id}.fastp.html
    """
}

process ALIGN_AND_SORT {
    tag "$sample_id"
    label 'medium'
    publishDir "${params.outdir}/alignment", mode: params.publish_dir_mode, pattern: '*.bam'
    publishDir "${params.outdir}/alignment", mode: params.publish_dir_mode, pattern: '*.bai'

    input:
    tuple val(sample_id), path(read1), path(read2)
    path bwa_index

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai")

    script:
    def rg = "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA\\tLB:${sample_id}\\tPU:${sample_id}"
    """
    bwa mem -M -t ${task.cpus} -R '${rg}' ${bwa_index} ${read1} ${read2} | \
        samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam -
    samtools index -@ ${task.cpus} ${sample_id}.sorted.bam
    """
}

process MARK_DUPLICATES {
    tag "$sample_id"
    label 'medium'
    publishDir "${params.outdir}/alignment", mode: params.publish_dir_mode, pattern: '*.bam'
    publishDir "${params.outdir}/alignment", mode: params.publish_dir_mode, pattern: '*.bai'
    publishDir "${params.outdir}/qc/picard", mode: params.publish_dir_mode, pattern: '*.dup_metrics.txt'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.md.bam"), path("${sample_id}.md.bam.bai"), emit: bam
    path("${sample_id}.dup_metrics.txt"), emit: metrics

    script:
    """
    gatk MarkDuplicates \
        -I ${bam} \
        -O ${sample_id}.md.bam \
        -M ${sample_id}.dup_metrics.txt \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY LENIENT
    """
}

process ENSURE_BAM_INDEX {
    tag "$sample_id"
    label 'small'
    publishDir "${params.outdir}/alignment", mode: params.publish_dir_mode, pattern: '*.bai'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path(bam), path("${bam.simpleName}.bai")

    script:
    """
    samtools index -@ ${task.cpus} ${bam}
    """
}

process BASE_RECALIBRATOR {
    tag "$sample_id"
    label 'medium'
    publishDir "${params.outdir}/gatk/bqsr", mode: params.publish_dir_mode

    input:
    tuple val(sample_id), path(bam), path(bai)
    path fasta
    path known_sites
    path intervals

    output:
    tuple val(sample_id), path(bam), path(bai), path("${sample_id}.recal.table")

    script:
    """
    gatk BaseRecalibrator \
        -R ${fasta} \
        -I ${bam} \
        --known-sites ${known_sites} \
        -L ${intervals} \
        -O ${sample_id}.recal.table
    """
}

process APPLY_BQSR {
    tag "$sample_id"
    label 'medium'
    publishDir "${params.outdir}/alignment", mode: params.publish_dir_mode, pattern: '*.bam'
    publishDir "${params.outdir}/alignment", mode: params.publish_dir_mode, pattern: '*.bai'

    input:
    tuple val(sample_id), path(bam), path(bai), path(recal_table)
    path fasta

    output:
    tuple val(sample_id), path("${sample_id}.bqsr.bam"), path("${sample_id}.bqsr.bam.bai")

    script:
    """
    gatk ApplyBQSR \
        -R ${fasta} \
        -I ${bam} \
        --bqsr-recal-file ${recal_table} \
        -O ${sample_id}.bqsr.bam
    samtools index -@ ${task.cpus} ${sample_id}.bqsr.bam
    """
}

process COLLECT_ALIGNMENT_QC {
    tag "$sample_id"
    label 'medium'
    publishDir "${params.outdir}/qc/alignment", mode: params.publish_dir_mode

    input:
    tuple val(sample_id), path(bam), path(bai)
    path fasta
    path intervals

    output:
    path("${sample_id}.flagstat.txt"), emit: flagstat
    path("${sample_id}.stats.txt"), emit: stats
    path("${sample_id}.idxstats.txt"), emit: idxstats
    path("${sample_id}.insert_size_metrics.txt"), emit: insert_metrics
    path("${sample_id}.insert_size_histogram.pdf"), emit: insert_histogram
    path("${sample_id}.hs_metrics.txt"), emit: hs_metrics

    script:
    """
    samtools flagstat ${bam} > ${sample_id}.flagstat.txt
    samtools stats ${bam} > ${sample_id}.stats.txt
    samtools idxstats ${bam} > ${sample_id}.idxstats.txt
    gatk CollectInsertSizeMetrics \
        -I ${bam} \
        -O ${sample_id}.insert_size_metrics.txt \
        -H ${sample_id}.insert_size_histogram.pdf
    gatk CollectHsMetrics \
        -I ${bam} \
        -O ${sample_id}.hs_metrics.txt \
        -R ${fasta} \
        -BI ${intervals} \
        -TI ${intervals}
    """
}

process HAPLOTYPE_CALLER {
    tag "$sample_id"
    label 'large'
    publishDir "${params.outdir}/variants", mode: params.publish_dir_mode, pattern: '*.vcf.gz'
    publishDir "${params.outdir}/variants", mode: params.publish_dir_mode, pattern: '*.vcf.gz.tbi'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path fasta
    path intervals

    output:
    tuple val(sample_id), path("${sample_id}.g.vcf.gz"), path("${sample_id}.g.vcf.gz.tbi")

    script:
    """
    gatk HaplotypeCaller \
        -R ${fasta} \
        -I ${bam} \
        -L ${intervals} \
        -ERC GVCF \
        -O ${sample_id}.g.vcf.gz
    """
}

process GENOTYPE_GVCF {
    tag "$sample_id"
    label 'medium'
    publishDir "${params.outdir}/variants", mode: params.publish_dir_mode, pattern: '*.vcf.gz'
    publishDir "${params.outdir}/variants", mode: params.publish_dir_mode, pattern: '*.vcf.gz.tbi'

    input:
    tuple val(sample_id), path(gvcf), path(tbi)
    path fasta
    path intervals

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi")

    script:
    """
    gatk GenotypeGVCFs \
        -R ${fasta} \
        -V ${gvcf} \
        -L ${intervals} \
        -O ${sample_id}.vcf.gz
    """
}

process VARIANT_QC {
    tag "$sample_id"
    label 'small'
    publishDir "${params.outdir}/qc/variants", mode: params.publish_dir_mode

    input:
    tuple val(sample_id), path(vcf), path(tbi)

    output:
    path("${sample_id}.bcftools_stats.txt")

    script:
    """
    bcftools stats ${vcf} > ${sample_id}.bcftools_stats.txt
    """
}

process CLINCNV_SAMPLE_COVERAGE {
    tag "$sample_id"
    label 'small'
    publishDir "${params.outdir}/clincnv/coverage", mode: params.publish_dir_mode

    input:
    tuple val(sample_id), path(bam), path(bai)
    path targets

    output:
    tuple val(sample_id), path("${sample_id}.depth.tsv")

    script:
    """
    awk 'BEGIN{OFS="\\t"} !/^#/ {print \$1,\$2,\$3}' ${targets} > targets.3col.bed
    samtools bedcov -Q ${params.clincnv_min_mapq} targets.3col.bed ${bam} | \
        awk 'BEGIN{OFS="\\t"} { len=\$3-\$2; cov=(len>0 ? \$4/len : 0); print cov }' > ${sample_id}.depth.tsv
    """
}

process CLINCNV_BUILD_COHORT_MATRIX {
    tag 'clincnv_matrix'
    label 'small'
    publishDir "${params.outdir}/clincnv", mode: params.publish_dir_mode

    input:
    path depth_files
    path targets

    output:
    path('clincnv_panel_of_normals.cov')

    script:
    """
    awk 'BEGIN{OFS="\\t"} !/^#/ {print \$1,\$2,\$3}' ${targets} > regions.tsv
    files=(\$(find . -maxdepth 1 -name '*.depth.tsv' | sort))
    if [[ \${#files[@]} -eq 0 ]]; then
        echo 'Nenhum arquivo de cobertura ClinCNV foi recebido' >&2
        exit 1
    fi
    header=\$'chr\tstart\tend'
    for f in "\${files[@]}"; do
        sample=\$(basename "\${f}" .depth.tsv)
        header+=\$'\t'"\${sample}"
    done
    printf '%s\n' "\${header}" > clincnv_panel_of_normals.cov
    paste regions.tsv "\${files[@]}" >> clincnv_panel_of_normals.cov
    """
}

process CLINCNV_PREPARE_SINGLE_MATRIX {
    tag "$sample_id"
    label 'small'
    publishDir "${params.outdir}/clincnv", mode: params.publish_dir_mode

    input:
    tuple val(sample_id), path(depth_file)
    path normal_cov

    output:
    tuple val(sample_id), path("${sample_id}.clincnv.cov")

    script:
    """
    matrix_lines=\$(( \$(wc -l < ${normal_cov}) - 1 ))
    sample_lines=\$(wc -l < ${depth_file})
    if [[ \${matrix_lines} -ne \${sample_lines} ]]; then
        echo "ClinCNV matrix/target mismatch for ${sample_id}: normal_cov=\${matrix_lines} rows, sample=\${sample_lines} rows" >&2
        exit 1
    fi
    {
        head -n 1 ${normal_cov} | awk -v sample='${sample_id}' 'BEGIN{OFS="\\t"} {print \$0, sample}'
        paste <(tail -n +2 ${normal_cov}) ${depth_file}
    } > ${sample_id}.clincnv.cov
    """
}

process CLINCNV_RUN_GERMLINE {
    tag "$sample_id"
    label 'large'
    publishDir params.clincnv_outdir, mode: params.publish_dir_mode

    input:
    tuple val(sample_id), path(clincnv_cov)
    path targets
    path clincnv_script

    output:
    path("${sample_id}.clincnv.files.txt")

    script:
    def hg38Arg = normalizeBool(params.clincnv_hg38, true) ? '--hg38' : ''
    def scoreArg = params.clincnv_score_g ? "--scoreG ${params.clincnv_score_g}" : ''
    def lengthArg = params.clincnv_length_g != null ? "--lengthG ${params.clincnv_length_g}" : ''
    def reanalyseArg = params.clincnv_reanalyse_cohort ? "--reanalyseCohort ${params.clincnv_reanalyse_cohort}" : ''
    def extraArgs = params.clincnv_extra_args ?: ''
    """
    mkdir -p clincnv_${sample_id}
    Rscript ${clincnv_script} \
        --normal ${clincnv_cov} \
        --bed ${targets} \
        --out clincnv_${sample_id} \
        --folderWithScript \
        ${clincnv_script.parent} \
        ${hg38Arg} \
        ${scoreArg} \
        ${lengthArg} \
        ${reanalyseArg} \
        ${extraArgs}
    find clincnv_${sample_id} -maxdepth 3 -type f | sort > ${sample_id}.clincnv.files.txt
    """
}

process MULTIQC {
    tag 'multiqc'
    label 'small'
    publishDir "${params.outdir}/multiqc", mode: params.publish_dir_mode

    input:
    path qc_files

    output:
    path('multiqc_report.html')
    path('multiqc_data')

    script:
    """
    multiqc . --force --outdir .
    """
}

process WRITE_RUN_SUMMARY {
    tag 'run_summary'
    label 'small'
    publishDir params.outdir, mode: params.publish_dir_mode

    input:
    val summary_map

    output:
    path('run_summary.json')

    script:
    def json = JsonOutput.prettyPrint(JsonOutput.toJson(summary_map))
    """
    cat > run_summary.json <<'JSON'
${json}
JSON
    """
}

workflow {
    log.info "Iniciando ${workflow.manifest.name ?: 'pipeline'} | modo=${params.mode} | genome=${params.genome}"

    requireParam('fasta')
    requireParam('panel_bed')

    def fasta = requireFilePath('fasta', params.fasta)
    def intervals = requireFilePath('panel_bed/intervals', params.intervals ?: params.panel_bed)
    def clincnvTargets = requireFilePath('clincnv_targets', params.clincnv_targets ?: params.panel_bed)

    def summary = [
        run_name           : workflow.runName,
        mode               : params.mode,
        input_format       : params.input_format,
        outdir             : params.outdir,
        genome             : params.genome,
        run_fastp          : normalizeBool(params.run_fastp, true),
        run_bqsr           : normalizeBool(params.run_bqsr, false),
        run_clincnv        : normalizeBool(params.run_clincnv, true),
        clincnv_script_dir : params.clincnv_script_dir,
        clincnv_normal_cov : params.clincnv_normal_cov
    ]

    if (params.mode == 'single_sample') {
        requireParam('sample_sheet')

        def finalBams
        def qcFiles = Channel.empty()

        if (params.input_format == 'fastq') {
            requireParam('bwa_index')

            def inputReads = loadRows(params.sample_sheet).map { row -> toFastqTuple(row) }
            def readsForAlignment
            def fastpJson = Channel.empty()
            def fastpHtml = Channel.empty()

            if (normalizeBool(params.run_fastp, true)) {
                def fastpResult = FASTP(inputReads)
                readsForAlignment = fastpResult.reads
                fastpJson = fastpResult.json
                fastpHtml = fastpResult.html
            } else {
                readsForAlignment = inputReads
            }

            def alignedBams = ALIGN_AND_SORT(readsForAlignment, requireFilePath('bwa_index', params.bwa_index))
            def dedup = MARK_DUPLICATES(alignedBams)
            finalBams = dedup.bam

            if (normalizeBool(params.run_bqsr, false)) {
                requireParam('known_sites')
                def recal = BASE_RECALIBRATOR(finalBams, fasta, requireFilePath('known_sites', params.known_sites), intervals)
                finalBams = APPLY_BQSR(recal, fasta)
            }

            def alignQc = COLLECT_ALIGNMENT_QC(finalBams, fasta, intervals)
            def gvcf = HAPLOTYPE_CALLER(finalBams, fasta, intervals)
            def vcf = GENOTYPE_GVCF(gvcf, fasta, intervals)
            def variantQc = VARIANT_QC(vcf)

            qcFiles = Channel.empty()
                .mix(fastpJson)
                .mix(fastpHtml)
                .mix(dedup.metrics)
                .mix(alignQc.flagstat)
                .mix(alignQc.stats)
                .mix(alignQc.idxstats)
                .mix(alignQc.insert_metrics)
                .mix(alignQc.insert_histogram)
                .mix(alignQc.hs_metrics)
                .mix(variantQc)
        } else if (params.input_format == 'bam') {
            def bamInput = loadRows(params.sample_sheet).map { row -> toBamTuple(row) }
            finalBams = ENSURE_BAM_INDEX(bamInput)
            def alignQc = COLLECT_ALIGNMENT_QC(finalBams, fasta, intervals)
            def gvcf = HAPLOTYPE_CALLER(finalBams, fasta, intervals)
            def vcf = GENOTYPE_GVCF(gvcf, fasta, intervals)
            def variantQc = VARIANT_QC(vcf)

            qcFiles = Channel.empty()
                .mix(alignQc.flagstat)
                .mix(alignQc.stats)
                .mix(alignQc.idxstats)
                .mix(alignQc.insert_metrics)
                .mix(alignQc.insert_histogram)
                .mix(alignQc.hs_metrics)
                .mix(variantQc)
        } else {
            error "input_format inválido: ${params.input_format}. Use 'fastq' ou 'bam'."
        }

        MULTIQC(qcFiles.collect())

        if (normalizeBool(params.run_clincnv, true)) {
            requireParam('clincnv_script_dir')
            requireParam('clincnv_normal_cov')
            def clincnvScript = requireFilePath('clincnv_script_dir/clinCNV.R', "${params.clincnv_script_dir}/clinCNV.R")
            def normalCov = requireFilePath('clincnv_normal_cov', params.clincnv_normal_cov)
            def sampleDepth = CLINCNV_SAMPLE_COVERAGE(finalBams, clincnvTargets)
            def singleMatrix = CLINCNV_PREPARE_SINGLE_MATRIX(sampleDepth, normalCov)
            CLINCNV_RUN_GERMLINE(singleMatrix, clincnvTargets, clincnvScript)
        }
    } else if (params.mode == 'pon_build') {
        requireParam('normal_bamsheet')
        def normalBams = loadRows(params.normal_bamsheet).map { row -> toNormalBamTuple(row) }
        def indexedNormals = ENSURE_BAM_INDEX(normalBams)
        def depthFiles = CLINCNV_SAMPLE_COVERAGE(indexedNormals, clincnvTargets)
            .map { sample_id, depth_file -> depth_file }
            .collect()
        CLINCNV_BUILD_COHORT_MATRIX(depthFiles, clincnvTargets)
    } else {
        error "Modo inválido: ${params.mode}. Use 'single_sample' ou 'pon_build'."
    }

    WRITE_RUN_SUMMARY(summary)
}
