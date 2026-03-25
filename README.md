# idp-germline-panel-nf

Pipeline Nextflow para análise germinativa de painel em `hg38`, com:

- QC de FASTQ e BAM
- alinhamento com `bwa mem`
- marcação de duplicatas com GATK
- BQSR opcional
- chamada de SNV/indel com `HaplotypeCaller` + `GenotypeGVCFs`
- CNV com ClinCNV
- agregação de QC com MultiQC

## Estado atual

O pipeline já suporta dois modos:

- `single_sample`: FASTQ ou BAM para análise germinativa completa
- `pon_build`: geração da matriz `.cov` do ClinCNV a partir de BAMs normais

A integração com ClinCNV foi ajustada para o fluxo real da ferramenta: ela **não** roda em BAM diretamente. O pipeline agora:

1. extrai cobertura média por alvo com `samtools bedcov`
2. gera uma matriz `.cov` para o cohort de normais
3. anexa a amostra analisada a essa matriz
4. executa `Rscript clinCNV.R --normal ... --bed ... --out ...`

## Estrutura mínima de entrada

### FASTQ

Use [`samplesheet.example.csv`](/home/bioinfo/idp-germline-panel-nf/samplesheet.example.csv) como modelo:

```csv
sample_id,fastq_r1,fastq_r2
SAMPLE_001,/data/fastq/SAMPLE_001_R1.fastq.gz,/data/fastq/SAMPLE_001_R2.fastq.gz
```

### BAM

Use [`samplesheet_bam.example.csv`](/home/bioinfo/idp-germline-panel-nf/samplesheet_bam.example.csv) como modelo:

```csv
sample_id,bam
SAMPLE_001,/data/bam/SAMPLE_001.bam
```

### Normais para ClinCNV

Use [`normal_bamsheet.example.csv`](/home/bioinfo/idp-germline-panel-nf/normal_bamsheet.example.csv) como modelo:

```csv
sample_id,bam
NORMAL_001,/data/bam/NORMAL_001.bam
NORMAL_002,/data/bam/NORMAL_002.bam
```

## Requisitos de referência

### Para SNV/indel

- `--fasta`: referência FASTA
- `--panel_bed`: BED dos alvos do painel
- `--intervals`: opcional; se omitido, usa `--panel_bed`
- `--bwa_index`: prefixo do índice do BWA quando `input_format=fastq`
- `--known_sites`: obrigatório apenas se `--run_bqsr true`

### Para ClinCNV

- `--clincnv_targets`: BED **anotado com GC no 4º campo**
- `--clincnv_script_dir`: diretório do checkout do ClinCNV contendo `clinCNV.R`
- `--clincnv_normal_cov`: matriz `.cov` com os normais, necessária em `single_sample`

Importante:

- O `ClinCNV` espera uma matriz `.cov` com cabeçalho e colunas `chr start end sample1 sample2 ...`.
- O BED do ClinCNV deve estar anotado com GC-content. O README oficial informa que o 4º campo deve conter GC de `0` a `1`.
- Para `hg38`, o pipeline passa `--hg38` automaticamente quando `clincnv_hg38=true`.

## Setup recomendado do ClinCNV

A forma mais direta para este pipeline é usar o código oficial do ClinCNV localmente.

### 1. Clonar o ClinCNV

```bash
git clone https://github.com/imgag/ClinCNV /opt/ClinCNV
```

### 2. Instalar dependências R

Segundo o README oficial do projeto, instale pelo menos:

```r
install.packages("optparse")
install.packages("robustbase")
install.packages("MASS")
install.packages("data.table")
install.packages("foreach")
install.packages("doParallel")
install.packages("mclust")
install.packages("R.utils")
install.packages("RColorBrewer")
install.packages("party")
install.packages("dbscan")
install.packages("umap")
install.packages("Rcpp")
```

### 3. Sintaxe ClinCNV usada pelo pipeline

O pipeline chama o ClinCNV neste formato:

```bash
Rscript /opt/ClinCNV/clinCNV.R \
  --normal results/clincnv/SAMPLE_001.clincnv.cov \
  --bed /refs/panel.gc_annotated.bed \
  --out results/clincnv/clincnv_SAMPLE_001 \
  --folderWithScript /opt/ClinCNV \
  --hg38 \
  --scoreG 50 \
  --lengthG 2 \
  --reanalyseCohort T
```

Esse formato segue o README oficial do ClinCNV para contexto germinativo, com `--normal`, `--bed`, `--out` e `--folderWithScript`.

### Alternativa com Singularity/Apptainer

O README oficial também informa um container SIF e a chamada:

```bash
singularity exec -B /bind/paths ClinCNV_v1.18.3.sif clinCNV.R parameters
```

Neste repositório, a integração automatizada está preparada para o modo com checkout local do ClinCNV via `Rscript`.

## Comandos padrão

### 1. Construir a matriz de normais do ClinCNV

```bash
nextflow run . -profile docker \
  --mode pon_build \
  --normal_bamsheet normal_bamsheet.example.csv \
  --fasta /refs/hg38.fa \
  --panel_bed /refs/panel_targets.bed \
  --clincnv_targets /refs/panel.gc_annotated.bed \
  --outdir results_pon
```

Saída principal esperada:

- `results_pon/clincnv/clincnv_panel_of_normals.cov`

### 2. Rodar uma amostra a partir de FASTQ

```bash
nextflow run . -profile docker \
  --mode single_sample \
  --input_format fastq \
  --sample_sheet samplesheet.example.csv \
  --fasta /refs/hg38.fa \
  --panel_bed /refs/panel_targets.bed \
  --intervals /refs/panel_targets.bed \
  --bwa_index /refs/bwa/hg38.fa \
  --known_sites /refs/dbsnp_hg38.vcf.gz \
  --run_bqsr true \
  --clincnv_targets /refs/panel.gc_annotated.bed \
  --clincnv_script_dir /opt/ClinCNV \
  --clincnv_normal_cov results_pon/clincnv/clincnv_panel_of_normals.cov \
  --outdir results_sample
```

### 3. Rodar uma amostra a partir de BAM

```bash
nextflow run . -profile docker \
  --mode single_sample \
  --input_format bam \
  --sample_sheet samplesheet_bam.example.csv \
  --fasta /refs/hg38.fa \
  --panel_bed /refs/panel_targets.bed \
  --clincnv_targets /refs/panel.gc_annotated.bed \
  --clincnv_script_dir /opt/ClinCNV \
  --clincnv_normal_cov results_pon/clincnv/clincnv_panel_of_normals.cov \
  --outdir results_sample
```

## Principais parâmetros

### Gerais

- `--mode`: `single_sample` ou `pon_build`
- `--input_format`: `fastq` ou `bam`
- `--outdir`: diretório final de resultados
- `--publish_dir_mode`: padrão `copy`

### SNV/indel

- `--run_fastp`: padrão `true`
- `--run_bqsr`: padrão `false`
- `--known_sites`: exigido se `run_bqsr=true`

### ClinCNV

- `--run_clincnv`: padrão `true`
- `--clincnv_hg38`: padrão `true`
- `--clincnv_min_mapq`: padrão `3`
- `--clincnv_score_g`: padrão `50`
- `--clincnv_length_g`: padrão `2`
- `--clincnv_reanalyse_cohort`: padrão `T`
- `--clincnv_extra_args`: argumentos extras repassados ao `clinCNV.R`

## Saídas principais

### `single_sample`

- `results/alignment/`: BAMs finais
- `results/variants/`: `g.vcf.gz` e `vcf.gz`
- `results/qc/`: métricas de alinhamento, PICARD/GATK e `bcftools stats`
- `results/multiqc/multiqc_report.html`
- `results/clincnv/`: arquivos do ClinCNV
- `results/run_summary.json`

### `pon_build`

- `results/clincnv/clincnv_panel_of_normals.cov`
- `results/run_summary.json`

## Limitações atuais

- O ClinCNV fica dependente de um checkout local funcional com R e pacotes instalados.
- O pipeline ainda não faz anotação automática do BED com GC-content; isso deve ser preparado antes.
- O fluxo de CNV está focado em germline por profundidade on-target; BAF/off-target ainda não foram integrados.

## Referências usadas para a integração do ClinCNV

- Repositório oficial: `https://github.com/imgag/ClinCNV`
- Checkout validado localmente em `/tmp/ClinCNV`
- Sintaxe oficial documentada no README do projeto:
  - `Rscript clinCNV.R --normal normal.cov --out outputFolder --bed annotatedBedFile.bed --folderWithScript $PWD`
  - `singularity exec -B /bind/paths ClinCNV_v1.18.3.sif clinCNV.R parameters`
