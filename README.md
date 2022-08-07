# nextflow-notes
Nextflow and the Half Blooded Prince

## Indexing

This is my informal learning notes of nextflow usage. For NGS analysis, we first prepare the reference therefore the script for that is:
We use `//` to comment a line.

```bash
// Declaration of DSL version. By default its dsl=2 anyways
nextflow.enable.dsl=2

// Input parameters
params.ref = "resources/sequence.fasta"
params.outdir="results/00_indexes"

// Definition of commands
process BWAINDEX{
	publishDir "$params.outdir"   //make output directory
	input:
	path (fasta)
	output:
	path "*"
	script:
	"""
	mkdir bwaidx
	cp $fasta bwaidx/
	bwa index bwaidx/$fasta
	samtools faidx bwaidx/$fasta
	"""
}

fa_ch=Channel.fromPath(params.ref, checkIfExists: true)
workflow {
  BWAINDEX(fa_ch)
}
```

`Channel.fromPath` is used to input files sequentially along with their path.

## Fastqc 

```bash
nextflow.enable.dsl=2
params.raw = "data/*{1,2}.fastq.gz"
params.outdir="results/01_rawfastqc"
process FASTQC {
	publishDir "$params.outdir"
	input:
	tuple val(sample_id), path(reads)
	output:
	file("*.{html,zip}")
	script:
	"""
	fastqc -t 10 ${reads[0]} ${reads[1]}
	"""
}


reads_ch = Channel.fromFilePairs(params.raw, checkIfExists: true )

workflow {
  FASTQC(reads_ch)
}


```
## Trimming


```bash
nextflow.enable.dsl=2
params.raw = "data/*{1,2}.fastq.gz"
params.adapter = "$baseDir/resources/universal_adapter.fasta"
process TRIM {
	publishDir "results/02_trimmed"
	input:
	tuple val(sample_id), file(reads)
	output:
          tuple val(sample_id), file(fq_1_paired), file(fq_2_paired)
	script:
    fq_1_paired = sample_id + '_R1_P.fastq.gz'
    fq_1_unpaired = sample_id + '_R1_UP.fastq.gz'
    fq_2_paired = sample_id + '_R2_P.fastq.gz'
    fq_2_unpaired = sample_id + '_R2_UP.fastq.gz'
	"""
	##mkdir -p trimmed/trim_paired trimmed/trim_unpaired
	trimmomatic PE -threads 20 ${reads[0]} ${reads[1]} \
	$fq_1_paired \
	$fq_1_unpaired \
	$fq_2_paired \
	$fq_2_unpaired \
	ILLUMINACLIP:${params.adapter}:2:30:10 \
	MINLEN:30
	"""
}

reads_ch = Channel.fromFilePairs(params.raw, checkIfExists: true )
workflow {
  TRIM(reads_ch)
}
```
