configfile: 'config.yaml'
workdir: config['workdir']


rule all:
    input:
        'results/counts.csv',
        expand('qc/{accession}.html', accession=config['samples'].keys())


rule get_fastq_se:
    output:
        'data/{accession}.fastq'
    wrapper:
        '0.50.4/bio/sra-tools/fasterq-dump'


rule cutadapt:
    input:
        'data/{accession}.fastq'
    output:
        fastq = 'trimmed/{accession}.trimmed.fastq',
        qc = 'trimmed/{accession}.qc.txt'
    params:
        lambda wildcards: '-a "A{100}" -q 20'
    log:
        'logs/cutadapt/{accession}.log'
    wrapper:
        '0.50.4/bio/cutadapt/se'


rule fastqc:
    input:
        'trimmed/{accession}.trimmed.fastq'
    output:
        html = 'qc/{accession}.html',
        zip = 'qc/{accession}_fastqc.zip'
    log:
        'logs/fastqc/{accession}.log'
    wrapper:
        '0.50.4/bio/fastqc'


rule kallisto_index:
    input:
        fasta = srcdir(config['reference'])
    output:
        index = 'reference/reference.idx'
    log:
        'logs/kallisto_index_reference.log'
    wrapper:
        '0.50.4/bio/kallisto/index'


rule kallisto_quant:
    input:
        fastq = 'trimmed/{accession}.trimmed.fastq',
        index = 'reference/reference.idx'
    output:
        directory('quantification/{accession}')
    params:
        extra = '--single -l 200 -s 20'
    log:
        'logs/kallisto_quant_{accession}.log'
    wrapper:
        '0.50.4/bio/kallisto/quant'


rule aggregate_results:
    input:
        quant_dirs = expand(
            'quantification/{accession}',
            accession=config['samples'].keys())
    output:
        fname = 'results/counts.csv'
    script:
        'scripts/aggregate.py'
