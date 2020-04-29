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
        lambda wildcards: '-a "A{100}" -q 20 -m 20'
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


rule star_index:
    input:
        fasta = srcdir(config['reference'])
    output:
        directory('reference/')
    log:
        'logs/star_index_reference.log'
    wrapper:
        '0.50.4/bio/star/index'


rule star_se:
    input:
        fq1 = 'trimmed/{accession}.trimmed.fastq',
        ref = 'reference/'
    output:
        'alignment/{accession}/Aligned.sortedByCoord.out.bam'
    log:
        'logs/star/{accession}.log'
    params:
        index = 'reference/',
        extra = '--outSAMtype BAM SortedByCoordinate'
    wrapper:
        '0.50.4/bio/star/align'


rule samtools_index:
    input:
        'alignment/{accession}/Aligned.sortedByCoord.out.bam'
    output:
        'alignment/{accession}/Aligned.sortedByCoord.out.bam.bai'
    wrapper:
        '0.50.4/bio/samtools/index'


rule aggregate_results:
    input:
        fastq_files = expand(
            'trimmed/{accession}.trimmed.fastq',
            accession=config['samples'].keys()),
        bam_files = expand(
            'alignment/{accession}/Aligned.sortedByCoord.out.bam',
            accession=config['samples'].keys()),
        index_files = expand(
            'alignment/{accession}/Aligned.sortedByCoord.out.bam.bai',
            accession=config['samples'].keys()),
        fname_ref = srcdir(config['reference'])
    output:
        fname = 'results/counts.csv',
        plot_dir = directory('plots')
    script:
        'scripts/aggregate.py'
