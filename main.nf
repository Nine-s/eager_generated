nextflow.enable.dsl = 2


include { FASTQC as FASTQC_PRE_PREPROCESSING; FASTQC as FASTQC_POST_PREPROCESSING } from '../../../ninon/description_prototype/modules_eager/fastqc.nf'
include { FASTP  } from '../../../ninon/description_prototype/modules_eager/fastp_aDNA.nf'
include { ADAPTER_REMOVAL  } from '../../../ninon/description_prototype/modules_eager/adapter_removal.nf'
include { BOWTIE2_INDEX  } from '../../../ninon/description_prototype/modules_eager/bowtie2_index.nf'
include { BOWTIE2  } from '../../../ninon/description_prototype/modules_eager/bowtie2.nf'
include { SAMTOOLS_FILTER  } from '../../../ninon/description_prototype/modules_eager/samtools_filter.nf'
include { BBDUK  } from '../../../ninon/description_prototype/modules_eager/bbduk.nf'
include { KRAKEN ; KRAKEN_PARSE ; KRAKEN_MERGE  } from '../../../ninon/description_prototype/modules_eager/kraken.nf'

workflow{
        read_pairs_ch = Channel
            .fromPath( params.csv_input )
            .splitCsv(header: true, sep: ',')
            .map {row -> tuple(row.sample, [row.path_r1, row.path_r2], row.condition)}
            .view()
        
FASTQC_PRE_PREPROCESSING(read_pairs_ch)
FASTP(read_pairs_ch)
//ADAPTER_REMOVAL(FASTP.out.reads)
BOWTIE2_INDEX(params.genome)
FASTQC_POST_PREPROCESSING(FASTP.out.reads)
BOWTIE2(FASTP.out.reads, params.genome, BOWTIE2_INDEX.out.index)
SAMTOOLS_FILTER(BOWTIE2.out.bam)
BBDUK(SAMTOOLS_FILTER.out.unmapped_fastq)
KRAKEN(BBDUK.out.fastq)
KRAKEN_PARSE(KRAKEN.out.kraken_report)
KRAKEN_MERGE(KRAKEN_PARSE.out.kraken_parsed)

}