/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { NANOFILT               } from '../modules/nf-core/nanofilt/main'
include { FASTQC as FASTQC_2     } from '../modules/nf-core/fastqc/main'
include { MINIMAP2_ALIGN         } from '../modules/nf-core/minimap2/align/main'  
include { MOSDEPTH               } from '../modules/nf-core/mosdepth/main'
include { SAMTOOLS_FLAGSTAT      } from '../modules/nf-core/samtools/flagstat/main'
include { GUNZIP                 } from '../modules/nf-core/gunzip/main' 
include { SAMTOOLS_FAIDX         } from '../modules/nf-core/samtools/faidx/main'  
include { MEDAKA_CONSENSUS       } from '../modules/local/medaka/main'
include { QUAST                  } from '../modules/nf-core/quast/main'
include { BAKTA_BAKTADBDOWNLOAD  } from '../modules/nf-core/bakta/baktadbdownload/main' 
include { BAKTA_BAKTA            } from '../modules/nf-core/bakta/bakta/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_microbiome_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MICROBIOME {

    take:
    reads          // channel: samplesheet read in from --input
    reference          // channel: fasta reference in from --input
    
    main:
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()
    //
    // MODULE: Run FastQC
    //
    FASTQC (
        reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    //
    // MODULE: Run NanoFilt
    //
    NANOFILT (
        reads,
        []
    )
    ch_versions = ch_versions.mix(NANOFILT.out.versions.first())
    ch_trimmed_reads = NANOFILT.out.filtreads
    //
    // MODULE: Run FastQC after NanoFilt
    //
    FASTQC_2 (
        ch_trimmed_reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_2.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC_2.out.versions.first())
    //
    // MODULE: Run Minimap2 Align
    //
    reads_fasta = ch_trimmed_reads
        .join(reference)
        .multiMap {
            meta, rds, fa ->
            reads: [ meta, rds ]
            reference: [ meta, fa ]
        }
    ch_trimmed_reads = reads_fasta.reads
    ch_reference     = reads_fasta.reference
    MINIMAP2_ALIGN (
        ch_trimmed_reads,
        ch_reference,
        true,
        'bai',
        false,
        false
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())
    ch_bam = MINIMAP2_ALIGN.out.bam
    ch_bai = MINIMAP2_ALIGN.out.index
    //
    // MODULE: Run Mosdepth
    //
    ch_align    = ch_bam.join(ch_bai)
    ch_mosdepth = ch_align.join(ch_reference).multiMap {
        meta, bam, bai, fa ->
        align: [ meta, bam, bai, [] ]
        reference: [ meta, fa ]
    }
    ch_mosdepth_input   = ch_mosdepth.align
    ch_reference        = ch_mosdepth.reference
    MOSDEPTH (
        ch_mosdepth_input,
        ch_reference
    )
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.global_txt.collect{it[1]})
    //        
    // MODULE: Run Flagstat
    //
    SAMTOOLS_FLAGSTAT (
        ch_bam.join(ch_bai)
    )
    ch_versions      = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_FLAGSTAT.out.flagstat.collect{it[1]})
    //    
    // MODULE: Run Gunzip
    //
    GUNZIP (
        ch_reference
    )
    ch_reference = GUNZIP.out.gunzip
    //    
    // MODULE: Run Faidx
    //
    SAMTOOLS_FAIDX (
        ch_reference,
        [[],[]],
        false
    )
    ch_fai = SAMTOOLS_FAIDX.out.fai
    // 
    // MODULE: Run Medaka Consensus
    //
    ch_reference_fai        = ch_reference.join(ch_fai)
    ch_align_reference_fai  = ch_align
        .join(ch_reference_fai)
        .multiMap {
            meta, bam, bai, fa, fai ->
                align: [ meta, bam, bai ]
        reference_fai: [ meta, fa, fai ]
        }
    ch_align        = ch_align_reference_fai.align
    ch_reference_fai    = ch_align_reference_fai.reference_fai
    MEDAKA_CONSENSUS (
        ch_align,
        ch_reference_fai
    )
        ch_contigs = MEDAKA_CONSENSUS.out.fa
    // 
    // MODULE: Run Bacta DB Download
    //
    BAKTA_BAKTADBDOWNLOAD ()
    ch_versions = ch_versions.mix(BAKTA_BAKTADBDOWNLOAD.out.versions)
    ch_bacta_db = BAKTA_BAKTADBDOWNLOAD.out.db
    // 
    // MODULE: Run Bacta Annotation
    //
    BAKTA_BAKTA (
        ch_contigs,
        ch_bacta_db,
        [],
        [],
        [],
        [],
    )
    ch_gff              = BAKTA_BAKTA.out.gff
    ch_versions         = ch_versions.mix(BAKTA_BAKTA.out.versions.first())
    ch_multiqc_files    = ch_multiqc_files.mix(BAKTA_BAKTA.out.txt.map{it[1]})
    // 
    // MODULE: Run Quast
    //
    ch_reference_contigs_gff  = ch_reference
        .join(ch_contigs)
        .join(ch_gff)
        .multiMap {
            meta, ref, contigs, gff ->
                reference: [ meta, ref ]
                  contigs: [ meta, contigs ]
                      gff: [ meta, gff ]
        }
    ch_reference    = ch_reference_contigs_gff.reference
    ch_contigs      = ch_reference_contigs_gff.contigs
    ch_gff          = ch_reference_contigs_gff.gff
    QUAST (
        ch_reference,
        ch_contigs,
        ch_gff,
    )
    ch_versions = ch_versions.mix(BAKTA_BAKTA.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(BAKTA_BAKTA.out.txt.map{it[1]})
    // 
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'microbiome_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
