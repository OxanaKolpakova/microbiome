process MEDAKA_CONSENSUS {
    tag "${meta.id}"
    conda 'bioconda::medaka'
    container 'biocontainers/medaka:2.0.0--py311hfd2b166_0'
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(bam), path(bai) 
    tuple val(meta2), path(reference), path(faindex)  

    output:
    tuple val(meta), path("*.fa"), emit: fa
    tuple val("${task.process}"), val('medaka'), eval("medaka --version 2>&1 | sed '1!d;s/.* //'"), topic: versions, emit: versions_medaka
    
    script:
    """
    medaka_consensus -i $bam -d $reference -o ${meta.id}_out -t ${task.cpus}
    mv ${meta.id}_out/consensus.fasta ${meta.id}.fa 
    """
}