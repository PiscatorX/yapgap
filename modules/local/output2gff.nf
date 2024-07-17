process TOGFF{

    tag "$meta.id"
    label 'process_low'

    // conda "bioconda::biopython=1.70"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'biocontainers/biopython:1.70--np112py27_1' }"

    input:
    val meta
    path output
    val tool

    output:
    tuple val(meta), path("*.gff"), emit: gff
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    
    gff3generatorx.py \
        ${output} \
	--outfile \
        ${prefix}_${tool}.gff

   cat <<-END_VERSIONS > versions.yml
   "${task.process}":
	gff3generatox.py: v1
   END_VERSIONS
   """
}

