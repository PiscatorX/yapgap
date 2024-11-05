process TRNASCANSE {

    tag "$meta.id"
    label 'process_low'

    conda "bioconda::trnascan-se=2.0.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trnascan-se:2.0.12--pl5321h031d066_0':
        'biocontainers/trnascan-se:2.0.12--pl5321h031d066_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*"), emit: trnascan
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${fasta.baseName}"

    """
    
    tRNAscan-SE \
        --threads $task.cpus \
        -o ${prefix}.out \
 	-f ${prefix}.stats \
        -m ${prefix}.log \
	-l ${prefix}.ss \
        -r ${prefix}.fpass_out \
        -F ${prefix}.fos \
	${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(tRNAscan-SE --version 2>&1) | grep tRNAscan-SE  /tmp/tr | cut -d " " -f 2- ))
    END_VERSIONS
     
    """

}




