process SPALN_ALIGN {

    tag "$meta.id | $meta_p.id"
    label 'process_low'
    
    conda (params.enable_conda ? "bioconda::spaln=2.4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/spaln:2.4.7--pl5321h9a82719_1':
        'quay.io/biocontainers/spaln:2.4.7--pl5321h9a82719_1' }"

    input:
    tuple val(meta), path(spaln_index)
    tuple val(meta_p), path(proteins)
    val spaln_q
    val spaln_options

    output:
    tuple val(meta_p),path("${chunk_name}.*rd"), emit: align
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    chunk_name = proteins.getBaseName()

    """
    spaln \
        -o $chunk_name \
        -O 12 \
        -Q ${spaln_q} \
        -t ${task.cpus} \
        -d genome_spaln \
	 $proteins 
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spaln: \$(echo \$( spaln 2>&1 | head -n3 | tail -n1 | cut -f4 -d " " ))
    END_VERSIONS
   """
}



// spaln \
//         -o $chunk_name ${spaln_options} \
//         -Q ${spaln_q} \
//         -O 12 -t${task.cpus} \
//         -d genome_spaln $proteins
