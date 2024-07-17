process INFERNAL_SEARCH {

    tag "$meta2.id"
    label 'process_high'
    
    conda (params.enable_conda ? "bioconda::infernal=1.1.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/infernal:1.1.4--h779adbc_0':
        'quay.io/biocontainers/infernal:1.1.4--h779adbc_0' }"

    input:
    val meta
    val db_size_mb
    each path(fasta)
    tuple val(meta2), path(cm_file)
    path cm_index


    output:
    tuple val(meta), path("*.tbl"), emit: tbl
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    rfam_tbl = fasta.getBaseName() + ".rfam.tbl"
    rfam_txt = fasta.getBaseName() + ".rfam.out"
    """
    
    cmsearch \
       -Z ${db_size_mb} \
       --rfam \
       --cpu ${task.cpus} \
       --cut_tc \
       --tblout \
      $rfam_tbl \
      -o $rfam_txt \
      $cm_file \
      $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cmsearch: \$(echo \$(cmsearch -h) | head -n2 | tail -n1  | cut -f3 -d " " ))
    END_VERSIONS

    """
}
