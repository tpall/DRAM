process SPLIT_POOLED_HITS {
    label 'process_single'

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_hmmer_mmseqs2_pruned:d2c88b719ab1322c"

    tag { db_name }

    input:
    tuple val( pool_name ), path( pooled_hits_csv )
    val( db_name )

    output:
    path( "*___mmseqs_${db_name}_formatted.csv" ), emit: per_genome_hits, optional: true

    script:
    """
    split_pooled_hits.py ${pooled_hits_csv} ${db_name}
    """
}
