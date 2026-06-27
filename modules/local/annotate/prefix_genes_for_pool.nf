process PREFIX_GENES_FOR_POOL {
    label 'process_single'

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_hmmer_mmseqs2_pruned:d2c88b719ab1322c"

    tag { input_fasta }

    input:
    tuple val( input_fasta ), path( proteins ), path( gene_locs )

    output:
    path( "${input_fasta}.pooled_genes.faa" ),     emit: faa
    path( "${input_fasta}.pooled_gene_locs.tsv" ), emit: locs

    script:
    // Gene ids from Prodigal are {scaffold}_{gene} and only unique WITHIN a genome
    // (megahit k141_* scaffold names collide across bins). Before pooling, make ids
    // globally unique by prefixing the genome name + the pipeline's "___" delimiter,
    // in both the FASTA headers and the gene-locs query_id column. SPLIT_POOLED_HITS
    // recovers the genome from this prefix and strips it to restore the original id.
    """
    sed 's/^>/>${input_fasta}___/' ${proteins} > ${input_fasta}.pooled_genes.faa
    awk 'BEGIN{FS=OFS="\\t"} NR==1{print; next} {\$1="${input_fasta}___"\$1; print}' ${gene_locs} > ${input_fasta}.pooled_gene_locs.tsv
    """
}
