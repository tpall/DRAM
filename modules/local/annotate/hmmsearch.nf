process HMM_SEARCH {
    label 'process_small'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_hmmer_mmseqs2_pruned:d2c88b719ab1322c"

    tag { input_fasta }

    input:
    tuple val( input_fasta ), path( fasta ), path( prodigal_locs_tsv )
    val ( e_value )
    path( database_loc )
    path( hmm_info_path )
    val (ec_from_info )
    val (db_name)

    output:
    tuple val( input_fasta ), path ( "${input_fasta}___formatted_${db_name}_hits.csv" ), emit: formatted_hits, optional: true

    script:
    def args = task.ext.args ?: ""
    def ec_flag = ec_from_info ? "--ec_from_info" : ""

    """
    # Bit-score prefilter (-T/--domT) instead of the e-value prefilter (-E): bit scores
    # are independent of search-space size (Z), so a pooled multi-genome search reports
    # the SAME hits as per-genome (hmm_parser then applies the per-profile bit-score
    # thresholds). Floor 10 sits below the lowest kofam threshold (15.07) so nothing a
    # real threshold would keep is pre-cut. NB: this is the perf/hmm-pool-bitscore
    # experiment — it changes results vs the -E ${e_value} prefilter.
    hmmsearch \\
    -T 10 --domT 10 \\
    --domtblout ${input_fasta}_hmmsearch.out \\
    --cpu ${task.cpus} \\
    ${database_loc}/*.hmm \\
    ${fasta} > /dev/null

    hmm_parser.py \\
        --hmm_domtbl ${input_fasta}_hmmsearch.out \\
        --hmm_info_path ${hmm_info_path} \\
        ${ec_flag} \\
        --gene_locs ${prodigal_locs_tsv} \\
        --db_name ${db_name} \\
        --output "${input_fasta}___formatted_${db_name}_hits.csv"
    """
}
