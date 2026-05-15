process ADJECTIVES {
    label 'process_small'

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/click_lark_networkx_numpy_pruned:4f796d4dd0e33580"

    input:
    path( ch_combined_annotations, stageAs: "raw-annotations.tsv" )
    path( rules_tsv )

    output:
    path("traits.xlsx"), emit: adjectives_ch

    script:
    def args = task.ext.args ?: ""

    """
    # export constants for script
    export FASTA_COLUMN="${params.CONSTANTS.FASTA_COLUMN}"

    adjectives.py --annotations ${ch_combined_annotations} --output traits.xlsx --rules_tsv '${rules_tsv}' ${args}
    """
}
