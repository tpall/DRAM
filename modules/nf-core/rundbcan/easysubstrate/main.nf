process RUNDBCAN_EASYSUBSTRATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dbcan:5.2.6--pyhdfd78af_0' :
        'biocontainers/dbcan:5.2.6--pyhdfd78af_0' }"

    input:
    tuple val(meta),  path(input_raw_data)
    tuple val(meta2), path(input_gff), val(gff_type)
    path  dbcan_db

    output:
    tuple val(meta), path("${prefix}_overview.tsv")            , emit: cazyme_annotation
    tuple val(meta), path("${prefix}_dbCAN_hmm_results.tsv")   , emit: dbcanhmm_results
    tuple val(meta), path("${prefix}_dbCANsub_hmm_results.tsv"), emit: dbcansub_results
    tuple val(meta), path("${prefix}_diamond.out")             , emit: dbcandiamond_results
    tuple val(meta), path("${prefix}_cgc.gff")                 , emit: cgc_gff
    tuple val(meta), path("${prefix}_cgc_standard_out.tsv")    , emit: cgc_standard_out
    tuple val(meta), path("${prefix}_diamond.out.tc")          , emit: diamond_out_tc
    tuple val(meta), path("${prefix}_TF_hmm_results.tsv")      , emit: tf_hmm_results, optional: true
    tuple val(meta), path("${prefix}_STP_hmm_results.tsv")     , emit: stp_hmm_results
    tuple val(meta), path("${prefix}_total_cgc_info.tsv")      , emit: total_cgc_info
    tuple val(meta), path("${prefix}_substrate_prediction.tsv"), emit: substrate_prediction
    tuple val(meta), path("${prefix}_synteny_pdf/")            , emit: synteny_pdf
    path  "versions.yml"                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # run_dbcan 5.2.6's TCDB DIAMOND search ignores easy_substrate's --threads
    # and lets DIAMOND grab every core on the node (256 on the HPC compute
    # nodes), oversubscribing the ${task.cpus}-CPU allocation -> SIGBUS / multi-
    # hour hangs / NODE_FAIL. Shadow `diamond` on PATH so every DIAMOND call
    # run_dbcan spawns is pinned to task.cpus regardless of what it requests.
    mkdir -p .dbcan_shim
    export DBCAN_REAL_DIAMOND="\$(command -v diamond)"
    cat > .dbcan_shim/diamond <<'DBCAN_SHIM'
#!/usr/bin/env bash
args=()
while [ \$# -gt 0 ]; do
    case "\$1" in
        --threads|-p)     shift; shift ;;
        --threads=*|-p=*) shift ;;
        *)                args+=("\$1"); shift ;;
    esac
done
exec "\$DBCAN_REAL_DIAMOND" "\${args[@]}" --threads ${task.cpus}
DBCAN_SHIM
    chmod +x .dbcan_shim/diamond
    export PATH="\$PWD/.dbcan_shim:\$PATH"

    run_dbcan easy_substrate \\
        --mode protein \\
        --db_dir ${dbcan_db} \\
        --input_raw_data ${input_raw_data} \\
        --output_dir . \\
        --input_gff ${input_gff} \\
        --gff_type ${gff_type} \\
        --threads ${task.cpus} \\
        ${args}

    mv overview.tsv             ${prefix}_overview.tsv
    mv dbCAN_hmm_results.tsv    ${prefix}_dbCAN_hmm_results.tsv
    mv dbCANsub_hmm_results.tsv ${prefix}_dbCANsub_hmm_results.tsv
    mv diamond.out              ${prefix}_diamond.out
    mv cgc.gff                  ${prefix}_cgc.gff
    mv cgc_standard_out.tsv     ${prefix}_cgc_standard_out.tsv
    mv diamond.out.tc           ${prefix}_diamond.out.tc
    mv STP_hmm_results.tsv      ${prefix}_STP_hmm_results.tsv
    mv total_cgc_info.tsv       ${prefix}_total_cgc_info.tsv
    mv CGC.faa                  ${prefix}_CGC.faa
    mv PUL_blast.out            ${prefix}_PUL_blast.out
    mv substrate_prediction.tsv ${prefix}_substrate_prediction.tsv
    mv synteny_pdf/             ${prefix}_synteny_pdf/
    if [ -f TF_hmm_results.tsv ]; then
        mv TF_hmm_results.tsv   ${prefix}_TF_hmm_results.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: \$(echo \$(run_dbcan version) | cut -f2 -d':' | cut -f2 -d' ')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_overview.tsv
    touch ${prefix}_dbCAN_hmm_results.tsv
    touch ${prefix}_dbCANsub_hmm_results.tsv
    touch ${prefix}_diamond.out
    touch ${prefix}_cgc.gff
    touch ${prefix}_cgc_standard_out.tsv
    touch ${prefix}_diamond.out.tc
    touch ${prefix}_TF_hmm_results.tsv
    touch ${prefix}_STP_hmm_results.tsv
    touch ${prefix}_total_cgc_info.tsv
    touch ${prefix}_CGC.faa
    touch ${prefix}_PUL_blast.out
    touch ${prefix}_substrate_prediction.tsv
    mkdir -p ${prefix}_synteny_pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: \$(echo \$(run_dbcan version) | cut -f2 -d':' | cut -f2 -d' ')
    END_VERSIONS
    """
}
