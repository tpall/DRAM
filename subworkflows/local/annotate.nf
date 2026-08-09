include { RENAME_FASTA           } from "../../modules/local/rename/rename_fasta.nf"
include { RENAME_PROTEINS        } from "../../modules/local/rename/rename_proteins.nf"
include { CALL                   } from "../../subworkflows/local/call.nf"
include { QC                     } from "../../subworkflows/local/qc.nf"
include { DB_SEARCH              } from "../../subworkflows/local/db_search.nf"
include { DRAMV_FLAGS            } from "../../modules/local/dramv/dramv_flags.nf"
include { GENE_LOCS              } from "../../modules/local/annotate/gene_locs.nf"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO ANNOTATE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ANNOTATE {
    take:
    ch_fasta  // channel: [ val(input_fasta name), path(fasta) ]
    default_sheet // Path to dummy sheet
    call     // boolean: whether gene calling flag is set
    use_kegg
    use_kofam
    use_dbcan
    use_camper
    use_fegenie
    use_methyl
    use_canthyd
    use_sulfur
    use_pfam
    use_merops
    use_uniref
    use_metals
    use_antismash
    use_rgi
    use_card
    use_tcdb
    use_dram_db
    use_vog

    main:
    // n_fastas = 0
    ch_rrna_collected = default_sheet
    ch_trna_collected = default_sheet
    ch_combined_annotations = default_sheet

    ch_quast_stats = default_sheet
    ch_collected_fna = default_sheet
    ch_gene_gff = default_sheet
    ch_filtered_fasta = default_sheet
    ch_called_genes = default_sheet

    if (call){
        fasta_name = ch_fasta.map { it[0] }
        fasta_files = ch_fasta.map { it[1] }

        // n_fastas = file("$params.input_fasta/${params.fasta_fmt}").size()

        if(params.rename) {
            // We need to use collect so that we pass all the fasta files to the rename process at once
            // Otherwise, it will try to rename each fasta file one at a time
            // Which since rename is so fast, will clog up job queues
            // so it is faster to rename all at once
            RENAME_FASTA( fasta_name.toList(), fasta_files.toList() )
            // we use flatten here to turn a list back into a channel
            renamed_fasta_paths = RENAME_FASTA.out.renamed_fasta_paths.flatten()
            // we need to recreate the fasta channel with the renamed fasta files
            ch_fasta = renamed_fasta_paths.map {
                fasta_name = it.getBaseName()
                tuple(fasta_name, it)
            }
        }

        CALL( ch_fasta )
        ch_quast_stats = CALL.out.ch_quast_stats
        ch_gene_locs = CALL.out.ch_gene_locs
        ch_called_proteins = CALL.out.ch_called_proteins
        ch_collected_fna = CALL.out.ch_collected_fna
        ch_gene_gff = CALL.out.ch_gene_gff
        ch_filtered_fasta = CALL.out.ch_filtered_fasta
        ch_called_genes = CALL.out.ch_called_genes

    }
    else {
        ch_called_proteins = channel
            .fromPath(file(params.input_genes) / params.genes_fmt, checkIfExists: true)
            .ifEmpty { exit 1, "If you specify --annotate without --call, you must provide a fasta file of called genes using --input_genes. Cannot find any called gene fasta files matching: ${params.input_genes}\n Path needs to follow pattern: path/to/directory/" }
            .map { file ->
                def input_fastaName = file.getBaseName()
                tuple(input_fastaName, file)
            }

        fasta_name = ch_called_proteins.map { tup -> tup[0] }
        fasta_files = ch_called_proteins.map { tup -> tup[1] }

        if(params.rename) {
            // We need to use collect so that we pass all the fasta files to the rename process at once
            // Otherwise, it will try to rename each fasta file one at a time
            // Which since rename is so fast, will clog up job queues
            // so it is faster to rename all at once
            RENAME_PROTEINS( fasta_name.toList(), fasta_files.toList() )
            // we use flatten here to turn a list back into a channel
            renamed_paths = RENAME_PROTEINS.out.renamed_paths.flatten()
            // we need to recreate the fasta channel with the renamed fasta files
            ch_called_proteins = renamed_paths.map {
                fasta_name = it.getBaseName()
                tuple(fasta_name, it)
            }
        }

        GENE_LOCS( ch_called_proteins)
        ch_gene_locs = GENE_LOCS.out.prodigal_locs_tsv
        // n_fastas = file("$params.input_genes/${params.genes_fmt}").size()
    }
    ch_antismash_map = ch_filtered_fasta
        .map { file ->
            def meta = [:]
            meta.id = file.getBaseName()
            tuple(meta, file)
        }

    ch_fna_map = ch_called_genes
        .map {
                file_name, file ->
                def meta = [:]
                meta.id = file_name
                tuple(meta, file)
            }

    ch_faa_map = ch_called_proteins
        .map {
                file_name, file ->
                def meta = [:]
                meta.id = file_name
                tuple(meta, file)
            }

    ch_gff_map = ch_gene_gff
        .map {
            file ->
            def meta = [:]
            meta.id = file.getBaseName()
            tuple(meta, file, "prodigal")
        }

    if (params.annotate){
        DB_SEARCH(
            ch_gene_locs,
            ch_called_proteins,
            ch_antismash_map,
            ch_fna_map,
            ch_faa_map,
            ch_gff_map,
            ch_gene_gff,
            default_sheet,
            use_kegg,
            use_kofam,
            use_dbcan,
            use_camper,
            use_fegenie,
            use_methyl,
            use_canthyd,
            use_sulfur,
            use_pfam,
            use_merops,
            use_uniref,
            use_metals,
            use_antismash,
            use_rgi,
            use_card,
            use_tcdb,
            use_dram_db,
            use_vog
            )
        ch_combined_annotations = DB_SEARCH.out.ch_combined_annotations
    }

    if (params.qc){
        QC( ch_fasta, default_sheet, ch_combined_annotations, ch_collected_fna, call )
        ch_rrna_collected = QC.out.ch_rrna_collected
        ch_trna_collected = QC.out.ch_trna_collected
        ch_combined_annotations = QC.out.ch_final_annots
    }

    // DRAM-v: add amg_flags + is_transposon to the combined annotations.
    if (params.use_dramv) {
        // use_dramv force-enables use_vog in the workflow, so vog_list is always
        // available here; the V flag reads VOGdb functional categories from it.
        ch_vog_list = channel.fromPath(params.vog_list, checkIfExists: true)
        DRAMV_FLAGS(
            ch_combined_annotations,
            ch_fasta.map { it -> it[1] }.collect(),
            ch_vog_list
        )
        ch_combined_annotations = DRAMV_FLAGS.out.combined_annotations_with_flags
    }

    emit:
    ch_rrna_collected
    ch_trna_collected
    ch_combined_annotations
    ch_quast_stats

}
