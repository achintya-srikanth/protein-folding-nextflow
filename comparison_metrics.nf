#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.fasta       = "data/mutations.fasta"
params.staged_json = "staged_pairs.json"
params.outdir      = "comparison_reports"

process fetch_uniprot_annotations {
    tag "${pdb_id}"
    
    input:
    val pdb_id

    output:
    tuple val(pdb_id), stdout

    script:
    """
    get_active_sites.py ${pdb_id} | tr -d '\\n'
    """
}

process compute_structural_metrics {
    tag "${meta.id}"
    publishDir {"${params.outdir}/${meta.id}/data"}, mode: 'copy', pattern: "*_metrics_report.json"
    
    input:
    tuple val(meta), path(wt_pdb), path(mt_pdb), path(fasta), val(active_residues)

    output:
    tuple val(meta), path(wt_pdb), path(mt_pdb), path("${meta.id}_metrics_report.json")

    script:
    def active_sites_arg = active_residues.trim() ? "--active_sites ${active_residues.trim()}" : ""
    """
    analyze_mutations.py \
        --wt_pdb "${wt_pdb}" \
        --mt_pdb "${mt_pdb}" \
        --fasta "${fasta}" \
        --wt_id "${meta.wt_id}" \
        --mt_id "${meta.mt_id}" \
        ${active_sites_arg} > "${meta.id}_metrics_report.json"
    """
}

process render_structural_visuals {
    tag "${meta.id}"
    publishDir {"${params.outdir}/${meta.id}/images"}, mode: 'copy'

    input:
    tuple val(meta), path(wt_pdb), path(mt_pdb), path(metrics_json)

    output:
    path "${meta.id}_mutation_overlay.png"

    script:
    """
    render_mutations.py \
        --wt_pdb "${wt_pdb}" \
        --mt_pdb "${mt_pdb}" \
        --json_report "${metrics_json}" \
        --id "${meta.id}"
    """
}

workflow {
    fasta_ref = file(params.fasta)
    
    // 1. Ingest the 12 pairs from the JSON manifest
    staged_ch = Channel.fromPath(params.staged_json)
        .splitJson()
        .map { entry ->
            def meta = [
                id: entry.meta_id,
                wt_id: entry.wt_id,
                mt_id: entry.mt_id,
                model: entry.model
            ]
            return [ entry.pdb_id, meta, file(entry.wt_path), file(entry.mt_path) ]
        }

    // 2. Query annotations uniquely per PDB ID (Emits 3 items total)
    unique_pdb_ids_ch = staged_ch.map { pdb_id, meta, wt, mt -> pdb_id }.unique()
    annotation_stream = fetch_uniprot_annotations(unique_pdb_ids_ch)

    // 3. EXHAUSTIVE MULTIPLEX FIX: Use .reduce to build a native, stable Groovy dictionary
    // This explicitly maps each [pdb_id, residues] into a single, permanent broadcast map
    annotation_map_ch = annotation_stream
        .reduce([:]) { map, item -> 
            def pdb_id = item[0]
            def residues = item[1]
            map[pdb_id] = residues
            return map
        }

    // 4. Combine the 12 staged items with the 1 global lookup map Value Channel
    processing_ch = staged_ch
        .combine(annotation_map_ch)
        .map { pdb_id, meta, wt_file, mt_file, global_map ->
            
            // Clean, syntax-safe lookup from our native broadcast dictionary
            def active_residues = global_map[pdb_id] ?: ""
            
            return [ meta, wt_file, mt_file, fasta_ref, active_residues ] 
        }

    // 5. Run calculation engines and visual image generation across all 12 pairs in parallel
    metrics_ch = compute_structural_metrics(processing_ch)
    render_structural_visuals(metrics_ch)
}