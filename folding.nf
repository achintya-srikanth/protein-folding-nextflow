nextflow.enable.dsl=2

// Default parameters
params.fasta = 'data/test.fasta'
params.outdir = 'results'

// Processes

process get_msa {

    tag "${id}"

    input:
    tuple val(id), val(sequence)

    output:
    tuple val(id), path("${id}_sequence_msa.a3m")

    script:
    """
    get_msa.py --id ${id} --sequence "${sequence}"
    """
}

process run_alphafold {
    label 'gpu'
    tag "${id}"

    input:
    tuple val(id), path(msa), val(sequence)

    output:
    tuple val(id), path("${id}_alphafold_result.pdb")

    script:
    """
    # 1. Sanitize the MSA file: remove Windows carriage returns (\r) and any leading blank lines/spaces
    tr -d '\\r' < "${msa}" | sed '/^\$/d' > sanitized_msa.a3m

    # 2. Re-materialize the clean query fasta 
    echo ">${id}" > temp.fasta
    echo "${sequence.trim()}" >> temp.fasta
    
    # 3. Pass the sanitized versions to your script
    run_alphafold.py --id "${id}" --msa sanitized_msa.a3m --fasta temp.fasta
    """
}

process run_esmfold {
    label 'gpu'
    
    tag "${id}"
    
    input:
    tuple val(id), val(sequence)

    output:
    tuple val(id), path("${id}_esmfold_result.pdb")

    script:
    """
    run_esmfold.py --id ${id} --sequence "${sequence.trim()}"
    """
}

process run_omegafold {
    label 'gpu'

    input:
    tuple val(id), val(sequence)

    output:
    tuple val(id), path("${id}_omegafold_result.pdb")

    script:
    """
    run_omegafold.py --id ${id} --sequence "${sequence.trim()}"
    """
}

process run_rosettafold {
    label 'gpu'

    tag "${id}"

    input:
    tuple val(id), val(sequence)

    output:
    tuple val(id), path("${id}_rosettafold_result.pdb")

    script:
    """
    run_rosettafold.py --id ${id} --sequence "${sequence.trim()}"
    """
}

process publish_results {

    tag "${id}"

    publishDir {"${params.outdir}/${id}"}, mode: 'copy'

    input:
    tuple val(id), path('inputs/*', stageAs: 'inputs/*')

    output:
    path "${id}_summary.txt"
    path "*.pdb"

    script:
    """
    echo "Summary of models for ${id}:" > ${id}_summary.txt

    for f in inputs/*; do
        echo "- \$(basename \$f)" >> ${id}_summary.txt
        cp \$f ./ # Re-materialize the file in the task directory for publishing
    done

    echo "Total models generated: \$(ls *.pdb | wc -l)" >> ${id}_summary.txt
    
    """

}

// Define the workflow
workflow {

    // Take fasta file as input and split into records (id and sequence)

    input_ch = Channel.fromPath(params.fasta)
                      .splitFasta(record: ['id': true, 'sequence': true])
                      .map { record -> [record.id, record.sequence] }
    
    // Send fasta file to MMseqs2 for MSA search

    msa_channel = get_msa(input_ch)

    // Join a3m with fasta on {id} for af2

    af2_input_ch = msa_channel.join(input_ch)

    // Parallel execution across models (Broadcasting)

    alphafold_pdb = run_alphafold(af2_input_ch)
    esmfold_pdb = run_esmfold(input_ch)
    omegafold_pdb = run_omegafold(input_ch)
    rosettafold_pdb = run_rosettafold(input_ch)

    // Gather results and publish to 'results'

    all_models = alphafold_pdb.mix(omegafold_pdb, esmfold_pdb, rosettafold_pdb).groupTuple()

    // Publish results
    publish_results(all_models)

}