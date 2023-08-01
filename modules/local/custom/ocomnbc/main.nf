process OCOMNBC {
    tag "$prefix"
    label 'process_high'
    container 'adbennett/ocom_nbc:v1.3'

    input:
    tuple val(prefix), path(table)

    output:
    tuple val(prefix), path("*.tsv")          , emit: nbc_output
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${prefix}"
    def table = "${table}"
    """
    #!/usr/bin/env python3

    # Setup
    import csv
    import sys
    import os
    import joblib
    output_header = ['NBC_predicted_species', 'NBC_predicted_scores']
    seq_column = '${prefix}'.upper() + '_sequence'
    py_version = sys.version
    py_version = py_version.split(" ")
    py_version = py_version[0]

    #@st.cache_data
    def get_prediction_and_scores(seqs):
        ''' loads the scikit-learn pipeline, returns predictions and score'''
        p = species_pipeline.predict(seqs)
        all_scores = species_pipeline.predict_proba(seqs)
        all_classes = list(species_pipeline.classes_)
        predicted_scores = []
 
        for prediction, scores in zip(p, all_scores):
            predicted_scores.append(scores[all_classes.index(prediction)])
        return p, predicted_scores

    # Load the pipeline binary from the docker image
    species_pipeline = joblib.load('/mnt/scratch/pipeline.pkl') 

    # Get a list of sequences from input .tsv file
    seqs=[]
    with open('${table}') as input:
        reader = csv.DictReader(input, delimiter="\\t", quotechar='"')

        for row in reader:
            seqs.append(row[seq_column])

    # Get list of NBC species and NBC scores
    species_list, scores_list = get_prediction_and_scores(seqs)

    # Create output .tsv file
    out_file = '${prefix}' + '_nbc_output.tsv'
    with open(out_file, 'w', newline='') as output:
        writer = csv.writer(output, delimiter='\\t')
        writer.writerow(output_header)
        writer.writerows(zip(species_list, scores_list))

    # Create version .yml file
    with open('versions.yml', 'w') as version_file:
        version_file.write(f"\\"{'${task.process}'}\\":\\n")
        version_file.write(f"    python: {py_version}\\n")
        version_file.write(f"    joblib: {joblib.__version__}\\n")
    """
}
