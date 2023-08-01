process POSTDEMUX_SAMPSHEET {
    tag "$samplesheet"
    label 'process_low'
    container 'quay.io/biocontainers/pandas:1.5.2'

    input:
    path samplesheet
    tuple val(prefix1), path(reads)
    tuple val(prefix2), path(raw_data)
    val(obi3_demux)

    output:
    path '*.csv'       , emit: samplesheet
    path "*.txt"       , emit: missing_samples
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def samplesheet   = "'$samplesheet'"
    def reads         = "'$reads'"
    def raw_data      = "'$raw_data'"
    def obi3_demux    = "'$obi3_demux'"
    """
    #!/usr/bin/env python3

    # Setup
    import pandas as pd
    import os
    import sys
    py_version = sys.version
    py_version = py_version.split(" ")
    py_version = py_version[0]

    if (${obi3_demux}.upper() == "TRUE"):
        obi3_demux = True
    else:
        obi3_demux = False

    sampsheet = pd.read_csv(${samplesheet})
    raw_data = ${raw_data}.split(" ")
    missing_samples = []
    drop_rows = []

    for row in range(len(sampsheet.index)):
        # If the size is less than 25, the file is likely empty
        if (os.stat(str(sampsheet.at[row, "sample"]) + ".R1.fq.gz").st_size <= 25):
            missing_samples.append(sampsheet.at[row, "sample"])
            drop_rows.append(row)
        else:
            sampsheet.at[row, "fastq_1"] = str(os.getcwd() + "/" + sampsheet.at[row, "sample"] + ".R1.fq.gz")
        
            if len(raw_data) != 2 or obi3_demux:
                sampsheet.at[row, "fastq_2"] = None
            else:
                sampsheet.at[row, "fastq_2"] = str(os.getcwd() + "/" + sampsheet.at[row, "sample"] + ".R2.fq.gz")

    # Drop missing samples from df so they don't affect downstream modules
    if (len(drop_rows) > 0):
        sampsheet = sampsheet.drop(drop_rows)

    sampsheet.to_csv("new_" + ${samplesheet}, index=False)

    with open("missing_samples.txt", 'w') as file:
        if (len(missing_samples) > 0):
            for i in missing_samples:
                file.write(str(i) + '\\n')
        else:
            file.write("No samples were missing after demultiplexing")

    # Create version .yml file
    with open('versions.yml', 'w') as version_file:
        version_file.write(f"\\"{'${task.process}'}\\":\\n")
        version_file.write(f"    python: {py_version}\\n")
        version_file.write(f"    pandas: {pd.__version__}\\n")
    """
}
