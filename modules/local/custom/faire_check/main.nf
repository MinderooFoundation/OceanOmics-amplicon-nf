process FAIRE_CHECK {
    tag "$metadata"
    label 'process_low'
    container 'quay.io/biocontainers/pandas:1.5.2'

    input:
    path(metadata)

    output:
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    #!/usr/bin/env python3

    # Setup
    import pandas as pd
    import os
    import sys
    py_version = sys.version
    py_version = py_version.split(" ")
    py_version = py_version[0]


    #sampsheet = pd.read_csv(${metadata})
    #raw_data = ${metadata}.split(" ")
    #missing_samples = []
    #discarded_samples = []
    #drop_rows = []

    #for row in range(len(sampsheet.index)):
    #    if "discarded" in sampsheet.columns:
    #        if (sampsheet.at[row, "discarded"] == True):
    #            discarded_samples.append(sampsheet.at[row, "sample"])

        # If the size is less than 25, the file is likely empty
    #    if not os.path.isfile(str(sampsheet.at[row, "sample"]) + ".R1.fq.gz"):
    #        missing_samples.append(sampsheet.at[row, "sample"])
    #        drop_rows.append(row)
    #    elif (os.stat(str(sampsheet.at[row, "sample"]) + ".R1.fq.gz").st_size <= 25):
    #        missing_samples.append(sampsheet.at[row, "sample"])
    #        drop_rows.append(row)
    #    else:
    #        sampsheet.at[row, "fastq_1"] = str(os.getcwd() + "/" + sampsheet.at[row, "sample"] + ".R1.fq.gz")

    #        if len(raw_data) != 2 or obi3_demux:
    #            sampsheet.at[row, "fastq_2"] = None
    #        else:
    #            sampsheet.at[row, "fastq_2"] = str(os.getcwd() + "/" + sampsheet.at[row, "sample"] + ".R2.fq.gz")

    #if "discarded" in sampsheet.columns:
    #    disc_miss_df = pd.DataFrame(data={"samples_discarded_post_qPCR_QC": discarded_samples})
    #    disc_miss_df["samples_with_no_reads_assigned"] = ""
    #    for row in range(len(disc_miss_df.index)):
    #        if (disc_miss_df.at[row, "samples_discarded_post_qPCR_QC"]) in missing_samples:
    #            disc_miss_df.at[row, "samples_with_no_reads_assigned"] = disc_miss_df.at[row, "samples_discarded_post_qPCR_QC"]

    #    for sam in missing_samples:
    #        if (sam not in discarded_samples):
    #            new_row = pd.DataFrame(data={"samples_discarded_post_qPCR_QC": [""], "samples_with_no_reads_assigned": [sam]})
    #            disc_miss_df = pd.concat([disc_miss_df, new_row])
    #    if len(disc_miss_df) > 0:
    #        disc_miss_df.to_csv("missing_samples.csv", index=False)
    #    else:
    #        with open("missing_samples.csv", 'w') as file:
    #            file.write("No_discarded_samples_and_no_samples_missing_after_demultiplexing")
    #else:
    #    with open("missing_samples.csv", 'w') as file:
    #        if (len(missing_samples) > 0):
    #            file.write('Missing_Samples\\n')
    #            for i in missing_samples:
    #                file.write(str(i) + '\\n')
    #        else:
    #            file.write("No_samples_were_missing_after_demultiplexing")

    # Drop missing samples from df so they don't affect downstream modules
    #if (len(drop_rows) > 0):
    #    sampsheet = sampsheet.drop(drop_rows)

    #sampsheet.to_csv("new_" + ${metadata}, index=False)

    # Create version .yml file
    with open('versions.yml', 'w') as version_file:
        version_file.write(f"\\"{'${task.process}'}\\":\\n")
        version_file.write(f"    python: {py_version}\\n")
        version_file.write(f"    pandas: {pd.__version__}\\n")
    """
}
