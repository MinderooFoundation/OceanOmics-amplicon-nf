process TRIM_AND_CONCAT {
    tag ""
    label 'process_medium'
    container 'quay.io/biocontainers/csvtk:0.26.0--h9ee0642_0'

    input:
    tuple val(meta), path(reads)
    path metadata

    output:
    tuple val(meta), path("*.fq.gz"), emit: reads
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Function to trim and concatenate files
    process_files() {
        file1=\$1
        file2=\$2
        suffix=\$3
        forward_trim=\$4
        reverse_trim=\$5
        prefix=\$(basename "\$file2" | rev | cut -d '_' -f 2- | rev)

        if [[ "\$suffix" == "1.fq" ]]; then
            # Trim reads
            zcat "\$file1" | sed "n;s/^.\\{\$forward_trim\\}//" > "\${prefix}_trimmed_forward.\${suffix}"
            zcat "\$file2" | sed "n;s/^.\\{\$reverse_trim\\}//" > "\${prefix}_trimmed_reverse.\${suffix}"
        else
            # Trim reads
            zcat "\$file1" | sed "n;s/^.\\{\$reverse_trim\\}//" > "\${prefix}_trimmed_forward.\${suffix}"
            zcat "\$file2" | sed "n;s/^.\\{\$forward_trim\\}//" > "\${prefix}_trimmed_reverse.\${suffix}"
        fi

        # Concatenate files
        cat "\${prefix}_trimmed_forward.\${suffix}" "\${prefix}_trimmed_reverse.\${suffix}" | gzip > "\${prefix}.R\${suffix}.gz"

        # Remove files that aren't needed any more
        rm "\${prefix}_trimmed_forward.\${suffix}" "\${prefix}_trimmed_reverse.\${suffix}"
    }
    export -f process_files

    # Function to trim a single end file
    trim_file() {
        file1=\$1
        forward_trim=\$2
        prefix=\$3

        # Trim reads
        zcat "\$file1" | sed "n;s/^.\\{\$forward_trim\\}//" | gzip > "\${prefix}.R1.fq.gz"
    }
    export -f process_files

    # Loop through each row of the metadata, getting just the relevant columns
    cat ${metadata} | csvtk cut -f sample,fw_primer,rv_primer | {
        skip_first=true  # Flag variable to skip the first iteration

        while IFS= read -r line; do
            if \$skip_first; then
                skip_first=false
                continue  
            fi

            # Create a temporary file descriptor
            exec 3<<<"\$line"
            # Split the line
            IFS=',' read -r -u 3 sample fw_primer rv_primer

            # If single end
            if [ -z "\$rv_primer" ]; then
                fw_length=\${#fw_primer}

                trim_file \${sample}.1.fq.gz \$fw_length \${sample} &

            # paired end
            else
                fw_length=\${#fw_primer}
                rv_length=\${#rv_primer}

                # Trim and concatenate the files
                process_files \${sample}_forward.1.fq.gz \${sample}_reverse.1.fq.gz "1.fq" \$fw_length \$rv_length &

                # Trim and concatenate the files
                process_files \${sample}_forward.2.fq.gz \${sample}_reverse.2.fq.gz "2.fq" \$fw_length \$rv_length &
            fi
            
            # Close the temporary file descriptor
            exec 3>&-
        done
        wait
    }

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: 4.0
        csvtk: 0.26.0
    END_VERSIONS
    """
}