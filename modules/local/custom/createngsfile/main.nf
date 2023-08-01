process CREATE_NGSFILE {
    tag "$samplesheet"
    label 'process_medium'
    container 'quay.io/biocontainers/csvtk:0.26.0--h9ee0642_0'

    input:
    path samplesheet

    output:
    path "*.txt"       , emit: ngsfile
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    touch ngsfile.txt
    cat ${samplesheet} | csvtk cut -f sample,fw_index,rv_index,fw_primer,rv_primer | {
        skip_first=true  # Flag variable to skip the first iteration

        while IFS= read -r line; do
            if \$skip_first; then
                skip_first=false
                continue  
            fi

            # Create a temporary file descriptor
            exec 3<<<"\$line"
            # Split the line
            IFS=',' read -r -u 3 sample fw_index rv_index fw_primer rv_primer

            if [ -z "\$rv_primer" ] || [ -z "\$rv_index" ]; then
                echo "rv_index and rv_primer are needed for obitools3 workflow, try cutadapt workflow instead"
                exit 1

            else
                fw_row="ngsfile\\t\${sample}\\t\${fw_index}:\${rv_index}\\t\${fw_primer}\\t\${rv_primer}"
                rv_row="ngsfile\\t\${sample}\\t\${rv_index}:\${fw_index}\\t\${rv_primer}\\t\${fw_primer}"
                echo -e "\${fw_row}" >> ngsfile.txt
                echo -e "\${rv_row}" >> ngsfile.txt
            fi

            # Close the temporary file descriptor
            exec 3>&-
        done
    }
 
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: 0.26.0
    END_VERSIONS
    """
}
