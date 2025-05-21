process INPUTFILE_INFO {
    tag "$meta"
    label 'process_low'
    container 'ubuntu:20.04'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.csv"), emit: csv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    sample=${meta.id}
    if [[ "${meta.single_end}" == "true" ]]; then
        md5sum="\$(md5sum ${reads})"
        checksum="\${md5sum%% *}"
        echo \${sample%_*},\${sample%_*},${reads},,\${checksum}, > \${sample}_file_info.csv
    else
        md5sum="\$(md5sum ${reads[0]})"
        md5sum2="\$(md5sum ${reads[1]})"
        checksum="\${md5sum%% *}"
        checksum2="\${md5sum2%% *}"
        echo \${sample%_*},\${sample%_*},${reads[0]},${reads[1]},\${checksum},\${checksum2} > \${sample}_file_info.csv
    fi
    """
}
