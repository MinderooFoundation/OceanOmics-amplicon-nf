process CUTADAPT_TRIM_LEFTOVER_PRIMERS {
    tag "$prefix"
    label 'process_medium'
    container 'quay.io/biocontainers/cutadapt:4.7--py310h4b81fae_1'

    input:
    tuple val(prefix), path(reads)
    path primer_file_5end
    path primer_file_3end
    val ulimit

    output:
    tuple val(prefix), path("trimmed_*.gz"), emit: reads
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args1 = task.ext.args1 ?: ''
    def args2 = task.ext.args2 ?: ''
    """
    #!/bin/bash

    # Too avoid too many open files error:
    ulimit -S -n ${ulimit}

    files=($reads)
    if [ "\${#files[@]}" -eq 1 ]; then
        last_line=""
        first_primer=true
        while IFS= read -r line; do
            if [[ ! \$line == ">"* ]]; then
                if [ \$first_primer = true ]; then
                    cutadapt -j ${task.cpus} \
                        -g "\$line" \
                        ${args1} \
                        -o \${line}_\$(basename \${files[0]}) \
                        ${reads}
                    last_line=\${line}
                    first_primer=false
                else
                    cutadapt -j ${task.cpus} \
                        -g "\$line" \
                        ${args1} \
                        -o \${line}_\$(basename \${files[0]}) \
                        \${last_line}_\$(basename \${files[0]})
                    last_line=\${line}
                fi
            fi
        done < ${primer_file_5end}

        while IFS= read -r line; do
            if [[ ! \$line == ">"* ]]; then
                cutadapt -j ${task.cpus} \
                    -a "\$line" \
                    ${args1} \
                    -o \${line}_\$(basename \${files[0]}) \
                    \${last_line}_\$(basename \${files[0]})
                last_line=\${line}
            fi
        done < ${primer_file_3end}
        mv \${last_line}_\$(basename \${files[0]}) trimmed_\$(basename \${files[0]})

    else
        last_line=""
        first_primer=true
        while IFS= read -r line; do
            if [[ ! \$line == ">"* ]]; then
                if [ \$first_primer = true ]; then
                    cutadapt -j ${task.cpus} \
                        -g "\$line" \
                        -G "\$line" \
                        ${args2} \
                        -o \${line}_\$(basename \${files[0]}) -p \${line}_\$(basename \${files[1]}) \
                        ${reads}
                    last_line=\${line}
                    first_primer=false
                else
                    cutadapt -j ${task.cpus} \
                        -g "\$line" \
                        -G "\$line" \
                        ${args2} \
                        -o \${line}_\$(basename \${files[0]}) -p \${line}_\$(basename \${files[1]}) \
                        \${last_line}_\$(basename \${files[0]}) \${last_line}_\$(basename \${files[1]})
                    last_line=\${line}
                fi
            fi
        done < ${primer_file_5end}

        while IFS= read -r line; do
            if [[ ! \$line == ">"* ]]; then
                cutadapt -j ${task.cpus} \
                    -a "\$line" \
                    -A "\$line" \
                    ${args2} \
                    -o \${line}_\$(basename \${files[0]}) -p \${line}_\$(basename \${files[1]}) \
                    \${last_line}_\$(basename \${files[0]}) \${last_line}_\$(basename \${files[1]})
                last_line=\${line}
            fi
        done < ${primer_file_3end}

        mv \${last_line}_\$(basename \${files[0]}) trimmed_\$(basename \${files[0]})
        mv \${last_line}_\$(basename \${files[1]}) trimmed_\$(basename \${files[1]})
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}
