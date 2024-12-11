process PRIMER_CONTAM_STATS {
    tag "$prefix"
    label 'process_low'

    container 'ubuntu:20.04'

    input:
    tuple val(prefix), path(fasta)
    val fw_primer
    val rv_primer

    output:
    tuple val(prefix), path("*.txt"), emit: txt
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    touch ${prefix}_primer_contam_stats.txt
    fw_prim="${fw_primer}"
    rv_prim="${rv_primer}"

    if [[ \${fw_prim} == *";"* ]]; then
        fw_primer_arr=(\${fw_prim//;/ })
        fw_arraylength=\${#fw_primer_arr[@]}

        for (( i=0; i<\${fw_arraylength}; i++ ));
        do
            fw_prim_firsthalf=\${fw_primer_arr[\$i]:0:\${#fw_primer_arr[\$i]}/2}
            fw_prim_secondhalf=\${fw_primer_arr[\$i]:\${#fw_primer_arr[\$i]}/2}
            fw_prim_rev=\$(echo \${fw_primer_arr[\$i]} | rev | sed 's/C/g/g' | sed 's/G/c/g' | sed 's/A/t/g' | sed 's/T/a/g' | sed 's/c/C/g' | sed 's/g/G/g' | sed 's/t/T/g' | sed 's/a/A/g')
            fw_prim_firsthalf_rev=\$(echo \${fw_prim_firsthalf} | rev | sed 's/C/g/g' | sed 's/G/c/g' | sed 's/A/t/g' | sed 's/T/a/g' | sed 's/c/C/g' | sed 's/g/G/g' | sed 's/t/T/g' | sed 's/a/A/g')
            fw_prim_secondhalf_rev=\$(echo \${fw_prim_secondhalf} | rev | sed 's/C/g/g' | sed 's/G/c/g' | sed 's/A/t/g' | sed 's/T/a/g' | sed 's/c/C/g' | sed 's/g/G/g' | sed 's/t/T/g' | sed 's/a/A/g')

            echo "fw prim \${fw_primer_arr[\$i]} found \$(grep \${fw_primer_arr[\$i]} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
            if [ \$(grep \${fw_primer_arr[\$i]} ${fasta} | wc -l) -ne 0 ]; then
                echo "<details>" >> ${prefix}_primer_contam_stats.txt
                echo "  <summary>ASVs with the fw primer</summary>" >> ${prefix}_primer_contam_stats.txt
                echo "  <p>\$(grep \${fw_primer_arr[\$i]} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
                echo "</details>" >> ${prefix}_primer_contam_stats.txt
            fi

            echo "<br>fw prim first half \${fw_prim_firsthalf} found \$(grep \${fw_prim_firsthalf} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
            if [ \$(grep \${fw_prim_firsthalf} ${fasta} | wc -l) -ne 0 ]; then
                echo "<details>" >> ${prefix}_primer_contam_stats.txt
                echo "  <summary>ASVs with the first half of the fw primer</summary>" >> ${prefix}_primer_contam_stats.txt
                echo "  <p>\$(grep \${fw_prim_firsthalf} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
                echo "</details>" >> ${prefix}_primer_contam_stats.txt
            fi

            echo "<br>fw prim second half \${fw_prim_secondhalf} found \$(grep \${fw_prim_secondhalf} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
            if [ \$(grep \${fw_prim_secondhalf} ${fasta} | wc -l) -ne 0 ]; then
                echo "<details>" >> ${prefix}_primer_contam_stats.txt
                echo "  <summary>ASVs with the second half of the fw primer</summary>" >> ${prefix}_primer_contam_stats.txt
                echo "  <p>\$(grep \${fw_prim_secondhalf} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
                echo "</details>" >> ${prefix}_primer_contam_stats.txt
            fi

            echo "<br>fw prim reverse-complemented \${fw_prim_rev} found \$(grep \${fw_prim_rev} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
            if [ \$(grep \${fw_prim_rev} ${fasta} | wc -l) -ne 0 ]; then
                echo "<details>" >> ${prefix}_primer_contam_stats.txt
                echo "  <summary>ASVs with fw primer reversed</summary>" >> ${prefix}_primer_contam_stats.txt
                echo "  <p>\$(grep \${fw_prim_rev} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
                echo "</details>" >> ${prefix}_primer_contam_stats.txt
            fi

            echo "<br>fw prim first half reverse-complemented \${fw_prim_firsthalf_rev} found \$(grep \${fw_prim_firsthalf_rev} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
            if [ \$(grep \${fw_prim_firsthalf_rev} ${fasta} | wc -l) -ne 0 ]; then
                echo "<details>" >> ${prefix}_primer_contam_stats.txt
                echo "  <summary>ASVs with the first half of the fw primer reversed</summary>" >> ${prefix}_primer_contam_stats.txt
                echo "  <p>\$(grep \${fw_prim_firsthalf_rev} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
                echo "</details>" >> ${prefix}_primer_contam_stats.txt
            fi

            echo "<br>fw prim second half reverse-complemented \${fw_prim_secondhalf_rev} found \$(grep \${fw_prim_secondhalf_rev} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
            if [ \$(grep \${fw_prim_secondhalf_rev} ${fasta} | wc -l) -ne 0 ]; then
                echo "<details>" >> ${prefix}_primer_contam_stats.txt
                echo "  <summary>ASVs with the second half of the fw primer reversed</summary>" >> ${prefix}_primer_contam_stats.txt
                echo "  <p>\$(grep \${fw_prim_secondhalf_rev} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
                echo "</details>" >> ${prefix}_primer_contam_stats.txt
            fi
        done
    else
        fw_prim_firsthalf=\${fw_prim:0:\${#fw_prim}/2}
        fw_prim_secondhalf=\${fw_prim:\${#fw_prim}/2}
        fw_prim_rev=\$(echo \${fw_prim} | rev | sed 's/C/g/g' | sed 's/G/c/g' | sed 's/A/t/g' | sed 's/T/a/g' | sed 's/c/C/g' | sed 's/g/G/g' | sed 's/t/T/g' | sed 's/a/A/g')
        fw_prim_firsthalf_rev=\$(echo \${fw_prim_firsthalf} | rev | sed 's/C/g/g' | sed 's/G/c/g' | sed 's/A/t/g' | sed 's/T/a/g' | sed 's/c/C/g' | sed 's/g/G/g' | sed 's/t/T/g' | sed 's/a/A/g')
        fw_prim_secondhalf_rev=\$(echo \${fw_prim_secondhalf} | rev | sed 's/C/g/g' | sed 's/G/c/g' | sed 's/A/t/g' | sed 's/T/a/g' | sed 's/c/C/g' | sed 's/g/G/g' | sed 's/t/T/g' | sed 's/a/A/g')

        echo "fw prim \${fw_prim} found \$(grep \${fw_prim} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
        if [ \$(grep \${fw_prim} ${fasta} | wc -l) -ne 0 ]; then
            echo "<details>" >> ${prefix}_primer_contam_stats.txt
            echo "  <summary>ASVs with the fw primer</summary>" >> ${prefix}_primer_contam_stats.txt
            echo "  <p>\$(grep \${fw_prim} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
            echo "</details>" >> ${prefix}_primer_contam_stats.txt
        fi

        echo "<br>fw prim first half \${fw_prim_firsthalf} found \$(grep \${fw_prim_firsthalf} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
        if [ \$(grep \${fw_prim_firsthalf} ${fasta} | wc -l) -ne 0 ]; then
            echo "<details>" >> ${prefix}_primer_contam_stats.txt
            echo "  <summary>ASVs with the first half of the fw primer</summary>" >> ${prefix}_primer_contam_stats.txt
            echo "  <p>\$(grep \${fw_prim_firsthalf} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
            echo "</details>" >> ${prefix}_primer_contam_stats.txt
        fi

        echo "<br>fw prim second half \${fw_prim_secondhalf} found \$(grep \${fw_prim_secondhalf} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
        if [ \$(grep \${fw_prim_secondhalf} ${fasta} | wc -l) -ne 0 ]; then
            echo "<details>" >> ${prefix}_primer_contam_stats.txt
            echo "  <summary>ASVs with the second half of the fw primer</summary>" >> ${prefix}_primer_contam_stats.txt
            echo "  <p>\$(grep \${fw_prim_secondhalf} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
            echo "</details>" >> ${prefix}_primer_contam_stats.txt
        fi

        echo "<br>fw prim reverse-complemented \${fw_prim_rev} found \$(grep \${fw_prim_rev} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
        if [ \$(grep \${fw_prim_rev} ${fasta} | wc -l) -ne 0 ]; then
            echo "<details>" >> ${prefix}_primer_contam_stats.txt
            echo "  <summary>ASVs with fw primer reversed</summary>" >> ${prefix}_primer_contam_stats.txt
            echo "  <p>\$(grep \${fw_prim_rev} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
            echo "</details>" >> ${prefix}_primer_contam_stats.txt
        fi

        echo "<br>fw prim first half reverse-complemented \${fw_prim_firsthalf_rev} found \$(grep \${fw_prim_firsthalf_rev} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
        if [ \$(grep \${fw_prim_firsthalf_rev} ${fasta} | wc -l) -ne 0 ]; then
            echo "<details>" >> ${prefix}_primer_contam_stats.txt
            echo "  <summary>ASVs with the first half of the fw primer reversed</summary>" >> ${prefix}_primer_contam_stats.txt
            echo "  <p>\$(grep \${fw_prim_firsthalf_rev} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
            echo "</details>" >> ${prefix}_primer_contam_stats.txt
        fi

        echo "<br>fw prim second half reverse-complemented \${fw_prim_secondhalf_rev} found \$(grep \${fw_prim_secondhalf_rev} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
        if [ \$(grep \${fw_prim_secondhalf_rev} ${fasta} | wc -l) -ne 0 ]; then
            echo "<details>" >> ${prefix}_primer_contam_stats.txt
            echo "  <summary>ASVs with the second half of the fw primer reversed</summary>" >> ${prefix}_primer_contam_stats.txt
            echo "  <p>\$(grep \${fw_prim_secondhalf_rev} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
            echo "</details>" >> ${prefix}_primer_contam_stats.txt
        fi
    fi

    if [[ \${rv_prim} == *";"* ]]; then
        rv_primer_arr=(\${rv_prim//;/ })
        rv_arraylength=\${#rv_primer_arr[@]}

        for (( i=0; i<\${rv_arraylength}; i++ ));
        do
            rv_prim_firsthalf=\${rv_primer_arr[\$i]:0:\${#rv_primer_arr[\$i]}/2}
            rv_prim_secondhalf=\${rv_primer_arr[\$i]:\${#rv_primer_arr[\$i]}/2}

            rv_prim_rev=\$(echo \${rv_primer_arr[\$i]} | rev | sed 's/C/g/g' | sed 's/G/c/g' | sed 's/A/t/g' | sed 's/T/a/g' | sed 's/c/C/g' | sed 's/g/G/g' | sed 's/t/T/g' | sed 's/a/A/g')
            rv_prim_firsthalf_rev=\$(echo \${rv_prim_firsthalf} | rev | sed 's/C/g/g' | sed 's/G/c/g' | sed 's/A/t/g' | sed 's/T/a/g' | sed 's/c/C/g' | sed 's/g/G/g' | sed 's/t/T/g' | sed 's/a/A/g')
            rv_prim_secondhalf_rev=\$(echo \${rv_prim_secondhalf} | rev | sed 's/C/g/g' | sed 's/G/c/g' | sed 's/A/t/g' | sed 's/T/a/g' | sed 's/c/C/g' | sed 's/g/G/g' | sed 's/t/T/g' | sed 's/a/A/g')

            echo "<br>rv prim \${rv_primer_arr[\$i]} found \$(grep \${rv_primer_arr[\$i]} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
            if [ \$(grep \${rv_primer_arr[\$i]} ${fasta} | wc -l) -ne 0 ]; then
                echo "<details>" >> ${prefix}_primer_contam_stats.txt
                echo "  <summary>ASVs with the rv primer</summary>" >> ${prefix}_primer_contam_stats.txt
                echo "  <p>\$(grep \${rv_primer_arr[\$i]} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
                echo "</details>" >> ${prefix}_primer_contam_stats.txt
            fi

            echo "<br>rv prim first half \${rv_prim_firsthalf} found \$(grep \${rv_prim_firsthalf} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
            if [ \$(grep \${rv_prim_firsthalf} ${fasta} | wc -l) -ne 0 ]; then
                echo "<details>" >> ${prefix}_primer_contam_stats.txt
                echo "  <summary>ASVs with the first half of the rv primer</summary>" >> ${prefix}_primer_contam_stats.txt
                echo "  <p>\$(grep \${rv_prim_firsthalf} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
                echo "</details>" >> ${prefix}_primer_contam_stats.txt
            fi

            echo "<br>rv prim second half \${rv_prim_secondhalf} found \$(grep \${rv_prim_secondhalf} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
            if [ \$(grep \${rv_prim_secondhalf} ${fasta} | wc -l) -ne 0 ]; then
                echo "<details>" >> ${prefix}_primer_contam_stats.txt
                echo "  <summary>ASVs with the second half of the rv primer</summary>" >> ${prefix}_primer_contam_stats.txt
                echo "  <p>\$(grep \${rv_prim_secondhalf} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
                echo "</details>" >> ${prefix}_primer_contam_stats.txt
            fi

            echo "<br>rv prim reverse-complemented \${rv_prim_rev} found \$(grep \${rv_prim_rev} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
            if [ \$(grep \${rv_prim_rev} ${fasta} | wc -l) -ne 0 ]; then
                echo "<details>" >> ${prefix}_primer_contam_stats.txt
                echo "  <summary>ASVs with rv primer reversed</summary>" >> ${prefix}_primer_contam_stats.txt
                echo "  <p>\$(grep \${rv_prim_rev} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
                echo "</details>" >> ${prefix}_primer_contam_stats.txt
            fi

            echo "<br>rv prim first half reverse-complemented \${rv_prim_firsthalf_rev} found \$(grep \${rv_prim_firsthalf_rev} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
            if [ \$(grep \${rv_prim_firsthalf_rev} ${fasta} | wc -l) -ne 0 ]; then
                echo "<details>" >> ${prefix}_primer_contam_stats.txt
                echo "  <summary>ASVs with the first half of the rv primer reversed</summary>" >> ${prefix}_primer_contam_stats.txt
                echo "  <p>\$(grep \${rv_prim_firsthalf_rev} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
                echo "</details>" >> ${prefix}_primer_contam_stats.txt
            fi

            echo "<br>rv prim second half reverse-complemented \${rv_prim_secondhalf_rev} found \$(grep \${rv_prim_secondhalf_rev} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
            if [ \$(grep \${rv_prim_secondhalf_rev} ${fasta} | wc -l) -ne 0 ]; then
                echo "<details>" >> ${prefix}_primer_contam_stats.txt
                echo "  <summary>ASVs with the second half of the rv primer reversed</summary>" >> ${prefix}_primer_contam_stats.txt
                echo "  <p>\$(grep \${rv_prim_secondhalf_rev} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
                echo "</details>" >> ${prefix}_primer_contam_stats.txt
            fi
        done
    else
        rv_prim_firsthalf=\${rv_prim:0:\${#rv_prim}/2}
        rv_prim_secondhalf=\${rv_prim:\${#rv_prim}/2}

        rv_prim_rev=\$(echo \${rv_prim} | rev | sed 's/C/g/g' | sed 's/G/c/g' | sed 's/A/t/g' | sed 's/T/a/g' | sed 's/c/C/g' | sed 's/g/G/g' | sed 's/t/T/g' | sed 's/a/A/g')
        rv_prim_firsthalf_rev=\$(echo \${rv_prim_firsthalf} | rev | sed 's/C/g/g' | sed 's/G/c/g' | sed 's/A/t/g' | sed 's/T/a/g' | sed 's/c/C/g' | sed 's/g/G/g' | sed 's/t/T/g' | sed 's/a/A/g')
        rv_prim_secondhalf_rev=\$(echo \${rv_prim_secondhalf} | rev | sed 's/C/g/g' | sed 's/G/c/g' | sed 's/A/t/g' | sed 's/T/a/g' | sed 's/c/C/g' | sed 's/g/G/g' | sed 's/t/T/g' | sed 's/a/A/g')

        echo "<br>rv prim \${rv_prim} found \$(grep \${rv_prim} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
        if [ \$(grep \${rv_prim} ${fasta} | wc -l) -ne 0 ]; then
            echo "<details>" >> ${prefix}_primer_contam_stats.txt
            echo "  <summary>ASVs with the rv primer</summary>" >> ${prefix}_primer_contam_stats.txt
            echo "  <p>\$(grep \${rv_prim} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
            echo "</details>" >> ${prefix}_primer_contam_stats.txt
        fi

        echo "<br>rv prim first half \${rv_prim_firsthalf} found \$(grep \${rv_prim_firsthalf} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
        if [ \$(grep \${rv_prim_firsthalf} ${fasta} | wc -l) -ne 0 ]; then
            echo "<details>" >> ${prefix}_primer_contam_stats.txt
            echo "  <summary>ASVs with the first half of the rv primer</summary>" >> ${prefix}_primer_contam_stats.txt
            echo "  <p>\$(grep \${rv_prim_firsthalf} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
            echo "</details>" >> ${prefix}_primer_contam_stats.txt
        fi

        echo "<br>rv prim second half \${rv_prim_secondhalf} found \$(grep \${rv_prim_secondhalf} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
        if [ \$(grep \${rv_prim_secondhalf} ${fasta} | wc -l) -ne 0 ]; then
            echo "<details>" >> ${prefix}_primer_contam_stats.txt
            echo "  <summary>ASVs with the second half of the rv primer</summary>" >> ${prefix}_primer_contam_stats.txt
            echo "  <p>\$(grep \${rv_prim_secondhalf} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
            echo "</details>" >> ${prefix}_primer_contam_stats.txt
        fi

        echo "<br>rv prim reverse-complemented \${rv_prim_rev} found \$(grep \${rv_prim_rev} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
        if [ \$(grep \${rv_prim_rev} ${fasta} | wc -l) -ne 0 ]; then
            echo "<details>" >> ${prefix}_primer_contam_stats.txt
            echo "  <summary>ASVs with rv primer reversed</summary>" >> ${prefix}_primer_contam_stats.txt
            echo "  <p>\$(grep \${rv_prim_rev} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
            echo "</details>" >> ${prefix}_primer_contam_stats.txt
        fi

        echo "<br>rv prim first half reverse-complemented \${rv_prim_firsthalf_rev} found \$(grep \${rv_prim_firsthalf_rev} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
        if [ \$(grep \${rv_prim_firsthalf_rev} ${fasta} | wc -l) -ne 0 ]; then
            echo "<details>" >> ${prefix}_primer_contam_stats.txt
            echo "  <summary>ASVs with the first half of the rv primer reversed</summary>" >> ${prefix}_primer_contam_stats.txt
            echo "  <p>\$(grep \${rv_prim_firsthalf_rev} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
            echo "</details>" >> ${prefix}_primer_contam_stats.txt
        fi

        echo "<br>rv prim second half reverse-complemented \${rv_prim_secondhalf_rev} found \$(grep \${rv_prim_secondhalf_rev} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
        if [ \$(grep \${rv_prim_secondhalf_rev} ${fasta} | wc -l) -ne 0 ]; then
            echo "<details>" >> ${prefix}_primer_contam_stats.txt
            echo "  <summary>ASVs with the second half of the rv primer reversed</summary>" >> ${prefix}_primer_contam_stats.txt
            echo "  <p>\$(grep \${rv_prim_secondhalf_rev} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
            echo "</details>" >> ${prefix}_primer_contam_stats.txt
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grep: \$(grep -V | head -n 1 | sed 's/grep (GNU grep) //g')
    END_VERSIONS
    """
}
