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
    fw_primer=${fw_primer}
    rv_primer=${rv_primer}
    fw_primer_firsthalf=\${fw_primer:0:\${#fw_primer}/2}
    fw_primer_secondhalf=\${fw_primer:\${#fw_primer}/2}
    rv_primer_firsthalf=\${rv_primer:0:\${#rv_primer}/2}
    rv_primer_secondhalf=\${rv_primer:\${#rv_primer}/2}

    fw_primer_rev=\$(echo \${fw_primer} | rev)
    rv_primer_rev=\$(echo \${rv_primer} | rev)
    fw_primer_firsthalf_rev=\$(echo \${fw_primer_firsthalf} | rev)
    fw_primer_secondhalf_rev=\$(echo \${fw_primer_secondhalf} | rev)
    rv_primer_firsthalf_rev=\$(echo \${rv_primer_firsthalf} | rev)
    rv_primer_secondhalf_rev=\$(echo \${rv_primer_secondhalf} | rev)

    touch ${prefix}_primer_contam_stats.txt
    echo "fw primer \${fw_primer} found \$(grep \${fw_primer} ${fasta} | wc -l) times" > ${prefix}_primer_contam_stats.txt
    if [ \$(grep \${fw_primer} ${fasta} | wc -l) -ne 0 ]; then
        echo "<details>" >> ${prefix}_primer_contam_stats.txt
        echo "  <summary>ASVs with the fw primer</summary>" >> ${prefix}_primer_contam_stats.txt
        echo "  <p>\$(grep \${fw_primer} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
        echo "</details>" >> ${prefix}_primer_contam_stats.txt
    fi

    echo "<br>rv primer \${rv_primer} found \$(grep \${rv_primer} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
    if [ \$(grep \${rv_primer} ${fasta} | wc -l) -ne 0 ]; then
        echo "<details>" >> ${prefix}_primer_contam_stats.txt
        echo "  <summary>ASVs with the rv primer</summary>" >> ${prefix}_primer_contam_stats.txt
        echo "  <p>\$(grep \${rv_primer} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
        echo "</details>" >> ${prefix}_primer_contam_stats.txt
    fi

    echo "<br>fw primer first half \${fw_primer_firsthalf} found \$(grep \${fw_primer_firsthalf} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
    if [ \$(grep \${fw_primer_firsthalf} ${fasta} | wc -l) -ne 0 ]; then
        echo "<details>" >> ${prefix}_primer_contam_stats.txt
        echo "  <summary>ASVs with the first half of the fw primer</summary>" >> ${prefix}_primer_contam_stats.txt
        echo "  <p>\$(grep \${fw_primer_firsthalf} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
        echo "</details>" >> ${prefix}_primer_contam_stats.txt
    fi

    echo "<br>fw primer second half \${fw_primer_secondhalf} found \$(grep \${fw_primer_secondhalf} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
    if [ \$(grep \${fw_primer_secondhalf} ${fasta} | wc -l) -ne 0 ]; then
        echo "<details>" >> ${prefix}_primer_contam_stats.txt
        echo "  <summary>ASVs with the second half of the fw primer</summary>" >> ${prefix}_primer_contam_stats.txt
        echo "  <p>\$(grep \${fw_primer_secondhalf} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
        echo "</details>" >> ${prefix}_primer_contam_stats.txt
    fi

    echo "<br>rv primer first half \${rv_primer_firsthalf} found \$(grep \${rv_primer_firsthalf} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
    if [ \$(grep \${rv_primer_firsthalf} ${fasta} | wc -l) -ne 0 ]; then
        echo "<details>" >> ${prefix}_primer_contam_stats.txt
        echo "  <summary>ASVs with the first half of the rv primer</summary>" >> ${prefix}_primer_contam_stats.txt
        echo "  <p>\$(grep \${rv_primer_firsthalf} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
        echo "</details>" >> ${prefix}_primer_contam_stats.txt
    fi

    echo "<br>rv primer second half \${rv_primer_secondhalf} found \$(grep \${rv_primer_secondhalf} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
    if [ \$(grep \${rv_primer_secondhalf} ${fasta} | wc -l) -ne 0 ]; then
        echo "<details>" >> ${prefix}_primer_contam_stats.txt
        echo "  <summary>ASVs with the second half of the rv primer</summary>" >> ${prefix}_primer_contam_stats.txt
        echo "  <p>\$(grep \${rv_primer_secondhalf} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
        echo "</details>" >> ${prefix}_primer_contam_stats.txt
    fi

    echo "<br>fw primer reversed \${fw_primer_rev} found \$(grep \${fw_primer_rev} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
    if [ \$(grep \${fw_primer_rev} ${fasta} | wc -l) -ne 0 ]; then
        echo "<details>" >> ${prefix}_primer_contam_stats.txt
        echo "  <summary>ASVs with fw primer reversed</summary>" >> ${prefix}_primer_contam_stats.txt
        echo "  <p>\$(grep \${fw_primer_rev} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
        echo "</details>" >> ${prefix}_primer_contam_stats.txt
    fi

    echo "<br>rv primer reversed \${rv_primer_rev} found \$(grep \${rv_primer_rev} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
    if [ \$(grep \${rv_primer_rev} ${fasta} | wc -l) -ne 0 ]; then
        echo "<details>" >> ${prefix}_primer_contam_stats.txt
        echo "  <summary>ASVs with rv primer reversed</summary>" >> ${prefix}_primer_contam_stats.txt
        echo "  <p>\$(grep \${rv_primer_rev} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
        echo "</details>" >> ${prefix}_primer_contam_stats.txt
    fi

    echo "<br>fw primer first half reversed \${fw_primer_firsthalf_rev} found \$(grep \${fw_primer_firsthalf_rev} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
    if [ \$(grep \${fw_primer_firsthalf_rev} ${fasta} | wc -l) -ne 0 ]; then
        echo "<details>" >> ${prefix}_primer_contam_stats.txt
        echo "  <summary>ASVs with the first half of the fw primer reversed</summary>" >> ${prefix}_primer_contam_stats.txt
        echo "  <p>\$(grep \${fw_primer_firsthalf_rev} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
        echo "</details>" >> ${prefix}_primer_contam_stats.txt
    fi

    echo "<br>fw primer second half reversed \${fw_primer_secondhalf_rev} found \$(grep \${fw_primer_secondhalf_rev} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
    if [ \$(grep \${fw_primer_secondhalf_rev} ${fasta} | wc -l) -ne 0 ]; then
        echo "<details>" >> ${prefix}_primer_contam_stats.txt
        echo "  <summary>ASVs with the second half of the fw primer reversed</summary>" >> ${prefix}_primer_contam_stats.txt
        echo "  <p>\$(grep \${fw_primer_secondhalf_rev} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
        echo "</details>" >> ${prefix}_primer_contam_stats.txt
    fi

    echo "<br>rv primer first half reversed \${rv_primer_firsthalf_rev} found \$(grep \${rv_primer_firsthalf_rev} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
    if [ \$(grep \${rv_primer_firsthalf_rev} ${fasta} | wc -l) -ne 0 ]; then
        echo "<details>" >> ${prefix}_primer_contam_stats.txt
        echo "  <summary>ASVs with the first half of the rv primer reversed</summary>" >> ${prefix}_primer_contam_stats.txt
        echo "  <p>\$(grep \${rv_primer_firsthalf_rev} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
        echo "</details>" >> ${prefix}_primer_contam_stats.txt
    fi

    echo "<br>rv primer second half reversed \${rv_primer_secondhalf_rev} found \$(grep \${rv_primer_secondhalf_rev} ${fasta} | wc -l) times" >> ${prefix}_primer_contam_stats.txt
    if [ \$(grep \${rv_primer_secondhalf_rev} ${fasta} | wc -l) -ne 0 ]; then
        echo "<details>" >> ${prefix}_primer_contam_stats.txt
        echo "  <summary>ASVs with the second half of the rv primer reversed</summary>" >> ${prefix}_primer_contam_stats.txt
        echo "  <p>\$(grep \${rv_primer_secondhalf_rev} ${fasta} -B 1 | grep -v "\\-\\-" | sed 's/>//g') </p>" >> ${prefix}_primer_contam_stats.txt
        echo "</details>" >> ${prefix}_primer_contam_stats.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grep: \$(grep -V | head -n 1 | sed 's/grep (GNU grep) //g')
    END_VERSIONS
    """
}
