process GET_PRIMERFILES {
    label 'process_low'

    container 'ubuntu:20.04'

    input:
    val fw_primer
    val rv_primer

    output:
    path "primers_5end.fa", emit: fasta_5end
    path "primers_3end.fa", emit: fasta_3end

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    touch primers_5end.fa
    touch primers_3end.fa

    if [[ ${fw_primer} == *";"* ]]; then
        fw_primer_arr=(\${${fw_primer}//;/ })
        i=0

        for fw_primer in "\${fw_primer_arr[@]}"
        do
            fw_primer_rvcomp=\$(echo \${fw_primer} | rev | sed 's/C/g/g' | sed 's/G/c/g' | sed 's/T/a/g' | sed 's/A/t/g' | sed 's/g/G/g' | sed 's/c/C/g' | sed 's/a/A/g' | sed 's/t/T/g')

            echo '>fw_primer_\${i};rightmost' >> primers_5end.fa
            echo '\${fw_primer};rightmost' >> primers_5end.fa

            echo '>fw_primer_\${i}_rvcomp' >> primers_3end.fa
            echo \${fw_primer_rvcomp} >> primers_3end.fa
            ((i++))
        done
    else
        fw_primer_rvcomp=\$(echo ${fw_primer} | rev | sed 's/C/g/g' | sed 's/G/c/g' | sed 's/T/a/g' | sed 's/A/t/g' | sed 's/g/G/g' | sed 's/c/C/g' | sed 's/a/A/g' | sed 's/t/T/g')

        echo '>fw_primer;rightmost' >> primers_5end.fa
        echo '${fw_primer};rightmost' >> primers_5end.fa

        echo '>fw_primer_rvcomp' >> primers_3end.fa
        echo \${fw_primer_rvcomp} >> primers_3end.fa
    fi

    if [[ ${rv_primer} == *";"* ]]; then
        rv_primer_arr=(\${${rv_primer}//;/ })
        i=0

        for rv_primer in "\${rv_primer_arr[@]}"
        do
            rv_primer_rvcomp=\$(echo \${rv_primer} | rev | sed 's/C/g/g' | sed 's/G/c/g' | sed 's/T/a/g' | sed 's/A/t/g' | sed 's/g/G/g' | sed 's/c/C/g' | sed 's/a/A/g' | sed 's/t/T/g')

            echo '>rv_primer_\${i};rightmost' >> primers_5end.fa
            echo '\${rv_primer};rightmost' >> primers_5end.fa

            echo '>rv_primer_\${i}_rvcomp' >> primers_3end.fa
            echo \${rv_primer_rvcomp} >> primers_3end.fa
            ((i++))
        done
    else
        rv_primer_rvcomp=\$(echo ${rv_primer} | rev | sed 's/C/g/g' | sed 's/G/c/g' | sed 's/T/a/g' | sed 's/A/t/g' | sed 's/g/G/g' | sed 's/c/C/g' | sed 's/a/A/g' | sed 's/t/T/g')

        echo '>rv_primer;rightmost' >> primers_5end.fa
        echo '${rv_primer};rightmost' >> primers_5end.fa

        echo '>rv_primer_rvcomp' >> primers_3end.fa
        echo \${rv_primer_rvcomp} >> primers_3end.fa
    fi
    """
}
