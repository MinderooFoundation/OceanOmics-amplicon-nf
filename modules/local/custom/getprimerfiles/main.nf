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

    fw_prim="${fw_primer}"
    rv_prim="${rv_primer}"

    if [[ \${fw_prim} == *";"* ]]; then
        fw_primer_arr=(\${fw_prim//;/ })
        fw_arraylength=\${#fw_primer_arr[@]}

        for (( i=0; i<\${fw_arraylength}; i++ ));
        do
            fw_primer_rvcomp=\$(echo \${fw_primer_arr[\$i]} | rev | sed 's/C/g/g' | sed 's/G/c/g' | sed 's/T/a/g' | sed 's/A/t/g' | sed 's/g/G/g' | sed 's/c/C/g' | sed 's/a/A/g' | sed 's/t/T/g' | sed 's/R/y/g' | sed 's/Y/r/g' | sed 's/K/m/g' | sed 's/M/k/g' | sed 's/B/v/g' | sed 's/D/h/g' | sed 's/H/d/g' | sed 's/V/b/g' | sed 's/y/Y/g' | sed 's/r/R/g' | sed 's/m/M/g' | sed 's/k/K/g' | sed 's/v/V/g' | sed 's/h/H/g' | sed 's/d/D/g' | sed 's/b/B/g')

            echo '>fw_primer_'\${i}';rightmost' >> primers_5end.fa
            echo \${fw_primer_arr[\$i]}';rightmost' >> primers_5end.fa

            echo '>fw_primer_'\${i}'_rvcomp' >> primers_3end.fa
            echo \${fw_primer_rvcomp} >> primers_3end.fa
        done
    else
        fw_primer_rvcomp=\$(echo ${fw_primer} | rev | sed 's/C/g/g' | sed 's/G/c/g' | sed 's/T/a/g' | sed 's/A/t/g' | sed 's/g/G/g' | sed 's/c/C/g' | sed 's/a/A/g' | sed 's/t/T/g' | sed 's/R/y/g' | sed 's/Y/r/g' | sed 's/K/m/g' | sed 's/M/k/g' | sed 's/B/v/g' | sed 's/D/h/g' | sed 's/H/d/g' | sed 's/V/b/g' | sed 's/y/Y/g' | sed 's/r/R/g' | sed 's/m/M/g' | sed 's/k/K/g' | sed 's/v/V/g' | sed 's/h/H/g' | sed 's/d/D/g' | sed 's/b/B/g')

        echo '>fw_primer;rightmost' >> primers_5end.fa
        echo '${fw_primer};rightmost' >> primers_5end.fa

        echo '>fw_primer_rvcomp' >> primers_3end.fa
        echo \${fw_primer_rvcomp} >> primers_3end.fa
    fi

    if [[ \${rv_prim} == *";"* ]]; then
        rv_primer_arr=(\${rv_prim//;/ })
        rv_arraylength=\${#rv_primer_arr[@]}

        for (( i=0; i<\${rv_arraylength}; i++ ));
        do
            rv_primer_rvcomp=\$(echo \${rv_primer_arr[\$i]} | rev | sed 's/C/g/g' | sed 's/G/c/g' | sed 's/T/a/g' | sed 's/A/t/g' | sed 's/g/G/g' | sed 's/c/C/g' | sed 's/a/A/g' | sed 's/t/T/g' | sed 's/R/y/g' | sed 's/Y/r/g' | sed 's/K/m/g' | sed 's/M/k/g' | sed 's/B/v/g' | sed 's/D/h/g' | sed 's/H/d/g' | sed 's/V/b/g' | sed 's/y/Y/g' | sed 's/r/R/g' | sed 's/m/M/g' | sed 's/k/K/g' | sed 's/v/V/g' | sed 's/h/H/g' | sed 's/d/D/g' | sed 's/b/B/g')

            echo '>rv_primer_'\${i}';rightmost' >> primers_5end.fa
            echo \${rv_primer_arr[\$i]}';rightmost' >> primers_5end.fa

            echo '>rv_primer_'\${i}'_rvcomp' >> primers_3end.fa
            echo \${rv_primer_rvcomp} >> primers_3end.fa
        done
    else
        rv_primer_rvcomp=\$(echo ${rv_primer} | rev | sed 's/C/g/g' | sed 's/G/c/g' | sed 's/T/a/g' | sed 's/A/t/g' | sed 's/g/G/g' | sed 's/c/C/g' | sed 's/a/A/g' | sed 's/t/T/g' | sed 's/R/y/g' | sed 's/Y/r/g' | sed 's/K/m/g' | sed 's/M/k/g' | sed 's/B/v/g' | sed 's/D/h/g' | sed 's/H/d/g' | sed 's/V/b/g' | sed 's/y/Y/g' | sed 's/r/R/g' | sed 's/m/M/g' | sed 's/k/K/g' | sed 's/v/V/g' | sed 's/h/H/g' | sed 's/d/D/g' | sed 's/b/B/g')

        echo '>rv_primer;rightmost' >> primers_5end.fa
        echo '${rv_primer};rightmost' >> primers_5end.fa

        echo '>rv_primer_rvcomp' >> primers_3end.fa
        echo \${rv_primer_rvcomp} >> primers_3end.fa
    fi
    """
}
