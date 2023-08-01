process COUNTOFSEQ {
    tag "$reads"
    label 'process_medium'

    container 'ubuntu:20.04'

    input:
    val prefix
    path reads

    output:
    path "*.txt"         , emit: counts
    path "versions.yml"  , emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    """
    if echo "${reads}" | grep -q "\\.gz\$"; then
        zcat ${reads} | awk -F'[>;]' '/^>/ {if (sample) { count[sample]++ }; sample = substr(\$0, 2); sub(/;.*/, "", sample)} !/^>/ {sequence[sample] = sequence[sample] \$0} END {count[sample]++; for (sample in count) { printf "%s\\n%d\\n", sample, count[sample] }}' > ${prefix}_counts.txt
    else
        cat ${reads} | awk -F'[>;]' '/^>/ {if (sample) { count[sample]++ }; sample = substr(\$0, 2); sub(/;.*/, "", sample)} !/^>/ {sequence[sample] = sequence[sample] \$0} END {count[sample]++; for (sample in count) { printf "%s\\n%d\\n", sample, count[sample] }}' > ${prefix}_counts.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: 1.3.4
    END_VERSIONS
    """
}