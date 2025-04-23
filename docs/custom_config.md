# MinderooFoundation/OceanOmics-amplicon-nf: assay_params

## Introduction

If you provide the `--assay` parameter, the pipeline will automate choosing certain other parameters in the pipeline. The pipeline will choose these parameters by pulling information from `conf/assay_params.config`.

### Custom parameters

If you would like to use a custom config file, you can provide this with the `-c` option. This is an example of a custom config file.

```
params {
    assay_params {
        '16SFish' {
            fw_primer = "GACCCTATGGAGCTTTAGAC"
            rv_primer = "CGCTGTTATCCCTADRGTAACT"
            asv_min_overlap = 24
            read_maxlen     = 1000
            seqtk_trim      = false
            seqtk_length    = 1000
        }
        'MiFish' {
            fw_primer = "GTCGGTAAAACTCGTGCCAGC"
            rv_primer = "CATAGTGGGGTATCTAATCCCAGTTTG"
            asv_min_overlap = 12
            read_maxlen     = 160
            seqtk_trim      = true
            seqtk_length    = 120
        }
    }
}
```

Note: currently the only parameters that can be automated per assay are `--fw_primer`, `--rv_primer`, `asv_min_overlap`, `read_maxlen`, `seqtk_trim`, and `seqtk_length`. Other parameters can still be added to your custom config file. E.g.,
```
params {
    blast_pid = 97
    skip_lulu = true
}
```
