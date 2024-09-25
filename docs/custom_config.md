# MinderooFoundation/OceanOmics-amplicon-nf: assay_params

## Introduction

If you provide the `--assay` parameter, the pipeline will automate choosing certain other parameters in the pipeline. The pipeline will choose these parameters by pulling information from `conf/assay_params.config`.

### Custom parameters

If you would like to use a custom config file, you can provide this with the `-c` option. Your custom file should look like this

```
params {
    assay_params {
        '16SFish' {
            fw_primer = "GACCCTATGGAGCTTTAGAC"
            rv_primer = "CGCTGTTATCCCTADRGTAACT"
        }
        'MiFish' {
            fw_primer = "GTCGGTAAAACTCGTGCCAGC"
            rv_primer = "CATAGTGGGGTATCTAATCCCAGTTTG"
        }
    }
}
```

Note: currently the only parameters that can be automated per assay/read length are `--fw_primer`, and `--rv_primer`. Other parameters can still be added to your custom config file. E.g.,
```
params {
    blast_pid = 97
    skip_lulu = true
}
```
