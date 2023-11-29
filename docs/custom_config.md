# MinderooFoundation/OceanOmics-amplicon-nf: assay_params

## Introduction

If you provide the `--assay_readlength` parameter, the pipeline will automate choosing certain other parameters in the pipeline. The pipeline will choose these parameters by pulling information from `conf/assay_params.config`.

### Custom parameters

If you would like to use a custom config file, you can provide this with the `-c` option. Your custom file should look like this

```
params {
    assay_params {
        '12S_150' {
            seqtk_trim       = false
            seqtk_length     = null
            asv_min_overlap  = 12
            lca_qcov         = 100
        }
        '12S_300' {
            seqtk_trim       = true
            seqtk_length     = 150
            asv_min_overlap  = 12
            lca_qcov         = 100
        }
    }
}
```

Note: currently the only parameters that can be automated per assay/read length are `--seqtk_trim`, `--seqtk_length`, `--asv_min_overlap`, and `--lca_qcov`. Other parameters can still be added to your custom config file. E.g.,
```
params {
    blast_pid = 97
    skip_lulu = true
}
```
