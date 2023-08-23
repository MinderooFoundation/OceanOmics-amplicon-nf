#!/usr/bin/env python3

import joblib
import numpy as np
from Bio import SeqIO
import sys
from argparse import ArgumentParser
import logging

CUTOFF = 0.9
TPROCESS = "OCOM_NBC"

def main():
    parser = ArgumentParser(
        prog = 'Fish Naive Bayes Classifier',
        description = 'Classifies species and genus for a given fasta file of 12S, 16S, or CO1 genes.',
        epilog = 'Prediction probabilities range from 0 to 1, the default cutoff is 0.9.\nModels were trained using scikit-learn v1.2.1.')

    parser.add_argument('-i', '--input',
            help='Name of the input fasta file', required=True)
    parser.add_argument('-o', '--output',
            help='Name of the output file - TSV format', required=True)
    parser.add_argument('-t', '--task_process',
            help='Task Process, if using in Nextflow module', required=False,
            default=TPROCESS, type=str)
    parser.add_argument('-c', '--cutoff',
            help=f'Score cutoff that sets predictions to NA. Default is {CUTOFF}.',
            default=CUTOFF, type=float)
    parser.add_argument('-v', help='Verbose mode (quiet by default)',
            action='store_true')

    args = parser.parse_args()

    to_parse = args.input
    output = args.output
    task_process = args.task_process
    cutoff = args.cutoff

    py_version = sys.version
    py_version = py_version.split(" ")
    py_version = py_version[0]

    verbose = args.v
    if verbose:
        logging.basicConfig(format='%(asctime)s: %(message)s', level=logging.INFO,
                datefmt="%Y-%m-%d %H:%M:%S")

    logging.info(f'Loading FASTA file "{to_parse}".')
    # load the entire file into memory so we can give all sequences to the classifiers at once
    seqs = [[x.description, str(x.seq)] for x in SeqIO.parse(to_parse, 'fasta')]
    y, X = list(), list()
    for a in seqs:
        y.append(a[0])
        X.append(a[1])
    logging.info('Finished loading.')

    logging.info('Deciding on the class of genes you have.')
    geneclass_pipeline = joblib.load('/mnt/scratch/C01_12S_16S_pipeline.joblib')
    geneclasses = geneclass_pipeline.classes_

    gene_dict = {x:[] for x in geneclasses}
    gene_dict['NA'] = [] # those genes where we're not sure what they are!

    gene_prob = geneclass_pipeline.predict_proba(X)
    for gene_probas, input_gene, label in zip(gene_prob, X, y):
        gene_prediction_label = geneclasses[np.argmax(gene_probas)]
        best_score = max(gene_probas)
        best_label = gene_prediction_label
        if best_score < 0.9:
            best_label = 'NA' # we don't know what this gene is
        gene_dict[best_label].append( (input_gene, label) )

    logging.info(f'We have the following gene classes: { {x:len(gene_dict[x]) for x in gene_dict.keys()} }')
    logging.info(f'Proceeding with species/genus prediction.')

    t12s_genus = joblib.load('/mnt/scratch/12S_genus_pipeline.joblib')
    t12s_species = joblib.load('/mnt/scratch/12S_species_pipeline.joblib')

    t16s_genus = joblib.load('/mnt/scratch/16S_genus_pipeline.joblib')
    t16s_species = joblib.load('/mnt/scratch/16S_species_pipeline.joblib')

    CO1_genus = joblib.load('/mnt/scratch/CO1_genus_pipeline.joblib')
    CO1_species = joblib.load('/mnt/scratch/CO1_species_pipeline.joblib')

    classifier_dict = {'12S': (t12s_genus, t12s_species), '16S': (t16s_genus, t16s_species),\
            'CO1': (CO1_genus, CO1_species)}

    logging.info('Finished loading models.')

    logging.info(f'Now predicting. Prediction cutoffs are set to {cutoff}.')
    with open(output, 'w') as out:
        out.write('Gene\tLabel\tGenus prediction\tGenus score\tSpecies prediction\tSpecies score\n')
        for gene in gene_dict: # gene: 12S, 16S, CO1
            if len(gene_dict[gene]) == 0: continue
            to_classify, labels = zip(*gene_dict[gene])
            if gene == 'NA':
                # these are the genes with no class (12S/16S/CO1)
                for label in labels:
                    out.write(f'{gene}\t{label}\tNA\tNA\tNA\tNA\n')
                continue

            this_models = classifier_dict[gene]
            genus_model, species_model = this_models
            # ('ACCGCGGTTATACGAGAGGCCCAAGTTGATAGACTCCGGCGTAAAGAGTGGTTAAGATAAATTTTAAACTAAAGCCGAACGCCCTCAAAGCTGTTATACGCTCCCGAGGGTAAGAAGCCCAATCACGAAAGTGGCTTTATACCAGCTGAACCCACGAAAGCTATGACA', '1749_1916_KC136482.1 334880 Oplegnathus woodwardi KC136482.1 Oplegnathus woodwardi isolate NZ_001 collection-date 29-Nov-2005 cytochrome b (CYTB) gene, partial cds; tRNA-Thr and tRNA-Pro genes, D-loop, tRNA-Phe, 12S ribosomal RNA, and tRNA-Val genes, complete sequence; and 16S ribosomal RNA gene, partial sequence; mitochondrial') 12S
            genus_preds = genus_model.predict_proba(to_classify)
            species_preds = species_model.predict_proba(to_classify)
            for label, genus_pred, species_pred in zip(labels, genus_preds, species_preds):
                max_genus = max(genus_pred) # that's a score
                if max_genus < cutoff:
                    genus_pred_label = 'NA'
                else:
                    genus_pred_label = genus_model.classes_[np.argmax(genus_pred)]

                max_species = max(species_pred)
                if max_species < cutoff:
                    species_pred_label = 'NA'
                else:
                    species_pred_label = species_model.classes_[np.argmax(species_pred)]
                out.write(f'{gene}\t{label}\t{genus_pred_label}\t{max_genus:.2f}\t{species_pred_label}\t{max_species:.2f}\n')

    logging.info(f'Done! Wrote results to "{output}".')

    # Create version .yml file
    with open('versions.yml', 'w') as version_file:
        version_file.write(f"\"{task_process}\":\n")
        version_file.write(f"    python: {py_version}\n")
        version_file.write(f"    joblib: {joblib.__version__}\n")
        version_file.write(f"    numpy: {np.__version__}\n")

main()
