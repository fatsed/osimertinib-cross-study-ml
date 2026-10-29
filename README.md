
# Confounding-Aware Cross-Study Evaluation of Osimertinib Response States

## About the Project

Osimertinib is widely used to treat EGFR-mutant non-small cell lung cancer (NSCLC), but some cancer cells can survive treatment and eventually develop resistance. Understanding the difference between drug-tolerant persister (DTP) cells and cells with acquired resistance remains an important challenge.

This project explores whether gene expression data can be used to distinguish these two states across independent studies.

One of the main concerns is that datasets often differ in their cell lines, treatment conditions, and experimental procedures. As a result, a machine learning model might learn these differences instead of the biological patterns related to drug response.

The goal is to investigate this problem and evaluate how reliably machine learning models perform when tested on data from different studies.

## Research Question

Can machine learning distinguish osimertinib-tolerant cells from acquired-resistant cells across independent studies without relying on study-specific differences?

## Planned Approach

The project will involve:

- Reviewing publicly available transcriptomic datasets and selecting suitable samples.
- Preparing and comparing gene expression data from different studies.
- Training baseline machine learning models.
- Evaluating model performance on independent datasets.
- Investigating how study origin, cell lines, and treatment conditions affect predictions.
- Examining whether important gene expression features remain consistent across studies.

## Data

The study will use publicly available transcriptomic datasets, primarily from the [NCBI Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/).

Dataset selection and preprocessing will be documented as the project progresses.

## Current Progress

The project is currently in its initial development stage. The research question and relevant literature have been reviewed, and the next step is to examine the available datasets before starting model development.

Code, analysis results, and figures will be added as the experiments are carried out.
