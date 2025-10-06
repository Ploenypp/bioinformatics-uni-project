# Bioinformatics University Project
This repository contains the work done for course LU3IN402 : Projet Bioformatique during the final year of my double bachlor's in biology and computer science.

_note: I have decided to rewrite, reorganize, and reannotate the code for readability and relevance to my background and current education._

## Project Objective 
to identify biologically significant motifs in DNA sequences

### Biological Context
While open reading frames (ORFs) - motifs defined by a start and a stop codon - may potentially code for proteins, this is not always the case. Transcription factors, proteins that regulate DNA transcription and, by extension, protein expression, often bind to sites close to the coding sequence. Detecting these sites helps confirm the biological relevance of a sequence.

By analyzing ChIP-Seq data (Chromatin Immunopreciptation Sequence), a technique that identifies DNA regions that interact with specific proteins, the project applies computational algorithms to identify frequently occuring, biologically significant candidate motifs. Statistical analyses can then be used to validate these findings.

## General Organization
`datasets` : data to be analyzed in the form of FASTA and JASPAR files
`LU3IN402` : work done for the course, using the JupyterNotebooks and coded provided by the professors

Files outside these folders are my rewritten works. Code taken from the original professor-provided files will be indicated in the files.

*currently working on* : `motif_detection.ipynb`
- generating candidate motifs and detecting their presence in the given data using aforementioend computational algorithms
