# NGSgeno - Next Generation Sequencing Genotyping pipeline

Integrates robotic sample preparation and Illumin Miseq DNA sequencing of amplicons with the APF's Musterer (and later Rodentity) database.

## Dependencies

BBtools https://jgi.doe.gov/data-and-tools/bbtools/
Java https://www.java.com/en/download/manual.jsp
Python 3.8 or newer (3.10+ recommended)
Microsoft Edge is the default browser but this can be changed in ngsrun.bat

Python modules: (see requirements.txt)
openpyxl - for reading/writing excel files
biopython - sequence matching 
requests - connecting to DBs
jsonpickle - save/load experiment info
streamlit - web interface
st_aggrid - interactive web tables
pandas - dataframes to support web tables?

## Changelog:

See changelog.txt

## Robotics operation guide

TBD

## Amplicon sequence analysis description

1. Sequences are copied back from the sequencing instrument (Miseq) to the appropriate folder on the pipeline workstation.
2. The script ngsmatch.py performs the following steps:
2.1. Reads belonging to each well are trimmed for poor quality bases using BBduk (BBtools)
2.2. Paired read ends are merged using BBmerge (BBtools)
2.3. The counts for all unique merged sequences are collected by the script
2.4. We attempt to match each unique sequence against all possible target amplicon sequences, first using exact matching, and then using inexact matching (75% identity) if exact matching fails.
Sequences with a high proportion of counts in any given well, or with at least 50 counts, are reported using their unique sequence. This lets us identify new, real amplicon sequences. Unique sequences that do not make up at least 20% of the reads in a well, or have at least 50 reads, do not get matched and are simply counted as "other" to save time.
2.5. A report named results.csv is generated at this stage with all required metadata to identify where the sequences and their read counts came from.
3. The gt_mice.py script reads the results.csv file and attempts to call Musterer-recognised genotypes (observables) for each tested assay. It writes a report in the standard per-well format named results_gt.csv and a per-mouse genotyping report named mice_gts.csv.

## Definitions of terms

Efficiency - found in the final reports, this is the total number of sequences that contributed to a called genotype divided by the total merged reads counted in the given wells for an assay. It is reported as a percentage.
Allele ratios - found in the final reports, this is the proportion of sequences that matched to a particular haplotype call e.g. wt/mut may give an allele ratio of 0.61/0.39.
Other - reported by ngsmatch.py, this is a collective name for all other unique sequences that were in too small a proportion to meaningfully affect a final genotype call. It is comprised of all unique sequences with less than 20% of the total sequences in a well, or less than 50 reads, whichever is larger. This is to keep the processing time reasonable - roughly 1 hour per 1000 wells vs 12 hours per 1000 wells.

## Credits
The NGS Genotyping Pipeline and NGSXplorer are the creations of The ANU Bioinformatics Consultancy (ABC), The Biomolecular Resource Facility (BRF), The John Curtin School of Medical Research (JCSMR), and The Australian National University (ANU).
The project was initiated at the end of 2018 by the BRF and the construction was undertaken by the ABC (primarily Bob Buckley assisted by Cameron Jack, and later Cameron Jack assisted by Gabi (Gabrielle) Ryan). 
Some additional components were built by the Informatics Team at the Australian Phenomics Facility (led by Philip Wu). Laboratory processes were constructed by the BRF Genotyping team, initially led by Sorelle Bowman, and later by Simone Kuelzer and Peter Milburn.
