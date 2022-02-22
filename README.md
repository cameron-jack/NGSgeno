# NGSgeno - Next Generation Sequencing Genotyping pipeline

Integrates robotic sample preparation and Illumin Miseq DNA sequencing of amplicons with the APF's Musterer (and later Rodentity) database.

## Dependencies

BBtools https://jgi.doe.gov/data-and-tools/bbtools/
Java https://www.java.com/en/download/manual.jsp
Python 3.8 or newer
Microsoft Edge is the default browser but this can be changed in ngsrun.bat

Python modules:
Openpyxl
Biopython
Requests

## To do
* Add column to per-well GT outputs for excess alleles in sequencing
* Combine custom and Musterer pipelines

## Changelog:

Pipeline version: 0.16
* Initial plan was to allow custom workflows to be added to a mouse pipeline run
* Now is a complete reworking into an online "application"
* Replace multiple pipeline pages with a single "Dashboard" interface, with expanding/collapsing sections
* Never leave the main screen - everything is displayed
* Incorporates CSS and Javascript for the first time to allow dynamic browser features
* Refreshes screen on each change
* Select a "run" folder
* Compose an "inventory" of samples and plates
* Plan and track the whole "experiment"
* Keep a "run log" which allows complete re-running of an experiment
* There are no longer any "standard" primer plates - all primer plates are built for purpose, after the Nimbus stage concludes

Pipeline version: 0.15
* Initial attempt to support Musterer mice in custom manifests (running on the mouse pipeline)
* Incomplete due to changes planned for version 0.16
* New error handling/reporting code throughout

Pipeline version: 0.14
* Replaced primercheck.py with validate_assays.py - Only checks against the assay list 
  (not primer plate), and generates list of required primers
* Additional protection aganst unguarded Nimbus output barcodes (creates guarded Echo_COC files)
* Creates guarded copy of custom manifest if it isn't already guarded
* Much improved analysis interface
* Caches moved to run folder

Pipeline version: 0.13
* New analysis (stage 3) interface for running ngsmatch.py and gt_mice.py
* barcode protection for Musterer mnnnm, Rodentity MnnnM, custom cnnnc, plate pnnnp
* Separate out file interfaces to improve code maintainability
* Eliminate excess white-space/empty rows in files as an issue?

Pipeline version: 0.12
* Separated Tm1a,b,c,d into their own cases for genotyping - genotyping *works*
* Removed --fast option from ngsmatch.py
* i7i5 barcodes now applied in rotating fashion (like primer pairs) in echo.py
* Sanity checking incorporated

Pipeline version: 0.11
* Much improved genotyping code
* Hard coded Stage3.bat to run custom pipeline matching
* ngsmatch can now work with custom assay inputs

Pipeline version: 0.10
* Added per-mouse genotyping report (much more human readable)
* Replaced genotyping algorithm to simplify and clean up process - needs another reboot sadly
* Bug fixes in lab pipeline to remove white spaces from names and empty rows from CSV files
* Improved matching performance by 1000x (9x algorithm, 8x multiprocessing, 15x caching)

Pipeline version: 0.9
* Sequence to target matching now matches against all potential target sequences
* Inexact matching has been fixed
* The matching stage is now significantly slower than before
* Sequences with low coverage are reported as "other"
* Per-well "calling" has been removed as unnecessary
* Absence of single outcome assays (e.g. CRE) now results in a NEG genotype rather than a failure to report 

Pipeline version: 0.8
* Extensive UI improvements to clearly delineate mouse runs from custom runs
* Run folders are generated and/or require mouse_ or custom_ to be prefixed - this will help keep the main project clean after being in use a long time
* Multiple taq/water plates are now used if one is insufficient
* Plate dead volumes have been increased to 50 and 700 nl for 384 and 6 well plates respectively
* Numerous code changes to improve modularity

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

