# NGSgeno - Next Generation Sequencing Genotyping pipeline

Pipeline version: 0.8
version_comments:
* Extensive UI improvements to clearly delineate mouse runs from custom runs
* Run folders are generated and/or require mouse_ or custom_ to be prefixed - this will help keep the main project clean after being in use a long time
* Multiple taq/water plates are now used if one is insufficient
* Plate dead volumes have been increased to 50 and 700 nl for 384 and 6 well plates respectively
* Numerous code changes to improve modularity

Integrates robotic sample preparation and Illumin Miseq DNA sequencing of amplicons with the APF's Musterer (and later Rodentity) database.
