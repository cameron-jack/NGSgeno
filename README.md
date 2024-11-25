# NGSgeno - Next Generation Sequencing Genotyping pipeline

Integrates robotic sample preparation and Illumina Miseq DNA sequencing of amplicons with the APF's Rodentity database.

# To do
* Documentation! Sphinx docs of code. 
* Modular indexing code
* New file generation/transaction code
* Secure access through SSL
* Multi-user logins and user tracking (attach user to logs and user to uploads and transacted steps)
* Access levels -> Superuser (can do everything), Admin (can perform various configurable actions), User - just run the pipeline and not change anything
* Integrate inputs and outputs with automatic uploads to Lab Archives

## Dependencies

BBtools https://jgi.doe.gov/data-and-tools/bbtools/
Java https://www.java.com/en/download/manual.jsp
Python 3.8 or newer (3.10+ recommended)
Microsoft Edge is the default browser but this can be changed in ngsrun.bat

Python modules: (see requirements.txt)
openpyxl - for reading/writing excel files
biopython==1.81 - sequence matching 
requests - connecting to DBs
jsonpickle - save/load experiment info
streamlit==1.26 - web interface
st_aggrid - interactive web tables
extra_streamlit_components - advanced GUI widgets
pandas - dataframes to support web tables

## Changelog:

See changelog.txt

## Robotics operation guide

TBD

## Credits
The NGS Genotyping Pipeline and NGSXplorer are the creations of The ANU Bioinformatics Consultancy (ABC), 
The Biomolecular Resource Facility (BRF), The John Curtin School of Medical Research (JCSMR), and The Australian National University (ANU).
The project was initiated at the end of 2018 by the BRF and the construction was undertaken by the ABC - primarily 
Bob Buckley assisted by Cameron Jack, and later Cameron Jack assisted by Gabi (Gabrielle) Ryan. 
Some additional components were built by the Informatics Team at the Australian Phenomics Facility (led by Philip Wu). 
Laboratory processes were constructed by the BRF Genotyping team, initially led by Sorelle Bowman, and later by Simone Kuelzer and Peter Milburn.

## License
This software (product) is issued with the MIT license

Copyright 2018 The Australian National University

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated 
documentation files (the “Software”), to deal in the Software without restriction, including without 
limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE 
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS 
OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT 
OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.