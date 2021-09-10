# Pipeline version 0.10

python ..\bin\ngsmatch.py %1
python ..\bin\makehtml.py Results.csv
python ..\bin\gt_mice.py -r ..\library\reference_sequences_20201111.txt -c ..\library\NGS_assay_conversions_20201125.xlsx -k ..\bin\config.ini -o results_gt.csv Results.csv

