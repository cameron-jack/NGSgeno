rem Pipeline version 0.10
rem %1 name of Stage3.csv file
rem %2 path to reference sequence file (.txt or .csv)
rem %3 path to NGS assay conversions file (.xlsx)

rem python ..\bin\ngsmatch.py %1
del ..\library\mouse_match_cache.jsn
del ..\library\mouse_miss_cache.jsn
python ..\bin\ngsmatch.py -t ..\library\reference_sequences_20211006.txt -o Results.csv -c ..\library\mouse_match_cache.jsn -m ..\library\mouse_miss_cache.jsn -k 1 --custom Stage3.csv
rem python ..\bin\makehtml.py Results.csv
rem python ..\bin\gt_mice.py -r ..\library\reference_sequences_20201111.txt -c ..\library\NGS_assay_conversions_20201125.xlsx -k ..\bin\config.ini -o results_gt.csv Results.csv
rem python ..\bin\gt_mice.py -r %2 -c %3 -k ..\bin\config.ini -o results_gt.csv Results.csv
