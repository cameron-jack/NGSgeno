python ..\bin\ngsmatch.py %1
python ..\bin\makehtml.py Results.csv
python ..\bin\genotype.py -g ..\library\groups_SB_20201111.csv -o Results_gt.csv Results.csv 
python ..\bin\sanity.py -o Results_sanity.csv Results_gt.csv
