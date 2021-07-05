@echo Run the NGS Geno Pipeline

set NGS="%HOMEDRIVE%%HOMEPATH%\NGSgeno"
set port=9000

cd %NGS%

start /d %NGS% python -m http.server --cgi %port% 

MicrosoftEdge http://localhost:9000/cgi-bin/cgi-nimbus1.py
