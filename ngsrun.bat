@echo Run the NGS Geno Pipeline v1.6

# set NGS="%HOMEDRIVE%%HOMEPATH%\NGSgeno"
#set port=9123

# cd %NGS%

# start /d %NGS% python -m http.server --cgi %port% 

#MicrosoftEdge http://localhost:9123/cgi-bin/cgi-nimbus.py

streamlit run ngsgeno_gr.py