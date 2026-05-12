REM Don't auto update 
REM cmd /k "venv_ngsgeno\Scripts\activate & git pull & streamlit run --server.enableXsrfProtection false src\ngsgeno.py"
cmd /k "venv_ngsgeno\Scripts\activate & streamlit run --server.enableXsrfProtection false src\ngsgeno.py"