$VER = "v2.04.000"
$DATE = "Date: {0}" -f (Get-Date)

$COMMENT = @"
* New: Help launched in separate tab, with more detailed instructions and screenshots
* New: HTML documentation and PDF generation with a doc folder
* New: NGSG Retype and NGSG Reference now added to core repository
* New: pydocgen.py script to generate HTML documentation from docstrings (API documentation)
* Changed: NGSG Version Select now runs on port 9222
* Changed: NGSG Reference now runs on port 9224
* Changed: NGSG Retype now runs on port 9225
* Changed: NGSG Retype now displays a message in the right-hand panel if no assays are found to need retyping
* Changed: all code is now in the src folder, instead of bin, with a single entry point in the root folder
* Changed: each pipeline stage is now in its own file in the src folder
* Changed: add_css function moved from display_components to stutil
* Fixed: simplified jsonpickle experiment loading
* Fixed: typo in code that collates already reported wells at the Report stage
* NOTE: Do not try to use multithreading in Streamlit, you cannot hold thread handles in event driven code
* TODO: coloured text in PDF report
* TODO: changing included variants in the alignment table does not redo the alignment
"@

Move-Item -Path "changelog.txt" -Destination "changelog_old.txt"
Set-Content -Path "changelog.txt" -Value $VER
Add-Content -Path "changelog.txt" -Value $DATE
Add-Content -Path "changelog.txt" -Value $COMMENT
Add-Content -Path "changelog.txt" -Value ""
Get-Content -Path "changelog_old.txt" | Add-Content -Path "changelog.txt"
Remove-Item -Path "changelog_old.txt"

git add -u
$MSG = $COMMENT
git commit -m $MSG
git tag -a $VER -m $MSG

git push
git push origin $VER
