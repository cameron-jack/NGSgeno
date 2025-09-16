$VER = "v2.02.000"
$DATE = "Date: {0}" -f (Get-Date)

$COMMENT = @"
* Fixed: retained parameter in match calling function caused crashed subprocesses
* Fixed: amplicon matching was looking for primers
* Fixed: debugging parameters weren't guarded against being absent
* Fixed: amplicon inexact matching was using incorrectly named variable
* Fixed: mixed types in view columns causing pyarrow error
* Fixed: genotyping error messages for unknown sequences
* Fixed: borked call to generate primer information
* New: unreported/reported amplicon lists
* New: tables being migrated to new Streamlit 1.49.0 style tables
* Changed: identity limits now restricted to 0.0-1.0 with warnings
* Changed: Streamlit updated to 1.49.1
* Changed: streamlit.aggrid updated 1.1.8
* Changed: Now reporting amplicons in PDF rather than rich text with fpdf2 v2.8.4
* Changed: Report page removes expanding sections for clarity
* TODO: coloured text in PDF report
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
