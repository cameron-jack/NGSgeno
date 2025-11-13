$VER = "v2.03.008"
$DATE = "Date: {0}" -f (Get-Date)

$COMMENT = @"
* Fixed: forgot to update version variable in ngsgeno.py
* New: updated README to NGSgeno2 and added acknowledgement of BRF genotyping team members
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
