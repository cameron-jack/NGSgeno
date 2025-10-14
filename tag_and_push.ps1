$VER = "v2.03.000"
$DATE = "Date: {0}" -f (Get-Date)

$COMMENT = @"
* Fixed: obsolete code in FASTA parser only observed when files had errors
* New: multiple sequence alignment added to genotyping report PDF
* Changed: default font changed to Inconsolata monospace for better alignment viewing
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
