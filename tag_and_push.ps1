$VER = "v2.02.003"
$DATE = "Date: {0}" -f (Get-Date)

$COMMENT = @"
* Fixed: Unknown sequences not being reported in genotyping report (accidentally skipped in v2.02.002)
* TODO: coloured text in PDF report
* TODO: add multiple sequence alignmnent to genotyping report screen
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
