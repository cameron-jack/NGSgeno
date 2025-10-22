$VER = "v2.03.002"
$DATE = "Date: {0}" -f (Get-Date)

$COMMENT = @"
* Fixed: non-blocking matching processes correctly closed and cleaned up
* Fixed: incorrect path to style.css
* TODO: User reports that amplicon matching progress was not being shown, grey screen only
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
