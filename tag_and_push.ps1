$VER = "v0.28.009"
$DATE = "Date: {0}" -f (Get-Date)

$COMMENT = @"
* FIXED: Primer display now correctly reports primer use
* CHANGED: primers without any uses and without any volume are not reported
* NEW: Primer display now flags red any primer that without sufficient availability
* CHANGED: warnings are no longer reported to console or debugging by default
* BUG: Autosaving feature is not working
* TODO: Put autosaving on file upload and generate functions
* TODO: require the assaylist file for PCR1 output file tracking
* TODO: Save buttons to CSV and Excel for primer display
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
