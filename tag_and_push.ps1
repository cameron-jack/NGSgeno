$VER = "v0.28.011"
$DATE = "Date: {0}" -f (Get-Date)

$COMMENT = @"
* NEW: Primer display on Primer tab 1 now has a download button to save table as CSV
* CHANGED: Many unnecessary screen messages now silent
* FIXED: ngsmatch.py now exits softly when no samples are present
* FIXED: table colouring reverts correctly for primers
* NEW: Primer table now displays full well and dead volumes
* TODO: require the assaylist file for PCR1 output file tracking
* TODO: Echo COC files should be invalidated by changed Nimbus files
* TODO: Excel format save button for primer list 
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
