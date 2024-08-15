$VER = "v0.28.010"
$DATE = "Date: {0}" -f (Get-Date)

$COMMENT = @"
* FIXED: autosave now works and runs after parsing uploads and after generating files
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
