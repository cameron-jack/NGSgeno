$VER = "v0.28.016"
$DATE = "Date: {0}" -f (Get-Date)

$COMMENT = @"
* Final preproduction fixes
* FIXED: Typo in efficiency column in the results.csv from the allele calling 
* CHANGED: interface to show FASTQ files already present on the Allele Calling page
* REMOVED: FASTQ upload
* NEW: auto saving the primer report file to primer_list.csv
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
