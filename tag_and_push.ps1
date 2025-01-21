$VER = "v1.01.00"
$DATE = "Date: {0}" -f (Get-Date)

$COMMENT = @"
* NEW: Illumina FASTQ sequence count added (takes about 4 seconds for a full run)
* FIXED: set/list index issue when checking for new Illumina FASTQ file names
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
