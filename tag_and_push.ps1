$VER = "v1.00.000"
$DATE = "Date: {0}" -f (Get-Date)

$COMMENT = @"
Production release
* NEW: supports new Illumina i100 sequencer files
* NEW: adds a Force option for generating primer picklists in spite of insufficient primer volmes
* NEW: proper messaging support extended through primer file generation functions
* CHANGED: push returned to tag_and_push.ps1
* FIXED: inccorect indentation in ngsmatch.py giving false lock file errors
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
