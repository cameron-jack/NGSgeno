$VER = "v0.28.014"
$DATE = "Date: {0}" -f (Get-Date)

$COMMENT = @"
* CHANGED: primer volumes must add up to 2000nl, primers now at 250nl
* NEW: Warning message about primer volumes
* CHANGED: Auto update pipeline version disabled
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