$VER = "v0.28.018"
$DATE = "Date: {0}" -f (Get-Date)

$COMMENT = @"
* NEW: Added 384-well DNA-plate upload for reruns of failed assays
* NEW: sample viewer display includes new 384-well plates
* NEW: 384-plates without matching Echo COC files now accepted
* UPDATED: various messages regarding DNA plates
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
