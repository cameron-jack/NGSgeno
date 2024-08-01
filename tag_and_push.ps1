$VER = "v0.28.006"
$DATE = "Date: {0}" -f (Get-Date)

$COMMENT = @"
* FIXED: util.py dependencies in ngsmatch.py removed to fix path issue
* FIXED: removed debugging line from display of report stage (7)
* NEW: this script!
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
