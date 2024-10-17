$VER = "v0.28.017"
$DATE = "Date: {0}" -f (Get-Date)

$COMMENT = @"
* Final preproduction fixes - allele calling
* FIXED: Splitting inexact matching by in-group and out-group caused some sequences to be incorrectly assigned
* NEW: Various parameters for matching have been exposed in the user interface
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
