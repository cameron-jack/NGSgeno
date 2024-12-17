$VER = "v0.28.019"
$DATE = "Date: {0}" -f (Get-Date)

$COMMENT = @"
* FIXED: Sample viewer code was incorrectly looking up barcodes in sample wells
* FIXED: All DNA plates were being added again to the list of sample plates, not just the separately loaded 384-well ones
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
