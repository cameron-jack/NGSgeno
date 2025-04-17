$VER = "v1.02.001"
$DATE = "Date: {0}" -f (Get-Date)

$COMMENT = @"
* CHANGED: Replaced biopython dependencey for FASTA parsing with Eslam's code
* CHANGED: dual bracket refseqs now treated as single bracket types
* FIXED: Fixed exact matching with variable regions
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
