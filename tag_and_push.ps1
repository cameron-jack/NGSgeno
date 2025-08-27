$VER = "v2.01.000"
$DATE = "Date: {0}" -f (Get-Date)

$COMMENT = @"
* Incomplete: amplicon reporting, stable but unfinished
* New: added amplicon matching - exact matching and variable region matching
* New: multiple reference files are now supported, and users can select which ones they want to use for given samples
* New: updates older reference sequences format to new format exp.reference_sequences[(source_fn, purpose)]=[(id, seq), ...]
* New: amplicon reports in rich text (RTF) using rich
* New: multiple sequence alignment in variant viewer using cogent3
* Changed: info viewer is now enabled by default
* Changed: read counts now include comma separator for thousands
* Removed: miss-cache user options
* Removed: miss-cache
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
