$VER = "v2.03.009"
$DATE = "Date: {0}" -f (Get-Date)

$COMMENT = @"
* Fixed: progress bar randomly restarts from zero
* Fixed: order of parameters for amplicon calling reordered to match allele calling
* Fixed: Progress bar only displays during Alleles stage
* Fixed: Progess bar is no longer displayed after completion
* Fixed: info viewers called before file uploads are processed
* Fixed: alignment table has spinner while waiting for completion
* NOTE: Do not try to use multithreading in Streamlit, you cannot hold thread handles in event driven code
* TODO: coloured text in PDF report
* TODO: changing included variants in the alignment table does not redo the alignment
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
