$VER = "v2.03.003"
$DATE = "Date: {0}" -f (Get-Date)

$COMMENT = @"
* Fixed: shifted infinite progress checking loop to end of code, as recommended by Streamlit
* New: ngsmatch.py extensively uses asyncio, aiofiles to improve concurrency
* New: requires aiofiles python module
* Fixed: inexact matching no longer matches against bracket characters
* Fixed: brackets also removed from variant sequence creation
* Changed: variants.fa no longer reports hits to other targets
* Changed: Amplicon matching has its own inexact matching code path
* Changed: Reverted to sans serif fonts now that sequence viewing is no longer required
* NOTE: Do not try to use multithreading in Streamlit, you cannot hold thread handles in event driven code
* TODO: User reports that amplicon matching progress was not being shown, grey screen only
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
