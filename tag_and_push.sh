#!/usr/bin/env bash

VER="v0.26.015"
DATE="Date: $(date)"

read -r -d '' COMMENT << EOM
* Added: Update_modules.bat for updating tools automatically.
* Added: Save message at the bottom of each page.
* Added: Pinned all tools to a version.
EOM

mv changelog.txt changelog_old.txt
echo "$VER" > changelog.txt
echo "$DATE" >> changelog.txt
echo "$COMMENT" >> changelog.txt
echo "" >> changelog.txt
cat changelog_old.txt >> changelog.txt
rm changelog_old.txt

git add -u
MSG="$(echo "$COMMENT")"
git commit -m "$MSG"
git tag -a "$VER" -m "$MSG"

#git push
#git push origin $VER
