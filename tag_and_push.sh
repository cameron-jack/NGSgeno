#!/usr/bin/env bash

VER="v0.21.020"
DATE="Date: $(date)"

read -r -d '' COMMENT << EOM
* Fixed: pcrplate and pcrwell should be pcrPlate and pcrWell respectively
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
git commit --amend -m "$MSG"
git tag -a "$VER" -m "$MSG"

git push
git push origin $VER
