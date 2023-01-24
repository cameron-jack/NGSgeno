#!/usr/bin/env bash

VER="v0.21.008"
DATE="Date: $(date)"

read -r -d '' COMMENT << EOM
* CHANGED: display viewer now controlled by checkbox
* NEW: batch file for starting startup.py
* NEW: first try at managing all uploads
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
