#!/usr/bin/env bash

VER="v0.28.002"
DATE="Date: $(date)"

read -r -d '' COMMENT << EOM
* FIXED: Rodentity and custom upload functions adapted to new file handling
* FIXED: auto deleting filenames/queues starting with underscore
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
