#!/usr/bin/env bash

VER="v0.28.000"
DATE="Date: $(date)"

read -r -d '' COMMENT << EOM
* NEW: all uploads are queued and processed one at a time, separating the GUI from the parsing
* NEW: multi-file plate uploads now get a file-plate-purpose mapping to solve conflicts
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
