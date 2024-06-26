#!/usr/bin/env bash

VER="v0.27.000"
DATE="Date: $(date)"

read -r -d '' COMMENT << EOM
* NEW: log has reordered columns and colour
* FIXED: sample display column spacing
* NEW: on-screen message displays from anywhere in code
* FIXED: file record synching on re-uploaded files
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
