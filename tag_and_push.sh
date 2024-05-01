#!/usr/bin/env bash

VER="v0.26.006"
DATE="Date: $(date)"

read -r -d '' COMMENT << EOM
* New: primer and index file uploads now protect against incorrect file upload
* Fixed: When files are deleted, associated plates are soft deleted too
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
