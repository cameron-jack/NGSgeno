#!/usr/bin/env bash

VER="v0.26.011"
DATE="Date: $(date)"

read -r -d '' COMMENT << EOM
* Updated jsonpickle
* Fixed plate, sample, and file display selections
* Fixed experiment saving/loading.
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
