#!/usr/bin/env bash

VER="v0.26.007"
DATE="Date: $(date)"

read -r -d '' COMMENT << EOM
* Fixed: transaction handling
* Fixed: Slowdowns from prints and saves
* Fixed: robustness of assay list uploads
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
