#!/usr/bin/env bash

VER="v0.26.012"
DATE="Date: $(date)"

read -r -d '' COMMENT << EOM
* Fixed: requirements file mistake
* Fixed: primer and index live displays
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
