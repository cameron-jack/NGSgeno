#!/usr/bin/env bash

VER="v0.23.001"
DATE="Date: $(date)"

read -r -d '' COMMENT << EOM
* Fixed: bug in parsing of aligner output
* Fixed: process handling
* Fixed: launch and match progress display
* Fixed: integration with browser
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

#git push
#git push origin $VER
