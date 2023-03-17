#!/usr/bin/env bash

VER="v0.21.010"
DATE="Date: $(date)"

read -r -d '' COMMENT << EOM
* FIXED: Not deleting obsolete pending_steps
* FIXED: Displaying obsolete pending files
* CHANGED: All parsing of text files ignores non-ASCII characters
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
