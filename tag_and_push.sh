#!/usr/bin/env bash

VER="v0.23.000"
DATE="Date: $(date)"

read -r -d '' COMMENT << EOM
* Fixed: clear old reference sequence entries on loading new references
* Fixed: sets html encoding for Windows plate viewer
* New: rewrite of matching code for correctness and speed
* TODO: one bug in matching code, likely of inexact matches
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
