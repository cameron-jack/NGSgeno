#!/usr/bin/env bash

VER="v0.21.006"
DATE="Date: $(date)"

read -r -d '' COMMENT << EOM
* NEW: Protection against non-ascii characters in reference sequence file
* NEW: Added checks for open files
* NEW: Increased reliability and reporting in matching code
EOM

mv changelog.txt changelog_old.txt
echo "$VER" > changelog.txt
echo "$DATE" >> changelog.txt
echo "$COMMENT" >> changelog.txt
echo "" >> changelog.txt
cat changelog_old.txt >> changelog.txt
rm changelog_old.txt

#MSG="$(printf "${COMMENT}")"
#sed -i "$MSG" changelog.txt
#echo "${COMMENT}" | sed '1s/^/1i /' | sed -i -f- changelog.txt
#echo "$DATE" | sed '1s/^/1i /' | sed -i -f- changelog.txt
#echo "$VER" | sed '1s/^/1i /' | sed -i -f- changelog.txt

git add -u
MSG="$(echo "$COMMENT")"
git commit --amend -m "$MSG"
git tag -a "$VER" -m "$MSG"

git push
git push origin $VER
