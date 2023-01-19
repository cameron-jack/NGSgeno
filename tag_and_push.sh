#!/usr/bin/env bash

VER="v0.21.007"
DATE="Date: $(date)"

read -r -d '' COMMENT << EOM
* NEW: ngsmatch.py creates debug.log in run directory if using -d
* FIXED: one too few parameters passed to process_well()
* FIXED: corrected hard coded paths to bbmap components
* CHANGED: improved debugging message handling
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
