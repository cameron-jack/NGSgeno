#!/usr/bin/env bash

VER="v0.21.013"
DATE="Date: $(date)"

read -r -d '' COMMENT << EOM
* Fixed: missing assays that did not match primer names
* Changed: assay-primer rather than primer-assay mapping
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
