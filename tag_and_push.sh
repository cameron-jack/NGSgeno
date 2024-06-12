#!/usr/bin/env bash

VER="v0.26.016"
DATE="Date: $(date)"

read -r -d '' COMMENT << EOM
* CHANGED: parsers now overwrite existing plate entries by default unless they have different purpose
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
