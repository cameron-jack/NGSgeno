#!/usr/bin/env bash
# set -x

VER="v0.21.003"
DATE="Date: $(date)"

read -r -d '' COMMENT << EOM
* FIXED: skipping first line of Stage3 output. Confused about lack of header\n
* FIXED: skipping some result outputs due to change of fields\n
* CHANGED: Now opens results at the end of matching, so no empty files\n
EOM

echo "$COMMENT" | sed '1s/^/1i /' | sed -i -f- changelog.txt
echo "$DATE" | sed '1s/^/1i /' | sed -i -f- changelog.txt
echo "$VER" | sed '1s/^/1i /' | sed -i -f- changelog.txt

git add -u
MSG="$(printf "${COMMENT}")"
git commit -m "$MSG"
git tag -a "$VER" -m "$MSG"

git push
git push origin $VER
