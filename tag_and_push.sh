#!/usr/bin/env bash
# set -x

VER="v0.21.004"
DATE="Date: $(date)"

read -r -d '' COMMENT << EOM
* CHANGED: custom manifest samples are now always guarded c

EOM

echo "$COMMENT" | sed '1s/^/1i /' | sed -i -f- changelog.txt
echo "$DATE" | sed '1s/^/1i /' | sed -i -f- changelog.txt
echo "$VER" | sed '1s/^/1i /' | sed -i -f- changelog.txt

git add -u
MSG="$(printf "${COMMENT}")"
git commit -m "$MSG"
#git tag -a "$VER" -m "$MSG"

git push
git push origin $VER
