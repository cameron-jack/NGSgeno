#!/usr/bin/env bash
# set -x

VER="v0.21.002"
DATE="Date: $(date)"

read -r -d '' COMMENT << EOM
* FIXED: removed the sampleNo field from FASTQ file names\n
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
