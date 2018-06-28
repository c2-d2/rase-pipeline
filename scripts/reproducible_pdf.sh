#! /usr/bin/env bash

set -e
set -o pipefail
set -u
#set -f

tmp="$(mktemp).pdf"

cat "$1" \
    | perl -pe 's/^\/CreationDate.*$/\/CreationDate (D:20180101000000)/g' \
    | perl -pe 's/^\/ModDate.*$/\/ModDate (D:20180101000000)/g' \
    > "$tmp"
mv "$tmp" "$1"

