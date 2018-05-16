#! /usr/bin/env bash


set -e
set -o pipefail

# $1 - file
# $2 - expected error code
test () {
    file="$1"
    experr="$2"
    echo "Testing $x"
    tmp=$(mktemp)
    ((./$file 2>&1)> "$tmp" || (code=$? && exit $(($code - $experr))) ) || ( echo "...failed ($code)" && cat "$tmp"  && exit 1)
}

DIR=`dirname $0`

# argparse should return exit status 2
for x in "$DIR"/*.py; do
    test "$x" 2
done

# R should exit with 1
for x in *.R; do
    test "$x" 1
done
