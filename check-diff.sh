#!/bin/bash -e
src="${1:src}"
for fic in $(ls ${src}/*)
do
    clang-format -style=file "${fic}" | diff -U5 "${fic}" -
done
