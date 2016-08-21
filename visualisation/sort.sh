#!/bin/bash

myArray=( "$@" )
directory="${myArray[0]}"
SOURCE="${myArray[0]}"
DEST='new'

I=0
for FILE in "$SOURCE"; do
    (( D = 1 + ++I / 1000 ))
    [[ -d $DEST/$D ]] || mkdir -p "$DEST/$D"  ## You can just skip dir checking but that would be slow.
    cp -v "$FILE" "$DEST/$D/m_$I"
done