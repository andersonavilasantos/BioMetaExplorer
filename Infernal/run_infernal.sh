#!/bin/bash

grep -v "^#" $1 | awk '{print $1"/"$8"-"$9, $8, $9, $1}' | esl-sfetch -Cf $2 -> $3