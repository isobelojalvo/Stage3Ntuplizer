#!/bin/bash
while IFS='' read -r line || [[ -n "$line" ]]; do
    echo "$line"
    ./dascli.py --query="lumi file=$line"
done < "$1"