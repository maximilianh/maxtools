#!/bin/bash

directory="/hive/data/genomes/wuhCor1/bed/uniprot/"
file="covid-19.new.xml"
uri="ftp://ftp.uniprot.org/pub/databases/uniprot/pre_release/covid-19.xml"

cd "$directory"
if [ ! -f "$file" ]
then
    touch $file --date="2000-01-01"
fi

touch "$file.time" --reference="$file"
curl -s -o "$file" -z "$file.time" "$uri"

if [[ "$file" -nt "$file.time" ]]
then
    printf "File $file has been updated.\n"
    /hive/data/genomes/wuhCor1/bed/uniprot/doUpdate.sh
fi
