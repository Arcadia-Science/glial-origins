#!/bin/bash

# This script queries the UniProt ID Mapping API to convert
# one list of gene names into another (usually UniProtKB)

# This script takes 3 arguments:
# $1: a newline-separated txt file with a list of gene names
#     with a filename ending in '_ids.txt'
# $2: a 'from' value which states the starting name type
# $3: a 'to' value which states the return name type

# Store the ids as a variable:
ids=$(cat $1 | tr '\n' ',')

# Write the curl forms to a file:
echo "from=$2&to=$3&ids=$ids" > $1-form.txt

# And submit, then clean up
curl -d @$1-form.txt https://rest.uniprot.org/idmapping/run > $1.txt
jobID=$(cat $1.txt | sed "s/.*://g" | sed "s/}//g" | sed 's/"//g')
rm $1.txt
rm $1-form.txt
waiting=1
finished={"jobStatus":"FINISHED"}
while [ $waiting == 1 ]
do
    # Wait a few seconds
    sleep 10
    
    #Check if it's finished. 
    curl -i "https://rest.uniprot.org/idmapping/status/$jobID" | tail -n1 > $1.status    
    status=$(cat $1.status)

    # Test whether the file storing the UniProt results is empty or not. If so we'll keep going, otherwise we append to the list of mapped IDs and end. 
    if [[ "$status" == "$finished" ]]
    then
        waiting=1
    else
        waiting=0
        curl -s "https://rest.uniprot.org/idmapping/uniprotkb/results/stream/$jobID?format=tsv" > $1.txt
        
        newname=$(echo $1 | sed 's/_ids.txt//g')
        mv $1.txt ${newname}_UniProtIDs.txt
        rm $1.status
    fi
done