#! /bin/bash
# locate journal.txt. Proceed if and only if one match is found
# TODO if $1 is missing, open the last one
cmd='grep -l "$1" .stages/*/journal.txt -r'
numf=$(eval $cmd | wc -l)
if (( numf > 1)); then
    echo "ERROR : matched more than 1 journal.txt"
    eval $cmd
    exit 0
elif (( numf == 0 )); then
    echo "ERROR : no matched journal.txt."
    echo "NOTE filenames are not matched against."
    exit 0
fi
journal=$(eval $cmd)

jobid=$(grep "JobID := " $journal -r)
jobid=${jobid:$(expr length "JobID := ")} # substring ${substr:position:length}

datadir=$(grep "DataDir := " $journal -r)
datadir=${datadir:$(expr length "DataDir := ")} # substring ${substr:position:length}

# check file existence
vim $(grep "$jobid" oe/*) $datadir/"vitals.txt"
