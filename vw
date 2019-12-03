#! /bin/bash
journal=$(../pick $1)
ret=$?
if [[ $ret != 0 ]]; then # check return code
    exit $ret
fi

jobid=$(grep "JobID := " $journal -r)
jobid=${jobid:$(expr length "JobID := ")} # substring ${substr:position:length}

# extract the pure number ID, which is presumably the longest number in the string. Copied from https://unix.stackexchange.com/questions/127009/how-can-i-print-the-longest-number-in-a-string
jobid=$(echo $jobid |
            awk '{gsub("[^0-9]+","\n"); print;}' |
            awk '{ if (length($0) > max) {max = length($0); maxline = $0} }
  END { print maxline }')

datadir=$(grep "DataDir := " $journal -r)
datadir=${datadir:$(expr length "DataDir := ")} # substring ${substr:position:length}

# open files in different tabs in ReadOnly mode
vim -R -O $(find oe/ -type f -name "*${jobid}*") -c "tabnew $datadir/vitals.txt" -c "tabr" -c "command! -nargs=1 Log tabnew $datadir/logs/rank<args>.log"
