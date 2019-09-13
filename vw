#! /bin/bash
journal=$(../pick $1)
ret=$?
if [[ $ret != 0 ]]; then # check return code
    exit $ret
fi

jobid=$(grep "JobID := " $journal -r)
jobid=${jobid:$(expr length "JobID := ")} # substring ${substr:position:length}

datadir=$(grep "DataDir := " $journal -r)
datadir=${datadir:$(expr length "DataDir := ")} # substring ${substr:position:length}

# open files in different tabs in ReadOnly mode
vim -R -O $(find oe/ -type f -name "*${jobid}*") -c "tabnew $datadir/vitals.txt" -c "tabr" -c "command! -nargs=1 Log tabnew $datadir/logs/rank<args>.log"
