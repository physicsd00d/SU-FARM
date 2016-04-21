#!/bin/bash
# This script is to run the aiaaSpace file from the command line because there is a massive memory leak in FACET so it must be completely shut down periodically.

#java -Djava.library.path=/Users/tcolvin1/Documents/testWorkspace/testFACET/bin -cp /Users/tcolvin1/Documents/testWorkspace/testFACET/target/classes:/Users/tcolvin1/Documents/testWorkspace/testFACET/bin/jlib/*:./ aiaaSpace.AiaaSpace

# Here is how to run the script with arguments
FACET=/Users/tcolvin1/Documents/testWorkspace/testFACET
run="java -Djava.library.path=$FACET/bin -Xmx768m -cp $FACET/target/classes:$FACET/bin/jlib/*:./ aiaaSpace.AiaaSpace"

#dataFolderName="/Users/tcolvin1/Desktop/aiaaSpace/PrunedTRX/runTheseFirstMonth/"
#outFolderName="/Users/tcolvin1/Desktop/aiaaSpace/MissionResults/BaselineOneMonth/"

dataFolderName="/Users/tcolvin1/Desktop/aiaaSpace/PrunedTRX/runTheseFebMarch/"
outFolderName="/Users/tcolvin1/Desktop/aiaaSpace/MissionResults/BaselineFebMarch/"

dataList=`ls $dataFolderName*.trx`

# Don't need the folder prefixes, just strip out the basenames
dataList=`basename $dataList`


for curData in $dataList; do
  dataParts=(${curData//_/ })
  dataLocation=${dataParts[1]}
  dataDate=${dataParts[2]}
  dataDate=${dataDate::${#dataDate}-4}  # need to get rid of the .trx at the end of the date.

#    # Skip over data files before this date
#    if [[ $dataDate < 20130122 ]]; then
#      continue
#    fi

#  echo "$outFolderName" "$dataFolderName$curData"
  eval $run "$outFolderName" "$dataFolderName$curData"
done











