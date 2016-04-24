#!/bin/bash
# This script is to run the aiaaSpace file from the command line because there is a massive memory leak in FACET so it must be completely shut down periodically.

#java -Djava.library.path=/Users/tcolvin1/Documents/testWorkspace/testFACET/bin -cp /Users/tcolvin1/Documents/testWorkspace/testFACET/target/classes:/Users/tcolvin1/Documents/testWorkspace/testFACET/bin/jlib/*:./ aiaaSpace.AiaaSpace

# Here is how to run the script with arguments
FACET=/Users/tcolvin1/Documents/testWorkspace/testFACET
run="java -Djava.library.path=$FACET/bin -Xmx768m -cp $FACET/target/classes:$FACET/bin/jlib/*:./ aiaaSpace.AiaaSpace"
#eval $run data1 data2

dataFolderName="/Users/tcolvin1/Desktop/aiaaSpace/PrunedTRX/runTheseFebMarch/"

#missionSuaFolderName="/Users/tcolvin1/Desktop/aiaaSpace/SuaFilesBrokenOut/FacetEnvelopes/"
#outFolderName="/Users/tcolvin1/Desktop/aiaaSpace/MissionResults/EnvelopesFebMarch/"

missionSuaFolderName="/Users/tcolvin1/Desktop/aiaaSpace/SuaFilesBrokenOut/TraditionalSUAs/"
outFolderName="/Users/tcolvin1/Desktop/aiaaSpace/MissionResults/TraditionalFebMarch/"

dataList=`ls $dataFolderName*.trx`
missionList=`ls $missionSuaFolderName\SUA_*`

# Don't need the folder prefixes, just strip out the basenames
dataList=`basename $dataList`


# For each SUA_ file, find the trx files that go with it
for curMission in $missionList; do
  curSUA=`basename $curMission`
  missionParts=(${curSUA//_/ })  # I think this is replacing all underscores with spaces, which bash will then parse as separate strings
  launchLocation=${missionParts[3]}
  echo $launchLocation

#  # Only run this mission file.
#  if [[ $curSUA != 'SUA_2018H_SS2_America' ]]; then
#    continue
#  fi

  for curData in $dataList; do
    dataParts=(${curData//_/ })
    dataLocation=${dataParts[1]}
    dataDate=${dataParts[2]}
    dataDate=${dataDate::${#dataDate}-4}  # need to get rid of the .trx at the end of the date.

#    # Skip over data files before this date
#    if [[ $dataDate < 20130122 ]]; then
#      continue
#    fi

    if grep -q $launchLocation <<<$dataLocation; then
      eval $run "$outFolderName" "$missionSuaFolderName$curSUA" "$dataFolderName$curData"
#      echo "$curData $dataDate"
 #     echo "$curSUA $curData"
    fi
  done
done











