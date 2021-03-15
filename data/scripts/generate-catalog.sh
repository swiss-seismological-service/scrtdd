#!/bin/bash

seiscomp_exec="/usr/local/bin/seiscomp exec"

# Database to read the events from, e.g.:
#   mysql://user:password@host/database
#   postgresql://user:password@host/database
#   sqlite3:///home/sysop/database.sqlite
CATALOG_DB="dbtype://user:password@host/database"

# Database to import the relocated catalog to (empty -> no import)
DESTINATION_DB=""

#
# Date range used for fetching events from the database
# The double-difference relocation will be peformend on those events
#
START_DATE=$(date -Iminutes -d "2009-01-01")  # set your starting date here
END_DATE=$(date -Iminutes)                    # now

# 
# Create a working directory for this relocation
#
workingdir="./$(date +%Y%m%d-%H%M -d $END_DATE)"

echo "Creating working directory $workingdir"

if [ -d $workingdir ]; then
  echo "directory $workingdir already exists: stop here"
  exit 1
fi

mkdir -p $workingdir

if [ ! -d $workingdir ]; then
  echo "Cannot create directory $workingdir: stop here"
  exit 1
fi

cd $workingdir

#
# Downloads the events
#
echo "Downloading events from $START_DATE to $END_DATE (database $CATALOG_DB)..."

ID_FILE=catalog-ids.csv

$seiscomp_exec sclistorg -d $CATALOG_DB --begin "$(date "+%Y-%m-%d %H:%M:00" -d $START_DATE)" \
          --end "$(date "+%Y-%m-%d %H:%M:00" -d $END_DATE)" \
          --org-type preferred \
   > $ID_FILE

if [ $? -ne 0 ] || [ ! -f $ID_FILE ]; then
  echo "Catalog not downloaded: stop here"
  exit 1
fi

#
# Relocate catalog
#
echo "Relocating events..."

#depending on the size of the logs, many files will be generated in the form scrtdd.log scrtdd.log.1 scrtdd.log.2 ...
RTDDLOG_FILE=scrtdd.log

$seiscomp_exec scrtdd -d $CATALOG_DB --reloc-catalog $ID_FILE --profile myProfile \
       --verbosity=3 --log-file $RTDDLOG_FILE

if [ $? -ne 0 ] || [ ! -f reloc-event.csv ] || [ ! -f reloc-phase.csv ] || [ ! -f reloc-station.csv ]; then
  echo "Catalog not relocated: stop here"
  exit 1
fi

#
# Convert catalog to Seiscomp XML
#
XMLRELOC_FILE=relocated.xml

echo "Converting relocated events to XML file $XMLRELOC_FILE..."

$seiscomp_exec scrtdd --dump-catalog-xml reloc-station.csv,reloc-event.csv,reloc-phase.csv > $XMLRELOC_FILE

if [ $? -ne 0 ] || [ ! -f $RELOC_FILE ]; then
  echo "Cannot convert the relocated catalog to XML: stop here"
  exit 1
fi

#
# Import the relocated catalog in a database (if destination db is defined)
#
if [ -n "${DESTINATION_DB}" ]; then

  echo "Importing $XMLRELOC_FILE to $DESTINATION_DB ..."

  $seiscomp_exec scdb -i  $XMLRELOC_FILE -d $DESTINATION_DB

  if [ $? -ne 0 ]; then
    echo "Cannot import the relocated catalog: stop here"
    exit 1
  fi

fi

#
# clean-up of useless files (no need to save the xml since we can easily re-create it)
# and keep a tar files of the results
#
rm -f $XMLRELOC_FILE
cd ..

echo "Creating backup file $workingdir.tar.bz2..."
tar cjSvf $workingdir.tar.bz2 $workingdir

echo "Deleting folder $workingdir..."
rm -rf $workingdir

