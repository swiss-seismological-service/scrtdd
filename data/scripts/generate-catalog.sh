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
# Create a working directory for this relocation
#
NOW=$(date -Iminutes)

WORKINGDIR="./$(date +%Y%m%d-%H%M -d $NOW)"

echo "Creating working directory $WORKINGDIR"

if [ -d $WORKINGDIR ]; then
  echo "directory $WORKINGDIR already exists: stop here"
  exit 1
fi

mkdir -p $WORKINGDIR

if [ ! -d $WORKINGDIR ]; then
  echo "Cannot create directory $WORKINGDIR: stop here"
  exit 1
fi

cd $WORKINGDIR

#
# Downloads the events
#
echo "Downloading events from $CATALOG_DB)..."

ID_FILE=catalog-ids.csv

$seiscomp_exec sclistorg -d $CATALOG_DB \
          --begin "2000-01-01 00:00:00" \
          --end "$(date "+%Y-%m-%d %H:%M:00" -d $NOW)" \
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
XMLRELOC_FILE=relocated.xml

$seiscomp_exec scrtdd -d $CATALOG_DB --reloc-catalog $ID_FILE \
       --profile myMultiEventProfile \
       --verbosity=3 --log-file $RTDDLOG_FILE \
       --xmlout > $XMLRELOC_FILE

if [ $? -ne 0 ] || [ ! -f reloc-event.csv ] || [ ! -f reloc-phase.csv ] || [ ! -f reloc-station.csv ]; then
  echo "Catalog not relocated: stop here"
  exit 1
fi

#
# Copy relocated catalog to real-time rtDD configuration
#
RTDD_BGCAT_DIR=~/.seiscomp/scrtdd/myRealTimeProfile/bgCat

echo "Copying reloc-station.csv reloc-phase.csv reloc-event.csv to $RTDD_BGCAT_DIR"

cp -b -S .old -f reloc-station.csv reloc-phase.csv reloc-event.csv $RTDD_BGCAT_DIR/
if [ $? -ne 0 ]; then
  echo "Cannot copy the relocated catalog to $RTDD_BGCAT_DIR: stop here"
  exit 1
fi
# Force scrtdd to reload the profile background catalog
$seiscomp_exec scrtdd --send-reload-profile-msg myRealTimeProfile --user rtddBgCatUpdated

#
# Import the relocated catalog in a database (if destination db is defined)
#
if [ -n "${DESTINATION_DB}" ]; then

  if [ $? -ne 0 ] || [ ! -f $XMLRELOC_FILE ]; then
    echo "Cannot find the relocated XML catalog: stop here"
    exit 1
  fi

  echo "Importing $XMLRELOC_FILE to $DESTINATION_DB ..."

  $seiscomp_exec scdb -i  $XMLRELOC_FILE -d $DESTINATION_DB

  if [ $? -ne 0 ]; then
    echo "Cannot import the relocated catalog: stop here"
    exit 1
  fi

  # No need to keep the xml file since we can easily re-create it
  rm -f $XMLRELOC_FILE 
fi

#
# Back-up the results in a compressed tar file and then delete the working directory
#
cd ..

echo "Creating backup file $WORKINGDIR.tar.bz2..."
tar cjSvf $WORKINGDIR.tar.bz2 $WORKINGDIR

echo "Deleting folder $WORKINGDIR..."
rm -rf $WORKINGDIR

