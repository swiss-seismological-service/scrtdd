#!/bin/bash

seiscomp_exec="/usr/local/bin/seiscomp exec"

# Database to read the events from, e.g.:
#   mysql://user:password@host/database
#   postgresql://user:password@host/database
#   sqlite3:///home/sysop/database.sqlite
CATALOG_DB="dbtype://user:password@host/database"

# Database to which import the relocated catalog
# (empty string -> do nothing)
DESTINATION_DB=""

# The real-time rtDD configuration folder to which copy the relocated catalog
# The real-time rtDD profile should point to this catalog to get automatic updates
# (empty string -> do nothing)
RTDD_BGCAT_DIR=""

# Folder to which copy the data for the web page interactive map
# You need a web server to visualize the page on a different machine or you
# can run a browser on the same machine and open the relocation-map.html file
# (empty string -> do nothing)
WEB_DIR=""

RTDD_PROFILE="myProfile"

LOGFLAG="--console=1 --verbosity=3"

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
# Downloads the events (configure sclistorg options as you please)
#
echo "Downloading events from $CATALOG_DB..."

ID_FILE=catalog-ids.csv

$seiscomp_exec sclistorg -d $CATALOG_DB \
          --begin "2000-01-01 00:00:00" \
          --end "$(date "+%Y-%m-%d %H:%M:00" -d $NOW)" \
          --org-type preferred \
          $LOGFLAG \
   > $ID_FILE

if [ $? -ne 0 ] || [ ! -f $ID_FILE ]; then
  echo "Catalog not downloaded: stop here"
  exit 1
fi

echo "Done: created file $ID_FILE"

#
# Relocate catalog
#
echo "Relocating events with profile $RTDD_PROFILE..."

XMLRELOC_FILE=relocated.xml

$seiscomp_exec scrtdd -d $CATALOG_DB \
       --reloc-catalog $ID_FILE \
       --profile $RTDD_PROFILE \
       --verbosity=3 --log-file scrtdd.log \
       --xmlout > $XMLRELOC_FILE

if [ $? -ne 0 ] || [ ! -f reloc-event.csv ] || [ ! -f reloc-phase.csv ] || [ ! -f reloc-station.csv ]; then
  echo "Catalog not relocated: stop here"
  exit 1
fi

if [ "$( wc --lines < reloc-event.csv )" -lt 2 ]; then
  echo "Relocated catalog contains no events: stop here"
  exit 1
fi

echo "Done: created files $XMLRELOC_FILE and reloc-event.csv,reloc-phase.csv,reloc-station.csv"

#
# Copy relocated catalog to real-time rtDD configuration
#
if [ -n "${RTDD_BGCAT_DIR}" ]; then

  echo "Copying reloc-station.csv reloc-phase.csv reloc-event.csv to $RTDD_BGCAT_DIR"
  
  cp -b -S .old -f reloc-station.csv reloc-phase.csv reloc-event.csv $RTDD_BGCAT_DIR/
  if [ $? -ne 0 ]; then
    echo "Cannot copy the relocated catalog to $RTDD_BGCAT_DIR: stop here"
    exit 1
  fi

  # Force scrtdd to reload the profile background catalog
  echo "Send message to scrtdd forcing a reloading of the background catalog for profile $RTDD_PROFILE"
  $seiscomp_exec scrtdd --send-reload-profile-msg $RTDD_PROFILE --user rtddBgCatUpdate $LOGFLAG
fi

#
# Copy the relocated catalog into the web folder
#
if [ -n "${WEB_DIR}" ]; then
  cp -f event.csv  $WEB_DIR/event.csv
  cp -f station.csv  $WEB_DIR/station.csv
  cp -f reloc-event.csv $WEB_DIR/me-dd-event.csv
  cp -f ../relocation-map.html $WEB_DIR/
  date > $WEB_DIR/LAST_RUN
fi

#
# Import the relocated catalog into the destination database
#
if [ -n "${DESTINATION_DB}" ]; then

  if [ ! -f $XMLRELOC_FILE ]; then
    echo "Cannot find the relocated XML catalog $XMLRELOC_FILE: stop here"
    exit 1
  fi

  # Optional: run scevent on $XMLRELOC_FILE
  # $seiscomp_exec scevent --ep $XMLRELOC_FILE  > relocated-with-events.xml

  # Optional: clean the database from the previous data
  # $seiscomp_exec scdbstrip -d $DESTINATION_DB $LOGFLAG --days 0
  # $seiscomp_exec scdbstrip -d $DESTINATION_DB $LOGFLAG --days 0 --check --clean-unused

  echo "Importing $XMLRELOC_FILE into $DESTINATION_DB ..."

  $seiscomp_exec scdb -i $XMLRELOC_FILE \
                      -d $DESTINATION_DB \
                      --verbosity=3 --log-file scdb.log

  if [ $? -ne 0 ]; then
    echo "Errors while importing the relocated catalog into database"
  fi

  echo "Done: imported $XMLRELOC_FILE into $DESTINATION_DB"
fi

#
# Back-up the results in a compressed tar file and then delete the working directory
#
cd ..

echo "Creating backup file $WORKINGDIR.tar.bz2..."
tar cjSvf $WORKINGDIR.tar.bz2 $WORKINGDIR

echo "Deleting folder $WORKINGDIR..."
rm -rf $WORKINGDIR

