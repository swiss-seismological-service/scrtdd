#!/bin/bash

seiscomp_exec="/usr/local/bin/seiscomp exec"

# Database to read the events from, e.g.:
#   mysql://user:password@host/database
#   postgresql://user:password@host/database
#   sqlite3:///home/sysop/database.sqlite
CATALOG_DB="dbtype://user:password@host/database"

echo "Downloading events from $CATALOG_DB..."

$seiscomp_exec sclistorg -d $CATALOG_DB \
          --begin "2009-01-01 00:00:00" \
          --end "2019-01-01 00:00:00" \
          --org-type preferred \
          --manual-only \
          --inc-agency "MyAgency" \
          --ev-type "earthquake,quarry blast,explosion,induced earthquake" \
          --area 46.58,8.32,46.67,8.46 \
          --verbosity=3 --log-file sclistorg.log \
   > catalog-ids.csv

echo "Relocating events..."

$seiscomp_exec scrtdd -d $CATALOG_DB \
       --reloc-catalog catalog-ids.csv \
       --profile myProfile \
       --verbosity=3 --log-file scrtdd.log \
       --xmlout > relocated.xml

echo "event.csv, phase.csv and stations.csv contain the input catalog in csv format"
echo "reloc-event.csv, reloc-phase.csv and reloc-stations.csv contain the relocated catalog in csv format"
echo "relocated.xml contains the relocated catalog in SCML format"

