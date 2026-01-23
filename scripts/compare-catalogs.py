#!/usr/bin/python3

import sys
import math
import pandas as pd
import pyproj

geod = pyproj.Geod(ellps='WGS84')


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


if len(sys.argv) != 4:
    eprint("compare-catalogs.py event1.csv event2.csv ev-search-time")
    eprint("ev-search-time is input to pd.Timedelta e.g. '0.5 sec'")
    exit(0)

loc1 = pd.read_csv(sys.argv[1], parse_dates=['isotime'])
loc2 = pd.read_csv(sys.argv[2], parse_dates=['isotime'])

EV_RANGE = pd.Timedelta(sys.argv[3])  # e.g. '0.5 sec'


def locationDifference(lat1, lon1, kmdep1, lat2, lon2, kmdep2):
    # horizontal distance in meters
    azimuth, backaz, hDistance = geod.inv(
        lon1, lat1, lon2, lat2, radians=False)
    vDistance = abs(kmdep1 * 1000. - kmdep2 * 1000.)  # km -> meter
    distance = (hDistance**2 + vDistance**2)**(1. / 2)
    return (hDistance, vDistance, distance)


print("id,isotime,latitude,longitude,depth,magnitude,time-diff(sec),epi-diff(m),depth-diff(m),hypo-diff(m)")

for ev1 in loc1.itertuples():

    int_start = ev1.isotime - EV_RANGE
    int_end = ev1.isotime + EV_RANGE

    ev2 = loc2[(loc2['isotime'] > int_start) & (loc2['isotime'] < int_end)]

    if ev2.empty:
        eprint("Only in %s,%s,%.6f,%.6f,%.3f,%.2f" %
               (sys.argv[1], ev1.isotime, ev1.latitude, ev1.longitude, ev1.depth, ev1.magnitude))
        continue

    if ev2.shape[0] != 1:
        eprint(
            'Multiple events found in %s around time %s' %
            (sys.argv[2], ev1.isotime))
        eprint(ev2)
        continue

    hDistance, vDistance, distance = locationDifference(
        ev1.latitude, ev1.longitude, ev1.depth,
        ev2.iloc[0].latitude, ev2.iloc[0].longitude, ev2.iloc[0].depth
    )
    print("%s,%s,%.6f,%.6f,%.3f,%.2f,%.3f,%.1f,%.1f,%.1f" %
          (ev1.id, ev1.isotime, ev1.latitude, ev1.longitude, ev1.depth, ev1.magnitude,
           (ev1.isotime - ev2.iloc[0].isotime).total_seconds(), hDistance, vDistance, distance)
          )
