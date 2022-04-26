#!/usr/bin/python3

import sys
import math
import pandas as pd
import geopy.distance


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


if len(sys.argv) != 3:
    eprint("compare-catalogs.py event1.csv event2.csv")
    exit(0)

EV_RANGE = pd.Timedelta('0.5 sec')

loc1 = pd.read_csv(sys.argv[1], parse_dates=['isotime'])
loc2 = pd.read_csv(sys.argv[2], parse_dates=['isotime'])

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

    coords_1 = (ev1.latitude, ev1.longitude)
    coords_2 = (ev2.iloc[0].latitude, ev2.iloc[0].longitude)
    coord_dist = geopy.distance.distance(coords_1, coords_2).km * 1000
    depth_diff = abs(ev1.depth - ev2.iloc[0].depth) * 1000
    distance = (coord_dist**2 + depth_diff**2)**(1. / 2)

    print("%s,%s,%.6f,%.6f,%.3f,%.2f,%.3f,%.1f,%.1f,%.1f" %
          (ev1.id, ev1.isotime, ev1.latitude, ev1.longitude, ev1.depth, ev1.magnitude,
           (ev1.isotime - ev2.iloc[0].isotime).total_seconds(), coord_dist, depth_diff, distance)
          )
