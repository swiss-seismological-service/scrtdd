/***************************************************************************
 * Copyright (C) gempa GmbH                                                *
 * All rights reserved.                                                    *
 * Contact: gempa GmbH (seiscomp-dev@gempa.de)                             *
 *                                                                         *
 * GNU Affero General Public License Usage                                 *
 * This file may be used under the terms of the GNU Affero                 *
 * Public License version 3.0 as published by the Free Software Foundation *
 * and appearing in the file LICENSE included in the packaging of this     *
 * file. Please review the following information to ensure the GNU Affero  *
 * Public License version 3.0 requirements will be met:                    *
 * https://www.gnu.org/licenses/agpl-3.0.html.                             *
 *                                                                         *
 * Other Usage                                                             *
 * Alternatively, this file may be used in accordance with the terms and   *
 * conditions contained in a signed written agreement between you and      *
 * gempa GmbH.                                                             *
 ***************************************************************************/

#include "geo.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>

#define WGS84_SEMI_MAJOR_AXIS 6378137.0
#define WGS84_FLATTENING (1.0/298.2572235630)


namespace
{
inline double rad2deg(double r) { return 180.0 * r / M_PI; }

inline double deg2rad(double d) { return M_PI * d / 180.0; }
}

namespace HDD
{

static int _delazi(double lat1,  double lon1, double lat2, double lon2,
                   double *dist, double *azi, double *baz) {
	double a,b,gam,cosa,cosb,cosc,sina,sinb,sinc,azi1,azi2,delta;
	if (lat1==lat2 && lon1==lon2) {
		*dist = *azi = *baz = 0.;
		return 0;
	}
	a     = M_PI_2 - lat2;
	b     = M_PI_2 - lat1;
	gam   = lon1 - lon2;
	cosa  = cos(a);
	cosb  = cos(b);
	sina  = sin(a);
	sinb  = sin(b);
	cosc  = cosa*cosb + sinb*sina*cos(gam);
	cosc = std::min(1.0, std::max(-1.0, cosc));
	delta = acos(cosc);
	sinc  = sin(delta);
	azi1  = acos((cosa-cosb*cosc)/(sinb*sinc));
	azi2  = acos((cosb-cosa*cosc)/(sina*sinc));

	if ((std::isnan(azi1) || std::isnan(azi2)) && fabs(lon2-lon1)<0.000001) {
		if (lat1>lat2) {
			azi1 = M_PI;
			azi2 = 0.;
		}
		else {
			azi1 = 0.;
			azi2 = M_PI;
		}
	}
	else {
		if (sin(gam)<0.)
			azi2 = M_PI+M_PI - azi2;
		else
			azi1 = M_PI+M_PI - azi1;
	}
	*dist = delta;
	*azi  = azi1;
	*baz  = azi2;
	return 0;
}

void delazi(double lat1, double lon1, double lat2, double lon2,
            double *dist, double *azi1, double *azi2) {
	double pi180 = M_PI/180., ip180 = 180./M_PI;
	_delazi(lat1*pi180, lon1*pi180, lat2*pi180, lon2*pi180,
	        dist, azi1, azi2);
// XXX	if (err) return err; // FIXME throw an exception
	*dist *= ip180;
	*azi1 *= ip180;
	*azi2 *= ip180;
}


static void mb_geocr( double lon, double lat, double *a, double *b, double *c ) {
	double blbda, bphi, ep, ug, vg;

	blbda = deg2rad(lon);
	bphi = deg2rad(lat);
	ep = 1.0 - WGS84_FLATTENING;
	ug = ep*ep*tan(bphi);
	vg = 1.0/sqrt(1.0+ug*ug);
	*a = vg*cos(blbda);
	*b = vg*sin(blbda);
	*c = ug*vg;
}


static double mb_azm(double x, double y) {
	double th;

	if ( x == 0.0 ) {
		if ( y > 0.0 ) return 90.0;
		if ( y < 0.0 ) return 270.0;
		return 0.0;
	}

	th = rad2deg(atan(fabs(y/x)));
	if ( x > 0.0 ) {
		if ( y < 0.0 )
			return 360.0-th;
		return th;
	}
	else {
		if ( y >= 0.0 )
			return 180.0-th;
		return 180.0+th;
	}
}


void delazi_wgs84(double elat, double elon, double slat, double slon,
                  double *distance, double *azim, double *bazim) {
	/* returns distance and azimuth in degrees of two locations on
	 * earth
	 *
	 * parameters of routine
	 * double     slat, slon;    input; latitude and longitude of station
	 * double     elat, elon;    input; latitude and longitude of epicentre
	 * double     *distance;     output; distance in degrees
	 * double     *azim;         output; azimuth in degrees
	 * double     *bazim;        output; back-azimuth in degrees
	 */
	double as, bs, cs, ds;
	double ae, be, ce, de;
	double bls, cbls, sbls, ble;
	double codel, bgdel;
	double xi, xj, xk;
	double sindt, cosz, sinz;

	/* check for equality */
	as = fabs(slat-elat) + fabs(slon-elon);
	if ( as < 1.0e-5 ) {
		*distance = 0.0;
		*azim = 0.0;
		*bazim = 0.0;
		return;
	}

	mb_geocr(slon, slat, &as, &bs, &cs);
	ds = sqrt(1.0 - cs*cs);
	mb_geocr(elon, elat, &ae, &be, &ce);
	de = sqrt(1.0 - ce*ce);

	bls = deg2rad(slon);
	cbls = cos(bls);
	sbls = sin(bls);
	codel = ae*as + be*bs + ce*cs;

	sindt = sqrt(1.0-codel*codel);
	if ( codel == 0.0 )
		bgdel = M_PI/2.0;
	else {
		bgdel = atan(fabs(sindt/codel));
		if ( codel <= 0.0 )
			bgdel = M_PI - bgdel;
	}

	*distance = rad2deg(bgdel);

	/* azimuths */
	xi = bs*ce - be*cs;
	xj = as*ce - ae*cs;
	xk = as*be - ae*bs;
	cosz = (xi*sbls + xj*cbls)/sindt;
	sinz = xk/(ds*sindt);
	*bazim = mb_azm(cosz, sinz);
	ble = deg2rad(elon);
	cosz = -(xi*sin(ble) + xj*cos(ble))/sindt;
	sinz = -xk/(de*sindt);
	*azim = mb_azm(cosz, sinz);
}


static int _delandaz2coord(double d, double az, double lat0, double lon0,
                           double *lat, double *lon) {
	double a, b, gam, cosa, cosb, cosd, sina, sinb, sind;

	if ( d > M_PI ) {
		d = M_PI*2 - d;
		az += M_PI;
	}

	b    = (M_PI_2 - lat0);
	cosb = cos(b);  sinb = sin(b);
	cosd = cos(d);  sind = sin(d);
	cosa = cosb*cosd + sinb*sind*cos(az);
	a    = acos(cosa);
	sina = sin(a);
	gam = (cosd - cosa*cosb)/(sina*sinb);
	if (gam > 1.) gam = 1.;
	if (gam <-1.) gam =-1.;
	gam = acos(gam);
	if (sin(az) < 0.) gam = -gam;

	*lat = M_PI_2 - a;
	*lon = fmod (lon0+gam+M_PI, 2*M_PI) - M_PI;

	return 0;
}

void delandaz2coord(double dist, double azi, double lat0, double lon0,
                    double *lat, double *lon) {
	double pi180 = M_PI/180., ip180 = 180./M_PI;
	_delandaz2coord(dist*pi180, azi*pi180, lat0*pi180, lon0*pi180, lat, lon);
// XXX	if (err) return err;
	*lat *= ip180;
	*lon *= ip180;
}

}
