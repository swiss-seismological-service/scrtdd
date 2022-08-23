
/*--------------------------------------------------------------------
    WARNING - this file is not gmt_map.c, it contains a subset of gmt_map.c

    subset of functions selected by A. Lomax, June 1998

    changes or additions indicated by                                 */

/*AJL ... AJL*/

/*AJL*/ /*END AJL*/

/*--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 *    The GMT-system:	@(#)gmt_map.c	2.56  09 Aug 1995
 *
 *    Copyright (c) 1991-1995 by P. Wessel and W. H. F. Smith
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*
 *
 *			G M T _ M A P . C
 *
 *- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * gmt_map.c contains code related to coordinate transformation
 *
 * 	Map Transformation Setup Routines
 *	These routines initializes the selected map transformation
 *	The names and main function are listed below
 *	NB! Note that the transformation function does not check that they are
 *	passed valid lon,lat numbers. I.e asking for log10 scaling using values
 *	<= 0 results in problems.
 *
 * Map_projections include functions that will set up the transformation
 * between xy and latlon for several map projections.
 *
 * A few of the core coordinate transformation functions are based on similar
 * FORTRAN routines written by Pat Manley, Doug Shearer, and Bill Haxby, and
 * have been rewritten in C and subsequently streamlined.  The rest is based
 * on P. Snyder, "Map Projections - a working manual", USGS Prof paper 1395.
 *
 * Transformations supported (both forward and inverse):
 *	Linear x/y[/z] scaling
 *	Polar (theta, r) scaling
 *	Mercator
 *	Stereographic
 *	Albers Equal-Area Conic
 *	Lambert Conformal Conic
 *	Cassini Cylindrical
 *	Oblique Mercator
 *	TM Transverse Mercator
 *	UTM Universal Transverse Mercator
 *	Lambert Azimuthal Equal-Area
 *	Mollweide Equal-Area
 *	Hammer-Aitoff Equal-Area
 *	Sinusoidal Equal-Area
 *	Winkel Tripel
 *	Orthographic
 *	Azimuthal Equidistant
 *	Robinson
 *	Eckert IV
 *	Cylindrical Equal-area (e.g., Peters, Gall, Behrmann)
 *	Cylindrical Equidistant (Plate Carree)
 *
 * The ellipsoid used is selectable by editing the .gmtdefaults in your
 * home directory.  If no such file, create one by running gmtdefaults.
 *
 * Usage: Initialize system by calling map_setup (separate module), and
 * then just use geo_to_xy() and xy_to_geo() functions.
 *
 * Author:	Paul Wessel
 * Date:	27-JUL-1992
 * Version:	v2.1
 *
 *
 * Functions include:
 *
...

 */

/*AJL*/

/*#include "gmt.h"*/

#include <stdexcept>
#include <string>

#include <cmath>
#include <cstring>

#include "map_project.h"

namespace HDD {
namespace GMT {

#define PI_2 (2.0 * M_PI)
#define D2R (M_PI / 180.0)
#define R2D (180.0 / M_PI)

#define SMALL 1.0e-10

/*#define d_log(x) ((x) <= 0.0 ? gmt_NaN : log (x))*/
#define d_log(x) ((x) <= 0.0 ? -1.0e10 : log(x))
#define d_sqrt(x) ((x) < 0.0 ? 0.0 : sqrt(x))

/*double gmt_NaN;*/

struct ELLIPSOID
{
  char name[20];
  int date;
  double eq_radius;
  double pol_radius;
  double flattening;
};

/* Information about a particular ellipsoid */
/* Table taken from Snyder "Map projection - a working manual",
                p 12 Table 1 */

#define N_ELLIPSOIDS 15

struct ELLIPSOID ellipse[N_ELLIPSOIDS] = {
    // name, date, eq_radius, pol_radius, flattening
    {"WGS-84", 1984, 6378137.0, 6356752.1, 1.0 / 298.254},
    {"GRS-80", 1980, 6378137.0, 6356752.3, 1.0 / 298.257},
    {"WGS-72", 1972, 6378135.0, 6356750.5, 1.0 / 298.26},
    {"Australian", 1965, 6378160.0, 6356774.7, 1.0 / 298.25},
    {"Krasovsky", 1940, 6378245.0, 6356863.0, 1.0 / 298.3},
    {"International", 1924, 6378388.0, 6356911.9, 1.0 / 297.0},
    {"Hayford-1909", 1909, 6378388.0, 6356911.9, 1.0 / 297.0},
    {"Clarke-1880", 1880, 6378249.1, 6356514.9, 1.0 / 293.46},
    {"Clarke-1866", 1866, 6378206.4, 6356583.8, 1.0 / 294.98},
    {"Airy", 1830, 6377563.4, 6356256.9, 1.0 / 299.32},
    {"Bessel", 1841, 6377397.2, 6356079.0, 1.0 / 299.15},
    {"Hayford-1830", 1830, 6377276.3, 6356075.4, 1.0 / 300.80},
    {"Sphere", 1980, 6371008.7714, 6371008.7714, 0.0},
    /* https://gisgeography.com/geodetic-datums-nad27-nad83-wgs84/
    NAD27 Datum vs NAD83 Datum
    The NAD27 datum was based on the Clarke Ellipsoid of 1866:
    Semi-major axis: 6,378,206.4 m
    Semi-minor axis: 6,356,583.8 m
    Inverse flattening: 294.98
    The NAD83 datum was based on the Geodetic Reference System (GRS80)
    Ellipsoid: Semi-major axis: 6,378,137.0 m Semi-minor axis: 6,356,752.3 m
    Inverse flattening: 298.26
     */
    {"NAD-27", 1927, 6378206.4, 6356583.8, 1.0 / 294.98},
    {"NAD-83", 1983, 6378137.4, 6356752.3, 1.0 / 294.26}

};

/* set constants that were set in GMT function map_setup */

/* use values from gmt_defaults.h: */

Ellip map_setup_proxy(const char *ellipsoid_name)
{
  int num_ellipsoid;

  /* determine ellipsoid */
  for (num_ellipsoid = 0; num_ellipsoid < N_ELLIPSOIDS &&
                          strcmp(ellipsoid_name, ellipse[num_ellipsoid].name);
       num_ellipsoid++)
    ;
  if (num_ellipsoid == N_ELLIPSOIDS)
    throw std::runtime_error("Invalid ellipsoid " +
                             std::string(ellipsoid_name));

  Ellip e{};
  e.EQ_RAD = ellipse[num_ellipsoid].eq_radius;
  double f = ellipse[num_ellipsoid].flattening;
  e.ECC2   = 2 * f - f * f;
  e.ECC4   = e.ECC2 * e.ECC2;
  e.ECC6   = e.ECC2 * e.ECC4;
  e.ECC    = d_sqrt(e.ECC2);
  return e;
}

/*END AJL*/

/*
 *	TRANSFORMATION ROUTINES FOR THE LAMBERT CONFORMAL CONIC PROJECTION
 */

/*** function to set up a Lambert Conformal Conic projection */

LAMBERT vlamb(const char *ellipsoid_name,
              double rlong0,
              double rlat0,
              double pha,
              double phb)
{
  LAMBERT proj{};
  proj.e = map_setup_proxy(ellipsoid_name);
  proj.pha = pha;
  proj.phb = phb;

  double t_pha, m_pha, t_phb, m_phb, t_rlat0;

  proj.NorthPole = (rlat0 > 0.0);
  // AJL 20090812 BUG FIX
  proj.Pole = (proj.NorthPole) ? 90.0 : -90.0;
  // proj.Pole = (NorthPole) ? 90.0 : -90.0;
  pha *= D2R;
  phb *= D2R;

  t_pha = tan(45.0 * D2R - 0.5 * pha) /
          pow((1.0 - proj.e.ECC * sin(pha)) / (1.0 + proj.e.ECC * sin(pha)),
              0.5 * proj.e.ECC);
  m_pha = cos(pha) / d_sqrt(1.0 - proj.e.ECC2 * pow(sin(pha), 2.0));
  t_phb = tan(45.0 * D2R - 0.5 * phb) /
          pow((1.0 - proj.e.ECC * sin(phb)) / (1.0 + proj.e.ECC * sin(phb)),
              0.5 * proj.e.ECC);
  m_phb   = cos(phb) / d_sqrt(1.0 - proj.e.ECC2 * pow(sin(phb), 2.0));
  t_rlat0 = tan(45.0 * D2R - 0.5 * rlat0 * D2R) /
            pow((1.0 - proj.e.ECC * sin(rlat0 * D2R)) /
                    (1.0 + proj.e.ECC * sin(rlat0 * D2R)),
                0.5 * proj.e.ECC);

  if (pha != phb)
    proj.LambertConfConic_N =
        (d_log(m_pha) - d_log(m_phb)) / (d_log(t_pha) - d_log(t_phb));
  else
    proj.LambertConfConic_N = sin(pha);

  proj.LambertConfConic_F =
      m_pha / (proj.LambertConfConic_N * pow(t_pha, proj.LambertConfConic_N));
  proj.CentralMeridian       = rlong0;
  proj.LambertConfConic_rho0 = proj.e.EQ_RAD * proj.LambertConfConic_F *
                               pow(t_rlat0, proj.LambertConfConic_N);

  return proj;
}

/*** function to do x,y to lat,long Lambert Conformal Conic projection */

void lamb(const LAMBERT &proj, double lon, double lat, double *x, double *y)
{
  double rho, theta, hold1, hold2, hold3;

  while ((lon - proj.CentralMeridian) < -180.0) lon += 360.0;
  while ((lon - proj.CentralMeridian) > 180.0) lon -= 360.0;
  lat *= D2R;

  hold2 = pow(((1.0 - proj.e.ECC * sin(lat)) / (1.0 + proj.e.ECC * sin(lat))),
              0.5 * proj.e.ECC);
  hold3 = tan(45.0 * D2R - 0.5 * lat);
  if (fabs(hold3) < SMALL) hold3 = 0.0;
  hold1 = (hold3 == 0.0) ? 0.0 : pow(hold3 / hold2, proj.LambertConfConic_N);
  rho   = proj.e.EQ_RAD * proj.LambertConfConic_F * hold1;
  theta = proj.LambertConfConic_N * (lon - proj.CentralMeridian) * D2R;

  *x = rho * sin(theta);
  *y = proj.LambertConfConic_rho0 - rho * cos(theta);
}

/*** function to do lat,long to x,y inverse
                        Lambert Conformal Conic projection */

void ilamb(const LAMBERT &proj, double *lon, double *lat, double x, double y)
{
  int i;
  double theta, temp, rho, t, tphi, phi, delta;

  theta = atan(x / (proj.LambertConfConic_rho0 - y));
  *lon  = (theta / proj.LambertConfConic_N) * R2D + proj.CentralMeridian;

  temp = x * x +
         (proj.LambertConfConic_rho0 - y) * (proj.LambertConfConic_rho0 - y);
  rho   = copysign(d_sqrt(temp), proj.LambertConfConic_N);
  t     = pow((rho / (proj.e.EQ_RAD * proj.LambertConfConic_F)),
          1. / proj.LambertConfConic_N);
  tphi  = 90.0 * D2R - 2.0 * atan(t);
  delta = 1.0;
  for (i = 0; i < 100 && delta > 1.0e-8; i++)
  {
    temp  = (1.0 - proj.e.ECC * sin(tphi)) / (1.0 + proj.e.ECC * sin(tphi));
    phi   = 90.0 * D2R - 2.0 * atan(t * pow(temp, 0.5 * proj.e.ECC));
    delta = fabs(fabs(tphi) - fabs(phi));
    tphi  = phi;
  }
  *lat = phi * R2D;
}

/* Transverse Mercator Projection (TM) */

TRANS_MERCATOR vtm(const char *ellipsoid_name,
                   double lon0,
                   double lat0,
                   bool use_false_easting,
                   long false_easting,
                   double map_scale_factor)
{
  TRANS_MERCATOR proj{};
  proj.e = map_setup_proxy(ellipsoid_name);

  /* Set up an TM projection */
  double e1;

  e1 = (1.0 - d_sqrt(1.0 - proj.e.ECC2)) / (1.0 + d_sqrt(1.0 - proj.e.ECC2));
  proj.t_e2  = proj.e.ECC2 / (1.0 - proj.e.ECC2);
  proj.t_c1  = (1.0 - 0.25 * proj.e.ECC2 - 3.0 * proj.e.ECC4 / 64.0 -
               5.0 * proj.e.ECC6 / 256.0);
  proj.t_c2  = (3.0 * proj.e.ECC2 / 8.0 + 3.0 * proj.e.ECC4 / 32.0 +
               45.0 * proj.e.ECC6 / 1024.0);
  proj.t_c3  = (15.0 * proj.e.ECC4 / 256.0 + 45.0 * proj.e.ECC6 / 1024.0);
  proj.t_c4  = (35.0 * proj.e.ECC6 / 3072.0);
  proj.t_ic1 = (1.5 * e1 - 27.0 * pow(e1, 3.0) / 32.0);
  proj.t_ic2 = (21.0 * e1 * e1 / 16.0 - 55.0 * pow(e1, 4.0) / 32.0);
  proj.t_ic3 = (151.0 * pow(e1, 3.0) / 96.0);
  proj.t_ic4 = (1097.0 * pow(e1, 4.0) / 512.0);
  proj.central_meridian  = lon0;
  proj.use_false_easting = use_false_easting;
  proj.false_easting = false_easting;
  proj.map_scale_factor = map_scale_factor;

  // get y offset of central parallel (use lat0 = 0.0 for standard TM w/o
  // offset)
  double lon, lat, x, y;
  lon = lon0;
  lat = lat0;
  proj.y_central_parralel = 0;
  tm(proj, lon, lat, &x, &y);
  proj.y_central_parralel = y;
  return proj;
}

void tm(
    const TRANS_MERCATOR &proj, double lon, double lat, double *x, double *y)
{
  /* Convert lon/lat to TM x/y */
  double N, T, T2, C, A, M, dlon, tan_lat, cos_lat, A2, A3, A5;

  dlon = lon - proj.central_meridian;
  if (fabs(dlon) > 360.0) dlon += copysign(360.0, -dlon);
  if (fabs(dlon) > 180.0) dlon = copysign(360.0 - fabs(dlon), -dlon);
  lat *= D2R;
  M = proj.e.EQ_RAD * (proj.t_c1 * lat - proj.t_c2 * sin(2.0 * lat) +
                       proj.t_c3 * sin(4.0 * lat) - proj.t_c4 * sin(6.0 * lat));
  if (fabs(lat) == M_PI_2)
  {
    *x = 0.0;
    *y = proj.map_scale_factor * M;
  }
  else
  {
    N       = proj.e.EQ_RAD / d_sqrt(1.0 - proj.e.ECC2 * pow(sin(lat), 2.0));
    tan_lat = tan(lat);
    cos_lat = cos(lat);
    T       = tan_lat * tan_lat;
    T2      = T * T;
    C       = proj.t_e2 * cos_lat * cos_lat;
    A       = dlon * D2R * cos_lat;
    A2      = A * A;
    A3      = A2 * A;
    A5      = A3 * A2;
    *x      = proj.map_scale_factor * N *
         (A + (1.0 - T + C) * (A3 * 0.16666666666666666667) +
          (5.0 - 18.0 * T + T2 + 72.0 * C - 58.0 * proj.t_e2) *
              (A5 * 0.00833333333333333333));
    A3 *= A;
    A5 *= A;
    *y = proj.map_scale_factor *
         (M + N * tan(lat) *
                  (0.5 * A2 +
                   (5.0 - T + 9.0 * C + 4.0 * C * C) *
                       (A3 * 0.04166666666666666667) +
                   (61.0 - 58.0 * T + T2 + 600.0 * C - 330.0 * proj.t_e2) *
                       (A5 * 0.00138888888888888889)));
  }

  // correct for TM x offset of central parallel
  *y -= proj.y_central_parralel;

  // correct for false easting
  if (proj.use_false_easting)
  {
    *x += proj.false_easting;
  }
}

void itm(
    const TRANS_MERCATOR &proj, double *lon, double *lat, double x, double y)
{

  // correct for TM x offset of central parallel
  y += proj.y_central_parralel;

  // correct for false easting
  if (proj.use_false_easting)
  {
    x -= proj.false_easting;
  }

  /* Convert TM x/y to lon/lat */
  double M, mu, phi1, C1, C12, T1, T12, tmp, tmp2, N1, R1, D, D2, D3, D5,
      cos_phi1, tan_phi1;

  M    = y / proj.map_scale_factor;
  mu   = M / (proj.e.EQ_RAD * proj.t_c1);
  phi1 = mu + proj.t_ic1 * sin(2.0 * mu) + proj.t_ic2 * sin(4.0 * mu) +
         proj.t_ic3 * sin(6.0 * mu) + proj.t_ic4 * sin(8.0 * mu);
  cos_phi1 = cos(phi1);
  tan_phi1 = tan(phi1);
  C1       = proj.t_e2 * cos_phi1 * cos_phi1;
  C12      = C1 * C1;
  T1       = tan_phi1 * tan_phi1;
  T12      = T1 * T1;
  tmp      = 1.0 - proj.e.ECC2 * (1.0 - cos_phi1 * cos_phi1);
  tmp2     = d_sqrt(tmp);
  N1       = proj.e.EQ_RAD / tmp2;
  R1       = proj.e.EQ_RAD * (1.0 - proj.e.ECC2) / (tmp * tmp2);
  D        = x / (N1 * proj.map_scale_factor);
  D2       = D * D;
  D3       = D2 * D;
  D5       = D3 * D2;

  *lon = proj.central_meridian +
         R2D *
             (D - (1.0 + 2.0 * T1 + C1) * (D3 * 0.16666666666666666667) +
              (5.0 - 2.0 * C1 + 28.0 * T1 - 3.0 * C12 + 8.0 * proj.t_e2 +
               24.0 * T12) *
                  (D5 * 0.00833333333333333333)) /
             cos(phi1);
  D3 *= D;
  D5 *= D;
  *lat =
      phi1 - (N1 * tan(phi1) / R1) *
                 (0.5 * D2 -
                  (5.0 + 3.0 * T1 + 10.0 * C1 - 4.0 * C12 - 9.0 * proj.t_e2) *
                      (D3 * 0.04166666666666666667) +
                  (61.0 + 90.0 * T1 + 298 * C1 + 45.0 * T12 -
                   252.0 * proj.t_e2 - 3.0 * C12) *
                      (D5 * 0.00138888888888888889));
  (*lat) *= R2D;
}

/* Azimuthal Equidistant Projection (AE) */

AZIMUTHAL_EQUIDIST
vazeqdist(const char *ellipsoid_name, double lon0, double lat0)
{
  AZIMUTHAL_EQUIDIST proj{};
  proj.e = map_setup_proxy(ellipsoid_name);

  /* Set up azimuthal equidistant projection */
  proj.central_meridian = lon0;
  proj.pole             = lat0;
  proj.sinp             = sin(lat0 * D2R);
  proj.cosp             = cos(lat0 * D2R);

  return proj;
}

void azeqdist(const AZIMUTHAL_EQUIDIST &proj,
              double lon,
              double lat,
              double *x,
              double *y)
{
  /* Convert lon/lat to azimuthal equidistant x/y */
  double k, dlon, cc, c, clat, clon, slat;

  while ((lon - proj.central_meridian) < -180.0) lon += 360.0;
  while ((lon - proj.central_meridian) > 180.0) lon -= 360.0;
  dlon = (lon - proj.central_meridian) * D2R;
  lat *= D2R;
  slat = sin(lat);
  clat = cos(lat);
  clon = cos(dlon);

  cc = proj.sinp * slat + proj.cosp * clat * clon;
  if (fabs(cc) >= 1.0)
    *x = *y = 0.0;
  else
  {
    c  = acos(cc);
    k  = proj.e.EQ_RAD * c / sin(c);
    *x = k * clat * sin(dlon);
    *y = k * (proj.cosp * slat - proj.sinp * clat * clon);
  }
}

void iazeqdist(const AZIMUTHAL_EQUIDIST &proj,
               double *lon,
               double *lat,
               double x,
               double y)
{
  /* Convert azimuthal equidistant x/yto lon/lat */
  double rho, c, sin_c, cos_c;

  rho = hypot(x, y);

  if (rho == 0.0)
  {
    *lat = proj.pole;
    *lon = proj.central_meridian;
  }
  else
  {
    c     = rho / proj.e.EQ_RAD;
    sin_c = sin(c);
    cos_c = cos(c);
    *lat  = asin(cos_c * proj.sinp + (y * sin_c * proj.cosp / rho)) * R2D;
    if (proj.pole == 90.0)
      *lon = proj.central_meridian + R2D * atan2(x, -y);
    else if (proj.pole == -90.0)
      *lon = proj.central_meridian + R2D * atan2(x, y);
    else
      *lon = proj.central_meridian +
             R2D * atan2(x * sin_c,
                         (rho * proj.cosp * cos_c - y * proj.sinp * sin_c));
    if ((*lon) <= -180) (*lon) += 360.0;
  }
}

} // namespace GMT
} // namespace HDD
