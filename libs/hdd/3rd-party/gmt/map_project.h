#ifndef __HDD_GMT_H__
#define __HDD_GMT_H__

namespace HDD {
namespace GMT {

struct Ellip
{
  double EQ_RAD;
  double ECC;
  double ECC2;
  double ECC4;
  double ECC6;
};

// Lambert Conformal Conic projection (e.g. LAMBERT)

struct LAMBERT
{
  Ellip e;
  double pha;
  double phb;
  bool NorthPole; /* TRUE if projection is on northern hermisphere, FALSE on
                     southern */
  double CentralMeridian; /* Central meridian for projection */
  double Pole;            /* +90 pr -90, depending on which pole */

  /* Lambert conformal conic parameters.
                  (See Snyder for details on all parameters) */
  double LambertConfConic_N;
  double LambertConfConic_F;
  double LambertConfConic_rho0;
};

LAMBERT vlamb(const char *ellipsoid_name,
              double rlong0,
              double rlat0,
              double pha,
              double phb);
void lamb(const LAMBERT &proj, double lon, double lat, double *x, double *y);
void ilamb(const LAMBERT &proj, double *lon, double *lat, double x, double y);

// Transverse Mercator Projection (TM) (e.g. TRANS_MERC)

struct TRANS_MERCATOR
{
  Ellip e;
  bool use_false_easting;  // flag to apply false easting (500km added to X when
                           // converting geog->UTM, and v.v.)
  long false_easting;  // false easting to use (e.g. many UTM systems use 500km added to X when converting geog->UTM, and v.v.)
  double map_scale_factor;
  double central_meridian; // Central meridian for projection
  double y_central_parralel; // y offset of central parallel
  double t_e2;
  double t_c1, t_c2, t_c3, t_c4;
  double t_ic1, t_ic2, t_ic3, t_ic4;
};

TRANS_MERCATOR vtm(const char *ellipsoid_name,
                   double lon0,
                   double lat0,
                   bool use_false_easting,
                   long false_easting,
                   double map_scale_factor);
void tm(
    const TRANS_MERCATOR &proj, double lon, double lat, double *x, double *y);
void itm(
    const TRANS_MERCATOR &proj, double *lon, double *lat, double x, double y);

// Azimuthal Equidistant Projection (AE) (e.g. AZIMUTHAL_EQUIDIST)

struct AZIMUTHAL_EQUIDIST
{
  Ellip e;
  bool north_pole; /* TRUE if projection is on northern hemisphere, FALSE on
                      southern */

  double central_meridian; // Central meridian (longitude) for projection
  double pole;             // Central latitude for projection
  double sinp;
  double cosp;
};

AZIMUTHAL_EQUIDIST
vazeqdist(const char *ellipsoid_name, double lon0, double lat0);

void azeqdist(const AZIMUTHAL_EQUIDIST &proj,
              double lon,
              double lat,
              double *x,
              double *y);

void iazeqdist(const AZIMUTHAL_EQUIDIST &proj,
               double *lon,
               double *lat,
               double x,
               double y);
} // namespace GMT
} // namespace HDD

#endif
