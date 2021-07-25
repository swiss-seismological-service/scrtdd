
int map_setup_proxy(int n_proj, char* ellipsoid_name);

// Lambert Conformal Conic projection (e.g. LAMBERT)
int vlamb(int n_proj, double rlong0, double rlat0, double pha, double phb);
int lamb(int n_proj, double lon, double lat, double *x, double *y);
int ilamb(int n_proj, double *lon, double *lat, double x, double y);

// Transverse Mercator Projection (TM) (e.g. TRANS_MERC) 
// Changed function names (vtm -> vtvm, tm -> tvm, itm -> itvm) to avoid ambiguity with "tm"
void vtvm(int n_proj, double lon0, double lat0, int use_false_easting);
void tvm(int n_proj, double lon, double lat, double *x, double *y);
void itvm(int n_proj, double *lon, double *lat, double x, double y);

// Universal Transverse Mercator Projection (UTM) (not used..?)
void vutm(int n_proj, double lon0, int lat0, int use_false_easting);
void utm(int n_proj, double lon, double lat, double *x, double *y);
void iutm(int n_proj, double *lon, double *lat, double x, double y);

// Azimuthal Equidistant Projection (AE) (e.g. AZIMUTHAL_EQUIDIST)
void vazeqdist(int n_proj, double lon0, double lat0);
void azeqdist(int n_proj, double lon, double lat, double *x, double *y);
void iazeqdist(int n_proj, double *lon, double *lat, double x, double y);
