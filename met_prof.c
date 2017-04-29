/*! 
  \file
  Extract vertical profile from meteorological data.
*/

#include "libtrac.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/* Maximum number of altitudes. */
#define NZ 1000

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  met_t *met;

  FILE *in, *out;

  static double timem[NZ], z, z0, z1, dz, lon, lon0, lon1, dlon, lonm[NZ],
    lat, lat0, lat1, dlat, latm[NZ], t, tm[NZ],
    u, um[NZ], v, vm[NZ], w, wm[NZ];

  static int epol, i, iz, np[NZ], redx, redy, redp;

  /* Allocate... */
  ALLOC(met, met_t, 1);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <prof.tab> <met0> [ <met1> ... ]");

  /* Read control parameters... */
  read_ctl(NULL, NULL, argc, argv, &ctl);
  z0 = scan_ctl(NULL, NULL, argc, argv, "Z0", -1, "0", NULL);
  z1 = scan_ctl(NULL, NULL, argc, argv, "Z1", -1, "60", NULL);
  dz = scan_ctl(NULL, NULL, argc, argv, "DZ", -1, "1", NULL);
  lon0 = scan_ctl(NULL, NULL, argc, argv, "LON0", -1, "0", NULL);
  lon1 = scan_ctl(NULL, NULL, argc, argv, "LON1", -1, "0", NULL);
  dlon = scan_ctl(NULL, NULL, argc, argv, "DLON", -1, "1", NULL);
  lat0 = scan_ctl(NULL, NULL, argc, argv, "LAT0", -1, "0", NULL);
  lat1 = scan_ctl(NULL, NULL, argc, argv, "LAT1", -1, "0", NULL);
  dlat = scan_ctl(NULL, NULL, argc, argv, "DLAT", -1, "1", NULL);
  epol = (int) scan_ctl(NULL, NULL, argc, argv, "EPOL", -1, "0", NULL);
  redx = (int) scan_ctl(NULL, NULL, argc, argv, "REDX", -1, "1", NULL);
  redy = (int) scan_ctl(NULL, NULL, argc, argv, "REDY", -1, "1", NULL);
  redp = (int) scan_ctl(NULL, NULL, argc, argv, "REDP", -1, "1", NULL);

  /* Loop over input files... */
  for (i = 3; i < argc; i++) {

    /* Read meteorological data... */
    if (!(in = fopen(argv[i], "r")))
      continue;
    else
      fclose(in);
    read_met(argv[i], met);

    /* Extrapolate... */
    if (epol)
      extrapolate_met(met);

    /* Reduce resolution... */
    reduce_met(met, redx, redy, redp);

    /* Average... */
    for (z = z0; z <= z1; z += dz) {
      iz = (int) ((z - z0) / dz);
      if (iz < 0 || iz > NZ)
	ERRMSG("Too many altitudes!");
      for (lon = lon0; lon <= lon1; lon += dlon)
	for (lat = lat0; lat <= lat1; lat += dlat) {
	  intpol_met_space(met, P(z), lon, lat, &t, &u, &v, &w);
	  timem[iz] += met->time;
	  lonm[iz] += lon;
	  latm[iz] += lat;
	  tm[iz] += t;
	  um[iz] += u;
	  vm[iz] += v;
	  wm[iz] += w;
	  np[iz]++;
	}
    }
  }

  /* Create output file... */
  printf("Write meteorological data file: %s\n", argv[2]);
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1  = time [s]\n"
	  "# $2  = altitude [km]\n"
	  "# $3  = longitude [deg]\n"
	  "# $4  = latitude [deg]\n"
	  "# $5  = pressure [hPa]\n"
	  "# $6  = temperature [K]\n"
	  "# $7  = zonal wind [m/s]\n"
	  "# $8  = meridional wind [m/s]\n"
	  "# $9  = vertical wind [hPa/s]\n\n");

  /* Write data... */
  for (z = z0; z <= z1; z += dz) {
    iz = (int) ((z - z0) / dz);
    fprintf(out, "%.2f %g %g %g %g %g %g %g %g\n",
	    timem[iz] / np[iz], z, lonm[iz] / np[iz], latm[iz] / np[iz], P(z),
	    tm[iz] / np[iz], um[iz] / np[iz], vm[iz] / np[iz],
	    wm[iz] / np[iz]);
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(met);

  return EXIT_SUCCESS;
}
