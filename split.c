/*! 
  \file
  Split air parcels into a larger number of parcels.
*/

#include "libtrac.h"

int main(
  int argc,
  char *argv[]) {

  static atm_t *atm, *atm2;

  static ctl_t ctl;

  gsl_rng *rng;

  double m, mtot = 0, dx, dz, mmax = 0, z0, z1, lon0, lon1, lat0, lat1;

  int i, ip, iq, n;

  /* Allocate... */
  ALLOC(atm, atm_t, 1);
  ALLOC(atm2, atm_t, 1);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <atm_in> <atm_out>");

  /* Read control parameters... */
  read_ctl(NULL, NULL, argc, argv, &ctl);
  n = (int) scan_ctl(NULL, NULL, argc, argv, "N", -1, "", NULL);
  m = scan_ctl(NULL, NULL, argc, argv, "M", -1, "-999", NULL);
  dz = scan_ctl(NULL, NULL, argc, argv, "DZ", -1, "0", NULL);
  z0 = scan_ctl(NULL, NULL, argc, argv, "Z0", -1, "0", NULL);
  z1 = scan_ctl(NULL, NULL, argc, argv, "Z1", -1, "0", NULL);
  dx = scan_ctl(NULL, NULL, argc, argv, "DX", -1, "0", NULL);
  lon0 = scan_ctl(NULL, NULL, argc, argv, "LON0", -1, "0", NULL);
  lon1 = scan_ctl(NULL, NULL, argc, argv, "LON1", -1, "0", NULL);
  lat0 = scan_ctl(NULL, NULL, argc, argv, "LAT0", -1, "0", NULL);
  lat1 = scan_ctl(NULL, NULL, argc, argv, "LAT1", -1, "0", NULL);

  /* Init random number generator... */
  gsl_rng_env_setup();
  rng = gsl_rng_alloc(gsl_rng_default);

  /* Read atmospheric data... */
  read_atm(NULL, argv[2], atm, &ctl);

  /* Get total and maximum mass... */
  if (ctl.qnt_mass >= 0)
    for (ip = 0; ip < atm->np; ip++) {
      mtot += atm->q[ctl.qnt_mass][ip];
      mmax = GSL_MAX(mmax, atm->q[ctl.qnt_mass][ip]);
    }
  if (m > 0)
    mtot = m;

  /* Loop over air parcels... */
  for (i = 0; i < n; i++) {

    /* Select air parcel... */
    if (ctl.qnt_mass >= 0)
      do {
	ip = (int) gsl_rng_uniform_int(rng, (long unsigned int) atm->np);
      } while (gsl_rng_uniform(rng) > atm->q[ctl.qnt_mass][ip] / mmax);
    else
      ip = (int) gsl_rng_uniform_int(rng, (long unsigned int) atm->np);

    /* Set time... */
    atm2->time[atm2->np] = atm->time[ip];

    /* Set vertical position... */
    if (z1 > z0)
      atm2->p[atm2->np] = P(z0 + (z1 - z0) * gsl_rng_uniform_pos(rng));
    else
      atm2->p[atm2->np] = atm->p[ip]
	+ dz2dp(gsl_ran_gaussian_ziggurat(rng, dz / 2.3548), atm->p[ip]);

    /* Set horizontal position... */
    if (lon1 > lon0 && lat1 > lat0) {
      atm2->lon[atm2->np] = lon0 + (lon1 - lon0) * gsl_rng_uniform_pos(rng);
      atm2->lat[atm2->np] = lat0 + (lat1 - lat0) * gsl_rng_uniform_pos(rng);
    } else {
      atm2->lon[atm2->np] = atm->lon[ip]
	+ gsl_ran_gaussian_ziggurat(rng, dx2deg(dx, atm->lat[ip]) / 2.3548);
      atm2->lat[atm2->np] = atm->lat[ip]
	+ gsl_ran_gaussian_ziggurat(rng, dy2deg(dx) / 2.3548);
    }

    /* Copy quantities... */
    for (iq = 0; iq < ctl.nq; iq++)
      atm2->q[iq][atm2->np] = atm->q[iq][ip];

    /* Adjust mass... */
    if (ctl.qnt_mass >= 0)
      atm2->q[ctl.qnt_mass][atm2->np] = mtot / n;

    /* Increment particle counter... */
    if ((++atm2->np) >= NP)
      ERRMSG("Too many air parcels!");
  }

  /* Save data and close file... */
  write_atm(NULL, argv[3], atm2, &ctl);

  /* Free... */
  free(atm);
  free(atm2);

  return EXIT_SUCCESS;
}
