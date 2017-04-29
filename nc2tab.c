/*! 
  \file
  Convert atmospheric data files between ASCII and netCDF.
*/

#include "libtrac.h"

int main(
  int argc,
  char *argv[]) {

  static atm_t atm;
  static ctl_t ctl;
  char *dirname=".";

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <ctl> <atm_in> <atm_out> <grid_out>");

  /* Read control parameter */
  read_ctl(dirname, argv[1], argc, argv, &ctl);

  /* Read atmospheric data... */
  read_atm(dirname, argv[2], &atm, &ctl);

  /* Write atmospheric data... */
  if (argv[3][0] != '-')
    write_atm(dirname, argv[3], &atm, &ctl);

  /* Write gridded data... */
  if (argv[4][0] != '-')
    write_grid(dirname, argv[4], &atm, &ctl, 0, 1e100);

  return EXIT_SUCCESS;
}
