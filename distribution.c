/*! 
  \file
  Calculate accuracy of trajectories based on reference calculations.
*/

#include "libtrac.h"
#include <gsl/gsl_sort.h>
#define MAX_CLASSES 1000

int main(
  int argc,
  char *argv[]) {

  atm_t *atm1, *atm2;

  char *name, atmStartName[LEN];
  char *year, *mon, *day, *hour, *min;

  ctl_t ctl;

  FILE *out;

  double maxd=0, d, x0[3], x1[3], ahtd, avtd, ahtd2, avtd2, rhtd, rvtd, t;

  double dh, dv, pdfh[MAX_CLASSES], pdfv[MAX_CLASSES], cdfv[MAX_CLASSES], cdfh[MAX_CLASSES];


  static double distv[NP];
  int distiv=0;
  
  static double disth[NP];
  int distih=0;
    
  int ip, i=0;

  int class, classes;

  /* Allocate... */
  ALLOC(atm1, atm_t, 1);
  ALLOC(atm2, atm_t, 1);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG
      ("Give parameters: <outfile> <atm1a> <atm1b> ");

  strcpy(atmStartName, argv[2] );
  if (argv[2][strlen(argv[2]) - 1] == 'c') {
    ctl.atm_iformat = 1;
  }

  /* Write info... */
  printf("Write trajectory analysis data: %s\n", argv[1]);

  /* Create output file... */
  if (!(out = fopen(argv[1], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = horizontal divergence (class center)\n"
	  "# $2 = horizontal CDF\n"
	  "# $3 = horizontal PDF\n"
	  "# $4 = vertical divergence (class center)\n"
	  "# $5 = vertical CDF\n"
	  "# $6 = vertical PDF\n");

  /* Read atmopheric data... */
  read_atm(NULL, argv[2], atm1, &ctl);
  read_atm(NULL, argv[3], atm2, &ctl);

  /* Check if structs match... */
  if (atm1->np != atm2->np)
    ERRMSG("Different numbers of parcels!");
  for (ip = 0; ip < atm1->np; ip++)
    if (fabs(atm1->time[ip] - atm2->time[ip]) > 1)
      printf("WARNING! Times do not match! dt=%f\n",
	       fabs(atm1->time[ip] - atm2->time[ip]));

  /* Init... */
  ahtd = avtd = rhtd = rvtd = 0;
  ahtd2 = avtd2 = 0;
  distih = distiv = 0;

  /* Loop over air parcels... */
  for (ip = 0; ip < atm1->np; ip++) {

      /* Sum up the vertical error... */
      d=fabs(Z(atm1->p[ip]) - Z(atm2->p[ip]));
      avtd +=  d;
      avtd2 +=  d*d;
      distv[distiv++] = d;

      /* Sum up the horizontal error... */
      geo2cart(0, atm1->lon[ip], atm1->lat[ip], x0);
      geo2cart(0, atm2->lon[ip], atm2->lat[ip], x1);
      d = DIST(x0, x1);
      ahtd += d;
      ahtd2 += d*d;
      disth[distih++] = d;
      if(d>maxd){
	maxd=d;
	i=ip;
      }
    }
    
  gsl_sort(distv, (size_t)1, (size_t)distiv);
  gsl_sort(disth, (size_t)1, (size_t)distih);
  
  
  classes = (atm1->np/5);
  if(classes>MAX_CLASSES) classes=MAX_CLASSES;

  dh = disth[distih-1] / classes;
  dv = distv[distiv-1] / classes;

  for (ip = 0; ip < distih; ip++) {
    class = (int)(disth[ip] / dh);
    pdfh[class]++;        

    class = (int)(distv[ip] / dv);
    pdfv[class]++;
  }

  cdfh[0]=pdfh[0];
  cdfv[0]=pdfv[0];
  for (ip=0; ip<classes-1; ip++) {
    cdfh[ip+1]=cdfh[ip]+pdfh[ip+1];
    cdfv[ip+1]=cdfv[ip]+pdfv[ip+1];
  }
  
  d = (double) distih;
  fprintf(out, "%g %g %g %g %g %g\n", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  for (ip=0; ip<classes-1; ip++) {     
     fprintf(out, "%g %g %g %g %g %g\n",
        dh * (ip + 0.5), pdfh[ip] / d, cdfh[ip] / d, 
        dv * (ip + 0.5), pdfv[ip] / d, cdfv[ip] / d);
  }
  
  /* Close file... */
  fclose(out);

  /* Free... */
  free(atm1);
  free(atm2);

  return EXIT_SUCCESS;
}
