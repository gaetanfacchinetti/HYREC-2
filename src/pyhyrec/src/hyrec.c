/*************************************************************************/
/*                    HYREC-2 MAIN FUNCTION                              */
/*************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "history.h"

int main(void) {
  remove("output_xe.dat");
  
  //HYREC_DATA rec_data;
  HYREC_DATA * data = malloc(sizeof(*data));

  double zmax = 8000.;
  double zmin = 0.;

  data->path_to_hyrec = "";
  hyrec_allocate(data, zmax, zmin);

  rec_get_cosmoparam(stdin, stderr, data->cosmo);
  hyrec_compute(data, MODEL);

  if (data->error == 1) fprintf(stderr,"%s\n",data->error_message);
  else {
    double z = zmax;
    char file[200];
    FILE *fout;
    fout = fopen("output_xe.dat", "a");
    while (z > zmin) {
      fprintf(fout, "%f %1.10E %1.10E\n",z,hyrec_xe(z, data),hyrec_Tm(z, data));
      z -= 1.;
      }
    fclose(fout);
  }
  
  hyrec_free(data);

  return 0;

}
