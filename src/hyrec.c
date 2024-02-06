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

  /*
  rec_get_cosmoparam(stdin, stderr, ptr_data->cosmo);
  hyrec_compute(ptr_data, MODEL);
  if (ptr_data->error == 1) fprintf(stderr,"%s\n",ptr_data->error_message);
  else {
    double z = zmax;
    while (z > zmin) {
      printf("%f %1.10E %1.10E\n",z,hyrec_xe(z, ptr_data),hyrec_Tm(z, ptr_data));
      z -= 1.;
      }
  }
  */
  
  hyrec_free(data);
  free(data);

  return 0;

}
