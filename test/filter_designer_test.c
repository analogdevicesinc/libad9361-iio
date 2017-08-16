

#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include "ad9361.h"
#include <stdlib.h>

int check_result(short *taps)
{
  FILE* fp;
  char buffer[255];

  fp = fopen("correct_taps.txt", "r");

  uint k = 0;
  while(fgets(buffer, 255, (FILE*) fp)) {
      int tap = atoi(buffer);
      printf("|%i|%i|\n", tap, taps[k]);
      if (tap != taps[k])
      {
        fclose(fp);
        return 1;
      }
      k++;
  }

  fclose(fp);
  return 0;
}

int main(void)
{
  struct filter_design_parameters fdp;

  // Initialize input arguments
  fdp.Rdata = 7680000;
  fdp.Fpass = 2250000;
  fdp.Fstop = 2750000;
  fdp.caldiv = 32;
  fdp.FIR = 2;
  fdp.HB1 = 2;
  fdp.DAC_div = 1;

  fdp.Type = "Lowpass";
  fdp.RxTx = "Tx";

  fdp.RFbw = 4236204;
  fdp.converter_rate = 122880000;
  fdp.PLL_rate = 983040000;
  fdp.Fcenter = 0;
  fdp.wnom = 4400000;
  fdp.FIRdBmin = 0;
  fdp.int_FIR = 1;
  fdp.PLL_mult = 8;
  fdp.Apass = 0.1250;
  fdp.Astop = 85;
  fdp.phEQ = -1;
  fdp.HB2 = 2;
  fdp.HB3 = 3;

  // Initialize taps
  short outputTaps[128];

  // Call designer
  ad9361_generate_fir_taps(&fdp, outputTaps);
  printf("Generated Taps\n");

  // Check
  return check_result(outputTaps);

}
