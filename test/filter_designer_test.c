

#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include "ad9361.h"
#include <stdlib.h>

int check_result(short *taps)
{
  FILE * fp;
  char * line = NULL;
  size_t len = 0;
  ssize_t read;
  fp = fopen("correct_taps.txt", "r");
  if (fp == NULL)
    exit(EXIT_FAILURE);

  uint k = 0;
  while ((read = getline(&line, &len, fp)) != -1) {
    int tap = atoi(line);
    printf("%i %i\n", tap, taps[k]);
    if (tap != taps[k])
    {
      fclose(fp);
      if (line)
        free(line);
      return 1;
    }
    k++;
  }
  fclose(fp);
  if (line)
    free(line);
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
