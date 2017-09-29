

#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include "ad9361.h"
#include <stdlib.h>
#include <stdint.h>

int check_result(short *taps)
{
  FILE* fp;
  char buffer[255];

  fp = fopen("bb_test_taps.txt", "r");

  uint8_t k = 0;

  while(fgets(buffer, 255, (FILE*) fp)) {
      int tap = atoi(buffer);
      printf("Tap %d |%i|%i\n", k+1, tap, taps[k]);
      if (abs(tap - taps[k])>2)
      {
        fclose(fp);
        return 2;
      }
      k++;
  }

  fclose(fp);
  return 0;
}

int main(void)
{
  struct filter_design_parameters fdp;
  int ret = ad9361_filter_config_from_rate(&fdp, 1000000, true);
  if (ret<0)
    exit(1);
  printf("Here\n");
  // Initialize taps
  short outputTaps[128];
  int num_taps;
  // Call designer
  ad9361_generate_fir_taps(&fdp, outputTaps, &num_taps);
  // Check
  return check_result(outputTaps);
}
