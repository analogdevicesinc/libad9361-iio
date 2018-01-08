

#include "ad9361.h"
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int check_result(short *taps, char tap_filename[100])
{
    FILE *fp;
    char buffer[255];
    printf("FN: %s\n",tap_filename);
    fp = fopen(tap_filename, "r");

    uint8_t k = 0;

    while (fgets(buffer, 255, (FILE *)fp)) {
        int tap = atoi(buffer);
        // printf("|%i|%i|\n", tap, taps[k]);
        // Taps need to be within 2 samples (remez had randomness in it)
        if (abs(tap - taps[k])>2) {
            fclose(fp);
            printf("Error: Tap missmatch\n");
            return -2;
        }
        k++;
    }
    fclose(fp);
    return 0;
}

int main(void)
{
    struct filter_design_parameters fdpTX;
    struct filter_design_parameters fdpRX;
    short outputTaps[128];
    int num_taps, ret, gain, k;
    char filename[100];

    unsigned long rates[] = {1000000, 10000000, 20000000, 60000000};

    for (k = 0; k < 4; k++) {

        // Generate rates in filter struct
        ret = ad9361_calculate_rf_clock_chain_fdp(&fdpTX, &fdpRX, rates[k]);

        // Test RX side
        ret = ad9361_generate_fir_taps(&fdpRX, outputTaps, &num_taps, &gain);
        if (ret < 0)
            return ret;
        sprintf(filename,"rateRX_%lu.taps",rates[k]);
        ret = check_result(outputTaps,filename);
        if (ret < 0)
            return ret;

        // Test TX side
        ret = ad9361_generate_fir_taps(&fdpTX, outputTaps, &num_taps, &gain);
        if (ret < 0)
            return ret;
        sprintf(filename,"rateTX_%lu.taps",rates[k]);
        ret = check_result(outputTaps,filename);
        if (ret < 0)
            return ret;
    }
    return 0;
}
