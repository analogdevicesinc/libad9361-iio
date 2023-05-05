

#include "ad9361.h"
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main(void)
{
    int ret, k;
    unsigned long Fpass, Fstop, wnomTX, wnomRX;

    unsigned long rates[] = {1000000, 10000000, 20000000, 60000000};

    unsigned long tx_path_clks[6];
    unsigned long rx_path_clks[6];
    unsigned int DAC_div;
    unsigned int tx_gain;
    unsigned int rx_gain;
    short tx_fir[128];
    short rx_fir[128];
    int tx_fir_size;
    int rx_fir_size;

    for (k = 0; k < 4; k++) {

        printf("Testing rate: %lu\n",rates[k]);

        Fpass = rates[k] / 3.0;
        Fstop = Fpass * 1.25;
        wnomTX = 1.6 * Fstop;
        wnomRX = 1.4 * Fstop;

        ret = ad9361_set_bb_rate_custom_filter_manual_bm(rates[k], Fpass, Fstop,
                wnomTX, wnomRX, tx_path_clks, rx_path_clks, &DAC_div, &tx_gain, 
                &rx_gain, tx_fir, rx_fir, &tx_fir_size, &rx_fir_size);
        if (ret<0)
            return ret;

        printf("TX FIR size: %i\n",tx_fir_size);
        printf("RX FIR size: %i\n",rx_fir_size);
        for (int i=0; i<tx_fir_size; i++)
            printf("%i,",tx_fir[i]);
        printf("\n");
        for (int i=0; i<rx_fir_size; i++)
            printf("%i,",rx_fir[i]);
        printf("\n");

        printf("TX path clocks: ");
        for (int i=0; i<6; i++)
            printf("%lu,",tx_path_clks[i]);
        printf("\n");
        printf("RX path clocks: ");
        for (int i=0; i<6; i++)
            printf("%lu,",rx_path_clks[i]);
        printf("\n");

        printf("DAC div: %i\n",DAC_div);
        printf("TX gain: %i\n",tx_gain);
        printf("RX gain: %i\n",rx_gain);

    }
    return 0;
}
