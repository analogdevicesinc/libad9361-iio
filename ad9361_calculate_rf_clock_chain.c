
#include "ad9361.h"
#include <stdbool.h>
#include <stdio.h>
#include <errno.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MIN_ADC_CLK 25000000UL /* 25 MHz */
#define MAX_ADC_CLK 640000000UL /* 640 MHz */
#define MIN_DAC_CLK 25000000UL /* 25 MHz */
#define MAX_DAC_CLK	(MAX_ADC_CLK / 2)
#define MAX_BBPLL_FREF 70007000UL /* 70 MHz + 100ppm */
#define MIN_BBPLL_FREQ 714928500UL /* 715 MHz - 100ppm */
#define MAX_BBPLL_FREQ 1430143000UL /* 1430 MHz + 100ppm */
#define MAX_BBPLL_DIV 64
#define MIN_BBPLL_DIV	2
#define MAX_DATA_RATE 61440000UL // Output of FIR (RX)
#define MIN_DATA_RATE MIN_ADC_CLK / 48
// Output of FIR on RX
#define MAX_FIR MAX_DATA_RATE * 2
// Associated with outputs of stage
#define MAX_RX_HB1 122880000UL
#define MAX_RX_HB2 245760000UL
#define MAX_RX_HB3 320000000UL
// Associated with inputs of stage
#define MAX_TX_HB1 122880000UL
#define MAX_TX_HB2 245760000UL
#define MAX_TX_HB3 320000000UL

const unsigned long RX_MAX_PATH_RATES[] = {MAX_ADC_CLK, MAX_RX_HB3, MAX_RX_HB2, MAX_RX_HB1, MAX_FIR};
const unsigned long TX_MAX_PATH_RATES[] = {MAX_DAC_CLK, MAX_TX_HB3, MAX_TX_HB2, MAX_TX_HB1, MAX_FIR};
const unsigned long RX_MIN_PATH_RATES[] = {MIN_ADC_CLK, 0, 0, 0, 0};
const unsigned long TX_MIN_PATH_RATES[] = {MIN_DAC_CLK, 0, 0, 0, 0};

#define check(val,min,max) ( (val)<=(max) ? (val)>=(min) : false )

unsigned long max_rate_found;


bool check_rates(int FIR, const int *HB_configs, unsigned long samp_rate,
                 unsigned long *rates)
{
    int j;
    bool c = true;

    rates[5] = samp_rate;
    for (j=4; j>0; j--) {
        rates[j] = rates[j+1]*HB_configs[j-1];
        if (j>1) {
            c &= check(rates[j],TX_MIN_PATH_RATES[j-1],TX_MAX_PATH_RATES[j-1]);
            c &= check(rates[j],RX_MIN_PATH_RATES[j-1],RX_MAX_PATH_RATES[j-1]);
        }
    }
    return c;
}

int determine_pll_div(unsigned long *rates)
{
    // Determine necessary PLL multiplier
    unsigned long long tmp;
    int PLL_mult = MAX_BBPLL_DIV;
    while (PLL_mult>1) {
        tmp = (unsigned long long) rates[1]*PLL_mult;
        if (check(tmp, MIN_BBPLL_FREQ, MAX_BBPLL_FREQ)) {
            rates[0] = (unsigned long) rates[1]*PLL_mult;
            return PLL_mult;
        }
        PLL_mult >>= 1;
    }
    return -1;
}

int check_dac_adc_config(unsigned long pll_bb, int PLL_mult,
                         int dec_table_index)
{
    // Need to determine if DAC divider is required and if ADC and DAC rates
    // can be satisfied
    unsigned long with_dd, without_dd;
    bool a, b, c;

    with_dd = pll_bb/PLL_mult/2;
    without_dd = pll_bb/PLL_mult/1;

    a = check(with_dd, MIN_DAC_CLK, MAX_DAC_CLK);
    b = check(without_dd, MIN_ADC_CLK, MAX_ADC_CLK);
    c = check(without_dd, MIN_DAC_CLK, MAX_DAC_CLK);

    if (c && b)
        return 1; // Run without dac div
    else if (a && b && (dec_table_index<5))
        return 2; // Run with dac div
    else
        return -1;
}

void set_rates(unsigned long *rx_path_clks,
               unsigned long *tx_path_clks, int DAC_div, unsigned long *rates,
               int dec_table_index)
{
    int k;

    // Check if ADC will run faster in config
    if (rates[1]>max_rate_found)
        max_rate_found = rates[1];
    else
        return;

    for (k=0; k<6; k++) {
        rx_path_clks[k] = rates[k];
        tx_path_clks[k] = rates[k];

        if (k>0) { // Adjust HB's for DAC divider setting
            if ((dec_table_index<2) && (k<4))
                tx_path_clks[k] = rates[k]/DAC_div;
            else if ((dec_table_index<4) && (k<3))
                tx_path_clks[k] = rates[k]/DAC_div;
        }
    }
}


int determine_path_rates_with_fir(unsigned long sample_rate,
                                  unsigned long rate_gov,
                                  unsigned long *rx_path_clks,
                                  unsigned long *tx_path_clks,
                                  int FIR)
{
    unsigned long rates[6];
    int PLL_mult, k;

    max_rate_found = 0UL;

    const int HB_configs[][4] = {
        {3,2,2,FIR}, //12
        {2,2,2,FIR}, //8
        {3,2,1,FIR}, //6
        {2,2,1,FIR}, //4
        {2,1,1,FIR}, //2
        {3,1,1,FIR}, //3
        {1,1,1,FIR}, //1
    };

    // RX Path:
    // BBPLL -> /PLL_div -> /HB3 -> /HB2 -> /HB1 -> /FIR
    // TX Path:
    // BBPLL -> /(PLL_div*DAC_div) -> /HB3 -> /HB2 -> /HB1 -> /FIR

    // Cycle through possible decimations from highest to lowest
    for (k=0; k<7; k++) {
        // HB3 cannot be 3 if rate_gov enabled
        if ((rate_gov>0) && HB_configs[k][0]==3)
            continue;
        // Check if HB and FIR rates are valid
        if ( check_rates(FIR, HB_configs[k], sample_rate, rates) ) {
            // Calculate PLL divider for configuration
            PLL_mult = determine_pll_div(rates);
            if (PLL_mult>0) {
                // Determine DAC divider setting and check ADC/DAC settings
                int dac_div = check_dac_adc_config(rates[0], PLL_mult, k);
                // printf("dac_div: %d\n",dac_div);
                if (dac_div>0)
                    set_rates(rx_path_clks, tx_path_clks, dac_div, rates, k);
            }
        }
    }

    if (max_rate_found==0UL)
        return -EINVAL;
    else
        return 0;
}

int ad9361_calculate_rf_clock_chain(unsigned long sample_rate,
                                    unsigned long rate_gov,
                                    unsigned long *rx_path_clks,
                                    unsigned long *tx_path_clks)
{
    int ret, k;
    int FIR[] = {4,2,1};

    // Check desired rate within bounds
    if (!check(sample_rate, MIN_DATA_RATE, MAX_DATA_RATE))
        return -EINVAL;

    // Rate selection will try to:
    // 1. Utilize the maximum decimation in the FIR
    // 2. Run the ADC/DAC as fast as possible
    // 3. Use the most decimation possible starting with HB3(closest to ADC)->HB1

    // Cycle through available FIR settings
    for (k=0; k<3; k++) {
        ret = determine_path_rates_with_fir(sample_rate, rate_gov, rx_path_clks,
                                            tx_path_clks, FIR[k]);
        if (ret==0)
            break;
    }
    return ret;
}
