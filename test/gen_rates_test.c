


#include "ad9361.h"
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>

/******************************Linux Reference*********************************/

struct PHY {
    bool bypass_rx_fir;
    int rx_fir_dec;
    bool bypass_tx_fir;
    int tx_fir_int;
//int tx_intdec;
    bool rx_eq_2tx;
};

#define MIN_ADC_CLK			25000000UL /* 25 MHz */
//#define MIN_ADC_CLK			(MIN_BBPLL_FREQ / MAX_BBPLL_DIV) /* 11.17MHz */
#define MAX_ADC_CLK			640000000UL /* 640 MHz */
#define MAX_DAC_CLK	(MAX_ADC_CLK / 2)

#define MAX_BBPLL_FREF			70007000UL /* 70 MHz + 100ppm */
#define MIN_BBPLL_FREQ			714928500UL /* 715 MHz - 100ppm */
#define MAX_BBPLL_FREQ			1430143000UL /* 1430 MHz + 100ppm */
#define MAX_BBPLL_DIV			64
#define MIN_BBPLL_DIV	2

enum ad9361_pdata_rx_freq {
    BBPLL_FREQ,
    ADC_FREQ,
    R2_FREQ,
    R1_FREQ,
    CLKRF_FREQ,
    RX_SAMPL_FREQ,
    NUM_RX_CLOCKS,
};

enum ad9361_pdata_tx_freq {
    IGNORE,
    DAC_FREQ,
    T2_FREQ,
    T1_FREQ,
    CLKTF_FREQ,
    TX_SAMPL_FREQ,
    NUM_TX_CLOCKS,
};

int linux_calculate_rf_clock_chain(
    unsigned long tx_sample_rate,
    uint32_t rate_gov,
    unsigned long *rx_path_clks,
    unsigned long *tx_path_clks)
{
    // PHY SETTINGS
    struct PHY phy;
    phy.bypass_rx_fir = false;
    phy.bypass_tx_fir = false;
    phy.rx_fir_dec = 4;
    phy.tx_fir_int = 4;
    //phy.tx_intdec = 4;
    phy.rx_eq_2tx = false;

    unsigned long clktf, clkrf, adc_rate = 0, dac_rate = 0;
    uint64_t bbpll_rate;
    int i, index_rx = -1, index_tx = -1, tmp;
    uint32_t div, tx_intdec, rx_intdec, recursion = 1;
    const char clk_dividers[][4] = {
        {12,3,2,2},
        {8,2,2,2},
        {6,3,2,1},//Updated to give correct orientation
        {4,2,2,1},
        {3,3,1,1},
        {2,2,1,1},
        {1,1,1,1},
    };


    if (phy.bypass_rx_fir)
        rx_intdec = 1;
    else
        rx_intdec = phy.rx_fir_dec;

    if (phy.bypass_tx_fir)
        tx_intdec = 1;
    else
        tx_intdec = phy.tx_fir_int;

    if ((rate_gov == 1) && ((rx_intdec * tx_sample_rate * 8) < MIN_ADC_CLK)) {
        recursion = 0;
        rate_gov = 0;
    }

    // dev_dbg(&phy.spi.dev, "%s: requested rate %lu TXFIR int %d RXFIR dec %d mode %s",
    // 	__func__, tx_sample_rate, tx_intdec, rx_intdec,
    // 	rate_gov ? "Nominal" : "Highest OSR");

    if (tx_sample_rate > 61440000UL)
        return -EINVAL;

    clktf = tx_sample_rate * tx_intdec;
    clkrf = tx_sample_rate * rx_intdec * (phy.rx_eq_2tx ? 2 : 1);

    for (i = rate_gov; i < 7; i++) {
        adc_rate = clkrf * clk_dividers[i][0];
        dac_rate = clktf * clk_dividers[i][0];

        if ((adc_rate <= MAX_ADC_CLK) && (adc_rate >= MIN_ADC_CLK)) {


            if (dac_rate > adc_rate)
                tmp = (dac_rate / adc_rate) * -1;
            else
                tmp = adc_rate / dac_rate;

            if (adc_rate <= MAX_DAC_CLK) {
                index_rx = i;
                index_tx = i - ((tmp == 1) ? 0 : tmp);
                dac_rate = adc_rate; /* ADC_CLK */
                break;
            } else {
                dac_rate = adc_rate / 2;  /* ADC_CLK/2 */
                index_rx = i;

                if (i == 4 && tmp >= 0)
                    index_tx = 7; /* STOP: 3/2 != 1 */
                else
                    index_tx = i + ((i == 5 && tmp >= 0) ? 1 : 2) -
                               ((tmp == 1) ? 0 : tmp);

                break;
            }
        }
    }

    if ((index_tx < 0 || index_tx > 6 || index_rx < 0 || index_rx > 6)
        && rate_gov < 7 && recursion) {
        return linux_calculate_rf_clock_chain(tx_sample_rate,
                                              ++rate_gov, rx_path_clks, tx_path_clks);
    } else if ((index_tx < 0 || index_tx > 6 || index_rx < 0 || index_rx > 6)) {
        // dev_err(&phy.spi.dev, "%s: Failed to find suitable dividers: %s",
        // __func__, (adc_rate < MIN_ADC_CLK) ? "ADC clock below limit" : "BBPLL rate above limit");
        printf("Failed to find clock rate\n");
        return -EINVAL;
    }

    /* Calculate target BBPLL rate */
    div = MAX_BBPLL_DIV;

    do {
        bbpll_rate = (uint64_t)adc_rate * div;
        div >>= 1;

    } while ((bbpll_rate > MAX_BBPLL_FREQ) && (div >= MIN_BBPLL_DIV));
    // printf("div: %d\n",div<<=1);

    rx_path_clks[BBPLL_FREQ] = bbpll_rate;
    rx_path_clks[ADC_FREQ] = adc_rate;
    rx_path_clks[R2_FREQ] = rx_path_clks[ADC_FREQ] / clk_dividers[index_rx][1];
    rx_path_clks[R1_FREQ] = rx_path_clks[R2_FREQ] / clk_dividers[index_rx][2];
    rx_path_clks[CLKRF_FREQ] = rx_path_clks[R1_FREQ] / clk_dividers[index_rx][3];
    rx_path_clks[RX_SAMPL_FREQ] = rx_path_clks[CLKRF_FREQ] / 	rx_intdec;

    tx_path_clks[BBPLL_FREQ] = bbpll_rate;
    tx_path_clks[DAC_FREQ] = dac_rate;
    tx_path_clks[T2_FREQ] = tx_path_clks[DAC_FREQ] / clk_dividers[index_tx][1];
    tx_path_clks[T1_FREQ] =tx_path_clks[T2_FREQ] / clk_dividers[index_tx][2];
    tx_path_clks[CLKTF_FREQ] = tx_path_clks[T1_FREQ] / clk_dividers[index_tx][3];
    tx_path_clks[TX_SAMPL_FREQ] = tx_path_clks[CLKTF_FREQ] / 	tx_intdec;

    return 0;
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

bool check_result(unsigned long *rx1, unsigned long *rx2, unsigned long *tx1,
                  unsigned long *tx2)
{
    bool r = false;
    int o=0;
    for (o=0; o<6; o++)
        r |= ( (rx1[o] != rx2[o]) || (tx1[o] != tx2[o]) ) ;

    if (r) {
        printf("LINUX\n");
        printf("RX %lu | TX %lu\n",rx1[o],tx1[o]);
        for (o=1; o<6; o++)
            printf("RX %lu (%lu)| TX %lu (%lu)\n",rx1[o],rx1[o-1]/rx1[o],tx1[o],
                   tx1[o-1]/tx1[o]);
        printf("----------\n");
        printf("LIBAD9361\n");
        printf("RX %lu | TX %lu\n",rx2[0],tx2[0]);
        for (o=1; o<6; o++)
            printf("RX %lu (%lu)| TX %lu (%lu)\n",rx2[o],rx2[o-1]/rx2[o],tx2[o],
                   tx2[o-1]/tx2[o]);
    }
    return r;
}

int main(void)
{
    int ret, k;
    unsigned long rx1[6], tx1[6];
    unsigned long rx2[6], tx2[6];
    unsigned long samples_rates;

    for (k=1; k<35000; k++) {

        samples_rates = 1000*k + 520000UL;

        uint32_t rate_governor = 0; //phy->rate_governor ? 1500000U : 1000000U;

        // Baseline from linux driver
        ret = linux_calculate_rf_clock_chain(samples_rates, rate_governor, rx1, tx1);
        if (ret<0)
            return ret;

        ret = ad9361_calculate_rf_clock_chain(samples_rates, rate_governor, rx2, tx2);
        if (ret<0)
            return ret;

        if (check_result(rx1,rx2,tx1,tx2))
            return -1;
    }
    return 0;
}
