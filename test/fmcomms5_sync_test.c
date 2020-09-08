
#include "ad9361.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#ifdef __APPLE__
#include <iio/iio.h>
#else
#include <iio.h>
#endif
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#include <time.h>
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MHZ(x) ((long long)(x*1000000.0 + .5))
#define GHZ(x) ((long long)(x*1000000000.0 + .5))

#define DEV_RX_NAME "cf-ad9361-A"
#define DEV_RX_SLAVE_NAME "cf-ad9361-B"
#define DEV_TX_NAME "cf-ad9361-dds-core-lpc"
#define DEV_TX_SLAVE_NAME "cf-ad9361-dds-core-B"
#define DEV_PHY_NAME "ad9361-phy"
#define DEV_PHY_SLAVE_NAME "ad9361-phy-B"

#define SAMPLES 32768
#define M_2PI 2 * M_PI
#define STALE_BUFFERS 20

#define ENABLE_PERFORMANCE_TESTS 0

bool same_chip = 0;

#define CHECK(expr)                                                            \
  if (expr < 0) {                                                              \
    return expr;                                                               \
  }

static struct iio_device *dev_phy, *dev_phy_slave;
static struct iio_device *dev_rx, *dev_rx_slave;
static struct iio_device *dev_tx, *dev_tx_slave;
static struct iio_buffer *rxbuf;
static struct iio_channel *rxa_chan_real, *rxa_chan_imag;
static struct iio_channel *rxb_chan_real, *rxb_chan_imag;

static void ad9361_sleep_ms(void)
{
#ifdef _WIN32
    Sleep(1); /* milliseconds */
#else
    struct timespec time;

    time.tv_sec = 0;
    time.tv_nsec = 1000 * 1000;
    nanosleep(&time, NULL);
#endif
}


int check_fmcomms5_connected(struct iio_context *ctx)
{
    dev_rx = iio_context_find_device(ctx, DEV_RX_NAME);
    dev_rx_slave = iio_context_find_device(ctx, DEV_RX_SLAVE_NAME);
    return (dev_rx && dev_rx_slave);
}

double calculate_phase(int16_t *a, int16_t *b, int16_t *c, int16_t *d,
                       int samples)
{
    int k = 0;
    double real = 0, imag = 0;
    for (; k < samples; k++) {
        real += ((double)a[k] * (double)c[k] + (double)b[k] * (double)d[k]);
        imag += ((double)a[k] * (double)d[k] - (double)b[k] * (double)c[k]);
    }
    return atan2(imag, real);
}

int streaming_interfaces(bool enable)
{
    if (enable) {
        rxa_chan_real = iio_device_find_channel(dev_rx, "voltage0", false);
        rxa_chan_imag = iio_device_find_channel(dev_rx, "voltage1", false);
        if (same_chip) {
            rxb_chan_real = iio_device_find_channel(dev_rx, "voltage2", false);
            rxb_chan_imag = iio_device_find_channel(dev_rx, "voltage3", false);
        } else {
            rxb_chan_real = iio_device_find_channel(dev_rx, "voltage4", false);
            rxb_chan_imag = iio_device_find_channel(dev_rx, "voltage5", false);
        }
        if (!(rxa_chan_real && rxa_chan_imag && rxb_chan_real && rxb_chan_imag))
            streaming_interfaces(false);

        iio_channel_enable(rxa_chan_real);
        iio_channel_enable(rxa_chan_imag);
        iio_channel_enable(rxb_chan_real);
        iio_channel_enable(rxb_chan_imag);
        rxbuf = iio_device_create_buffer(dev_rx, SAMPLES, false);
        if (!rxbuf)
            streaming_interfaces(false);
    } else {
        if (rxbuf) {
            iio_buffer_destroy(rxbuf);
        }
        if (rxa_chan_real) {
            iio_channel_disable(rxa_chan_real);
        }
        if (rxa_chan_imag) {
            iio_channel_disable(rxa_chan_imag);
        }
        if (rxb_chan_real) {
            iio_channel_disable(rxb_chan_real);
        }
        if (rxb_chan_imag) {
            iio_channel_disable(rxb_chan_imag);
        }
        return -1;
    }
    return 0;
}

void read_buffer_data(struct iio_channel *chn, struct iio_buffer *buf,
                      void *dst, size_t len)
{
    uintptr_t src_ptr, dst_ptr = (uintptr_t)dst, end = dst_ptr + len;
    unsigned int bytes = iio_channel_get_data_format(chn)->length / 8;
    uintptr_t buf_end = (uintptr_t)iio_buffer_end(buf);
    ptrdiff_t buf_step = iio_buffer_step(buf);

    for (src_ptr = (uintptr_t)iio_buffer_first(buf, chn);
         src_ptr < buf_end && dst_ptr + bytes <= end;
         src_ptr += buf_step, dst_ptr += bytes)
        iio_channel_convert(chn, (void *)dst_ptr, (const void *)src_ptr);
}

double estimate_phase_diff(double *estimate)
{
    ssize_t nbytes_rx = iio_buffer_refill(rxbuf);
    if (!nbytes_rx)
        return nbytes_rx;

    int16_t myData0_i[SAMPLES], myData0_q[SAMPLES];
    int16_t myData2_i[SAMPLES], myData2_q[SAMPLES];

    // Read data from all channels
    read_buffer_data(rxa_chan_real, rxbuf, myData0_i, SAMPLES * sizeof(int16_t));
    read_buffer_data(rxa_chan_imag, rxbuf, myData0_q, SAMPLES * sizeof(int16_t));
    read_buffer_data(rxb_chan_real, rxbuf, myData2_i, SAMPLES * sizeof(int16_t));
    read_buffer_data(rxb_chan_imag, rxbuf, myData2_q, SAMPLES * sizeof(int16_t));

    ad9361_sleep_ms();

    *estimate =
        calculate_phase(myData0_i, myData0_q, myData2_i, myData2_q, SAMPLES) *
        180 / M_PI;

    return 0;
}

int setup_iio_devices(struct iio_context *ctx)
{
    dev_rx = iio_context_find_device(ctx, DEV_RX_NAME);
    dev_rx_slave = iio_context_find_device(ctx, DEV_RX_SLAVE_NAME);
    dev_phy = iio_context_find_device(ctx, DEV_PHY_NAME);
    dev_phy_slave = iio_context_find_device(ctx, DEV_PHY_SLAVE_NAME);
    dev_tx = iio_context_find_device(ctx, DEV_TX_NAME);
    dev_tx_slave = iio_context_find_device(ctx, DEV_TX_SLAVE_NAME);
    return (dev_rx && dev_rx_slave && dev_phy && dev_phy_slave && dev_tx &&
            dev_tx_slave);
}

int check_phase_sma(struct iio_context *ctx, double *est)
{
    int ret, g;

    // Set up devices
    if (!setup_iio_devices(ctx))
        return -ENODEV;

    // Enable channels
    if (streaming_interfaces(true) < 0)
        return -ENODEV;

    for (g=0; g<STALE_BUFFERS; g++)
        ret = estimate_phase_diff(est);
    CHECK(ret);

    // Disable channels
    streaming_interfaces(false);

    return 0;
}

int main(void)
{
    // Set up context
    struct iio_context *ctx;
    const char* uri = getenv("URI_FMCOMMS5");
    if (uri == NULL)
        exit(0);// Cant find anything don't run tests
    ctx = iio_create_context_from_uri(uri);
    if (ctx==NULL)
        exit(0);// Cant find anything don't run tests
    if (!check_fmcomms5_connected(ctx))
        exit(0);// Cant find anything don't run tests

    // Test sync
    long long freq = 0;
    int freqGHZ, ret=0;
    for (freqGHZ = 1; freqGHZ<9; freqGHZ++) {

        printf("#### Calibrating FMComms5 at LO %f GHz ####\n",(double)freqGHZ/10);
        freq = 100000000*freqGHZ;

        ret = ad9361_fmcomms5_phase_sync(ctx, freq);
        if (ret<0) {
            printf("Error: %d\n",ret);
            break;
        }
    }

#if (ENABLE_PERFORMANCE_TESTS > 0)

    // Test sync performance
    double phase_tolerance = 3; // Degrees
    // These tests assume you have a signal splitter from TX1A_A to
    // RX1A_A, RX2A_A, RX1A_B, and RX2A_B using matched length cables
    double est = 0, est2 = 0;
    for (freqGHZ = 1; freqGHZ<9; freqGHZ++) {

        printf("#### Calibrating FMComms5 at LO %f GHz ####\n",(double)freqGHZ/10);
        freq = 100000000*freqGHZ;

        ret = ad9361_fmcomms5_phase_sync(ctx, freq);
        if (ret<0) {
            printf("FS Error: %d\n",ret);
            break;
        }
        // Check results
        same_chip = 1;
        ret = check_phase_sma(ctx, &est);
        if (ret<0) {
            printf("CP1 Error: %d\n",ret);
            break;
        }

        same_chip = 0;
        est2 = 0;
        ret = check_phase_sma(ctx, &est2);
        if (ret<0) {
            printf("CP2 Error: %d\n",ret);
            break;
        }
        printf("Same Chips Phase: %f | ",est);
        printf("Accross Chips Phase: %f\n",est2);

        if ((fabs(est)>phase_tolerance) || (fabs(est2)>phase_tolerance)) {
            printf("Phase calibration not within tolerance\n");
            ret = -2;
            break;
        }
    }

#endif

    // Cleanup
    if (ctx) {
        iio_context_destroy(ctx);
    }
    exit(ret);
}
