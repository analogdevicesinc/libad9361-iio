
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

#define ENABLE_PERFORMANCE_TESTS 1

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



int main(void)
{
    // Set up context
    struct iio_context *ctx;
    struct iio_device *dev;
    char * uri = "ip:analog";
    ctx = iio_create_context_from_uri(uri);
    if (ctx==NULL)
        exit(0);// Cant find anything don't run tests
    if (!check_fmcomms5_connected(ctx))
        exit(0);// Cant find anything don't run tests

    // Set up digital loopback
    dev = iio_context_find_device(ctx, "ad9361-phy");
    iio_device_debug_attr_write_bool(dev, "loopback", 1);
    dev = iio_context_find_device(ctx, "ad9361-phy-B");
    iio_device_debug_attr_write_bool(dev, "loopback", 1);

    // Test sync
    long long freq = 0;
    int freqGHZ, ret=0;
    for (freqGHZ = 1; freqGHZ<30; freqGHZ++) {

        printf("#### Calibrating FMComms5 at LO %f GHz ####\n",(double)1/10);
        // freq = 100000000*freqGHZ;
        freq = 100000000;

        printf("######\nBuffer 1\n");
        ret = ad9361_fmcomms5_phase_sync(ctx, freq);
        if (ret<0) {
            printf("Error: %d\n",ret);
            break;
        }

        printf("######\nBuffer 2\n");
        ret = ad9361_fmcomms5_phase_sync(ctx, freq);
        if (ret<0) {
            printf("Error: %d\n",ret);
            break;
        }
    }

    // Cleanup
    if (ctx) {
        iio_context_destroy(ctx);
    }
    exit(ret);
}
