

#include "ad9361.h"
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#ifdef __APPLE__
#include <iio/iio.h>
#else
#include <iio.h>
#endif

int main(void)
{
    int ret, k, g;
    struct iio_context *ctx;
    struct iio_device *dev;
    struct iio_channel *chan;
    long long current_rate;

    unsigned long rates[] = {520888, 600000, 1000000, 10000000, 20000000, 46000000};

    char * uri = "ip:192.168.2.1";
    ctx = iio_create_context_from_uri(uri);
    if (ctx == NULL)
        exit(0);// Cant find anything don't run tests
    dev = iio_context_find_device(ctx, "ad9361-phy");

    for (g = 0; g < 2; g++) {
        printf("########################\n");
        printf("Loop %d\n",g);
        printf("########################\n");

        for (k = 5; k >= 0; k--) {
            printf("Testing rate: %lu\n",rates[k]);

            ret = ad9361_set_bb_rate(dev, rates[k]);
            if (ret<0) {
                printf("ad9361_set_bb_rate Failed: %d\n",ret);
                return ret;
            }
            // Checks
            chan = iio_device_find_channel(dev, "voltage0", true);
            if (chan == NULL)
                return -ENODEV;
            ret = iio_channel_attr_read_longlong(chan, "sampling_frequency", &current_rate);
            if (ret < 0)
                return ret;
            if (abs(current_rate != (long long) rates[k])> 2)
                return -1;

            printf("FIR rate check passed\n");

        }
    }

    ////////////////////////////////////////////////////////////////////////////
    printf("########################\n");
    printf("Running tests backwards\n");
    printf("########################\n");

    for (g = 0; g < 2; g++) {
        printf("########################\n");
        printf("Loop %d\n",g);
        printf("########################\n");

        for (k = 0; k < 6; k++) {
            printf("Testing rate: %lu\n",rates[k]);

            ret = ad9361_set_bb_rate(dev, rates[k]);
            if (ret<0) {
                printf("ad9361_set_bb_rate Failed: %d\n",ret);
                return ret;
            }
            // Checks
            chan = iio_device_find_channel(dev, "voltage0", true);
            if (chan == NULL)
                return -ENODEV;
            ret = iio_channel_attr_read_longlong(chan, "sampling_frequency", &current_rate);
            if (ret < 0)
                return ret;
            if (abs(current_rate != (long long) rates[k])> 2)
                return -1;

            printf("FIR rate check passed\n");

        }
    }
    printf("Overall Passed\n");


    return 0;
}
