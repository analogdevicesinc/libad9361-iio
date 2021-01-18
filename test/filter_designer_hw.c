

#include "ad9361.h"
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __APPLE__
#include <iio/iio.h>
#else
#include <iio.h>
#endif

int main(void)
{
    int ret, k;
    unsigned long Fpass, Fstop, wnomTX, wnomRX;
    struct iio_context *ctx;
    struct iio_device *dev;

    unsigned long rates[] = {1000000, 10000000, 20000000, 60000000};

    const char* uri = getenv("URI_AD9361");
    if (uri == NULL)
        exit(0);// Cant find anything don't run tests
    ctx = iio_create_context_from_uri(uri);
    if (ctx == NULL)
        exit(0);// Cant find anything don't run tests
    dev = iio_context_find_device(ctx, "ad9361-phy");

    for (k = 0; k < 4; k++) {

        printf("Testing rate: %lu\n",rates[k]);

        ret = ad9361_set_bb_rate_custom_filter_auto(dev, rates[k]);
        if (ret<0)
            return ret;

        Fpass = rates[k] / 3.0;
        Fstop = Fpass * 1.25;
        wnomTX = 1.6 * Fstop;
        wnomRX = 1.4 * Fstop;

        ret = ad9361_set_bb_rate_custom_filter_manual(dev, rates[k], Fpass, Fstop,
                wnomTX, wnomRX);
        if (ret<0)
            return ret;

    }
    return 0;
}
