

#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#ifdef __APPLE__
#include <iio/iio.h>
#else
#include <iio.h>
#endif
#include "ad9361.h"
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define MHZ(x) ((long long)(x*1000000.0 + .5))
#define GHZ(x) ((long long)(x*1000000000.0 + .5))

#define DEV_RX_NAME "cf-ad9361-A"
#define DEV_RX_SLAVE_NAME "cf-ad9361-B"

static struct iio_device *dev_rx, *dev_rx_slave;

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
  char * uri = "ip:192.168.3.2";
  ctx = iio_create_context_from_uri(uri);
  if (ctx==NULL)
    exit(0);// Cant find anything don't run tests
  if (!check_fmcomms5_connected(ctx))
    exit(0);// Cant find anything don't run tests
  // Sync
  int freqMHZ, ret=0;
  for (freqMHZ = 3; freqMHZ<60; freqMHZ++) {
      ret = ad9361_fmcomms5_phase_sync(ctx, MHZ(freqMHZ), GHZ(0.90846));
      if (ret<0)
          break;
  }
  // Cleanup
  if (ctx) { iio_context_destroy(ctx); }
  exit(ret);
}
