/*
 * Copyright (C) 2017 Analog Devices, Inc.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 */

#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#ifdef __APPLE__
#include <iio/iio.h>
#else
#include <iio.h>
#endif
#include "ad9361.h"
#include <stdint.h>
#include <stdlib.h>
#include <errno.h>
#include <getopt.h>

#define DEV_RX_NAME "cf-ad9361-A"
#define DEV_RX_SLAVE_NAME "cf-ad9361-B"

static struct iio_device *dev_rx, *dev_rx_slave;

int check_fmcomms5_connected(struct iio_context *ctx) {
  dev_rx = iio_context_find_device(ctx, DEV_RX_NAME);
  dev_rx_slave = iio_context_find_device(ctx, DEV_RX_SLAVE_NAME);
  return (dev_rx && dev_rx_slave);
}

#define MY_NAME "ad9361_fm5_phase_sync"

enum backend {
  LOCAL,
  NETWORK,
  AUTO,
};

static const struct option options[] = {
    {"help", no_argument, 0, 'h'},
    {"network", required_argument, 0, 'n'},
    {"uri", required_argument, 0, 'u'},
    {"sample-rate", required_argument, 0, 's'},
    {"lo", required_argument, 0, 'l'},
    {0, 0, 0, 0},
};

static const char *options_descriptions[] = {
    "Show this help and quit.",
    "Use the network backend with the provided hostname.",
    "Use the context at the provided URI.",
    "Set desired sample rate (Hz).",
    "Set desired lo frequency (Hz).",
};

static void usage(void) {
  unsigned int i;

  printf("Usage:\n\t" MY_NAME " [-s <sample-rate>]\n\t" MY_NAME
         " [-l <lo>]\n\t" MY_NAME " [-n <hostname>]\n\t" MY_NAME
         " [-u <uri>]\n\nOptions:\n");
  for (i = 0; options[i].name; i++)
    printf("\t-%c, --%s\n\t\t\t%s\n", options[i].val, options[i].name,
           options_descriptions[i]);
}

int main(int argc, char **argv) {
  struct iio_context *ctx;
  int c, option_index = 0, arg_index = 0, ip_index = 0, uri_index = 0;
  enum backend backend = LOCAL;
  bool detect_context = false, set_lo = false, set_fs = false;
  int ret;
  long long lo;
  long long fs;

  while ((c = getopt_long(argc, argv, "+hn:u:s:l:", options, &option_index)) !=
         -1) {
    switch (c) {
    case 'h':
      usage();
      return EXIT_SUCCESS;
    case 'n':
      if (backend != LOCAL) {
        fprintf(stderr, "-n and -u are mutually exclusive\n");
        return EXIT_FAILURE;
      }
      backend = NETWORK;
      arg_index += 2;
      ip_index = arg_index;
      break;
    case 'u':
      if (backend != LOCAL) {
        fprintf(stderr, "-n and -u are mutually exclusive\n");
        return EXIT_FAILURE;
      }
      backend = AUTO;
      arg_index += 2;
      uri_index = arg_index;
      break;
    case 's':
      arg_index += 2;
      fs = atoll(argv[arg_index]);
      set_fs = true;
      break;
    case 'l':
      arg_index += 2;
      lo = atoll(argv[arg_index]);
      set_lo = true;
      break;
    case '?':
      return EXIT_FAILURE;
    }
  }
  if (!set_lo){
    fprintf(stderr, "Must set lo frequency -l.\n\n");
    return EXIT_FAILURE;
  }
  if (!set_fs){
    fprintf(stderr, "Must set sample-rate -s.\n\n");
    return EXIT_FAILURE;
  }

  if (arg_index >= argc) {
    fprintf(stderr, "Incorrect number of arguments.\n\n");
    usage();
    return EXIT_FAILURE;
  }

  if (backend == NETWORK)
    ctx = iio_create_network_context(argv[ip_index]);
  else if (backend == AUTO)
    ctx = iio_create_context_from_uri(argv[uri_index]);
  else
    ctx = iio_create_default_context(); // FIXME

  if (!ctx) {
    if (!detect_context) {
      char buf[1024];

      iio_strerror(errno, buf, sizeof(buf));
      fprintf(stderr, "Unable to create IIO context: %s\n", buf);
    }

    return EXIT_FAILURE;
  }

  printf("IIO context created with %s backend.\n", iio_context_get_name(ctx));

  if (!check_fmcomms5_connected(ctx)) {
    fprintf(stderr, "Not connected to FMComms5 board.\n");
    return EXIT_FAILURE;
  }

  ret = ad9361_fmcomms5_phase_sync(ctx, fs, lo);
  if (ret < 0){
    fprintf(stderr, "Failed to sync FMComms5 board.\n");
  }
  else
    printf("FMComms5 RX/TX channels have been synchronized.\n");

  if (ctx) {
    iio_context_destroy(ctx);
  }

  exit(0);
}
