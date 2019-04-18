/*
 * libad9361
 *
 * Copyright (C) 2014 Analog Devices, Inc.
 * Author: Travis Collins <travis.collins@analog.com>
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
 *
 * */

#include <errno.h>
#include <getopt.h>
#include <iio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ad9361.h>

#define MY_NAME "ad9361_filter"

#ifdef _WIN32
#define snprintf sprintf_s
#endif


static const struct option options[] = {
    {"help", no_argument, 0, 'h'},
    {"uri", required_argument, 0, 'u'},
    {"filename", required_argument, 0, 'f'},
    {"Fpass", required_argument, 0, 'p'},
    {"Fstop", no_argument, 0, 's'},
    {"TxRFAnalogCutoff", no_argument, 0, 'x'},
    {"RxRFAnalogCutoff", no_argument, 0, 'c'},
    {"DataRate", no_argument, 0, 'r'},
    {0, 0, 0, 0},
};

static const char *options_descriptions[] = {
    "Show this help and quit.",
    "Use the context at the provided URI.",
    "Filename of output file",
    "Frequency in Hz of end of passband",
    "Frequency in Hz of start of stopband",
    "Frequency in Hz of 3dB corner of transmitter analog filter",
    "Frequency in Hz of 3dB corner of receiver analog filter",
    "Data rate of TX and RX",
};

static void usage(void)
{
    unsigned int i;

    printf("Usage:\n\t" MY_NAME " [-r <DataRate>]\n\t"
           "\n\nOptions:\n");
    for (i = 0; options[i].name; i++)
        printf("\t-%c, --%s\n\t\t\t%s\n",
               options[i].val, options[i].name,
               options_descriptions[i]);
}

int main(int argc, char **argv)
{
    FILE *fp;
    struct iio_context *ctx;
    int c, option_index = 0;
    const char *arg_uri = NULL;
    const char *filename = NULL;
    int ret;
    double ApassTx, AstopTx, ApassRx, AstopRx, maxInputFS, maxInputDB;

    unsigned long rate = 0;
    unsigned long Fpass = -1;
    unsigned long Fstop = -1;
    unsigned long wnomTX = -1;
    unsigned long wnomRX = -1;
    char *filter_data;
    bool custom_settings = false;

    while ((c = getopt_long(argc, argv, "+hu:f:p:s:at:ar:r",
                            options, &option_index)) != -1) {
        switch (c) {
        case 'h':
            usage();
            return EXIT_SUCCESS;
        case 'u':
            arg_uri = optarg;
            break;
        case 'f':
            filename = optarg;
            break;
        case 'p':
            custom_settings = true;
            Fpass = atol(optarg);
            break;
        case 's':
            custom_settings = true;
            Fstop = atol(optarg);
            break;
        case 'x':
            custom_settings = true;
            wnomTX = atol(optarg);
            break;
        case 'c':
            custom_settings = true;
            wnomRX = atol(optarg);
            break;
        case 'r':
            rate = atol(optarg);
            break;
        case '?':
            return EXIT_FAILURE;
        }
    }

    if (optind != argc) {
        fprintf(stderr, "Incorrect number of arguments.\n\n");
        usage();
        return EXIT_FAILURE;
    }

    if (rate == 0)
    {
      fprintf(stderr, "DataRate (-r) must be set.\n\n");
      usage();
      return EXIT_FAILURE;
    }

    if (custom_settings) {
        if ((Fpass == -1) || (Fstop == -1) || (wnomTX == -1) || (wnomRX == -1)) {
            fprintf(stderr, "All custom filter settings must be set.\n\n");
            usage();
            return EXIT_FAILURE;
        }
    } else {
        Fpass = rate / 3.0;
        Fstop = Fpass * 1.25;
        wnomTX = 1.6 * Fstop;
        wnomRX = 1.4 * Fstop;
    }

    // Design filter
    ret = ad9361_set_bb_rate_custom_filter_manual_file(rate, Fpass, Fstop,
            wnomTX, wnomRX, &filter_data, &ApassTx, &AstopTx, &ApassRx, &AstopRx,
             &maxInputFS, &maxInputDB);
    if (ret<0)
    {
        fprintf(stderr, "Filter generation failed %d.\n\n", ret);
        return EXIT_FAILURE;
    }
    printf("Generated Filter Results\n");
    printf(" Passband Ripple: %f dB (Tx) | %f dB (Rx)\n", ApassTx, ApassRx);
    printf(" Stopband Attenuation: %f dB (Tx) | %f dB (Rx)\n", AstopTx, AstopRx);
    printf(" Max input: %f FS (%f dB)\n", maxInputFS, maxInputDB);

    // Export filter
    if (filename!=NULL) {
        fp = fopen(filename,"w+");
        if (fp)
            fputs(filter_data,fp);
        else {
            fprintf(stderr, "Opening %s failed.\n\n", filename);
            return EXIT_FAILURE;
        }
        fclose(fp);
    } else
        printf("%s\n",filter_data);

}
