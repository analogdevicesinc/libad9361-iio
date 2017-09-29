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

/* Include Files */
#include "rt_nonfinite.h"
#include "internal_design_filter_cg.h"

/* Type Definitions */
#include <stdio.h>
#ifndef struct_emxArray__common
#define struct_emxArray__common

struct emxArray__common
{
  void *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray__common*/

#ifndef typedef_emxArray__common
#define typedef_emxArray__common

typedef struct emxArray__common emxArray__common;

#endif                                 /*typedef_emxArray__common*/

#ifndef struct_emxArray_boolean_T
#define struct_emxArray_boolean_T

struct emxArray_boolean_T
{
  boolean_T *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_boolean_T*/

#ifndef typedef_emxArray_boolean_T
#define typedef_emxArray_boolean_T

typedef struct emxArray_boolean_T emxArray_boolean_T;

#endif                                 /*typedef_emxArray_boolean_T*/

#ifndef struct_emxArray_cint8_T
#define struct_emxArray_cint8_T

struct emxArray_cint8_T
{
  cint8_T *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_cint8_T*/

#ifndef typedef_emxArray_cint8_T
#define typedef_emxArray_cint8_T

typedef struct emxArray_cint8_T emxArray_cint8_T;

#endif                                 /*typedef_emxArray_cint8_T*/

#ifndef struct_emxArray_creal_T
#define struct_emxArray_creal_T

struct emxArray_creal_T
{
  creal_T *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_creal_T*/

#ifndef typedef_emxArray_creal_T
#define typedef_emxArray_creal_T

typedef struct emxArray_creal_T emxArray_creal_T;

#endif                                 /*typedef_emxArray_creal_T*/

#ifndef struct_emxArray_int32_T
#define struct_emxArray_int32_T

struct emxArray_int32_T
{
  int *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_int32_T*/

#ifndef typedef_emxArray_int32_T
#define typedef_emxArray_int32_T

typedef struct emxArray_int32_T emxArray_int32_T;

#endif                                 /*typedef_emxArray_int32_T*/

#ifndef struct_emxArray_int8_T
#define struct_emxArray_int8_T

struct emxArray_int8_T
{
  signed char *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_int8_T*/

#ifndef typedef_emxArray_int8_T
#define typedef_emxArray_int8_T

typedef struct emxArray_int8_T emxArray_int8_T;

#endif                                 /*typedef_emxArray_int8_T*/

#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
  double *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_real_T*/

#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T

typedef struct emxArray_real_T emxArray_real_T;

#endif                                 /*typedef_emxArray_real_T*/

#ifndef struct_emxArray_uint32_T
#define struct_emxArray_uint32_T

struct emxArray_uint32_T
{
  unsigned int *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_uint32_T*/

#ifndef typedef_emxArray_uint32_T
#define typedef_emxArray_uint32_T

typedef struct emxArray_uint32_T emxArray_uint32_T;

#endif                                 /*typedef_emxArray_uint32_T*/

#ifndef typedef_struct_T
#define typedef_struct_T

typedef struct {
  double nfft[2048];
  double Fs;
  double w[2048];
  double fvflag;
  char range[8];
  double centerdc;
  char configlevel[7];
} struct_T;

#endif                                 /*typedef_struct_T*/

/* Function Declarations */
static void analogresp(const char type[2], const double f[2048], double
  Fconverter, const double b1_data[], const int b1_size[2], const
  emxArray_creal_T *a1, const double b2_data[], const int b2_size[2], const
  emxArray_creal_T *a2, creal_T abc[2048]);
static void b_abs(const emxArray_creal_T *x, emxArray_real_T *y);
static void b_analogresp(const char type[2], const emxArray_real_T *f, double
  Fconverter, const double b1_data[], const int b1_size[2], const
  emxArray_creal_T *a1, const double b2_data[], const int b2_size[2], const
  emxArray_creal_T *a2, emxArray_creal_T *abc);
static void b_butter_cg(double Wn, double num[4], emxArray_creal_T *den);
static void b_cos(emxArray_real_T *x);
static void b_exp(creal_T x[2048]);
static void b_firfreqz(const double b[15], const struct_T *options, creal_T h
  [2048], double w[2048]);
static void b_firpm_cg(double order, const double ff[4], const emxArray_real_T
  *amplitudes, const emxArray_real_T *frequencies, const emxArray_real_T
  *weights, emxArray_real_T *h, boolean_T *valid, double *err);
static void b_fix(double *x);
static void b_freqs_cg(const double b_data[], const int b_size[2], const
  emxArray_creal_T *a, const emxArray_real_T *w, emxArray_creal_T *h);
static void b_freqz_cg(const double b[15], const double w[2048], double Fs,
  creal_T hh[2048]);
static void b_generateCascadedResponseRx(const char enables[4], const
  emxArray_real_T *w, double Fs, const double hb1_coeff[15], const double
  hb2_coeff[7], const double hb3_coeff_data[], const int hb3_coeff_size[2],
  const double dec_int3_coeff_data[], const int dec_int3_coeff_size[2],
  emxArray_creal_T *combinedResponse);
static double b_log2(double x);
static void b_mldivide(const double A[4], const double B[2], double Y[2]);
static void b_polyval(const double p[7], const creal_T x[2048], creal_T y[2048]);
static void b_power(const emxArray_real_T *a, emxArray_real_T *y);
static void b_rdivide(const emxArray_creal_T *x, const emxArray_creal_T *y,
                      emxArray_creal_T *z);
static void b_sin(double x[2048]);
static void b_sinc(emxArray_real_T *x, emxArray_real_T *y);
static void b_sqrt(creal_T *x);
static boolean_T b_strcmp(const char a[2]);
static double b_sum(const emxArray_real_T *x);
static void b_us(const double o[7], double u[7]);
static void b_xscal(int n, const creal_T a, emxArray_creal_T *x, int ix0, int
                    incx);
static void b_xzlartg(const creal_T f, const creal_T g, double *cs, creal_T *sn);
static void butter_cg(double Wn, double num[2], creal_T den_data[], int
                      den_size[2]);
static void c_abs(const emxArray_real_T *x, emxArray_real_T *y);
static void c_analogresp(const emxArray_real_T *f, double Fconverter, const
  double b1_data[], const int b1_size[2], const emxArray_creal_T *a1, const
  double b2_data[], const int b2_size[2], const emxArray_creal_T *a2,
  emxArray_creal_T *abc);
static void c_exp(emxArray_creal_T *x);
static void c_firfreqz(const double b[7], const struct_T *options, creal_T h
  [2048], double w[2048]);
static void c_fix(emxArray_real_T *x);
static void c_freqz_cg(const double b[7], const double w[2048], double Fs,
  creal_T hh[2048]);
static void c_generateCascadedResponseRx(const char enables[4], const
  emxArray_real_T *w, double Fs, const double hb1_coeff[15], const double
  hb2_coeff[7], const double hb3_coeff_data[], const int hb3_coeff_size[2],
  const double dec_int3_coeff_data[], const int dec_int3_coeff_size[2], const
  double extraTaps_data[], const int extraTaps_size[2], emxArray_creal_T
  *combinedResponse);
static void c_polyval(const double p_data[], const int p_size[2], const creal_T
                      x[2048], creal_T y[2048]);
static void c_rdivide(const emxArray_real_T *x, const emxArray_real_T *y,
                      emxArray_real_T *z);
static boolean_T c_strcmp(const char a[2]);
static void c_us(const double o_data[], const int o_size[2], double u_data[],
                 int u_size[2]);
static int cfprintf(const char * varargin_1);
static double dBinv(double dBinput);
static void d_analogresp(const emxArray_real_T *f, double Fconverter, const
  double b1_data[], const int b1_size[2], const emxArray_creal_T *a1, const
  double b2_data[], const int b2_size[2], const emxArray_creal_T *a2,
  emxArray_creal_T *abc);
static void d_firfreqz(const double b_data[], const int b_size[2], const
  struct_T *options, creal_T h[2048], double w[2048]);
static void d_freqz_cg(const double b_data[], const int b_size[2], const double
  w[2048], double Fs, creal_T hh[2048]);
static void d_polyval(const double p[30], const creal_T x[2048], creal_T y[2048]);
static void d_us(const double o[15], double u[30]);
static int div_s32_floor(int numerator, int denominator);
static void e_firfreqz(const double b[30], const struct_T *options, creal_T h
  [2048], double w[2048]);
static void e_freqz_cg(const double b[30], const double w[2048], double Fs,
  creal_T hh[2048]);
static void e_polyval(const double p[60], const creal_T x[2048], creal_T y[2048]);
static void e_us(const double o[15], double u[60]);
static int eml_zlahqr(emxArray_creal_T *h);
static void emxEnsureCapacity(emxArray__common *emxArray, int oldNumel, int
  elementSize);
static void emxFree_boolean_T(emxArray_boolean_T **pEmxArray);
static void emxFree_cint8_T(emxArray_cint8_T **pEmxArray);
static void emxFree_creal_T(emxArray_creal_T **pEmxArray);
static void emxFree_int32_T(emxArray_int32_T **pEmxArray);
static void emxFree_int8_T(emxArray_int8_T **pEmxArray);
static void emxFree_real_T(emxArray_real_T **pEmxArray);
static void emxFree_uint32_T(emxArray_uint32_T **pEmxArray);
static void emxInit_boolean_T(emxArray_boolean_T **pEmxArray, int numDimensions);
static void emxInit_cint8_T(emxArray_cint8_T **pEmxArray, int numDimensions);
static void emxInit_creal_T(emxArray_creal_T **pEmxArray, int numDimensions);
static void emxInit_creal_T1(emxArray_creal_T **pEmxArray, int numDimensions);
static void emxInit_int32_T(emxArray_int32_T **pEmxArray, int numDimensions);
static void emxInit_int32_T1(emxArray_int32_T **pEmxArray, int numDimensions);
static void emxInit_int8_T(emxArray_int8_T **pEmxArray, int numDimensions);
static void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions);
static void emxInit_real_T1(emxArray_real_T **pEmxArray, int numDimensions);
static void emxInit_uint32_T(emxArray_uint32_T **pEmxArray, int numDimensions);
static void f_firfreqz(const double b[60], const struct_T *options, creal_T h
  [2048], double w[2048]);
static void f_freqz_cg(const double b[60], const double w[2048], double Fs,
  creal_T hh[2048]);
static void f_polyval(const double p[14], const creal_T x[2048], creal_T y[2048]);
static void f_us(const double o[7], double u[14]);
static void fileManager(FILE * *f, boolean_T *a);
static void firfreqz(const struct_T *options, creal_T h[2048], double w[2048]);
static void firpm_cg(double order, const double ff[4], const emxArray_real_T
                     *amplitudes, const emxArray_real_T *frequencies, const
                     emxArray_real_T *weights, emxArray_real_T *h);
static void firpmgrid_cg(double nfilt, const double ff[4], emxArray_real_T
  *gridactual);
static void freqs_cg(const double b_data[], const int b_size[2], const
                     emxArray_creal_T *a, const double w[2048], creal_T h[2048]);
static void freqz_cg(const double w[2048], double Fs, creal_T hh[2048]);
static void g_firfreqz(const double b[14], const struct_T *options, creal_T h
  [2048], double w[2048]);
static void g_freqz_cg(const double b[14], const double w[2048], double Fs,
  creal_T hh[2048]);
static void g_polyval(const double p_data[], const int p_size[2], const
                      emxArray_creal_T *x, emxArray_creal_T *y);
static void genWeights(double filterOrder, const double bands[4], const
  emxArray_real_T *frequencies, const emxArray_real_T *weights, const
  emxArray_real_T *amplitudes, emxArray_real_T *grid, emxArray_real_T *des,
  emxArray_real_T *wt);
static void generateCascadedResponseRx(const char enables[4], const double w
  [2048], double Fs, const double hb1_coeff[15], const double hb2_coeff[7],
  const double hb3_coeff_data[], const int hb3_coeff_size[2], const double
  dec_int3_coeff_data[], const int dec_int3_coeff_size[2], creal_T
  combinedResponse[2048]);
static void h_freqz_cg(const emxArray_real_T *w, double Fs, emxArray_creal_T *hh);
static void i_freqz_cg(const double b[15], const emxArray_real_T *w, double Fs,
  emxArray_creal_T *hh);
static void j_freqz_cg(const double b[7], const emxArray_real_T *w, double Fs,
  emxArray_creal_T *hh);
static void k_freqz_cg(const double b_data[], const int b_size[2], const
  emxArray_real_T *w, double Fs, emxArray_creal_T *hh);
static void l_freqz_cg(const double b[30], const emxArray_real_T *w, double Fs,
  emxArray_creal_T *hh);
static void lp2lp_cg(const emxArray_creal_T *a, const emxArray_real_T *b, double
                     d, double wo, emxArray_creal_T *at, emxArray_real_T *bt,
                     double *dt);
static void m_freqz_cg(const double b[60], const emxArray_real_T *w, double Fs,
  emxArray_creal_T *hh);
static double mag2db(double y);
static void mldivide(const double A[4], const double B[4], double Y[4]);
static void n_freqz_cg(const double b[14], const emxArray_real_T *w, double Fs,
  emxArray_creal_T *hh);
static void poly(const emxArray_creal_T *x, emxArray_creal_T *c);
static void polyval(const double p[15], const creal_T x[2048], creal_T y[2048]);
static void power(const double a[2048], double y[2048]);
static void rdivide(const emxArray_real_T *x, double y, emxArray_real_T *z);
static creal_T recip(const creal_T y);
static double remezdd(double k, double n, double m, const emxArray_real_T *x);
static void remezm(double nfilt, const double edge[4], const emxArray_real_T
                   *grid, emxArray_real_T *des, emxArray_real_T *wt,
                   emxArray_real_T *h, double *dev, boolean_T *valid);
static void removeTrailingZero(const double b_data[], const int b_size[2], const
  emxArray_creal_T *a, double bR_data[], int bR_size[2], emxArray_creal_T *aR);
static double rt_atan2d_snf(double u0, double u1);
static double rt_hypotd_snf(double u0, double u1);
static double rt_powd_snf(double u0, double u1);
static double rt_remd_snf(double u0, double u1);
static double rt_roundd_snf(double u);
static void sinc(double x[2048], double y[2048]);
static double sum(const double x[2048]);
static void us(const double o[15], double u[15]);
static void vector_poly(const emxArray_creal_T *x, emxArray_creal_T *c);
static double xdlapy3(double x1, double x2, double x3);
static void xgehrd(emxArray_creal_T *a);
static double xnrm2(int n, const emxArray_creal_T *x, int ix0);
static void xscal(int n, const creal_T a, emxArray_creal_T *x, int ix0);
static void xzgeev(const emxArray_creal_T *A, int *info, emxArray_creal_T
                   *alpha1, emxArray_creal_T *beta1);
static void xzhgeqz(const emxArray_creal_T *A, int ilo, int ihi, int *info,
                    emxArray_creal_T *alpha1, emxArray_creal_T *beta1);
static void xzlartg(const creal_T f, const creal_T g, double *cs, creal_T *sn,
                    creal_T *r);
static void zp2ss_cg(emxArray_creal_T *a, emxArray_real_T *b, emxArray_creal_T
                     *c, double *d);

/* Function Definitions */

/*
 * UNTITLED Summary of this function goes here
 *    Detailed explanation goes here
 * Arguments    : const char type[2]
 *                const double f[2048]
 *                double Fconverter
 *                const double b1_data[]
 *                const int b1_size[2]
 *                const emxArray_creal_T *a1
 *                const double b2_data[]
 *                const int b2_size[2]
 *                const emxArray_creal_T *a2
 *                creal_T abc[2048]
 * Return Type  : void
 */
static void analogresp(const char type[2], const double f[2048], double
  Fconverter, const double b1_data[], const int b1_size[2], const
  emxArray_creal_T *a1, const double b2_data[], const int b2_size[2], const
  emxArray_creal_T *a2, creal_T abc[2048])
{
  boolean_T b_bool;
  int kstr;
  int exitg2;
  static const char cv27[2] = { 'T', 'x' };

  int exitg1;
  double b_f[2048];
  static const char cv28[2] = { 'R', 'x' };

  double a[2048];
  double dv15[2048];
  static creal_T dcv1[2048];
  double abc_im;
  double a_im;
  b_bool = false;
  kstr = 0;
  do {
    exitg2 = 0;
    if (kstr + 1 < 3) {
      if (type[kstr] != cv27[kstr]) {
        exitg2 = 1;
      } else {
        kstr++;
      }
    } else {
      b_bool = true;
      exitg2 = 1;
    }
  } while (exitg2 == 0);

  if (b_bool) {
    kstr = 0;
  } else {
    b_bool = false;
    kstr = 0;
    do {
      exitg1 = 0;
      if (kstr + 1 < 3) {
        if (type[kstr] != cv28[kstr]) {
          exitg1 = 1;
        } else {
          kstr++;
        }
      } else {
        b_bool = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);

    if (b_bool) {
      kstr = 1;
    } else {
      kstr = -1;
    }
  }

  switch (kstr) {
   case 0:
    for (kstr = 0; kstr < 2048; kstr++) {
      b_f[kstr] = f[kstr] / Fconverter;
    }

    sinc(b_f, a);
    for (kstr = 0; kstr < 2048; kstr++) {
      dv15[kstr] = 6.2831853071795862 * f[kstr];
    }

    freqs_cg(b1_data, b1_size, a1, dv15, abc);
    for (kstr = 0; kstr < 2048; kstr++) {
      dv15[kstr] = 6.2831853071795862 * f[kstr];
    }

    freqs_cg(b2_data, b2_size, a2, dv15, dcv1);
    for (kstr = 0; kstr < 2048; kstr++) {
      abc_im = a[kstr] * abc[kstr].re;
      a_im = a[kstr] * abc[kstr].im;
      abc[kstr].re = abc_im * dcv1[kstr].re - a_im * dcv1[kstr].im;
      abc[kstr].im = abc_im * dcv1[kstr].im + a_im * dcv1[kstr].re;
    }
    break;

   case 1:
    for (kstr = 0; kstr < 2048; kstr++) {
      b_f[kstr] = f[kstr] / Fconverter;
    }

    sinc(b_f, a);
    for (kstr = 0; kstr < 2048; kstr++) {
      b_f[kstr] = rt_powd_snf(a[kstr], 3.0);
      dv15[kstr] = 6.2831853071795862 * f[kstr];
    }

    freqs_cg(b1_data, b1_size, a1, dv15, abc);
    for (kstr = 0; kstr < 2048; kstr++) {
      dv15[kstr] = 6.2831853071795862 * f[kstr];
    }

    freqs_cg(b2_data, b2_size, a2, dv15, dcv1);
    for (kstr = 0; kstr < 2048; kstr++) {
      abc_im = abc[kstr].re * dcv1[kstr].im + abc[kstr].im * dcv1[kstr].re;
      abc[kstr].re = b_f[kstr] * (abc[kstr].re * dcv1[kstr].re - abc[kstr].im *
        dcv1[kstr].im);
      abc[kstr].im = b_f[kstr] * abc_im;
    }
    break;

   default:
    /*  Default to Rx */
    for (kstr = 0; kstr < 2048; kstr++) {
      b_f[kstr] = f[kstr] / Fconverter;
    }

    sinc(b_f, a);
    for (kstr = 0; kstr < 2048; kstr++) {
      b_f[kstr] = rt_powd_snf(a[kstr], 3.0);
      dv15[kstr] = 6.2831853071795862 * f[kstr];
    }

    freqs_cg(b1_data, b1_size, a1, dv15, abc);
    for (kstr = 0; kstr < 2048; kstr++) {
      dv15[kstr] = 6.2831853071795862 * f[kstr];
    }

    freqs_cg(b2_data, b2_size, a2, dv15, dcv1);
    for (kstr = 0; kstr < 2048; kstr++) {
      abc_im = abc[kstr].re * dcv1[kstr].im + abc[kstr].im * dcv1[kstr].re;
      abc[kstr].re = b_f[kstr] * (abc[kstr].re * dcv1[kstr].re - abc[kstr].im *
        dcv1[kstr].im);
      abc[kstr].im = b_f[kstr] * abc_im;
    }
    break;
  }
}

/*
 * Arguments    : const emxArray_creal_T *x
 *                emxArray_real_T *y
 * Return Type  : void
 */
static void b_abs(const emxArray_creal_T *x, emxArray_real_T *y)
{
  int k;
  k = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)y, k, (int)sizeof(double));
  for (k = 0; k + 1 <= x->size[1]; k++) {
    y->data[k] = rt_hypotd_snf(x->data[k].re, x->data[k].im);
  }
}

/*
 * UNTITLED Summary of this function goes here
 *    Detailed explanation goes here
 * Arguments    : const char type[2]
 *                const emxArray_real_T *f
 *                double Fconverter
 *                const double b1_data[]
 *                const int b1_size[2]
 *                const emxArray_creal_T *a1
 *                const double b2_data[]
 *                const int b2_size[2]
 *                const emxArray_creal_T *a2
 *                emxArray_creal_T *abc
 * Return Type  : void
 */
static void b_analogresp(const char type[2], const emxArray_real_T *f, double
  Fconverter, const double b1_data[], const int b1_size[2], const
  emxArray_creal_T *a1, const double b2_data[], const int b2_size[2], const
  emxArray_creal_T *a2, emxArray_creal_T *abc)
{
  boolean_T b_bool;
  int kstr;
  int exitg2;
  static const char cv40[2] = { 'T', 'x' };

  emxArray_creal_T *r14;
  int exitg1;
  emxArray_real_T *r15;
  emxArray_real_T *r16;
  static const char cv41[2] = { 'R', 'x' };

  emxArray_real_T *r17;
  int loop_ub;
  emxArray_real_T *r18;
  double abc_re;
  double abc_im;
  double re;
  double im;
  b_bool = false;
  kstr = 0;
  do {
    exitg2 = 0;
    if (kstr + 1 < 3) {
      if (type[kstr] != cv40[kstr]) {
        exitg2 = 1;
      } else {
        kstr++;
      }
    } else {
      b_bool = true;
      exitg2 = 1;
    }
  } while (exitg2 == 0);

  if (b_bool) {
    kstr = 0;
  } else {
    b_bool = false;
    kstr = 0;
    do {
      exitg1 = 0;
      if (kstr + 1 < 3) {
        if (type[kstr] != cv41[kstr]) {
          exitg1 = 1;
        } else {
          kstr++;
        }
      } else {
        b_bool = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);

    if (b_bool) {
      kstr = 1;
    } else {
      kstr = -1;
    }
  }

  emxInit_creal_T(&r14, 2);
  emxInit_real_T(&r15, 2);
  emxInit_real_T(&r16, 2);
  switch (kstr) {
   case 0:
    emxInit_real_T(&r17, 2);
    rdivide(f, Fconverter, r15);
    b_sinc(r15, r16);
    kstr = r17->size[0] * r17->size[1];
    r17->size[0] = 1;
    r17->size[1] = f->size[1];
    emxEnsureCapacity((emxArray__common *)r17, kstr, (int)sizeof(double));
    loop_ub = f->size[0] * f->size[1];
    for (kstr = 0; kstr < loop_ub; kstr++) {
      r17->data[kstr] = 6.2831853071795862 * f->data[kstr];
    }

    emxInit_real_T(&r18, 2);
    b_freqs_cg(b1_data, b1_size, a1, r17, abc);
    kstr = r18->size[0] * r18->size[1];
    r18->size[0] = 1;
    r18->size[1] = f->size[1];
    emxEnsureCapacity((emxArray__common *)r18, kstr, (int)sizeof(double));
    loop_ub = f->size[0] * f->size[1];
    emxFree_real_T(&r17);
    for (kstr = 0; kstr < loop_ub; kstr++) {
      r18->data[kstr] = 6.2831853071795862 * f->data[kstr];
    }

    b_freqs_cg(b2_data, b2_size, a2, r18, r14);
    kstr = abc->size[0] * abc->size[1];
    abc->size[0] = 1;
    abc->size[1] = r16->size[1];
    emxEnsureCapacity((emxArray__common *)abc, kstr, (int)sizeof(creal_T));
    loop_ub = r16->size[0] * r16->size[1];
    emxFree_real_T(&r18);
    for (kstr = 0; kstr < loop_ub; kstr++) {
      abc_re = r16->data[kstr] * abc->data[kstr].re;
      abc_im = r16->data[kstr] * abc->data[kstr].im;
      re = r14->data[kstr].re;
      im = r14->data[kstr].im;
      abc->data[kstr].re = abc_re * re - abc_im * im;
      abc->data[kstr].im = abc_re * im + abc_im * re;
    }
    break;

   case 1:
    emxInit_real_T(&r17, 2);
    kstr = r17->size[0] * r17->size[1];
    r17->size[0] = 1;
    r17->size[1] = f->size[1];
    emxEnsureCapacity((emxArray__common *)r17, kstr, (int)sizeof(double));
    loop_ub = f->size[0] * f->size[1];
    for (kstr = 0; kstr < loop_ub; kstr++) {
      r17->data[kstr] = 6.2831853071795862 * f->data[kstr];
    }

    emxInit_real_T(&r18, 2);
    b_freqs_cg(b1_data, b1_size, a1, r17, abc);
    kstr = r18->size[0] * r18->size[1];
    r18->size[0] = 1;
    r18->size[1] = f->size[1];
    emxEnsureCapacity((emxArray__common *)r18, kstr, (int)sizeof(double));
    loop_ub = f->size[0] * f->size[1];
    emxFree_real_T(&r17);
    for (kstr = 0; kstr < loop_ub; kstr++) {
      r18->data[kstr] = 6.2831853071795862 * f->data[kstr];
    }

    b_freqs_cg(b2_data, b2_size, a2, r18, r14);
    rdivide(f, Fconverter, r15);
    b_sinc(r15, r16);
    b_power(r16, r15);
    kstr = abc->size[0] * abc->size[1];
    abc->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)abc, kstr, (int)sizeof(creal_T));
    kstr = abc->size[0];
    loop_ub = abc->size[1];
    loop_ub *= kstr;
    emxFree_real_T(&r18);
    for (kstr = 0; kstr < loop_ub; kstr++) {
      abc_re = abc->data[kstr].re * r14->data[kstr].re - abc->data[kstr].im *
        r14->data[kstr].im;
      abc_im = abc->data[kstr].re * r14->data[kstr].im + abc->data[kstr].im *
        r14->data[kstr].re;
      abc->data[kstr].re = r15->data[kstr] * abc_re;
      abc->data[kstr].im = r15->data[kstr] * abc_im;
    }
    break;

   default:
    emxInit_real_T(&r17, 2);

    /*  Default to Rx */
    kstr = r17->size[0] * r17->size[1];
    r17->size[0] = 1;
    r17->size[1] = f->size[1];
    emxEnsureCapacity((emxArray__common *)r17, kstr, (int)sizeof(double));
    loop_ub = f->size[0] * f->size[1];
    for (kstr = 0; kstr < loop_ub; kstr++) {
      r17->data[kstr] = 6.2831853071795862 * f->data[kstr];
    }

    emxInit_real_T(&r18, 2);
    b_freqs_cg(b1_data, b1_size, a1, r17, abc);
    kstr = r18->size[0] * r18->size[1];
    r18->size[0] = 1;
    r18->size[1] = f->size[1];
    emxEnsureCapacity((emxArray__common *)r18, kstr, (int)sizeof(double));
    loop_ub = f->size[0] * f->size[1];
    emxFree_real_T(&r17);
    for (kstr = 0; kstr < loop_ub; kstr++) {
      r18->data[kstr] = 6.2831853071795862 * f->data[kstr];
    }

    b_freqs_cg(b2_data, b2_size, a2, r18, r14);
    rdivide(f, Fconverter, r15);
    b_sinc(r15, r16);
    b_power(r16, r15);
    kstr = abc->size[0] * abc->size[1];
    abc->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)abc, kstr, (int)sizeof(creal_T));
    kstr = abc->size[0];
    loop_ub = abc->size[1];
    loop_ub *= kstr;
    emxFree_real_T(&r18);
    for (kstr = 0; kstr < loop_ub; kstr++) {
      abc_re = abc->data[kstr].re * r14->data[kstr].re - abc->data[kstr].im *
        r14->data[kstr].im;
      abc_im = abc->data[kstr].re * r14->data[kstr].im + abc->data[kstr].im *
        r14->data[kstr].re;
      abc->data[kstr].re = r15->data[kstr] * abc_re;
      abc->data[kstr].im = r15->data[kstr] * abc_im;
    }
    break;
  }

  emxFree_real_T(&r16);
  emxFree_real_T(&r15);
  emxFree_creal_T(&r14);
}

/*
 * BUTTER_CG Butterworth digital and analog filter design.  Codegen support
 *
 *  This function is based on 'butter' by The MathWorks Inc.
 * Arguments    : double Wn
 *                double num[4]
 *                emxArray_creal_T *den
 * Return Type  : void
 */
static void b_butter_cg(double Wn, double num[4], emxArray_creal_T *den)
{
  emxArray_creal_T *a;
  emxArray_real_T *b;
  emxArray_creal_T *c;
  emxArray_creal_T *b_a;
  emxArray_real_T *b_b;
  double d;
  double b_d;
  double y_re;
  double y_im;
  int k;
  static const double c_a[4] = { 0.0, 0.0, 0.0, 0.037037037037037035 };

  emxInit_creal_T(&a, 2);
  emxInit_real_T1(&b, 1);
  emxInit_creal_T(&c, 2);
  emxInit_creal_T(&b_a, 2);
  emxInit_real_T1(&b_b, 1);

  /*  Set types we will only use */
  /*  Cast to enforce precision rules */
  /*  step 1: get analog, pre-warped frequencies */
  /*  step 2: convert to low-pass prototype estimate */
  /*  lowpass */
  /*  Cast to enforce precision rules */
  /*  step 3: Get N-th order Butterworth analog lowpass prototype */
  /*  Transform to state-space */
  zp2ss_cg(a, b, c, &d);

  /*  step 4: Transform to lowpass, bandpass, highpass, or bandstop of desired Wn */
  /*  Lowpass */
  lp2lp_cg(a, b, d, Wn, b_a, b_b, &b_d);

  /*  step 5: Use Bilinear transformation to find discrete equivalent: */
  /*  nargout <= 3 */
  /*  Transform to zero-pole-gain and polynomial forms: */
  /*  nargout <= 2 */
  poly(b_a, den);

  /*  This internal function returns more exact numerator vectors */
  /*  for the num/den case. */
  /*  Wn input is two element band edge vector */
  /* --------------------------------- */
  /*  lowpass */
  y_re = den->data[0].re;
  y_im = den->data[0].im;
  k = 0;
  emxFree_real_T(&b_b);
  emxFree_creal_T(&b_a);
  emxFree_creal_T(&c);
  emxFree_real_T(&b);
  emxFree_creal_T(&a);
  while (k <= den->size[1] - 2) {
    d = 0.0 * y_im + 0.0 * y_re;
    y_re = (0.0 * y_re - 0.0 * y_im) + den->data[k + 1].re;
    y_im = d + den->data[k + 1].im;
    k++;
  }

  for (k = 0; k < 4; k++) {
    d = c_a[k] * y_re;
    b_d = c_a[k] * y_im;
    if (b_d == 0.0) {
      d /= 0.037037037037037035;
    } else if (d == 0.0) {
      d = 0.0;
    } else {
      d /= 0.037037037037037035;
    }

    num[k] = d;
  }

  /*  num = poly(a-b*c)+(d-1)*den; */
}

/*
 * Arguments    : emxArray_real_T *x
 * Return Type  : void
 */
static void b_cos(emxArray_real_T *x)
{
  int nx;
  int k;
  nx = x->size[1];
  for (k = 0; k + 1 <= nx; k++) {
    x->data[k] = cos(x->data[k]);
  }
}

/*
 * Arguments    : creal_T x[2048]
 * Return Type  : void
 */
static void b_exp(creal_T x[2048])
{
  int k;
  double x_re;
  double r;
  for (k = 0; k < 2048; k++) {
    if (x[k].im == 0.0) {
      x_re = exp(x[k].re);
      r = 0.0;
    } else if (rtIsInf(x[k].im) && rtIsInf(x[k].re) && (x[k].re < 0.0)) {
      x_re = 0.0;
      r = 0.0;
    } else {
      r = exp(x[k].re / 2.0);
      x_re = r * (r * cos(x[k].im));
      r *= r * sin(x[k].im);
    }

    x[k].re = x_re;
    x[k].im = r;
  }
}

/*
 * Make b a row
 * Arguments    : const double b[15]
 *                const struct_T *options
 *                creal_T h[2048]
 *                double w[2048]
 * Return Type  : void
 */
static void b_firfreqz(const double b[15], const struct_T *options, creal_T h
  [2048], double w[2048])
{
  double digw[2048];
  creal_T dcv2[2048];
  int i55;
  double bim;
  double h_re;
  double brm;
  double d;

  /* -------------------------------------------------------------------------- */
  /*  Actual Frequency Response Computation */
  /* if fvflag, */
  /*    Frequency vector specified.  Use Horner's method of polynomial */
  /*    evaluation at the frequency points and divide the numerator */
  /*    by the denominator. */
  /*  */
  /*    Note: we use positive i here because of the relationship */
  /*             polyval(a,exp(1i*w)) = fft(a).*exp(1i*w*(length(a)-1)) */
  /*                ( assuming w = 2*pi*(0:length(a)-1)/length(a) ) */
  /*  */
  /*  Fs was specified, freq. vector is in Hz */
  /*  Convert from Hz to rad/sample for computational purposes */
  /*  Digital frequency must be used for this calculation */
  for (i55 = 0; i55 < 2048; i55++) {
    w[i55] = options->w[i55];
    bim = 6.2831853071795862 * options->w[i55] / options->Fs;
    dcv2[i55].re = bim * 0.0;
    dcv2[i55].im = bim;
    digw[i55] = bim;
  }

  b_exp(dcv2);
  polyval(b, dcv2, h);
  for (i55 = 0; i55 < 2048; i55++) {
    dcv2[i55].re = 14.0 * (digw[i55] * 0.0);
    dcv2[i55].im = 14.0 * digw[i55];
  }

  b_exp(dcv2);
  for (i55 = 0; i55 < 2048; i55++) {
    h_re = h[i55].re;
    if (dcv2[i55].im == 0.0) {
      if (h[i55].im == 0.0) {
        h[i55].re /= dcv2[i55].re;
        h[i55].im = 0.0;
      } else if (h[i55].re == 0.0) {
        h[i55].re = 0.0;
        h[i55].im /= dcv2[i55].re;
      } else {
        h[i55].re /= dcv2[i55].re;
        h[i55].im /= dcv2[i55].re;
      }
    } else if (dcv2[i55].re == 0.0) {
      if (h[i55].re == 0.0) {
        h[i55].re = h[i55].im / dcv2[i55].im;
        h[i55].im = 0.0;
      } else if (h[i55].im == 0.0) {
        h[i55].re = 0.0;
        h[i55].im = -(h_re / dcv2[i55].im);
      } else {
        h[i55].re = h[i55].im / dcv2[i55].im;
        h[i55].im = -(h_re / dcv2[i55].im);
      }
    } else {
      brm = fabs(dcv2[i55].re);
      bim = fabs(dcv2[i55].im);
      if (brm > bim) {
        bim = dcv2[i55].im / dcv2[i55].re;
        d = dcv2[i55].re + bim * dcv2[i55].im;
        h[i55].re = (h[i55].re + bim * h[i55].im) / d;
        h[i55].im = (h[i55].im - bim * h_re) / d;
      } else if (bim == brm) {
        if (dcv2[i55].re > 0.0) {
          bim = 0.5;
        } else {
          bim = -0.5;
        }

        if (dcv2[i55].im > 0.0) {
          d = 0.5;
        } else {
          d = -0.5;
        }

        h[i55].re = (h[i55].re * bim + h[i55].im * d) / brm;
        h[i55].im = (h[i55].im * bim - h_re * d) / brm;
      } else {
        bim = dcv2[i55].re / dcv2[i55].im;
        d = dcv2[i55].im + bim * dcv2[i55].re;
        h[i55].re = (bim * h[i55].re + h[i55].im) / d;
        h[i55].im = (bim * h[i55].im - h_re) / d;
      }
    }
  }
}

/*
 * FIRPM Parks-McClellan optimal equiripple FIR filter design.
 *
 *  This function is based on 'firpm' by The MathWorks Inc.
 * Arguments    : double order
 *                const double ff[4]
 *                const emxArray_real_T *amplitudes
 *                const emxArray_real_T *frequencies
 *                const emxArray_real_T *weights
 *                emxArray_real_T *h
 *                boolean_T *valid
 *                double *err
 * Return Type  : void
 */
static void b_firpm_cg(double order, const double ff[4], const emxArray_real_T
  *amplitudes, const emxArray_real_T *frequencies, const emxArray_real_T
  *weights, emxArray_real_T *h, boolean_T *valid, double *err)
{
  emxArray_real_T *grid;
  emxArray_real_T *des;
  emxArray_real_T *wt;
  emxArray_real_T *r26;
  double b_ff[4];
  int i46;
  emxArray_real_T *b_h;
  double x;
  int h_idx_0;
  int i47;
  emxArray_real_T *c_h;
  int i48;
  int loop_ub;
  emxArray_real_T *d_h;
  emxInit_real_T(&grid, 2);
  emxInit_real_T(&des, 2);
  emxInit_real_T(&wt, 2);
  emxInit_real_T(&r26, 2);
  genWeights(order, ff, frequencies, weights, amplitudes, grid, des, wt);

  /*  Workaround */
  /* ftype = 2; */
  /* sign_val = 1; */
  /*  Always bandpass designs */
  /*  cast to enforce precision rules */
  /*  Call actual design algorithm */
  rdivide(grid, 2.0, r26);
  emxFree_real_T(&grid);
  for (i46 = 0; i46 < 4; i46++) {
    b_ff[i46] = ff[i46] / 2.0;
  }

  emxInit_real_T(&b_h, 2);
  remezm(order + 1.0, b_ff, r26, des, wt, b_h, &x, valid);
  *err = fabs(x);
  h_idx_0 = b_h->size[0] * b_h->size[1];
  i46 = h->size[0] * h->size[1];
  h->size[0] = 1;
  h->size[1] = h_idx_0;
  emxEnsureCapacity((emxArray__common *)h, i46, (int)sizeof(double));
  emxFree_real_T(&r26);
  emxFree_real_T(&wt);
  emxFree_real_T(&des);
  for (i46 = 0; i46 < h_idx_0; i46++) {
    h->data[h->size[0] * i46] = b_h->data[i46];
  }

  emxFree_real_T(&b_h);

  /*  make it a row */
  x = (double)h->size[1] - rt_remd_snf(order + 1.0, 2.0);
  if (1.0 > x) {
    i46 = 1;
    h_idx_0 = 1;
    i47 = 0;
  } else {
    i46 = (int)x;
    h_idx_0 = -1;
    i47 = 1;
  }

  emxInit_real_T(&c_h, 2);
  i48 = c_h->size[0] * c_h->size[1];
  c_h->size[0] = 1;
  c_h->size[1] = (h->size[1] + div_s32_floor(i47 - i46, h_idx_0)) + 1;
  emxEnsureCapacity((emxArray__common *)c_h, i48, (int)sizeof(double));
  loop_ub = h->size[1];
  for (i48 = 0; i48 < loop_ub; i48++) {
    c_h->data[c_h->size[0] * i48] = h->data[h->size[0] * i48];
  }

  loop_ub = div_s32_floor(i47 - i46, h_idx_0);
  for (i47 = 0; i47 <= loop_ub; i47++) {
    c_h->data[c_h->size[0] * (i47 + h->size[1])] = h->data[(i46 + h_idx_0 * i47)
      - 1];
  }

  i46 = h->size[0] * h->size[1];
  h->size[0] = 1;
  h->size[1] = c_h->size[1];
  emxEnsureCapacity((emxArray__common *)h, i46, (int)sizeof(double));
  loop_ub = c_h->size[1];
  for (i46 = 0; i46 < loop_ub; i46++) {
    h->data[h->size[0] * i46] = c_h->data[c_h->size[0] * i46];
  }

  emxFree_real_T(&c_h);
  if (1 > h->size[1]) {
    i46 = 1;
    h_idx_0 = 1;
    i47 = 0;
  } else {
    i46 = h->size[1];
    h_idx_0 = -1;
    i47 = 1;
  }

  emxInit_real_T(&d_h, 2);
  i48 = d_h->size[0] * d_h->size[1];
  d_h->size[0] = 1;
  d_h->size[1] = div_s32_floor(i47 - i46, h_idx_0) + 1;
  emxEnsureCapacity((emxArray__common *)d_h, i48, (int)sizeof(double));
  loop_ub = div_s32_floor(i47 - i46, h_idx_0);
  for (i47 = 0; i47 <= loop_ub; i47++) {
    d_h->data[d_h->size[0] * i47] = h->data[(i46 + h_idx_0 * i47) - 1];
  }

  i46 = h->size[0] * h->size[1];
  h->size[0] = 1;
  h->size[1] = d_h->size[1];
  emxEnsureCapacity((emxArray__common *)h, i46, (int)sizeof(double));
  loop_ub = d_h->size[1];
  for (i46 = 0; i46 < loop_ub; i46++) {
    h->data[h->size[0] * i46] = d_h->data[d_h->size[0] * i46];
  }

  emxFree_real_T(&d_h);
}

/*
 * Arguments    : double *x
 * Return Type  : void
 */
static void b_fix(double *x)
{
  if (*x < 0.0) {
    *x = ceil(*x);
  } else {
    *x = floor(*x);
  }
}

/*
 * FREQS Laplace-transform (s-domain) frequency response with codegen support
 *
 *  This function is based on 'freqs' by The MathWorks Inc.
 * Arguments    : const double b_data[]
 *                const int b_size[2]
 *                const emxArray_creal_T *a
 *                const emxArray_real_T *w
 *                emxArray_creal_T *h
 * Return Type  : void
 */
static void b_freqs_cg(const double b_data[], const int b_size[2], const
  emxArray_creal_T *a, const emxArray_real_T *w, emxArray_creal_T *h)
{
  emxArray_creal_T *s;
  emxArray_creal_T *b_a;
  double b_b_data[4];
  int b_b_size[2];
  int i35;
  int loop_ub;
  emxArray_creal_T *y;
  boolean_T b8;
  int k;
  double a_re;
  double a_im;
  double s_re;
  double s_im;
  emxInit_creal_T(&s, 2);
  emxInit_creal_T(&b_a, 2);
  removeTrailingZero(b_data, b_size, a, b_b_data, b_b_size, b_a);
  i35 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = w->size[1];
  emxEnsureCapacity((emxArray__common *)s, i35, (int)sizeof(creal_T));
  loop_ub = w->size[0] * w->size[1];
  for (i35 = 0; i35 < loop_ub; i35++) {
    s->data[i35].re = w->data[i35] * 0.0;
    s->data[i35].im = w->data[i35];
  }

  emxInit_creal_T(&y, 2);
  i35 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = s->size[1];
  emxEnsureCapacity((emxArray__common *)y, i35, (int)sizeof(creal_T));
  if ((y->size[1] == 0) || (b_a->size[1] == 0)) {
    b8 = true;
  } else {
    b8 = false;
  }

  if (!b8) {
    i35 = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)y, i35, (int)sizeof(creal_T));
    loop_ub = y->size[1];
    for (i35 = 0; i35 < loop_ub; i35++) {
      y->data[y->size[0] * i35] = b_a->data[0];
    }

    for (k = 0; k <= b_a->size[1] - 2; k++) {
      i35 = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = s->size[1];
      emxEnsureCapacity((emxArray__common *)y, i35, (int)sizeof(creal_T));
      a_re = b_a->data[k + 1].re;
      a_im = b_a->data[k + 1].im;
      loop_ub = s->size[0] * s->size[1];
      for (i35 = 0; i35 < loop_ub; i35++) {
        s_re = s->data[i35].re * y->data[i35].re - s->data[i35].im * y->data[i35]
          .im;
        s_im = s->data[i35].re * y->data[i35].im + s->data[i35].im * y->data[i35]
          .re;
        y->data[i35].re = s_re + a_re;
        y->data[i35].im = s_im + a_im;
      }
    }
  }

  g_polyval(b_b_data, b_b_size, s, b_a);
  b_rdivide(b_a, y, h);
  emxFree_creal_T(&y);
  emxFree_creal_T(&b_a);
  emxFree_creal_T(&s);
}

/*
 * FREQZ_CG Frequency response of digital filter with codegen support
 *
 *  This function is based on 'freqz' by The MathWorks Inc.
 * Arguments    : const double b[15]
 *                const double w[2048]
 *                double Fs
 *                creal_T hh[2048]
 * Return Type  : void
 */
static void b_freqz_cg(const double b[15], const double w[2048], double Fs,
  creal_T hh[2048])
{
  static struct_T options;
  int i9;
  static const char cv14[8] = { 'o', 'n', 'e', 's', 'i', 'd', 'e', 'd' };

  double b_w[2048];
  static const char cv15[7] = { 'o', 'm', 'i', 't', 't', 'e', 'd' };

  /*  Cast to enforce precision rules */
  options.Fs = Fs;
  memcpy(&options.nfft[0], &w[0], sizeof(double) << 11);
  memcpy(&options.w[0], &w[0], sizeof(double) << 11);

  /*  Remaining are default or for advanced use */
  options.fvflag = 1.0;
  for (i9 = 0; i9 < 8; i9++) {
    options.range[i9] = cv14[i9];
  }

  options.centerdc = 0.0;
  for (i9 = 0; i9 < 7; i9++) {
    options.configlevel[i9] = cv15[i9];
  }

  b_firfreqz(b, &options, hh, b_w);

  /*  Generate the default structure to pass to freqzplot */
  /*  If rad/sample, Fs is empty */
}

/*
 * Arguments    : const char enables[4]
 *                const emxArray_real_T *w
 *                double Fs
 *                const double hb1_coeff[15]
 *                const double hb2_coeff[7]
 *                const double hb3_coeff_data[]
 *                const int hb3_coeff_size[2]
 *                const double dec_int3_coeff_data[]
 *                const int dec_int3_coeff_size[2]
 *                emxArray_creal_T *combinedResponse
 * Return Type  : void
 */
static void b_generateCascadedResponseRx(const char enables[4], const
  emxArray_real_T *w, double Fs, const double hb1_coeff[15], const double
  hb2_coeff[7], const double hb3_coeff_data[], const int hb3_coeff_size[2],
  const double dec_int3_coeff_data[], const int dec_int3_coeff_size[2],
  emxArray_creal_T *combinedResponse)
{
  boolean_T b_bool;
  int ix;
  int exitg12;
  emxArray_creal_T *d2;
  int exitg11;
  emxArray_creal_T *d3;
  static const char cv29[4] = { '2', '1', '1', '1' };

  double u[15];
  double b_u[7];
  double tmp_data[29];
  int tmp_size[2];
  double c_u[30];
  double d_u[60];
  double e_u[14];
  int iy;
  int k;
  int exitg10;
  static const char cv30[4] = { '1', '2', '1', '1' };

  int exitg9;
  double combinedResponse_re;
  double combinedResponse_im;
  static const char cv31[4] = { '1', '1', '2', '1' };

  double d2_re;
  double d2_im;
  int exitg8;
  static const char cv32[4] = { '2', '2', '1', '1' };

  int exitg7;
  static const char cv33[4] = { '2', '1', '2', '1' };

  int exitg6;
  static const char cv34[4] = { '1', '2', '2', '1' };

  int exitg5;
  static const char cv35[4] = { '2', '2', '2', '1' };

  int exitg4;
  static const char cv36[4] = { '1', '1', '1', '3' };

  int exitg3;
  static const char cv37[4] = { '2', '1', '1', '3' };

  int exitg2;
  static const char cv38[4] = { '1', '2', '1', '3' };

  int exitg1;
  static const char cv39[4] = { '2', '2', '1', '3' };

  /*  Cast */
  b_bool = false;
  ix = 1;
  do {
    exitg12 = 0;
    if (ix < 5) {
      if (enables[ix - 1] != '1') {
        exitg12 = 1;
      } else {
        ix++;
      }
    } else {
      b_bool = true;
      exitg12 = 1;
    }
  } while (exitg12 == 0);

  if (b_bool) {
    ix = 0;
  } else {
    b_bool = false;
    ix = 0;
    do {
      exitg11 = 0;
      if (ix + 1 < 5) {
        if (enables[ix] != cv29[ix]) {
          exitg11 = 1;
        } else {
          ix++;
        }
      } else {
        b_bool = true;
        exitg11 = 1;
      }
    } while (exitg11 == 0);

    if (b_bool) {
      ix = 1;
    } else {
      b_bool = false;
      ix = 0;
      do {
        exitg10 = 0;
        if (ix + 1 < 5) {
          if (enables[ix] != cv30[ix]) {
            exitg10 = 1;
          } else {
            ix++;
          }
        } else {
          b_bool = true;
          exitg10 = 1;
        }
      } while (exitg10 == 0);

      if (b_bool) {
        ix = 2;
      } else {
        b_bool = false;
        ix = 0;
        do {
          exitg9 = 0;
          if (ix + 1 < 5) {
            if (enables[ix] != cv31[ix]) {
              exitg9 = 1;
            } else {
              ix++;
            }
          } else {
            b_bool = true;
            exitg9 = 1;
          }
        } while (exitg9 == 0);

        if (b_bool) {
          ix = 3;
        } else {
          b_bool = false;
          ix = 0;
          do {
            exitg8 = 0;
            if (ix + 1 < 5) {
              if (enables[ix] != cv32[ix]) {
                exitg8 = 1;
              } else {
                ix++;
              }
            } else {
              b_bool = true;
              exitg8 = 1;
            }
          } while (exitg8 == 0);

          if (b_bool) {
            ix = 4;
          } else {
            b_bool = false;
            ix = 0;
            do {
              exitg7 = 0;
              if (ix + 1 < 5) {
                if (enables[ix] != cv33[ix]) {
                  exitg7 = 1;
                } else {
                  ix++;
                }
              } else {
                b_bool = true;
                exitg7 = 1;
              }
            } while (exitg7 == 0);

            if (b_bool) {
              ix = 5;
            } else {
              b_bool = false;
              ix = 0;
              do {
                exitg6 = 0;
                if (ix + 1 < 5) {
                  if (enables[ix] != cv34[ix]) {
                    exitg6 = 1;
                  } else {
                    ix++;
                  }
                } else {
                  b_bool = true;
                  exitg6 = 1;
                }
              } while (exitg6 == 0);

              if (b_bool) {
                ix = 6;
              } else {
                b_bool = false;
                ix = 0;
                do {
                  exitg5 = 0;
                  if (ix + 1 < 5) {
                    if (enables[ix] != cv35[ix]) {
                      exitg5 = 1;
                    } else {
                      ix++;
                    }
                  } else {
                    b_bool = true;
                    exitg5 = 1;
                  }
                } while (exitg5 == 0);

                if (b_bool) {
                  ix = 7;
                } else {
                  b_bool = false;
                  ix = 0;
                  do {
                    exitg4 = 0;
                    if (ix + 1 < 5) {
                      if (enables[ix] != cv36[ix]) {
                        exitg4 = 1;
                      } else {
                        ix++;
                      }
                    } else {
                      b_bool = true;
                      exitg4 = 1;
                    }
                  } while (exitg4 == 0);

                  if (b_bool) {
                    ix = 8;
                  } else {
                    b_bool = false;
                    ix = 0;
                    do {
                      exitg3 = 0;
                      if (ix + 1 < 5) {
                        if (enables[ix] != cv37[ix]) {
                          exitg3 = 1;
                        } else {
                          ix++;
                        }
                      } else {
                        b_bool = true;
                        exitg3 = 1;
                      }
                    } while (exitg3 == 0);

                    if (b_bool) {
                      ix = 9;
                    } else {
                      b_bool = false;
                      ix = 0;
                      do {
                        exitg2 = 0;
                        if (ix + 1 < 5) {
                          if (enables[ix] != cv38[ix]) {
                            exitg2 = 1;
                          } else {
                            ix++;
                          }
                        } else {
                          b_bool = true;
                          exitg2 = 1;
                        }
                      } while (exitg2 == 0);

                      if (b_bool) {
                        ix = 10;
                      } else {
                        b_bool = false;
                        ix = 0;
                        do {
                          exitg1 = 0;
                          if (ix + 1 < 5) {
                            if (enables[ix] != cv39[ix]) {
                              exitg1 = 1;
                            } else {
                              ix++;
                            }
                          } else {
                            b_bool = true;
                            exitg1 = 1;
                          }
                        } while (exitg1 == 0);

                        if (b_bool) {
                          ix = 11;
                        } else {
                          ix = -1;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  emxInit_creal_T(&d2, 2);
  emxInit_creal_T(&d3, 2);
  switch (ix) {
   case 0:
    /*  only FIR */
    h_freqz_cg(w, Fs, combinedResponse);
    break;

   case 1:
    /*  Hb1 */
    memset(&u[0], 0, 15U * sizeof(double));
    ix = 0;
    iy = 0;
    for (k = 0; k < 15; k++) {
      u[iy] = hb1_coeff[ix];
      ix++;
      iy++;
    }

    i_freqz_cg(u, w, Fs, combinedResponse);
    break;

   case 2:
    /*  Hb2 */
    for (ix = 0; ix < 7; ix++) {
      b_u[ix] = 0.0;
    }

    ix = 0;
    iy = 0;
    for (k = 0; k < 7; k++) {
      b_u[iy] = hb2_coeff[ix];
      ix++;
      iy++;
    }

    j_freqz_cg(b_u, w, Fs, combinedResponse);
    break;

   case 3:
    /*  Hb3 */
    c_us(hb3_coeff_data, hb3_coeff_size, tmp_data, tmp_size);
    k_freqz_cg(tmp_data, tmp_size, w, Fs, combinedResponse);
    break;

   case 4:
    /*  Hb2,Hb1 */
    memset(&c_u[0], 0, 30U * sizeof(double));
    ix = 0;
    iy = 0;
    for (k = 0; k < 15; k++) {
      c_u[iy] = hb1_coeff[ix];
      ix++;
      iy += 2;
    }

    l_freqz_cg(c_u, w, Fs, combinedResponse);
    for (ix = 0; ix < 7; ix++) {
      b_u[ix] = 0.0;
    }

    ix = 0;
    iy = 0;
    for (k = 0; k < 7; k++) {
      b_u[iy] = hb2_coeff[ix];
      ix++;
      iy++;
    }

    j_freqz_cg(b_u, w, Fs, d2);
    ix = combinedResponse->size[0] * combinedResponse->size[1];
    combinedResponse->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)combinedResponse, ix, (int)sizeof
                      (creal_T));
    ix = combinedResponse->size[0];
    iy = combinedResponse->size[1];
    iy *= ix;
    for (ix = 0; ix < iy; ix++) {
      combinedResponse_re = combinedResponse->data[ix].re;
      combinedResponse_im = combinedResponse->data[ix].im;
      d2_re = d2->data[ix].re;
      d2_im = d2->data[ix].im;
      combinedResponse->data[ix].re = combinedResponse_re * d2_re -
        combinedResponse_im * d2_im;
      combinedResponse->data[ix].im = combinedResponse_re * d2_im +
        combinedResponse_im * d2_re;
    }
    break;

   case 5:
    /*  Hb3,Hb1 */
    memset(&c_u[0], 0, 30U * sizeof(double));
    ix = 0;
    iy = 0;
    for (k = 0; k < 15; k++) {
      c_u[iy] = hb1_coeff[ix];
      ix++;
      iy += 2;
    }

    l_freqz_cg(c_u, w, Fs, combinedResponse);
    for (ix = 0; ix < 7; ix++) {
      b_u[ix] = 0.0;
    }

    ix = 0;
    iy = 0;
    for (k = 0; k < 7; k++) {
      b_u[iy] = hb2_coeff[ix];
      ix++;
      iy++;
    }

    j_freqz_cg(b_u, w, Fs, d2);
    ix = combinedResponse->size[0] * combinedResponse->size[1];
    combinedResponse->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)combinedResponse, ix, (int)sizeof
                      (creal_T));
    ix = combinedResponse->size[0];
    iy = combinedResponse->size[1];
    iy *= ix;
    for (ix = 0; ix < iy; ix++) {
      combinedResponse_re = combinedResponse->data[ix].re;
      combinedResponse_im = combinedResponse->data[ix].im;
      d2_re = d2->data[ix].re;
      d2_im = d2->data[ix].im;
      combinedResponse->data[ix].re = combinedResponse_re * d2_re -
        combinedResponse_im * d2_im;
      combinedResponse->data[ix].im = combinedResponse_re * d2_im +
        combinedResponse_im * d2_re;
    }
    break;

   case 6:
    /*  Hb3,Hb2 */
    memset(&c_u[0], 0, 30U * sizeof(double));
    ix = 0;
    iy = 0;
    for (k = 0; k < 15; k++) {
      c_u[iy] = hb1_coeff[ix];
      ix++;
      iy += 2;
    }

    l_freqz_cg(c_u, w, Fs, combinedResponse);
    c_us(hb3_coeff_data, hb3_coeff_size, tmp_data, tmp_size);
    k_freqz_cg(tmp_data, tmp_size, w, Fs, d2);
    ix = combinedResponse->size[0] * combinedResponse->size[1];
    combinedResponse->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)combinedResponse, ix, (int)sizeof
                      (creal_T));
    ix = combinedResponse->size[0];
    iy = combinedResponse->size[1];
    iy *= ix;
    for (ix = 0; ix < iy; ix++) {
      combinedResponse_re = combinedResponse->data[ix].re;
      combinedResponse_im = combinedResponse->data[ix].im;
      d2_re = d2->data[ix].re;
      d2_im = d2->data[ix].im;
      combinedResponse->data[ix].re = combinedResponse_re * d2_re -
        combinedResponse_im * d2_im;
      combinedResponse->data[ix].im = combinedResponse_re * d2_im +
        combinedResponse_im * d2_re;
    }
    break;

   case 7:
    /*  Hb3,Hb2,Hb1 */
    memset(&d_u[0], 0, 60U * sizeof(double));
    ix = 0;
    iy = 0;
    for (k = 0; k < 15; k++) {
      d_u[iy] = hb1_coeff[ix];
      ix++;
      iy += 4;
    }

    m_freqz_cg(d_u, w, Fs, combinedResponse);
    memset(&e_u[0], 0, 14U * sizeof(double));
    ix = 0;
    iy = 0;
    for (k = 0; k < 7; k++) {
      e_u[iy] = hb2_coeff[ix];
      ix++;
      iy += 2;
    }

    n_freqz_cg(e_u, w, Fs, d2);
    c_us(hb3_coeff_data, hb3_coeff_size, tmp_data, tmp_size);
    k_freqz_cg(tmp_data, tmp_size, w, Fs, d3);
    ix = combinedResponse->size[0] * combinedResponse->size[1];
    combinedResponse->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)combinedResponse, ix, (int)sizeof
                      (creal_T));
    ix = combinedResponse->size[0];
    iy = combinedResponse->size[1];
    iy *= ix;
    for (ix = 0; ix < iy; ix++) {
      combinedResponse_re = combinedResponse->data[ix].re * d2->data[ix].re -
        combinedResponse->data[ix].im * d2->data[ix].im;
      combinedResponse_im = combinedResponse->data[ix].re * d2->data[ix].im +
        combinedResponse->data[ix].im * d2->data[ix].re;
      d2_re = d3->data[ix].re;
      d2_im = d3->data[ix].im;
      combinedResponse->data[ix].re = combinedResponse_re * d2_re -
        combinedResponse_im * d2_im;
      combinedResponse->data[ix].im = combinedResponse_re * d2_im +
        combinedResponse_im * d2_re;
    }
    break;

   case 8:
    /*  Dec/Int3 */
    c_us(dec_int3_coeff_data, dec_int3_coeff_size, tmp_data, tmp_size);
    k_freqz_cg(tmp_data, tmp_size, w, Fs, combinedResponse);
    break;

   case 9:
    /*  Dec/Int3,Hb1 */
    memset(&c_u[0], 0, 30U * sizeof(double));
    ix = 0;
    iy = 0;
    for (k = 0; k < 15; k++) {
      c_u[iy] = hb1_coeff[ix];
      ix++;
      iy += 2;
    }

    l_freqz_cg(c_u, w, Fs, combinedResponse);
    c_us(dec_int3_coeff_data, dec_int3_coeff_size, tmp_data, tmp_size);
    k_freqz_cg(tmp_data, tmp_size, w, Fs, d2);
    ix = combinedResponse->size[0] * combinedResponse->size[1];
    combinedResponse->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)combinedResponse, ix, (int)sizeof
                      (creal_T));
    ix = combinedResponse->size[0];
    iy = combinedResponse->size[1];
    iy *= ix;
    for (ix = 0; ix < iy; ix++) {
      combinedResponse_re = combinedResponse->data[ix].re;
      combinedResponse_im = combinedResponse->data[ix].im;
      d2_re = d2->data[ix].re;
      d2_im = d2->data[ix].im;
      combinedResponse->data[ix].re = combinedResponse_re * d2_re -
        combinedResponse_im * d2_im;
      combinedResponse->data[ix].im = combinedResponse_re * d2_im +
        combinedResponse_im * d2_re;
    }
    break;

   case 10:
    /*  Dec/Int3,Hb2 */
    memset(&e_u[0], 0, 14U * sizeof(double));
    ix = 0;
    iy = 0;
    for (k = 0; k < 7; k++) {
      e_u[iy] = hb2_coeff[ix];
      ix++;
      iy += 2;
    }

    n_freqz_cg(e_u, w, Fs, combinedResponse);
    c_us(dec_int3_coeff_data, dec_int3_coeff_size, tmp_data, tmp_size);
    k_freqz_cg(tmp_data, tmp_size, w, Fs, d2);
    ix = combinedResponse->size[0] * combinedResponse->size[1];
    combinedResponse->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)combinedResponse, ix, (int)sizeof
                      (creal_T));
    ix = combinedResponse->size[0];
    iy = combinedResponse->size[1];
    iy *= ix;
    for (ix = 0; ix < iy; ix++) {
      combinedResponse_re = combinedResponse->data[ix].re;
      combinedResponse_im = combinedResponse->data[ix].im;
      d2_re = d2->data[ix].re;
      d2_im = d2->data[ix].im;
      combinedResponse->data[ix].re = combinedResponse_re * d2_re -
        combinedResponse_im * d2_im;
      combinedResponse->data[ix].im = combinedResponse_re * d2_im +
        combinedResponse_im * d2_re;
    }
    break;

   case 11:
    /*  Dec/Int3,Hb2,Hb1 */
    memset(&d_u[0], 0, 60U * sizeof(double));
    ix = 0;
    iy = 0;
    for (k = 0; k < 15; k++) {
      d_u[iy] = hb1_coeff[ix];
      ix++;
      iy += 4;
    }

    m_freqz_cg(d_u, w, Fs, combinedResponse);
    memset(&e_u[0], 0, 14U * sizeof(double));
    ix = 0;
    iy = 0;
    for (k = 0; k < 7; k++) {
      e_u[iy] = hb2_coeff[ix];
      ix++;
      iy += 2;
    }

    n_freqz_cg(e_u, w, Fs, d2);
    c_us(dec_int3_coeff_data, dec_int3_coeff_size, tmp_data, tmp_size);
    k_freqz_cg(tmp_data, tmp_size, w, Fs, d3);
    ix = combinedResponse->size[0] * combinedResponse->size[1];
    combinedResponse->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)combinedResponse, ix, (int)sizeof
                      (creal_T));
    ix = combinedResponse->size[0];
    iy = combinedResponse->size[1];
    iy *= ix;
    for (ix = 0; ix < iy; ix++) {
      combinedResponse_re = combinedResponse->data[ix].re * d2->data[ix].re -
        combinedResponse->data[ix].im * d2->data[ix].im;
      combinedResponse_im = combinedResponse->data[ix].re * d2->data[ix].im +
        combinedResponse->data[ix].im * d2->data[ix].re;
      d2_re = d3->data[ix].re;
      d2_im = d3->data[ix].im;
      combinedResponse->data[ix].re = combinedResponse_re * d2_re -
        combinedResponse_im * d2_im;
      combinedResponse->data[ix].im = combinedResponse_re * d2_im +
        combinedResponse_im * d2_re;
    }
    break;
  }

  emxFree_creal_T(&d3);
  emxFree_creal_T(&d2);

  /*  Add filter extra filter to end of cascade */
}

/*
 * Arguments    : double x
 * Return Type  : double
 */
static double b_log2(double x)
{
  double f;
  double t;
  int eint;
  if (x == 0.0) {
    f = rtMinusInf;
  } else if (x < 0.0) {
    f = rtNaN;
  } else if ((!rtIsInf(x)) && (!rtIsNaN(x))) {
    t = frexp(x, &eint);
    if (t == 0.5) {
      f = (double)eint - 1.0;
    } else {
      f = log(t) / 0.69314718055994529 + (double)eint;
    }
  } else {
    f = x;
  }

  return f;
}

/*
 * Arguments    : const double A[4]
 *                const double B[2]
 *                double Y[2]
 * Return Type  : void
 */
static void b_mldivide(const double A[4], const double B[2], double Y[2])
{
  int r1;
  int r2;
  double a21;
  if (fabs(A[1]) > fabs(A[0])) {
    r1 = 1;
    r2 = 0;
  } else {
    r1 = 0;
    r2 = 1;
  }

  a21 = A[r2] / A[r1];
  Y[1] = (B[r2] - B[r1] * a21) / (A[2 + r2] - a21 * A[2 + r1]);
  Y[0] = (B[r1] - Y[1] * A[2 + r1]) / A[r1];
}

/*
 * Arguments    : const double p[7]
 *                const creal_T x[2048]
 *                creal_T y[2048]
 * Return Type  : void
 */
static void b_polyval(const double p[7], const creal_T x[2048], creal_T y[2048])
{
  int i10;
  int k;
  double x_im;
  for (i10 = 0; i10 < 2048; i10++) {
    y[i10].re = p[0];
    y[i10].im = 0.0;
  }

  for (k = 0; k < 6; k++) {
    for (i10 = 0; i10 < 2048; i10++) {
      x_im = x[i10].re * y[i10].im + x[i10].im * y[i10].re;
      y[i10].re = (x[i10].re * y[i10].re - x[i10].im * y[i10].im) + p[k + 1];
      y[i10].im = x_im;
    }
  }
}

/*
 * Arguments    : const emxArray_real_T *a
 *                emxArray_real_T *y
 * Return Type  : void
 */
static void b_power(const emxArray_real_T *a, emxArray_real_T *y)
{
  int k;
  k = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = a->size[1];
  emxEnsureCapacity((emxArray__common *)y, k, (int)sizeof(double));
  for (k = 0; k + 1 <= a->size[1]; k++) {
    y->data[k] = rt_powd_snf(a->data[k], 3.0);
  }
}

/*
 * Arguments    : const emxArray_creal_T *x
 *                const emxArray_creal_T *y
 *                emxArray_creal_T *z
 * Return Type  : void
 */
static void b_rdivide(const emxArray_creal_T *x, const emxArray_creal_T *y,
                      emxArray_creal_T *z)
{
  int i25;
  int loop_ub;
  double x_re;
  double x_im;
  double y_re;
  double y_im;
  double brm;
  double bim;
  double s;
  i25 = z->size[0] * z->size[1];
  z->size[0] = 1;
  z->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)z, i25, (int)sizeof(creal_T));
  loop_ub = x->size[0] * x->size[1];
  for (i25 = 0; i25 < loop_ub; i25++) {
    x_re = x->data[i25].re;
    x_im = x->data[i25].im;
    y_re = y->data[i25].re;
    y_im = y->data[i25].im;
    if (y_im == 0.0) {
      if (x_im == 0.0) {
        z->data[i25].re = x_re / y_re;
        z->data[i25].im = 0.0;
      } else if (x_re == 0.0) {
        z->data[i25].re = 0.0;
        z->data[i25].im = x_im / y_re;
      } else {
        z->data[i25].re = x_re / y_re;
        z->data[i25].im = x_im / y_re;
      }
    } else if (y_re == 0.0) {
      if (x_re == 0.0) {
        z->data[i25].re = x_im / y_im;
        z->data[i25].im = 0.0;
      } else if (x_im == 0.0) {
        z->data[i25].re = 0.0;
        z->data[i25].im = -(x_re / y_im);
      } else {
        z->data[i25].re = x_im / y_im;
        z->data[i25].im = -(x_re / y_im);
      }
    } else {
      brm = fabs(y_re);
      bim = fabs(y_im);
      if (brm > bim) {
        s = y_im / y_re;
        bim = y_re + s * y_im;
        z->data[i25].re = (x_re + s * x_im) / bim;
        z->data[i25].im = (x_im - s * x_re) / bim;
      } else if (bim == brm) {
        if (y_re > 0.0) {
          s = 0.5;
        } else {
          s = -0.5;
        }

        if (y_im > 0.0) {
          bim = 0.5;
        } else {
          bim = -0.5;
        }

        z->data[i25].re = (x_re * s + x_im * bim) / brm;
        z->data[i25].im = (x_im * s - x_re * bim) / brm;
      } else {
        s = y_re / y_im;
        bim = y_im + s * y_re;
        z->data[i25].re = (s * x_re + x_im) / bim;
        z->data[i25].im = (s * x_im - x_re) / bim;
      }
    }
  }
}

/*
 * Arguments    : double x[2048]
 * Return Type  : void
 */
static void b_sin(double x[2048])
{
  int k;
  for (k = 0; k < 2048; k++) {
    x[k] = sin(x[k]);
  }
}

/*
 * Arguments    : emxArray_real_T *x
 *                emxArray_real_T *y
 * Return Type  : void
 */
static void b_sinc(emxArray_real_T *x, emxArray_real_T *y)
{
  emxArray_boolean_T *b_x;
  int i33;
  int loop_ub;
  emxArray_int32_T *ii;
  int nx;
  int idx;
  boolean_T exitg1;
  boolean_T guard1 = false;
  emxArray_int32_T *i;
  emxArray_real_T *c_x;
  emxArray_real_T *r19;
  emxInit_boolean_T(&b_x, 2);
  i33 = b_x->size[0] * b_x->size[1];
  b_x->size[0] = 1;
  b_x->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)b_x, i33, (int)sizeof(boolean_T));
  loop_ub = x->size[0] * x->size[1];
  for (i33 = 0; i33 < loop_ub; i33++) {
    b_x->data[i33] = (x->data[i33] == 0.0);
  }

  emxInit_int32_T(&ii, 2);
  nx = b_x->size[1];
  idx = 0;
  i33 = ii->size[0] * ii->size[1];
  ii->size[0] = 1;
  ii->size[1] = b_x->size[1];
  emxEnsureCapacity((emxArray__common *)ii, i33, (int)sizeof(int));
  loop_ub = 1;
  exitg1 = false;
  while ((!exitg1) && (loop_ub <= nx)) {
    guard1 = false;
    if (b_x->data[loop_ub - 1]) {
      idx++;
      ii->data[idx - 1] = loop_ub;
      if (idx >= nx) {
        exitg1 = true;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      loop_ub++;
    }
  }

  if (b_x->size[1] == 1) {
    if (idx == 0) {
      i33 = ii->size[0] * ii->size[1];
      ii->size[0] = 1;
      ii->size[1] = 0;
      emxEnsureCapacity((emxArray__common *)ii, i33, (int)sizeof(int));
    }
  } else {
    i33 = ii->size[0] * ii->size[1];
    if (1 > idx) {
      ii->size[1] = 0;
    } else {
      ii->size[1] = idx;
    }

    emxEnsureCapacity((emxArray__common *)ii, i33, (int)sizeof(int));
  }

  emxFree_boolean_T(&b_x);
  emxInit_int32_T(&i, 2);
  i33 = i->size[0] * i->size[1];
  i->size[0] = 1;
  i->size[1] = ii->size[1];
  emxEnsureCapacity((emxArray__common *)i, i33, (int)sizeof(int));
  loop_ub = ii->size[0] * ii->size[1];
  for (i33 = 0; i33 < loop_ub; i33++) {
    i->data[i33] = ii->data[i33];
  }

  i33 = ii->size[0] * ii->size[1];
  ii->size[0] = 1;
  ii->size[1] = i->size[1];
  emxEnsureCapacity((emxArray__common *)ii, i33, (int)sizeof(int));
  loop_ub = i->size[0] * i->size[1];
  for (i33 = 0; i33 < loop_ub; i33++) {
    ii->data[i33] = i->data[i33];
  }

  loop_ub = ii->size[0] * ii->size[1];
  for (i33 = 0; i33 < loop_ub; i33++) {
    x->data[ii->data[i33] - 1] = 1.0;
  }

  i33 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)y, i33, (int)sizeof(double));
  loop_ub = x->size[0] * x->size[1];
  for (i33 = 0; i33 < loop_ub; i33++) {
    y->data[i33] = 3.1415926535897931 * x->data[i33];
  }

  emxInit_real_T(&c_x, 2);
  i33 = c_x->size[0] * c_x->size[1];
  c_x->size[0] = 1;
  c_x->size[1] = y->size[1];
  emxEnsureCapacity((emxArray__common *)c_x, i33, (int)sizeof(double));
  loop_ub = y->size[0] * y->size[1];
  for (i33 = 0; i33 < loop_ub; i33++) {
    c_x->data[i33] = y->data[i33];
  }

  for (loop_ub = 0; loop_ub + 1 <= y->size[1]; loop_ub++) {
    c_x->data[loop_ub] = sin(c_x->data[loop_ub]);
  }

  emxInit_real_T(&r19, 2);
  i33 = r19->size[0] * r19->size[1];
  r19->size[0] = 1;
  r19->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)r19, i33, (int)sizeof(double));
  loop_ub = x->size[0] * x->size[1];
  for (i33 = 0; i33 < loop_ub; i33++) {
    r19->data[i33] = 3.1415926535897931 * x->data[i33];
  }

  c_rdivide(c_x, r19, y);
  i33 = ii->size[0] * ii->size[1];
  ii->size[0] = 1;
  ii->size[1] = i->size[1];
  emxEnsureCapacity((emxArray__common *)ii, i33, (int)sizeof(int));
  loop_ub = i->size[0] * i->size[1];
  emxFree_real_T(&r19);
  emxFree_real_T(&c_x);
  for (i33 = 0; i33 < loop_ub; i33++) {
    ii->data[i33] = i->data[i33];
  }

  emxFree_int32_T(&i);
  loop_ub = ii->size[0] * ii->size[1];
  for (i33 = 0; i33 < loop_ub; i33++) {
    y->data[ii->data[i33] - 1] = 1.0;
  }

  emxFree_int32_T(&ii);
}

/*
 * Arguments    : creal_T *x
 * Return Type  : void
 */
static void b_sqrt(creal_T *x)
{
  double absxi;
  double absxr;
  if (x->im == 0.0) {
    if (x->re < 0.0) {
      absxi = 0.0;
      absxr = sqrt(fabs(x->re));
    } else {
      absxi = sqrt(x->re);
      absxr = 0.0;
    }
  } else if (x->re == 0.0) {
    if (x->im < 0.0) {
      absxi = sqrt(-x->im / 2.0);
      absxr = -absxi;
    } else {
      absxi = sqrt(x->im / 2.0);
      absxr = absxi;
    }
  } else if (rtIsNaN(x->re) || rtIsNaN(x->im)) {
    absxi = rtNaN;
    absxr = rtNaN;
  } else if (rtIsInf(x->im)) {
    absxi = rtInf;
    absxr = x->im;
  } else if (rtIsInf(x->re)) {
    if (x->re < 0.0) {
      absxi = 0.0;
      absxr = rtInf;
    } else {
      absxi = rtInf;
      absxr = 0.0;
    }
  } else {
    absxr = fabs(x->re);
    absxi = fabs(x->im);
    if ((absxr > 4.4942328371557893E+307) || (absxi > 4.4942328371557893E+307))
    {
      absxr *= 0.5;
      absxi *= 0.5;
      absxi = rt_hypotd_snf(absxr, absxi);
      if (absxi > absxr) {
        absxi = sqrt(absxi) * sqrt(1.0 + absxr / absxi);
      } else {
        absxi = sqrt(absxi) * 1.4142135623730951;
      }
    } else {
      absxi = sqrt((rt_hypotd_snf(absxr, absxi) + absxr) * 0.5);
    }

    if (x->re > 0.0) {
      absxr = 0.5 * (x->im / absxi);
    } else {
      if (x->im < 0.0) {
        absxr = -absxi;
      } else {
        absxr = absxi;
      }

      absxi = 0.5 * (x->im / absxr);
    }
  }

  x->re = absxi;
  x->im = absxr;
}

/*
 * Arguments    : const char a[2]
 * Return Type  : boolean_T
 */
static boolean_T b_strcmp(const char a[2])
{
  boolean_T b_bool;
  int kstr;
  int exitg1;
  static const char cv0[2] = { 'R', 'x' };

  b_bool = false;
  kstr = 0;
  do {
    exitg1 = 0;
    if (kstr + 1 < 3) {
      if (a[kstr] != cv0[kstr]) {
        exitg1 = 1;
      } else {
        kstr++;
      }
    } else {
      b_bool = true;
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  return b_bool;
}

/*
 * Arguments    : const emxArray_real_T *x
 * Return Type  : double
 */
static double b_sum(const emxArray_real_T *x)
{
  double y;
  int k;
  if (x->size[1] == 0) {
    y = 0.0;
  } else {
    y = x->data[0];
    for (k = 2; k <= x->size[1]; k++) {
      y += x->data[k - 1];
    }
  }

  return y;
}

/*
 * Arguments    : const double o[7]
 *                double u[7]
 * Return Type  : void
 */
static void b_us(const double o[7], double u[7])
{
  int ix;
  int iy;
  int k;
  for (ix = 0; ix < 7; ix++) {
    u[ix] = 0.0;
  }

  ix = 0;
  iy = 0;
  for (k = 0; k < 7; k++) {
    u[iy] = o[ix];
    ix++;
    iy++;
  }
}

/*
 * Arguments    : int n
 *                const creal_T a
 *                emxArray_creal_T *x
 *                int ix0
 *                int incx
 * Return Type  : void
 */
static void b_xscal(int n, const creal_T a, emxArray_creal_T *x, int ix0, int
                    incx)
{
  double a_re;
  double a_im;
  int i53;
  int k;
  double x_re;
  double x_im;
  a_re = a.re;
  a_im = a.im;
  if (!(incx < 1)) {
    i53 = ix0 + incx * (n - 1);
    for (k = ix0; k <= i53; k += incx) {
      x_re = x->data[k - 1].re;
      x_im = x->data[k - 1].im;
      x->data[k - 1].re = a_re * x_re - a_im * x_im;
      x->data[k - 1].im = a_re * x_im + a_im * x_re;
    }
  }
}

/*
 * Arguments    : const creal_T f
 *                const creal_T g
 *                double *cs
 *                creal_T *sn
 * Return Type  : void
 */
static void b_xzlartg(const creal_T f, const creal_T g, double *cs, creal_T *sn)
{
  double scale;
  double g2;
  double f2s;
  double fs_re;
  double fs_im;
  double gs_re;
  double gs_im;
  boolean_T guard1 = false;
  double g2s;
  double d;
  scale = fabs(f.re);
  g2 = fabs(f.im);
  if (g2 > scale) {
    scale = g2;
  }

  f2s = fabs(g.re);
  g2 = fabs(g.im);
  if (g2 > f2s) {
    f2s = g2;
  }

  if (f2s > scale) {
    scale = f2s;
  }

  fs_re = f.re;
  fs_im = f.im;
  gs_re = g.re;
  gs_im = g.im;
  guard1 = false;
  if (scale >= 7.4428285367870146E+137) {
    do {
      fs_re *= 1.3435752215134178E-138;
      fs_im *= 1.3435752215134178E-138;
      gs_re *= 1.3435752215134178E-138;
      gs_im *= 1.3435752215134178E-138;
      scale *= 1.3435752215134178E-138;
    } while (!(scale < 7.4428285367870146E+137));

    guard1 = true;
  } else if (scale <= 1.3435752215134178E-138) {
    if ((g.re == 0.0) && (g.im == 0.0)) {
      *cs = 1.0;
      sn->re = 0.0;
      sn->im = 0.0;
    } else {
      do {
        fs_re *= 7.4428285367870146E+137;
        fs_im *= 7.4428285367870146E+137;
        gs_re *= 7.4428285367870146E+137;
        gs_im *= 7.4428285367870146E+137;
        scale *= 7.4428285367870146E+137;
      } while (!(scale > 1.3435752215134178E-138));

      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1) {
    scale = fs_re * fs_re + fs_im * fs_im;
    g2 = gs_re * gs_re + gs_im * gs_im;
    f2s = g2;
    if (1.0 > g2) {
      f2s = 1.0;
    }

    if (scale <= f2s * 2.0041683600089728E-292) {
      if ((f.re == 0.0) && (f.im == 0.0)) {
        *cs = 0.0;
        d = rt_hypotd_snf(gs_re, gs_im);
        sn->re = gs_re / d;
        sn->im = -gs_im / d;
      } else {
        g2s = sqrt(g2);
        *cs = rt_hypotd_snf(fs_re, fs_im) / g2s;
        f2s = fabs(f.re);
        g2 = fabs(f.im);
        if (g2 > f2s) {
          f2s = g2;
        }

        if (f2s > 1.0) {
          d = rt_hypotd_snf(f.re, f.im);
          fs_re = f.re / d;
          fs_im = f.im / d;
        } else {
          g2 = 7.4428285367870146E+137 * f.re;
          scale = 7.4428285367870146E+137 * f.im;
          d = rt_hypotd_snf(g2, scale);
          fs_re = g2 / d;
          fs_im = scale / d;
        }

        gs_re /= g2s;
        gs_im = -gs_im / g2s;
        sn->re = fs_re * gs_re - fs_im * gs_im;
        sn->im = fs_re * gs_im + fs_im * gs_re;
      }
    } else {
      f2s = sqrt(1.0 + g2 / scale);
      *cs = 1.0 / f2s;
      d = scale + g2;
      fs_re = f2s * fs_re / d;
      fs_im = f2s * fs_im / d;
      sn->re = fs_re * gs_re - fs_im * -gs_im;
      sn->im = fs_re * -gs_im + fs_im * gs_re;
    }
  }
}

/*
 * BUTTER_CG Butterworth digital and analog filter design.  Codegen support
 *
 *  This function is based on 'butter' by The MathWorks Inc.
 * Arguments    : double Wn
 *                double num[2]
 *                creal_T den_data[]
 *                int den_size[2]
 * Return Type  : void
 */
static void butter_cg(double Wn, double num[2], creal_T den_data[], int
                      den_size[2])
{
  emxArray_creal_T *c;
  int k;
  creal_T tmp_data[1];
  int tmp_size[2];
  double b_tmp_data[1];
  int b_tmp_size[1];
  emxArray_creal_T *a;
  emxArray_real_T *b;
  emxArray_creal_T c_tmp_data;
  emxArray_real_T d_tmp_data;
  double d;
  int loop_ub;
  double y_re;
  double y_im;
  double ai;
  emxInit_creal_T(&c, 2);

  /*  Set types we will only use */
  /*  Cast to enforce precision rules */
  /*  step 1: get analog, pre-warped frequencies */
  /*  step 2: convert to low-pass prototype estimate */
  /*  lowpass */
  /*  Cast to enforce precision rules */
  /*  step 3: Get N-th order Butterworth analog lowpass prototype */
  /*  Transform to state-space */
  /* ZP2SS  Zero-pole to state-space conversion. Codegen support */
  /*  */
  /*  This function is based on 'zp2ss' by The MathWorks Inc. */
  /*  Strip infinities and throw away. */
  /*  Group into complex pairs */
  /*  try */
  /*      % z and p should have real elements and exact complex conjugate pair. */
  /*      z = cplxpair(zF,0); */
  /*      p = cplxpair(pF,0); */
  /*  catch */
  /*      % If fail, revert to use the old default tolerance. */
  /*      % The use of tolerance in checking for real entries and conjugate pairs */
  /*      % may result in misinterpretation for edge cases. Please review the */
  /*      % process of how z and p are generated. */
  /*      z = cplxpair(zF,1e6*nz*norm(zF)*eps + eps); */
  /*      p = cplxpair(pF,1e6*np*norm(pF)*eps + eps); */
  /*  end */
  /*  Initialize state-space matrices for running series */
  /*  If odd number of poles AND zeros, convert the pole and zero */
  /*  at the end into state-space. */
  /*    H(s) = (s-z1)/(s-p1) = (s + num(2)) / (s + den(2)) */
  /*  If odd number of poles only, convert the pole at the */
  /*  end into state-space. */
  /*   H(s) = 1/(s-p1) = 1/(s + den(2)) */
  /*  If odd number of zeros only, convert the zero at the */
  /*  end, along with a pole-pair into state-space. */
  /*    H(s) = (s+num(2))/(s^2+den(2)s+den(3)) */
  /*  Now we have an even number of poles and zeros, although not */
  /*  necessarily the same number - there may be more poles. */
  /*    H(s) = (s^2+num(2)s+num(3))/(s^2+den(2)s+den(3)) */
  /*  Loop through rest of pairs, connecting in series to build the model. */
  /*  Take care of any left over unmatched pole pairs. */
  /*    H(s) = 1/(s^2+den(2)s+den(3)) */
  /*  Apply gain k: */
  /*  step 4: Transform to lowpass, bandpass, highpass, or bandstop of desired Wn */
  /*  Lowpass */
  k = c->size[0] * c->size[1];
  c->size[0] = 1;
  c->size[1] = 1;
  emxEnsureCapacity((emxArray__common *)c, k, (int)sizeof(creal_T));
  c->data[0].re = 1.0;
  c->data[0].im = 0.0;
  tmp_size[0] = 1;
  tmp_size[1] = 1;
  tmp_data[0].re = -1.0;
  tmp_data[0].im = 0.0;
  b_tmp_size[0] = 1;
  b_tmp_data[0] = 1.0;
  emxInit_creal_T(&a, 2);
  emxInit_real_T1(&b, 1);
  c_tmp_data.data = (creal_T *)&tmp_data;
  c_tmp_data.size = (int *)&tmp_size;
  c_tmp_data.allocatedSize = 1;
  c_tmp_data.numDimensions = 2;
  c_tmp_data.canFreeData = false;
  d_tmp_data.data = (double *)&b_tmp_data;
  d_tmp_data.size = (int *)&b_tmp_size;
  d_tmp_data.allocatedSize = 1;
  d_tmp_data.numDimensions = 1;
  d_tmp_data.canFreeData = false;
  lp2lp_cg(&c_tmp_data, &d_tmp_data, 0.0, Wn, a, b, &d);

  /*  step 5: Use Bilinear transformation to find discrete equivalent: */
  /*  nargout <= 3 */
  /*  Transform to zero-pole-gain and polynomial forms: */
  /*  nargout <= 2 */
  poly(a, c);
  den_size[0] = 1;
  den_size[1] = c->size[1];
  loop_ub = c->size[0] * c->size[1];
  emxFree_real_T(&b);
  emxFree_creal_T(&a);
  for (k = 0; k < loop_ub; k++) {
    den_data[k] = c->data[k];
  }

  emxFree_creal_T(&c);

  /*  This internal function returns more exact numerator vectors */
  /*  for the num/den case. */
  /*  Wn input is two element band edge vector */
  /* --------------------------------- */
  /*  lowpass */
  y_re = den_data[0].re;
  y_im = den_data[0].im;
  for (k = 0; k <= den_size[1] - 2; k++) {
    d = 0.0 * y_im + 0.0 * y_re;
    y_re = (0.0 * y_re - 0.0 * y_im) + den_data[k + 1].re;
    y_im = d + den_data[k + 1].im;
  }

  for (k = 0; k < 2; k++) {
    d = (double)k * y_re;
    ai = (double)k * y_im;
    if ((!(ai == 0.0)) && (d == 0.0)) {
      d = 0.0;
    }

    num[k] = d;
  }

  /*  num = poly(a-b*c)+(d-1)*den; */
}

/*
 * Arguments    : const emxArray_real_T *x
 *                emxArray_real_T *y
 * Return Type  : void
 */
static void c_abs(const emxArray_real_T *x, emxArray_real_T *y)
{
  int k;
  k = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)y, k, (int)sizeof(double));
  for (k = 0; k + 1 <= x->size[1]; k++) {
    y->data[k] = fabs(x->data[k]);
  }
}

/*
 * UNTITLED Summary of this function goes here
 *    Detailed explanation goes here
 * Arguments    : const emxArray_real_T *f
 *                double Fconverter
 *                const double b1_data[]
 *                const int b1_size[2]
 *                const emxArray_creal_T *a1
 *                const double b2_data[]
 *                const int b2_size[2]
 *                const emxArray_creal_T *a2
 *                emxArray_creal_T *abc
 * Return Type  : void
 */
static void c_analogresp(const emxArray_real_T *f, double Fconverter, const
  double b1_data[], const int b1_size[2], const emxArray_creal_T *a1, const
  double b2_data[], const int b2_size[2], const emxArray_creal_T *a2,
  emxArray_creal_T *abc)
{
  emxArray_real_T *r27;
  int b_abc;
  int loop_ub;
  emxArray_real_T *r28;
  emxArray_creal_T *r29;
  emxArray_real_T *r30;
  emxArray_real_T *r31;
  double abc_re;
  double abc_im;
  emxInit_real_T(&r27, 2);
  b_abc = r27->size[0] * r27->size[1];
  r27->size[0] = 1;
  r27->size[1] = f->size[1];
  emxEnsureCapacity((emxArray__common *)r27, b_abc, (int)sizeof(double));
  loop_ub = f->size[0] * f->size[1];
  for (b_abc = 0; b_abc < loop_ub; b_abc++) {
    r27->data[b_abc] = 6.2831853071795862 * f->data[b_abc];
  }

  emxInit_real_T(&r28, 2);
  b_freqs_cg(b1_data, b1_size, a1, r27, abc);
  b_abc = r28->size[0] * r28->size[1];
  r28->size[0] = 1;
  r28->size[1] = f->size[1];
  emxEnsureCapacity((emxArray__common *)r28, b_abc, (int)sizeof(double));
  loop_ub = f->size[0] * f->size[1];
  emxFree_real_T(&r27);
  for (b_abc = 0; b_abc < loop_ub; b_abc++) {
    r28->data[b_abc] = 6.2831853071795862 * f->data[b_abc];
  }

  emxInit_creal_T(&r29, 2);
  emxInit_real_T(&r30, 2);
  emxInit_real_T(&r31, 2);
  b_freqs_cg(b2_data, b2_size, a2, r28, r29);
  rdivide(f, Fconverter, r30);
  b_sinc(r30, r31);
  b_power(r31, r30);
  b_abc = abc->size[0] * abc->size[1];
  abc->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)abc, b_abc, (int)sizeof(creal_T));
  b_abc = abc->size[0];
  loop_ub = abc->size[1];
  loop_ub *= b_abc;
  emxFree_real_T(&r28);
  emxFree_real_T(&r31);
  for (b_abc = 0; b_abc < loop_ub; b_abc++) {
    abc_re = abc->data[b_abc].re * r29->data[b_abc].re - abc->data[b_abc].im *
      r29->data[b_abc].im;
    abc_im = abc->data[b_abc].re * r29->data[b_abc].im + abc->data[b_abc].im *
      r29->data[b_abc].re;
    abc->data[b_abc].re = r30->data[b_abc] * abc_re;
    abc->data[b_abc].im = r30->data[b_abc] * abc_im;
  }

  emxFree_real_T(&r30);
  emxFree_creal_T(&r29);
}

/*
 * Arguments    : emxArray_creal_T *x
 * Return Type  : void
 */
static void c_exp(emxArray_creal_T *x)
{
  int nx;
  int k;
  double x_re;
  double r;
  nx = x->size[1];
  for (k = 0; k + 1 <= nx; k++) {
    if (x->data[k].im == 0.0) {
      x_re = exp(x->data[k].re);
      r = 0.0;
    } else if (rtIsInf(x->data[k].im) && rtIsInf(x->data[k].re) && (x->data[k].
                re < 0.0)) {
      x_re = 0.0;
      r = 0.0;
    } else {
      r = exp(x->data[k].re / 2.0);
      x_re = r * (r * cos(x->data[k].im));
      r *= r * sin(x->data[k].im);
    }

    x->data[k].re = x_re;
    x->data[k].im = r;
  }
}

/*
 * Make b a row
 * Arguments    : const double b[7]
 *                const struct_T *options
 *                creal_T h[2048]
 *                double w[2048]
 * Return Type  : void
 */
static void c_firfreqz(const double b[7], const struct_T *options, creal_T h
  [2048], double w[2048])
{
  double digw[2048];
  creal_T dcv3[2048];
  int i56;
  double bim;
  double h_re;
  double brm;
  double d;

  /* -------------------------------------------------------------------------- */
  /*  Actual Frequency Response Computation */
  /* if fvflag, */
  /*    Frequency vector specified.  Use Horner's method of polynomial */
  /*    evaluation at the frequency points and divide the numerator */
  /*    by the denominator. */
  /*  */
  /*    Note: we use positive i here because of the relationship */
  /*             polyval(a,exp(1i*w)) = fft(a).*exp(1i*w*(length(a)-1)) */
  /*                ( assuming w = 2*pi*(0:length(a)-1)/length(a) ) */
  /*  */
  /*  Fs was specified, freq. vector is in Hz */
  /*  Convert from Hz to rad/sample for computational purposes */
  /*  Digital frequency must be used for this calculation */
  for (i56 = 0; i56 < 2048; i56++) {
    w[i56] = options->w[i56];
    bim = 6.2831853071795862 * options->w[i56] / options->Fs;
    dcv3[i56].re = bim * 0.0;
    dcv3[i56].im = bim;
    digw[i56] = bim;
  }

  b_exp(dcv3);
  b_polyval(b, dcv3, h);
  for (i56 = 0; i56 < 2048; i56++) {
    dcv3[i56].re = 6.0 * (digw[i56] * 0.0);
    dcv3[i56].im = 6.0 * digw[i56];
  }

  b_exp(dcv3);
  for (i56 = 0; i56 < 2048; i56++) {
    h_re = h[i56].re;
    if (dcv3[i56].im == 0.0) {
      if (h[i56].im == 0.0) {
        h[i56].re /= dcv3[i56].re;
        h[i56].im = 0.0;
      } else if (h[i56].re == 0.0) {
        h[i56].re = 0.0;
        h[i56].im /= dcv3[i56].re;
      } else {
        h[i56].re /= dcv3[i56].re;
        h[i56].im /= dcv3[i56].re;
      }
    } else if (dcv3[i56].re == 0.0) {
      if (h[i56].re == 0.0) {
        h[i56].re = h[i56].im / dcv3[i56].im;
        h[i56].im = 0.0;
      } else if (h[i56].im == 0.0) {
        h[i56].re = 0.0;
        h[i56].im = -(h_re / dcv3[i56].im);
      } else {
        h[i56].re = h[i56].im / dcv3[i56].im;
        h[i56].im = -(h_re / dcv3[i56].im);
      }
    } else {
      brm = fabs(dcv3[i56].re);
      bim = fabs(dcv3[i56].im);
      if (brm > bim) {
        bim = dcv3[i56].im / dcv3[i56].re;
        d = dcv3[i56].re + bim * dcv3[i56].im;
        h[i56].re = (h[i56].re + bim * h[i56].im) / d;
        h[i56].im = (h[i56].im - bim * h_re) / d;
      } else if (bim == brm) {
        if (dcv3[i56].re > 0.0) {
          bim = 0.5;
        } else {
          bim = -0.5;
        }

        if (dcv3[i56].im > 0.0) {
          d = 0.5;
        } else {
          d = -0.5;
        }

        h[i56].re = (h[i56].re * bim + h[i56].im * d) / brm;
        h[i56].im = (h[i56].im * bim - h_re * d) / brm;
      } else {
        bim = dcv3[i56].re / dcv3[i56].im;
        d = dcv3[i56].im + bim * dcv3[i56].re;
        h[i56].re = (bim * h[i56].re + h[i56].im) / d;
        h[i56].im = (bim * h[i56].im - h_re) / d;
      }
    }
  }
}

/*
 * Arguments    : emxArray_real_T *x
 * Return Type  : void
 */
static void c_fix(emxArray_real_T *x)
{
  int nx;
  int k;
  double b_x;
  nx = x->size[1];
  for (k = 0; k + 1 <= nx; k++) {
    if (x->data[k] < 0.0) {
      b_x = ceil(x->data[k]);
    } else {
      b_x = floor(x->data[k]);
    }

    x->data[k] = b_x;
  }
}

/*
 * FREQZ_CG Frequency response of digital filter with codegen support
 *
 *  This function is based on 'freqz' by The MathWorks Inc.
 * Arguments    : const double b[7]
 *                const double w[2048]
 *                double Fs
 *                creal_T hh[2048]
 * Return Type  : void
 */
static void c_freqz_cg(const double b[7], const double w[2048], double Fs,
  creal_T hh[2048])
{
  static struct_T options;
  int i11;
  static const char cv16[8] = { 'o', 'n', 'e', 's', 'i', 'd', 'e', 'd' };

  double b_w[2048];
  static const char cv17[7] = { 'o', 'm', 'i', 't', 't', 'e', 'd' };

  /*  Cast to enforce precision rules */
  options.Fs = Fs;
  memcpy(&options.nfft[0], &w[0], sizeof(double) << 11);
  memcpy(&options.w[0], &w[0], sizeof(double) << 11);

  /*  Remaining are default or for advanced use */
  options.fvflag = 1.0;
  for (i11 = 0; i11 < 8; i11++) {
    options.range[i11] = cv16[i11];
  }

  options.centerdc = 0.0;
  for (i11 = 0; i11 < 7; i11++) {
    options.configlevel[i11] = cv17[i11];
  }

  c_firfreqz(b, &options, hh, b_w);

  /*  Generate the default structure to pass to freqzplot */
  /*  If rad/sample, Fs is empty */
}

/*
 * Arguments    : const char enables[4]
 *                const emxArray_real_T *w
 *                double Fs
 *                const double hb1_coeff[15]
 *                const double hb2_coeff[7]
 *                const double hb3_coeff_data[]
 *                const int hb3_coeff_size[2]
 *                const double dec_int3_coeff_data[]
 *                const int dec_int3_coeff_size[2]
 *                const double extraTaps_data[]
 *                const int extraTaps_size[2]
 *                emxArray_creal_T *combinedResponse
 * Return Type  : void
 */
static void c_generateCascadedResponseRx(const char enables[4], const
  emxArray_real_T *w, double Fs, const double hb1_coeff[15], const double
  hb2_coeff[7], const double hb3_coeff_data[], const int hb3_coeff_size[2],
  const double dec_int3_coeff_data[], const int dec_int3_coeff_size[2], const
  double extraTaps_data[], const int extraTaps_size[2], emxArray_creal_T
  *combinedResponse)
{
  boolean_T b_bool;
  int ix;
  int exitg12;
  emxArray_creal_T *d2;
  int exitg11;
  emxArray_creal_T *d3;
  static const char cv42[4] = { '2', '1', '1', '1' };

  double u[15];
  double b_u[7];
  double tmp_data[29];
  int tmp_size[2];
  double c_u[30];
  double d_u[60];
  double e_u[14];
  int stages;
  int c;
  int k;
  int u_size[2];
  int exitg10;
  static const char cv43[4] = { '1', '2', '1', '1' };

  double u_data[1024];
  int exitg9;
  double combinedResponse_re;
  double combinedResponse_im;
  static const char cv44[4] = { '1', '1', '2', '1' };

  double d2_re;
  double d2_im;
  int exitg8;
  static const char cv45[4] = { '2', '2', '1', '1' };

  int exitg7;
  static const char cv46[4] = { '2', '1', '2', '1' };

  int exitg6;
  static const char cv47[4] = { '1', '2', '2', '1' };

  int exitg5;
  static const char cv48[4] = { '2', '2', '2', '1' };

  int exitg4;
  static const char cv49[4] = { '1', '1', '1', '3' };

  int exitg3;
  static const char cv50[4] = { '2', '1', '1', '3' };

  int exitg2;
  static const char cv51[4] = { '1', '2', '1', '3' };

  int exitg1;
  static const char cv52[4] = { '2', '2', '1', '3' };

  /*  Cast */
  b_bool = false;
  ix = 1;
  do {
    exitg12 = 0;
    if (ix < 5) {
      if (enables[ix - 1] != '1') {
        exitg12 = 1;
      } else {
        ix++;
      }
    } else {
      b_bool = true;
      exitg12 = 1;
    }
  } while (exitg12 == 0);

  if (b_bool) {
    ix = 0;
  } else {
    b_bool = false;
    ix = 0;
    do {
      exitg11 = 0;
      if (ix + 1 < 5) {
        if (enables[ix] != cv42[ix]) {
          exitg11 = 1;
        } else {
          ix++;
        }
      } else {
        b_bool = true;
        exitg11 = 1;
      }
    } while (exitg11 == 0);

    if (b_bool) {
      ix = 1;
    } else {
      b_bool = false;
      ix = 0;
      do {
        exitg10 = 0;
        if (ix + 1 < 5) {
          if (enables[ix] != cv43[ix]) {
            exitg10 = 1;
          } else {
            ix++;
          }
        } else {
          b_bool = true;
          exitg10 = 1;
        }
      } while (exitg10 == 0);

      if (b_bool) {
        ix = 2;
      } else {
        b_bool = false;
        ix = 0;
        do {
          exitg9 = 0;
          if (ix + 1 < 5) {
            if (enables[ix] != cv44[ix]) {
              exitg9 = 1;
            } else {
              ix++;
            }
          } else {
            b_bool = true;
            exitg9 = 1;
          }
        } while (exitg9 == 0);

        if (b_bool) {
          ix = 3;
        } else {
          b_bool = false;
          ix = 0;
          do {
            exitg8 = 0;
            if (ix + 1 < 5) {
              if (enables[ix] != cv45[ix]) {
                exitg8 = 1;
              } else {
                ix++;
              }
            } else {
              b_bool = true;
              exitg8 = 1;
            }
          } while (exitg8 == 0);

          if (b_bool) {
            ix = 4;
          } else {
            b_bool = false;
            ix = 0;
            do {
              exitg7 = 0;
              if (ix + 1 < 5) {
                if (enables[ix] != cv46[ix]) {
                  exitg7 = 1;
                } else {
                  ix++;
                }
              } else {
                b_bool = true;
                exitg7 = 1;
              }
            } while (exitg7 == 0);

            if (b_bool) {
              ix = 5;
            } else {
              b_bool = false;
              ix = 0;
              do {
                exitg6 = 0;
                if (ix + 1 < 5) {
                  if (enables[ix] != cv47[ix]) {
                    exitg6 = 1;
                  } else {
                    ix++;
                  }
                } else {
                  b_bool = true;
                  exitg6 = 1;
                }
              } while (exitg6 == 0);

              if (b_bool) {
                ix = 6;
              } else {
                b_bool = false;
                ix = 0;
                do {
                  exitg5 = 0;
                  if (ix + 1 < 5) {
                    if (enables[ix] != cv48[ix]) {
                      exitg5 = 1;
                    } else {
                      ix++;
                    }
                  } else {
                    b_bool = true;
                    exitg5 = 1;
                  }
                } while (exitg5 == 0);

                if (b_bool) {
                  ix = 7;
                } else {
                  b_bool = false;
                  ix = 0;
                  do {
                    exitg4 = 0;
                    if (ix + 1 < 5) {
                      if (enables[ix] != cv49[ix]) {
                        exitg4 = 1;
                      } else {
                        ix++;
                      }
                    } else {
                      b_bool = true;
                      exitg4 = 1;
                    }
                  } while (exitg4 == 0);

                  if (b_bool) {
                    ix = 8;
                  } else {
                    b_bool = false;
                    ix = 0;
                    do {
                      exitg3 = 0;
                      if (ix + 1 < 5) {
                        if (enables[ix] != cv50[ix]) {
                          exitg3 = 1;
                        } else {
                          ix++;
                        }
                      } else {
                        b_bool = true;
                        exitg3 = 1;
                      }
                    } while (exitg3 == 0);

                    if (b_bool) {
                      ix = 9;
                    } else {
                      b_bool = false;
                      ix = 0;
                      do {
                        exitg2 = 0;
                        if (ix + 1 < 5) {
                          if (enables[ix] != cv51[ix]) {
                            exitg2 = 1;
                          } else {
                            ix++;
                          }
                        } else {
                          b_bool = true;
                          exitg2 = 1;
                        }
                      } while (exitg2 == 0);

                      if (b_bool) {
                        ix = 10;
                      } else {
                        b_bool = false;
                        ix = 0;
                        do {
                          exitg1 = 0;
                          if (ix + 1 < 5) {
                            if (enables[ix] != cv52[ix]) {
                              exitg1 = 1;
                            } else {
                              ix++;
                            }
                          } else {
                            b_bool = true;
                            exitg1 = 1;
                          }
                        } while (exitg1 == 0);

                        if (b_bool) {
                          ix = 11;
                        } else {
                          ix = -1;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  emxInit_creal_T(&d2, 2);
  emxInit_creal_T(&d3, 2);
  switch (ix) {
   case 0:
    /*  only FIR */
    h_freqz_cg(w, Fs, combinedResponse);
    stages = 1;
    break;

   case 1:
    /*  Hb1 */
    memset(&u[0], 0, 15U * sizeof(double));
    ix = 0;
    stages = 0;
    for (k = 0; k < 15; k++) {
      u[stages] = hb1_coeff[ix];
      ix++;
      stages++;
    }

    i_freqz_cg(u, w, Fs, combinedResponse);
    stages = 1;
    break;

   case 2:
    /*  Hb2 */
    for (ix = 0; ix < 7; ix++) {
      b_u[ix] = 0.0;
    }

    ix = 0;
    stages = 0;
    for (k = 0; k < 7; k++) {
      b_u[stages] = hb2_coeff[ix];
      ix++;
      stages++;
    }

    j_freqz_cg(b_u, w, Fs, combinedResponse);
    stages = 1;
    break;

   case 3:
    /*  Hb3 */
    c_us(hb3_coeff_data, hb3_coeff_size, tmp_data, tmp_size);
    k_freqz_cg(tmp_data, tmp_size, w, Fs, combinedResponse);
    stages = 1;
    break;

   case 4:
    /*  Hb2,Hb1 */
    memset(&c_u[0], 0, 30U * sizeof(double));
    ix = 0;
    stages = 0;
    for (k = 0; k < 15; k++) {
      c_u[stages] = hb1_coeff[ix];
      ix++;
      stages += 2;
    }

    l_freqz_cg(c_u, w, Fs, combinedResponse);
    for (ix = 0; ix < 7; ix++) {
      b_u[ix] = 0.0;
    }

    ix = 0;
    stages = 0;
    for (k = 0; k < 7; k++) {
      b_u[stages] = hb2_coeff[ix];
      ix++;
      stages++;
    }

    j_freqz_cg(b_u, w, Fs, d2);
    ix = combinedResponse->size[0] * combinedResponse->size[1];
    combinedResponse->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)combinedResponse, ix, (int)sizeof
                      (creal_T));
    ix = combinedResponse->size[0];
    stages = combinedResponse->size[1];
    stages *= ix;
    for (ix = 0; ix < stages; ix++) {
      combinedResponse_re = combinedResponse->data[ix].re;
      combinedResponse_im = combinedResponse->data[ix].im;
      d2_re = d2->data[ix].re;
      d2_im = d2->data[ix].im;
      combinedResponse->data[ix].re = combinedResponse_re * d2_re -
        combinedResponse_im * d2_im;
      combinedResponse->data[ix].im = combinedResponse_re * d2_im +
        combinedResponse_im * d2_re;
    }

    stages = 2;
    break;

   case 5:
    /*  Hb3,Hb1 */
    memset(&c_u[0], 0, 30U * sizeof(double));
    ix = 0;
    stages = 0;
    for (k = 0; k < 15; k++) {
      c_u[stages] = hb1_coeff[ix];
      ix++;
      stages += 2;
    }

    l_freqz_cg(c_u, w, Fs, combinedResponse);
    for (ix = 0; ix < 7; ix++) {
      b_u[ix] = 0.0;
    }

    ix = 0;
    stages = 0;
    for (k = 0; k < 7; k++) {
      b_u[stages] = hb2_coeff[ix];
      ix++;
      stages++;
    }

    j_freqz_cg(b_u, w, Fs, d2);
    ix = combinedResponse->size[0] * combinedResponse->size[1];
    combinedResponse->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)combinedResponse, ix, (int)sizeof
                      (creal_T));
    ix = combinedResponse->size[0];
    stages = combinedResponse->size[1];
    stages *= ix;
    for (ix = 0; ix < stages; ix++) {
      combinedResponse_re = combinedResponse->data[ix].re;
      combinedResponse_im = combinedResponse->data[ix].im;
      d2_re = d2->data[ix].re;
      d2_im = d2->data[ix].im;
      combinedResponse->data[ix].re = combinedResponse_re * d2_re -
        combinedResponse_im * d2_im;
      combinedResponse->data[ix].im = combinedResponse_re * d2_im +
        combinedResponse_im * d2_re;
    }

    stages = 2;
    break;

   case 6:
    /*  Hb3,Hb2 */
    memset(&c_u[0], 0, 30U * sizeof(double));
    ix = 0;
    stages = 0;
    for (k = 0; k < 15; k++) {
      c_u[stages] = hb1_coeff[ix];
      ix++;
      stages += 2;
    }

    l_freqz_cg(c_u, w, Fs, combinedResponse);
    c_us(hb3_coeff_data, hb3_coeff_size, tmp_data, tmp_size);
    k_freqz_cg(tmp_data, tmp_size, w, Fs, d2);
    ix = combinedResponse->size[0] * combinedResponse->size[1];
    combinedResponse->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)combinedResponse, ix, (int)sizeof
                      (creal_T));
    ix = combinedResponse->size[0];
    stages = combinedResponse->size[1];
    stages *= ix;
    for (ix = 0; ix < stages; ix++) {
      combinedResponse_re = combinedResponse->data[ix].re;
      combinedResponse_im = combinedResponse->data[ix].im;
      d2_re = d2->data[ix].re;
      d2_im = d2->data[ix].im;
      combinedResponse->data[ix].re = combinedResponse_re * d2_re -
        combinedResponse_im * d2_im;
      combinedResponse->data[ix].im = combinedResponse_re * d2_im +
        combinedResponse_im * d2_re;
    }

    stages = 2;
    break;

   case 7:
    /*  Hb3,Hb2,Hb1 */
    memset(&d_u[0], 0, 60U * sizeof(double));
    ix = 0;
    stages = 0;
    for (k = 0; k < 15; k++) {
      d_u[stages] = hb1_coeff[ix];
      ix++;
      stages += 4;
    }

    m_freqz_cg(d_u, w, Fs, combinedResponse);
    memset(&e_u[0], 0, 14U * sizeof(double));
    ix = 0;
    stages = 0;
    for (k = 0; k < 7; k++) {
      e_u[stages] = hb2_coeff[ix];
      ix++;
      stages += 2;
    }

    n_freqz_cg(e_u, w, Fs, d2);
    c_us(hb3_coeff_data, hb3_coeff_size, tmp_data, tmp_size);
    k_freqz_cg(tmp_data, tmp_size, w, Fs, d3);
    ix = combinedResponse->size[0] * combinedResponse->size[1];
    combinedResponse->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)combinedResponse, ix, (int)sizeof
                      (creal_T));
    ix = combinedResponse->size[0];
    stages = combinedResponse->size[1];
    stages *= ix;
    for (ix = 0; ix < stages; ix++) {
      combinedResponse_re = combinedResponse->data[ix].re * d2->data[ix].re -
        combinedResponse->data[ix].im * d2->data[ix].im;
      combinedResponse_im = combinedResponse->data[ix].re * d2->data[ix].im +
        combinedResponse->data[ix].im * d2->data[ix].re;
      d2_re = d3->data[ix].re;
      d2_im = d3->data[ix].im;
      combinedResponse->data[ix].re = combinedResponse_re * d2_re -
        combinedResponse_im * d2_im;
      combinedResponse->data[ix].im = combinedResponse_re * d2_im +
        combinedResponse_im * d2_re;
    }

    stages = 3;
    break;

   case 8:
    /*  Dec/Int3 */
    c_us(dec_int3_coeff_data, dec_int3_coeff_size, tmp_data, tmp_size);
    k_freqz_cg(tmp_data, tmp_size, w, Fs, combinedResponse);
    stages = 1;
    break;

   case 9:
    /*  Dec/Int3,Hb1 */
    memset(&c_u[0], 0, 30U * sizeof(double));
    ix = 0;
    stages = 0;
    for (k = 0; k < 15; k++) {
      c_u[stages] = hb1_coeff[ix];
      ix++;
      stages += 2;
    }

    l_freqz_cg(c_u, w, Fs, combinedResponse);
    c_us(dec_int3_coeff_data, dec_int3_coeff_size, tmp_data, tmp_size);
    k_freqz_cg(tmp_data, tmp_size, w, Fs, d2);
    ix = combinedResponse->size[0] * combinedResponse->size[1];
    combinedResponse->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)combinedResponse, ix, (int)sizeof
                      (creal_T));
    ix = combinedResponse->size[0];
    stages = combinedResponse->size[1];
    stages *= ix;
    for (ix = 0; ix < stages; ix++) {
      combinedResponse_re = combinedResponse->data[ix].re;
      combinedResponse_im = combinedResponse->data[ix].im;
      d2_re = d2->data[ix].re;
      d2_im = d2->data[ix].im;
      combinedResponse->data[ix].re = combinedResponse_re * d2_re -
        combinedResponse_im * d2_im;
      combinedResponse->data[ix].im = combinedResponse_re * d2_im +
        combinedResponse_im * d2_re;
    }

    stages = 2;
    break;

   case 10:
    /*  Dec/Int3,Hb2 */
    memset(&e_u[0], 0, 14U * sizeof(double));
    ix = 0;
    stages = 0;
    for (k = 0; k < 7; k++) {
      e_u[stages] = hb2_coeff[ix];
      ix++;
      stages += 2;
    }

    n_freqz_cg(e_u, w, Fs, combinedResponse);
    c_us(dec_int3_coeff_data, dec_int3_coeff_size, tmp_data, tmp_size);
    k_freqz_cg(tmp_data, tmp_size, w, Fs, d2);
    ix = combinedResponse->size[0] * combinedResponse->size[1];
    combinedResponse->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)combinedResponse, ix, (int)sizeof
                      (creal_T));
    ix = combinedResponse->size[0];
    stages = combinedResponse->size[1];
    stages *= ix;
    for (ix = 0; ix < stages; ix++) {
      combinedResponse_re = combinedResponse->data[ix].re;
      combinedResponse_im = combinedResponse->data[ix].im;
      d2_re = d2->data[ix].re;
      d2_im = d2->data[ix].im;
      combinedResponse->data[ix].re = combinedResponse_re * d2_re -
        combinedResponse_im * d2_im;
      combinedResponse->data[ix].im = combinedResponse_re * d2_im +
        combinedResponse_im * d2_re;
    }

    stages = 2;
    break;

   case 11:
    /*  Dec/Int3,Hb2,Hb1 */
    memset(&d_u[0], 0, 60U * sizeof(double));
    ix = 0;
    stages = 0;
    for (k = 0; k < 15; k++) {
      d_u[stages] = hb1_coeff[ix];
      ix++;
      stages += 4;
    }

    m_freqz_cg(d_u, w, Fs, combinedResponse);
    memset(&e_u[0], 0, 14U * sizeof(double));
    ix = 0;
    stages = 0;
    for (k = 0; k < 7; k++) {
      e_u[stages] = hb2_coeff[ix];
      ix++;
      stages += 2;
    }

    n_freqz_cg(e_u, w, Fs, d2);
    c_us(dec_int3_coeff_data, dec_int3_coeff_size, tmp_data, tmp_size);
    k_freqz_cg(tmp_data, tmp_size, w, Fs, d3);
    ix = combinedResponse->size[0] * combinedResponse->size[1];
    combinedResponse->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)combinedResponse, ix, (int)sizeof
                      (creal_T));
    ix = combinedResponse->size[0];
    stages = combinedResponse->size[1];
    stages *= ix;
    for (ix = 0; ix < stages; ix++) {
      combinedResponse_re = combinedResponse->data[ix].re * d2->data[ix].re -
        combinedResponse->data[ix].im * d2->data[ix].im;
      combinedResponse_im = combinedResponse->data[ix].re * d2->data[ix].im +
        combinedResponse->data[ix].im * d2->data[ix].re;
      d2_re = d3->data[ix].re;
      d2_im = d3->data[ix].im;
      combinedResponse->data[ix].re = combinedResponse_re * d2_re -
        combinedResponse_im * d2_im;
      combinedResponse->data[ix].im = combinedResponse_re * d2_im +
        combinedResponse_im * d2_re;
    }

    stages = 3;
    break;
  }

  emxFree_creal_T(&d3);

  /*  Add filter extra filter to end of cascade */
  if (!(extraTaps_size[1] == 0)) {
    c = (int)rt_powd_snf(2.0, stages);
    ix = c * extraTaps_size[1];
    u_size[0] = 1;
    u_size[1] = (short)ix;
    stages = (short)ix;
    for (ix = 0; ix < stages; ix++) {
      u_data[ix] = 0.0;
    }

    ix = 0;
    stages = 0;
    for (k = 1; k <= extraTaps_size[1]; k++) {
      u_data[stages] = extraTaps_data[ix];
      ix++;
      stages += c;
    }

    k_freqz_cg(u_data, u_size, w, Fs, d2);
    ix = combinedResponse->size[0] * combinedResponse->size[1];
    combinedResponse->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)combinedResponse, ix, (int)sizeof
                      (creal_T));
    ix = combinedResponse->size[0];
    stages = combinedResponse->size[1];
    stages *= ix;
    for (ix = 0; ix < stages; ix++) {
      combinedResponse_re = combinedResponse->data[ix].re;
      combinedResponse_im = combinedResponse->data[ix].im;
      d2_re = d2->data[ix].re;
      d2_im = d2->data[ix].im;
      combinedResponse->data[ix].re = combinedResponse_re * d2_re -
        combinedResponse_im * d2_im;
      combinedResponse->data[ix].im = combinedResponse_re * d2_im +
        combinedResponse_im * d2_re;
    }
  }

  emxFree_creal_T(&d2);
}

/*
 * Arguments    : const double p_data[]
 *                const int p_size[2]
 *                const creal_T x[2048]
 *                creal_T y[2048]
 * Return Type  : void
 */
static void c_polyval(const double p_data[], const int p_size[2], const creal_T
                      x[2048], creal_T y[2048])
{
  int i12;
  int k;
  double x_im;
  if (!(p_size[1] == 0)) {
    for (i12 = 0; i12 < 2048; i12++) {
      y[i12].re = p_data[0];
      y[i12].im = 0.0;
    }

    for (k = 0; k <= p_size[1] - 2; k++) {
      for (i12 = 0; i12 < 2048; i12++) {
        x_im = x[i12].re * y[i12].im + x[i12].im * y[i12].re;
        y[i12].re = (x[i12].re * y[i12].re - x[i12].im * y[i12].im) + p_data[k +
          1];
        y[i12].im = x_im;
      }
    }
  }
}

/*
 * Arguments    : const emxArray_real_T *x
 *                const emxArray_real_T *y
 *                emxArray_real_T *z
 * Return Type  : void
 */
static void c_rdivide(const emxArray_real_T *x, const emxArray_real_T *y,
                      emxArray_real_T *z)
{
  int i34;
  int loop_ub;
  i34 = z->size[0] * z->size[1];
  z->size[0] = 1;
  z->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)z, i34, (int)sizeof(double));
  loop_ub = x->size[0] * x->size[1];
  for (i34 = 0; i34 < loop_ub; i34++) {
    z->data[i34] = x->data[i34] / y->data[i34];
  }
}

/*
 * Arguments    : const char a[2]
 * Return Type  : boolean_T
 */
static boolean_T c_strcmp(const char a[2])
{
  boolean_T b_bool;
  int kstr;
  int exitg1;
  static const char cv26[2] = { 'T', 'x' };

  b_bool = false;
  kstr = 0;
  do {
    exitg1 = 0;
    if (kstr + 1 < 3) {
      if (a[kstr] != cv26[kstr]) {
        exitg1 = 1;
      } else {
        kstr++;
      }
    } else {
      b_bool = true;
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  return b_bool;
}

/*
 * Arguments    : const double o_data[]
 *                const int o_size[2]
 *                double u_data[]
 *                int u_size[2]
 * Return Type  : void
 */
static void c_us(const double o_data[], const int o_size[2], double u_data[],
                 int u_size[2])
{
  int ix;
  int iy;
  int k;
  u_size[0] = 1;
  u_size[1] = (signed char)o_size[1];
  ix = (signed char)o_size[1];
  for (iy = 0; iy < ix; iy++) {
    u_data[iy] = 0.0;
  }

  ix = 0;
  iy = 0;
  for (k = 1; k <= o_size[1]; k++) {
    u_data[iy] = o_data[ix];
    ix++;
    iy++;
  }
}

/*
 * Arguments    : const char * varargin_1
 * Return Type  : int
 */
static int cfprintf(const char * varargin_1)
{
  int nbytesint;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T autoflush;
  static const char cfmt[4] = { '%', 's', '\x0a', '\x00' };

  b_NULL = NULL;
  nbytesint = 0;
  fileManager(&filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    nbytesint = fprintf(filestar, cfmt, varargin_1);
    fflush(filestar);
  }

  return nbytesint;
}

/*
 * Arguments    : double dBinput
 * Return Type  : double
 */
static double dBinv(double dBinput)
{
  double dBoutput;
  if (dBinput > -150.0) {
    dBoutput = rt_powd_snf(10.0, dBinput / 20.0);
  } else {
    dBoutput = 0.0;
  }

  return dBoutput;
}

/*
 * UNTITLED Summary of this function goes here
 *    Detailed explanation goes here
 * Arguments    : const emxArray_real_T *f
 *                double Fconverter
 *                const double b1_data[]
 *                const int b1_size[2]
 *                const emxArray_creal_T *a1
 *                const double b2_data[]
 *                const int b2_size[2]
 *                const emxArray_creal_T *a2
 *                emxArray_creal_T *abc
 * Return Type  : void
 */
static void d_analogresp(const emxArray_real_T *f, double Fconverter, const
  double b1_data[], const int b1_size[2], const emxArray_creal_T *a1, const
  double b2_data[], const int b2_size[2], const emxArray_creal_T *a2,
  emxArray_creal_T *abc)
{
  emxArray_real_T *r32;
  emxArray_real_T *r33;
  emxArray_real_T *r34;
  int i49;
  int loop_ub;
  emxArray_real_T *r35;
  emxArray_creal_T *r36;
  double re;
  double im;
  double b_re;
  double b_im;
  emxInit_real_T(&r32, 2);
  emxInit_real_T(&r33, 2);
  emxInit_real_T(&r34, 2);
  rdivide(f, Fconverter, r33);
  b_sinc(r33, r32);
  i49 = r34->size[0] * r34->size[1];
  r34->size[0] = 1;
  r34->size[1] = f->size[1];
  emxEnsureCapacity((emxArray__common *)r34, i49, (int)sizeof(double));
  loop_ub = f->size[0] * f->size[1];
  emxFree_real_T(&r33);
  for (i49 = 0; i49 < loop_ub; i49++) {
    r34->data[i49] = 6.2831853071795862 * f->data[i49];
  }

  emxInit_real_T(&r35, 2);
  b_freqs_cg(b1_data, b1_size, a1, r34, abc);
  i49 = r35->size[0] * r35->size[1];
  r35->size[0] = 1;
  r35->size[1] = f->size[1];
  emxEnsureCapacity((emxArray__common *)r35, i49, (int)sizeof(double));
  loop_ub = f->size[0] * f->size[1];
  emxFree_real_T(&r34);
  for (i49 = 0; i49 < loop_ub; i49++) {
    r35->data[i49] = 6.2831853071795862 * f->data[i49];
  }

  emxInit_creal_T(&r36, 2);
  b_freqs_cg(b2_data, b2_size, a2, r35, r36);
  i49 = abc->size[0] * abc->size[1];
  abc->size[0] = 1;
  abc->size[1] = r32->size[1];
  emxEnsureCapacity((emxArray__common *)abc, i49, (int)sizeof(creal_T));
  loop_ub = r32->size[0] * r32->size[1];
  emxFree_real_T(&r35);
  for (i49 = 0; i49 < loop_ub; i49++) {
    re = r32->data[i49] * abc->data[i49].re;
    im = r32->data[i49] * abc->data[i49].im;
    b_re = r36->data[i49].re;
    b_im = r36->data[i49].im;
    abc->data[i49].re = re * b_re - im * b_im;
    abc->data[i49].im = re * b_im + im * b_re;
  }

  emxFree_creal_T(&r36);
  emxFree_real_T(&r32);
}

/*
 * Make b a row
 * Arguments    : const double b_data[]
 *                const int b_size[2]
 *                const struct_T *options
 *                creal_T h[2048]
 *                double w[2048]
 * Return Type  : void
 */
static void d_firfreqz(const double b_data[], const int b_size[2], const
  struct_T *options, creal_T h[2048], double w[2048])
{
  int b_b_size[2];
  int loop_ub;
  int i57;
  double b_b_data[29];
  double digw[2048];
  creal_T dcv4[2048];
  double bim;
  double h_re;
  double brm;
  double d;

  /* -------------------------------------------------------------------------- */
  b_b_size[0] = 1;
  b_b_size[1] = b_size[1];
  loop_ub = b_size[1];
  for (i57 = 0; i57 < loop_ub; i57++) {
    b_b_data[i57] = b_data[i57];
  }

  /*  Actual Frequency Response Computation */
  /* if fvflag, */
  /*    Frequency vector specified.  Use Horner's method of polynomial */
  /*    evaluation at the frequency points and divide the numerator */
  /*    by the denominator. */
  /*  */
  /*    Note: we use positive i here because of the relationship */
  /*             polyval(a,exp(1i*w)) = fft(a).*exp(1i*w*(length(a)-1)) */
  /*                ( assuming w = 2*pi*(0:length(a)-1)/length(a) ) */
  /*  */
  /*  Fs was specified, freq. vector is in Hz */
  /*  Convert from Hz to rad/sample for computational purposes */
  /*  Digital frequency must be used for this calculation */
  for (i57 = 0; i57 < 2048; i57++) {
    w[i57] = options->w[i57];
    bim = 6.2831853071795862 * options->w[i57] / options->Fs;
    dcv4[i57].re = bim * 0.0;
    dcv4[i57].im = bim;
    digw[i57] = bim;
  }

  b_exp(dcv4);
  c_polyval(b_b_data, b_b_size, dcv4, h);
  for (i57 = 0; i57 < 2048; i57++) {
    dcv4[i57].re = ((double)b_b_size[1] - 1.0) * (digw[i57] * 0.0);
    dcv4[i57].im = ((double)b_b_size[1] - 1.0) * digw[i57];
  }

  b_exp(dcv4);
  for (i57 = 0; i57 < 2048; i57++) {
    h_re = h[i57].re;
    if (dcv4[i57].im == 0.0) {
      if (h[i57].im == 0.0) {
        h[i57].re /= dcv4[i57].re;
        h[i57].im = 0.0;
      } else if (h[i57].re == 0.0) {
        h[i57].re = 0.0;
        h[i57].im /= dcv4[i57].re;
      } else {
        h[i57].re /= dcv4[i57].re;
        h[i57].im /= dcv4[i57].re;
      }
    } else if (dcv4[i57].re == 0.0) {
      if (h[i57].re == 0.0) {
        h[i57].re = h[i57].im / dcv4[i57].im;
        h[i57].im = 0.0;
      } else if (h[i57].im == 0.0) {
        h[i57].re = 0.0;
        h[i57].im = -(h_re / dcv4[i57].im);
      } else {
        h[i57].re = h[i57].im / dcv4[i57].im;
        h[i57].im = -(h_re / dcv4[i57].im);
      }
    } else {
      brm = fabs(dcv4[i57].re);
      bim = fabs(dcv4[i57].im);
      if (brm > bim) {
        bim = dcv4[i57].im / dcv4[i57].re;
        d = dcv4[i57].re + bim * dcv4[i57].im;
        h[i57].re = (h[i57].re + bim * h[i57].im) / d;
        h[i57].im = (h[i57].im - bim * h_re) / d;
      } else if (bim == brm) {
        if (dcv4[i57].re > 0.0) {
          bim = 0.5;
        } else {
          bim = -0.5;
        }

        if (dcv4[i57].im > 0.0) {
          d = 0.5;
        } else {
          d = -0.5;
        }

        h[i57].re = (h[i57].re * bim + h[i57].im * d) / brm;
        h[i57].im = (h[i57].im * bim - h_re * d) / brm;
      } else {
        bim = dcv4[i57].re / dcv4[i57].im;
        d = dcv4[i57].im + bim * dcv4[i57].re;
        h[i57].re = (bim * h[i57].re + h[i57].im) / d;
        h[i57].im = (bim * h[i57].im - h_re) / d;
      }
    }
  }
}

/*
 * FREQZ_CG Frequency response of digital filter with codegen support
 *
 *  This function is based on 'freqz' by The MathWorks Inc.
 * Arguments    : const double b_data[]
 *                const int b_size[2]
 *                const double w[2048]
 *                double Fs
 *                creal_T hh[2048]
 * Return Type  : void
 */
static void d_freqz_cg(const double b_data[], const int b_size[2], const double
  w[2048], double Fs, creal_T hh[2048])
{
  static struct_T options;
  int i13;
  static const char cv18[8] = { 'o', 'n', 'e', 's', 'i', 'd', 'e', 'd' };

  double b_w[2048];
  static const char cv19[7] = { 'o', 'm', 'i', 't', 't', 'e', 'd' };

  /*  Cast to enforce precision rules */
  options.Fs = Fs;
  memcpy(&options.nfft[0], &w[0], sizeof(double) << 11);
  memcpy(&options.w[0], &w[0], sizeof(double) << 11);

  /*  Remaining are default or for advanced use */
  options.fvflag = 1.0;
  for (i13 = 0; i13 < 8; i13++) {
    options.range[i13] = cv18[i13];
  }

  options.centerdc = 0.0;
  for (i13 = 0; i13 < 7; i13++) {
    options.configlevel[i13] = cv19[i13];
  }

  d_firfreqz(b_data, b_size, &options, hh, b_w);

  /*  Generate the default structure to pass to freqzplot */
  /*  If rad/sample, Fs is empty */
}

/*
 * Arguments    : const double p[30]
 *                const creal_T x[2048]
 *                creal_T y[2048]
 * Return Type  : void
 */
static void d_polyval(const double p[30], const creal_T x[2048], creal_T y[2048])
{
  int i14;
  int k;
  double x_im;
  for (i14 = 0; i14 < 2048; i14++) {
    y[i14].re = p[0];
    y[i14].im = 0.0;
  }

  for (k = 0; k < 29; k++) {
    for (i14 = 0; i14 < 2048; i14++) {
      x_im = x[i14].re * y[i14].im + x[i14].im * y[i14].re;
      y[i14].re = (x[i14].re * y[i14].re - x[i14].im * y[i14].im) + p[k + 1];
      y[i14].im = x_im;
    }
  }
}

/*
 * Arguments    : const double o[15]
 *                double u[30]
 * Return Type  : void
 */
static void d_us(const double o[15], double u[30])
{
  int ix;
  int iy;
  int k;
  memset(&u[0], 0, 30U * sizeof(double));
  ix = 0;
  iy = 0;
  for (k = 0; k < 15; k++) {
    u[iy] = o[ix];
    ix++;
    iy += 2;
  }
}

/*
 * Arguments    : int numerator
 *                int denominator
 * Return Type  : int
 */
static int div_s32_floor(int numerator, int denominator)
{
  int quotient;
  unsigned int absNumerator;
  unsigned int absDenominator;
  boolean_T quotientNeedsNegation;
  unsigned int tempAbsQuotient;
  if (denominator == 0) {
    if (numerator >= 0) {
      quotient = MAX_int32_T;
    } else {
      quotient = MIN_int32_T;
    }
  } else {
    if (numerator < 0) {
      absNumerator = ~(unsigned int)numerator + 1U;
    } else {
      absNumerator = (unsigned int)numerator;
    }

    if (denominator < 0) {
      absDenominator = ~(unsigned int)denominator + 1U;
    } else {
      absDenominator = (unsigned int)denominator;
    }

    quotientNeedsNegation = ((numerator < 0) != (denominator < 0));
    tempAbsQuotient = absNumerator / absDenominator;
    if (quotientNeedsNegation) {
      absNumerator %= absDenominator;
      if (absNumerator > 0U) {
        tempAbsQuotient++;
      }

      quotient = -(int)tempAbsQuotient;
    } else {
      quotient = (int)tempAbsQuotient;
    }
  }

  return quotient;
}

/*
 * Make b a row
 * Arguments    : const double b[30]
 *                const struct_T *options
 *                creal_T h[2048]
 *                double w[2048]
 * Return Type  : void
 */
static void e_firfreqz(const double b[30], const struct_T *options, creal_T h
  [2048], double w[2048])
{
  double digw[2048];
  creal_T dcv5[2048];
  int i58;
  double bim;
  double h_re;
  double brm;
  double d;

  /* -------------------------------------------------------------------------- */
  /*  Actual Frequency Response Computation */
  /* if fvflag, */
  /*    Frequency vector specified.  Use Horner's method of polynomial */
  /*    evaluation at the frequency points and divide the numerator */
  /*    by the denominator. */
  /*  */
  /*    Note: we use positive i here because of the relationship */
  /*             polyval(a,exp(1i*w)) = fft(a).*exp(1i*w*(length(a)-1)) */
  /*                ( assuming w = 2*pi*(0:length(a)-1)/length(a) ) */
  /*  */
  /*  Fs was specified, freq. vector is in Hz */
  /*  Convert from Hz to rad/sample for computational purposes */
  /*  Digital frequency must be used for this calculation */
  for (i58 = 0; i58 < 2048; i58++) {
    w[i58] = options->w[i58];
    bim = 6.2831853071795862 * options->w[i58] / options->Fs;
    dcv5[i58].re = bim * 0.0;
    dcv5[i58].im = bim;
    digw[i58] = bim;
  }

  b_exp(dcv5);
  d_polyval(b, dcv5, h);
  for (i58 = 0; i58 < 2048; i58++) {
    dcv5[i58].re = 29.0 * (digw[i58] * 0.0);
    dcv5[i58].im = 29.0 * digw[i58];
  }

  b_exp(dcv5);
  for (i58 = 0; i58 < 2048; i58++) {
    h_re = h[i58].re;
    if (dcv5[i58].im == 0.0) {
      if (h[i58].im == 0.0) {
        h[i58].re /= dcv5[i58].re;
        h[i58].im = 0.0;
      } else if (h[i58].re == 0.0) {
        h[i58].re = 0.0;
        h[i58].im /= dcv5[i58].re;
      } else {
        h[i58].re /= dcv5[i58].re;
        h[i58].im /= dcv5[i58].re;
      }
    } else if (dcv5[i58].re == 0.0) {
      if (h[i58].re == 0.0) {
        h[i58].re = h[i58].im / dcv5[i58].im;
        h[i58].im = 0.0;
      } else if (h[i58].im == 0.0) {
        h[i58].re = 0.0;
        h[i58].im = -(h_re / dcv5[i58].im);
      } else {
        h[i58].re = h[i58].im / dcv5[i58].im;
        h[i58].im = -(h_re / dcv5[i58].im);
      }
    } else {
      brm = fabs(dcv5[i58].re);
      bim = fabs(dcv5[i58].im);
      if (brm > bim) {
        bim = dcv5[i58].im / dcv5[i58].re;
        d = dcv5[i58].re + bim * dcv5[i58].im;
        h[i58].re = (h[i58].re + bim * h[i58].im) / d;
        h[i58].im = (h[i58].im - bim * h_re) / d;
      } else if (bim == brm) {
        if (dcv5[i58].re > 0.0) {
          bim = 0.5;
        } else {
          bim = -0.5;
        }

        if (dcv5[i58].im > 0.0) {
          d = 0.5;
        } else {
          d = -0.5;
        }

        h[i58].re = (h[i58].re * bim + h[i58].im * d) / brm;
        h[i58].im = (h[i58].im * bim - h_re * d) / brm;
      } else {
        bim = dcv5[i58].re / dcv5[i58].im;
        d = dcv5[i58].im + bim * dcv5[i58].re;
        h[i58].re = (bim * h[i58].re + h[i58].im) / d;
        h[i58].im = (bim * h[i58].im - h_re) / d;
      }
    }
  }
}

/*
 * FREQZ_CG Frequency response of digital filter with codegen support
 *
 *  This function is based on 'freqz' by The MathWorks Inc.
 * Arguments    : const double b[30]
 *                const double w[2048]
 *                double Fs
 *                creal_T hh[2048]
 * Return Type  : void
 */
static void e_freqz_cg(const double b[30], const double w[2048], double Fs,
  creal_T hh[2048])
{
  static struct_T options;
  int i15;
  static const char cv20[8] = { 'o', 'n', 'e', 's', 'i', 'd', 'e', 'd' };

  double b_w[2048];
  static const char cv21[7] = { 'o', 'm', 'i', 't', 't', 'e', 'd' };

  /*  Cast to enforce precision rules */
  options.Fs = Fs;
  memcpy(&options.nfft[0], &w[0], sizeof(double) << 11);
  memcpy(&options.w[0], &w[0], sizeof(double) << 11);

  /*  Remaining are default or for advanced use */
  options.fvflag = 1.0;
  for (i15 = 0; i15 < 8; i15++) {
    options.range[i15] = cv20[i15];
  }

  options.centerdc = 0.0;
  for (i15 = 0; i15 < 7; i15++) {
    options.configlevel[i15] = cv21[i15];
  }

  e_firfreqz(b, &options, hh, b_w);

  /*  Generate the default structure to pass to freqzplot */
  /*  If rad/sample, Fs is empty */
}

/*
 * Arguments    : const double p[60]
 *                const creal_T x[2048]
 *                creal_T y[2048]
 * Return Type  : void
 */
static void e_polyval(const double p[60], const creal_T x[2048], creal_T y[2048])
{
  int i16;
  int k;
  double x_im;
  for (i16 = 0; i16 < 2048; i16++) {
    y[i16].re = p[0];
    y[i16].im = 0.0;
  }

  for (k = 0; k < 59; k++) {
    for (i16 = 0; i16 < 2048; i16++) {
      x_im = x[i16].re * y[i16].im + x[i16].im * y[i16].re;
      y[i16].re = (x[i16].re * y[i16].re - x[i16].im * y[i16].im) + p[k + 1];
      y[i16].im = x_im;
    }
  }
}

/*
 * Arguments    : const double o[15]
 *                double u[60]
 * Return Type  : void
 */
static void e_us(const double o[15], double u[60])
{
  int ix;
  int iy;
  int k;
  memset(&u[0], 0, 60U * sizeof(double));
  ix = 0;
  iy = 0;
  for (k = 0; k < 15; k++) {
    u[iy] = o[ix];
    ix++;
    iy += 4;
  }
}

/*
 * Arguments    : emxArray_creal_T *h
 * Return Type  : int
 */
static int eml_zlahqr(emxArray_creal_T *h)
{
  int info;
  int n;
  int ldh;
  creal_T v[2];
  int knt;
  int i;
  double SMLNUM;
  double htmp1;
  boolean_T exitg1;
  double tst;
  double aa;
  int L;
  creal_T u2;
  boolean_T goto140;
  int its;
  boolean_T exitg2;
  int k;
  boolean_T exitg4;
  boolean_T guard1 = false;
  creal_T b_u2;
  double t_re;
  double t_im;
  int b_n;
  creal_T x2;
  boolean_T goto70;
  double ab;
  int m;
  double ba;
  boolean_T exitg3;
  double u_re;
  double u_im;
  double s;
  int b_k;
  double b_SMLNUM;
  creal_T c_u2;
  creal_T b_v;
  int c_k;
  creal_T d_u2;
  creal_T e_u2;
  n = h->size[0];
  ldh = h->size[0];
  info = 0;
  if ((h->size[0] != 0) && (1 != h->size[0])) {
    for (knt = 0; knt + 1 <= n - 3; knt++) {
      h->data[(knt + h->size[0] * knt) + 2].re = 0.0;
      h->data[(knt + h->size[0] * knt) + 2].im = 0.0;
      h->data[(knt + h->size[0] * knt) + 3].re = 0.0;
      h->data[(knt + h->size[0] * knt) + 3].im = 0.0;
    }

    if (1 <= n - 2) {
      h->data[(n + h->size[0] * (n - 3)) - 1].re = 0.0;
      h->data[(n + h->size[0] * (n - 3)) - 1].im = 0.0;
    }

    for (i = 1; i + 1 <= n; i++) {
      if (h->data[i + h->size[0] * (i - 1)].im != 0.0) {
        htmp1 = h->data[i + h->size[0] * (i - 1)].re;
        tst = h->data[i + h->size[0] * (i - 1)].im;
        aa = fabs(h->data[i + h->size[0] * (i - 1)].re) + fabs(h->data[i +
          h->size[0] * (i - 1)].im);
        if (tst == 0.0) {
          u2.re = htmp1 / aa;
          u2.im = 0.0;
        } else if (htmp1 == 0.0) {
          u2.re = 0.0;
          u2.im = tst / aa;
        } else {
          u2.re = htmp1 / aa;
          u2.im = tst / aa;
        }

        aa = rt_hypotd_snf(u2.re, u2.im);
        if (-u2.im == 0.0) {
          u2.re /= aa;
          u2.im = 0.0;
        } else if (u2.re == 0.0) {
          u2.re = 0.0;
          u2.im = -u2.im / aa;
        } else {
          u2.re /= aa;
          u2.im = -u2.im / aa;
        }

        tst = h->data[i + h->size[0] * (i - 1)].re;
        aa = h->data[i + h->size[0] * (i - 1)].im;
        h->data[i + h->size[0] * (i - 1)].re = rt_hypotd_snf(tst, aa);
        h->data[i + h->size[0] * (i - 1)].im = 0.0;
        b_xscal(n - i, u2, h, (i + i * ldh) + 1, ldh);
        b_u2.re = u2.re;
        b_u2.im = -u2.im;
        if (n <= i + 2) {
          b_n = n;
        } else {
          b_n = i + 2;
        }

        xscal(b_n, b_u2, h, 1 + i * ldh);
      }
    }

    SMLNUM = 2.2250738585072014E-308 * ((double)n / 2.2204460492503131E-16);
    i = n - 1;
    exitg1 = false;
    while ((!exitg1) && (i + 1 >= 1)) {
      L = -1;
      goto140 = false;
      its = 0;
      exitg2 = false;
      while ((!exitg2) && (its < 31)) {
        k = i;
        exitg4 = false;
        while ((!exitg4) && ((k + 1 > L + 2) && (!(fabs(h->data[k + h->size[0] *
                   (k - 1)].re) + fabs(h->data[k + h->size[0] * (k - 1)].im) <=
                  SMLNUM)))) {
          tst = (fabs(h->data[(k + h->size[0] * (k - 1)) - 1].re) + fabs(h->
                  data[(k + h->size[0] * (k - 1)) - 1].im)) + (fabs(h->data[k +
            h->size[0] * k].re) + fabs(h->data[k + h->size[0] * k].im));
          if (tst == 0.0) {
            if (k - 1 >= 1) {
              tst = fabs(h->data[(k + h->size[0] * (k - 2)) - 1].re);
            }

            if (k + 2 <= n) {
              tst += fabs(h->data[(k + h->size[0] * k) + 1].re);
            }
          }

          guard1 = false;
          if (fabs(h->data[k + h->size[0] * (k - 1)].re) <=
              2.2204460492503131E-16 * tst) {
            htmp1 = fabs(h->data[k + h->size[0] * (k - 1)].re) + fabs(h->data[k
              + h->size[0] * (k - 1)].im);
            tst = fabs(h->data[(k + h->size[0] * k) - 1].re) + fabs(h->data[(k +
              h->size[0] * k) - 1].im);
            if (htmp1 > tst) {
              ab = htmp1;
              ba = tst;
            } else {
              ab = tst;
              ba = htmp1;
            }

            htmp1 = fabs(h->data[k + h->size[0] * k].re) + fabs(h->data[k +
              h->size[0] * k].im);
            t_re = h->data[(k + h->size[0] * (k - 1)) - 1].re - h->data[k +
              h->size[0] * k].re;
            t_im = h->data[(k + h->size[0] * (k - 1)) - 1].im - h->data[k +
              h->size[0] * k].im;
            tst = fabs(t_re) + fabs(t_im);
            if (htmp1 > tst) {
              aa = htmp1;
              htmp1 = tst;
            } else {
              aa = tst;
            }

            s = aa + ab;
            tst = 2.2204460492503131E-16 * (htmp1 * (aa / s));
            if ((SMLNUM >= tst) || rtIsNaN(tst)) {
              b_SMLNUM = SMLNUM;
            } else {
              b_SMLNUM = tst;
            }

            if (ba * (ab / s) <= b_SMLNUM) {
              exitg4 = true;
            } else {
              guard1 = true;
            }
          } else {
            guard1 = true;
          }

          if (guard1) {
            k--;
          }
        }

        L = k - 1;
        if (k + 1 > 1) {
          h->data[k + h->size[0] * (k - 1)].re = 0.0;
          h->data[k + h->size[0] * (k - 1)].im = 0.0;
        }

        if (k + 1 >= i + 1) {
          goto140 = true;
          exitg2 = true;
        } else {
          if (its == 10) {
            t_re = 0.75 * fabs(h->data[(k + h->size[0] * k) + 1].re) + h->data[k
              + h->size[0] * k].re;
            t_im = h->data[k + h->size[0] * k].im;
          } else if (its == 20) {
            t_re = 0.75 * fabs(h->data[i + h->size[0] * (i - 1)].re) + h->data[i
              + h->size[0] * i].re;
            t_im = h->data[i + h->size[0] * i].im;
          } else {
            t_re = h->data[i + h->size[0] * i].re;
            t_im = h->data[i + h->size[0] * i].im;
            x2 = h->data[(i + h->size[0] * i) - 1];
            b_sqrt(&x2);
            u2 = h->data[i + h->size[0] * (i - 1)];
            b_sqrt(&u2);
            u_re = x2.re * u2.re - x2.im * u2.im;
            u_im = x2.re * u2.im + x2.im * u2.re;
            s = fabs(u_re) + fabs(u_im);
            if (s != 0.0) {
              tst = h->data[(i + h->size[0] * (i - 1)) - 1].re - h->data[i +
                h->size[0] * i].re;
              aa = h->data[(i + h->size[0] * (i - 1)) - 1].im - h->data[i +
                h->size[0] * i].im;
              t_re = 0.5 * tst;
              t_im = 0.5 * aa;
              aa = fabs(t_re) + fabs(t_im);
              tst = fabs(t_re) + fabs(t_im);
              if (!((s >= tst) || rtIsNaN(tst))) {
                s = tst;
              }

              if (t_im == 0.0) {
                x2.re = t_re / s;
                x2.im = 0.0;
              } else if (t_re == 0.0) {
                x2.re = 0.0;
                x2.im = t_im / s;
              } else {
                x2.re = t_re / s;
                x2.im = t_im / s;
              }

              htmp1 = x2.re;
              ab = x2.re;
              x2.re = x2.re * x2.re - x2.im * x2.im;
              x2.im = htmp1 * x2.im + x2.im * ab;
              if (u_im == 0.0) {
                u2.re = u_re / s;
                u2.im = 0.0;
              } else if (u_re == 0.0) {
                u2.re = 0.0;
                u2.im = u_im / s;
              } else {
                u2.re = u_re / s;
                u2.im = u_im / s;
              }

              x2.re += u2.re * u2.re - u2.im * u2.im;
              x2.im += u2.re * u2.im + u2.im * u2.re;
              b_sqrt(&x2);
              u2.re = s * x2.re;
              u2.im = s * x2.im;
              if (aa > 0.0) {
                if (t_im == 0.0) {
                  x2.re = t_re / aa;
                  x2.im = 0.0;
                } else if (t_re == 0.0) {
                  x2.re = 0.0;
                  x2.im = t_im / aa;
                } else {
                  x2.re = t_re / aa;
                  x2.im = t_im / aa;
                }

                if (x2.re * u2.re + x2.im * u2.im < 0.0) {
                  u2.re = -u2.re;
                  u2.im = -u2.im;
                }
              }

              aa = t_re + u2.re;
              htmp1 = t_im + u2.im;
              if (htmp1 == 0.0) {
                if (u_im == 0.0) {
                  ba = u_re / aa;
                  tst = 0.0;
                } else if (u_re == 0.0) {
                  ba = 0.0;
                  tst = u_im / aa;
                } else {
                  ba = u_re / aa;
                  tst = u_im / aa;
                }
              } else if (aa == 0.0) {
                if (u_re == 0.0) {
                  ba = u_im / htmp1;
                  tst = 0.0;
                } else if (u_im == 0.0) {
                  ba = 0.0;
                  tst = -(u_re / htmp1);
                } else {
                  ba = u_im / htmp1;
                  tst = -(u_re / htmp1);
                }
              } else {
                ab = fabs(aa);
                tst = fabs(htmp1);
                if (ab > tst) {
                  s = htmp1 / aa;
                  tst = aa + s * htmp1;
                  ba = (u_re + s * u_im) / tst;
                  tst = (u_im - s * u_re) / tst;
                } else if (tst == ab) {
                  if (aa > 0.0) {
                    aa = 0.5;
                  } else {
                    aa = -0.5;
                  }

                  if (htmp1 > 0.0) {
                    tst = 0.5;
                  } else {
                    tst = -0.5;
                  }

                  ba = (u_re * aa + u_im * tst) / ab;
                  tst = (u_im * aa - u_re * tst) / ab;
                } else {
                  s = aa / htmp1;
                  tst = htmp1 + s * aa;
                  ba = (s * u_re + u_im) / tst;
                  tst = (s * u_im - u_re) / tst;
                }
              }

              t_re = h->data[i + h->size[0] * i].re - (u_re * ba - u_im * tst);
              t_im = h->data[i + h->size[0] * i].im - (u_re * tst + u_im * ba);
            }
          }

          goto70 = false;
          m = i;
          exitg3 = false;
          while ((!exitg3) && (m > k + 1)) {
            u2.re = h->data[(m + h->size[0] * (m - 1)) - 1].re - t_re;
            u2.im = h->data[(m + h->size[0] * (m - 1)) - 1].im - t_im;
            tst = h->data[m + h->size[0] * (m - 1)].re;
            s = (fabs(u2.re) + fabs(u2.im)) + fabs(tst);
            if (u2.im == 0.0) {
              u2.re /= s;
              u2.im = 0.0;
            } else if (u2.re == 0.0) {
              u2.re = 0.0;
              u2.im /= s;
            } else {
              u2.re /= s;
              u2.im /= s;
            }

            tst /= s;
            v[0] = u2;
            v[1].re = tst;
            v[1].im = 0.0;
            if (fabs(h->data[(m + h->size[0] * (m - 2)) - 1].re) * fabs(tst) <=
                2.2204460492503131E-16 * ((fabs(u2.re) + fabs(u2.im)) * ((fabs
                   (h->data[(m + h->size[0] * (m - 1)) - 1].re) + fabs(h->data
                    [(m + h->size[0] * (m - 1)) - 1].im)) + (fabs(h->data[m +
                    h->size[0] * m].re) + fabs(h->data[m + h->size[0] * m].im)))))
            {
              goto70 = true;
              exitg3 = true;
            } else {
              m--;
            }
          }

          if (!goto70) {
            u2.re = h->data[k + h->size[0] * k].re - t_re;
            u2.im = h->data[k + h->size[0] * k].im - t_im;
            tst = h->data[(k + h->size[0] * k) + 1].re;
            s = (fabs(u2.re) + fabs(u2.im)) + fabs(tst);
            if (u2.im == 0.0) {
              u2.re /= s;
              u2.im = 0.0;
            } else if (u2.re == 0.0) {
              u2.re = 0.0;
              u2.im /= s;
            } else {
              u2.re /= s;
              u2.im /= s;
            }

            tst /= s;
            v[0] = u2;
            v[1].re = tst;
            v[1].im = 0.0;
          }

          for (b_k = m; b_k <= i; b_k++) {
            if (b_k > m) {
              v[0] = h->data[(b_k + h->size[0] * (b_k - 2)) - 1];
              v[1] = h->data[b_k + h->size[0] * (b_k - 2)];
            }

            t_re = v[1].re;
            t_im = v[1].im;
            u2 = v[0];
            x2.re = 0.0;
            x2.im = 0.0;
            tst = rt_hypotd_snf(v[1].re, v[1].im);
            if ((tst != 0.0) || (v[0].im != 0.0)) {
              aa = xdlapy3(v[0].re, v[0].im, tst);
              if (v[0].re >= 0.0) {
                aa = -aa;
              }

              if (fabs(aa) < 1.0020841800044864E-292) {
                knt = 0;
                do {
                  knt++;
                  t_re *= 9.9792015476736E+291;
                  t_im *= 9.9792015476736E+291;
                  aa *= 9.9792015476736E+291;
                  u2.re *= 9.9792015476736E+291;
                  u2.im *= 9.9792015476736E+291;
                } while (!(fabs(aa) >= 1.0020841800044864E-292));

                aa = xdlapy3(u2.re, u2.im, rt_hypotd_snf(t_re, t_im));
                if (u2.re >= 0.0) {
                  aa = -aa;
                }

                htmp1 = aa - u2.re;
                if (0.0 - u2.im == 0.0) {
                  x2.re = htmp1 / aa;
                  x2.im = 0.0;
                } else if (htmp1 == 0.0) {
                  x2.re = 0.0;
                  x2.im = (0.0 - u2.im) / aa;
                } else {
                  x2.re = htmp1 / aa;
                  x2.im = (0.0 - u2.im) / aa;
                }

                d_u2.re = u2.re - aa;
                d_u2.im = u2.im;
                u2 = recip(d_u2);
                tst = t_re;
                t_re = u2.re * t_re - u2.im * t_im;
                t_im = u2.re * t_im + u2.im * tst;
                for (c_k = 1; c_k <= knt; c_k++) {
                  aa *= 1.0020841800044864E-292;
                }

                u2.re = aa;
                u2.im = 0.0;
              } else {
                htmp1 = aa - v[0].re;
                if (0.0 - v[0].im == 0.0) {
                  x2.re = htmp1 / aa;
                  x2.im = 0.0;
                } else if (htmp1 == 0.0) {
                  x2.re = 0.0;
                  x2.im = (0.0 - v[0].im) / aa;
                } else {
                  x2.re = htmp1 / aa;
                  x2.im = (0.0 - v[0].im) / aa;
                }

                b_v.re = v[0].re - aa;
                b_v.im = v[0].im;
                b_u2 = recip(b_v);
                t_re = b_u2.re * v[1].re - b_u2.im * v[1].im;
                t_im = b_u2.re * v[1].im + b_u2.im * v[1].re;
                u2.re = aa;
                u2.im = 0.0;
              }
            }

            v[0] = u2;
            v[1].re = t_re;
            v[1].im = t_im;
            if (b_k > m) {
              h->data[(b_k + h->size[0] * (b_k - 2)) - 1] = u2;
              h->data[b_k + h->size[0] * (b_k - 2)].re = 0.0;
              h->data[b_k + h->size[0] * (b_k - 2)].im = 0.0;
            }

            htmp1 = x2.re * t_re - x2.im * t_im;
            for (knt = b_k - 1; knt + 1 <= n; knt++) {
              ab = x2.re * h->data[(b_k + h->size[0] * knt) - 1].re - -x2.im *
                h->data[(b_k + h->size[0] * knt) - 1].im;
              tst = x2.re * h->data[(b_k + h->size[0] * knt) - 1].im + -x2.im *
                h->data[(b_k + h->size[0] * knt) - 1].re;
              u2.re = ab + htmp1 * h->data[b_k + h->size[0] * knt].re;
              u2.im = tst + htmp1 * h->data[b_k + h->size[0] * knt].im;
              h->data[(b_k + h->size[0] * knt) - 1].re -= u2.re;
              h->data[(b_k + h->size[0] * knt) - 1].im -= u2.im;
              h->data[b_k + h->size[0] * knt].re -= u2.re * t_re - u2.im * t_im;
              h->data[b_k + h->size[0] * knt].im -= u2.re * t_im + u2.im * t_re;
            }

            if (b_k + 2 <= i + 1) {
              c_k = b_k;
            } else {
              c_k = i - 1;
            }

            for (knt = 0; knt + 1 <= c_k + 2; knt++) {
              ab = x2.re * h->data[knt + h->size[0] * (b_k - 1)].re - x2.im *
                h->data[knt + h->size[0] * (b_k - 1)].im;
              tst = x2.re * h->data[knt + h->size[0] * (b_k - 1)].im + x2.im *
                h->data[knt + h->size[0] * (b_k - 1)].re;
              u2.re = ab + htmp1 * h->data[knt + h->size[0] * b_k].re;
              u2.im = tst + htmp1 * h->data[knt + h->size[0] * b_k].im;
              h->data[knt + h->size[0] * (b_k - 1)].re -= u2.re;
              h->data[knt + h->size[0] * (b_k - 1)].im -= u2.im;
              h->data[knt + h->size[0] * b_k].re -= u2.re * t_re - u2.im * -t_im;
              h->data[knt + h->size[0] * b_k].im -= u2.re * -t_im + u2.im * t_re;
            }

            if ((b_k == m) && (m > k + 1)) {
              u2.re = 1.0 - x2.re;
              u2.im = 0.0 - x2.im;
              aa = rt_hypotd_snf(u2.re, u2.im);
              if (u2.im == 0.0) {
                u2.re /= aa;
                u2.im = 0.0;
              } else if (u2.re == 0.0) {
                u2.re = 0.0;
                u2.im /= aa;
              } else {
                u2.re /= aa;
                u2.im /= aa;
              }

              tst = h->data[m + h->size[0] * (m - 1)].re;
              aa = h->data[m + h->size[0] * (m - 1)].im;
              h->data[m + h->size[0] * (m - 1)].re = tst * u2.re - aa * -u2.im;
              h->data[m + h->size[0] * (m - 1)].im = tst * -u2.im + aa * u2.re;
              if (m + 2 <= i + 1) {
                tst = h->data[(m + h->size[0] * m) + 1].re;
                aa = h->data[(m + h->size[0] * m) + 1].im;
                h->data[(m + h->size[0] * m) + 1].re = tst * u2.re - aa * u2.im;
                h->data[(m + h->size[0] * m) + 1].im = tst * u2.im + aa * u2.re;
              }

              for (knt = m; knt <= i + 1; knt++) {
                if (knt != m + 1) {
                  if (n > knt) {
                    b_xscal(n - knt, u2, h, knt + knt * ldh, ldh);
                  }

                  e_u2.re = u2.re;
                  e_u2.im = -u2.im;
                  xscal(knt - 1, e_u2, h, 1 + (knt - 1) * ldh);
                }
              }
            }
          }

          u2 = h->data[i + h->size[0] * (i - 1)];
          if (h->data[i + h->size[0] * (i - 1)].im != 0.0) {
            tst = rt_hypotd_snf(h->data[i + h->size[0] * (i - 1)].re, h->data[i
                                + h->size[0] * (i - 1)].im);
            h->data[i + h->size[0] * (i - 1)].re = tst;
            h->data[i + h->size[0] * (i - 1)].im = 0.0;
            if (u2.im == 0.0) {
              u2.re /= tst;
              u2.im = 0.0;
            } else if (u2.re == 0.0) {
              u2.re = 0.0;
              u2.im /= tst;
            } else {
              u2.re /= tst;
              u2.im /= tst;
            }

            if (n > i + 1) {
              c_u2.re = u2.re;
              c_u2.im = -u2.im;
              b_xscal((n - i) - 1, c_u2, h, (i + (i + 1) * ldh) + 1, ldh);
            }

            xscal(i, u2, h, 1 + i * ldh);
          }

          its++;
        }
      }

      if (!goto140) {
        info = i + 1;
        exitg1 = true;
      } else {
        i = L;
      }
    }
  }

  return info;
}

/*
 * Arguments    : emxArray__common *emxArray
 *                int oldNumel
 *                int elementSize
 * Return Type  : void
 */
static void emxEnsureCapacity(emxArray__common *emxArray, int oldNumel, int
  elementSize)
{
  int newNumel;
  int i;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }

  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }

  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }

    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i <<= 1;
      }
    }

    newData = calloc((unsigned int)i, (unsigned int)elementSize);
    if (emxArray->data != NULL) {
      memcpy(newData, emxArray->data, (unsigned int)(elementSize * oldNumel));
      if (emxArray->canFreeData) {
        free(emxArray->data);
      }
    }

    emxArray->data = newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

/*
 * Arguments    : emxArray_boolean_T **pEmxArray
 * Return Type  : void
 */
static void emxFree_boolean_T(emxArray_boolean_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_boolean_T *)NULL) {
    if (((*pEmxArray)->data != (boolean_T *)NULL) && (*pEmxArray)->canFreeData)
    {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_boolean_T *)NULL;
  }
}

/*
 * Arguments    : emxArray_cint8_T **pEmxArray
 * Return Type  : void
 */
static void emxFree_cint8_T(emxArray_cint8_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_cint8_T *)NULL) {
    if (((*pEmxArray)->data != (cint8_T *)NULL) && (*pEmxArray)->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_cint8_T *)NULL;
  }
}

/*
 * Arguments    : emxArray_creal_T **pEmxArray
 * Return Type  : void
 */
static void emxFree_creal_T(emxArray_creal_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_creal_T *)NULL) {
    if (((*pEmxArray)->data != (creal_T *)NULL) && (*pEmxArray)->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_creal_T *)NULL;
  }
}

/*
 * Arguments    : emxArray_int32_T **pEmxArray
 * Return Type  : void
 */
static void emxFree_int32_T(emxArray_int32_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_int32_T *)NULL) {
    if (((*pEmxArray)->data != (int *)NULL) && (*pEmxArray)->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_int32_T *)NULL;
  }
}

/*
 * Arguments    : emxArray_int8_T **pEmxArray
 * Return Type  : void
 */
static void emxFree_int8_T(emxArray_int8_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_int8_T *)NULL) {
    if (((*pEmxArray)->data != (signed char *)NULL) && (*pEmxArray)->canFreeData)
    {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_int8_T *)NULL;
  }
}

/*
 * Arguments    : emxArray_real_T **pEmxArray
 * Return Type  : void
 */
static void emxFree_real_T(emxArray_real_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_real_T *)NULL) {
    if (((*pEmxArray)->data != (double *)NULL) && (*pEmxArray)->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_real_T *)NULL;
  }
}

/*
 * Arguments    : emxArray_uint32_T **pEmxArray
 * Return Type  : void
 */
static void emxFree_uint32_T(emxArray_uint32_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_uint32_T *)NULL) {
    if (((*pEmxArray)->data != (unsigned int *)NULL) && (*pEmxArray)
        ->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_uint32_T *)NULL;
  }
}

/*
 * Arguments    : emxArray_boolean_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
static void emxInit_boolean_T(emxArray_boolean_T **pEmxArray, int numDimensions)
{
  emxArray_boolean_T *emxArray;
  int i;
  *pEmxArray = (emxArray_boolean_T *)malloc(sizeof(emxArray_boolean_T));
  emxArray = *pEmxArray;
  emxArray->data = (boolean_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 * Arguments    : emxArray_cint8_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
static void emxInit_cint8_T(emxArray_cint8_T **pEmxArray, int numDimensions)
{
  emxArray_cint8_T *emxArray;
  int i;
  *pEmxArray = (emxArray_cint8_T *)malloc(sizeof(emxArray_cint8_T));
  emxArray = *pEmxArray;
  emxArray->data = (cint8_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 * Arguments    : emxArray_creal_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
static void emxInit_creal_T(emxArray_creal_T **pEmxArray, int numDimensions)
{
  emxArray_creal_T *emxArray;
  int i;
  *pEmxArray = (emxArray_creal_T *)malloc(sizeof(emxArray_creal_T));
  emxArray = *pEmxArray;
  emxArray->data = (creal_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 * Arguments    : emxArray_creal_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
static void emxInit_creal_T1(emxArray_creal_T **pEmxArray, int numDimensions)
{
  emxArray_creal_T *emxArray;
  int i;
  *pEmxArray = (emxArray_creal_T *)malloc(sizeof(emxArray_creal_T));
  emxArray = *pEmxArray;
  emxArray->data = (creal_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 * Arguments    : emxArray_int32_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
static void emxInit_int32_T(emxArray_int32_T **pEmxArray, int numDimensions)
{
  emxArray_int32_T *emxArray;
  int i;
  *pEmxArray = (emxArray_int32_T *)malloc(sizeof(emxArray_int32_T));
  emxArray = *pEmxArray;
  emxArray->data = (int *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 * Arguments    : emxArray_int32_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
static void emxInit_int32_T1(emxArray_int32_T **pEmxArray, int numDimensions)
{
  emxArray_int32_T *emxArray;
  int i;
  *pEmxArray = (emxArray_int32_T *)malloc(sizeof(emxArray_int32_T));
  emxArray = *pEmxArray;
  emxArray->data = (int *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 * Arguments    : emxArray_int8_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
static void emxInit_int8_T(emxArray_int8_T **pEmxArray, int numDimensions)
{
  emxArray_int8_T *emxArray;
  int i;
  *pEmxArray = (emxArray_int8_T *)malloc(sizeof(emxArray_int8_T));
  emxArray = *pEmxArray;
  emxArray->data = (signed char *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 * Arguments    : emxArray_real_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
static void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions)
{
  emxArray_real_T *emxArray;
  int i;
  *pEmxArray = (emxArray_real_T *)malloc(sizeof(emxArray_real_T));
  emxArray = *pEmxArray;
  emxArray->data = (double *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 * Arguments    : emxArray_real_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
static void emxInit_real_T1(emxArray_real_T **pEmxArray, int numDimensions)
{
  emxArray_real_T *emxArray;
  int i;
  *pEmxArray = (emxArray_real_T *)malloc(sizeof(emxArray_real_T));
  emxArray = *pEmxArray;
  emxArray->data = (double *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 * Arguments    : emxArray_uint32_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
static void emxInit_uint32_T(emxArray_uint32_T **pEmxArray, int numDimensions)
{
  emxArray_uint32_T *emxArray;
  int i;
  *pEmxArray = (emxArray_uint32_T *)malloc(sizeof(emxArray_uint32_T));
  emxArray = *pEmxArray;
  emxArray->data = (unsigned int *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 * Make b a row
 * Arguments    : const double b[60]
 *                const struct_T *options
 *                creal_T h[2048]
 *                double w[2048]
 * Return Type  : void
 */
static void f_firfreqz(const double b[60], const struct_T *options, creal_T h
  [2048], double w[2048])
{
  double digw[2048];
  creal_T dcv6[2048];
  int i59;
  double bim;
  double h_re;
  double brm;
  double d;

  /* -------------------------------------------------------------------------- */
  /*  Actual Frequency Response Computation */
  /* if fvflag, */
  /*    Frequency vector specified.  Use Horner's method of polynomial */
  /*    evaluation at the frequency points and divide the numerator */
  /*    by the denominator. */
  /*  */
  /*    Note: we use positive i here because of the relationship */
  /*             polyval(a,exp(1i*w)) = fft(a).*exp(1i*w*(length(a)-1)) */
  /*                ( assuming w = 2*pi*(0:length(a)-1)/length(a) ) */
  /*  */
  /*  Fs was specified, freq. vector is in Hz */
  /*  Convert from Hz to rad/sample for computational purposes */
  /*  Digital frequency must be used for this calculation */
  for (i59 = 0; i59 < 2048; i59++) {
    w[i59] = options->w[i59];
    bim = 6.2831853071795862 * options->w[i59] / options->Fs;
    dcv6[i59].re = bim * 0.0;
    dcv6[i59].im = bim;
    digw[i59] = bim;
  }

  b_exp(dcv6);
  e_polyval(b, dcv6, h);
  for (i59 = 0; i59 < 2048; i59++) {
    dcv6[i59].re = 59.0 * (digw[i59] * 0.0);
    dcv6[i59].im = 59.0 * digw[i59];
  }

  b_exp(dcv6);
  for (i59 = 0; i59 < 2048; i59++) {
    h_re = h[i59].re;
    if (dcv6[i59].im == 0.0) {
      if (h[i59].im == 0.0) {
        h[i59].re /= dcv6[i59].re;
        h[i59].im = 0.0;
      } else if (h[i59].re == 0.0) {
        h[i59].re = 0.0;
        h[i59].im /= dcv6[i59].re;
      } else {
        h[i59].re /= dcv6[i59].re;
        h[i59].im /= dcv6[i59].re;
      }
    } else if (dcv6[i59].re == 0.0) {
      if (h[i59].re == 0.0) {
        h[i59].re = h[i59].im / dcv6[i59].im;
        h[i59].im = 0.0;
      } else if (h[i59].im == 0.0) {
        h[i59].re = 0.0;
        h[i59].im = -(h_re / dcv6[i59].im);
      } else {
        h[i59].re = h[i59].im / dcv6[i59].im;
        h[i59].im = -(h_re / dcv6[i59].im);
      }
    } else {
      brm = fabs(dcv6[i59].re);
      bim = fabs(dcv6[i59].im);
      if (brm > bim) {
        bim = dcv6[i59].im / dcv6[i59].re;
        d = dcv6[i59].re + bim * dcv6[i59].im;
        h[i59].re = (h[i59].re + bim * h[i59].im) / d;
        h[i59].im = (h[i59].im - bim * h_re) / d;
      } else if (bim == brm) {
        if (dcv6[i59].re > 0.0) {
          bim = 0.5;
        } else {
          bim = -0.5;
        }

        if (dcv6[i59].im > 0.0) {
          d = 0.5;
        } else {
          d = -0.5;
        }

        h[i59].re = (h[i59].re * bim + h[i59].im * d) / brm;
        h[i59].im = (h[i59].im * bim - h_re * d) / brm;
      } else {
        bim = dcv6[i59].re / dcv6[i59].im;
        d = dcv6[i59].im + bim * dcv6[i59].re;
        h[i59].re = (bim * h[i59].re + h[i59].im) / d;
        h[i59].im = (bim * h[i59].im - h_re) / d;
      }
    }
  }
}

/*
 * FREQZ_CG Frequency response of digital filter with codegen support
 *
 *  This function is based on 'freqz' by The MathWorks Inc.
 * Arguments    : const double b[60]
 *                const double w[2048]
 *                double Fs
 *                creal_T hh[2048]
 * Return Type  : void
 */
static void f_freqz_cg(const double b[60], const double w[2048], double Fs,
  creal_T hh[2048])
{
  static struct_T options;
  int i17;
  static const char cv22[8] = { 'o', 'n', 'e', 's', 'i', 'd', 'e', 'd' };

  double b_w[2048];
  static const char cv23[7] = { 'o', 'm', 'i', 't', 't', 'e', 'd' };

  /*  Cast to enforce precision rules */
  options.Fs = Fs;
  memcpy(&options.nfft[0], &w[0], sizeof(double) << 11);
  memcpy(&options.w[0], &w[0], sizeof(double) << 11);

  /*  Remaining are default or for advanced use */
  options.fvflag = 1.0;
  for (i17 = 0; i17 < 8; i17++) {
    options.range[i17] = cv22[i17];
  }

  options.centerdc = 0.0;
  for (i17 = 0; i17 < 7; i17++) {
    options.configlevel[i17] = cv23[i17];
  }

  f_firfreqz(b, &options, hh, b_w);

  /*  Generate the default structure to pass to freqzplot */
  /*  If rad/sample, Fs is empty */
}

/*
 * Arguments    : const double p[14]
 *                const creal_T x[2048]
 *                creal_T y[2048]
 * Return Type  : void
 */
static void f_polyval(const double p[14], const creal_T x[2048], creal_T y[2048])
{
  int i18;
  int k;
  double x_im;
  for (i18 = 0; i18 < 2048; i18++) {
    y[i18].re = p[0];
    y[i18].im = 0.0;
  }

  for (k = 0; k < 13; k++) {
    for (i18 = 0; i18 < 2048; i18++) {
      x_im = x[i18].re * y[i18].im + x[i18].im * y[i18].re;
      y[i18].re = (x[i18].re * y[i18].re - x[i18].im * y[i18].im) + p[k + 1];
      y[i18].im = x_im;
    }
  }
}

/*
 * Arguments    : const double o[7]
 *                double u[14]
 * Return Type  : void
 */
static void f_us(const double o[7], double u[14])
{
  int ix;
  int iy;
  int k;
  memset(&u[0], 0, 14U * sizeof(double));
  ix = 0;
  iy = 0;
  for (k = 0; k < 7; k++) {
    u[iy] = o[ix];
    ix++;
    iy += 2;
  }
}

/*
 * Arguments    : FILE * *f
 *                boolean_T *a
 * Return Type  : void
 */
static void fileManager(FILE * *f, boolean_T *a)
{
  *f = stdout;
  *a = true;
}

/*
 * Make b a row
 * Arguments    : const struct_T *options
 *                creal_T h[2048]
 *                double w[2048]
 * Return Type  : void
 */
static void firfreqz(const struct_T *options, creal_T h[2048], double w[2048])
{
  int i54;
  double brm;
  double bim;
  double d;

  /* -------------------------------------------------------------------------- */
  /*  Actual Frequency Response Computation */
  /* if fvflag, */
  /*    Frequency vector specified.  Use Horner's method of polynomial */
  /*    evaluation at the frequency points and divide the numerator */
  /*    by the denominator. */
  /*  */
  /*    Note: we use positive i here because of the relationship */
  /*             polyval(a,exp(1i*w)) = fft(a).*exp(1i*w*(length(a)-1)) */
  /*                ( assuming w = 2*pi*(0:length(a)-1)/length(a) ) */
  /*  */
  /*  Fs was specified, freq. vector is in Hz */
  /*  Convert from Hz to rad/sample for computational purposes */
  /*  Digital frequency must be used for this calculation */
  for (i54 = 0; i54 < 2048; i54++) {
    w[i54] = options->w[i54];
    h[i54].re = 0.0 * (6.2831853071795862 * options->w[i54] / options->Fs * 0.0);
    h[i54].im = 0.0 * (6.2831853071795862 * options->w[i54] / options->Fs);
  }

  b_exp(h);
  for (i54 = 0; i54 < 2048; i54++) {
    if (h[i54].im == 0.0) {
      h[i54].re = 1.0 / h[i54].re;
      h[i54].im = 0.0;
    } else if (h[i54].re == 0.0) {
      h[i54].re = 0.0;
      h[i54].im = -(1.0 / h[i54].im);
    } else {
      brm = fabs(h[i54].re);
      bim = fabs(h[i54].im);
      if (brm > bim) {
        bim = h[i54].im / h[i54].re;
        d = h[i54].re + bim * h[i54].im;
        h[i54].re = (1.0 + bim * 0.0) / d;
        h[i54].im = (0.0 - bim) / d;
      } else if (bim == brm) {
        if (h[i54].re > 0.0) {
          bim = 0.5;
        } else {
          bim = -0.5;
        }

        if (h[i54].im > 0.0) {
          d = 0.5;
        } else {
          d = -0.5;
        }

        h[i54].re = (bim + 0.0 * d) / brm;
        h[i54].im = (0.0 * bim - d) / brm;
      } else {
        bim = h[i54].re / h[i54].im;
        d = h[i54].im + bim * h[i54].re;
        h[i54].re = bim / d;
        h[i54].im = (bim * 0.0 - 1.0) / d;
      }
    }
  }
}

/*
 * FIRPM Parks-McClellan optimal equiripple FIR filter design.
 *
 *  This function is based on 'firpm' by The MathWorks Inc.
 * Arguments    : double order
 *                const double ff[4]
 *                const emxArray_real_T *amplitudes
 *                const emxArray_real_T *frequencies
 *                const emxArray_real_T *weights
 *                emxArray_real_T *h
 * Return Type  : void
 */
static void firpm_cg(double order, const double ff[4], const emxArray_real_T
                     *amplitudes, const emxArray_real_T *frequencies, const
                     emxArray_real_T *weights, emxArray_real_T *h)
{
  emxArray_real_T *grid;
  emxArray_real_T *des;
  emxArray_real_T *wt;
  emxArray_real_T *r20;
  double b_ff[4];
  int i36;
  emxArray_real_T *b_h;
  double err;
  boolean_T valid;
  int h_idx_0;
  int i37;
  emxArray_real_T *c_h;
  int i38;
  int loop_ub;
  emxArray_real_T *d_h;
  emxInit_real_T(&grid, 2);
  emxInit_real_T(&des, 2);
  emxInit_real_T(&wt, 2);
  emxInit_real_T(&r20, 2);
  genWeights(order, ff, frequencies, weights, amplitudes, grid, des, wt);

  /*  Workaround */
  /* ftype = 2; */
  /* sign_val = 1; */
  /*  Always bandpass designs */
  /*  cast to enforce precision rules */
  /*  Call actual design algorithm */
  rdivide(grid, 2.0, r20);
  emxFree_real_T(&grid);
  for (i36 = 0; i36 < 4; i36++) {
    b_ff[i36] = ff[i36] / 2.0;
  }

  emxInit_real_T(&b_h, 2);
  remezm(order + 1.0, b_ff, r20, des, wt, b_h, &err, &valid);
  h_idx_0 = b_h->size[0] * b_h->size[1];
  i36 = h->size[0] * h->size[1];
  h->size[0] = 1;
  h->size[1] = h_idx_0;
  emxEnsureCapacity((emxArray__common *)h, i36, (int)sizeof(double));
  emxFree_real_T(&r20);
  emxFree_real_T(&wt);
  emxFree_real_T(&des);
  for (i36 = 0; i36 < h_idx_0; i36++) {
    h->data[h->size[0] * i36] = b_h->data[i36];
  }

  emxFree_real_T(&b_h);

  /*  make it a row */
  err = (double)h->size[1] - rt_remd_snf(order + 1.0, 2.0);
  if (1.0 > err) {
    i36 = 1;
    h_idx_0 = 1;
    i37 = 0;
  } else {
    i36 = (int)err;
    h_idx_0 = -1;
    i37 = 1;
  }

  emxInit_real_T(&c_h, 2);
  i38 = c_h->size[0] * c_h->size[1];
  c_h->size[0] = 1;
  c_h->size[1] = (h->size[1] + div_s32_floor(i37 - i36, h_idx_0)) + 1;
  emxEnsureCapacity((emxArray__common *)c_h, i38, (int)sizeof(double));
  loop_ub = h->size[1];
  for (i38 = 0; i38 < loop_ub; i38++) {
    c_h->data[c_h->size[0] * i38] = h->data[h->size[0] * i38];
  }

  loop_ub = div_s32_floor(i37 - i36, h_idx_0);
  for (i37 = 0; i37 <= loop_ub; i37++) {
    c_h->data[c_h->size[0] * (i37 + h->size[1])] = h->data[(i36 + h_idx_0 * i37)
      - 1];
  }

  i36 = h->size[0] * h->size[1];
  h->size[0] = 1;
  h->size[1] = c_h->size[1];
  emxEnsureCapacity((emxArray__common *)h, i36, (int)sizeof(double));
  loop_ub = c_h->size[1];
  for (i36 = 0; i36 < loop_ub; i36++) {
    h->data[h->size[0] * i36] = c_h->data[c_h->size[0] * i36];
  }

  emxFree_real_T(&c_h);
  if (1 > h->size[1]) {
    i36 = 1;
    h_idx_0 = 1;
    i37 = 0;
  } else {
    i36 = h->size[1];
    h_idx_0 = -1;
    i37 = 1;
  }

  emxInit_real_T(&d_h, 2);
  i38 = d_h->size[0] * d_h->size[1];
  d_h->size[0] = 1;
  d_h->size[1] = div_s32_floor(i37 - i36, h_idx_0) + 1;
  emxEnsureCapacity((emxArray__common *)d_h, i38, (int)sizeof(double));
  loop_ub = div_s32_floor(i37 - i36, h_idx_0);
  for (i37 = 0; i37 <= loop_ub; i37++) {
    d_h->data[d_h->size[0] * i37] = h->data[(i36 + h_idx_0 * i37) - 1];
  }

  i36 = h->size[0] * h->size[1];
  h->size[0] = 1;
  h->size[1] = d_h->size[1];
  emxEnsureCapacity((emxArray__common *)h, i36, (int)sizeof(double));
  loop_ub = d_h->size[1];
  for (i36 = 0; i36 < loop_ub; i36++) {
    h->data[h->size[0] * i36] = d_h->data[d_h->size[0] * i36];
  }

  emxFree_real_T(&d_h);
}

/*
 * firpmgrid
 * Arguments    : double nfilt
 *                const double ff[4]
 *                emxArray_real_T *gridactual
 * Return Type  : void
 */
static void firpmgrid_cg(double nfilt, const double ff[4], emxArray_real_T
  *gridactual)
{
  static double grid[100000];
  double ngrid;
  double b_ngrid;
  double delf;
  double j;
  double l;
  double gridSize;
  emxArray_real_T *newgrid;
  emxArray_int32_T *r21;
  double a;
  int k;
  int nm1d2;
  double ndbl;
  double apnd;
  double cdiff;
  double delf1;
  double absa;
  double absb;
  int n;

  /* -------------------------------------------------------------------------- */
  memset(&grid[0], 0, 100000U * sizeof(double));

  /*  Make large initial memory */
  /*     Generate frequency grid */
  ngrid = nfilt / 2.0;
  grid[0] = ff[0];
  if (ngrid < 0.0) {
    b_ngrid = ceil(ngrid);
  } else {
    b_ngrid = floor(ngrid);
  }

  delf = 1.0 / (16.0 * b_ngrid);

  /*  If value at frequency 0 is constrained, make sure first grid point */
  /*  is not too small: */
  j = 1.0;
  l = 1.0;
  gridSize = 1.0;
  emxInit_real_T(&newgrid, 2);
  emxInit_int32_T(&r21, 2);
  while (l + 1.0 <= 4.0) {
    a = grid[(int)j - 1] + delf;
    ngrid = ff[(int)(l + 1.0) - 1] + delf;
    if (rtIsNaN(a) || rtIsNaN(delf) || rtIsNaN(ngrid)) {
      k = newgrid->size[0] * newgrid->size[1];
      newgrid->size[0] = 1;
      newgrid->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)newgrid, k, (int)sizeof(double));
      newgrid->data[0] = rtNaN;
    } else if ((delf == 0.0) || ((a < ngrid) && (delf < 0.0)) || ((ngrid < a) &&
                (delf > 0.0))) {
      k = newgrid->size[0] * newgrid->size[1];
      newgrid->size[0] = 1;
      newgrid->size[1] = 0;
      emxEnsureCapacity((emxArray__common *)newgrid, k, (int)sizeof(double));
    } else if ((rtIsInf(a) || rtIsInf(ngrid)) && (rtIsInf(delf) || (a == ngrid)))
    {
      k = newgrid->size[0] * newgrid->size[1];
      newgrid->size[0] = 1;
      newgrid->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)newgrid, k, (int)sizeof(double));
      newgrid->data[0] = rtNaN;
    } else if (rtIsInf(delf)) {
      k = newgrid->size[0] * newgrid->size[1];
      newgrid->size[0] = 1;
      newgrid->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)newgrid, k, (int)sizeof(double));
      newgrid->data[0] = a;
    } else if ((floor(a) == a) && (floor(delf) == delf)) {
      k = newgrid->size[0] * newgrid->size[1];
      newgrid->size[0] = 1;
      newgrid->size[1] = (int)floor((ngrid - a) / delf) + 1;
      emxEnsureCapacity((emxArray__common *)newgrid, k, (int)sizeof(double));
      nm1d2 = (int)floor((ngrid - a) / delf);
      for (k = 0; k <= nm1d2; k++) {
        newgrid->data[newgrid->size[0] * k] = a + delf * (double)k;
      }
    } else {
      ndbl = floor((ngrid - a) / delf + 0.5);
      apnd = a + ndbl * delf;
      if (delf > 0.0) {
        cdiff = apnd - ngrid;
      } else {
        cdiff = ngrid - apnd;
      }

      absa = fabs(a);
      absb = fabs(ngrid);
      if ((absa >= absb) || rtIsNaN(absb)) {
        absb = absa;
      }

      if (fabs(cdiff) < 4.4408920985006262E-16 * absb) {
        ndbl++;
        apnd = ngrid;
      } else if (cdiff > 0.0) {
        apnd = a + (ndbl - 1.0) * delf;
      } else {
        ndbl++;
      }

      if (ndbl >= 0.0) {
        n = (int)ndbl;
      } else {
        n = 0;
      }

      k = newgrid->size[0] * newgrid->size[1];
      newgrid->size[0] = 1;
      newgrid->size[1] = n;
      emxEnsureCapacity((emxArray__common *)newgrid, k, (int)sizeof(double));
      if (n > 0) {
        newgrid->data[0] = a;
        if (n > 1) {
          newgrid->data[n - 1] = apnd;
          nm1d2 = (n - 1) / 2;
          for (k = 1; k < nm1d2; k++) {
            ngrid = (double)k * delf;
            newgrid->data[k] = a + ngrid;
            newgrid->data[(n - k) - 1] = apnd - ngrid;
          }

          if (nm1d2 << 1 == n - 1) {
            newgrid->data[nm1d2] = (a + apnd) / 2.0;
          } else {
            ngrid = (double)nm1d2 * delf;
            newgrid->data[nm1d2] = a + ngrid;
            newgrid->data[nm1d2 + 1] = apnd - ngrid;
          }
        }
      }
    }

    if (newgrid->size[1] < 11) {
      delf1 = ((ff[(int)(l + 1.0) - 1] + delf) - (grid[(int)j - 1] + delf)) /
        10.0;
      a = grid[(int)j - 1] + delf1;
      ngrid = ff[(int)(l + 1.0) - 1] + delf1;
      if (rtIsNaN(a) || rtIsNaN(delf1) || rtIsNaN(ngrid)) {
        k = newgrid->size[0] * newgrid->size[1];
        newgrid->size[0] = 1;
        newgrid->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)newgrid, k, (int)sizeof(double));
        newgrid->data[0] = rtNaN;
      } else if ((delf1 == 0.0) || ((a < ngrid) && (delf1 < 0.0)) || ((ngrid < a)
                  && (delf1 > 0.0))) {
        k = newgrid->size[0] * newgrid->size[1];
        newgrid->size[0] = 1;
        newgrid->size[1] = 0;
        emxEnsureCapacity((emxArray__common *)newgrid, k, (int)sizeof(double));
      } else if ((rtIsInf(a) || rtIsInf(ngrid)) && (rtIsInf(delf1) || (a ==
                   ngrid))) {
        k = newgrid->size[0] * newgrid->size[1];
        newgrid->size[0] = 1;
        newgrid->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)newgrid, k, (int)sizeof(double));
        newgrid->data[0] = rtNaN;
      } else if (rtIsInf(delf1)) {
        k = newgrid->size[0] * newgrid->size[1];
        newgrid->size[0] = 1;
        newgrid->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)newgrid, k, (int)sizeof(double));
        newgrid->data[0] = a;
      } else if ((floor(a) == a) && (floor(delf1) == delf1)) {
        k = newgrid->size[0] * newgrid->size[1];
        newgrid->size[0] = 1;
        newgrid->size[1] = (int)floor((ngrid - a) / delf1) + 1;
        emxEnsureCapacity((emxArray__common *)newgrid, k, (int)sizeof(double));
        nm1d2 = (int)floor((ngrid - a) / delf1);
        for (k = 0; k <= nm1d2; k++) {
          newgrid->data[newgrid->size[0] * k] = a + delf1 * (double)k;
        }
      } else {
        ndbl = floor((ngrid - a) / delf1 + 0.5);
        apnd = a + ndbl * delf1;
        if (delf1 > 0.0) {
          cdiff = apnd - ngrid;
        } else {
          cdiff = ngrid - apnd;
        }

        absa = fabs(a);
        absb = fabs(ngrid);
        if ((absa >= absb) || rtIsNaN(absb)) {
          absb = absa;
        }

        if (fabs(cdiff) < 4.4408920985006262E-16 * absb) {
          ndbl++;
          apnd = ngrid;
        } else if (cdiff > 0.0) {
          apnd = a + (ndbl - 1.0) * delf1;
        } else {
          ndbl++;
        }

        if (ndbl >= 0.0) {
          n = (int)ndbl;
        } else {
          n = 0;
        }

        k = newgrid->size[0] * newgrid->size[1];
        newgrid->size[0] = 1;
        newgrid->size[1] = n;
        emxEnsureCapacity((emxArray__common *)newgrid, k, (int)sizeof(double));
        if (n > 0) {
          newgrid->data[0] = a;
          if (n > 1) {
            newgrid->data[n - 1] = apnd;
            nm1d2 = (n - 1) / 2;
            for (k = 1; k < nm1d2; k++) {
              ngrid = (double)k * delf1;
              newgrid->data[k] = a + ngrid;
              newgrid->data[(n - k) - 1] = apnd - ngrid;
            }

            if (nm1d2 << 1 == n - 1) {
              newgrid->data[nm1d2] = (a + apnd) / 2.0;
            } else {
              ngrid = (double)nm1d2 * delf1;
              newgrid->data[nm1d2] = a + ngrid;
              newgrid->data[nm1d2 + 1] = apnd - ngrid;
            }
          }
        }
      }
    }

    /* grid = [grid newgrid]; */
    k = newgrid->size[1];
    nm1d2 = r21->size[0] * r21->size[1];
    r21->size[0] = 1;
    r21->size[1] = (int)((double)k - 1.0) + 1;
    emxEnsureCapacity((emxArray__common *)r21, nm1d2, (int)sizeof(int));
    nm1d2 = (int)((double)k - 1.0);
    for (k = 0; k <= nm1d2; k++) {
      r21->data[r21->size[0] * k] = (int)(gridSize + (1.0 + (double)k));
    }

    nm1d2 = newgrid->size[0] * newgrid->size[1];
    for (k = 0; k < nm1d2; k++) {
      grid[r21->data[k] - 1] = newgrid->data[k];
    }

    gridSize += (double)newgrid->size[1];

    /* jend = length(grid); */
    /* length(grid); */
    if (gridSize > 1.0) {
      grid[(int)(gridSize - 1.0) - 1] = ff[(int)(l + 1.0) - 1];
      j = gridSize;
    } else {
      j = 2.0;
    }

    l += 2.0;
    if (l + 1.0 <= 4.0) {
      grid[(int)j - 1] = ff[(int)l - 1];
    }
  }

  emxFree_int32_T(&r21);
  emxFree_real_T(&newgrid);
  ngrid = j - 1.0;

  /*  If value at frequency 1 is constrained, remove that grid point: */
  if (grid[(int)(j - 1.0) - 1] > 1.0 - delf) {
    if (ff[2] < 1.0 - delf) {
      ngrid = (j - 1.0) - 1.0;
    } else {
      grid[(int)(j - 1.0) - 1] = ff[2];
    }
  }

  if (1.0 > ngrid) {
    nm1d2 = 0;
  } else {
    nm1d2 = (int)ngrid;
  }

  k = gridactual->size[0] * gridactual->size[1];
  gridactual->size[0] = 1;
  gridactual->size[1] = nm1d2;
  emxEnsureCapacity((emxArray__common *)gridactual, k, (int)sizeof(double));
  for (k = 0; k < nm1d2; k++) {
    gridactual->data[gridactual->size[0] * k] = grid[k];
  }
}

/*
 * FREQS Laplace-transform (s-domain) frequency response with codegen support
 *
 *  This function is based on 'freqs' by The MathWorks Inc.
 * Arguments    : const double b_data[]
 *                const int b_size[2]
 *                const emxArray_creal_T *a
 *                const double w[2048]
 *                creal_T h[2048]
 * Return Type  : void
 */
static void freqs_cg(const double b_data[], const int b_size[2], const
                     emxArray_creal_T *a, const double w[2048], creal_T h[2048])
{
  emxArray_creal_T *b_a;
  double b_b_data[4];
  int b_b_size[2];
  static creal_T s[2048];
  int i21;
  creal_T y[2048];
  boolean_T b0;
  int k;
  double h_re;
  double bim;
  double d;
  double brm;
  emxInit_creal_T(&b_a, 2);
  removeTrailingZero(b_data, b_size, a, b_b_data, b_b_size, b_a);
  for (i21 = 0; i21 < 2048; i21++) {
    s[i21].re = w[i21] * 0.0;
    s[i21].im = w[i21];
  }

  b0 = (b_a->size[1] == 0);
  if (!b0) {
    for (i21 = 0; i21 < 2048; i21++) {
      y[i21] = b_a->data[0];
    }

    for (k = 0; k <= b_a->size[1] - 2; k++) {
      bim = b_a->data[k + 1].re;
      d = b_a->data[k + 1].im;
      for (i21 = 0; i21 < 2048; i21++) {
        brm = s[i21].re * y[i21].im + s[i21].im * y[i21].re;
        y[i21].re = (s[i21].re * y[i21].re - s[i21].im * y[i21].im) + bim;
        y[i21].im = brm + d;
      }
    }
  }

  emxFree_creal_T(&b_a);
  c_polyval(b_b_data, b_b_size, s, h);
  for (i21 = 0; i21 < 2048; i21++) {
    h_re = h[i21].re;
    if (y[i21].im == 0.0) {
      if (h[i21].im == 0.0) {
        h[i21].re /= y[i21].re;
        h[i21].im = 0.0;
      } else if (h[i21].re == 0.0) {
        h[i21].re = 0.0;
        h[i21].im /= y[i21].re;
      } else {
        h[i21].re /= y[i21].re;
        h[i21].im /= y[i21].re;
      }
    } else if (y[i21].re == 0.0) {
      if (h[i21].re == 0.0) {
        h[i21].re = h[i21].im / y[i21].im;
        h[i21].im = 0.0;
      } else if (h[i21].im == 0.0) {
        h[i21].re = 0.0;
        h[i21].im = -(h_re / y[i21].im);
      } else {
        h[i21].re = h[i21].im / y[i21].im;
        h[i21].im = -(h_re / y[i21].im);
      }
    } else {
      brm = fabs(y[i21].re);
      bim = fabs(y[i21].im);
      if (brm > bim) {
        bim = y[i21].im / y[i21].re;
        d = y[i21].re + bim * y[i21].im;
        h[i21].re = (h[i21].re + bim * h[i21].im) / d;
        h[i21].im = (h[i21].im - bim * h_re) / d;
      } else if (bim == brm) {
        if (y[i21].re > 0.0) {
          bim = 0.5;
        } else {
          bim = -0.5;
        }

        if (y[i21].im > 0.0) {
          d = 0.5;
        } else {
          d = -0.5;
        }

        h[i21].re = (h[i21].re * bim + h[i21].im * d) / brm;
        h[i21].im = (h[i21].im * bim - h_re * d) / brm;
      } else {
        bim = y[i21].re / y[i21].im;
        d = y[i21].im + bim * y[i21].re;
        h[i21].re = (bim * h[i21].re + h[i21].im) / d;
        h[i21].im = (bim * h[i21].im - h_re) / d;
      }
    }
  }
}

/*
 * FREQZ_CG Frequency response of digital filter with codegen support
 *
 *  This function is based on 'freqz' by The MathWorks Inc.
 * Arguments    : const double w[2048]
 *                double Fs
 *                creal_T hh[2048]
 * Return Type  : void
 */
static void freqz_cg(const double w[2048], double Fs, creal_T hh[2048])
{
  struct_T options;
  int i7;
  static const char cv12[8] = { 'o', 'n', 'e', 's', 'i', 'd', 'e', 'd' };

  double b_w[2048];
  static const char cv13[7] = { 'o', 'm', 'i', 't', 't', 'e', 'd' };

  /*  Cast to enforce precision rules */
  options.Fs = Fs;
  memcpy(&options.nfft[0], &w[0], sizeof(double) << 11);
  memcpy(&options.w[0], &w[0], sizeof(double) << 11);

  /*  Remaining are default or for advanced use */
  options.fvflag = 1.0;
  for (i7 = 0; i7 < 8; i7++) {
    options.range[i7] = cv12[i7];
  }

  options.centerdc = 0.0;
  for (i7 = 0; i7 < 7; i7++) {
    options.configlevel[i7] = cv13[i7];
  }

  firfreqz(&options, hh, b_w);

  /*  Generate the default structure to pass to freqzplot */
  /*  If rad/sample, Fs is empty */
}

/*
 * Make b a row
 * Arguments    : const double b[14]
 *                const struct_T *options
 *                creal_T h[2048]
 *                double w[2048]
 * Return Type  : void
 */
static void g_firfreqz(const double b[14], const struct_T *options, creal_T h
  [2048], double w[2048])
{
  double digw[2048];
  creal_T dcv7[2048];
  int i60;
  double bim;
  double h_re;
  double brm;
  double d;

  /* -------------------------------------------------------------------------- */
  /*  Actual Frequency Response Computation */
  /* if fvflag, */
  /*    Frequency vector specified.  Use Horner's method of polynomial */
  /*    evaluation at the frequency points and divide the numerator */
  /*    by the denominator. */
  /*  */
  /*    Note: we use positive i here because of the relationship */
  /*             polyval(a,exp(1i*w)) = fft(a).*exp(1i*w*(length(a)-1)) */
  /*                ( assuming w = 2*pi*(0:length(a)-1)/length(a) ) */
  /*  */
  /*  Fs was specified, freq. vector is in Hz */
  /*  Convert from Hz to rad/sample for computational purposes */
  /*  Digital frequency must be used for this calculation */
  for (i60 = 0; i60 < 2048; i60++) {
    w[i60] = options->w[i60];
    bim = 6.2831853071795862 * options->w[i60] / options->Fs;
    dcv7[i60].re = bim * 0.0;
    dcv7[i60].im = bim;
    digw[i60] = bim;
  }

  b_exp(dcv7);
  f_polyval(b, dcv7, h);
  for (i60 = 0; i60 < 2048; i60++) {
    dcv7[i60].re = 13.0 * (digw[i60] * 0.0);
    dcv7[i60].im = 13.0 * digw[i60];
  }

  b_exp(dcv7);
  for (i60 = 0; i60 < 2048; i60++) {
    h_re = h[i60].re;
    if (dcv7[i60].im == 0.0) {
      if (h[i60].im == 0.0) {
        h[i60].re /= dcv7[i60].re;
        h[i60].im = 0.0;
      } else if (h[i60].re == 0.0) {
        h[i60].re = 0.0;
        h[i60].im /= dcv7[i60].re;
      } else {
        h[i60].re /= dcv7[i60].re;
        h[i60].im /= dcv7[i60].re;
      }
    } else if (dcv7[i60].re == 0.0) {
      if (h[i60].re == 0.0) {
        h[i60].re = h[i60].im / dcv7[i60].im;
        h[i60].im = 0.0;
      } else if (h[i60].im == 0.0) {
        h[i60].re = 0.0;
        h[i60].im = -(h_re / dcv7[i60].im);
      } else {
        h[i60].re = h[i60].im / dcv7[i60].im;
        h[i60].im = -(h_re / dcv7[i60].im);
      }
    } else {
      brm = fabs(dcv7[i60].re);
      bim = fabs(dcv7[i60].im);
      if (brm > bim) {
        bim = dcv7[i60].im / dcv7[i60].re;
        d = dcv7[i60].re + bim * dcv7[i60].im;
        h[i60].re = (h[i60].re + bim * h[i60].im) / d;
        h[i60].im = (h[i60].im - bim * h_re) / d;
      } else if (bim == brm) {
        if (dcv7[i60].re > 0.0) {
          bim = 0.5;
        } else {
          bim = -0.5;
        }

        if (dcv7[i60].im > 0.0) {
          d = 0.5;
        } else {
          d = -0.5;
        }

        h[i60].re = (h[i60].re * bim + h[i60].im * d) / brm;
        h[i60].im = (h[i60].im * bim - h_re * d) / brm;
      } else {
        bim = dcv7[i60].re / dcv7[i60].im;
        d = dcv7[i60].im + bim * dcv7[i60].re;
        h[i60].re = (bim * h[i60].re + h[i60].im) / d;
        h[i60].im = (bim * h[i60].im - h_re) / d;
      }
    }
  }
}

/*
 * FREQZ_CG Frequency response of digital filter with codegen support
 *
 *  This function is based on 'freqz' by The MathWorks Inc.
 * Arguments    : const double b[14]
 *                const double w[2048]
 *                double Fs
 *                creal_T hh[2048]
 * Return Type  : void
 */
static void g_freqz_cg(const double b[14], const double w[2048], double Fs,
  creal_T hh[2048])
{
  static struct_T options;
  int i19;
  static const char cv24[8] = { 'o', 'n', 'e', 's', 'i', 'd', 'e', 'd' };

  double b_w[2048];
  static const char cv25[7] = { 'o', 'm', 'i', 't', 't', 'e', 'd' };

  /*  Cast to enforce precision rules */
  options.Fs = Fs;
  memcpy(&options.nfft[0], &w[0], sizeof(double) << 11);
  memcpy(&options.w[0], &w[0], sizeof(double) << 11);

  /*  Remaining are default or for advanced use */
  options.fvflag = 1.0;
  for (i19 = 0; i19 < 8; i19++) {
    options.range[i19] = cv24[i19];
  }

  options.centerdc = 0.0;
  for (i19 = 0; i19 < 7; i19++) {
    options.configlevel[i19] = cv25[i19];
  }

  g_firfreqz(b, &options, hh, b_w);

  /*  Generate the default structure to pass to freqzplot */
  /*  If rad/sample, Fs is empty */
}

/*
 * Arguments    : const double p_data[]
 *                const int p_size[2]
 *                const emxArray_creal_T *x
 *                emxArray_creal_T *y
 * Return Type  : void
 */
static void g_polyval(const double p_data[], const int p_size[2], const
                      emxArray_creal_T *x, emxArray_creal_T *y)
{
  int i29;
  boolean_T b4;
  int loop_ub;
  int k;
  double x_re;
  double x_im;
  i29 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)y, i29, (int)sizeof(creal_T));
  if ((y->size[1] == 0) || (p_size[1] == 0)) {
    b4 = true;
  } else {
    b4 = false;
  }

  if (!b4) {
    i29 = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)y, i29, (int)sizeof(creal_T));
    loop_ub = y->size[1];
    for (i29 = 0; i29 < loop_ub; i29++) {
      y->data[y->size[0] * i29].re = p_data[0];
      y->data[y->size[0] * i29].im = 0.0;
    }

    for (k = 0; k <= p_size[1] - 2; k++) {
      i29 = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = x->size[1];
      emxEnsureCapacity((emxArray__common *)y, i29, (int)sizeof(creal_T));
      loop_ub = x->size[0] * x->size[1];
      for (i29 = 0; i29 < loop_ub; i29++) {
        x_re = x->data[i29].re * y->data[i29].re - x->data[i29].im * y->data[i29]
          .im;
        x_im = x->data[i29].re * y->data[i29].im + x->data[i29].im * y->data[i29]
          .re;
        y->data[i29].re = x_re + p_data[k + 1];
        y->data[i29].im = x_im;
      }
    }
  }
}

/*
 * Arguments    : double filterOrder
 *                const double bands[4]
 *                const emxArray_real_T *frequencies
 *                const emxArray_real_T *weights
 *                const emxArray_real_T *amplitudes
 *                emxArray_real_T *grid
 *                emxArray_real_T *des
 *                emxArray_real_T *wt
 * Return Type  : void
 */
static void genWeights(double filterOrder, const double bands[4], const
  emxArray_real_T *frequencies, const emxArray_real_T *weights, const
  emxArray_real_T *amplitudes, emxArray_real_T *grid, emxArray_real_T *des,
  emxArray_real_T *wt)
{
  int iv1[2];
  int ixstart;
  emxArray_uint32_T *positionsOfNewFreqIndx;
  int n;
  int ind;
  emxArray_real_T *varargin_2;
  emxArray_real_T *b_frequencies;
  double b_grid;
  int itmp;
  int ix;
  boolean_T exitg1;

  /*  */
  firpmgrid_cg(filterOrder + 1.0, bands, grid);
  for (ixstart = 0; ixstart < 2; ixstart++) {
    iv1[ixstart] = grid->size[ixstart];
  }

  emxInit_uint32_T(&positionsOfNewFreqIndx, 2);
  ixstart = positionsOfNewFreqIndx->size[0] * positionsOfNewFreqIndx->size[1];
  positionsOfNewFreqIndx->size[0] = 1;
  positionsOfNewFreqIndx->size[1] = iv1[1];
  emxEnsureCapacity((emxArray__common *)positionsOfNewFreqIndx, ixstart, (int)
                    sizeof(unsigned int));
  n = iv1[1];
  for (ixstart = 0; ixstart < n; ixstart++) {
    positionsOfNewFreqIndx->data[ixstart] = 0U;
  }

  ind = 0;
  emxInit_real_T(&varargin_2, 2);
  emxInit_real_T(&b_frequencies, 2);
  while (ind <= iv1[1] - 1) {
    ixstart = b_frequencies->size[0] * b_frequencies->size[1];
    b_frequencies->size[0] = 1;
    b_frequencies->size[1] = frequencies->size[1];
    emxEnsureCapacity((emxArray__common *)b_frequencies, ixstart, (int)sizeof
                      (double));
    b_grid = grid->data[ind];
    n = frequencies->size[0] * frequencies->size[1];
    for (ixstart = 0; ixstart < n; ixstart++) {
      b_frequencies->data[ixstart] = frequencies->data[ixstart] - b_grid;
    }

    c_abs(b_frequencies, varargin_2);
    ixstart = 1;
    n = varargin_2->size[1];
    b_grid = varargin_2->data[0];
    itmp = 1;
    if (varargin_2->size[1] > 1) {
      if (rtIsNaN(varargin_2->data[0])) {
        ix = 2;
        exitg1 = false;
        while ((!exitg1) && (ix <= n)) {
          ixstart = ix;
          if (!rtIsNaN(varargin_2->data[ix - 1])) {
            b_grid = varargin_2->data[ix - 1];
            itmp = ix;
            exitg1 = true;
          } else {
            ix++;
          }
        }
      }

      if (ixstart < varargin_2->size[1]) {
        while (ixstart + 1 <= n) {
          if (varargin_2->data[ixstart] < b_grid) {
            b_grid = varargin_2->data[ixstart];
            itmp = ixstart + 1;
          }

          ixstart++;
        }
      }
    }

    positionsOfNewFreqIndx->data[ind] = (unsigned int)itmp;
    ind++;
  }

  emxFree_real_T(&b_frequencies);
  emxFree_real_T(&varargin_2);
  ixstart = wt->size[0] * wt->size[1];
  wt->size[0] = 1;
  wt->size[1] = positionsOfNewFreqIndx->size[1];
  emxEnsureCapacity((emxArray__common *)wt, ixstart, (int)sizeof(double));
  n = positionsOfNewFreqIndx->size[0] * positionsOfNewFreqIndx->size[1];
  for (ixstart = 0; ixstart < n; ixstart++) {
    wt->data[ixstart] = weights->data[(int)positionsOfNewFreqIndx->data[ixstart]
      - 1];
  }

  ixstart = des->size[0] * des->size[1];
  des->size[0] = 1;
  des->size[1] = positionsOfNewFreqIndx->size[1];
  emxEnsureCapacity((emxArray__common *)des, ixstart, (int)sizeof(double));
  n = positionsOfNewFreqIndx->size[0] * positionsOfNewFreqIndx->size[1];
  for (ixstart = 0; ixstart < n; ixstart++) {
    des->data[ixstart] = amplitudes->data[(int)positionsOfNewFreqIndx->
      data[ixstart] - 1];
  }

  emxFree_uint32_T(&positionsOfNewFreqIndx);
}

/*
 * Arguments    : const char enables[4]
 *                const double w[2048]
 *                double Fs
 *                const double hb1_coeff[15]
 *                const double hb2_coeff[7]
 *                const double hb3_coeff_data[]
 *                const int hb3_coeff_size[2]
 *                const double dec_int3_coeff_data[]
 *                const int dec_int3_coeff_size[2]
 *                creal_T combinedResponse[2048]
 * Return Type  : void
 */
static void generateCascadedResponseRx(const char enables[4], const double w
  [2048], double Fs, const double hb1_coeff[15], const double hb2_coeff[7],
  const double hb3_coeff_data[], const int hb3_coeff_size[2], const double
  dec_int3_coeff_data[], const int dec_int3_coeff_size[2], creal_T
  combinedResponse[2048])
{
  boolean_T b_bool;
  int kstr;
  int exitg12;
  int exitg11;
  double dv2[15];
  double dv3[7];
  double tmp_data[29];
  int tmp_size[2];
  double dv4[30];
  double dv5[30];
  double dv6[30];
  double dv7[60];
  double dv8[30];
  double dv9[14];
  double dv10[60];
  static const char cv1[4] = { '2', '1', '1', '1' };

  double dv11[7];
  double dv12[7];
  double dv13[14];
  double dv14[14];
  static creal_T d2[2048];
  double combinedResponse_re;
  static creal_T d3[2048];
  int exitg10;
  static const char cv2[4] = { '1', '2', '1', '1' };

  double combinedResponse_im;
  int exitg9;
  static const char cv3[4] = { '1', '1', '2', '1' };

  int exitg8;
  static const char cv4[4] = { '2', '2', '1', '1' };

  int exitg7;
  static const char cv5[4] = { '2', '1', '2', '1' };

  int exitg6;
  static const char cv6[4] = { '1', '2', '2', '1' };

  int exitg5;
  static const char cv7[4] = { '2', '2', '2', '1' };

  int exitg4;
  static const char cv8[4] = { '1', '1', '1', '3' };

  int exitg3;
  static const char cv9[4] = { '2', '1', '1', '3' };

  int exitg2;
  static const char cv10[4] = { '1', '2', '1', '3' };

  int exitg1;
  static const char cv11[4] = { '2', '2', '1', '3' };

  /*  Cast */
  b_bool = false;
  kstr = 1;
  do {
    exitg12 = 0;
    if (kstr < 5) {
      if (enables[kstr - 1] != '1') {
        exitg12 = 1;
      } else {
        kstr++;
      }
    } else {
      b_bool = true;
      exitg12 = 1;
    }
  } while (exitg12 == 0);

  if (b_bool) {
    kstr = 0;
  } else {
    b_bool = false;
    kstr = 0;
    do {
      exitg11 = 0;
      if (kstr + 1 < 5) {
        if (enables[kstr] != cv1[kstr]) {
          exitg11 = 1;
        } else {
          kstr++;
        }
      } else {
        b_bool = true;
        exitg11 = 1;
      }
    } while (exitg11 == 0);

    if (b_bool) {
      kstr = 1;
    } else {
      b_bool = false;
      kstr = 0;
      do {
        exitg10 = 0;
        if (kstr + 1 < 5) {
          if (enables[kstr] != cv2[kstr]) {
            exitg10 = 1;
          } else {
            kstr++;
          }
        } else {
          b_bool = true;
          exitg10 = 1;
        }
      } while (exitg10 == 0);

      if (b_bool) {
        kstr = 2;
      } else {
        b_bool = false;
        kstr = 0;
        do {
          exitg9 = 0;
          if (kstr + 1 < 5) {
            if (enables[kstr] != cv3[kstr]) {
              exitg9 = 1;
            } else {
              kstr++;
            }
          } else {
            b_bool = true;
            exitg9 = 1;
          }
        } while (exitg9 == 0);

        if (b_bool) {
          kstr = 3;
        } else {
          b_bool = false;
          kstr = 0;
          do {
            exitg8 = 0;
            if (kstr + 1 < 5) {
              if (enables[kstr] != cv4[kstr]) {
                exitg8 = 1;
              } else {
                kstr++;
              }
            } else {
              b_bool = true;
              exitg8 = 1;
            }
          } while (exitg8 == 0);

          if (b_bool) {
            kstr = 4;
          } else {
            b_bool = false;
            kstr = 0;
            do {
              exitg7 = 0;
              if (kstr + 1 < 5) {
                if (enables[kstr] != cv5[kstr]) {
                  exitg7 = 1;
                } else {
                  kstr++;
                }
              } else {
                b_bool = true;
                exitg7 = 1;
              }
            } while (exitg7 == 0);

            if (b_bool) {
              kstr = 5;
            } else {
              b_bool = false;
              kstr = 0;
              do {
                exitg6 = 0;
                if (kstr + 1 < 5) {
                  if (enables[kstr] != cv6[kstr]) {
                    exitg6 = 1;
                  } else {
                    kstr++;
                  }
                } else {
                  b_bool = true;
                  exitg6 = 1;
                }
              } while (exitg6 == 0);

              if (b_bool) {
                kstr = 6;
              } else {
                b_bool = false;
                kstr = 0;
                do {
                  exitg5 = 0;
                  if (kstr + 1 < 5) {
                    if (enables[kstr] != cv7[kstr]) {
                      exitg5 = 1;
                    } else {
                      kstr++;
                    }
                  } else {
                    b_bool = true;
                    exitg5 = 1;
                  }
                } while (exitg5 == 0);

                if (b_bool) {
                  kstr = 7;
                } else {
                  b_bool = false;
                  kstr = 0;
                  do {
                    exitg4 = 0;
                    if (kstr + 1 < 5) {
                      if (enables[kstr] != cv8[kstr]) {
                        exitg4 = 1;
                      } else {
                        kstr++;
                      }
                    } else {
                      b_bool = true;
                      exitg4 = 1;
                    }
                  } while (exitg4 == 0);

                  if (b_bool) {
                    kstr = 8;
                  } else {
                    b_bool = false;
                    kstr = 0;
                    do {
                      exitg3 = 0;
                      if (kstr + 1 < 5) {
                        if (enables[kstr] != cv9[kstr]) {
                          exitg3 = 1;
                        } else {
                          kstr++;
                        }
                      } else {
                        b_bool = true;
                        exitg3 = 1;
                      }
                    } while (exitg3 == 0);

                    if (b_bool) {
                      kstr = 9;
                    } else {
                      b_bool = false;
                      kstr = 0;
                      do {
                        exitg2 = 0;
                        if (kstr + 1 < 5) {
                          if (enables[kstr] != cv10[kstr]) {
                            exitg2 = 1;
                          } else {
                            kstr++;
                          }
                        } else {
                          b_bool = true;
                          exitg2 = 1;
                        }
                      } while (exitg2 == 0);

                      if (b_bool) {
                        kstr = 10;
                      } else {
                        b_bool = false;
                        kstr = 0;
                        do {
                          exitg1 = 0;
                          if (kstr + 1 < 5) {
                            if (enables[kstr] != cv11[kstr]) {
                              exitg1 = 1;
                            } else {
                              kstr++;
                            }
                          } else {
                            b_bool = true;
                            exitg1 = 1;
                          }
                        } while (exitg1 == 0);

                        if (b_bool) {
                          kstr = 11;
                        } else {
                          kstr = -1;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  switch (kstr) {
   case 0:
    /*  only FIR */
    freqz_cg(w, Fs, combinedResponse);
    break;

   case 1:
    /*  Hb1 */
    us(hb1_coeff, dv2);
    b_freqz_cg(dv2, w, Fs, combinedResponse);
    break;

   case 2:
    /*  Hb2 */
    b_us(hb2_coeff, dv3);
    c_freqz_cg(dv3, w, Fs, combinedResponse);
    break;

   case 3:
    /*  Hb3 */
    c_us(hb3_coeff_data, hb3_coeff_size, tmp_data, tmp_size);
    d_freqz_cg(tmp_data, tmp_size, w, Fs, combinedResponse);
    break;

   case 4:
    /*  Hb2,Hb1 */
    d_us(hb1_coeff, dv4);
    e_freqz_cg(dv4, w, Fs, combinedResponse);
    b_us(hb2_coeff, dv11);
    c_freqz_cg(dv11, w, Fs, d2);
    for (kstr = 0; kstr < 2048; kstr++) {
      combinedResponse_re = combinedResponse[kstr].re;
      combinedResponse[kstr].re = combinedResponse[kstr].re * d2[kstr].re -
        combinedResponse[kstr].im * d2[kstr].im;
      combinedResponse[kstr].im = combinedResponse_re * d2[kstr].im +
        combinedResponse[kstr].im * d2[kstr].re;
    }
    break;

   case 5:
    /*  Hb3,Hb1 */
    d_us(hb1_coeff, dv5);
    e_freqz_cg(dv5, w, Fs, combinedResponse);
    b_us(hb2_coeff, dv12);
    c_freqz_cg(dv12, w, Fs, d2);
    for (kstr = 0; kstr < 2048; kstr++) {
      combinedResponse_re = combinedResponse[kstr].re;
      combinedResponse[kstr].re = combinedResponse[kstr].re * d2[kstr].re -
        combinedResponse[kstr].im * d2[kstr].im;
      combinedResponse[kstr].im = combinedResponse_re * d2[kstr].im +
        combinedResponse[kstr].im * d2[kstr].re;
    }
    break;

   case 6:
    /*  Hb3,Hb2 */
    d_us(hb1_coeff, dv6);
    e_freqz_cg(dv6, w, Fs, combinedResponse);
    c_us(hb3_coeff_data, hb3_coeff_size, tmp_data, tmp_size);
    d_freqz_cg(tmp_data, tmp_size, w, Fs, d2);
    for (kstr = 0; kstr < 2048; kstr++) {
      combinedResponse_re = combinedResponse[kstr].re;
      combinedResponse[kstr].re = combinedResponse[kstr].re * d2[kstr].re -
        combinedResponse[kstr].im * d2[kstr].im;
      combinedResponse[kstr].im = combinedResponse_re * d2[kstr].im +
        combinedResponse[kstr].im * d2[kstr].re;
    }
    break;

   case 7:
    /*  Hb3,Hb2,Hb1 */
    e_us(hb1_coeff, dv7);
    f_freqz_cg(dv7, w, Fs, combinedResponse);
    f_us(hb2_coeff, dv13);
    g_freqz_cg(dv13, w, Fs, d2);
    c_us(hb3_coeff_data, hb3_coeff_size, tmp_data, tmp_size);
    d_freqz_cg(tmp_data, tmp_size, w, Fs, d3);
    for (kstr = 0; kstr < 2048; kstr++) {
      combinedResponse_re = combinedResponse[kstr].re * d2[kstr].re -
        combinedResponse[kstr].im * d2[kstr].im;
      combinedResponse_im = combinedResponse[kstr].re * d2[kstr].im +
        combinedResponse[kstr].im * d2[kstr].re;
      combinedResponse[kstr].re = combinedResponse_re * d3[kstr].re -
        combinedResponse_im * d3[kstr].im;
      combinedResponse[kstr].im = combinedResponse_re * d3[kstr].im +
        combinedResponse_im * d3[kstr].re;
    }
    break;

   case 8:
    /*  Dec/Int3 */
    c_us(dec_int3_coeff_data, dec_int3_coeff_size, tmp_data, tmp_size);
    d_freqz_cg(tmp_data, tmp_size, w, Fs, combinedResponse);
    break;

   case 9:
    /*  Dec/Int3,Hb1 */
    d_us(hb1_coeff, dv8);
    e_freqz_cg(dv8, w, Fs, combinedResponse);
    c_us(dec_int3_coeff_data, dec_int3_coeff_size, tmp_data, tmp_size);
    d_freqz_cg(tmp_data, tmp_size, w, Fs, d2);
    for (kstr = 0; kstr < 2048; kstr++) {
      combinedResponse_re = combinedResponse[kstr].re;
      combinedResponse[kstr].re = combinedResponse[kstr].re * d2[kstr].re -
        combinedResponse[kstr].im * d2[kstr].im;
      combinedResponse[kstr].im = combinedResponse_re * d2[kstr].im +
        combinedResponse[kstr].im * d2[kstr].re;
    }
    break;

   case 10:
    /*  Dec/Int3,Hb2 */
    f_us(hb2_coeff, dv9);
    g_freqz_cg(dv9, w, Fs, combinedResponse);
    c_us(dec_int3_coeff_data, dec_int3_coeff_size, tmp_data, tmp_size);
    d_freqz_cg(tmp_data, tmp_size, w, Fs, d2);
    for (kstr = 0; kstr < 2048; kstr++) {
      combinedResponse_re = combinedResponse[kstr].re;
      combinedResponse[kstr].re = combinedResponse[kstr].re * d2[kstr].re -
        combinedResponse[kstr].im * d2[kstr].im;
      combinedResponse[kstr].im = combinedResponse_re * d2[kstr].im +
        combinedResponse[kstr].im * d2[kstr].re;
    }
    break;

   case 11:
    /*  Dec/Int3,Hb2,Hb1 */
    e_us(hb1_coeff, dv10);
    f_freqz_cg(dv10, w, Fs, combinedResponse);
    f_us(hb2_coeff, dv14);
    g_freqz_cg(dv14, w, Fs, d2);
    c_us(dec_int3_coeff_data, dec_int3_coeff_size, tmp_data, tmp_size);
    d_freqz_cg(tmp_data, tmp_size, w, Fs, d3);
    for (kstr = 0; kstr < 2048; kstr++) {
      combinedResponse_re = combinedResponse[kstr].re * d2[kstr].re -
        combinedResponse[kstr].im * d2[kstr].im;
      combinedResponse_im = combinedResponse[kstr].re * d2[kstr].im +
        combinedResponse[kstr].im * d2[kstr].re;
      combinedResponse[kstr].re = combinedResponse_re * d3[kstr].re -
        combinedResponse_im * d3[kstr].im;
      combinedResponse[kstr].im = combinedResponse_re * d3[kstr].im +
        combinedResponse_im * d3[kstr].re;
    }
    break;
  }

  /*  Add filter extra filter to end of cascade */
}

/*
 * FREQZ_CG Frequency response of digital filter with codegen support
 *
 *  This function is based on 'freqz' by The MathWorks Inc.
 * Arguments    : const emxArray_real_T *w
 *                double Fs
 *                emxArray_creal_T *hh
 * Return Type  : void
 */
static void h_freqz_cg(const emxArray_real_T *w, double Fs, emxArray_creal_T *hh)
{
  emxArray_real_T *r6;
  int i24;
  int loop_ub;
  emxArray_real_T *digw;
  emxArray_creal_T *s;
  emxArray_creal_T *y;
  boolean_T b1;
  double re;
  double im;
  emxInit_real_T(&r6, 2);

  /*  Cast to enforce precision rules */
  /*  Remaining are default or for advanced use */
  /*  Make b a row */
  /* -------------------------------------------------------------------------- */
  /*  Actual Frequency Response Computation */
  /* if fvflag, */
  /*    Frequency vector specified.  Use Horner's method of polynomial */
  /*    evaluation at the frequency points and divide the numerator */
  /*    by the denominator. */
  /*  */
  /*    Note: we use positive i here because of the relationship */
  /*             polyval(a,exp(1i*w)) = fft(a).*exp(1i*w*(length(a)-1)) */
  /*                ( assuming w = 2*pi*(0:length(a)-1)/length(a) ) */
  /*  */
  /*  Fs was specified, freq. vector is in Hz */
  i24 = r6->size[0] * r6->size[1];
  r6->size[0] = 1;
  r6->size[1] = w->size[1];
  emxEnsureCapacity((emxArray__common *)r6, i24, (int)sizeof(double));
  loop_ub = w->size[0] * w->size[1];
  for (i24 = 0; i24 < loop_ub; i24++) {
    r6->data[i24] = 6.2831853071795862 * w->data[i24];
  }

  emxInit_real_T(&digw, 2);
  emxInit_creal_T(&s, 2);
  rdivide(r6, Fs, digw);

  /*  Convert from Hz to rad/sample for computational purposes */
  i24 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = digw->size[1];
  emxEnsureCapacity((emxArray__common *)s, i24, (int)sizeof(creal_T));
  loop_ub = digw->size[0] * digw->size[1];
  emxFree_real_T(&r6);
  for (i24 = 0; i24 < loop_ub; i24++) {
    s->data[i24].re = digw->data[i24] * 0.0;
    s->data[i24].im = digw->data[i24];
  }

  emxInit_creal_T(&y, 2);
  c_exp(s);

  /*  Digital frequency must be used for this calculation */
  i24 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = s->size[1];
  emxEnsureCapacity((emxArray__common *)y, i24, (int)sizeof(creal_T));
  b1 = (y->size[1] == 0);
  if (!b1) {
    i24 = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)y, i24, (int)sizeof(creal_T));
    loop_ub = y->size[1];
    for (i24 = 0; i24 < loop_ub; i24++) {
      y->data[y->size[0] * i24].re = 1.0;
      y->data[y->size[0] * i24].im = 0.0;
    }
  }

  i24 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = digw->size[1];
  emxEnsureCapacity((emxArray__common *)s, i24, (int)sizeof(creal_T));
  loop_ub = digw->size[0] * digw->size[1];
  for (i24 = 0; i24 < loop_ub; i24++) {
    re = digw->data[i24] * 0.0;
    im = digw->data[i24];
    s->data[i24].re = 0.0 * re;
    s->data[i24].im = 0.0 * im;
  }

  emxFree_real_T(&digw);
  c_exp(s);
  b_rdivide(y, s, hh);

  /*  Generate the default structure to pass to freqzplot */
  /*  If rad/sample, Fs is empty */
  emxFree_creal_T(&y);
  emxFree_creal_T(&s);
}

/*
 * FREQZ_CG Frequency response of digital filter with codegen support
 *
 *  This function is based on 'freqz' by The MathWorks Inc.
 * Arguments    : const double b[15]
 *                const emxArray_real_T *w
 *                double Fs
 *                emxArray_creal_T *hh
 * Return Type  : void
 */
static void i_freqz_cg(const double b[15], const emxArray_real_T *w, double Fs,
  emxArray_creal_T *hh)
{
  emxArray_real_T *r7;
  int i26;
  int loop_ub;
  emxArray_real_T *digw;
  emxArray_creal_T *s;
  emxArray_creal_T *y;
  boolean_T b2;
  int k;
  double re;
  double im;
  emxInit_real_T(&r7, 2);

  /*  Cast to enforce precision rules */
  /*  Remaining are default or for advanced use */
  /*  Make b a row */
  /* -------------------------------------------------------------------------- */
  /*  Actual Frequency Response Computation */
  /* if fvflag, */
  /*    Frequency vector specified.  Use Horner's method of polynomial */
  /*    evaluation at the frequency points and divide the numerator */
  /*    by the denominator. */
  /*  */
  /*    Note: we use positive i here because of the relationship */
  /*             polyval(a,exp(1i*w)) = fft(a).*exp(1i*w*(length(a)-1)) */
  /*                ( assuming w = 2*pi*(0:length(a)-1)/length(a) ) */
  /*  */
  /*  Fs was specified, freq. vector is in Hz */
  i26 = r7->size[0] * r7->size[1];
  r7->size[0] = 1;
  r7->size[1] = w->size[1];
  emxEnsureCapacity((emxArray__common *)r7, i26, (int)sizeof(double));
  loop_ub = w->size[0] * w->size[1];
  for (i26 = 0; i26 < loop_ub; i26++) {
    r7->data[i26] = 6.2831853071795862 * w->data[i26];
  }

  emxInit_real_T(&digw, 2);
  emxInit_creal_T(&s, 2);
  rdivide(r7, Fs, digw);

  /*  Convert from Hz to rad/sample for computational purposes */
  i26 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = digw->size[1];
  emxEnsureCapacity((emxArray__common *)s, i26, (int)sizeof(creal_T));
  loop_ub = digw->size[0] * digw->size[1];
  emxFree_real_T(&r7);
  for (i26 = 0; i26 < loop_ub; i26++) {
    s->data[i26].re = digw->data[i26] * 0.0;
    s->data[i26].im = digw->data[i26];
  }

  emxInit_creal_T(&y, 2);
  c_exp(s);

  /*  Digital frequency must be used for this calculation */
  i26 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = s->size[1];
  emxEnsureCapacity((emxArray__common *)y, i26, (int)sizeof(creal_T));
  b2 = (y->size[1] == 0);
  if (!b2) {
    i26 = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)y, i26, (int)sizeof(creal_T));
    loop_ub = y->size[1];
    for (i26 = 0; i26 < loop_ub; i26++) {
      y->data[y->size[0] * i26].re = b[0];
      y->data[y->size[0] * i26].im = 0.0;
    }

    for (k = 0; k < 14; k++) {
      i26 = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = s->size[1];
      emxEnsureCapacity((emxArray__common *)y, i26, (int)sizeof(creal_T));
      loop_ub = s->size[0] * s->size[1];
      for (i26 = 0; i26 < loop_ub; i26++) {
        re = s->data[i26].re * y->data[i26].re - s->data[i26].im * y->data[i26].
          im;
        im = s->data[i26].re * y->data[i26].im + s->data[i26].im * y->data[i26].
          re;
        y->data[i26].re = re + b[k + 1];
        y->data[i26].im = im;
      }
    }
  }

  i26 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = digw->size[1];
  emxEnsureCapacity((emxArray__common *)s, i26, (int)sizeof(creal_T));
  loop_ub = digw->size[0] * digw->size[1];
  for (i26 = 0; i26 < loop_ub; i26++) {
    re = digw->data[i26] * 0.0;
    im = digw->data[i26];
    s->data[i26].re = 14.0 * re;
    s->data[i26].im = 14.0 * im;
  }

  emxFree_real_T(&digw);
  c_exp(s);
  b_rdivide(y, s, hh);

  /*  Generate the default structure to pass to freqzplot */
  /*  If rad/sample, Fs is empty */
  emxFree_creal_T(&y);
  emxFree_creal_T(&s);
}

/*
 * FREQZ_CG Frequency response of digital filter with codegen support
 *
 *  This function is based on 'freqz' by The MathWorks Inc.
 * Arguments    : const double b[7]
 *                const emxArray_real_T *w
 *                double Fs
 *                emxArray_creal_T *hh
 * Return Type  : void
 */
static void j_freqz_cg(const double b[7], const emxArray_real_T *w, double Fs,
  emxArray_creal_T *hh)
{
  emxArray_real_T *r8;
  int i27;
  int loop_ub;
  emxArray_real_T *digw;
  emxArray_creal_T *s;
  emxArray_creal_T *y;
  boolean_T b3;
  int k;
  double re;
  double im;
  emxInit_real_T(&r8, 2);

  /*  Cast to enforce precision rules */
  /*  Remaining are default or for advanced use */
  /*  Make b a row */
  /* -------------------------------------------------------------------------- */
  /*  Actual Frequency Response Computation */
  /* if fvflag, */
  /*    Frequency vector specified.  Use Horner's method of polynomial */
  /*    evaluation at the frequency points and divide the numerator */
  /*    by the denominator. */
  /*  */
  /*    Note: we use positive i here because of the relationship */
  /*             polyval(a,exp(1i*w)) = fft(a).*exp(1i*w*(length(a)-1)) */
  /*                ( assuming w = 2*pi*(0:length(a)-1)/length(a) ) */
  /*  */
  /*  Fs was specified, freq. vector is in Hz */
  i27 = r8->size[0] * r8->size[1];
  r8->size[0] = 1;
  r8->size[1] = w->size[1];
  emxEnsureCapacity((emxArray__common *)r8, i27, (int)sizeof(double));
  loop_ub = w->size[0] * w->size[1];
  for (i27 = 0; i27 < loop_ub; i27++) {
    r8->data[i27] = 6.2831853071795862 * w->data[i27];
  }

  emxInit_real_T(&digw, 2);
  emxInit_creal_T(&s, 2);
  rdivide(r8, Fs, digw);

  /*  Convert from Hz to rad/sample for computational purposes */
  i27 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = digw->size[1];
  emxEnsureCapacity((emxArray__common *)s, i27, (int)sizeof(creal_T));
  loop_ub = digw->size[0] * digw->size[1];
  emxFree_real_T(&r8);
  for (i27 = 0; i27 < loop_ub; i27++) {
    s->data[i27].re = digw->data[i27] * 0.0;
    s->data[i27].im = digw->data[i27];
  }

  emxInit_creal_T(&y, 2);
  c_exp(s);

  /*  Digital frequency must be used for this calculation */
  i27 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = s->size[1];
  emxEnsureCapacity((emxArray__common *)y, i27, (int)sizeof(creal_T));
  b3 = (y->size[1] == 0);
  if (!b3) {
    i27 = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)y, i27, (int)sizeof(creal_T));
    loop_ub = y->size[1];
    for (i27 = 0; i27 < loop_ub; i27++) {
      y->data[y->size[0] * i27].re = b[0];
      y->data[y->size[0] * i27].im = 0.0;
    }

    for (k = 0; k < 6; k++) {
      i27 = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = s->size[1];
      emxEnsureCapacity((emxArray__common *)y, i27, (int)sizeof(creal_T));
      loop_ub = s->size[0] * s->size[1];
      for (i27 = 0; i27 < loop_ub; i27++) {
        re = s->data[i27].re * y->data[i27].re - s->data[i27].im * y->data[i27].
          im;
        im = s->data[i27].re * y->data[i27].im + s->data[i27].im * y->data[i27].
          re;
        y->data[i27].re = re + b[k + 1];
        y->data[i27].im = im;
      }
    }
  }

  i27 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = digw->size[1];
  emxEnsureCapacity((emxArray__common *)s, i27, (int)sizeof(creal_T));
  loop_ub = digw->size[0] * digw->size[1];
  for (i27 = 0; i27 < loop_ub; i27++) {
    re = digw->data[i27] * 0.0;
    im = digw->data[i27];
    s->data[i27].re = 6.0 * re;
    s->data[i27].im = 6.0 * im;
  }

  emxFree_real_T(&digw);
  c_exp(s);
  b_rdivide(y, s, hh);

  /*  Generate the default structure to pass to freqzplot */
  /*  If rad/sample, Fs is empty */
  emxFree_creal_T(&y);
  emxFree_creal_T(&s);
}

/*
 * FREQZ_CG Frequency response of digital filter with codegen support
 *
 *  This function is based on 'freqz' by The MathWorks Inc.
 * Arguments    : const double b_data[]
 *                const int b_size[2]
 *                const emxArray_real_T *w
 *                double Fs
 *                emxArray_creal_T *hh
 * Return Type  : void
 */
static void k_freqz_cg(const double b_data[], const int b_size[2], const
  emxArray_real_T *w, double Fs, emxArray_creal_T *hh)
{
  int b_b_size[2];
  int loop_ub;
  int i28;
  emxArray_real_T *r9;
  double b_b_data[1024];
  emxArray_real_T *digw;
  emxArray_creal_T *s;
  emxArray_creal_T *r10;
  double re;
  double im;

  /*  Cast to enforce precision rules */
  /*  Remaining are default or for advanced use */
  /*  Make b a row */
  /* -------------------------------------------------------------------------- */
  b_b_size[0] = 1;
  b_b_size[1] = b_size[1];
  loop_ub = b_size[1];
  for (i28 = 0; i28 < loop_ub; i28++) {
    b_b_data[i28] = b_data[i28];
  }

  emxInit_real_T(&r9, 2);

  /*  Actual Frequency Response Computation */
  /* if fvflag, */
  /*    Frequency vector specified.  Use Horner's method of polynomial */
  /*    evaluation at the frequency points and divide the numerator */
  /*    by the denominator. */
  /*  */
  /*    Note: we use positive i here because of the relationship */
  /*             polyval(a,exp(1i*w)) = fft(a).*exp(1i*w*(length(a)-1)) */
  /*                ( assuming w = 2*pi*(0:length(a)-1)/length(a) ) */
  /*  */
  /*  Fs was specified, freq. vector is in Hz */
  i28 = r9->size[0] * r9->size[1];
  r9->size[0] = 1;
  r9->size[1] = w->size[1];
  emxEnsureCapacity((emxArray__common *)r9, i28, (int)sizeof(double));
  loop_ub = w->size[0] * w->size[1];
  for (i28 = 0; i28 < loop_ub; i28++) {
    r9->data[i28] = 6.2831853071795862 * w->data[i28];
  }

  emxInit_real_T(&digw, 2);
  emxInit_creal_T(&s, 2);
  rdivide(r9, Fs, digw);

  /*  Convert from Hz to rad/sample for computational purposes */
  i28 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = digw->size[1];
  emxEnsureCapacity((emxArray__common *)s, i28, (int)sizeof(creal_T));
  loop_ub = digw->size[0] * digw->size[1];
  emxFree_real_T(&r9);
  for (i28 = 0; i28 < loop_ub; i28++) {
    s->data[i28].re = digw->data[i28] * 0.0;
    s->data[i28].im = digw->data[i28];
  }

  emxInit_creal_T(&r10, 2);
  c_exp(s);

  /*  Digital frequency must be used for this calculation */
  g_polyval(b_b_data, b_b_size, s, r10);
  i28 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = digw->size[1];
  emxEnsureCapacity((emxArray__common *)s, i28, (int)sizeof(creal_T));
  loop_ub = digw->size[0] * digw->size[1];
  for (i28 = 0; i28 < loop_ub; i28++) {
    re = digw->data[i28] * 0.0;
    im = digw->data[i28];
    s->data[i28].re = ((double)b_b_size[1] - 1.0) * re;
    s->data[i28].im = ((double)b_b_size[1] - 1.0) * im;
  }

  emxFree_real_T(&digw);
  c_exp(s);
  b_rdivide(r10, s, hh);

  /*  Generate the default structure to pass to freqzplot */
  /*  If rad/sample, Fs is empty */
  emxFree_creal_T(&r10);
  emxFree_creal_T(&s);
}

/*
 * FREQZ_CG Frequency response of digital filter with codegen support
 *
 *  This function is based on 'freqz' by The MathWorks Inc.
 * Arguments    : const double b[30]
 *                const emxArray_real_T *w
 *                double Fs
 *                emxArray_creal_T *hh
 * Return Type  : void
 */
static void l_freqz_cg(const double b[30], const emxArray_real_T *w, double Fs,
  emxArray_creal_T *hh)
{
  emxArray_real_T *r11;
  int i30;
  int loop_ub;
  emxArray_real_T *digw;
  emxArray_creal_T *s;
  emxArray_creal_T *y;
  boolean_T b5;
  int k;
  double re;
  double im;
  emxInit_real_T(&r11, 2);

  /*  Cast to enforce precision rules */
  /*  Remaining are default or for advanced use */
  /*  Make b a row */
  /* -------------------------------------------------------------------------- */
  /*  Actual Frequency Response Computation */
  /* if fvflag, */
  /*    Frequency vector specified.  Use Horner's method of polynomial */
  /*    evaluation at the frequency points and divide the numerator */
  /*    by the denominator. */
  /*  */
  /*    Note: we use positive i here because of the relationship */
  /*             polyval(a,exp(1i*w)) = fft(a).*exp(1i*w*(length(a)-1)) */
  /*                ( assuming w = 2*pi*(0:length(a)-1)/length(a) ) */
  /*  */
  /*  Fs was specified, freq. vector is in Hz */
  i30 = r11->size[0] * r11->size[1];
  r11->size[0] = 1;
  r11->size[1] = w->size[1];
  emxEnsureCapacity((emxArray__common *)r11, i30, (int)sizeof(double));
  loop_ub = w->size[0] * w->size[1];
  for (i30 = 0; i30 < loop_ub; i30++) {
    r11->data[i30] = 6.2831853071795862 * w->data[i30];
  }

  emxInit_real_T(&digw, 2);
  emxInit_creal_T(&s, 2);
  rdivide(r11, Fs, digw);

  /*  Convert from Hz to rad/sample for computational purposes */
  i30 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = digw->size[1];
  emxEnsureCapacity((emxArray__common *)s, i30, (int)sizeof(creal_T));
  loop_ub = digw->size[0] * digw->size[1];
  emxFree_real_T(&r11);
  for (i30 = 0; i30 < loop_ub; i30++) {
    s->data[i30].re = digw->data[i30] * 0.0;
    s->data[i30].im = digw->data[i30];
  }

  emxInit_creal_T(&y, 2);
  c_exp(s);

  /*  Digital frequency must be used for this calculation */
  i30 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = s->size[1];
  emxEnsureCapacity((emxArray__common *)y, i30, (int)sizeof(creal_T));
  b5 = (y->size[1] == 0);
  if (!b5) {
    i30 = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)y, i30, (int)sizeof(creal_T));
    loop_ub = y->size[1];
    for (i30 = 0; i30 < loop_ub; i30++) {
      y->data[y->size[0] * i30].re = b[0];
      y->data[y->size[0] * i30].im = 0.0;
    }

    for (k = 0; k < 29; k++) {
      i30 = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = s->size[1];
      emxEnsureCapacity((emxArray__common *)y, i30, (int)sizeof(creal_T));
      loop_ub = s->size[0] * s->size[1];
      for (i30 = 0; i30 < loop_ub; i30++) {
        re = s->data[i30].re * y->data[i30].re - s->data[i30].im * y->data[i30].
          im;
        im = s->data[i30].re * y->data[i30].im + s->data[i30].im * y->data[i30].
          re;
        y->data[i30].re = re + b[k + 1];
        y->data[i30].im = im;
      }
    }
  }

  i30 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = digw->size[1];
  emxEnsureCapacity((emxArray__common *)s, i30, (int)sizeof(creal_T));
  loop_ub = digw->size[0] * digw->size[1];
  for (i30 = 0; i30 < loop_ub; i30++) {
    re = digw->data[i30] * 0.0;
    im = digw->data[i30];
    s->data[i30].re = 29.0 * re;
    s->data[i30].im = 29.0 * im;
  }

  emxFree_real_T(&digw);
  c_exp(s);
  b_rdivide(y, s, hh);

  /*  Generate the default structure to pass to freqzplot */
  /*  If rad/sample, Fs is empty */
  emxFree_creal_T(&y);
  emxFree_creal_T(&s);
}

/*
 * LP2LP Lowpass to lowpass analog filter transformation. Codegen support
 *
 *  This function is based on 'lp2lp' by The MathWorks Inc.
 * Arguments    : const emxArray_creal_T *a
 *                const emxArray_real_T *b
 *                double d
 *                double wo
 *                emxArray_creal_T *at
 *                emxArray_real_T *bt
 *                double *dt
 * Return Type  : void
 */
static void lp2lp_cg(const emxArray_creal_T *a, const emxArray_real_T *b, double
                     d, double wo, emxArray_creal_T *at, emxArray_real_T *bt,
                     double *dt)
{
  int i50;
  int loop_ub;

  /*  Transform lowpass to lowpass */
  i50 = at->size[0] * at->size[1];
  at->size[0] = a->size[0];
  at->size[1] = a->size[1];
  emxEnsureCapacity((emxArray__common *)at, i50, (int)sizeof(creal_T));
  loop_ub = a->size[0] * a->size[1];
  for (i50 = 0; i50 < loop_ub; i50++) {
    at->data[i50].re = wo * a->data[i50].re;
    at->data[i50].im = wo * a->data[i50].im;
  }

  i50 = bt->size[0];
  bt->size[0] = b->size[0];
  emxEnsureCapacity((emxArray__common *)bt, i50, (int)sizeof(double));
  loop_ub = b->size[0];
  for (i50 = 0; i50 < loop_ub; i50++) {
    bt->data[i50] = wo * b->data[i50];
  }

  *dt = d;
}

/*
 * FREQZ_CG Frequency response of digital filter with codegen support
 *
 *  This function is based on 'freqz' by The MathWorks Inc.
 * Arguments    : const double b[60]
 *                const emxArray_real_T *w
 *                double Fs
 *                emxArray_creal_T *hh
 * Return Type  : void
 */
static void m_freqz_cg(const double b[60], const emxArray_real_T *w, double Fs,
  emxArray_creal_T *hh)
{
  emxArray_real_T *r12;
  int i31;
  int loop_ub;
  emxArray_real_T *digw;
  emxArray_creal_T *s;
  emxArray_creal_T *y;
  boolean_T b6;
  int k;
  double re;
  double im;
  emxInit_real_T(&r12, 2);

  /*  Cast to enforce precision rules */
  /*  Remaining are default or for advanced use */
  /*  Make b a row */
  /* -------------------------------------------------------------------------- */
  /*  Actual Frequency Response Computation */
  /* if fvflag, */
  /*    Frequency vector specified.  Use Horner's method of polynomial */
  /*    evaluation at the frequency points and divide the numerator */
  /*    by the denominator. */
  /*  */
  /*    Note: we use positive i here because of the relationship */
  /*             polyval(a,exp(1i*w)) = fft(a).*exp(1i*w*(length(a)-1)) */
  /*                ( assuming w = 2*pi*(0:length(a)-1)/length(a) ) */
  /*  */
  /*  Fs was specified, freq. vector is in Hz */
  i31 = r12->size[0] * r12->size[1];
  r12->size[0] = 1;
  r12->size[1] = w->size[1];
  emxEnsureCapacity((emxArray__common *)r12, i31, (int)sizeof(double));
  loop_ub = w->size[0] * w->size[1];
  for (i31 = 0; i31 < loop_ub; i31++) {
    r12->data[i31] = 6.2831853071795862 * w->data[i31];
  }

  emxInit_real_T(&digw, 2);
  emxInit_creal_T(&s, 2);
  rdivide(r12, Fs, digw);

  /*  Convert from Hz to rad/sample for computational purposes */
  i31 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = digw->size[1];
  emxEnsureCapacity((emxArray__common *)s, i31, (int)sizeof(creal_T));
  loop_ub = digw->size[0] * digw->size[1];
  emxFree_real_T(&r12);
  for (i31 = 0; i31 < loop_ub; i31++) {
    s->data[i31].re = digw->data[i31] * 0.0;
    s->data[i31].im = digw->data[i31];
  }

  emxInit_creal_T(&y, 2);
  c_exp(s);

  /*  Digital frequency must be used for this calculation */
  i31 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = s->size[1];
  emxEnsureCapacity((emxArray__common *)y, i31, (int)sizeof(creal_T));
  b6 = (y->size[1] == 0);
  if (!b6) {
    i31 = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)y, i31, (int)sizeof(creal_T));
    loop_ub = y->size[1];
    for (i31 = 0; i31 < loop_ub; i31++) {
      y->data[y->size[0] * i31].re = b[0];
      y->data[y->size[0] * i31].im = 0.0;
    }

    for (k = 0; k < 59; k++) {
      i31 = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = s->size[1];
      emxEnsureCapacity((emxArray__common *)y, i31, (int)sizeof(creal_T));
      loop_ub = s->size[0] * s->size[1];
      for (i31 = 0; i31 < loop_ub; i31++) {
        re = s->data[i31].re * y->data[i31].re - s->data[i31].im * y->data[i31].
          im;
        im = s->data[i31].re * y->data[i31].im + s->data[i31].im * y->data[i31].
          re;
        y->data[i31].re = re + b[k + 1];
        y->data[i31].im = im;
      }
    }
  }

  i31 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = digw->size[1];
  emxEnsureCapacity((emxArray__common *)s, i31, (int)sizeof(creal_T));
  loop_ub = digw->size[0] * digw->size[1];
  for (i31 = 0; i31 < loop_ub; i31++) {
    re = digw->data[i31] * 0.0;
    im = digw->data[i31];
    s->data[i31].re = 59.0 * re;
    s->data[i31].im = 59.0 * im;
  }

  emxFree_real_T(&digw);
  c_exp(s);
  b_rdivide(y, s, hh);

  /*  Generate the default structure to pass to freqzplot */
  /*  If rad/sample, Fs is empty */
  emxFree_creal_T(&y);
  emxFree_creal_T(&s);
}

/*
 * Arguments    : double y
 * Return Type  : double
 */
static double mag2db(double y)
{
  return 20.0 * log10(y);
}

/*
 * Arguments    : const double A[4]
 *                const double B[4]
 *                double Y[4]
 * Return Type  : void
 */
static void mldivide(const double A[4], const double B[4], double Y[4])
{
  int r1;
  int r2;
  double a21;
  double a22;
  int k;
  if (fabs(A[1]) > fabs(A[0])) {
    r1 = 1;
    r2 = 0;
  } else {
    r1 = 0;
    r2 = 1;
  }

  a21 = A[r2] / A[r1];
  a22 = A[2 + r2] - a21 * A[2 + r1];
  for (k = 0; k < 2; k++) {
    Y[1 + (k << 1)] = (B[r2 + (k << 1)] - B[r1 + (k << 1)] * a21) / a22;
    Y[k << 1] = (B[r1 + (k << 1)] - Y[1 + (k << 1)] * A[2 + r1]) / A[r1];
  }
}

/*
 * FREQZ_CG Frequency response of digital filter with codegen support
 *
 *  This function is based on 'freqz' by The MathWorks Inc.
 * Arguments    : const double b[14]
 *                const emxArray_real_T *w
 *                double Fs
 *                emxArray_creal_T *hh
 * Return Type  : void
 */
static void n_freqz_cg(const double b[14], const emxArray_real_T *w, double Fs,
  emxArray_creal_T *hh)
{
  emxArray_real_T *r13;
  int i32;
  int loop_ub;
  emxArray_real_T *digw;
  emxArray_creal_T *s;
  emxArray_creal_T *y;
  boolean_T b7;
  int k;
  double re;
  double im;
  emxInit_real_T(&r13, 2);

  /*  Cast to enforce precision rules */
  /*  Remaining are default or for advanced use */
  /*  Make b a row */
  /* -------------------------------------------------------------------------- */
  /*  Actual Frequency Response Computation */
  /* if fvflag, */
  /*    Frequency vector specified.  Use Horner's method of polynomial */
  /*    evaluation at the frequency points and divide the numerator */
  /*    by the denominator. */
  /*  */
  /*    Note: we use positive i here because of the relationship */
  /*             polyval(a,exp(1i*w)) = fft(a).*exp(1i*w*(length(a)-1)) */
  /*                ( assuming w = 2*pi*(0:length(a)-1)/length(a) ) */
  /*  */
  /*  Fs was specified, freq. vector is in Hz */
  i32 = r13->size[0] * r13->size[1];
  r13->size[0] = 1;
  r13->size[1] = w->size[1];
  emxEnsureCapacity((emxArray__common *)r13, i32, (int)sizeof(double));
  loop_ub = w->size[0] * w->size[1];
  for (i32 = 0; i32 < loop_ub; i32++) {
    r13->data[i32] = 6.2831853071795862 * w->data[i32];
  }

  emxInit_real_T(&digw, 2);
  emxInit_creal_T(&s, 2);
  rdivide(r13, Fs, digw);

  /*  Convert from Hz to rad/sample for computational purposes */
  i32 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = digw->size[1];
  emxEnsureCapacity((emxArray__common *)s, i32, (int)sizeof(creal_T));
  loop_ub = digw->size[0] * digw->size[1];
  emxFree_real_T(&r13);
  for (i32 = 0; i32 < loop_ub; i32++) {
    s->data[i32].re = digw->data[i32] * 0.0;
    s->data[i32].im = digw->data[i32];
  }

  emxInit_creal_T(&y, 2);
  c_exp(s);

  /*  Digital frequency must be used for this calculation */
  i32 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = s->size[1];
  emxEnsureCapacity((emxArray__common *)y, i32, (int)sizeof(creal_T));
  b7 = (y->size[1] == 0);
  if (!b7) {
    i32 = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)y, i32, (int)sizeof(creal_T));
    loop_ub = y->size[1];
    for (i32 = 0; i32 < loop_ub; i32++) {
      y->data[y->size[0] * i32].re = b[0];
      y->data[y->size[0] * i32].im = 0.0;
    }

    for (k = 0; k < 13; k++) {
      i32 = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = s->size[1];
      emxEnsureCapacity((emxArray__common *)y, i32, (int)sizeof(creal_T));
      loop_ub = s->size[0] * s->size[1];
      for (i32 = 0; i32 < loop_ub; i32++) {
        re = s->data[i32].re * y->data[i32].re - s->data[i32].im * y->data[i32].
          im;
        im = s->data[i32].re * y->data[i32].im + s->data[i32].im * y->data[i32].
          re;
        y->data[i32].re = re + b[k + 1];
        y->data[i32].im = im;
      }
    }
  }

  i32 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = digw->size[1];
  emxEnsureCapacity((emxArray__common *)s, i32, (int)sizeof(creal_T));
  loop_ub = digw->size[0] * digw->size[1];
  for (i32 = 0; i32 < loop_ub; i32++) {
    re = digw->data[i32] * 0.0;
    im = digw->data[i32];
    s->data[i32].re = 13.0 * re;
    s->data[i32].im = 13.0 * im;
  }

  emxFree_real_T(&digw);
  c_exp(s);
  b_rdivide(y, s, hh);

  /*  Generate the default structure to pass to freqzplot */
  /*  If rad/sample, Fs is empty */
  emxFree_creal_T(&y);
  emxFree_creal_T(&s);
}

/*
 * Arguments    : const emxArray_creal_T *x
 *                emxArray_creal_T *c
 * Return Type  : void
 */
static void poly(const emxArray_creal_T *x, emxArray_creal_T *c)
{
  emxArray_creal_T *alpha1;
  int info;
  boolean_T p;
  int i;
  boolean_T exitg2;
  emxArray_creal_T *beta1;
  int exitg1;
  double x_re;
  double alpha1_re;
  double x_im;
  double alpha1_im;
  emxArray_creal_T *h;
  boolean_T b_x;
  double beta1_re;
  double beta1_im;
  double brm;
  int istart;
  int jend;
  emxInit_creal_T1(&alpha1, 1);
  if ((x->size[0] == 0) || (x->size[1] == 0)) {
    info = alpha1->size[0];
    alpha1->size[0] = x->size[0];
    emxEnsureCapacity((emxArray__common *)alpha1, info, (int)sizeof(creal_T));
    i = x->size[0];
    for (info = 0; info < i; info++) {
      alpha1->data[info].re = 0.0;
      alpha1->data[info].im = 0.0;
    }
  } else if ((x->size[0] == 1) && (x->size[1] == 1)) {
    info = alpha1->size[0];
    alpha1->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)alpha1, info, (int)sizeof(creal_T));
    alpha1->data[0] = x->data[0];
  } else {
    p = (x->size[0] == x->size[1]);
    if (p) {
      info = 0;
      exitg2 = false;
      while ((!exitg2) && (info <= x->size[1] - 1)) {
        i = 0;
        do {
          exitg1 = 0;
          if (i <= info) {
            x_re = x->data[info + x->size[0] * i].re;
            x_im = -x->data[info + x->size[0] * i].im;
            b_x = ((x->data[i + x->size[0] * info].re == x_re) && (x->data[i +
                    x->size[0] * info].im == x_im));
            if (!b_x) {
              p = false;
              exitg1 = 1;
            } else {
              i++;
            }
          } else {
            info++;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    }

    if (p) {
      emxInit_creal_T(&beta1, 2);
      info = beta1->size[0] * beta1->size[1];
      beta1->size[0] = x->size[0];
      beta1->size[1] = x->size[1];
      emxEnsureCapacity((emxArray__common *)beta1, info, (int)sizeof(creal_T));
      i = x->size[0] * x->size[1];
      for (info = 0; info < i; info++) {
        beta1->data[info] = x->data[info];
      }

      emxInit_creal_T(&h, 2);
      xgehrd(beta1);
      eml_zlahqr(beta1);
      info = h->size[0] * h->size[1];
      h->size[0] = beta1->size[0];
      h->size[1] = beta1->size[1];
      emxEnsureCapacity((emxArray__common *)h, info, (int)sizeof(creal_T));
      i = beta1->size[0] * beta1->size[1];
      for (info = 0; info < i; info++) {
        h->data[info] = beta1->data[info];
      }

      if ((beta1->size[0] == 0) || (beta1->size[1] == 0) || (3 >= beta1->size[0]))
      {
      } else {
        istart = 4;
        if (beta1->size[0] - 4 < beta1->size[1] - 1) {
          jend = beta1->size[0] - 3;
        } else {
          jend = beta1->size[1];
        }

        for (info = 1; info <= jend; info++) {
          for (i = istart; i <= beta1->size[0]; i++) {
            h->data[(i + h->size[0] * (info - 1)) - 1].re = 0.0;
            h->data[(i + h->size[0] * (info - 1)) - 1].im = 0.0;
          }

          istart++;
        }
      }

      emxFree_creal_T(&beta1);
      info = alpha1->size[0];
      alpha1->size[0] = h->size[0];
      emxEnsureCapacity((emxArray__common *)alpha1, info, (int)sizeof(creal_T));
      for (info = 0; info + 1 <= h->size[0]; info++) {
        alpha1->data[info] = h->data[info + h->size[0] * info];
      }

      emxFree_creal_T(&h);
    } else {
      emxInit_creal_T1(&beta1, 1);
      xzgeev(x, &info, alpha1, beta1);
      info = alpha1->size[0];
      emxEnsureCapacity((emxArray__common *)alpha1, info, (int)sizeof(creal_T));
      i = alpha1->size[0];
      for (info = 0; info < i; info++) {
        alpha1_re = alpha1->data[info].re;
        alpha1_im = alpha1->data[info].im;
        beta1_re = beta1->data[info].re;
        beta1_im = beta1->data[info].im;
        if (beta1_im == 0.0) {
          if (alpha1_im == 0.0) {
            alpha1->data[info].re = alpha1_re / beta1_re;
            alpha1->data[info].im = 0.0;
          } else if (alpha1_re == 0.0) {
            alpha1->data[info].re = 0.0;
            alpha1->data[info].im = alpha1_im / beta1_re;
          } else {
            alpha1->data[info].re = alpha1_re / beta1_re;
            alpha1->data[info].im = alpha1_im / beta1_re;
          }
        } else if (beta1_re == 0.0) {
          if (alpha1_re == 0.0) {
            alpha1->data[info].re = alpha1_im / beta1_im;
            alpha1->data[info].im = 0.0;
          } else if (alpha1_im == 0.0) {
            alpha1->data[info].re = 0.0;
            alpha1->data[info].im = -(alpha1_re / beta1_im);
          } else {
            alpha1->data[info].re = alpha1_im / beta1_im;
            alpha1->data[info].im = -(alpha1_re / beta1_im);
          }
        } else {
          brm = fabs(beta1_re);
          x_re = fabs(beta1_im);
          if (brm > x_re) {
            x_im = beta1_im / beta1_re;
            x_re = beta1_re + x_im * beta1_im;
            alpha1->data[info].re = (alpha1_re + x_im * alpha1_im) / x_re;
            alpha1->data[info].im = (alpha1_im - x_im * alpha1_re) / x_re;
          } else if (x_re == brm) {
            if (beta1_re > 0.0) {
              x_im = 0.5;
            } else {
              x_im = -0.5;
            }

            if (beta1_im > 0.0) {
              x_re = 0.5;
            } else {
              x_re = -0.5;
            }

            alpha1->data[info].re = (alpha1_re * x_im + alpha1_im * x_re) / brm;
            alpha1->data[info].im = (alpha1_im * x_im - alpha1_re * x_re) / brm;
          } else {
            x_im = beta1_re / beta1_im;
            x_re = beta1_im + x_im * beta1_re;
            alpha1->data[info].re = (x_im * alpha1_re + alpha1_im) / x_re;
            alpha1->data[info].im = (x_im * alpha1_im - alpha1_re) / x_re;
          }
        }
      }

      emxFree_creal_T(&beta1);
    }
  }

  vector_poly(alpha1, c);
  emxFree_creal_T(&alpha1);
}

/*
 * Arguments    : const double p[15]
 *                const creal_T x[2048]
 *                creal_T y[2048]
 * Return Type  : void
 */
static void polyval(const double p[15], const creal_T x[2048], creal_T y[2048])
{
  int i8;
  int k;
  double x_im;
  for (i8 = 0; i8 < 2048; i8++) {
    y[i8].re = p[0];
    y[i8].im = 0.0;
  }

  for (k = 0; k < 14; k++) {
    for (i8 = 0; i8 < 2048; i8++) {
      x_im = x[i8].re * y[i8].im + x[i8].im * y[i8].re;
      y[i8].re = (x[i8].re * y[i8].re - x[i8].im * y[i8].im) + p[k + 1];
      y[i8].im = x_im;
    }
  }
}

/*
 * Arguments    : const double a[2048]
 *                double y[2048]
 * Return Type  : void
 */
static void power(const double a[2048], double y[2048])
{
  int k;
  for (k = 0; k < 2048; k++) {
    y[k] = a[k] * a[k];
  }
}

/*
 * Arguments    : const emxArray_real_T *x
 *                double y
 *                emxArray_real_T *z
 * Return Type  : void
 */
static void rdivide(const emxArray_real_T *x, double y, emxArray_real_T *z)
{
  int i23;
  int loop_ub;
  i23 = z->size[0] * z->size[1];
  z->size[0] = 1;
  z->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)z, i23, (int)sizeof(double));
  loop_ub = x->size[0] * x->size[1];
  for (i23 = 0; i23 < loop_ub; i23++) {
    z->data[i23] = x->data[i23] / y;
  }
}

/*
 * Arguments    : const creal_T y
 * Return Type  : creal_T
 */
static creal_T recip(const creal_T y)
{
  creal_T z;
  double brm;
  double bim;
  double d;
  brm = fabs(y.re);
  bim = fabs(y.im);
  if (y.im == 0.0) {
    z.re = 1.0 / y.re;
    z.im = 0.0;
  } else if (y.re == 0.0) {
    z.re = 0.0;
    z.im = -1.0 / y.im;
  } else if (brm > bim) {
    bim = y.im / y.re;
    d = y.re + bim * y.im;
    z.re = 1.0 / d;
    z.im = -bim / d;
  } else if (brm == bim) {
    bim = 0.5;
    if (y.re < 0.0) {
      bim = -0.5;
    }

    d = 0.5;
    if (y.im < 0.0) {
      d = -0.5;
    }

    z.re = bim / brm;
    z.im = -d / brm;
  } else {
    bim = y.re / y.im;
    d = y.im + bim * y.re;
    z.re = bim / d;
    z.im = -1.0 / d;
  }

  return z;
}

/*
 * REMEZDD Lagrange interpolation coefficients.
 * Arguments    : double k
 *                double n
 *                double m
 *                const emxArray_real_T *x
 * Return Type  : double
 */
static double remezdd(double k, double n, double m, const emxArray_real_T *x)
{
  double y;
  int l;
  emxArray_real_T *xx;
  emxArray_int32_T *r25;
  int i45;
  int i;
  int end;
  double b_x;
  int loop_ub;

  /*  */
  /*    Author: T. Krauss 1993 */
  /*        Was Revision: 1.4, Date: 1994/01/25 17:59:44 */
  y = 1.0;
  l = 0;
  emxInit_real_T(&xx, 2);
  emxInit_int32_T(&r25, 2);
  while (l <= (int)m - 1) {
    if ((m == 0.0) || (((m > 0.0) && (1.0 + (double)l > n)) || ((0.0 > m) && (n >
           1.0 + (double)l)))) {
      i45 = 1;
      i = 1;
      end = 0;
    } else {
      i45 = l + 1;
      i = (int)m;
      end = (int)n;
    }

    b_x = x->data[(int)k - 1];
    loop_ub = xx->size[0] * xx->size[1];
    xx->size[0] = 1;
    xx->size[1] = div_s32_floor(end - i45, i) + 1;
    emxEnsureCapacity((emxArray__common *)xx, loop_ub, (int)sizeof(double));
    loop_ub = div_s32_floor(end - i45, i);
    for (end = 0; end <= loop_ub; end++) {
      xx->data[xx->size[0] * end] = 2.0 * (b_x - x->data[(i45 + i * end) - 1]);
    }

    end = xx->size[1] - 1;
    loop_ub = 0;
    for (i = 0; i <= end; i++) {
      if (xx->data[i] != 0.0) {
        loop_ub++;
      }
    }

    i45 = r25->size[0] * r25->size[1];
    r25->size[0] = 1;
    r25->size[1] = loop_ub;
    emxEnsureCapacity((emxArray__common *)r25, i45, (int)sizeof(int));
    loop_ub = 0;
    for (i = 0; i <= end; i++) {
      if (xx->data[i] != 0.0) {
        r25->data[loop_ub] = i + 1;
        loop_ub++;
      }
    }

    if (r25->size[1] == 0) {
      b_x = 1.0;
    } else {
      b_x = xx->data[r25->data[0] - 1];
      for (loop_ub = 2; loop_ub <= r25->size[1]; loop_ub++) {
        b_x *= xx->data[r25->data[r25->size[0] * (loop_ub - 1)] - 1];
      }
    }

    y *= b_x;
    l++;
  }

  emxFree_int32_T(&r25);
  emxFree_real_T(&xx);
  y = 1.0 / y;

  /*  EOF */
  return y;
}

/*
 * remezm function
 *  Inputs
 *      nfilt - filter length
 *      edge - vector of band edges (between 0 and .5)
 *      grid - frequency grid (between 0 and .5)
 *      des - desired function on frequency grid
 *      wt - weight function on frequency grid
 *      neg == 1 ==> antisymmetric imp resp,
 *          == 0 ==> symmetric imp resp
 *  Outputs
 *      h - coefficients of basis functions
 *      dev - computed error
 *      iext - indices of extremal frequencies
 * Arguments    : double nfilt
 *                const double edge[4]
 *                const emxArray_real_T *grid
 *                emxArray_real_T *des
 *                emxArray_real_T *wt
 *                emxArray_real_T *h
 *                double *dev
 *                boolean_T *valid
 * Return Type  : void
 */
static void remezm(double nfilt, const double edge[4], const emxArray_real_T
                   *grid, emxArray_real_T *des, emxArray_real_T *wt,
                   emxArray_real_T *h, double *dev, boolean_T *valid)
{
  double nodd;
  double nfcns;
  int varargin_2;
  int ngrid;
  emxArray_real_T *j;
  emxArray_real_T *x2;
  int i39;
  double temp;
  int loop_ub;
  emxArray_real_T *iext;
  emxArray_real_T *r22;
  emxArray_real_T *y;
  emxArray_real_T *x;
  double comp;
  double dtemp;
  double b_y1;
  int luck;
  int nut1;
  double err;
  int ixstart;
  int flag;
  double devl;
  emxArray_real_T *ad;
  int niter;
  double jchnge;
  double d1;
  emxArray_real_T *l;
  emxArray_real_T *a;
  emxArray_int32_T *r23;
  emxArray_real_T *b;
  emxArray_real_T *b_x;
  emxArray_real_T *b_a;
  emxArray_real_T *b_wt;
  emxArray_real_T *c_wt;
  emxArray_real_T *c_x;
  emxArray_real_T *d_x;
  emxArray_real_T *e_x;
  emxArray_real_T *f_x;
  emxArray_real_T *g_x;
  emxArray_real_T *h_x;
  emxArray_real_T *i_x;
  emxArray_real_T *j_x;
  emxArray_real_T *k_x;
  emxArray_real_T *l_x;
  emxArray_real_T *m_x;
  emxArray_real_T *n_x;
  emxArray_int8_T *r24;
  emxArray_int32_T *b_iext;
  emxArray_real_T *mtmp;
  emxArray_real_T *b_mtmp;
  emxArray_real_T *c_iext;
  emxArray_real_T *d_iext;
  emxArray_real_T *k1;
  emxArray_real_T *b_k1;
  emxArray_real_T *e_iext;
  emxArray_real_T *f_iext;
  emxArray_real_T *b_j;
  emxArray_real_T *c_j;
  emxArray_real_T *d_j;
  emxArray_real_T *e_j;
  boolean_T guard1 = false;
  int exitg1;
  double f_j;
  double c_k1;
  int nut;
  double knz;
  double b_l;
  int nu;
  double varargin_1[2];
  double dnum;
  int i40;
  int i41;
  int i42;
  int i43;
  int i44;
  boolean_T guard2 = false;
  int flag34;
  int exitg5;
  boolean_T exitg11;
  boolean_T exitg4;
  boolean_T exitg15;
  boolean_T exitg10;
  boolean_T exitg2;
  boolean_T exitg12;
  boolean_T exitg14;
  boolean_T exitg8;
  boolean_T exitg13;
  boolean_T exitg6;
  boolean_T exitg16;
  boolean_T exitg3;
  boolean_T exitg9;
  boolean_T exitg7;

  /*  */
  *valid = true;
  nodd = rt_remd_snf(nfilt, 2.0);

  /*  nodd == 1 ==> filter length is odd */
  /*  nodd == 0 ==> filter length is even */
  nfcns = nfilt / 2.0;
  b_fix(&nfcns);
  if (nodd == 1.0) {
    nfcns++;
  }

  varargin_2 = grid->size[1];
  ngrid = grid->size[1];
  emxInit_real_T(&j, 2);
  emxInit_real_T(&x2, 2);
  if (nodd != 1.0) {
    i39 = x2->size[0] * x2->size[1];
    x2->size[0] = 1;
    x2->size[1] = grid->size[1];
    emxEnsureCapacity((emxArray__common *)x2, i39, (int)sizeof(double));
    loop_ub = grid->size[0] * grid->size[1];
    for (i39 = 0; i39 < loop_ub; i39++) {
      x2->data[i39] = 3.1415926535897931 * grid->data[i39];
    }

    b_cos(x2);
    c_rdivide(des, x2, j);
    i39 = des->size[0] * des->size[1];
    des->size[0] = 1;
    des->size[1] = j->size[1];
    emxEnsureCapacity((emxArray__common *)des, i39, (int)sizeof(double));
    loop_ub = j->size[0] * j->size[1];
    for (i39 = 0; i39 < loop_ub; i39++) {
      des->data[i39] = j->data[i39];
    }

    i39 = x2->size[0] * x2->size[1];
    x2->size[0] = 1;
    x2->size[1] = grid->size[1];
    emxEnsureCapacity((emxArray__common *)x2, i39, (int)sizeof(double));
    loop_ub = grid->size[0] * grid->size[1];
    for (i39 = 0; i39 < loop_ub; i39++) {
      x2->data[i39] = 3.1415926535897931 * grid->data[i39];
    }

    emxInit_real_T(&r22, 2);
    b_cos(x2);
    i39 = r22->size[0] * r22->size[1];
    r22->size[0] = 1;
    r22->size[1] = x2->size[1];
    emxEnsureCapacity((emxArray__common *)r22, i39, (int)sizeof(double));
    loop_ub = x2->size[0] * x2->size[1];
    for (i39 = 0; i39 < loop_ub; i39++) {
      r22->data[i39] = x2->data[i39];
    }

    i39 = wt->size[0] * wt->size[1];
    wt->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)wt, i39, (int)sizeof(double));
    ixstart = wt->size[0];
    flag = wt->size[1];
    loop_ub = ixstart * flag;
    for (i39 = 0; i39 < loop_ub; i39++) {
      wt->data[i39] *= r22->data[i39];
    }

    emxFree_real_T(&r22);
  }

  temp = ((double)grid->size[1] - 1.0) / nfcns;
  if (rtIsNaN(nfcns)) {
    i39 = j->size[0] * j->size[1];
    j->size[0] = 1;
    j->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)j, i39, (int)sizeof(double));
    j->data[0] = rtNaN;
  } else if (nfcns < 1.0) {
    i39 = j->size[0] * j->size[1];
    j->size[0] = 1;
    j->size[1] = 0;
    emxEnsureCapacity((emxArray__common *)j, i39, (int)sizeof(double));
  } else if (rtIsInf(nfcns) && (1.0 == nfcns)) {
    i39 = j->size[0] * j->size[1];
    j->size[0] = 1;
    j->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)j, i39, (int)sizeof(double));
    j->data[0] = rtNaN;
  } else {
    i39 = j->size[0] * j->size[1];
    j->size[0] = 1;
    j->size[1] = (int)floor(nfcns - 1.0) + 1;
    emxEnsureCapacity((emxArray__common *)j, i39, (int)sizeof(double));
    loop_ub = (int)floor(nfcns - 1.0);
    for (i39 = 0; i39 <= loop_ub; i39++) {
      j->data[j->size[0] * i39] = 1.0 + (double)i39;
    }
  }

  i39 = x2->size[0] * x2->size[1];
  x2->size[0] = 1;
  x2->size[1] = j->size[1] + 1;
  emxEnsureCapacity((emxArray__common *)x2, i39, (int)sizeof(double));
  loop_ub = j->size[1];
  for (i39 = 0; i39 < loop_ub; i39++) {
    x2->data[x2->size[0] * i39] = temp * (j->data[j->size[0] * i39] - 1.0) + 1.0;
  }

  emxInit_real_T1(&iext, 1);
  x2->data[x2->size[0] * j->size[1]] = grid->size[1];
  c_fix(x2);
  i39 = iext->size[0];
  iext->size[0] = x2->size[1] + 1;
  emxEnsureCapacity((emxArray__common *)iext, i39, (int)sizeof(double));
  loop_ub = x2->size[1];
  for (i39 = 0; i39 < loop_ub; i39++) {
    iext->data[i39] = x2->data[x2->size[0] * i39];
  }

  emxInit_real_T(&y, 2);
  emxInit_real_T(&x, 2);
  iext->data[x2->size[1]] = 0.0;

  /*  Remez exchange loop */
  comp = -1.0;
  dtemp = -1.0;
  b_y1 = -1.0;
  luck = -1;
  nut1 = -1;
  err = -1.0;
  i39 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = 1;
  emxEnsureCapacity((emxArray__common *)y, i39, (int)sizeof(double));
  y->data[0] = -1.0;
  *dev = -1.0;
  devl = -1.0;
  i39 = x->size[0] * x->size[1];
  x->size[0] = 1;
  x->size[1] = (int)(nfcns + 1.0);
  emxEnsureCapacity((emxArray__common *)x, i39, (int)sizeof(double));
  loop_ub = (int)(nfcns + 1.0);
  for (i39 = 0; i39 < loop_ub; i39++) {
    x->data[i39] = 0.0;
  }

  emxInit_real_T(&ad, 2);
  niter = 0;
  jchnge = 1.0;
  d1 = (nfcns - 1.0) / 15.0;
  b_fix(&d1);
  i39 = ad->size[0] * ad->size[1];
  ad->size[0] = 1;
  ad->size[1] = (int)(nfcns + 1.0);
  emxEnsureCapacity((emxArray__common *)ad, i39, (int)sizeof(double));
  loop_ub = (int)(nfcns + 1.0);
  for (i39 = 0; i39 < loop_ub; i39++) {
    ad->data[i39] = 0.0;
  }

  /*  index manager(s) */
  emxInit_real_T(&l, 2);
  emxInit_real_T(&a, 2);
  emxInit_int32_T(&r23, 2);
  emxInit_real_T1(&b, 1);
  emxInit_real_T(&b_x, 2);
  emxInit_real_T(&b_a, 2);
  emxInit_real_T(&b_wt, 2);
  emxInit_real_T(&c_wt, 2);
  emxInit_real_T(&c_x, 2);
  emxInit_real_T(&d_x, 2);
  emxInit_real_T(&e_x, 2);
  emxInit_real_T(&f_x, 2);
  emxInit_real_T(&g_x, 2);
  emxInit_real_T(&h_x, 2);
  emxInit_real_T(&i_x, 2);
  emxInit_real_T(&j_x, 2);
  emxInit_real_T(&k_x, 2);
  emxInit_real_T(&l_x, 2);
  emxInit_real_T(&m_x, 2);
  emxInit_real_T(&n_x, 2);
  emxInit_int8_T(&r24, 2);
  emxInit_int32_T1(&b_iext, 1);
  emxInit_real_T(&mtmp, 2);
  emxInit_real_T1(&b_mtmp, 1);
  emxInit_real_T(&c_iext, 2);
  emxInit_real_T1(&d_iext, 1);
  emxInit_real_T(&k1, 2);
  emxInit_real_T1(&b_k1, 1);
  emxInit_real_T(&e_iext, 2);
  emxInit_real_T1(&f_iext, 1);
  emxInit_real_T1(&b_j, 1);
  emxInit_real_T(&c_j, 2);
  emxInit_real_T1(&d_j, 1);
  emxInit_real_T1(&e_j, 1);
  guard1 = false;
  do {
    exitg1 = 0;
    if (jchnge > 0.0) {
      iext->data[(int)((nfcns + 1.0) + 1.0) - 1] = (double)varargin_2 + 1.0;
      niter++;
      if (niter > 250) {
        guard1 = true;
        exitg1 = 1;
      } else {
        if (1.0 > nfcns + 1.0) {
          loop_ub = 0;
        } else {
          loop_ub = (int)(nfcns + 1.0);
        }

        i39 = l->size[0] * l->size[1];
        l->size[0] = 1;
        l->size[1] = loop_ub;
        emxEnsureCapacity((emxArray__common *)l, i39, (int)sizeof(double));
        for (i39 = 0; i39 < loop_ub; i39++) {
          l->data[l->size[0] * i39] = iext->data[i39];
        }

        i39 = x->size[0] * x->size[1];
        x->size[0] = 1;
        x->size[1] = l->size[1];
        emxEnsureCapacity((emxArray__common *)x, i39, (int)sizeof(double));
        nut = l->size[0] * l->size[1];
        for (i39 = 0; i39 < nut; i39++) {
          x->data[i39] = 6.2831853071795862 * grid->data[(int)l->data[i39] - 1];
        }

        b_cos(x);
        for (ixstart = 0; ixstart < (int)(nfcns + 1.0); ixstart++) {
          ad->data[ixstart] = remezdd(1.0 + (double)ixstart, nfcns + 1.0, d1 +
            1.0, x);
        }

        for (i39 = 0; i39 < 2; i39++) {
          varargin_1[i39] = ad->size[i39];
        }

        i39 = j->size[0] * j->size[1];
        j->size[0] = 1;
        j->size[1] = (int)varargin_1[1];
        emxEnsureCapacity((emxArray__common *)j, i39, (int)sizeof(double));
        nut = (int)varargin_1[1];
        for (i39 = 0; i39 < nut; i39++) {
          j->data[i39] = 1.0;
        }

        if (2.0 > nfcns + 1.0) {
          i39 = 0;
          i40 = 1;
          i41 = 0;
          i42 = 0;
          i43 = 1;
        } else {
          i39 = 1;
          i40 = 2;
          i41 = (int)(nfcns + 1.0);
          i42 = 1;
          i43 = 2;
        }

        i44 = r24->size[0] * r24->size[1];
        r24->size[0] = 1;
        r24->size[1] = (int)varargin_1[1];
        emxEnsureCapacity((emxArray__common *)r24, i44, (int)sizeof(signed char));
        nut = (int)varargin_1[1];
        for (i44 = 0; i44 < nut; i44++) {
          r24->data[r24->size[0] * i44] = 1;
        }

        nut = div_s32_floor((i41 - i39) - 1, i40);
        for (i41 = 0; i41 <= nut; i41++) {
          j->data[i42 + i43 * i41] = -(double)r24->data[i39 + i40 * i41];
        }

        i39 = b->size[0];
        b->size[0] = l->size[1];
        emxEnsureCapacity((emxArray__common *)b, i39, (int)sizeof(double));
        nut = l->size[1];
        for (i39 = 0; i39 < nut; i39++) {
          b->data[i39] = des->data[(int)l->data[l->size[0] * i39] - 1];
        }

        guard2 = false;
        if (ad->size[1] == 1) {
          guard2 = true;
        } else {
          i39 = b_iext->size[0];
          b_iext->size[0] = loop_ub;
          emxEnsureCapacity((emxArray__common *)b_iext, i39, (int)sizeof(int));
          for (i39 = 0; i39 < loop_ub; i39++) {
            b_iext->data[i39] = (int)iext->data[i39];
          }

          if (b_iext->size[0] == 1) {
            guard2 = true;
          } else {
            dnum = 0.0;
            for (i39 = 0; i39 < ad->size[1]; i39++) {
              dnum += ad->data[ad->size[0] * i39] * b->data[i39];
            }
          }
        }

        if (guard2) {
          dnum = 0.0;
          for (i39 = 0; i39 < ad->size[1]; i39++) {
            dnum += ad->data[ad->size[0] * i39] * b->data[i39];
          }
        }

        i39 = c_wt->size[0] * c_wt->size[1];
        c_wt->size[0] = 1;
        c_wt->size[1] = l->size[1];
        emxEnsureCapacity((emxArray__common *)c_wt, i39, (int)sizeof(double));
        loop_ub = l->size[0] * l->size[1];
        for (i39 = 0; i39 < loop_ub; i39++) {
          c_wt->data[i39] = wt->data[(int)l->data[i39] - 1];
        }

        c_rdivide(ad, c_wt, x2);
        i39 = b->size[0];
        b->size[0] = x2->size[1];
        emxEnsureCapacity((emxArray__common *)b, i39, (int)sizeof(double));
        loop_ub = x2->size[1];
        for (i39 = 0; i39 < loop_ub; i39++) {
          b->data[i39] = x2->data[x2->size[0] * i39];
        }

        if ((j->size[1] == 1) || (b->size[0] == 1)) {
          temp = 0.0;
          for (i39 = 0; i39 < j->size[1]; i39++) {
            temp += j->data[j->size[0] * i39] * b->data[i39];
          }
        } else {
          temp = 0.0;
          for (i39 = 0; i39 < j->size[1]; i39++) {
            temp += j->data[j->size[0] * i39] * b->data[i39];
          }
        }

        *dev = dnum / temp;
        nu = 1;
        if (*dev > 0.0) {
          nu = -1;
        }

        *dev *= -(double)nu;
        temp = (double)nu * *dev;
        i39 = b_a->size[0] * b_a->size[1];
        b_a->size[0] = 1;
        b_a->size[1] = j->size[1];
        emxEnsureCapacity((emxArray__common *)b_a, i39, (int)sizeof(double));
        loop_ub = j->size[0] * j->size[1];
        for (i39 = 0; i39 < loop_ub; i39++) {
          b_a->data[i39] = temp * j->data[i39];
        }

        i39 = b_wt->size[0] * b_wt->size[1];
        b_wt->size[0] = 1;
        b_wt->size[1] = l->size[1];
        emxEnsureCapacity((emxArray__common *)b_wt, i39, (int)sizeof(double));
        loop_ub = l->size[0] * l->size[1];
        for (i39 = 0; i39 < loop_ub; i39++) {
          b_wt->data[i39] = wt->data[(int)l->data[i39] - 1];
        }

        c_rdivide(b_a, b_wt, x2);
        i39 = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = l->size[1];
        emxEnsureCapacity((emxArray__common *)y, i39, (int)sizeof(double));
        loop_ub = l->size[0] * l->size[1];
        for (i39 = 0; i39 < loop_ub; i39++) {
          y->data[i39] = des->data[(int)l->data[i39] - 1] + x2->data[i39];
        }

        if (*dev <= devl) {
          /* warning(message('signal:firpm:DidNotConverge',niter)) */
          cfprintf("DidNotConverge");
          i39 = h->size[0] * h->size[1];
          h->size[0] = (int)nfilt;
          h->size[1] = 1;
          emxEnsureCapacity((emxArray__common *)h, i39, (int)sizeof(double));
          loop_ub = (int)nfilt;
          for (i39 = 0; i39 < loop_ub; i39++) {
            h->data[i39] = 0.0;
          }

          *dev = -1.0;

          /* iext */
          *valid = false;
          exitg1 = 1;
        } else {
          devl = *dev;
          jchnge = 0.0;
          c_k1 = iext->data[0];
          knz = iext->data[(int)(nfcns + 1.0) - 1];
          temp = 0.0;
          nut = -nu;
          f_j = 1.0;
          flag34 = 1;
          while (f_j < (nfcns + 1.0) + 1.0) {
            dnum = iext->data[(int)(unsigned int)f_j];
            b_l = iext->data[(int)f_j - 1] + 1.0;
            nut = -nut;
            if (f_j == 2.0) {
              b_y1 = comp;
            }

            comp = *dev;
            flag = 1;
            if (iext->data[(int)f_j - 1] + 1.0 < iext->data[(int)(f_j + 1.0) - 1])
            {
              /*  gee */
              err = cos(6.2831853071795862 * grid->data[(int)(iext->data[(int)
                         f_j - 1] + 1.0) - 1]);
              i39 = n_x->size[0] * n_x->size[1];
              n_x->size[0] = 1;
              n_x->size[1] = x->size[1];
              emxEnsureCapacity((emxArray__common *)n_x, i39, (int)sizeof(double));
              loop_ub = x->size[0] * x->size[1];
              for (i39 = 0; i39 < loop_ub; i39++) {
                n_x->data[i39] = err - x->data[i39];
              }

              c_rdivide(ad, n_x, j);
              i39 = b->size[0];
              b->size[0] = y->size[1];
              emxEnsureCapacity((emxArray__common *)b, i39, (int)sizeof(double));
              loop_ub = y->size[1];
              for (i39 = 0; i39 < loop_ub; i39++) {
                b->data[i39] = y->data[y->size[0] * i39];
              }

              if ((j->size[1] == 1) || (b->size[0] == 1)) {
                dtemp = 0.0;
                for (i39 = 0; i39 < j->size[1]; i39++) {
                  dtemp += j->data[j->size[0] * i39] * b->data[i39];
                }
              } else {
                dtemp = 0.0;
                for (i39 = 0; i39 < j->size[1]; i39++) {
                  dtemp += j->data[j->size[0] * i39] * b->data[i39];
                }
              }

              err = b_sum(j);
              err = (dtemp / err - des->data[(int)(iext->data[(int)f_j - 1] +
                      1.0) - 1]) * wt->data[(int)(iext->data[(int)f_j - 1] + 1.0)
                - 1];
              dtemp = (double)nut * err - *dev;
              if (dtemp > 0.0) {
                comp = (double)nut * err;
                b_l = (iext->data[(int)f_j - 1] + 1.0) + 1.0;
                exitg16 = false;
                while ((!exitg16) && (b_l < dnum)) {
                  /*  gee */
                  err = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                  i39 = m_x->size[0] * m_x->size[1];
                  m_x->size[0] = 1;
                  m_x->size[1] = x->size[1];
                  emxEnsureCapacity((emxArray__common *)m_x, i39, (int)sizeof
                                    (double));
                  loop_ub = x->size[0] * x->size[1];
                  for (i39 = 0; i39 < loop_ub; i39++) {
                    m_x->data[i39] = err - x->data[i39];
                  }

                  c_rdivide(ad, m_x, j);
                  i39 = b->size[0];
                  b->size[0] = y->size[1];
                  emxEnsureCapacity((emxArray__common *)b, i39, (int)sizeof
                                    (double));
                  loop_ub = y->size[1];
                  for (i39 = 0; i39 < loop_ub; i39++) {
                    b->data[i39] = y->data[y->size[0] * i39];
                  }

                  if ((j->size[1] == 1) || (b->size[0] == 1)) {
                    dtemp = 0.0;
                    for (i39 = 0; i39 < j->size[1]; i39++) {
                      dtemp += j->data[j->size[0] * i39] * b->data[i39];
                    }
                  } else {
                    dtemp = 0.0;
                    for (i39 = 0; i39 < j->size[1]; i39++) {
                      dtemp += j->data[j->size[0] * i39] * b->data[i39];
                    }
                  }

                  err = b_sum(j);
                  err = (dtemp / err - des->data[(int)b_l - 1]) * wt->data[(int)
                    b_l - 1];
                  dtemp = (double)nut * err - comp;
                  if (dtemp > 0.0) {
                    comp = (double)nut * err;
                    b_l++;
                  } else {
                    exitg16 = true;
                  }
                }

                iext->data[(int)f_j - 1] = b_l - 1.0;
                f_j++;
                temp = b_l - 1.0;
                jchnge++;
                flag = 0;
              }
            }

            if (flag != 0) {
              b_l -= 2.0;
              exitg15 = false;
              while ((!exitg15) && (b_l > temp)) {
                /*  gee */
                err = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                i39 = l_x->size[0] * l_x->size[1];
                l_x->size[0] = 1;
                l_x->size[1] = x->size[1];
                emxEnsureCapacity((emxArray__common *)l_x, i39, (int)sizeof
                                  (double));
                loop_ub = x->size[0] * x->size[1];
                for (i39 = 0; i39 < loop_ub; i39++) {
                  l_x->data[i39] = err - x->data[i39];
                }

                c_rdivide(ad, l_x, j);
                i39 = b->size[0];
                b->size[0] = y->size[1];
                emxEnsureCapacity((emxArray__common *)b, i39, (int)sizeof(double));
                loop_ub = y->size[1];
                for (i39 = 0; i39 < loop_ub; i39++) {
                  b->data[i39] = y->data[y->size[0] * i39];
                }

                if ((j->size[1] == 1) || (b->size[0] == 1)) {
                  dtemp = 0.0;
                  for (i39 = 0; i39 < j->size[1]; i39++) {
                    dtemp += j->data[j->size[0] * i39] * b->data[i39];
                  }
                } else {
                  dtemp = 0.0;
                  for (i39 = 0; i39 < j->size[1]; i39++) {
                    dtemp += j->data[j->size[0] * i39] * b->data[i39];
                  }
                }

                err = b_sum(j);
                err = (dtemp / err - des->data[(int)b_l - 1]) * wt->data[(int)
                  b_l - 1];
                dtemp = (double)nut * err - comp;
                if ((dtemp > 0.0) || (jchnge > 0.0)) {
                  exitg15 = true;
                } else {
                  b_l--;
                }
              }

              if (b_l <= temp) {
                b_l = iext->data[(int)f_j - 1] + 1.0;
                if (jchnge > 0.0) {
                  iext->data[(int)f_j - 1] = (iext->data[(int)f_j - 1] + 1.0) -
                    1.0;
                  f_j++;
                  temp = b_l - 1.0;
                  jchnge++;
                } else {
                  b_l = (iext->data[(int)f_j - 1] + 1.0) + 1.0;
                  exitg14 = false;
                  while ((!exitg14) && (b_l < dnum)) {
                    /*  gee */
                    err = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                    i39 = k_x->size[0] * k_x->size[1];
                    k_x->size[0] = 1;
                    k_x->size[1] = x->size[1];
                    emxEnsureCapacity((emxArray__common *)k_x, i39, (int)sizeof
                                      (double));
                    loop_ub = x->size[0] * x->size[1];
                    for (i39 = 0; i39 < loop_ub; i39++) {
                      k_x->data[i39] = err - x->data[i39];
                    }

                    c_rdivide(ad, k_x, j);
                    i39 = b->size[0];
                    b->size[0] = y->size[1];
                    emxEnsureCapacity((emxArray__common *)b, i39, (int)sizeof
                                      (double));
                    loop_ub = y->size[1];
                    for (i39 = 0; i39 < loop_ub; i39++) {
                      b->data[i39] = y->data[y->size[0] * i39];
                    }

                    if ((j->size[1] == 1) || (b->size[0] == 1)) {
                      dtemp = 0.0;
                      for (i39 = 0; i39 < j->size[1]; i39++) {
                        dtemp += j->data[j->size[0] * i39] * b->data[i39];
                      }
                    } else {
                      dtemp = 0.0;
                      for (i39 = 0; i39 < j->size[1]; i39++) {
                        dtemp += j->data[j->size[0] * i39] * b->data[i39];
                      }
                    }

                    err = b_sum(j);
                    err = (dtemp / err - des->data[(int)b_l - 1]) * wt->data
                      [(int)b_l - 1];
                    dtemp = (double)nut * err - comp;
                    if (dtemp > 0.0) {
                      exitg14 = true;
                    } else {
                      b_l++;
                    }
                  }

                  if ((b_l < dnum) && (dtemp > 0.0)) {
                    comp = (double)nut * err;
                    b_l++;
                    exitg13 = false;
                    while ((!exitg13) && (b_l < dnum)) {
                      /*  gee */
                      err = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                      i39 = j_x->size[0] * j_x->size[1];
                      j_x->size[0] = 1;
                      j_x->size[1] = x->size[1];
                      emxEnsureCapacity((emxArray__common *)j_x, i39, (int)
                                        sizeof(double));
                      loop_ub = x->size[0] * x->size[1];
                      for (i39 = 0; i39 < loop_ub; i39++) {
                        j_x->data[i39] = err - x->data[i39];
                      }

                      c_rdivide(ad, j_x, j);
                      i39 = b->size[0];
                      b->size[0] = y->size[1];
                      emxEnsureCapacity((emxArray__common *)b, i39, (int)sizeof
                                        (double));
                      loop_ub = y->size[1];
                      for (i39 = 0; i39 < loop_ub; i39++) {
                        b->data[i39] = y->data[y->size[0] * i39];
                      }

                      if ((j->size[1] == 1) || (b->size[0] == 1)) {
                        dtemp = 0.0;
                        for (i39 = 0; i39 < j->size[1]; i39++) {
                          dtemp += j->data[j->size[0] * i39] * b->data[i39];
                        }
                      } else {
                        dtemp = 0.0;
                        for (i39 = 0; i39 < j->size[1]; i39++) {
                          dtemp += j->data[j->size[0] * i39] * b->data[i39];
                        }
                      }

                      err = b_sum(j);
                      err = (dtemp / err - des->data[(int)b_l - 1]) * wt->data
                        [(int)b_l - 1];
                      dtemp = (double)nut * err - comp;
                      if (dtemp > 0.0) {
                        comp = (double)nut * err;
                        b_l++;
                      } else {
                        exitg13 = true;
                      }
                    }

                    iext->data[(int)f_j - 1] = b_l - 1.0;
                    f_j++;
                    temp = b_l - 1.0;
                    jchnge = 1.0;
                  } else {
                    temp = iext->data[(int)f_j - 1];
                    f_j++;
                  }
                }
              } else if (dtemp > 0.0) {
                comp = (double)nut * err;
                b_l--;
                exitg12 = false;
                while ((!exitg12) && (b_l > temp)) {
                  /*  gee */
                  err = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                  i39 = i_x->size[0] * i_x->size[1];
                  i_x->size[0] = 1;
                  i_x->size[1] = x->size[1];
                  emxEnsureCapacity((emxArray__common *)i_x, i39, (int)sizeof
                                    (double));
                  loop_ub = x->size[0] * x->size[1];
                  for (i39 = 0; i39 < loop_ub; i39++) {
                    i_x->data[i39] = err - x->data[i39];
                  }

                  c_rdivide(ad, i_x, j);
                  i39 = b->size[0];
                  b->size[0] = y->size[1];
                  emxEnsureCapacity((emxArray__common *)b, i39, (int)sizeof
                                    (double));
                  loop_ub = y->size[1];
                  for (i39 = 0; i39 < loop_ub; i39++) {
                    b->data[i39] = y->data[y->size[0] * i39];
                  }

                  if ((j->size[1] == 1) || (b->size[0] == 1)) {
                    dtemp = 0.0;
                    for (i39 = 0; i39 < j->size[1]; i39++) {
                      dtemp += j->data[j->size[0] * i39] * b->data[i39];
                    }
                  } else {
                    dtemp = 0.0;
                    for (i39 = 0; i39 < j->size[1]; i39++) {
                      dtemp += j->data[j->size[0] * i39] * b->data[i39];
                    }
                  }

                  err = b_sum(j);
                  err = (dtemp / err - des->data[(int)b_l - 1]) * wt->data[(int)
                    b_l - 1];
                  dtemp = (double)nut * err - comp;
                  if (dtemp > 0.0) {
                    comp = (double)nut * err;
                    b_l--;
                  } else {
                    exitg12 = true;
                  }
                }

                temp = iext->data[(int)f_j - 1];
                iext->data[(int)f_j - 1] = b_l + 1.0;
                f_j++;
                jchnge++;
              } else {
                temp = iext->data[(int)f_j - 1];
                f_j++;
              }
            }
          }

          do {
            exitg5 = 0;
            if (f_j == (nfcns + 1.0) + 1.0) {
              varargin_1[1] = iext->data[0];
              ixstart = 1;
              if (rtIsNaN(c_k1)) {
                flag = 2;
                exitg11 = false;
                while ((!exitg11) && (flag < 3)) {
                  ixstart = 2;
                  if (!rtIsNaN(varargin_1[1])) {
                    c_k1 = varargin_1[1];
                    exitg11 = true;
                  } else {
                    flag = 3;
                  }
                }
              }

              if ((ixstart < 2) && (varargin_1[1] < c_k1)) {
                c_k1 = varargin_1[1];
              }

              varargin_1[1] = iext->data[(int)(nfcns + 1.0) - 1];
              ixstart = 1;
              if (rtIsNaN(knz)) {
                flag = 2;
                exitg10 = false;
                while ((!exitg10) && (flag < 3)) {
                  ixstart = 2;
                  if (!rtIsNaN(varargin_1[1])) {
                    knz = varargin_1[1];
                    exitg10 = true;
                  } else {
                    flag = 3;
                  }
                }
              }

              if ((ixstart < 2) && (varargin_1[1] > knz)) {
                knz = varargin_1[1];
              }

              nut1 = nut;
              nut = -nu;
              comp *= 1.00001;
              luck = 1;
              flag = 1;
              b_l = 1.0;
              exitg8 = false;
              while ((!exitg8) && (b_l < c_k1)) {
                /*  gee */
                err = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                i39 = h_x->size[0] * h_x->size[1];
                h_x->size[0] = 1;
                h_x->size[1] = x->size[1];
                emxEnsureCapacity((emxArray__common *)h_x, i39, (int)sizeof
                                  (double));
                loop_ub = x->size[0] * x->size[1];
                for (i39 = 0; i39 < loop_ub; i39++) {
                  h_x->data[i39] = err - x->data[i39];
                }

                c_rdivide(ad, h_x, j);
                i39 = b->size[0];
                b->size[0] = y->size[1];
                emxEnsureCapacity((emxArray__common *)b, i39, (int)sizeof(double));
                loop_ub = y->size[1];
                for (i39 = 0; i39 < loop_ub; i39++) {
                  b->data[i39] = y->data[y->size[0] * i39];
                }

                if ((j->size[1] == 1) || (b->size[0] == 1)) {
                  dtemp = 0.0;
                  for (i39 = 0; i39 < j->size[1]; i39++) {
                    dtemp += j->data[j->size[0] * i39] * b->data[i39];
                  }
                } else {
                  dtemp = 0.0;
                  for (i39 = 0; i39 < j->size[1]; i39++) {
                    dtemp += j->data[j->size[0] * i39] * b->data[i39];
                  }
                }

                err = b_sum(j);
                err = (dtemp / err - des->data[(int)b_l - 1]) * wt->data[(int)
                  b_l - 1];
                dtemp = err * -(double)nu - comp;
                if (dtemp > 0.0) {
                  comp = -(double)nu * err;
                  b_l++;
                  exitg9 = false;
                  while ((!exitg9) && (b_l < c_k1)) {
                    /*  gee */
                    err = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                    i39 = g_x->size[0] * g_x->size[1];
                    g_x->size[0] = 1;
                    g_x->size[1] = x->size[1];
                    emxEnsureCapacity((emxArray__common *)g_x, i39, (int)sizeof
                                      (double));
                    loop_ub = x->size[0] * x->size[1];
                    for (i39 = 0; i39 < loop_ub; i39++) {
                      g_x->data[i39] = err - x->data[i39];
                    }

                    c_rdivide(ad, g_x, j);
                    i39 = b->size[0];
                    b->size[0] = y->size[1];
                    emxEnsureCapacity((emxArray__common *)b, i39, (int)sizeof
                                      (double));
                    loop_ub = y->size[1];
                    for (i39 = 0; i39 < loop_ub; i39++) {
                      b->data[i39] = y->data[y->size[0] * i39];
                    }

                    if ((j->size[1] == 1) || (b->size[0] == 1)) {
                      dtemp = 0.0;
                      for (i39 = 0; i39 < j->size[1]; i39++) {
                        dtemp += j->data[j->size[0] * i39] * b->data[i39];
                      }
                    } else {
                      dtemp = 0.0;
                      for (i39 = 0; i39 < j->size[1]; i39++) {
                        dtemp += j->data[j->size[0] * i39] * b->data[i39];
                      }
                    }

                    err = b_sum(j);
                    err = (dtemp / err - des->data[(int)b_l - 1]) * wt->data
                      [(int)b_l - 1];
                    dtemp = -(double)nu * err - comp;
                    if (dtemp > 0.0) {
                      comp = -(double)nu * err;
                      b_l++;
                    } else {
                      exitg9 = true;
                    }
                  }

                  iext->data[(int)((nfcns + 1.0) + 1.0) - 1] = b_l - 1.0;
                  f_j = ((nfcns + 1.0) + 1.0) + 1.0;
                  jchnge++;
                  flag = 0;
                  exitg8 = true;
                } else {
                  b_l++;
                }
              }

              if (flag != 0) {
                luck = 6;
                nut = -nut1;
                comp = b_y1 * 1.00001;
                b_l = ((double)ngrid + 1.0) - 1.0;
                exitg6 = false;
                while ((!exitg6) && (b_l > knz)) {
                  /*  gee */
                  err = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                  i39 = f_x->size[0] * f_x->size[1];
                  f_x->size[0] = 1;
                  f_x->size[1] = x->size[1];
                  emxEnsureCapacity((emxArray__common *)f_x, i39, (int)sizeof
                                    (double));
                  loop_ub = x->size[0] * x->size[1];
                  for (i39 = 0; i39 < loop_ub; i39++) {
                    f_x->data[i39] = err - x->data[i39];
                  }

                  c_rdivide(ad, f_x, j);
                  i39 = b->size[0];
                  b->size[0] = y->size[1];
                  emxEnsureCapacity((emxArray__common *)b, i39, (int)sizeof
                                    (double));
                  loop_ub = y->size[1];
                  for (i39 = 0; i39 < loop_ub; i39++) {
                    b->data[i39] = y->data[y->size[0] * i39];
                  }

                  if ((j->size[1] == 1) || (b->size[0] == 1)) {
                    dtemp = 0.0;
                    for (i39 = 0; i39 < j->size[1]; i39++) {
                      dtemp += j->data[j->size[0] * i39] * b->data[i39];
                    }
                  } else {
                    dtemp = 0.0;
                    for (i39 = 0; i39 < j->size[1]; i39++) {
                      dtemp += j->data[j->size[0] * i39] * b->data[i39];
                    }
                  }

                  err = b_sum(j);
                  err = (dtemp / err - des->data[(int)b_l - 1]) * wt->data[(int)
                    b_l - 1];
                  dtemp = err * -(double)nut1 - comp;
                  if (dtemp > 0.0) {
                    comp = -(double)nut1 * err;
                    luck = 16;
                    b_l--;
                    exitg7 = false;
                    while ((!exitg7) && (b_l > knz)) {
                      /*  gee */
                      err = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                      i39 = e_x->size[0] * e_x->size[1];
                      e_x->size[0] = 1;
                      e_x->size[1] = x->size[1];
                      emxEnsureCapacity((emxArray__common *)e_x, i39, (int)
                                        sizeof(double));
                      loop_ub = x->size[0] * x->size[1];
                      for (i39 = 0; i39 < loop_ub; i39++) {
                        e_x->data[i39] = err - x->data[i39];
                      }

                      c_rdivide(ad, e_x, j);
                      i39 = b->size[0];
                      b->size[0] = y->size[1];
                      emxEnsureCapacity((emxArray__common *)b, i39, (int)sizeof
                                        (double));
                      loop_ub = y->size[1];
                      for (i39 = 0; i39 < loop_ub; i39++) {
                        b->data[i39] = y->data[y->size[0] * i39];
                      }

                      if ((j->size[1] == 1) || (b->size[0] == 1)) {
                        dtemp = 0.0;
                        for (i39 = 0; i39 < j->size[1]; i39++) {
                          dtemp += j->data[j->size[0] * i39] * b->data[i39];
                        }
                      } else {
                        dtemp = 0.0;
                        for (i39 = 0; i39 < j->size[1]; i39++) {
                          dtemp += j->data[j->size[0] * i39] * b->data[i39];
                        }
                      }

                      err = b_sum(j);
                      err = (dtemp / err - des->data[(int)b_l - 1]) * wt->data
                        [(int)b_l - 1];
                      dtemp = -(double)nut1 * err - comp;
                      if (dtemp > 0.0) {
                        comp = -(double)nut1 * err;
                        b_l--;
                      } else {
                        exitg7 = true;
                      }
                    }

                    iext->data[(int)((nfcns + 1.0) + 1.0) - 1] = b_l + 1.0;
                    f_j = ((nfcns + 1.0) + 1.0) + 1.0;
                    jchnge++;
                    flag = 0;
                    exitg6 = true;
                  } else {
                    b_l--;
                  }
                }

                if (flag != 0) {
                  flag34 = 0;
                  if (luck != 6) {
                    temp = (nfcns + 1.0) - nfcns;
                    if (2.0 > temp) {
                      i39 = -2;
                      i40 = 0;
                    } else {
                      i39 = -1;
                      i40 = (int)temp;
                    }

                    temp = (nfcns + 1.0) - nfcns;
                    if (temp > (nfcns + 1.0) - 1.0) {
                      i41 = 1;
                      i42 = 0;
                    } else {
                      i41 = (int)temp;
                      i42 = (int)((nfcns + 1.0) - 1.0);
                    }

                    /*  Update index */
                    temp = (nfcns + 1.0) - nfcns;
                    if (2.0 > temp) {
                      i43 = -2;
                      i44 = 0;
                    } else {
                      i43 = -1;
                      i44 = (int)temp;
                    }

                    temp = (nfcns + 1.0) - nfcns;
                    if (temp > (nfcns + 1.0) - 1.0) {
                      ixstart = 1;
                      flag = 0;
                    } else {
                      ixstart = (int)temp;
                      flag = (int)((nfcns + 1.0) - 1.0);
                    }

                    nu = mtmp->size[0] * mtmp->size[1];
                    mtmp->size[0] = 1;
                    mtmp->size[1] = ((i40 - i39) + i42) - i41;
                    emxEnsureCapacity((emxArray__common *)mtmp, nu, (int)sizeof
                                      (double));
                    mtmp->data[0] = c_k1;
                    loop_ub = i40 - i39;
                    for (nu = 0; nu <= loop_ub - 3; nu++) {
                      mtmp->data[mtmp->size[0] * (nu + 1)] = iext->data[(i39 +
                        nu) + 2];
                    }

                    loop_ub = i42 - i41;
                    for (i42 = 0; i42 <= loop_ub; i42++) {
                      mtmp->data[mtmp->size[0] * (((i42 + i40) - i39) - 1)] =
                        iext->data[(i41 + i42) - 1];
                    }

                    loop_ub = mtmp->size[1];
                    i39 = r23->size[0] * r23->size[1];
                    r23->size[0] = 1;
                    r23->size[1] = loop_ub;
                    emxEnsureCapacity((emxArray__common *)r23, i39, (int)sizeof
                                      (int));
                    for (i39 = 0; i39 < loop_ub; i39++) {
                      r23->data[r23->size[0] * i39] = i39;
                    }

                    i39 = b_mtmp->size[0];
                    b_mtmp->size[0] = ((i44 - i43) + flag) - ixstart;
                    emxEnsureCapacity((emxArray__common *)b_mtmp, i39, (int)
                                      sizeof(double));
                    b_mtmp->data[0] = c_k1;
                    loop_ub = i44 - i43;
                    for (i39 = 0; i39 <= loop_ub - 3; i39++) {
                      b_mtmp->data[i39 + 1] = iext->data[(i43 + i39) + 2];
                    }

                    loop_ub = flag - ixstart;
                    for (i39 = 0; i39 <= loop_ub; i39++) {
                      b_mtmp->data[((i39 + i44) - i43) - 1] = iext->data
                        [(ixstart + i39) - 1];
                    }

                    loop_ub = r23->size[1];
                    for (i39 = 0; i39 < loop_ub; i39++) {
                      iext->data[r23->data[r23->size[0] * i39]] = b_mtmp->data[(*
                        (int (*)[2])r23->size)[0] * i39];
                    }

                    jchnge++;
                  }

                  exitg5 = 1;
                }
              }
            } else {
              exitg5 = 1;
            }
          } while (exitg5 == 0);

          if ((flag34 != 0) && (f_j > (nfcns + 1.0) + 1.0)) {
            if (luck > 9) {
              if (2.0 > nfcns + 1.0) {
                i39 = 0;
                i40 = 0;
              } else {
                i39 = 1;
                i40 = (int)(nfcns + 1.0);
              }

              if (nfcns + 1.0 > (nfcns + 1.0) - 1.0) {
                i41 = 0;
                i42 = 0;
              } else {
                i41 = (int)(nfcns + 1.0) - 1;
                i42 = (int)((nfcns + 1.0) - 1.0);
              }

              /*  Update index */
              if (2.0 > nfcns + 1.0) {
                i43 = 0;
                i44 = 0;
              } else {
                i43 = 1;
                i44 = (int)(nfcns + 1.0);
              }

              if (nfcns + 1.0 > (nfcns + 1.0) - 1.0) {
                ixstart = 0;
                flag = 0;
              } else {
                ixstart = (int)(nfcns + 1.0) - 1;
                flag = (int)((nfcns + 1.0) - 1.0);
              }

              nu = e_iext->size[0] * e_iext->size[1];
              e_iext->size[0] = 1;
              e_iext->size[1] = (((i40 - i39) + i42) - i41) + 2;
              emxEnsureCapacity((emxArray__common *)e_iext, nu, (int)sizeof
                                (double));
              loop_ub = i40 - i39;
              for (nu = 0; nu < loop_ub; nu++) {
                e_iext->data[e_iext->size[0] * nu] = iext->data[i39 + nu];
              }

              loop_ub = i42 - i41;
              for (nu = 0; nu < loop_ub; nu++) {
                e_iext->data[e_iext->size[0] * ((nu + i40) - i39)] = iext->
                  data[i41 + nu];
              }

              e_iext->data[e_iext->size[0] * (((i40 - i39) + i42) - i41)] =
                iext->data[(int)((nfcns + 1.0) + 1.0) - 1];
              e_iext->data[e_iext->size[0] * ((((i40 - i39) + i42) - i41) + 1)] =
                iext->data[(int)((nfcns + 1.0) + 1.0) - 1];
              loop_ub = e_iext->size[1];
              i39 = r23->size[0] * r23->size[1];
              r23->size[0] = 1;
              r23->size[1] = loop_ub;
              emxEnsureCapacity((emxArray__common *)r23, i39, (int)sizeof(int));
              for (i39 = 0; i39 < loop_ub; i39++) {
                r23->data[r23->size[0] * i39] = i39;
              }

              temp = iext->data[(int)((nfcns + 1.0) + 1.0) - 1];
              dnum = iext->data[(int)((nfcns + 1.0) + 1.0) - 1];
              i39 = f_iext->size[0];
              f_iext->size[0] = (((i44 - i43) + flag) - ixstart) + 2;
              emxEnsureCapacity((emxArray__common *)f_iext, i39, (int)sizeof
                                (double));
              loop_ub = i44 - i43;
              for (i39 = 0; i39 < loop_ub; i39++) {
                f_iext->data[i39] = iext->data[i43 + i39];
              }

              loop_ub = flag - ixstart;
              for (i39 = 0; i39 < loop_ub; i39++) {
                f_iext->data[(i39 + i44) - i43] = iext->data[ixstart + i39];
              }

              f_iext->data[((i44 - i43) + flag) - ixstart] = temp;
              f_iext->data[(((i44 - i43) + flag) - ixstart) + 1] = dnum;
              loop_ub = r23->size[1];
              for (i39 = 0; i39 < loop_ub; i39++) {
                iext->data[r23->data[r23->size[0] * i39]] = f_iext->data[(*(int
                  (*)[2])r23->size)[0] * i39];
              }

              jchnge++;
            } else {
              ixstart = 1;
              if (rtIsNaN(b_y1)) {
                flag = 2;
                exitg4 = false;
                while ((!exitg4) && (flag < 3)) {
                  ixstart = 2;
                  if (!rtIsNaN(comp)) {
                    b_y1 = comp;
                    exitg4 = true;
                  } else {
                    flag = 3;
                  }
                }
              }

              if ((ixstart < 2) && (comp > b_y1)) {
                b_y1 = comp;
              }

              c_k1 = iext->data[(int)((nfcns + 1.0) + 1.0) - 1];
              comp = b_y1 * 1.00001;
              b_l = ((double)ngrid + 1.0) - 1.0;
              exitg2 = false;
              while ((!exitg2) && (b_l > knz)) {
                /*  gee */
                err = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                i39 = d_x->size[0] * d_x->size[1];
                d_x->size[0] = 1;
                d_x->size[1] = x->size[1];
                emxEnsureCapacity((emxArray__common *)d_x, i39, (int)sizeof
                                  (double));
                loop_ub = x->size[0] * x->size[1];
                for (i39 = 0; i39 < loop_ub; i39++) {
                  d_x->data[i39] = err - x->data[i39];
                }

                c_rdivide(ad, d_x, j);
                i39 = b->size[0];
                b->size[0] = y->size[1];
                emxEnsureCapacity((emxArray__common *)b, i39, (int)sizeof(double));
                loop_ub = y->size[1];
                for (i39 = 0; i39 < loop_ub; i39++) {
                  b->data[i39] = y->data[y->size[0] * i39];
                }

                if ((j->size[1] == 1) || (b->size[0] == 1)) {
                  dtemp = 0.0;
                  for (i39 = 0; i39 < j->size[1]; i39++) {
                    dtemp += j->data[j->size[0] * i39] * b->data[i39];
                  }
                } else {
                  dtemp = 0.0;
                  for (i39 = 0; i39 < j->size[1]; i39++) {
                    dtemp += j->data[j->size[0] * i39] * b->data[i39];
                  }
                }

                err = b_sum(j);
                err = (dtemp / err - des->data[(int)b_l - 1]) * wt->data[(int)
                  b_l - 1];
                dtemp = err * -(double)nut1 - comp;
                if (dtemp > 0.0) {
                  comp = -(double)nut1 * err;
                  luck += 10;
                  b_l--;
                  exitg3 = false;
                  while ((!exitg3) && (b_l > knz)) {
                    /*  gee */
                    err = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                    i39 = c_x->size[0] * c_x->size[1];
                    c_x->size[0] = 1;
                    c_x->size[1] = x->size[1];
                    emxEnsureCapacity((emxArray__common *)c_x, i39, (int)sizeof
                                      (double));
                    loop_ub = x->size[0] * x->size[1];
                    for (i39 = 0; i39 < loop_ub; i39++) {
                      c_x->data[i39] = err - x->data[i39];
                    }

                    c_rdivide(ad, c_x, j);
                    i39 = b->size[0];
                    b->size[0] = y->size[1];
                    emxEnsureCapacity((emxArray__common *)b, i39, (int)sizeof
                                      (double));
                    loop_ub = y->size[1];
                    for (i39 = 0; i39 < loop_ub; i39++) {
                      b->data[i39] = y->data[y->size[0] * i39];
                    }

                    if ((j->size[1] == 1) || (b->size[0] == 1)) {
                      dtemp = 0.0;
                      for (i39 = 0; i39 < j->size[1]; i39++) {
                        dtemp += j->data[j->size[0] * i39] * b->data[i39];
                      }
                    } else {
                      dtemp = 0.0;
                      for (i39 = 0; i39 < j->size[1]; i39++) {
                        dtemp += j->data[j->size[0] * i39] * b->data[i39];
                      }
                    }

                    err = b_sum(j);
                    err = (dtemp / err - des->data[(int)b_l - 1]) * wt->data
                      [(int)b_l - 1];
                    dtemp = -(double)nut1 * err - comp;
                    if (dtemp > 0.0) {
                      comp = -(double)nut1 * err;
                      b_l--;
                    } else {
                      exitg3 = true;
                    }
                  }

                  iext->data[(int)((nfcns + 1.0) + 1.0) - 1] = b_l + 1.0;
                  jchnge++;
                  if (2.0 > nfcns + 1.0) {
                    i39 = 0;
                    i40 = 0;
                  } else {
                    i39 = 1;
                    i40 = (int)(nfcns + 1.0);
                  }

                  if (nfcns + 1.0 > (nfcns + 1.0) - 1.0) {
                    i41 = 0;
                    i42 = 0;
                  } else {
                    i41 = (int)(nfcns + 1.0) - 1;
                    i42 = (int)((nfcns + 1.0) - 1.0);
                  }

                  /*  Update index */
                  if (2.0 > nfcns + 1.0) {
                    i43 = 0;
                    i44 = 0;
                  } else {
                    i43 = 1;
                    i44 = (int)(nfcns + 1.0);
                  }

                  if (nfcns + 1.0 > (nfcns + 1.0) - 1.0) {
                    ixstart = 0;
                    flag = 0;
                  } else {
                    ixstart = (int)(nfcns + 1.0) - 1;
                    flag = (int)((nfcns + 1.0) - 1.0);
                  }

                  nu = c_iext->size[0] * c_iext->size[1];
                  c_iext->size[0] = 1;
                  c_iext->size[1] = (((i40 - i39) + i42) - i41) + 2;
                  emxEnsureCapacity((emxArray__common *)c_iext, nu, (int)sizeof
                                    (double));
                  loop_ub = i40 - i39;
                  for (nu = 0; nu < loop_ub; nu++) {
                    c_iext->data[c_iext->size[0] * nu] = iext->data[i39 + nu];
                  }

                  loop_ub = i42 - i41;
                  for (nu = 0; nu < loop_ub; nu++) {
                    c_iext->data[c_iext->size[0] * ((nu + i40) - i39)] =
                      iext->data[i41 + nu];
                  }

                  c_iext->data[c_iext->size[0] * (((i40 - i39) + i42) - i41)] =
                    iext->data[(int)((nfcns + 1.0) + 1.0) - 1];
                  c_iext->data[c_iext->size[0] * ((((i40 - i39) + i42) - i41) +
                    1)] = iext->data[(int)((nfcns + 1.0) + 1.0) - 1];
                  loop_ub = c_iext->size[1];
                  i39 = r23->size[0] * r23->size[1];
                  r23->size[0] = 1;
                  r23->size[1] = loop_ub;
                  emxEnsureCapacity((emxArray__common *)r23, i39, (int)sizeof
                                    (int));
                  for (i39 = 0; i39 < loop_ub; i39++) {
                    r23->data[r23->size[0] * i39] = i39;
                  }

                  temp = iext->data[(int)((nfcns + 1.0) + 1.0) - 1];
                  dnum = iext->data[(int)((nfcns + 1.0) + 1.0) - 1];
                  i39 = d_iext->size[0];
                  d_iext->size[0] = (((i44 - i43) + flag) - ixstart) + 2;
                  emxEnsureCapacity((emxArray__common *)d_iext, i39, (int)sizeof
                                    (double));
                  loop_ub = i44 - i43;
                  for (i39 = 0; i39 < loop_ub; i39++) {
                    d_iext->data[i39] = iext->data[i43 + i39];
                  }

                  loop_ub = flag - ixstart;
                  for (i39 = 0; i39 < loop_ub; i39++) {
                    d_iext->data[(i39 + i44) - i43] = iext->data[ixstart + i39];
                  }

                  d_iext->data[((i44 - i43) + flag) - ixstart] = temp;
                  d_iext->data[(((i44 - i43) + flag) - ixstart) + 1] = dnum;
                  loop_ub = r23->size[1];
                  for (i39 = 0; i39 < loop_ub; i39++) {
                    iext->data[r23->data[r23->size[0] * i39]] = d_iext->data
                      [(*(int (*)[2])r23->size)[0] * i39];
                  }

                  exitg2 = true;
                } else {
                  b_l--;
                }
              }

              if (luck != 6) {
                temp = (nfcns + 1.0) - nfcns;
                if (2.0 > temp) {
                  i39 = -2;
                  i40 = 0;
                } else {
                  i39 = -1;
                  i40 = (int)temp;
                }

                temp = (nfcns + 1.0) - nfcns;
                if (temp > (nfcns + 1.0) - 1.0) {
                  i41 = 1;
                  i42 = 0;
                } else {
                  i41 = (int)temp;
                  i42 = (int)((nfcns + 1.0) - 1.0);
                }

                /*  Update index */
                temp = (nfcns + 1.0) - nfcns;
                if (2.0 > temp) {
                  i43 = -2;
                  i44 = 0;
                } else {
                  i43 = -1;
                  i44 = (int)temp;
                }

                temp = (nfcns + 1.0) - nfcns;
                if (temp > (nfcns + 1.0) - 1.0) {
                  ixstart = 1;
                  flag = 0;
                } else {
                  ixstart = (int)temp;
                  flag = (int)((nfcns + 1.0) - 1.0);
                }

                nu = k1->size[0] * k1->size[1];
                k1->size[0] = 1;
                k1->size[1] = ((i40 - i39) + i42) - i41;
                emxEnsureCapacity((emxArray__common *)k1, nu, (int)sizeof(double));
                k1->data[0] = c_k1;
                loop_ub = i40 - i39;
                for (nu = 0; nu <= loop_ub - 3; nu++) {
                  k1->data[k1->size[0] * (nu + 1)] = iext->data[(i39 + nu) + 2];
                }

                loop_ub = i42 - i41;
                for (i42 = 0; i42 <= loop_ub; i42++) {
                  k1->data[k1->size[0] * (((i42 + i40) - i39) - 1)] = iext->
                    data[(i41 + i42) - 1];
                }

                loop_ub = k1->size[1];
                i39 = r23->size[0] * r23->size[1];
                r23->size[0] = 1;
                r23->size[1] = loop_ub;
                emxEnsureCapacity((emxArray__common *)r23, i39, (int)sizeof(int));
                for (i39 = 0; i39 < loop_ub; i39++) {
                  r23->data[r23->size[0] * i39] = i39;
                }

                i39 = b_k1->size[0];
                b_k1->size[0] = ((i44 - i43) + flag) - ixstart;
                emxEnsureCapacity((emxArray__common *)b_k1, i39, (int)sizeof
                                  (double));
                b_k1->data[0] = c_k1;
                loop_ub = i44 - i43;
                for (i39 = 0; i39 <= loop_ub - 3; i39++) {
                  b_k1->data[i39 + 1] = iext->data[(i43 + i39) + 2];
                }

                loop_ub = flag - ixstart;
                for (i39 = 0; i39 <= loop_ub; i39++) {
                  b_k1->data[((i39 + i44) - i43) - 1] = iext->data[(ixstart +
                    i39) - 1];
                }

                loop_ub = r23->size[1];
                for (i39 = 0; i39 < loop_ub; i39++) {
                  iext->data[r23->data[r23->size[0] * i39]] = b_k1->data[(*(int
                    (*)[2])r23->size)[0] * i39];
                }

                jchnge++;
              }
            }
          }

          guard1 = false;
        }
      }
    } else {
      guard1 = true;
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  if (guard1) {
    /*  Inverse Fourier transformation */
    f_j = -1.0;
    c_k1 = -1.0;

    /*  initialize memory */
    /* x(nzz) = -2; */
    i39 = x2->size[0] * x2->size[1];
    x2->size[0] = 1;
    x2->size[1] = x->size[1] + 1;
    emxEnsureCapacity((emxArray__common *)x2, i39, (int)sizeof(double));
    loop_ub = x->size[1];
    for (i39 = 0; i39 < loop_ub; i39++) {
      x2->data[x2->size[0] * i39] = x->data[x->size[0] * i39];
    }

    x2->data[x2->size[0] * x->size[1]] = -2.0;
    jchnge = 2.0 * nfcns - 1.0;
    knz = 1.0 / jchnge;
    b_l = 1.0;
    nu = 0;
    if (((edge[0] == 0.0) && (edge[3] == 0.5)) || (nfcns <= 3.0)) {
      nu = 1;
    }

    if (nu != 1) {
      dtemp = cos(6.2831853071795862 * grid->data[0]);
      dnum = cos(6.2831853071795862 * grid->data[varargin_2 - 1]);
      c_k1 = 2.0 / (dtemp - dnum);
      f_j = -(dtemp + dnum) / (dtemp - dnum);
    }

    i39 = a->size[0] * a->size[1];
    a->size[0] = 1;
    a->size[1] = (int)nfcns;
    emxEnsureCapacity((emxArray__common *)a, i39, (int)sizeof(double));
    for (nut = 0; nut < (int)nfcns; nut++) {
      temp = ((1.0 + (double)nut) - 1.0) * knz;
      dnum = cos(6.2831853071795862 * temp);
      if (nu != 1) {
        dnum = (dnum - f_j) / c_k1;
        temp = acos(dnum) / 6.2831853071795862;
      }

      err = x2->data[(int)b_l - 1];
      while ((dnum <= err) && (err - dnum >= 1.0E-6)) {
        b_l++;
        err = x2->data[(int)b_l - 1];
      }

      if (fabs(dnum - err) < 1.0E-6) {
        a->data[nut] = y->data[(int)b_l - 1];
      } else {
        /*  gee */
        loop_ub = (int)(nfcns + 1.0);
        err = cos(6.2831853071795862 * temp);
        i39 = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = (int)(nfcns + 1.0);
        emxEnsureCapacity((emxArray__common *)b_x, i39, (int)sizeof(double));
        for (i39 = 0; i39 < loop_ub; i39++) {
          b_x->data[b_x->size[0] * i39] = err - x2->data[i39];
        }

        c_rdivide(ad, b_x, j);
        i39 = b->size[0];
        b->size[0] = y->size[1];
        emxEnsureCapacity((emxArray__common *)b, i39, (int)sizeof(double));
        loop_ub = y->size[1];
        for (i39 = 0; i39 < loop_ub; i39++) {
          b->data[i39] = y->data[y->size[0] * i39];
        }

        if ((j->size[1] == 1) || (b->size[0] == 1)) {
          dtemp = 0.0;
          for (i39 = 0; i39 < j->size[1]; i39++) {
            dtemp += j->data[j->size[0] * i39] * b->data[i39];
          }
        } else {
          dtemp = 0.0;
          for (i39 = 0; i39 < j->size[1]; i39++) {
            dtemp += j->data[j->size[0] * i39] * b->data[i39];
          }
        }

        err = b_sum(j);
        a->data[nut] = dtemp / err;
      }

      ixstart = 1;
      if ((int)(b_l - 1.0) > 1) {
        ixstart = (int)(b_l - 1.0);
      }

      b_l = ixstart;
    }

    temp = 6.2831853071795862 / jchnge;
    i39 = j->size[0] * j->size[1];
    j->size[0] = 1;
    j->size[1] = (int)nfcns;
    emxEnsureCapacity((emxArray__common *)j, i39, (int)sizeof(double));
    for (nut = 0; nut < (int)nfcns; nut++) {
      dnum = ((1.0 + (double)nut) - 1.0) * temp;
      if (nfcns - 1.0 < 1.0) {
        j->data[nut] = a->data[0];
      } else {
        if (2.0 > nfcns) {
          i39 = 1;
          i40 = 1;
        } else {
          i39 = 2;
          i40 = (int)nfcns + 1;
        }

        if (nfcns - 1.0 < 1.0) {
          i41 = y->size[0] * y->size[1];
          y->size[0] = 1;
          y->size[1] = 0;
          emxEnsureCapacity((emxArray__common *)y, i41, (int)sizeof(double));
        } else {
          i41 = y->size[0] * y->size[1];
          y->size[0] = 1;
          y->size[1] = (int)floor((nfcns - 1.0) - 1.0) + 1;
          emxEnsureCapacity((emxArray__common *)y, i41, (int)sizeof(double));
          loop_ub = (int)floor((nfcns - 1.0) - 1.0);
          for (i41 = 0; i41 <= loop_ub; i41++) {
            y->data[y->size[0] * i41] = 1.0 + (double)i41;
          }
        }

        i41 = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)y, i41, (int)sizeof(double));
        ixstart = y->size[0];
        flag = y->size[1];
        loop_ub = ixstart * flag;
        for (i41 = 0; i41 < loop_ub; i41++) {
          y->data[i41] *= dnum;
        }

        b_cos(y);
        i41 = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)y, i41, (int)sizeof(double));
        ixstart = y->size[0];
        flag = y->size[1];
        loop_ub = ixstart * flag;
        for (i41 = 0; i41 < loop_ub; i41++) {
          y->data[i41] *= 2.0;
        }

        i41 = b->size[0];
        b->size[0] = i40 - i39;
        emxEnsureCapacity((emxArray__common *)b, i41, (int)sizeof(double));
        loop_ub = i40 - i39;
        for (i41 = 0; i41 < loop_ub; i41++) {
          b->data[i41] = a->data[(i39 + i41) - 1];
        }

        if ((y->size[1] == 1) || (i40 - i39 == 1)) {
          dtemp = 0.0;
          for (i39 = 0; i39 < y->size[1]; i39++) {
            dtemp += y->data[y->size[0] * i39] * b->data[i39];
          }
        } else {
          dtemp = 0.0;
          for (i39 = 0; i39 < y->size[1]; i39++) {
            dtemp += y->data[y->size[0] * i39] * b->data[i39];
          }
        }

        j->data[nut] = a->data[0] + dtemp;
      }
    }

    if (2.0 > nfcns) {
      i39 = -1;
      i40 = 0;
    } else {
      i39 = 0;
      i40 = (int)nfcns;
    }

    i41 = b_j->size[0];
    b_j->size[0] = i40 - i39;
    emxEnsureCapacity((emxArray__common *)b_j, i41, (int)sizeof(double));
    b_j->data[0] = j->data[0];
    loop_ub = i40 - i39;
    for (i40 = 0; i40 <= loop_ub - 2; i40++) {
      b_j->data[i40 + 1] = 2.0 * j->data[(i39 + i40) + 1];
    }

    i39 = iext->size[0];
    iext->size[0] = b_j->size[0];
    emxEnsureCapacity((emxArray__common *)iext, i39, (int)sizeof(double));
    loop_ub = b_j->size[0];
    for (i39 = 0; i39 < loop_ub; i39++) {
      iext->data[i39] = b_j->data[i39] / jchnge;
    }

    i39 = j->size[0] * j->size[1];
    j->size[0] = 1;
    j->size[1] = (int)nfcns;
    emxEnsureCapacity((emxArray__common *)j, i39, (int)sizeof(double));
    loop_ub = (int)nfcns;
    for (i39 = 0; i39 < loop_ub; i39++) {
      j->data[i39] = 0.0;
    }

    i39 = l->size[0] * l->size[1];
    l->size[0] = 1;
    l->size[1] = (int)nfcns - 1;
    emxEnsureCapacity((emxArray__common *)l, i39, (int)sizeof(double));
    loop_ub = (int)nfcns - 1;
    for (i39 = 0; i39 < loop_ub; i39++) {
      l->data[i39] = 0.0;
    }

    if (nu != 1) {
      j->data[0] = 2.0 * iext->data[(int)nfcns - 1] * f_j + iext->data[(int)
        nfcns - 2];
      j->data[1] = 2.0 * c_k1 * iext->data[(int)nfcns - 1];
      l->data[0] = iext->data[(int)nfcns - 3] - iext->data[(int)nfcns - 1];
      for (nut = 0; nut <= (int)nfcns - 3; nut++) {
        if (2 + nut == (int)nfcns - 1) {
          c_k1 /= 2.0;
          f_j /= 2.0;
        }

        j->data[2 + nut] = 0.0;
        i39 = x2->size[0] * x2->size[1];
        x2->size[0] = 1;
        x2->size[1] = (int)((2.0 + (double)nut) - 1.0) + 1;
        emxEnsureCapacity((emxArray__common *)x2, i39, (int)sizeof(double));
        loop_ub = (int)((2.0 + (double)nut) - 1.0);
        for (i39 = 0; i39 <= loop_ub; i39++) {
          x2->data[x2->size[0] * i39] = 1.0 + (double)i39;
        }

        loop_ub = x2->size[0] * x2->size[1];
        for (i39 = 0; i39 < loop_ub; i39++) {
          a->data[(int)x2->data[i39] - 1] = j->data[(int)x2->data[i39] - 1];
        }

        temp = 2.0 * f_j;
        loop_ub = x2->size[0] * x2->size[1];
        for (i39 = 0; i39 < loop_ub; i39++) {
          j->data[(int)x2->data[i39] - 1] = temp * a->data[(int)x2->data[i39] -
            1];
        }

        j->data[1] += 2.0 * a->data[0] * c_k1;
        i39 = x2->size[0] * x2->size[1];
        x2->size[0] = 1;
        x2->size[1] = (int)(((2.0 + (double)nut) - 1.0) - 1.0) + 1;
        emxEnsureCapacity((emxArray__common *)x2, i39, (int)sizeof(double));
        loop_ub = (int)(((2.0 + (double)nut) - 1.0) - 1.0);
        for (i39 = 0; i39 <= loop_ub; i39++) {
          x2->data[x2->size[0] * i39] = 1.0 + (double)i39;
        }

        i39 = x->size[0] * x->size[1];
        x->size[0] = 1;
        x->size[1] = x2->size[1];
        emxEnsureCapacity((emxArray__common *)x, i39, (int)sizeof(double));
        loop_ub = x2->size[0] * x2->size[1];
        for (i39 = 0; i39 < loop_ub; i39++) {
          x->data[i39] = a->data[(int)x2->data[i39]];
        }

        i39 = d_j->size[0];
        d_j->size[0] = x2->size[0] * x2->size[1];
        emxEnsureCapacity((emxArray__common *)d_j, i39, (int)sizeof(double));
        loop_ub = x2->size[0] * x2->size[1];
        for (i39 = 0; i39 < loop_ub; i39++) {
          d_j->data[i39] = (j->data[(int)x2->data[i39] - 1] + l->data[(int)
                            x2->data[i39] - 1]) + c_k1 * x->data[i39];
        }

        loop_ub = d_j->size[0];
        for (i39 = 0; i39 < loop_ub; i39++) {
          j->data[(int)x2->data[i39] - 1] = d_j->data[i39];
        }

        i39 = x2->size[0] * x2->size[1];
        x2->size[0] = 1;
        x2->size[1] = (int)(((2.0 + (double)nut) + 1.0) - 3.0) + 1;
        emxEnsureCapacity((emxArray__common *)x2, i39, (int)sizeof(double));
        loop_ub = (int)(((2.0 + (double)nut) + 1.0) - 3.0);
        for (i39 = 0; i39 <= loop_ub; i39++) {
          x2->data[x2->size[0] * i39] = 3.0 + (double)i39;
        }

        i39 = x->size[0] * x->size[1];
        x->size[0] = 1;
        x->size[1] = x2->size[1];
        emxEnsureCapacity((emxArray__common *)x, i39, (int)sizeof(double));
        loop_ub = x2->size[0] * x2->size[1];
        for (i39 = 0; i39 < loop_ub; i39++) {
          x->data[i39] = a->data[(int)x2->data[i39] - 2];
        }

        i39 = e_j->size[0];
        e_j->size[0] = x2->size[0] * x2->size[1];
        emxEnsureCapacity((emxArray__common *)e_j, i39, (int)sizeof(double));
        loop_ub = x2->size[0] * x2->size[1];
        for (i39 = 0; i39 < loop_ub; i39++) {
          e_j->data[i39] = j->data[(int)x2->data[i39] - 1] + c_k1 * x->data[i39];
        }

        loop_ub = e_j->size[0];
        for (i39 = 0; i39 < loop_ub; i39++) {
          j->data[(int)x2->data[i39] - 1] = e_j->data[i39];
        }

        if (2 + nut != (int)nfcns - 1) {
          i39 = x2->size[0] * x2->size[1];
          x2->size[0] = 1;
          x2->size[1] = (int)((2.0 + (double)nut) - 1.0) + 1;
          emxEnsureCapacity((emxArray__common *)x2, i39, (int)sizeof(double));
          loop_ub = (int)((2.0 + (double)nut) - 1.0);
          for (i39 = 0; i39 <= loop_ub; i39++) {
            x2->data[x2->size[0] * i39] = 1.0 + (double)i39;
          }

          loop_ub = x2->size[0] * x2->size[1];
          for (i39 = 0; i39 < loop_ub; i39++) {
            l->data[(int)x2->data[i39] - 1] = -a->data[(int)x2->data[i39] - 1];
          }

          l->data[0] += iext->data[((int)nfcns - nut) - 4];
        }
      }

      loop_ub = (int)nfcns;
      for (i39 = 0; i39 < loop_ub; i39++) {
        iext->data[i39] = j->data[i39];
      }
    }

    /*  alpha must be at lease >=3 */
    if (nfcns <= 3.0) {
      /* alpha(nfcns + 1) = 0; */
      /* alpha(nfcns + 2) = 0; */
      i39 = j->size[0] * j->size[1];
      j->size[0] = 1;
      j->size[1] = iext->size[0] + 2;
      emxEnsureCapacity((emxArray__common *)j, i39, (int)sizeof(double));
      loop_ub = iext->size[0];
      for (i39 = 0; i39 < loop_ub; i39++) {
        j->data[j->size[0] * i39] = iext->data[i39];
      }

      j->data[j->size[0] * iext->size[0]] = 0.0;
      j->data[j->size[0] * (iext->size[0] + 1)] = 0.0;
    } else {
      i39 = j->size[0] * j->size[1];
      j->size[0] = 1;
      j->size[1] = iext->size[0];
      emxEnsureCapacity((emxArray__common *)j, i39, (int)sizeof(double));
      loop_ub = iext->size[0];
      for (i39 = 0; i39 < loop_ub; i39++) {
        j->data[j->size[0] * i39] = iext->data[i39];
      }
    }

    /* alpha=alpha'; */
    /*  now that's done! */
    if (nodd != 0.0) {
      i39 = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = (int)floor(-((0.0 - (nfcns - 1.0)) - -1.0)) + 1;
      emxEnsureCapacity((emxArray__common *)x, i39, (int)sizeof(double));
      loop_ub = (int)floor(-((0.0 - (nfcns - 1.0)) - -1.0));
      for (i39 = 0; i39 <= loop_ub; i39++) {
        x->data[x->size[0] * i39] = j->data[(int)((nfcns + 1.0) + (-1.0 -
          (double)i39)) - 1];
      }

      i39 = h->size[0] * h->size[1];
      h->size[0] = 1;
      h->size[1] = x->size[1] + 1;
      emxEnsureCapacity((emxArray__common *)h, i39, (int)sizeof(double));
      loop_ub = x->size[1];
      for (i39 = 0; i39 < loop_ub; i39++) {
        h->data[h->size[0] * i39] = 0.5 * x->data[x->size[0] * i39];
      }

      h->data[h->size[0] * x->size[1]] = j->data[0];
    } else {
      if ((nfcns - (nfcns - 1.0)) + 2.0 > nfcns) {
        i39 = 0;
        i40 = 1;
      } else {
        i39 = (int)nfcns - 1;
        i40 = -1;
      }

      i41 = x2->size[0] * x2->size[1];
      x2->size[0] = 1;
      x2->size[1] = (int)floor(-((0.0 - (nfcns - 1.0)) - -2.0)) + 1;
      emxEnsureCapacity((emxArray__common *)x2, i41, (int)sizeof(double));
      loop_ub = (int)floor(-((0.0 - (nfcns - 1.0)) - -2.0));
      for (i41 = 0; i41 <= loop_ub; i41++) {
        x2->data[x2->size[0] * i41] = j->data[(int)((nfcns + 1.0) + (double)(int)
          (-2.0 - (double)i41)) - 1];
      }

      i41 = c_j->size[0] * c_j->size[1];
      c_j->size[0] = 1;
      c_j->size[1] = 2 + x2->size[1];
      emxEnsureCapacity((emxArray__common *)c_j, i41, (int)sizeof(double));
      c_j->data[0] = j->data[(int)nfcns - 1];
      loop_ub = x2->size[1];
      for (i41 = 0; i41 < loop_ub; i41++) {
        c_j->data[c_j->size[0] * (i41 + 1)] = x2->data[x2->size[0] * i41] +
          j->data[i39 + i40 * i41];
      }

      c_j->data[c_j->size[0] * (1 + x2->size[1])] = 2.0 * j->data[0] + j->data[1];
      i39 = h->size[0] * h->size[1];
      h->size[0] = 1;
      h->size[1] = c_j->size[1];
      emxEnsureCapacity((emxArray__common *)h, i39, (int)sizeof(double));
      loop_ub = c_j->size[1];
      for (i39 = 0; i39 < loop_ub; i39++) {
        h->data[h->size[0] * i39] = 0.25 * c_j->data[c_j->size[0] * i39];
      }
    }
  }

  emxFree_real_T(&e_j);
  emxFree_real_T(&d_j);
  emxFree_real_T(&c_j);
  emxFree_real_T(&b_j);
  emxFree_real_T(&f_iext);
  emxFree_real_T(&e_iext);
  emxFree_real_T(&b_k1);
  emxFree_real_T(&k1);
  emxFree_real_T(&d_iext);
  emxFree_real_T(&c_iext);
  emxFree_real_T(&b_mtmp);
  emxFree_real_T(&mtmp);
  emxFree_int32_T(&b_iext);
  emxFree_int8_T(&r24);
  emxFree_real_T(&n_x);
  emxFree_real_T(&m_x);
  emxFree_real_T(&l_x);
  emxFree_real_T(&k_x);
  emxFree_real_T(&j_x);
  emxFree_real_T(&i_x);
  emxFree_real_T(&h_x);
  emxFree_real_T(&g_x);
  emxFree_real_T(&f_x);
  emxFree_real_T(&e_x);
  emxFree_real_T(&d_x);
  emxFree_real_T(&c_x);
  emxFree_real_T(&c_wt);
  emxFree_real_T(&b_wt);
  emxFree_real_T(&b_a);
  emxFree_real_T(&b_x);
  emxFree_real_T(&b);
  emxFree_int32_T(&r23);
  emxFree_real_T(&a);
  emxFree_real_T(&x2);
  emxFree_real_T(&l);
  emxFree_real_T(&ad);
  emxFree_real_T(&x);
  emxFree_real_T(&y);
  emxFree_real_T(&iext);
  emxFree_real_T(&j);
}

/*
 * Arguments    : const double b_data[]
 *                const int b_size[2]
 *                const emxArray_creal_T *a
 *                double bR_data[]
 *                int bR_size[2]
 *                emxArray_creal_T *aR
 * Return Type  : void
 */
static void removeTrailingZero(const double b_data[], const int b_size[2], const
  emxArray_creal_T *a, double bR_data[], int bR_size[2], emxArray_creal_T *aR)
{
  int idx;
  int ii_data[1];
  signed char ii_size[2];
  int i22;
  int ii;
  boolean_T exitg2;
  int b_ii_data[1];
  signed char b_ii_size[2];
  boolean_T exitg1;
  boolean_T b_a;
  signed char varargin_1_data[1];
  int varargin_2_data[1];
  signed char csz_idx_1;
  int trz_len_data_idx_0;
  double d0;

  /*  end freqs */
  idx = 0;
  for (i22 = 0; i22 < 2; i22++) {
    ii_size[i22] = 1;
  }

  ii = b_size[1];
  exitg2 = false;
  while ((!exitg2) && (ii > 0)) {
    if (b_data[ii - 1] != 0.0) {
      idx = 1;
      ii_data[0] = ii;
      exitg2 = true;
    } else {
      ii--;
    }
  }

  if (idx == 0) {
    ii_size[1] = 0;
  }

  idx = 0;
  for (i22 = 0; i22 < 2; i22++) {
    b_ii_size[i22] = 1;
  }

  ii = a->size[1];
  exitg1 = false;
  while ((!exitg1) && (ii > 0)) {
    b_a = ((a->data[ii - 1].re != 0.0) || (a->data[ii - 1].im != 0.0));
    if (b_a) {
      idx = 1;
      b_ii_data[0] = ii;
      exitg1 = true;
    } else {
      ii--;
    }
  }

  if (idx == 0) {
    b_ii_size[1] = 0;
  }

  ii = ii_size[1];
  for (i22 = 0; i22 < ii; i22++) {
    varargin_1_data[i22] = (signed char)((signed char)b_size[1] - (signed char)
      ii_data[i22]);
  }

  idx = a->size[1];
  ii = b_ii_size[1];
  for (i22 = 0; i22 < ii; i22++) {
    varargin_2_data[i22] = idx - b_ii_data[i22];
  }

  if (ii_size[1] <= b_ii_size[1]) {
    csz_idx_1 = ii_size[1];
  } else {
    csz_idx_1 = 0;
  }

  ii = 1;
  while (ii <= csz_idx_1) {
    if (varargin_1_data[0] <= varargin_2_data[0]) {
      trz_len_data_idx_0 = varargin_1_data[0];
    } else {
      trz_len_data_idx_0 = varargin_2_data[0];
    }

    ii = 2;
  }

  if (trz_len_data_idx_0 > 0) {
    d0 = (double)b_size[1] - (double)trz_len_data_idx_0;
    if (1.0 > d0) {
      ii = 0;
    } else {
      ii = (int)d0;
    }

    bR_size[0] = 1;
    bR_size[1] = ii;
    for (i22 = 0; i22 < ii; i22++) {
      bR_data[bR_size[0] * i22] = b_data[i22];
    }

    d0 = (double)a->size[1] - (double)trz_len_data_idx_0;
    if (1.0 > d0) {
      ii = 0;
    } else {
      ii = (int)d0;
    }

    i22 = aR->size[0] * aR->size[1];
    aR->size[0] = 1;
    aR->size[1] = ii;
    emxEnsureCapacity((emxArray__common *)aR, i22, (int)sizeof(creal_T));
    for (i22 = 0; i22 < ii; i22++) {
      aR->data[aR->size[0] * i22] = a->data[i22];
    }
  } else {
    bR_size[0] = 1;
    bR_size[1] = b_size[1];
    ii = b_size[0] * b_size[1];
    for (i22 = 0; i22 < ii; i22++) {
      bR_data[i22] = b_data[i22];
    }

    i22 = aR->size[0] * aR->size[1];
    aR->size[0] = 1;
    aR->size[1] = a->size[1];
    emxEnsureCapacity((emxArray__common *)aR, i22, (int)sizeof(creal_T));
    ii = a->size[0] * a->size[1];
    for (i22 = 0; i22 < ii; i22++) {
      aR->data[i22] = a->data[i22];
    }
  }

  /*  end removeTrailingZero */
}

/*
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_atan2d_snf(double u0, double u1)
{
  double y;
  int b_u0;
  int b_u1;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    if (u0 > 0.0) {
      b_u0 = 1;
    } else {
      b_u0 = -1;
    }

    if (u1 > 0.0) {
      b_u1 = 1;
    } else {
      b_u1 = -1;
    }

    y = atan2(b_u0, b_u1);
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = atan2(u0, u1);
  }

  return y;
}

/*
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_hypotd_snf(double u0, double u1)
{
  double y;
  double a;
  double b;
  a = fabs(u0);
  b = fabs(u1);
  if (a < b) {
    a /= b;
    y = b * sqrt(a * a + 1.0);
  } else if (a > b) {
    b /= a;
    y = a * sqrt(b * b + 1.0);
  } else if (rtIsNaN(b)) {
    y = b;
  } else {
    y = a * 1.4142135623730951;
  }

  return y;
}

/*
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_powd_snf(double u0, double u1)
{
  double y;
  double d2;
  double d3;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    d2 = fabs(u0);
    d3 = fabs(u1);
    if (rtIsInf(u1)) {
      if (d2 == 1.0) {
        y = rtNaN;
      } else if (d2 > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d3 == 0.0) {
      y = 1.0;
    } else if (d3 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

/*
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_remd_snf(double u0, double u1)
{
  double y;
  double b_u1;
  double tr;
  if (!((!rtIsNaN(u0)) && (!rtIsInf(u0)) && ((!rtIsNaN(u1)) && (!rtIsInf(u1)))))
  {
    y = rtNaN;
  } else {
    if (u1 < 0.0) {
      b_u1 = ceil(u1);
    } else {
      b_u1 = floor(u1);
    }

    if ((u1 != 0.0) && (u1 != b_u1)) {
      tr = u0 / u1;
      if (fabs(tr - rt_roundd_snf(tr)) <= DBL_EPSILON * fabs(tr)) {
        y = 0.0;
      } else {
        y = fmod(u0, u1);
      }
    } else {
      y = fmod(u0, u1);
    }
  }

  return y;
}

/*
 * Arguments    : double u
 * Return Type  : double
 */
static double rt_roundd_snf(double u)
{
  double y;
  if (fabs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = floor(u + 0.5);
    } else if (u > -0.5) {
      y = u * 0.0;
    } else {
      y = ceil(u - 0.5);
    }
  } else {
    y = u;
  }

  return y;
}

/*
 * Arguments    : double x[2048]
 *                double y[2048]
 * Return Type  : void
 */
static void sinc(double x[2048], double y[2048])
{
  int idx;
  short ii_data[2048];
  short ii_size[2];
  int i20;
  int ii;
  static const short iv0[2] = { 1, 2048 };

  boolean_T exitg1;
  boolean_T guard1 = false;
  short i_data[2048];
  idx = 0;
  for (i20 = 0; i20 < 2; i20++) {
    ii_size[i20] = iv0[i20];
  }

  ii = 1;
  exitg1 = false;
  while ((!exitg1) && (ii < 2049)) {
    guard1 = false;
    if (x[ii - 1] == 0.0) {
      idx++;
      ii_data[idx - 1] = (short)ii;
      if (idx >= 2048) {
        exitg1 = true;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      ii++;
    }
  }

  if (1 > idx) {
    idx = 0;
  }

  ii = ii_size[0] * idx;
  for (i20 = 0; i20 < ii; i20++) {
    i_data[i20] = ii_data[i20];
  }

  for (i20 = 0; i20 < idx; i20++) {
    ii_data[i20] = i_data[i20];
  }

  for (i20 = 0; i20 < idx; i20++) {
    x[ii_data[i20] - 1] = 1.0;
  }

  for (i20 = 0; i20 < 2048; i20++) {
    y[i20] = 3.1415926535897931 * x[i20];
  }

  b_sin(y);
  for (i20 = 0; i20 < 2048; i20++) {
    y[i20] /= 3.1415926535897931 * x[i20];
  }

  for (i20 = 0; i20 < idx; i20++) {
    ii_data[i20] = i_data[i20];
  }

  for (i20 = 0; i20 < idx; i20++) {
    y[ii_data[i20] - 1] = 1.0;
  }
}

/*
 * Arguments    : const double x[2048]
 * Return Type  : double
 */
static double sum(const double x[2048])
{
  double y;
  int k;
  y = x[0];
  for (k = 0; k < 2047; k++) {
    y += x[k + 1];
  }

  return y;
}

/*
 * Arguments    : const double o[15]
 *                double u[15]
 * Return Type  : void
 */
static void us(const double o[15], double u[15])
{
  int ix;
  int iy;
  int k;
  memset(&u[0], 0, 15U * sizeof(double));
  ix = 0;
  iy = 0;
  for (k = 0; k < 15; k++) {
    u[iy] = o[ix];
    ix++;
    iy++;
  }
}

/*
 * Arguments    : const emxArray_creal_T *x
 *                emxArray_creal_T *c
 * Return Type  : void
 */
static void vector_poly(const emxArray_creal_T *x, emxArray_creal_T *c)
{
  int n;
  int unnamed_idx_1;
  int k;
  double x_re;
  double x_im;
  double c_re;
  double c_im;
  n = x->size[0];
  unnamed_idx_1 = x->size[0] + 1;
  k = c->size[0] * c->size[1];
  c->size[0] = 1;
  c->size[1] = unnamed_idx_1;
  emxEnsureCapacity((emxArray__common *)c, k, (int)sizeof(creal_T));
  c->data[0].re = 1.0;
  c->data[0].im = 0.0;
  for (unnamed_idx_1 = 0; unnamed_idx_1 + 1 <= n; unnamed_idx_1++) {
    x_re = -x->data[unnamed_idx_1].re;
    x_im = -x->data[unnamed_idx_1].im;
    c_re = c->data[unnamed_idx_1].re;
    c_im = c->data[unnamed_idx_1].im;
    c->data[unnamed_idx_1 + 1].re = x_re * c_re - x_im * c_im;
    c->data[unnamed_idx_1 + 1].im = x_re * c_im + x_im * c_re;
    for (k = unnamed_idx_1; k + 1 > 1; k--) {
      x_re = x->data[unnamed_idx_1].re * c->data[k - 1].re - x->
        data[unnamed_idx_1].im * c->data[k - 1].im;
      x_im = x->data[unnamed_idx_1].re * c->data[k - 1].im + x->
        data[unnamed_idx_1].im * c->data[k - 1].re;
      c->data[k].re -= x_re;
      c->data[k].im -= x_im;
    }
  }
}

/*
 * Arguments    : double x1
 *                double x2
 *                double x3
 * Return Type  : double
 */
static double xdlapy3(double x1, double x2, double x3)
{
  double y;
  double a;
  double b;
  double c;
  a = fabs(x1);
  b = fabs(x2);
  c = fabs(x3);
  if ((a >= b) || rtIsNaN(b)) {
    y = a;
  } else {
    y = b;
  }

  if (c > y) {
    y = c;
  }

  if ((y > 0.0) && (!rtIsInf(y))) {
    a /= y;
    b /= y;
    c /= y;
    y *= sqrt((a * a + c * c) + b * b);
  } else {
    y = (a + b) + c;
  }

  return y;
}

/*
 * Arguments    : emxArray_creal_T *a
 * Return Type  : void
 */
static void xgehrd(emxArray_creal_T *a)
{
  int n;
  int ntau;
  emxArray_creal_T *tau;
  emxArray_creal_T *work;
  int i51;
  int i;
  int im1n;
  int in;
  int b_i;
  int c;
  creal_T alpha1;
  double c_re;
  double c_im;
  double xnorm;
  double beta1;
  int jy;
  int knt;
  boolean_T b_tau;
  double ai;
  int lastv;
  int lastc;
  int k;
  boolean_T exitg6;
  creal_T b_a;
  boolean_T exitg5;
  int ix;
  double temp_im;
  int exitg4;
  boolean_T exitg3;
  creal_T b_alpha1;
  boolean_T exitg2;
  int exitg1;
  n = a->size[0];
  if (a->size[0] < 1) {
    ntau = 0;
  } else {
    ntau = a->size[0] - 1;
  }

  emxInit_creal_T1(&tau, 1);
  emxInit_creal_T1(&work, 1);
  i51 = tau->size[0];
  tau->size[0] = ntau;
  emxEnsureCapacity((emxArray__common *)tau, i51, (int)sizeof(creal_T));
  ntau = a->size[0];
  i51 = work->size[0];
  work->size[0] = ntau;
  emxEnsureCapacity((emxArray__common *)work, i51, (int)sizeof(creal_T));
  for (i51 = 0; i51 < ntau; i51++) {
    work->data[i51].re = 0.0;
    work->data[i51].im = 0.0;
  }

  for (i = 0; i + 1 < n; i++) {
    im1n = i * n + 2;
    in = (i + 1) * n;
    if (i + 3 <= n) {
      b_i = i + 3;
    } else {
      b_i = n;
    }

    ntau = b_i + i * n;
    c = (n - i) - 2;
    alpha1 = a->data[(i + a->size[0] * i) + 1];
    c_re = 0.0;
    c_im = 0.0;
    if (!(c + 1 <= 0)) {
      xnorm = xnrm2(c, a, ntau);
      if ((xnorm != 0.0) || (a->data[(i + a->size[0] * i) + 1].im != 0.0)) {
        beta1 = xdlapy3(a->data[(i + a->size[0] * i) + 1].re, a->data[(i +
          a->size[0] * i) + 1].im, xnorm);
        if (a->data[(i + a->size[0] * i) + 1].re >= 0.0) {
          beta1 = -beta1;
        }

        if (fabs(beta1) < 1.0020841800044864E-292) {
          knt = 0;
          do {
            knt++;
            i51 = (ntau + c) - 1;
            for (k = ntau; k <= i51; k++) {
              xnorm = a->data[k - 1].re;
              ai = a->data[k - 1].im;
              a->data[k - 1].re = 9.9792015476736E+291 * xnorm - 0.0 * ai;
              a->data[k - 1].im = 9.9792015476736E+291 * ai + 0.0 * xnorm;
            }

            beta1 *= 9.9792015476736E+291;
            alpha1.re *= 9.9792015476736E+291;
            alpha1.im *= 9.9792015476736E+291;
          } while (!(fabs(beta1) >= 1.0020841800044864E-292));

          xnorm = xnrm2(c, a, ntau);
          beta1 = xdlapy3(alpha1.re, alpha1.im, xnorm);
          if (alpha1.re >= 0.0) {
            beta1 = -beta1;
          }

          xnorm = beta1 - alpha1.re;
          if (0.0 - alpha1.im == 0.0) {
            c_re = xnorm / beta1;
          } else if (xnorm == 0.0) {
            c_im = (0.0 - alpha1.im) / beta1;
          } else {
            c_re = xnorm / beta1;
            c_im = (0.0 - alpha1.im) / beta1;
          }

          b_alpha1.re = alpha1.re - beta1;
          b_alpha1.im = alpha1.im;
          xscal(c, recip(b_alpha1), a, ntau);
          for (k = 1; k <= knt; k++) {
            beta1 *= 1.0020841800044864E-292;
          }

          alpha1.re = beta1;
          alpha1.im = 0.0;
        } else {
          xnorm = beta1 - a->data[(i + a->size[0] * i) + 1].re;
          ai = 0.0 - a->data[(i + a->size[0] * i) + 1].im;
          if (ai == 0.0) {
            c_re = xnorm / beta1;
          } else if (xnorm == 0.0) {
            c_im = ai / beta1;
          } else {
            c_re = xnorm / beta1;
            c_im = ai / beta1;
          }

          b_a.re = a->data[(i + a->size[0] * i) + 1].re - beta1;
          b_a.im = a->data[(i + a->size[0] * i) + 1].im;
          xscal(c, recip(b_a), a, ntau);
          alpha1.re = beta1;
          alpha1.im = 0.0;
        }
      }
    }

    tau->data[i].re = c_re;
    tau->data[i].im = c_im;
    a->data[(i + a->size[0] * i) + 1].re = 1.0;
    a->data[(i + a->size[0] * i) + 1].im = 0.0;
    c = (n - i) - 3;
    jy = (i + im1n) - 1;
    b_tau = ((tau->data[i].re != 0.0) || (tau->data[i].im != 0.0));
    if (b_tau) {
      lastv = c + 2;
      ntau = jy + c;
      exitg6 = false;
      while ((!exitg6) && (lastv > 0)) {
        b_tau = ((a->data[ntau + 1].re == 0.0) && (a->data[ntau + 1].im == 0.0));
        if (b_tau) {
          lastv--;
          ntau--;
        } else {
          exitg6 = true;
        }
      }

      lastc = n;
      exitg5 = false;
      while ((!exitg5) && (lastc > 0)) {
        ntau = in + lastc;
        c = ntau;
        do {
          exitg4 = 0;
          if ((n > 0) && (c <= ntau + (lastv - 1) * n)) {
            b_tau = ((a->data[c - 1].re != 0.0) || (a->data[c - 1].im != 0.0));
            if (b_tau) {
              exitg4 = 1;
            } else {
              c += n;
            }
          } else {
            lastc--;
            exitg4 = 2;
          }
        } while (exitg4 == 0);

        if (exitg4 == 1) {
          exitg5 = true;
        }
      }
    } else {
      lastv = 0;
      lastc = 0;
    }

    if (lastv > 0) {
      if (lastc != 0) {
        for (ntau = 1; ntau <= lastc; ntau++) {
          work->data[ntau - 1].re = 0.0;
          work->data[ntau - 1].im = 0.0;
        }

        ix = jy;
        i51 = (in + n * (lastv - 1)) + 1;
        knt = in + 1;
        while ((n > 0) && (knt <= i51)) {
          c_re = a->data[ix].re - 0.0 * a->data[ix].im;
          c_im = a->data[ix].im + 0.0 * a->data[ix].re;
          ntau = 0;
          k = (knt + lastc) - 1;
          for (c = knt; c <= k; c++) {
            xnorm = a->data[c - 1].re * c_re - a->data[c - 1].im * c_im;
            ai = a->data[c - 1].re * c_im + a->data[c - 1].im * c_re;
            work->data[ntau].re += xnorm;
            work->data[ntau].im += ai;
            ntau++;
          }

          ix++;
          knt += n;
        }
      }

      c_re = -tau->data[i].re;
      c_im = -tau->data[i].im;
      if (!((c_re == 0.0) && (c_im == 0.0))) {
        ntau = in;
        for (knt = 1; knt <= lastv; knt++) {
          b_tau = ((a->data[jy].re != 0.0) || (a->data[jy].im != 0.0));
          if (b_tau) {
            beta1 = a->data[jy].re * c_re + a->data[jy].im * c_im;
            temp_im = a->data[jy].re * c_im - a->data[jy].im * c_re;
            ix = 0;
            i51 = lastc + ntau;
            for (k = ntau; k + 1 <= i51; k++) {
              xnorm = work->data[ix].re * beta1 - work->data[ix].im * temp_im;
              ai = work->data[ix].re * temp_im + work->data[ix].im * beta1;
              a->data[k].re += xnorm;
              a->data[k].im += ai;
              ix++;
            }
          }

          jy++;
          ntau += n;
        }
      }
    }

    c = (n - i) - 3;
    im1n = (i + im1n) - 1;
    jy = (i + in) + 2;
    beta1 = tau->data[i].re;
    temp_im = -tau->data[i].im;
    if ((beta1 != 0.0) || (temp_im != 0.0)) {
      lastv = c + 2;
      ntau = im1n + c;
      exitg3 = false;
      while ((!exitg3) && (lastv > 0)) {
        b_tau = ((a->data[ntau + 1].re == 0.0) && (a->data[ntau + 1].im == 0.0));
        if (b_tau) {
          lastv--;
          ntau--;
        } else {
          exitg3 = true;
        }
      }

      lastc = (n - i) - 1;
      exitg2 = false;
      while ((!exitg2) && (lastc > 0)) {
        ntau = jy + (lastc - 1) * n;
        c = ntau;
        do {
          exitg1 = 0;
          if (c <= (ntau + lastv) - 1) {
            b_tau = ((a->data[c - 1].re != 0.0) || (a->data[c - 1].im != 0.0));
            if (b_tau) {
              exitg1 = 1;
            } else {
              c++;
            }
          } else {
            lastc--;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    } else {
      lastv = 0;
      lastc = 0;
    }

    if (lastv > 0) {
      if (lastc != 0) {
        for (ntau = 1; ntau <= lastc; ntau++) {
          work->data[ntau - 1].re = 0.0;
          work->data[ntau - 1].im = 0.0;
        }

        ntau = 0;
        i51 = jy + n * (lastc - 1);
        knt = jy;
        while ((n > 0) && (knt <= i51)) {
          ix = im1n;
          c_re = 0.0;
          c_im = 0.0;
          k = (knt + lastv) - 1;
          for (c = knt - 1; c + 1 <= k; c++) {
            c_re += a->data[c].re * a->data[ix].re + a->data[c].im * a->data[ix]
              .im;
            c_im += a->data[c].re * a->data[ix].im - a->data[c].im * a->data[ix]
              .re;
            ix++;
          }

          work->data[ntau].re += c_re - 0.0 * c_im;
          work->data[ntau].im += c_im + 0.0 * c_re;
          ntau++;
          knt += n;
        }
      }

      c_re = -beta1;
      c_im = -temp_im;
      if (!((c_re == 0.0) && (c_im == 0.0))) {
        ntau = jy - 1;
        jy = 0;
        for (knt = 1; knt <= lastc; knt++) {
          b_tau = ((work->data[jy].re != 0.0) || (work->data[jy].im != 0.0));
          if (b_tau) {
            beta1 = work->data[jy].re * c_re + work->data[jy].im * c_im;
            temp_im = work->data[jy].re * c_im - work->data[jy].im * c_re;
            ix = im1n;
            i51 = lastv + ntau;
            for (k = ntau; k + 1 <= i51; k++) {
              xnorm = a->data[ix].re * beta1 - a->data[ix].im * temp_im;
              ai = a->data[ix].re * temp_im + a->data[ix].im * beta1;
              a->data[k].re += xnorm;
              a->data[k].im += ai;
              ix++;
            }
          }

          jy++;
          ntau += n;
        }
      }
    }

    a->data[(i + a->size[0] * i) + 1] = alpha1;
  }

  emxFree_creal_T(&work);
  emxFree_creal_T(&tau);
}

/*
 * Arguments    : int n
 *                const emxArray_creal_T *x
 *                int ix0
 * Return Type  : double
 */
static double xnrm2(int n, const emxArray_creal_T *x, int ix0)
{
  double y;
  double scale;
  int kend;
  int k;
  double absxk;
  double t;
  y = 0.0;
  if (!(n < 1)) {
    if (n == 1) {
      y = rt_hypotd_snf(x->data[ix0 - 1].re, x->data[ix0 - 1].im);
    } else {
      scale = 2.2250738585072014E-308;
      kend = (ix0 + n) - 1;
      for (k = ix0; k <= kend; k++) {
        absxk = fabs(x->data[k - 1].re);
        if (absxk > scale) {
          t = scale / absxk;
          y = 1.0 + y * t * t;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }

        absxk = fabs(x->data[k - 1].im);
        if (absxk > scale) {
          t = scale / absxk;
          y = 1.0 + y * t * t;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * sqrt(y);
    }
  }

  return y;
}

/*
 * Arguments    : int n
 *                const creal_T a
 *                emxArray_creal_T *x
 *                int ix0
 * Return Type  : void
 */
static void xscal(int n, const creal_T a, emxArray_creal_T *x, int ix0)
{
  double a_re;
  double a_im;
  int i52;
  int k;
  double x_re;
  double x_im;
  a_re = a.re;
  a_im = a.im;
  i52 = (ix0 + n) - 1;
  for (k = ix0; k <= i52; k++) {
    x_re = x->data[k - 1].re;
    x_im = x->data[k - 1].im;
    x->data[k - 1].re = a_re * x_re - a_im * x_im;
    x->data[k - 1].im = a_re * x_im + a_im * x_re;
  }
}

/*
 * Arguments    : const emxArray_creal_T *A
 *                int *info
 *                emxArray_creal_T *alpha1
 *                emxArray_creal_T *beta1
 * Return Type  : void
 */
static void xzgeev(const emxArray_creal_T *A, int *info, emxArray_creal_T
                   *alpha1, emxArray_creal_T *beta1)
{
  emxArray_creal_T *At;
  int ii;
  int n;
  int nzcount;
  double anrm;
  boolean_T exitg7;
  double absxk;
  boolean_T ilascl;
  double anrmto;
  int ilo;
  double ctoc;
  int ihi;
  boolean_T notdone;
  emxArray_creal_T *b_A;
  double cfrom1;
  int exitg2;
  double cto1;
  int i;
  int j;
  double mul;
  boolean_T exitg5;
  creal_T b_At;
  creal_T c_At;
  double c;
  creal_T atmp;
  boolean_T exitg6;
  int exitg1;
  boolean_T d_At;
  boolean_T guard2 = false;
  double stemp_re;
  boolean_T exitg3;
  boolean_T exitg4;
  boolean_T guard1 = false;
  emxInit_creal_T(&At, 2);
  ii = At->size[0] * At->size[1];
  At->size[0] = A->size[0];
  At->size[1] = A->size[1];
  emxEnsureCapacity((emxArray__common *)At, ii, (int)sizeof(creal_T));
  n = A->size[0] * A->size[1];
  for (ii = 0; ii < n; ii++) {
    At->data[ii] = A->data[ii];
  }

  nzcount = 0;
  ii = alpha1->size[0];
  alpha1->size[0] = A->size[0];
  emxEnsureCapacity((emxArray__common *)alpha1, ii, (int)sizeof(creal_T));
  n = A->size[0];
  for (ii = 0; ii < n; ii++) {
    alpha1->data[ii].re = 0.0;
    alpha1->data[ii].im = 0.0;
  }

  ii = beta1->size[0];
  beta1->size[0] = A->size[0];
  emxEnsureCapacity((emxArray__common *)beta1, ii, (int)sizeof(creal_T));
  n = A->size[0];
  for (ii = 0; ii < n; ii++) {
    beta1->data[ii].re = 0.0;
    beta1->data[ii].im = 0.0;
  }

  if (!((A->size[0] == 0) || (A->size[1] == 0))) {
    anrm = 0.0;
    ii = 0;
    exitg7 = false;
    while ((!exitg7) && (ii <= A->size[0] * A->size[1] - 1)) {
      absxk = rt_hypotd_snf(A->data[ii].re, A->data[ii].im);
      if (rtIsNaN(absxk)) {
        anrm = rtNaN;
        exitg7 = true;
      } else {
        if (absxk > anrm) {
          anrm = absxk;
        }

        ii++;
      }
    }

    if (!((!rtIsInf(anrm)) && (!rtIsNaN(anrm)))) {
      ii = alpha1->size[0];
      alpha1->size[0] = A->size[0];
      emxEnsureCapacity((emxArray__common *)alpha1, ii, (int)sizeof(creal_T));
      n = A->size[0];
      for (ii = 0; ii < n; ii++) {
        alpha1->data[ii].re = rtNaN;
        alpha1->data[ii].im = 0.0;
      }

      ii = beta1->size[0];
      beta1->size[0] = A->size[0];
      emxEnsureCapacity((emxArray__common *)beta1, ii, (int)sizeof(creal_T));
      n = A->size[0];
      for (ii = 0; ii < n; ii++) {
        beta1->data[ii].re = rtNaN;
        beta1->data[ii].im = 0.0;
      }
    } else {
      ilascl = false;
      anrmto = anrm;
      if ((anrm > 0.0) && (anrm < 6.7178761075670888E-139)) {
        anrmto = 6.7178761075670888E-139;
        ilascl = true;
      } else {
        if (anrm > 1.4885657073574029E+138) {
          anrmto = 1.4885657073574029E+138;
          ilascl = true;
        }
      }

      if (ilascl) {
        absxk = anrm;
        ctoc = anrmto;
        notdone = true;
        while (notdone) {
          cfrom1 = absxk * 2.0041683600089728E-292;
          cto1 = ctoc / 4.9896007738368E+291;
          if ((cfrom1 > ctoc) && (ctoc != 0.0)) {
            mul = 2.0041683600089728E-292;
            absxk = cfrom1;
          } else if (cto1 > absxk) {
            mul = 4.9896007738368E+291;
            ctoc = cto1;
          } else {
            mul = ctoc / absxk;
            notdone = false;
          }

          ii = At->size[0] * At->size[1];
          emxEnsureCapacity((emxArray__common *)At, ii, (int)sizeof(creal_T));
          ii = At->size[0];
          nzcount = At->size[1];
          n = ii * nzcount;
          for (ii = 0; ii < n; ii++) {
            At->data[ii].re *= mul;
            At->data[ii].im *= mul;
          }
        }
      }

      ilo = 0;
      ihi = At->size[0];
      if (At->size[0] <= 1) {
        ihi = 1;
      } else {
        emxInit_creal_T(&b_A, 2);
        do {
          exitg2 = 0;
          i = 0;
          j = 0;
          notdone = false;
          ii = ihi;
          exitg5 = false;
          while ((!exitg5) && (ii > 0)) {
            nzcount = 0;
            i = ii;
            j = ihi;
            n = 1;
            exitg6 = false;
            while ((!exitg6) && (n <= ihi)) {
              d_At = ((At->data[(ii + At->size[0] * (n - 1)) - 1].re != 0.0) ||
                      (At->data[(ii + At->size[0] * (n - 1)) - 1].im != 0.0));
              guard2 = false;
              if (d_At || (ii == n)) {
                if (nzcount == 0) {
                  j = n;
                  nzcount = 1;
                  guard2 = true;
                } else {
                  nzcount = 2;
                  exitg6 = true;
                }
              } else {
                guard2 = true;
              }

              if (guard2) {
                n++;
              }
            }

            if (nzcount < 2) {
              notdone = true;
              exitg5 = true;
            } else {
              ii--;
            }
          }

          if (!notdone) {
            exitg2 = 2;
          } else {
            ii = b_A->size[0] * b_A->size[1];
            b_A->size[0] = At->size[0];
            b_A->size[1] = At->size[1];
            emxEnsureCapacity((emxArray__common *)b_A, ii, (int)sizeof(creal_T));
            n = At->size[0] * At->size[1];
            for (ii = 0; ii < n; ii++) {
              b_A->data[ii] = At->data[ii];
            }

            if (i != ihi) {
              for (ii = 0; ii + 1 <= At->size[0]; ii++) {
                atmp = b_A->data[(i + b_A->size[0] * ii) - 1];
                b_A->data[(i + b_A->size[0] * ii) - 1] = b_A->data[(ihi +
                  b_A->size[0] * ii) - 1];
                b_A->data[(ihi + b_A->size[0] * ii) - 1] = atmp;
              }
            }

            if (j != ihi) {
              for (ii = 0; ii + 1 <= ihi; ii++) {
                atmp = b_A->data[ii + b_A->size[0] * (j - 1)];
                b_A->data[ii + b_A->size[0] * (j - 1)] = b_A->data[ii +
                  b_A->size[0] * (ihi - 1)];
                b_A->data[ii + b_A->size[0] * (ihi - 1)] = atmp;
              }
            }

            ii = At->size[0] * At->size[1];
            At->size[0] = b_A->size[0];
            At->size[1] = b_A->size[1];
            emxEnsureCapacity((emxArray__common *)At, ii, (int)sizeof(creal_T));
            n = b_A->size[0] * b_A->size[1];
            for (ii = 0; ii < n; ii++) {
              At->data[ii] = b_A->data[ii];
            }

            ihi--;
            if (ihi == 1) {
              exitg2 = 1;
            }
          }
        } while (exitg2 == 0);

        if (exitg2 == 1) {
        } else {
          do {
            exitg1 = 0;
            i = 0;
            j = 0;
            notdone = false;
            n = ilo + 1;
            exitg3 = false;
            while ((!exitg3) && (n <= ihi)) {
              nzcount = 0;
              i = ihi;
              j = n;
              ii = ilo + 1;
              exitg4 = false;
              while ((!exitg4) && (ii <= ihi)) {
                d_At = ((At->data[(ii + At->size[0] * (n - 1)) - 1].re != 0.0) ||
                        (At->data[(ii + At->size[0] * (n - 1)) - 1].im != 0.0));
                guard1 = false;
                if (d_At || (ii == n)) {
                  if (nzcount == 0) {
                    i = ii;
                    nzcount = 1;
                    guard1 = true;
                  } else {
                    nzcount = 2;
                    exitg4 = true;
                  }
                } else {
                  guard1 = true;
                }

                if (guard1) {
                  ii++;
                }
              }

              if (nzcount < 2) {
                notdone = true;
                exitg3 = true;
              } else {
                n++;
              }
            }

            if (!notdone) {
              exitg1 = 1;
            } else {
              ii = b_A->size[0] * b_A->size[1];
              b_A->size[0] = At->size[0];
              b_A->size[1] = At->size[1];
              emxEnsureCapacity((emxArray__common *)b_A, ii, (int)sizeof(creal_T));
              n = At->size[0] * At->size[1];
              for (ii = 0; ii < n; ii++) {
                b_A->data[ii] = At->data[ii];
              }

              if (i != ilo + 1) {
                for (ii = ilo; ii + 1 <= At->size[0]; ii++) {
                  atmp = b_A->data[(i + b_A->size[0] * ii) - 1];
                  b_A->data[(i + b_A->size[0] * ii) - 1] = b_A->data[ilo +
                    b_A->size[0] * ii];
                  b_A->data[ilo + b_A->size[0] * ii] = atmp;
                }
              }

              if (j != ilo + 1) {
                for (ii = 0; ii + 1 <= ihi; ii++) {
                  atmp = b_A->data[ii + b_A->size[0] * (j - 1)];
                  b_A->data[ii + b_A->size[0] * (j - 1)] = b_A->data[ii +
                    b_A->size[0] * ilo];
                  b_A->data[ii + b_A->size[0] * ilo] = atmp;
                }
              }

              ii = At->size[0] * At->size[1];
              At->size[0] = b_A->size[0];
              At->size[1] = b_A->size[1];
              emxEnsureCapacity((emxArray__common *)At, ii, (int)sizeof(creal_T));
              n = b_A->size[0] * b_A->size[1];
              for (ii = 0; ii < n; ii++) {
                At->data[ii] = b_A->data[ii];
              }

              ilo++;
              if (ilo + 1 == ihi) {
                exitg1 = 1;
              }
            }
          } while (exitg1 == 0);
        }

        emxFree_creal_T(&b_A);
      }

      n = At->size[0];
      if ((!(At->size[0] <= 1)) && (!(ihi < ilo + 3))) {
        for (ii = ilo; ii + 1 < ihi - 1; ii++) {
          for (nzcount = ihi - 1; nzcount + 1 > ii + 2; nzcount--) {
            b_At = At->data[(nzcount + At->size[0] * ii) - 1];
            c_At = At->data[nzcount + At->size[0] * ii];
            xzlartg(b_At, c_At, &c, &atmp, &At->data[(nzcount + At->size[0] * ii)
                    - 1]);
            At->data[nzcount + At->size[0] * ii].re = 0.0;
            At->data[nzcount + At->size[0] * ii].im = 0.0;
            for (j = ii + 1; j + 1 <= n; j++) {
              absxk = atmp.re * At->data[nzcount + At->size[0] * j].re - atmp.im
                * At->data[nzcount + At->size[0] * j].im;
              ctoc = atmp.re * At->data[nzcount + At->size[0] * j].im + atmp.im *
                At->data[nzcount + At->size[0] * j].re;
              stemp_re = c * At->data[(nzcount + At->size[0] * j) - 1].re +
                absxk;
              absxk = c * At->data[(nzcount + At->size[0] * j) - 1].im + ctoc;
              ctoc = At->data[(nzcount + At->size[0] * j) - 1].re;
              cfrom1 = At->data[(nzcount + At->size[0] * j) - 1].im;
              cto1 = At->data[(nzcount + At->size[0] * j) - 1].im;
              mul = At->data[(nzcount + At->size[0] * j) - 1].re;
              At->data[nzcount + At->size[0] * j].re = c * At->data[nzcount +
                At->size[0] * j].re - (atmp.re * ctoc + atmp.im * cfrom1);
              At->data[nzcount + At->size[0] * j].im = c * At->data[nzcount +
                At->size[0] * j].im - (atmp.re * cto1 - atmp.im * mul);
              At->data[(nzcount + At->size[0] * j) - 1].re = stemp_re;
              At->data[(nzcount + At->size[0] * j) - 1].im = absxk;
            }

            atmp.re = -atmp.re;
            atmp.im = -atmp.im;
            for (i = 0; i + 1 <= ihi; i++) {
              absxk = atmp.re * At->data[i + At->size[0] * (nzcount - 1)].re -
                atmp.im * At->data[i + At->size[0] * (nzcount - 1)].im;
              ctoc = atmp.re * At->data[i + At->size[0] * (nzcount - 1)].im +
                atmp.im * At->data[i + At->size[0] * (nzcount - 1)].re;
              stemp_re = c * At->data[i + At->size[0] * nzcount].re + absxk;
              absxk = c * At->data[i + At->size[0] * nzcount].im + ctoc;
              ctoc = At->data[i + At->size[0] * nzcount].re;
              cfrom1 = At->data[i + At->size[0] * nzcount].im;
              cto1 = At->data[i + At->size[0] * nzcount].im;
              mul = At->data[i + At->size[0] * nzcount].re;
              At->data[i + At->size[0] * (nzcount - 1)].re = c * At->data[i +
                At->size[0] * (nzcount - 1)].re - (atmp.re * ctoc + atmp.im *
                cfrom1);
              At->data[i + At->size[0] * (nzcount - 1)].im = c * At->data[i +
                At->size[0] * (nzcount - 1)].im - (atmp.re * cto1 - atmp.im *
                mul);
              At->data[i + At->size[0] * nzcount].re = stemp_re;
              At->data[i + At->size[0] * nzcount].im = absxk;
            }
          }
        }
      }

      xzhgeqz(At, ilo + 1, ihi, &nzcount, alpha1, beta1);
      if ((nzcount == 0) && ilascl) {
        notdone = true;
        while (notdone) {
          cfrom1 = anrmto * 2.0041683600089728E-292;
          cto1 = anrm / 4.9896007738368E+291;
          if ((cfrom1 > anrm) && (anrm != 0.0)) {
            mul = 2.0041683600089728E-292;
            anrmto = cfrom1;
          } else if (cto1 > anrmto) {
            mul = 4.9896007738368E+291;
            anrm = cto1;
          } else {
            mul = anrm / anrmto;
            notdone = false;
          }

          ii = alpha1->size[0];
          emxEnsureCapacity((emxArray__common *)alpha1, ii, (int)sizeof(creal_T));
          n = alpha1->size[0];
          for (ii = 0; ii < n; ii++) {
            alpha1->data[ii].re *= mul;
            alpha1->data[ii].im *= mul;
          }
        }
      }
    }
  }

  emxFree_creal_T(&At);
  *info = nzcount;
}

/*
 * Arguments    : const emxArray_creal_T *A
 *                int ilo
 *                int ihi
 *                int *info
 *                emxArray_creal_T *alpha1
 *                emxArray_creal_T *beta1
 * Return Type  : void
 */
static void xzhgeqz(const emxArray_creal_T *A, int ilo, int ihi, int *info,
                    emxArray_creal_T *alpha1, emxArray_creal_T *beta1)
{
  emxArray_creal_T *b_A;
  int jm1;
  int jp1;
  double eshift_re;
  double eshift_im;
  creal_T ctemp;
  double anorm;
  double scale;
  double reAij;
  double sumsq;
  double b_atol;
  boolean_T firstNonZero;
  int j;
  double ascale;
  double bscale;
  int i;
  boolean_T failed;
  double imAij;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  int ifirst;
  int istart;
  double temp2;
  int ilast;
  int ilastm1;
  int ifrstm;
  int ilastm;
  int iiter;
  boolean_T goto60;
  boolean_T goto70;
  boolean_T goto90;
  int jiter;
  int exitg1;
  boolean_T exitg3;
  boolean_T b_guard1 = false;
  creal_T b_ascale;
  creal_T shift;
  creal_T c_A;
  creal_T d_A;
  double ad22_re;
  double ad22_im;
  boolean_T exitg2;
  double t1_im;
  emxInit_creal_T(&b_A, 2);
  jm1 = b_A->size[0] * b_A->size[1];
  b_A->size[0] = A->size[0];
  b_A->size[1] = A->size[1];
  emxEnsureCapacity((emxArray__common *)b_A, jm1, (int)sizeof(creal_T));
  jp1 = A->size[0] * A->size[1];
  for (jm1 = 0; jm1 < jp1; jm1++) {
    b_A->data[jm1] = A->data[jm1];
  }

  *info = 0;
  if ((A->size[0] == 1) && (A->size[1] == 1)) {
    ihi = 1;
  }

  jm1 = alpha1->size[0];
  alpha1->size[0] = A->size[0];
  emxEnsureCapacity((emxArray__common *)alpha1, jm1, (int)sizeof(creal_T));
  jp1 = A->size[0];
  for (jm1 = 0; jm1 < jp1; jm1++) {
    alpha1->data[jm1].re = 0.0;
    alpha1->data[jm1].im = 0.0;
  }

  jm1 = beta1->size[0];
  beta1->size[0] = A->size[0];
  emxEnsureCapacity((emxArray__common *)beta1, jm1, (int)sizeof(creal_T));
  jp1 = A->size[0];
  for (jm1 = 0; jm1 < jp1; jm1++) {
    beta1->data[jm1].re = 1.0;
    beta1->data[jm1].im = 0.0;
  }

  eshift_re = 0.0;
  eshift_im = 0.0;
  ctemp.re = 0.0;
  ctemp.im = 0.0;
  anorm = 0.0;
  if (!(ilo > ihi)) {
    scale = 0.0;
    sumsq = 0.0;
    firstNonZero = true;
    for (j = ilo; j <= ihi; j++) {
      jm1 = j + 1;
      if (ihi < j + 1) {
        jm1 = ihi;
      }

      for (i = ilo; i <= jm1; i++) {
        reAij = A->data[(i + A->size[0] * (j - 1)) - 1].re;
        imAij = A->data[(i + A->size[0] * (j - 1)) - 1].im;
        if (reAij != 0.0) {
          anorm = fabs(reAij);
          if (firstNonZero) {
            sumsq = 1.0;
            scale = anorm;
            firstNonZero = false;
          } else if (scale < anorm) {
            temp2 = scale / anorm;
            sumsq = 1.0 + sumsq * temp2 * temp2;
            scale = anorm;
          } else {
            temp2 = anorm / scale;
            sumsq += temp2 * temp2;
          }
        }

        if (imAij != 0.0) {
          anorm = fabs(imAij);
          if (firstNonZero) {
            sumsq = 1.0;
            scale = anorm;
            firstNonZero = false;
          } else if (scale < anorm) {
            temp2 = scale / anorm;
            sumsq = 1.0 + sumsq * temp2 * temp2;
            scale = anorm;
          } else {
            temp2 = anorm / scale;
            sumsq += temp2 * temp2;
          }
        }
      }
    }

    anorm = scale * sqrt(sumsq);
  }

  reAij = 2.2204460492503131E-16 * anorm;
  b_atol = 2.2250738585072014E-308;
  if (reAij > 2.2250738585072014E-308) {
    b_atol = reAij;
  }

  reAij = 2.2250738585072014E-308;
  if (anorm > 2.2250738585072014E-308) {
    reAij = anorm;
  }

  ascale = 1.0 / reAij;
  bscale = 1.0 / sqrt(A->size[0]);
  failed = true;
  for (j = ihi; j + 1 <= A->size[0]; j++) {
    alpha1->data[j] = A->data[j + A->size[0] * j];
  }

  guard1 = false;
  guard2 = false;
  if (ihi >= ilo) {
    ifirst = ilo;
    istart = ilo;
    ilast = ihi - 1;
    ilastm1 = ihi - 2;
    ifrstm = ilo;
    ilastm = ihi;
    iiter = 0;
    goto60 = false;
    goto70 = false;
    goto90 = false;
    jiter = 1;
    do {
      exitg1 = 0;
      if (jiter <= 30 * ((ihi - ilo) + 1)) {
        if (ilast + 1 == ilo) {
          goto60 = true;
        } else if (fabs(b_A->data[ilast + b_A->size[0] * ilastm1].re) + fabs
                   (b_A->data[ilast + b_A->size[0] * ilastm1].im) <= b_atol) {
          b_A->data[ilast + b_A->size[0] * ilastm1].re = 0.0;
          b_A->data[ilast + b_A->size[0] * ilastm1].im = 0.0;
          goto60 = true;
        } else {
          j = ilastm1;
          exitg3 = false;
          while ((!exitg3) && (j + 1 >= ilo)) {
            if (j + 1 == ilo) {
              firstNonZero = true;
            } else if (fabs(b_A->data[j + b_A->size[0] * (j - 1)].re) + fabs
                       (b_A->data[j + b_A->size[0] * (j - 1)].im) <= b_atol) {
              b_A->data[j + b_A->size[0] * (j - 1)].re = 0.0;
              b_A->data[j + b_A->size[0] * (j - 1)].im = 0.0;
              firstNonZero = true;
            } else {
              firstNonZero = false;
            }

            if (firstNonZero) {
              ifirst = j + 1;
              goto70 = true;
              exitg3 = true;
            } else {
              j--;
            }
          }
        }

        if (goto60 || goto70) {
          firstNonZero = true;
        } else {
          firstNonZero = false;
        }

        if (!firstNonZero) {
          jp1 = alpha1->size[0];
          jm1 = alpha1->size[0];
          alpha1->size[0] = jp1;
          emxEnsureCapacity((emxArray__common *)alpha1, jm1, (int)sizeof(creal_T));
          for (jm1 = 0; jm1 < jp1; jm1++) {
            alpha1->data[jm1].re = rtNaN;
            alpha1->data[jm1].im = 0.0;
          }

          jp1 = beta1->size[0];
          jm1 = beta1->size[0];
          beta1->size[0] = jp1;
          emxEnsureCapacity((emxArray__common *)beta1, jm1, (int)sizeof(creal_T));
          for (jm1 = 0; jm1 < jp1; jm1++) {
            beta1->data[jm1].re = rtNaN;
            beta1->data[jm1].im = 0.0;
          }

          *info = 1;
          exitg1 = 1;
        } else {
          b_guard1 = false;
          if (goto60) {
            goto60 = false;
            alpha1->data[ilast] = b_A->data[ilast + b_A->size[0] * ilast];
            ilast = ilastm1;
            ilastm1--;
            if (ilast + 1 < ilo) {
              failed = false;
              guard2 = true;
              exitg1 = 1;
            } else {
              iiter = 0;
              eshift_re = 0.0;
              eshift_im = 0.0;
              ilastm = ilast + 1;
              if (ifrstm > ilast + 1) {
                ifrstm = ilo;
              }

              b_guard1 = true;
            }
          } else {
            if (goto70) {
              goto70 = false;
              iiter++;
              ifrstm = ifirst;
              if (iiter - iiter / 10 * 10 != 0) {
                anorm = ascale * b_A->data[ilastm1 + b_A->size[0] * ilastm1].re;
                reAij = ascale * b_A->data[ilastm1 + b_A->size[0] * ilastm1].im;
                if (reAij == 0.0) {
                  shift.re = anorm / bscale;
                  shift.im = 0.0;
                } else if (anorm == 0.0) {
                  shift.re = 0.0;
                  shift.im = reAij / bscale;
                } else {
                  shift.re = anorm / bscale;
                  shift.im = reAij / bscale;
                }

                anorm = ascale * b_A->data[ilast + b_A->size[0] * ilast].re;
                reAij = ascale * b_A->data[ilast + b_A->size[0] * ilast].im;
                if (reAij == 0.0) {
                  ad22_re = anorm / bscale;
                  ad22_im = 0.0;
                } else if (anorm == 0.0) {
                  ad22_re = 0.0;
                  ad22_im = reAij / bscale;
                } else {
                  ad22_re = anorm / bscale;
                  ad22_im = reAij / bscale;
                }

                temp2 = 0.5 * (shift.re + ad22_re);
                t1_im = 0.5 * (shift.im + ad22_im);
                anorm = ascale * b_A->data[ilastm1 + b_A->size[0] * ilast].re;
                reAij = ascale * b_A->data[ilastm1 + b_A->size[0] * ilast].im;
                if (reAij == 0.0) {
                  sumsq = anorm / bscale;
                  imAij = 0.0;
                } else if (anorm == 0.0) {
                  sumsq = 0.0;
                  imAij = reAij / bscale;
                } else {
                  sumsq = anorm / bscale;
                  imAij = reAij / bscale;
                }

                anorm = ascale * b_A->data[ilast + b_A->size[0] * ilastm1].re;
                reAij = ascale * b_A->data[ilast + b_A->size[0] * ilastm1].im;
                if (reAij == 0.0) {
                  scale = anorm / bscale;
                  anorm = 0.0;
                } else if (anorm == 0.0) {
                  scale = 0.0;
                  anorm = reAij / bscale;
                } else {
                  scale = anorm / bscale;
                  anorm = reAij / bscale;
                }

                reAij = shift.re * ad22_im + shift.im * ad22_re;
                shift.re = ((temp2 * temp2 - t1_im * t1_im) + (sumsq * scale -
                  imAij * anorm)) - (shift.re * ad22_re - shift.im * ad22_im);
                shift.im = ((temp2 * t1_im + t1_im * temp2) + (sumsq * anorm +
                  imAij * scale)) - reAij;
                b_sqrt(&shift);
                if ((temp2 - ad22_re) * shift.re + (t1_im - ad22_im) * shift.im <=
                    0.0) {
                  shift.re += temp2;
                  shift.im += t1_im;
                } else {
                  shift.re = temp2 - shift.re;
                  shift.im = t1_im - shift.im;
                }
              } else {
                anorm = ascale * b_A->data[ilast + b_A->size[0] * ilastm1].re;
                reAij = ascale * b_A->data[ilast + b_A->size[0] * ilastm1].im;
                if (reAij == 0.0) {
                  sumsq = anorm / bscale;
                  imAij = 0.0;
                } else if (anorm == 0.0) {
                  sumsq = 0.0;
                  imAij = reAij / bscale;
                } else {
                  sumsq = anorm / bscale;
                  imAij = reAij / bscale;
                }

                eshift_re += sumsq;
                eshift_im += imAij;
                shift.re = eshift_re;
                shift.im = eshift_im;
              }

              j = ilastm1;
              jp1 = ilastm1 + 1;
              exitg2 = false;
              while ((!exitg2) && (j + 1 > ifirst)) {
                istart = j + 1;
                ctemp.re = ascale * b_A->data[j + b_A->size[0] * j].re -
                  shift.re * bscale;
                ctemp.im = ascale * b_A->data[j + b_A->size[0] * j].im -
                  shift.im * bscale;
                anorm = fabs(ctemp.re) + fabs(ctemp.im);
                temp2 = ascale * (fabs(b_A->data[jp1 + b_A->size[0] * j].re) +
                                  fabs(b_A->data[jp1 + b_A->size[0] * j].im));
                reAij = anorm;
                if (temp2 > anorm) {
                  reAij = temp2;
                }

                if ((reAij < 1.0) && (reAij != 0.0)) {
                  anorm /= reAij;
                  temp2 /= reAij;
                }

                if ((fabs(b_A->data[j + b_A->size[0] * (j - 1)].re) + fabs
                     (b_A->data[j + b_A->size[0] * (j - 1)].im)) * temp2 <=
                    anorm * b_atol) {
                  goto90 = true;
                  exitg2 = true;
                } else {
                  jp1 = j;
                  j--;
                }
              }

              if (!goto90) {
                istart = ifirst;
                ctemp.re = ascale * b_A->data[(ifirst + b_A->size[0] * (ifirst -
                  1)) - 1].re - shift.re * bscale;
                ctemp.im = ascale * b_A->data[(ifirst + b_A->size[0] * (ifirst -
                  1)) - 1].im - shift.im * bscale;
                goto90 = true;
              }
            }

            if (goto90) {
              goto90 = false;
              b_ascale.re = ascale * b_A->data[istart + b_A->size[0] * (istart -
                1)].re;
              b_ascale.im = ascale * b_A->data[istart + b_A->size[0] * (istart -
                1)].im;
              b_xzlartg(ctemp, b_ascale, &imAij, &shift);
              j = istart;
              jm1 = istart - 2;
              while (j < ilast + 1) {
                if (j > istart) {
                  c_A = b_A->data[(j + b_A->size[0] * jm1) - 1];
                  d_A = b_A->data[j + b_A->size[0] * jm1];
                  xzlartg(c_A, d_A, &imAij, &shift, &b_A->data[(j + b_A->size[0]
                           * jm1) - 1]);
                  b_A->data[j + b_A->size[0] * jm1].re = 0.0;
                  b_A->data[j + b_A->size[0] * jm1].im = 0.0;
                }

                for (jp1 = j - 1; jp1 + 1 <= ilastm; jp1++) {
                  anorm = shift.re * b_A->data[j + b_A->size[0] * jp1].re -
                    shift.im * b_A->data[j + b_A->size[0] * jp1].im;
                  reAij = shift.re * b_A->data[j + b_A->size[0] * jp1].im +
                    shift.im * b_A->data[j + b_A->size[0] * jp1].re;
                  ad22_re = imAij * b_A->data[(j + b_A->size[0] * jp1) - 1].re +
                    anorm;
                  ad22_im = imAij * b_A->data[(j + b_A->size[0] * jp1) - 1].im +
                    reAij;
                  anorm = b_A->data[(j + b_A->size[0] * jp1) - 1].re;
                  reAij = b_A->data[(j + b_A->size[0] * jp1) - 1].im;
                  scale = b_A->data[(j + b_A->size[0] * jp1) - 1].im;
                  sumsq = b_A->data[(j + b_A->size[0] * jp1) - 1].re;
                  b_A->data[j + b_A->size[0] * jp1].re = imAij * b_A->data[j +
                    b_A->size[0] * jp1].re - (shift.re * anorm + shift.im *
                    reAij);
                  b_A->data[j + b_A->size[0] * jp1].im = imAij * b_A->data[j +
                    b_A->size[0] * jp1].im - (shift.re * scale - shift.im *
                    sumsq);
                  b_A->data[(j + b_A->size[0] * jp1) - 1].re = ad22_re;
                  b_A->data[(j + b_A->size[0] * jp1) - 1].im = ad22_im;
                }

                shift.re = -shift.re;
                shift.im = -shift.im;
                jp1 = j;
                if (ilast + 1 < j + 2) {
                  jp1 = ilast - 1;
                }

                for (i = ifrstm - 1; i + 1 <= jp1 + 2; i++) {
                  anorm = shift.re * b_A->data[i + b_A->size[0] * (j - 1)].re -
                    shift.im * b_A->data[i + b_A->size[0] * (j - 1)].im;
                  reAij = shift.re * b_A->data[i + b_A->size[0] * (j - 1)].im +
                    shift.im * b_A->data[i + b_A->size[0] * (j - 1)].re;
                  ad22_re = imAij * b_A->data[i + b_A->size[0] * j].re + anorm;
                  ad22_im = imAij * b_A->data[i + b_A->size[0] * j].im + reAij;
                  anorm = b_A->data[i + b_A->size[0] * j].re;
                  reAij = b_A->data[i + b_A->size[0] * j].im;
                  scale = b_A->data[i + b_A->size[0] * j].im;
                  sumsq = b_A->data[i + b_A->size[0] * j].re;
                  b_A->data[i + b_A->size[0] * (j - 1)].re = imAij * b_A->data[i
                    + b_A->size[0] * (j - 1)].re - (shift.re * anorm + shift.im *
                    reAij);
                  b_A->data[i + b_A->size[0] * (j - 1)].im = imAij * b_A->data[i
                    + b_A->size[0] * (j - 1)].im - (shift.re * scale - shift.im *
                    sumsq);
                  b_A->data[i + b_A->size[0] * j].re = ad22_re;
                  b_A->data[i + b_A->size[0] * j].im = ad22_im;
                }

                jm1 = j - 1;
                j++;
              }
            }

            b_guard1 = true;
          }

          if (b_guard1) {
            jiter++;
          }
        }
      } else {
        guard2 = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  } else {
    guard1 = true;
  }

  if (guard2) {
    if (failed) {
      *info = ilast + 1;
      for (jp1 = 0; jp1 + 1 <= ilast + 1; jp1++) {
        alpha1->data[jp1].re = rtNaN;
        alpha1->data[jp1].im = 0.0;
        beta1->data[jp1].re = rtNaN;
        beta1->data[jp1].im = 0.0;
      }
    } else {
      guard1 = true;
    }
  }

  if (guard1) {
    for (j = 0; j + 1 < ilo; j++) {
      alpha1->data[j] = b_A->data[j + b_A->size[0] * j];
    }
  }

  emxFree_creal_T(&b_A);
}

/*
 * Arguments    : const creal_T f
 *                const creal_T g
 *                double *cs
 *                creal_T *sn
 *                creal_T *r
 * Return Type  : void
 */
static void xzlartg(const creal_T f, const creal_T g, double *cs, creal_T *sn,
                    creal_T *r)
{
  double scale;
  double f2s;
  double x;
  double fs_re;
  double fs_im;
  double gs_re;
  double gs_im;
  int count;
  int rescaledir;
  boolean_T guard1 = false;
  double g2;
  double g2s;
  scale = fabs(f.re);
  f2s = fabs(f.im);
  if (f2s > scale) {
    scale = f2s;
  }

  x = fabs(g.re);
  f2s = fabs(g.im);
  if (f2s > x) {
    x = f2s;
  }

  if (x > scale) {
    scale = x;
  }

  fs_re = f.re;
  fs_im = f.im;
  gs_re = g.re;
  gs_im = g.im;
  count = 0;
  rescaledir = 0;
  guard1 = false;
  if (scale >= 7.4428285367870146E+137) {
    do {
      count++;
      fs_re *= 1.3435752215134178E-138;
      fs_im *= 1.3435752215134178E-138;
      gs_re *= 1.3435752215134178E-138;
      gs_im *= 1.3435752215134178E-138;
      scale *= 1.3435752215134178E-138;
    } while (!(scale < 7.4428285367870146E+137));

    rescaledir = 1;
    guard1 = true;
  } else if (scale <= 1.3435752215134178E-138) {
    if ((g.re == 0.0) && (g.im == 0.0)) {
      *cs = 1.0;
      sn->re = 0.0;
      sn->im = 0.0;
      *r = f;
    } else {
      do {
        count++;
        fs_re *= 7.4428285367870146E+137;
        fs_im *= 7.4428285367870146E+137;
        gs_re *= 7.4428285367870146E+137;
        gs_im *= 7.4428285367870146E+137;
        scale *= 7.4428285367870146E+137;
      } while (!(scale > 1.3435752215134178E-138));

      rescaledir = -1;
      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1) {
    scale = fs_re * fs_re + fs_im * fs_im;
    g2 = gs_re * gs_re + gs_im * gs_im;
    x = g2;
    if (1.0 > g2) {
      x = 1.0;
    }

    if (scale <= x * 2.0041683600089728E-292) {
      if ((f.re == 0.0) && (f.im == 0.0)) {
        *cs = 0.0;
        r->re = rt_hypotd_snf(g.re, g.im);
        r->im = 0.0;
        g2 = rt_hypotd_snf(gs_re, gs_im);
        sn->re = gs_re / g2;
        sn->im = -gs_im / g2;
      } else {
        g2s = sqrt(g2);
        *cs = rt_hypotd_snf(fs_re, fs_im) / g2s;
        x = fabs(f.re);
        f2s = fabs(f.im);
        if (f2s > x) {
          x = f2s;
        }

        if (x > 1.0) {
          g2 = rt_hypotd_snf(f.re, f.im);
          fs_re = f.re / g2;
          fs_im = f.im / g2;
        } else {
          scale = 7.4428285367870146E+137 * f.re;
          f2s = 7.4428285367870146E+137 * f.im;
          g2 = rt_hypotd_snf(scale, f2s);
          fs_re = scale / g2;
          fs_im = f2s / g2;
        }

        gs_re /= g2s;
        gs_im = -gs_im / g2s;
        sn->re = fs_re * gs_re - fs_im * gs_im;
        sn->im = fs_re * gs_im + fs_im * gs_re;
        r->re = *cs * f.re + (sn->re * g.re - sn->im * g.im);
        r->im = *cs * f.im + (sn->re * g.im + sn->im * g.re);
      }
    } else {
      f2s = sqrt(1.0 + g2 / scale);
      r->re = f2s * fs_re;
      r->im = f2s * fs_im;
      *cs = 1.0 / f2s;
      g2 += scale;
      scale = r->re / g2;
      f2s = r->im / g2;
      sn->re = scale * gs_re - f2s * -gs_im;
      sn->im = scale * -gs_im + f2s * gs_re;
      if (rescaledir > 0) {
        for (rescaledir = 1; rescaledir <= count; rescaledir++) {
          r->re *= 7.4428285367870146E+137;
          r->im *= 7.4428285367870146E+137;
        }
      } else {
        if (rescaledir < 0) {
          for (rescaledir = 1; rescaledir <= count; rescaledir++) {
            r->re *= 1.3435752215134178E-138;
            r->im *= 1.3435752215134178E-138;
          }
        }
      }
    }
  }
}

/*
 * ZP2SS  Zero-pole to state-space conversion. Codegen support
 *
 *  This function is based on 'zp2ss' by The MathWorks Inc.
 * Arguments    : emxArray_creal_T *a
 *                emxArray_real_T *b
 *                emxArray_creal_T *c
 *                double *d
 * Return Type  : void
 */
static void zp2ss_cg(emxArray_creal_T *a, emxArray_real_T *b, emxArray_creal_T
                     *c, double *d)
{
  int i5;
  creal_T b_c[3];
  int j;
  double den[3];
  creal_T x[2];
  static const creal_T pF[3] = { { -0.49999999999999978,/* re */
      0.86602540378443871              /* im */
    }, { -0.49999999999999978,         /* re */
      -0.86602540378443871             /* im */
    }, { -1.0,                         /* re */
      0.0                              /* im */
    } };

  double wn;
  double b1[2];
  int k;
  double v[2];
  double t[4];
  double b_den[4];
  double A[4];
  emxArray_int8_T *reshapes_f2;
  emxArray_cint8_T *b_a;
  emxArray_cint8_T *result;
  int i6;
  signed char a_re;
  signed char a_im;
  emxArray_creal_T *b_b1;
  creal_T b_A[4];
  emxArray_creal_T *varargin_2;
  double c_im;
  int result_re;
  int result_im;
  emxArray_creal_T *r5;

  /*  Strip infinities and throw away. */
  /*  Group into complex pairs */
  /*  try */
  /*      % z and p should have real elements and exact complex conjugate pair. */
  /*      z = cplxpair(zF,0); */
  /*      p = cplxpair(pF,0); */
  /*  catch */
  /*      % If fail, revert to use the old default tolerance. */
  /*      % The use of tolerance in checking for real entries and conjugate pairs */
  /*      % may result in misinterpretation for edge cases. Please review the */
  /*      % process of how z and p are generated. */
  /*      z = cplxpair(zF,1e6*nz*norm(zF)*eps + eps); */
  /*      p = cplxpair(pF,1e6*np*norm(pF)*eps + eps); */
  /*  end */
  /*  Initialize state-space matrices for running series */
  /*  If odd number of poles AND zeros, convert the pole and zero */
  /*  at the end into state-space. */
  /*    H(s) = (s-z1)/(s-p1) = (s + num(2)) / (s + den(2)) */
  /*  If odd number of poles only, convert the pole at the */
  /*  end into state-space. */
  /*   H(s) = 1/(s-p1) = 1/(s + den(2)) */
  i5 = a->size[0] * a->size[1];
  a->size[0] = 1;
  a->size[1] = 1;
  emxEnsureCapacity((emxArray__common *)a, i5, (int)sizeof(creal_T));
  a->data[0].re = -1.0;
  a->data[0].im = 0.0;
  i5 = b->size[0];
  b->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)b, i5, (int)sizeof(double));
  b->data[0] = 1.0;
  i5 = c->size[0] * c->size[1];
  c->size[0] = 1;
  c->size[1] = 1;
  emxEnsureCapacity((emxArray__common *)c, i5, (int)sizeof(creal_T));
  c->data[0].re = 1.0;
  c->data[0].im = 0.0;

  /*  If odd number of zeros only, convert the zero at the */
  /*  end, along with a pole-pair into state-space. */
  /*    H(s) = (s+num(2))/(s^2+den(2)s+den(3)) */
  /*  Now we have an even number of poles and zeros, although not */
  /*  necessarily the same number - there may be more poles. */
  /*    H(s) = (s^2+num(2)s+num(3))/(s^2+den(2)s+den(3)) */
  /*  Loop through rest of pairs, connecting in series to build the model. */
  /*  Take care of any left over unmatched pole pairs. */
  /*    H(s) = 1/(s^2+den(2)s+den(3)) */
  b_c[0].re = 1.0;
  b_c[0].im = 0.0;
  for (j = 0; j < 2; j++) {
    x[j] = pF[j];
    wn = b_c[j].re;
    b_c[j + 1].re = -x[j].re * b_c[j].re - -x[j].im * b_c[j].im;
    b_c[j + 1].im = -x[j].re * b_c[j].im + -x[j].im * wn;
    k = j;
    while (k + 1 > 1) {
      b_c[1].re -= x[j].re * b_c[0].re - x[j].im * b_c[0].im;
      b_c[1].im -= x[j].re * b_c[0].im + x[j].im * b_c[0].re;
      k = 0;
    }
  }

  for (i5 = 0; i5 < 3; i5++) {
    den[i5] = b_c[i5].re;
  }

  for (k = 0; k < 2; k++) {
    x[k] = pF[k];
    b1[k] = rt_hypotd_snf(x[k].re, x[k].im);
  }

  wn = sqrt(b1[0] * b1[1]);
  v[0] = 1.0;
  v[1] = 1.0 / wn;
  for (i5 = 0; i5 < 4; i5++) {
    t[i5] = 0.0;
  }

  /*  Balancing transformation */
  b_den[0] = -den[1];
  b_den[2] = -den[2];
  for (j = 0; j < 2; j++) {
    t[j + (j << 1)] = v[j];
    b_den[1 + (j << 1)] = 1.0 - (double)j;
  }

  mldivide(t, b_den, A);
  for (i5 = 0; i5 < 2; i5++) {
    v[i5] = 1.0 - (double)i5;
  }

  emxInit_int8_T(&reshapes_f2, 2);
  b_mldivide(t, v, b1);

  /*  [a,b,c,d] = series(a,b,c,d,a1,b1,c1,d1); */
  /*  Next lines perform series connection */
  i5 = reshapes_f2->size[0] * reshapes_f2->size[1];
  reshapes_f2->size[0] = 1;
  reshapes_f2->size[1] = 2;
  emxEnsureCapacity((emxArray__common *)reshapes_f2, i5, (int)sizeof(signed char));
  for (i5 = 0; i5 < 2; i5++) {
    reshapes_f2->data[i5] = 0;
  }

  emxInit_cint8_T(&b_a, 2);
  i5 = b_a->size[0] * b_a->size[1];
  b_a->size[0] = a->size[0];
  b_a->size[1] = a->size[1];
  emxEnsureCapacity((emxArray__common *)b_a, i5, (int)sizeof(cint8_T));
  j = a->size[1];
  for (i5 = 0; i5 < j; i5++) {
    k = a->size[0];
    for (i6 = 0; i6 < k; i6++) {
      a_re = (signed char)a->data[i6 + a->size[0] * i5].re;
      a_im = (signed char)a->data[i6 + a->size[0] * i5].im;
      b_a->data[i6 + b_a->size[0] * i5].re = a_re;
      b_a->data[i6 + b_a->size[0] * i5].im = a_im;
    }
  }

  emxInit_cint8_T(&result, 2);
  i5 = result->size[0] * result->size[1];
  result->size[0] = 1;
  result->size[1] = 1 + reshapes_f2->size[1];
  emxEnsureCapacity((emxArray__common *)result, i5, (int)sizeof(cint8_T));
  for (i5 = 0; i5 < 1; i5++) {
    for (i6 = 0; i6 < 1; i6++) {
      result->data[0] = b_a->data[0];
    }
  }

  emxFree_cint8_T(&b_a);
  j = reshapes_f2->size[1];
  for (i5 = 0; i5 < j; i5++) {
    k = reshapes_f2->size[0];
    for (i6 = 0; i6 < k; i6++) {
      result->data[i6 + result->size[0] * (i5 + 1)].re = reshapes_f2->data[i6 +
        reshapes_f2->size[0] * i5];
      result->data[i6 + result->size[0] * (i5 + 1)].im = 0;
    }
  }

  emxFree_int8_T(&reshapes_f2);
  for (i5 = 0; i5 < 2; i5++) {
    x[i5].re = b1[i5];
    x[i5].im = 0.0;
  }

  emxInit_creal_T(&b_b1, 2);
  i5 = b_b1->size[0] * b_b1->size[1];
  b_b1->size[0] = 2;
  b_b1->size[1] = c->size[1];
  emxEnsureCapacity((emxArray__common *)b_b1, i5, (int)sizeof(creal_T));
  for (i5 = 0; i5 < 2; i5++) {
    j = c->size[1];
    for (i6 = 0; i6 < j; i6++) {
      wn = c->data[c->size[0] * i6].re;
      c_im = c->data[c->size[0] * i6].im;
      b_b1->data[i5 + b_b1->size[0] * i6].re = x[i5].re * wn - x[i5].im * c_im;
      b_b1->data[i5 + b_b1->size[0] * i6].im = x[i5].re * c_im + x[i5].im * wn;
    }
  }

  for (i5 = 0; i5 < 2; i5++) {
    for (i6 = 0; i6 < 2; i6++) {
      wn = 0.0;
      for (j = 0; j < 2; j++) {
        wn += A[i5 + (j << 1)] * t[j + (i6 << 1)];
      }

      b_A[i5 + (i6 << 1)].re = wn;
      b_A[i5 + (i6 << 1)].im = 0.0;
    }
  }

  emxInit_creal_T(&varargin_2, 2);
  i5 = varargin_2->size[0] * varargin_2->size[1];
  varargin_2->size[0] = 2;
  varargin_2->size[1] = b_b1->size[1] + 2;
  emxEnsureCapacity((emxArray__common *)varargin_2, i5, (int)sizeof(creal_T));
  j = b_b1->size[1];
  for (i5 = 0; i5 < j; i5++) {
    for (i6 = 0; i6 < 2; i6++) {
      varargin_2->data[i6 + varargin_2->size[0] * i5] = b_b1->data[i6 +
        b_b1->size[0] * i5];
    }
  }

  for (i5 = 0; i5 < 2; i5++) {
    for (i6 = 0; i6 < 2; i6++) {
      varargin_2->data[i6 + varargin_2->size[0] * (i5 + b_b1->size[1])] = b_A[i6
        + (i5 << 1)];
    }
  }

  emxFree_creal_T(&b_b1);
  if (!(result->size[1] == 0)) {
    j = result->size[1];
  } else {
    j = 3;
  }

  k = !(result->size[1] == 0);
  i5 = a->size[0] * a->size[1];
  a->size[0] = k + 2;
  a->size[1] = j;
  emxEnsureCapacity((emxArray__common *)a, i5, (int)sizeof(creal_T));
  for (i5 = 0; i5 < j; i5++) {
    for (i6 = 0; i6 < k; i6++) {
      result_re = result->data[i6 + k * i5].re;
      result_im = result->data[i6 + k * i5].im;
      a->data[i6 + a->size[0] * i5].re = result_re;
      a->data[i6 + a->size[0] * i5].im = result_im;
    }
  }

  emxFree_cint8_T(&result);
  for (i5 = 0; i5 < j; i5++) {
    for (i6 = 0; i6 < 2; i6++) {
      a->data[(i6 + k) + a->size[0] * i5] = varargin_2->data[i6 + (i5 << 1)];
    }
  }

  emxFree_creal_T(&varargin_2);
  j = b->size[0];
  i5 = b->size[0];
  b->size[0] = j + 2;
  emxEnsureCapacity((emxArray__common *)b, i5, (int)sizeof(double));
  for (i5 = 0; i5 < 2; i5++) {
    b->data[j + i5] = b1[i5] * 0.0;
  }

  for (i5 = 0; i5 < 2; i5++) {
    wn = 0.0;
    for (i6 = 0; i6 < 2; i6++) {
      wn += (double)i6 * t[i6 + (i5 << 1)];
    }

    x[i5].re = wn;
    x[i5].im = 0.0;
  }

  emxInit_creal_T(&r5, 2);
  i5 = r5->size[0] * r5->size[1];
  r5->size[0] = 1;
  r5->size[1] = c->size[1] + 2;
  emxEnsureCapacity((emxArray__common *)r5, i5, (int)sizeof(creal_T));
  j = c->size[1];
  for (i5 = 0; i5 < j; i5++) {
    r5->data[r5->size[0] * i5].re = 0.0 * c->data[c->size[0] * i5].re;
    r5->data[r5->size[0] * i5].im = 0.0 * c->data[c->size[0] * i5].im;
  }

  for (i5 = 0; i5 < 2; i5++) {
    r5->data[r5->size[0] * (i5 + c->size[1])] = x[i5];
  }

  i5 = c->size[0] * c->size[1];
  c->size[0] = 1;
  c->size[1] = r5->size[1];
  emxEnsureCapacity((emxArray__common *)c, i5, (int)sizeof(creal_T));
  j = r5->size[1];
  for (i5 = 0; i5 < j; i5++) {
    c->data[c->size[0] * i5] = r5->data[r5->size[0] * i5];
  }

  emxFree_creal_T(&r5);

  /*  Apply gain k: */
  i5 = c->size[0] * c->size[1];
  c->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)c, i5, (int)sizeof(creal_T));
  j = c->size[0];
  k = c->size[1];
  j *= k;
  for (i5 = 0; i5 < j; i5++) {
    c->data[i5].re *= 0.99999999999999989;
    c->data[i5].im *= 0.99999999999999989;
  }

  *d = 0.0;
}

/*
 * %%%%%%%%%%%%%%%%%%% Perform Input Checks
 * Arguments    : double Rdata
 *                double Fpass
 *                double Fstop
 *                double caldiv
 *                double FIR
 *                double HB1
 *                double PLL_mult
 *                double Apass
 *                double Astop
 *                double phEQ
 *                double HB2
 *                double HB3
 *                const char Type[7]
 *                const char RxTx[2]
 *                double RFbw
 *                double DAC_div
 *                double converter_rate
 *                double PLL_rate
 *                double Fcenter
 *                double wnom
 *                double FIRdBmin
 *                double int_FIR
 *                short outputTaps[128]
 * Return Type  : void
 */
void internal_design_filter_cg(double Rdata, double Fpass, double Fstop, double
  caldiv, double FIR, double HB1, double PLL_mult, double Apass, double Astop,
  double phEQ, double HB2, double HB3, const char Type[7], const char RxTx[2],
  double RFbw, double DAC_div, double converter_rate, double PLL_rate, double
  Fcenter, double wnom, double FIRdBmin, double int_FIR, short outputTaps[128])
{
  emxArray_creal_T *a1;
  emxArray_creal_T *a2;
  double b1[4];
  double b2[2];
  creal_T a2_data[3];
  int a2_size[2];
  int b1_size[2];
  double b1_data[4];
  int i0;
  int b2_size[2];
  int loop_ub;
  double b2_data[4];
  double hb1_coeff[15];
  static const double dv0[15] = { -0.00323486328125, 0.0, 0.01910400390625, 0.0,
    -0.07049560546875, 0.0, 0.30450439453125, 0.5, 0.30450439453125, 0.0,
    -0.07049560546875, 0.0, 0.01910400390625, 0.0, -0.00323486328125 };

  static const double dv1[15] = { -0.00390625, 0.0, 0.0205078125, 0.0,
    -0.07177734375, 0.0, 0.30224609375, 0.49462890625, 0.30224609375, 0.0,
    -0.07177734375, 0.0, 0.0205078125, 0.0, -0.00390625 };

  int hb3_coeff_size[2];
  double hb3_coeff_data[5];
  int dec_int3_coeff_size[2];
  static const double y[3] = { 0.25, 0.5, 0.25 };

  static const double b_y[5] = { 0.0625, 0.25, 0.375, 0.25, 0.0625 };

  double dec_int3_coeff_data[29];
  static const double c_y[29] = { 0.00146484375, -0.00077311197916666663, 0.0,
    -0.00634765625, -0.00048828125, 0.0, 0.019490559895833332,
    0.0090738932291666661, 0.0, -0.0494384765625, -0.0404052734375, 0.0,
    0.14522298177083331, 0.25541178385416663, 0.33333333333333331,
    0.25541178385416663, 0.14522298177083331, 0.0, -0.0404052734375,
    -0.0494384765625, 0.0, 0.0090738932291666661, 0.019490559895833332, 0.0,
    -0.00048828125, -0.00634765625, 0.0, -0.00077311197916666663, 0.00146484375
  };

  static const double d_y[17] = { 0.00335693359375, 0.00506591796875, 0.0,
    -0.02398681640625, -0.035400390625, 0.0, 0.1168212890625, 0.24664306640625,
    0.3125, 0.24664306640625, 0.1168212890625, 0.0, -0.035400390625,
    -0.02398681640625, 0.0, 0.00506591796875, 0.00335693359375 };

  int hb3;
  int dec_int3;
  double sigma;
  signed char i1;
  char enables[4];
  double w[2048];
  double phi[2048];
  int i;
  static const double hb2_coeff[7] = { -0.03515625, 0.0, 0.28515625, 0.5,
    0.28515625, 0.0, -0.03515625 };

  static creal_T combinedResponse[2048];
  static creal_T dcv0[2048];
  double b_combinedResponse[2048];
  double invariance[2048];
  double sigmax;
  double b_phi[2048];
  double b_w[2048];
  double clkFIR;
  double Gpass;
  double Gstop;
  emxArray_real_T *fg;
  emxArray_real_T *omega;
  emxArray_creal_T *c_combinedResponse;
  emxArray_creal_T *r0;
  int d_combinedResponse;
  emxArray_creal_T *rg2;
  double apnd;
  double absa;
  double im;
  emxArray_creal_T *rg;
  emxArray_real_T *F3;
  emxArray_real_T *sw;
  emxArray_real_T *fg2;
  int n;
  emxArray_real_T *omega2;
  emxArray_creal_T *rgN;
  int i2;
  emxArray_real_T *b_omega2;
  int i3;
  emxArray_real_T *c_omega2;
  emxArray_creal_T *e_combinedResponse;
  emxArray_real_T *x;
  emxArray_real_T *weight;
  emxArray_real_T *b_weight;
  boolean_T exitg7;
  emxArray_real_T *F1;
  emxArray_real_T *F2;
  emxArray_real_T *A1;
  emxArray_real_T *A2;
  emxArray_real_T *W1;
  emxArray_real_T *W2;
  int Nmax;
  double tap_store_data[1024];
  double Apass_actual_vector_data[8];
  double Astop_actual_vector_data[8];
  emxArray_real_T *ccoef;
  emxArray_real_T *F4;
  emxArray_real_T *h;
  emxArray_creal_T *b_rg;
  emxArray_real_T *d_omega2;
  emxArray_creal_T *b_rg2;
  emxArray_real_T *e_omega2;
  emxArray_real_T *f_omega2;
  emxArray_real_T *g_omega2;
  emxArray_creal_T *r1;
  emxArray_real_T *h_omega2;
  emxArray_creal_T *r2;
  emxArray_real_T *i_omega2;
  emxArray_real_T *j_omega2;
  emxArray_real_T *k_omega2;
  emxArray_real_T *b_omega;
  emxArray_real_T *b_F3;
  emxArray_real_T *b_sw;
  emxArray_real_T *b_A1;
  emxArray_real_T *b_F1;
  emxArray_real_T *b_W1;
  emxArray_real_T *c_A1;
  emxArray_real_T *c_F1;
  emxArray_real_T *c_W1;
  int exitg2;
  double d_F1[4];
  boolean_T exitg6;
  double e_F1[4];
  double b_tap_store_data[128];
  short i4;
  boolean_T valid;
  int tap_store_size[2];
  int b_tap_store_size[2];
  int c_tap_store_size[2];
  int d_tap_store_size[2];
  double c_F3[4];
  boolean_T exitg5;
  boolean_T exitg4;
  boolean_T exitg3;
  emxArray_real_T *r3;
  emxArray_real_T *r4;
  double firTapsPreScale[128];
  boolean_T exitg1;
  (void)caldiv;
  (void)PLL_mult;
  (void)Type;
  (void)RFbw;
  (void)PLL_rate;
  (void)Fcenter;

  /*  internal_design_filter_cg */
  /*  */
  /*  This implementation of the ADI ad936x-filter-wizard supports */
  /*  code generation from MATLAB Coder */
  /*  */
  /*  Edits by: Travis F. Collins <travisfcollins@gmail.com> */
  /*  */
  /*  When calling this function utilize the function "process_input" to make */
  /*  sure all the necessary fields exist */
  /*  */
  /*  Todo: */
  /*  - Set fractionalLength based on FSR of inputs (ML has no doc on this) */
  /*  */
  /*  Inputs */
  /*  ============================================ */
  /*  Rdata      = Input/output sample data rate (in Hz) */
  /*  Fpass      = Passband frequency (in Hz) */
  /*  Fstop      = Stopband frequency (in Hz) */
  /*  caldiv     = The actual discrete register value that describes the */
  /*               rolloff for the analog filters */
  /*  FIR        = FIR interpolation/decimation factor */
  /*  HB1        = HB1 interpolation/decimation rates */
  /*  PLL_mult   = PLL multiplication */
  /*  Apass      = Max ripple allowed in passband (in dB) */
  /*  Astop      = Min attenuation in stopband (in dB) */
  /*  phEQ       = Phase equalization on (not -1)/off (-1) */
  /*  HB2        = HB2 interpolation/decimation rates */
  /*  HB3        = HB3 interpolation/decimation rates */
  /*  Type       = The type of filter required. one of: */
  /*               'Lowpass' (default,and only support) */
  /*               'Bandpass' (Not implemented) */
  /*               'Equalize' (Not implemented) */
  /*               'Root Raised Cosine' (Not implemented) */
  /*  RxTx       = Is this 'Rx' or 'Tx'?. */
  /*  RFbw */
  /*  DAC_div    = The ADC/DAC ratio, for Rx channels, this is */
  /*               always '1', for Tx, it is either '1' or '2' */
  /*  converter_rate = Rate of converter */
  /*  PLL_rate   = the PLL rate in Hz */
  /*  Fcenter    = Center Frequency in Hz (only used for Bandpass), */
  /*               otherwise 0 */
  /*  wnom       = analog cutoff frequency (in Hz) */
  /*  FIRdBmin   = min rejection that FIR is required to have (in dB) */
  /*  int_FIR    = use AD9361 FIR on (1)/off (0) */
  /*  */
  /*  Outputs */
  /*  =============================================== */
  /*  firtaps          = fixed point FIR coefficients */
  /*  Scalar checks */
  /*  String checks */
  /*  This is unused */
  /* %%%%%%%%%%%%%%%%%%% Build processing struct */
  /*  Not used */
  /*  Design analog filters */
  emxInit_creal_T(&a1, 2);
  emxInit_creal_T(&a2, 2);
  if (b_strcmp(RxTx)) {
    /*  Define the analog filters (for design purpose) */
    butter_cg(6.2831853071795862 * (wnom * 1.7857142857142858), b2, a2_data,
              a2_size);
    b1_size[0] = 1;
    b1_size[1] = 2;
    for (i0 = 0; i0 < 2; i0++) {
      b1_data[i0] = b2[i0];
    }

    i0 = a1->size[0] * a1->size[1];
    a1->size[0] = 1;
    a1->size[1] = a2_size[1];
    emxEnsureCapacity((emxArray__common *)a1, i0, (int)sizeof(creal_T));
    loop_ub = a2_size[0] * a2_size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      a1->data[i0] = a2_data[i0];
    }

    /*  1st order */
    b_butter_cg(6.2831853071795862 * wnom, b1, a2);
    b2_size[0] = 1;
    b2_size[1] = 4;
    for (i0 = 0; i0 < 4; i0++) {
      b2_data[i0] = b1[i0];
    }

    /*  3rd order */
    /*  Define the digital filters with fixed coefficients */
    memcpy(&hb1_coeff[0], &dv1[0], 15U * sizeof(double));
    hb3_coeff_size[0] = 1;
    hb3_coeff_size[1] = 5;
    for (i0 = 0; i0 < 5; i0++) {
      hb3_coeff_data[i0] = b_y[i0];
    }

    dec_int3_coeff_size[0] = 1;
    dec_int3_coeff_size[1] = 17;
    memcpy(&dec_int3_coeff_data[0], &d_y[0], 17U * sizeof(double));
  } else {
    /*  Define the analog filters (for design purpose) */
    b_butter_cg(6.2831853071795862 * wnom, b1, a1);
    b1_size[0] = 1;
    b1_size[1] = 4;
    for (i0 = 0; i0 < 4; i0++) {
      b1_data[i0] = b1[i0];
    }

    /*  3rd order */
    butter_cg(6.2831853071795862 * (wnom * 3.125), b2, a2_data, a2_size);
    b2_size[0] = 1;
    b2_size[1] = 2;
    for (i0 = 0; i0 < 2; i0++) {
      b2_data[i0] = b2[i0];
    }

    i0 = a2->size[0] * a2->size[1];
    a2->size[0] = 1;
    a2->size[1] = a2_size[1];
    emxEnsureCapacity((emxArray__common *)a2, i0, (int)sizeof(creal_T));
    loop_ub = a2_size[0] * a2_size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      a2->data[i0] = a2_data[i0];
    }

    /*  1st order */
    /*  Define the digital filters with fixed coefficients */
    memcpy(&hb1_coeff[0], &dv0[0], 15U * sizeof(double));
    hb3_coeff_size[0] = 1;
    hb3_coeff_size[1] = 3;
    for (i0 = 0; i0 < 3; i0++) {
      hb3_coeff_data[i0] = y[i0];
    }

    dec_int3_coeff_size[0] = 1;
    dec_int3_coeff_size[1] = 29;
    memcpy(&dec_int3_coeff_data[0], &c_y[0], 29U * sizeof(double));
  }

  /*  Configure staging of filters */
  if (HB3 == 2.0) {
    hb3 = 2;
    dec_int3 = 49;
  } else if (HB3 == 3.0) {
    hb3 = 1;
    dec_int3 = 51;
  } else {
    hb3 = 1;
    dec_int3 = 49;
  }

  /*  convert the enables into a string */
  sigma = rt_roundd_snf(HB1);
  if (sigma < 128.0) {
    if (sigma >= -128.0) {
      i1 = (signed char)sigma;
    } else {
      i1 = MIN_int8_T;
    }
  } else if (sigma >= 128.0) {
    i1 = MAX_int8_T;
  } else {
    i1 = 0;
  }

  i0 = 48 + i1;
  if (i0 > 127) {
    i0 = 127;
  }

  enables[0] = (signed char)i0;
  sigma = rt_roundd_snf(HB2);
  if (sigma < 128.0) {
    if (sigma >= -128.0) {
      i1 = (signed char)sigma;
    } else {
      i1 = MIN_int8_T;
    }
  } else if (sigma >= 128.0) {
    i1 = MAX_int8_T;
  } else {
    i1 = 0;
  }

  i0 = 48 + i1;
  if (i0 > 127) {
    i0 = 127;
  }

  enables[1] = (signed char)i0;
  enables[2] = (signed char)(hb3 + 48);
  enables[3] = (signed char)dec_int3;

  /*  Find out the best fit delay on passband */
  memset(&w[0], 0, sizeof(double) << 11);
  memset(&phi[0], 0, sizeof(double) << 11);
  w[0] = -Fpass;
  for (i = 0; i < 2047; i++) {
    w[i + 1] = w[0] - 2.0 * w[0] * (2.0 + (double)i) / 2048.0;
  }

  /*  Generate target responses used in filter design phase */
  /*  Generate responses then convolve */
  generateCascadedResponseRx(enables, w, converter_rate, hb1_coeff, hb2_coeff,
    hb3_coeff_data, hb3_coeff_size, dec_int3_coeff_data, dec_int3_coeff_size,
    combinedResponse);

  /*  Determine overall response with analog filters inline */
  analogresp(RxTx, w, converter_rate, b1_data, b1_size, a1, b2_data, b2_size, a2,
             dcv0);
  for (i0 = 0; i0 < 2048; i0++) {
    sigmax = combinedResponse[i0].re;
    combinedResponse[i0].re = combinedResponse[i0].re * dcv0[i0].re -
      combinedResponse[i0].im * dcv0[i0].im;
    combinedResponse[i0].im = sigmax * dcv0[i0].im + combinedResponse[i0].im *
      dcv0[i0].re;
    b_combinedResponse[i0] = combinedResponse[i0].re;
  }

  power(b_combinedResponse, invariance);
  for (i0 = 0; i0 < 2048; i0++) {
    b_combinedResponse[i0] = combinedResponse[i0].im;
  }

  power(b_combinedResponse, b_phi);
  for (i0 = 0; i0 < 2048; i0++) {
    invariance[i0] += b_phi[i0];
  }

  phi[0] = rt_atan2d_snf(combinedResponse[0].im, combinedResponse[0].re);
  for (i = 0; i < 2047; i++) {
    sigma = rt_atan2d_snf(combinedResponse[i + 1].im, combinedResponse[i + 1].re)
      - phi[i];
    phi[i + 1] = phi[i] + (sigma - 6.2831853071795862 * floor(sigma /
      6.2831853071795862 + 0.5));
  }

  sigma = sum(invariance);
  for (i0 = 0; i0 < 2048; i0++) {
    b_combinedResponse[i0] = w[i0] * invariance[i0];
  }

  sigmax = sum(b_combinedResponse);
  if ((phEQ == 0.0) || (phEQ == -1.0)) {
    for (i0 = 0; i0 < 2048; i0++) {
      b_combinedResponse[i0] = w[i0] * phi[i0] * invariance[i0];
      b_phi[i0] = phi[i0] * invariance[i0];
      b_w[i0] = w[i0] * w[i0] * invariance[i0];
    }

    sigma = -((sigma * sum(b_combinedResponse) - sigmax * sum(b_phi)) / (sigma *
               sum(b_w) - sigmax * sigmax)) / 6.2831853071795862;
  } else {
    sigma = phEQ * 1.0E-9;
  }

  /*  Design the FIR */
  clkFIR = Rdata * FIR;
  Gpass = floor(16384.0 * Fpass / clkFIR);
  Gstop = ceil(16384.0 * Fstop / clkFIR);
  if (!((Gpass <= Gstop - 1.0) || rtIsNaN(Gstop - 1.0))) {
    Gpass = Gstop - 1.0;
  }

  emxInit_real_T(&fg, 2);
  i0 = fg->size[0] * fg->size[1];
  fg->size[0] = 1;
  fg->size[1] = (int)(Gpass + 1.0);
  emxEnsureCapacity((emxArray__common *)fg, i0, (int)sizeof(double));
  loop_ub = (int)(Gpass + 1.0);
  for (i0 = 0; i0 < loop_ub; i0++) {
    fg->data[i0] = 0.0;
  }

  emxInit_real_T(&omega, 2);
  i0 = omega->size[0] * omega->size[1];
  omega->size[0] = 1;
  omega->size[1] = (int)(Gpass + 1.0);
  emxEnsureCapacity((emxArray__common *)omega, i0, (int)sizeof(double));
  loop_ub = (int)(Gpass + 1.0);
  for (i0 = 0; i0 < loop_ub; i0++) {
    omega->data[i0] = 0.0;
  }

  /*  passband */
  for (i = 0; i < (int)(Gpass + 1.0); i++) {
    fg->data[i] = ((1.0 + (double)i) - 1.0) / 16384.0;
    omega->data[i] = fg->data[i] * clkFIR;
  }

  emxInit_creal_T(&c_combinedResponse, 2);
  emxInit_creal_T(&r0, 2);

  /*  Generate responses then convolve */
  b_generateCascadedResponseRx(enables, omega, converter_rate, hb1_coeff,
    hb2_coeff, hb3_coeff_data, hb3_coeff_size, dec_int3_coeff_data,
    dec_int3_coeff_size, c_combinedResponse);

  /*  Determine overall response with analog filters inline */
  b_analogresp(RxTx, omega, converter_rate, b1_data, b1_size, a1, b2_data,
               b2_size, a2, r0);
  i0 = c_combinedResponse->size[0] * c_combinedResponse->size[1];
  c_combinedResponse->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)c_combinedResponse, i0, (int)sizeof
                    (creal_T));
  dec_int3 = c_combinedResponse->size[0];
  d_combinedResponse = c_combinedResponse->size[1];
  loop_ub = dec_int3 * d_combinedResponse;
  for (i0 = 0; i0 < loop_ub; i0++) {
    sigmax = c_combinedResponse->data[i0].re;
    apnd = c_combinedResponse->data[i0].im;
    absa = r0->data[i0].re;
    im = r0->data[i0].im;
    c_combinedResponse->data[i0].re = sigmax * absa - apnd * im;
    c_combinedResponse->data[i0].im = sigmax * im + apnd * absa;
  }

  emxInit_creal_T(&rg2, 2);
  i0 = rg2->size[0] * rg2->size[1];
  rg2->size[0] = 1;
  rg2->size[1] = omega->size[1];
  emxEnsureCapacity((emxArray__common *)rg2, i0, (int)sizeof(creal_T));
  loop_ub = omega->size[0] * omega->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    absa = omega->data[i0] * -0.0;
    im = omega->data[i0] * -6.2831853071795862;
    rg2->data[i0].re = sigma * absa;
    rg2->data[i0].im = sigma * im;
  }

  emxInit_creal_T(&rg, 2);
  emxInit_real_T(&F3, 2);
  c_exp(rg2);
  b_rdivide(rg2, c_combinedResponse, rg);
  b_abs(c_combinedResponse, F3);
  sigma = Gpass + 1.0;

  /*  Expand memory correctly */
  emxInit_real_T(&sw, 2);
  if (rtIsNaN(Gstop)) {
    i0 = sw->size[0] * sw->size[1];
    sw->size[0] = 1;
    sw->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)sw, i0, (int)sizeof(double));
    sw->data[0] = rtNaN;
  } else if (8192.0 < Gstop) {
    i0 = sw->size[0] * sw->size[1];
    sw->size[0] = 1;
    sw->size[1] = 0;
    emxEnsureCapacity((emxArray__common *)sw, i0, (int)sizeof(double));
  } else if (rtIsInf(Gstop) && (Gstop == 8192.0)) {
    i0 = sw->size[0] * sw->size[1];
    sw->size[0] = 1;
    sw->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)sw, i0, (int)sizeof(double));
    sw->data[0] = rtNaN;
  } else if (Gstop == Gstop) {
    i0 = sw->size[0] * sw->size[1];
    sw->size[0] = 1;
    sw->size[1] = (int)(8192.0 - Gstop) + 1;
    emxEnsureCapacity((emxArray__common *)sw, i0, (int)sizeof(double));
    loop_ub = (int)(8192.0 - Gstop);
    for (i0 = 0; i0 <= loop_ub; i0++) {
      sw->data[sw->size[0] * i0] = Gstop + (double)i0;
    }
  } else {
    sigmax = floor((8192.0 - Gstop) + 0.5);
    apnd = Gstop + sigmax;
    absa = fabs(Gstop);
    if (!(absa >= 8192.0)) {
      absa = 8192.0;
    }

    if (fabs(apnd - 8192.0) < 4.4408920985006262E-16 * absa) {
      sigmax++;
      apnd = 8192.0;
    } else if (apnd - 8192.0 > 0.0) {
      apnd = Gstop + (sigmax - 1.0);
    } else {
      sigmax++;
    }

    if (sigmax >= 0.0) {
      n = (int)sigmax;
    } else {
      n = 0;
    }

    i0 = sw->size[0] * sw->size[1];
    sw->size[0] = 1;
    sw->size[1] = n;
    emxEnsureCapacity((emxArray__common *)sw, i0, (int)sizeof(double));
    if (n > 0) {
      sw->data[0] = Gstop;
      if (n > 1) {
        sw->data[n - 1] = apnd;
        dec_int3 = (n - 1) / 2;
        for (d_combinedResponse = 1; d_combinedResponse < dec_int3;
             d_combinedResponse++) {
          sw->data[d_combinedResponse] = Gstop + (double)d_combinedResponse;
          sw->data[(n - d_combinedResponse) - 1] = apnd - (double)
            d_combinedResponse;
        }

        if (dec_int3 << 1 == n - 1) {
          sw->data[dec_int3] = (Gstop + apnd) / 2.0;
        } else {
          sw->data[dec_int3] = Gstop + (double)dec_int3;
          sw->data[dec_int3 + 1] = apnd - (double)dec_int3;
        }
      }
    }
  }

  emxInit_real_T(&fg2, 2);
  i0 = fg2->size[0] * fg2->size[1];
  fg2->size[0] = 1;
  fg2->size[1] = (int)((unsigned int)sw->size[1] + fg->size[1]);
  emxEnsureCapacity((emxArray__common *)fg2, i0, (int)sizeof(double));
  loop_ub = (int)((unsigned int)sw->size[1] + fg->size[1]);
  for (i0 = 0; i0 < loop_ub; i0++) {
    fg2->data[i0] = 0.0;
  }

  loop_ub = fg->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    fg2->data[i0] = fg->data[fg->size[0] * i0];
  }

  if (rtIsNaN(Gstop)) {
    i0 = sw->size[0] * sw->size[1];
    sw->size[0] = 1;
    sw->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)sw, i0, (int)sizeof(double));
    sw->data[0] = rtNaN;
  } else if (8192.0 < Gstop) {
    i0 = sw->size[0] * sw->size[1];
    sw->size[0] = 1;
    sw->size[1] = 0;
    emxEnsureCapacity((emxArray__common *)sw, i0, (int)sizeof(double));
  } else if (rtIsInf(Gstop) && (Gstop == 8192.0)) {
    i0 = sw->size[0] * sw->size[1];
    sw->size[0] = 1;
    sw->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)sw, i0, (int)sizeof(double));
    sw->data[0] = rtNaN;
  } else if (Gstop == Gstop) {
    i0 = sw->size[0] * sw->size[1];
    sw->size[0] = 1;
    sw->size[1] = (int)(8192.0 - Gstop) + 1;
    emxEnsureCapacity((emxArray__common *)sw, i0, (int)sizeof(double));
    loop_ub = (int)(8192.0 - Gstop);
    for (i0 = 0; i0 <= loop_ub; i0++) {
      sw->data[sw->size[0] * i0] = Gstop + (double)i0;
    }
  } else {
    sigmax = floor((8192.0 - Gstop) + 0.5);
    apnd = Gstop + sigmax;
    absa = fabs(Gstop);
    if (!(absa >= 8192.0)) {
      absa = 8192.0;
    }

    if (fabs(apnd - 8192.0) < 4.4408920985006262E-16 * absa) {
      sigmax++;
      apnd = 8192.0;
    } else if (apnd - 8192.0 > 0.0) {
      apnd = Gstop + (sigmax - 1.0);
    } else {
      sigmax++;
    }

    if (sigmax >= 0.0) {
      n = (int)sigmax;
    } else {
      n = 0;
    }

    i0 = sw->size[0] * sw->size[1];
    sw->size[0] = 1;
    sw->size[1] = n;
    emxEnsureCapacity((emxArray__common *)sw, i0, (int)sizeof(double));
    if (n > 0) {
      sw->data[0] = Gstop;
      if (n > 1) {
        sw->data[n - 1] = apnd;
        dec_int3 = (n - 1) / 2;
        for (d_combinedResponse = 1; d_combinedResponse < dec_int3;
             d_combinedResponse++) {
          sw->data[d_combinedResponse] = Gstop + (double)d_combinedResponse;
          sw->data[(n - d_combinedResponse) - 1] = apnd - (double)
            d_combinedResponse;
        }

        if (dec_int3 << 1 == n - 1) {
          sw->data[dec_int3] = (Gstop + apnd) / 2.0;
        } else {
          sw->data[dec_int3] = Gstop + (double)dec_int3;
          sw->data[dec_int3 + 1] = apnd - (double)dec_int3;
        }
      }
    }
  }

  emxInit_real_T(&omega2, 2);
  i0 = omega2->size[0] * omega2->size[1];
  omega2->size[0] = 1;
  omega2->size[1] = (int)((unsigned int)sw->size[1] + omega->size[1]);
  emxEnsureCapacity((emxArray__common *)omega2, i0, (int)sizeof(double));
  loop_ub = (int)((unsigned int)sw->size[1] + omega->size[1]);
  for (i0 = 0; i0 < loop_ub; i0++) {
    omega2->data[i0] = 0.0;
  }

  loop_ub = omega->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    omega2->data[i0] = omega->data[omega->size[0] * i0];
  }

  if (rtIsNaN(Gstop)) {
    i0 = sw->size[0] * sw->size[1];
    sw->size[0] = 1;
    sw->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)sw, i0, (int)sizeof(double));
    sw->data[0] = rtNaN;
  } else if (8192.0 < Gstop) {
    i0 = sw->size[0] * sw->size[1];
    sw->size[0] = 1;
    sw->size[1] = 0;
    emxEnsureCapacity((emxArray__common *)sw, i0, (int)sizeof(double));
  } else if (rtIsInf(Gstop) && (Gstop == 8192.0)) {
    i0 = sw->size[0] * sw->size[1];
    sw->size[0] = 1;
    sw->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)sw, i0, (int)sizeof(double));
    sw->data[0] = rtNaN;
  } else if (Gstop == Gstop) {
    i0 = sw->size[0] * sw->size[1];
    sw->size[0] = 1;
    sw->size[1] = (int)(8192.0 - Gstop) + 1;
    emxEnsureCapacity((emxArray__common *)sw, i0, (int)sizeof(double));
    loop_ub = (int)(8192.0 - Gstop);
    for (i0 = 0; i0 <= loop_ub; i0++) {
      sw->data[sw->size[0] * i0] = Gstop + (double)i0;
    }
  } else {
    sigmax = floor((8192.0 - Gstop) + 0.5);
    apnd = Gstop + sigmax;
    absa = fabs(Gstop);
    if (!(absa >= 8192.0)) {
      absa = 8192.0;
    }

    if (fabs(apnd - 8192.0) < 4.4408920985006262E-16 * absa) {
      sigmax++;
      apnd = 8192.0;
    } else if (apnd - 8192.0 > 0.0) {
      apnd = Gstop + (sigmax - 1.0);
    } else {
      sigmax++;
    }

    if (sigmax >= 0.0) {
      n = (int)sigmax;
    } else {
      n = 0;
    }

    i0 = sw->size[0] * sw->size[1];
    sw->size[0] = 1;
    sw->size[1] = n;
    emxEnsureCapacity((emxArray__common *)sw, i0, (int)sizeof(double));
    if (n > 0) {
      sw->data[0] = Gstop;
      if (n > 1) {
        sw->data[n - 1] = apnd;
        dec_int3 = (n - 1) / 2;
        for (d_combinedResponse = 1; d_combinedResponse < dec_int3;
             d_combinedResponse++) {
          sw->data[d_combinedResponse] = Gstop + (double)d_combinedResponse;
          sw->data[(n - d_combinedResponse) - 1] = apnd - (double)
            d_combinedResponse;
        }

        if (dec_int3 << 1 == n - 1) {
          sw->data[dec_int3] = (Gstop + apnd) / 2.0;
        } else {
          sw->data[dec_int3] = Gstop + (double)dec_int3;
          sw->data[dec_int3 + 1] = apnd - (double)dec_int3;
        }
      }
    }
  }

  emxInit_creal_T(&rgN, 2);
  i0 = rgN->size[0] * rgN->size[1];
  rgN->size[0] = 1;
  rgN->size[1] = (int)((unsigned int)sw->size[1] + rg->size[1]);
  emxEnsureCapacity((emxArray__common *)rgN, i0, (int)sizeof(creal_T));
  loop_ub = (int)((unsigned int)sw->size[1] + rg->size[1]);
  for (i0 = 0; i0 < loop_ub; i0++) {
    rgN->data[i0].re = 0.0;
    rgN->data[i0].im = 0.0;
  }

  loop_ub = rg->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    rgN->data[i0] = rg->data[rg->size[0] * i0];
  }

  /*  stop band */
  for (dec_int3 = 0; dec_int3 < (int)(8192.0 + (1.0 - Gstop)); dec_int3++) {
    sigma++;
    fg2->data[(int)sigma - 1] = (Gstop + (double)dec_int3) / 16384.0;
    omega2->data[(int)sigma - 1] = fg2->data[(int)sigma - 1] * clkFIR;
    rgN->data[(int)sigma - 1].re = 0.0;
    rgN->data[(int)sigma - 1].im = 0.0;
  }

  /*  Generate responses then convolve */
  if (Gpass + 2.0 > omega2->size[1]) {
    i0 = 0;
    i2 = 0;
  } else {
    i0 = (int)(Gpass + 2.0) - 1;
    i2 = omega2->size[1];
  }

  emxInit_real_T(&b_omega2, 2);
  i3 = b_omega2->size[0] * b_omega2->size[1];
  b_omega2->size[0] = 1;
  b_omega2->size[1] = i2 - i0;
  emxEnsureCapacity((emxArray__common *)b_omega2, i3, (int)sizeof(double));
  loop_ub = i2 - i0;
  for (i2 = 0; i2 < loop_ub; i2++) {
    b_omega2->data[b_omega2->size[0] * i2] = omega2->data[i0 + i2];
  }

  b_generateCascadedResponseRx(enables, b_omega2, converter_rate, hb1_coeff,
    hb2_coeff, hb3_coeff_data, hb3_coeff_size, dec_int3_coeff_data,
    dec_int3_coeff_size, c_combinedResponse);
  emxFree_real_T(&b_omega2);
  if (Gpass + 2.0 > omega2->size[1]) {
    i0 = 0;
    i2 = 0;
  } else {
    i0 = (int)(Gpass + 2.0) - 1;
    i2 = omega2->size[1];
  }

  emxInit_real_T(&c_omega2, 2);
  i3 = c_omega2->size[0] * c_omega2->size[1];
  c_omega2->size[0] = 1;
  c_omega2->size[1] = i2 - i0;
  emxEnsureCapacity((emxArray__common *)c_omega2, i3, (int)sizeof(double));
  loop_ub = i2 - i0;
  for (i2 = 0; i2 < loop_ub; i2++) {
    c_omega2->data[c_omega2->size[0] * i2] = omega2->data[i0 + i2];
  }

  emxInit_creal_T(&e_combinedResponse, 2);
  b_analogresp(RxTx, c_omega2, converter_rate, b1_data, b1_size, a1, b2_data,
               b2_size, a2, r0);
  i0 = e_combinedResponse->size[0] * e_combinedResponse->size[1];
  e_combinedResponse->size[0] = 1;
  e_combinedResponse->size[1] = c_combinedResponse->size[1];
  emxEnsureCapacity((emxArray__common *)e_combinedResponse, i0, (int)sizeof
                    (creal_T));
  loop_ub = c_combinedResponse->size[0] * c_combinedResponse->size[1];
  emxFree_real_T(&c_omega2);
  for (i0 = 0; i0 < loop_ub; i0++) {
    sigmax = c_combinedResponse->data[i0].re;
    apnd = c_combinedResponse->data[i0].im;
    absa = r0->data[i0].re;
    im = r0->data[i0].im;
    e_combinedResponse->data[i0].re = sigmax * absa - apnd * im;
    e_combinedResponse->data[i0].im = sigmax * im + apnd * absa;
  }

  emxFree_creal_T(&c_combinedResponse);
  b_abs(e_combinedResponse, fg);
  emxFree_creal_T(&e_combinedResponse);
  if (b_strcmp(RxTx)) {
    rdivide(fg, dBinv(-Astop), omega);
  } else {
    emxInit_real_T(&x, 2);
    sigma = sqrt(FIR);
    i0 = x->size[0] * x->size[1];
    x->size[0] = 1;
    x->size[1] = fg->size[1];
    emxEnsureCapacity((emxArray__common *)x, i0, (int)sizeof(double));
    loop_ub = fg->size[0] * fg->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      x->data[i0] = sigma * fg->data[i0];
    }

    rdivide(x, dBinv(-Astop), omega);
    emxFree_real_T(&x);
  }

  sigma = dBinv(FIRdBmin);
  dec_int3 = omega->size[1];
  i0 = fg->size[0] * fg->size[1];
  fg->size[0] = 1;
  fg->size[1] = omega->size[1];
  emxEnsureCapacity((emxArray__common *)fg, i0, (int)sizeof(double));
  for (d_combinedResponse = 0; d_combinedResponse + 1 <= dec_int3;
       d_combinedResponse++) {
    if ((omega->data[d_combinedResponse] >= sigma) || rtIsNaN(sigma)) {
      sigmax = omega->data[d_combinedResponse];
    } else {
      sigmax = sigma;
    }

    fg->data[d_combinedResponse] = sigmax;
  }

  if (phEQ == -1.0) {
    b_abs(rgN, sw);
    i0 = rgN->size[0] * rgN->size[1];
    rgN->size[0] = 1;
    rgN->size[1] = sw->size[1];
    emxEnsureCapacity((emxArray__common *)rgN, i0, (int)sizeof(creal_T));
    loop_ub = sw->size[0] * sw->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      rgN->data[i0].re = sw->data[i0];
      rgN->data[i0].im = 0.0;
    }
  }

  emxInit_real_T(&weight, 2);
  rdivide(F3, dBinv(Apass / 2.0) - 1.0, sw);
  i0 = weight->size[0] * weight->size[1];
  weight->size[0] = 1;
  weight->size[1] = sw->size[1] + fg->size[1];
  emxEnsureCapacity((emxArray__common *)weight, i0, (int)sizeof(double));
  loop_ub = sw->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    weight->data[weight->size[0] * i0] = sw->data[sw->size[0] * i0];
  }

  loop_ub = fg->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    weight->data[weight->size[0] * (i0 + sw->size[1])] = fg->data[fg->size[0] *
      i0];
  }

  dec_int3 = 1;
  n = weight->size[1];
  sigma = weight->data[0];
  if (weight->size[1] > 1) {
    if (rtIsNaN(weight->data[0])) {
      d_combinedResponse = 2;
      exitg7 = false;
      while ((!exitg7) && (d_combinedResponse <= n)) {
        dec_int3 = d_combinedResponse;
        if (!rtIsNaN(weight->data[d_combinedResponse - 1])) {
          sigma = weight->data[d_combinedResponse - 1];
          exitg7 = true;
        } else {
          d_combinedResponse++;
        }
      }
    }

    if (dec_int3 < weight->size[1]) {
      while (dec_int3 + 1 <= n) {
        if (weight->data[dec_int3] > sigma) {
          sigma = weight->data[dec_int3];
        }

        dec_int3++;
      }
    }
  }

  emxInit_real_T(&b_weight, 2);
  i0 = b_weight->size[0] * b_weight->size[1];
  b_weight->size[0] = 1;
  b_weight->size[1] = weight->size[1];
  emxEnsureCapacity((emxArray__common *)b_weight, i0, (int)sizeof(double));
  loop_ub = weight->size[0] * weight->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    b_weight->data[i0] = weight->data[i0];
  }

  rdivide(b_weight, sigma, weight);

  /*  Set up design for FIR filter */
  i0 = fg->size[0] * fg->size[1];
  fg->size[0] = 1;
  fg->size[1] = rgN->size[1];
  emxEnsureCapacity((emxArray__common *)fg, i0, (int)sizeof(double));
  loop_ub = rgN->size[0] * rgN->size[1];
  emxFree_real_T(&b_weight);
  for (i0 = 0; i0 < loop_ub; i0++) {
    fg->data[i0] = rgN->data[i0].re;
  }

  if (1.0 > Gpass + 1.0) {
    loop_ub = 0;
  } else {
    loop_ub = (int)(Gpass + 1.0);
  }

  emxInit_real_T(&F1, 2);
  i0 = F1->size[0] * F1->size[1];
  F1->size[0] = 1;
  F1->size[1] = loop_ub;
  emxEnsureCapacity((emxArray__common *)F1, i0, (int)sizeof(double));
  for (i0 = 0; i0 < loop_ub; i0++) {
    F1->data[F1->size[0] * i0] = fg2->data[i0] * 2.0;
  }

  if (Gpass + 2.0 > fg2->size[1]) {
    i0 = 0;
    i2 = 0;
  } else {
    i0 = (int)(Gpass + 2.0) - 1;
    i2 = fg2->size[1];
  }

  emxInit_real_T(&F2, 2);
  i3 = F2->size[0] * F2->size[1];
  F2->size[0] = 1;
  F2->size[1] = i2 - i0;
  emxEnsureCapacity((emxArray__common *)F2, i3, (int)sizeof(double));
  loop_ub = i2 - i0;
  for (i2 = 0; i2 < loop_ub; i2++) {
    F2->data[F2->size[0] * i2] = fg2->data[i0 + i2] * 2.0;
  }

  if (1.0 > Gpass + 1.0) {
    loop_ub = 0;
  } else {
    loop_ub = (int)(Gpass + 1.0);
  }

  emxInit_real_T(&A1, 2);
  i0 = A1->size[0] * A1->size[1];
  A1->size[0] = 1;
  A1->size[1] = loop_ub;
  emxEnsureCapacity((emxArray__common *)A1, i0, (int)sizeof(double));
  for (i0 = 0; i0 < loop_ub; i0++) {
    A1->data[A1->size[0] * i0] = fg->data[i0];
  }

  if (Gpass + 2.0 > fg->size[1]) {
    i0 = 0;
    i2 = 0;
  } else {
    i0 = (int)(Gpass + 2.0) - 1;
    i2 = fg->size[1];
  }

  emxInit_real_T(&A2, 2);
  i3 = A2->size[0] * A2->size[1];
  A2->size[0] = 1;
  A2->size[1] = i2 - i0;
  emxEnsureCapacity((emxArray__common *)A2, i3, (int)sizeof(double));
  loop_ub = i2 - i0;
  for (i2 = 0; i2 < loop_ub; i2++) {
    A2->data[A2->size[0] * i2] = fg->data[i0 + i2];
  }

  if (1.0 > Gpass + 1.0) {
    loop_ub = 0;
  } else {
    loop_ub = (int)(Gpass + 1.0);
  }

  emxInit_real_T(&W1, 2);
  i0 = W1->size[0] * W1->size[1];
  W1->size[0] = 1;
  W1->size[1] = loop_ub;
  emxEnsureCapacity((emxArray__common *)W1, i0, (int)sizeof(double));
  for (i0 = 0; i0 < loop_ub; i0++) {
    W1->data[W1->size[0] * i0] = weight->data[i0];
  }

  if (Gpass + 2.0 > weight->size[1]) {
    i0 = 0;
    i2 = 0;
  } else {
    i0 = (int)(Gpass + 2.0) - 1;
    i2 = weight->size[1];
  }

  emxInit_real_T(&W2, 2);
  i3 = W2->size[0] * W2->size[1];
  W2->size[0] = 1;
  W2->size[1] = i2 - i0;
  emxEnsureCapacity((emxArray__common *)W2, i3, (int)sizeof(double));
  loop_ub = i2 - i0;
  for (i2 = 0; i2 < loop_ub; i2++) {
    W2->data[W2->size[0] * i2] = weight->data[i0 + i2];
  }

  /*  Determine the number of taps for FIR */
  if (b_strcmp(RxTx)) {
    if (hb3 == 1) {
      sigma = 16.0 * floor(converter_rate / Rdata);
      if (sigma <= 128.0) {
        apnd = sigma;
      } else {
        apnd = 128.0;
      }
    } else {
      sigma = 16.0 * floor(converter_rate / (2.0 * Rdata));
      if (sigma <= 128.0) {
        apnd = sigma;
      } else {
        apnd = 128.0;
      }
    }
  } else {
    switch ((int)FIR) {
     case 1:
      Nmax = 64;
      break;

     case 2:
      Nmax = 128;
      break;

     case 4:
      Nmax = 128;
      break;
    }

    sigma = 16.0 * floor(converter_rate * DAC_div / (2.0 * Rdata));
    if (sigma <= Nmax) {
      apnd = sigma;
    } else {
      apnd = Nmax;
    }
  }

  clkFIR = apnd / 16.0;
  loop_ub = (int)clkFIR * (int)apnd;
  for (i0 = 0; i0 < loop_ub; i0++) {
    tap_store_data[i0] = 0.0;
  }

  loop_ub = (int)(apnd / 16.0);
  for (i0 = 0; i0 < loop_ub; i0++) {
    Apass_actual_vector_data[i0] = 0.0;
  }

  loop_ub = (int)(apnd / 16.0);
  for (i0 = 0; i0 < loop_ub; i0++) {
    Astop_actual_vector_data[i0] = 0.0;
  }

  i = 0;

  /*  Design filter */
  emxInit_real_T(&ccoef, 2);
  emxInit_real_T(&F4, 2);
  emxInit_real_T(&h, 2);
  emxInit_creal_T(&b_rg, 2);
  emxInit_real_T(&d_omega2, 2);
  emxInit_creal_T(&b_rg2, 2);
  emxInit_real_T(&e_omega2, 2);
  emxInit_real_T(&f_omega2, 2);
  emxInit_real_T(&g_omega2, 2);
  emxInit_creal_T(&r1, 2);
  emxInit_real_T(&h_omega2, 2);
  emxInit_creal_T(&r2, 2);
  emxInit_real_T(&i_omega2, 2);
  emxInit_real_T(&j_omega2, 2);
  emxInit_real_T(&k_omega2, 2);
  emxInit_real_T(&b_omega, 2);
  emxInit_real_T(&b_F3, 2);
  emxInit_real_T(&b_sw, 2);
  emxInit_real_T(&b_A1, 2);
  emxInit_real_T(&b_F1, 2);
  emxInit_real_T(&b_W1, 2);
  emxInit_real_T(&c_A1, 2);
  emxInit_real_T(&c_F1, 2);
  emxInit_real_T(&c_W1, 2);
  do {
    exitg2 = 0;
    if (int_FIR != 0.0) {
      d_F1[0] = F1->data[0];
      d_F1[1] = F1->data[F1->size[1] - 1];
      d_F1[2] = F2->data[0];
      d_F1[3] = F2->data[F2->size[1] - 1];
      i0 = b_A1->size[0] * b_A1->size[1];
      b_A1->size[0] = 1;
      b_A1->size[1] = A1->size[1] + A2->size[1];
      emxEnsureCapacity((emxArray__common *)b_A1, i0, (int)sizeof(double));
      loop_ub = A1->size[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        b_A1->data[b_A1->size[0] * i0] = A1->data[A1->size[0] * i0];
      }

      loop_ub = A2->size[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        b_A1->data[b_A1->size[0] * (i0 + A1->size[1])] = A2->data[A2->size[0] *
          i0];
      }

      i0 = b_F1->size[0] * b_F1->size[1];
      b_F1->size[0] = 1;
      b_F1->size[1] = F1->size[1] + F2->size[1];
      emxEnsureCapacity((emxArray__common *)b_F1, i0, (int)sizeof(double));
      loop_ub = F1->size[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        b_F1->data[b_F1->size[0] * i0] = F1->data[F1->size[0] * i0];
      }

      loop_ub = F2->size[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        b_F1->data[b_F1->size[0] * (i0 + F1->size[1])] = F2->data[F2->size[0] *
          i0];
      }

      i0 = b_W1->size[0] * b_W1->size[1];
      b_W1->size[0] = 1;
      b_W1->size[1] = W1->size[1] + W2->size[1];
      emxEnsureCapacity((emxArray__common *)b_W1, i0, (int)sizeof(double));
      loop_ub = W1->size[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        b_W1->data[b_W1->size[0] * i0] = W1->data[W1->size[0] * i0];
      }

      loop_ub = W2->size[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        b_W1->data[b_W1->size[0] * (i0 + W1->size[1])] = W2->data[W2->size[0] *
          i0];
      }

      firpm_cg(apnd - 1.0, d_F1, b_A1, b_F1, b_W1, ccoef);
    } else {
      /*  Check different designs until we reach required ripple condition */
      sigma = rt_powd_snf(10.0, -Astop / 20.0);

      /*  Peak Ripple */
      i0 = ccoef->size[0] * ccoef->size[1];
      ccoef->size[0] = 1;
      ccoef->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)ccoef, i0, (int)sizeof(double));
      ccoef->data[0] = 0.0;

      /*  Predef type */
      d_combinedResponse = 0;
      exitg6 = false;
      while ((!exitg6) && (d_combinedResponse < 126)) {
        e_F1[0] = F1->data[0];
        e_F1[1] = F1->data[F1->size[1] - 1];
        e_F1[2] = F2->data[0];
        e_F1[3] = F2->data[F2->size[1] - 1];
        i0 = c_A1->size[0] * c_A1->size[1];
        c_A1->size[0] = 1;
        c_A1->size[1] = A1->size[1] + A2->size[1];
        emxEnsureCapacity((emxArray__common *)c_A1, i0, (int)sizeof(double));
        loop_ub = A1->size[1];
        for (i0 = 0; i0 < loop_ub; i0++) {
          c_A1->data[c_A1->size[0] * i0] = A1->data[A1->size[0] * i0];
        }

        loop_ub = A2->size[1];
        for (i0 = 0; i0 < loop_ub; i0++) {
          c_A1->data[c_A1->size[0] * (i0 + A1->size[1])] = A2->data[A2->size[0] *
            i0];
        }

        i0 = c_F1->size[0] * c_F1->size[1];
        c_F1->size[0] = 1;
        c_F1->size[1] = F1->size[1] + F2->size[1];
        emxEnsureCapacity((emxArray__common *)c_F1, i0, (int)sizeof(double));
        loop_ub = F1->size[1];
        for (i0 = 0; i0 < loop_ub; i0++) {
          c_F1->data[c_F1->size[0] * i0] = F1->data[F1->size[0] * i0];
        }

        loop_ub = F2->size[1];
        for (i0 = 0; i0 < loop_ub; i0++) {
          c_F1->data[c_F1->size[0] * (i0 + F1->size[1])] = F2->data[F2->size[0] *
            i0];
        }

        i0 = c_W1->size[0] * c_W1->size[1];
        c_W1->size[0] = 1;
        c_W1->size[1] = W1->size[1] + W2->size[1];
        emxEnsureCapacity((emxArray__common *)c_W1, i0, (int)sizeof(double));
        loop_ub = W1->size[1];
        for (i0 = 0; i0 < loop_ub; i0++) {
          c_W1->data[c_W1->size[0] * i0] = W1->data[W1->size[0] * i0];
        }

        loop_ub = W2->size[1];
        for (i0 = 0; i0 < loop_ub; i0++) {
          c_W1->data[c_W1->size[0] * (i0 + W1->size[1])] = W2->data[W2->size[0] *
            i0];
        }

        b_firpm_cg(3.0 + (double)d_combinedResponse, e_F1, c_A1, c_F1, c_W1,
                   ccoef, &valid, &sigmax);

        /*  Check if design meets specs */
        if ((sigmax < sigma) && valid) {
          exitg6 = true;
        } else {
          d_combinedResponse++;
        }
      }
    }

    /*  Enable phase equalization and apply update to taps */
    if (phEQ != -1.0) {
      if (1 > fg2->size[1]) {
        i0 = 1;
        i2 = 1;
        i3 = 0;
      } else {
        i0 = fg2->size[1];
        i2 = -1;
        i3 = 1;
      }

      d_combinedResponse = fg->size[0] * fg->size[1];
      fg->size[0] = 1;
      fg->size[1] = div_s32_floor(i3 - i0, i2) + 1;
      emxEnsureCapacity((emxArray__common *)fg, d_combinedResponse, (int)sizeof
                        (double));
      loop_ub = div_s32_floor(i3 - i0, i2);
      for (i3 = 0; i3 <= loop_ub; i3++) {
        fg->data[fg->size[0] * i3] = 0.5 - fg2->data[(i0 + i2 * i3) - 1];
      }

      if (1 > rgN->size[1]) {
        i0 = 1;
        i2 = 1;
        i3 = 0;
      } else {
        i0 = rgN->size[1];
        i2 = -1;
        i3 = 1;
      }

      d_combinedResponse = omega->size[0] * omega->size[1];
      omega->size[0] = 1;
      omega->size[1] = div_s32_floor(i3 - i0, i2) + 1;
      emxEnsureCapacity((emxArray__common *)omega, d_combinedResponse, (int)
                        sizeof(double));
      loop_ub = div_s32_floor(i3 - i0, i2);
      for (i3 = 0; i3 <= loop_ub; i3++) {
        omega->data[omega->size[0] * i3] = rgN->data[(i0 + i2 * i3) - 1].im;
      }

      if (1 > weight->size[1]) {
        i0 = 1;
        i2 = 1;
        i3 = 0;
      } else {
        i0 = weight->size[1];
        i2 = -1;
        i3 = 1;
      }

      d_combinedResponse = sw->size[0] * sw->size[1];
      sw->size[0] = 1;
      sw->size[1] = div_s32_floor(i3 - i0, i2) + 1;
      emxEnsureCapacity((emxArray__common *)sw, d_combinedResponse, (int)sizeof
                        (double));
      loop_ub = div_s32_floor(i3 - i0, i2);
      for (d_combinedResponse = 0; d_combinedResponse <= loop_ub;
           d_combinedResponse++) {
        sw->data[sw->size[0] * d_combinedResponse] = weight->data[(i0 + i2 *
          d_combinedResponse) - 1];
      }

      if (1.0 > (8192.0 - Gstop) + 1.0) {
        loop_ub = 0;
      } else {
        loop_ub = (int)((8192.0 - Gstop) + 1.0);
      }

      d_combinedResponse = F3->size[0] * F3->size[1];
      F3->size[0] = 1;
      F3->size[1] = loop_ub;
      emxEnsureCapacity((emxArray__common *)F3, d_combinedResponse, (int)sizeof
                        (double));
      for (d_combinedResponse = 0; d_combinedResponse < loop_ub;
           d_combinedResponse++) {
        F3->data[F3->size[0] * d_combinedResponse] = fg->data[d_combinedResponse]
          * 2.0;
      }

      if ((8192.0 - Gstop) + 2.0 > fg->size[1]) {
        d_combinedResponse = 0;
        dec_int3 = 0;
      } else {
        d_combinedResponse = (int)((8192.0 - Gstop) + 2.0) - 1;
        dec_int3 = fg->size[1];
      }

      n = F4->size[0] * F4->size[1];
      F4->size[0] = 1;
      F4->size[1] = dec_int3 - d_combinedResponse;
      emxEnsureCapacity((emxArray__common *)F4, n, (int)sizeof(double));
      loop_ub = dec_int3 - d_combinedResponse;
      for (dec_int3 = 0; dec_int3 < loop_ub; dec_int3++) {
        F4->data[F4->size[0] * dec_int3] = fg->data[d_combinedResponse +
          dec_int3] * 2.0;
      }

      if (1.0 > (8192.0 - Gstop) + 1.0) {
        loop_ub = 0;
      } else {
        loop_ub = (int)((8192.0 - Gstop) + 1.0);
      }

      if ((8192.0 - Gstop) + 2.0 > omega->size[1]) {
        d_combinedResponse = 0;
        dec_int3 = 0;
      } else {
        d_combinedResponse = (int)((8192.0 - Gstop) + 2.0) - 1;
        dec_int3 = omega->size[1];
      }

      if (1.0 > (8192.0 - Gstop) + 1.0) {
        hb3 = 0;
      } else {
        hb3 = (int)((8192.0 - Gstop) + 1.0);
      }

      if ((8192.0 - Gstop) + 2.0 > div_s32_floor(i3 - i0, i2) + 1) {
        n = 0;
        i0 = -1;
      } else {
        n = (int)((8192.0 - Gstop) + 2.0) - 1;
        i0 = div_s32_floor(i3 - i0, i2);
      }

      if (int_FIR != 0.0) {
        sigma = apnd - 1.0;
      } else {
        sigma = (double)ccoef->size[1] - 1.0;
      }

      c_F3[0] = F3->data[0];
      c_F3[1] = F3->data[F3->size[1] - 1];
      c_F3[2] = F4->data[0];
      c_F3[3] = F4->data[F4->size[1] - 1];
      i2 = b_omega->size[0] * b_omega->size[1];
      b_omega->size[0] = 1;
      b_omega->size[1] = (loop_ub + dec_int3) - d_combinedResponse;
      emxEnsureCapacity((emxArray__common *)b_omega, i2, (int)sizeof(double));
      for (i2 = 0; i2 < loop_ub; i2++) {
        b_omega->data[b_omega->size[0] * i2] = omega->data[i2];
      }

      dec_int3 -= d_combinedResponse;
      for (i2 = 0; i2 < dec_int3; i2++) {
        b_omega->data[b_omega->size[0] * (i2 + loop_ub)] = omega->
          data[d_combinedResponse + i2];
      }

      i2 = b_F3->size[0] * b_F3->size[1];
      b_F3->size[0] = 1;
      b_F3->size[1] = F3->size[1] + F4->size[1];
      emxEnsureCapacity((emxArray__common *)b_F3, i2, (int)sizeof(double));
      loop_ub = F3->size[1];
      for (i2 = 0; i2 < loop_ub; i2++) {
        b_F3->data[b_F3->size[0] * i2] = F3->data[F3->size[0] * i2];
      }

      loop_ub = F4->size[1];
      for (i2 = 0; i2 < loop_ub; i2++) {
        b_F3->data[b_F3->size[0] * (i2 + F3->size[1])] = F4->data[F4->size[0] *
          i2];
      }

      i2 = b_sw->size[0] * b_sw->size[1];
      b_sw->size[0] = 1;
      b_sw->size[1] = ((hb3 + i0) - n) + 1;
      emxEnsureCapacity((emxArray__common *)b_sw, i2, (int)sizeof(double));
      for (i2 = 0; i2 < hb3; i2++) {
        b_sw->data[b_sw->size[0] * i2] = sw->data[i2];
      }

      loop_ub = i0 - n;
      for (i0 = 0; i0 <= loop_ub; i0++) {
        b_sw->data[b_sw->size[0] * (i0 + hb3)] = sw->data[n + i0];
      }

      firpm_cg(sigma, c_F3, b_omega, b_F3, b_sw, fg);
      i0 = fg->size[1];
      for (d_combinedResponse = 0; d_combinedResponse < i0; d_combinedResponse++)
      {
        fg->data[d_combinedResponse] = -fg->data[d_combinedResponse] *
          rt_powd_snf(-1.0, (1.0 + (double)d_combinedResponse) - 1.0);
      }
    } else {
      for (i0 = 0; i0 < 2; i0++) {
        b2[i0] = ccoef->size[i0];
      }

      i0 = fg->size[0] * fg->size[1];
      fg->size[0] = 1;
      fg->size[1] = (int)b2[1];
      emxEnsureCapacity((emxArray__common *)fg, i0, (int)sizeof(double));
      loop_ub = (int)b2[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        fg->data[i0] = 0.0;
      }
    }

    loop_ub = ccoef->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      tap_store_data[i + (int)clkFIR * i0] = ccoef->data[ccoef->size[0] * i0] +
        fg->data[fg->size[0] * i0];
    }

    /*  scoef ==0 when no EQ */
    /*  TODO: Set fractionalLength based on FSR of input */
    if (1 > ccoef->size[1]) {
      loop_ub = 0;
    } else {
      loop_ub = ccoef->size[1];
    }

    for (i0 = 0; i0 < loop_ub; i0++) {
      sigma = tap_store_data[i + (int)clkFIR * i0] * 65536.0;
      sigmax = fabs(sigma);
      if (sigmax < 4.503599627370496E+15) {
        if (sigmax >= 0.5) {
          sigma = floor(sigma + 0.5);
        } else {
          sigma *= 0.0;
        }
      }

      if (sigma < 32768.0) {
        if (sigma >= -32768.0) {
          i4 = (short)sigma;
        } else {
          i4 = MIN_int16_T;
        }
      } else if (sigma >= 32768.0) {
        i4 = MAX_int16_T;
      } else {
        i4 = 0;
      }

      b_tap_store_data[i0] = (double)i4 * 1.52587890625E-5;
    }

    for (i0 = 0; i0 < loop_ub; i0++) {
      tap_store_data[i + (int)clkFIR * i0] = b_tap_store_data[i0];
    }

    if (b_strcmp(RxTx)) {
      if (1.0 > Gpass + 1.0) {
        loop_ub = 0;
      } else {
        loop_ub = (int)(Gpass + 1.0);
      }

      if (1 > ccoef->size[1]) {
        hb3 = 0;
      } else {
        hb3 = ccoef->size[1];
      }

      i0 = k_omega2->size[0] * k_omega2->size[1];
      k_omega2->size[0] = 1;
      k_omega2->size[1] = loop_ub;
      emxEnsureCapacity((emxArray__common *)k_omega2, i0, (int)sizeof(double));
      for (i0 = 0; i0 < loop_ub; i0++) {
        k_omega2->data[k_omega2->size[0] * i0] = omega2->data[i0];
      }

      tap_store_size[0] = 1;
      tap_store_size[1] = hb3;
      for (i0 = 0; i0 < hb3; i0++) {
        b_tap_store_data[i0] = tap_store_data[i + (int)clkFIR * i0];
      }

      c_generateCascadedResponseRx(enables, k_omega2, converter_rate, hb1_coeff,
        hb2_coeff, hb3_coeff_data, hb3_coeff_size, dec_int3_coeff_data,
        dec_int3_coeff_size, b_tap_store_data, tap_store_size, rg2);
      if (Gpass + 2.0 > omega2->size[1]) {
        i0 = 0;
        i2 = 0;
      } else {
        i0 = (int)(Gpass + 2.0) - 1;
        i2 = omega2->size[1];
      }

      if (1 > ccoef->size[1]) {
        loop_ub = 0;
      } else {
        loop_ub = ccoef->size[1];
      }

      i3 = j_omega2->size[0] * j_omega2->size[1];
      j_omega2->size[0] = 1;
      j_omega2->size[1] = i2 - i0;
      emxEnsureCapacity((emxArray__common *)j_omega2, i3, (int)sizeof(double));
      hb3 = i2 - i0;
      for (i2 = 0; i2 < hb3; i2++) {
        j_omega2->data[j_omega2->size[0] * i2] = omega2->data[i0 + i2];
      }

      c_tap_store_size[0] = 1;
      c_tap_store_size[1] = loop_ub;
      for (i0 = 0; i0 < loop_ub; i0++) {
        b_tap_store_data[i0] = tap_store_data[i + (int)clkFIR * i0];
      }

      c_generateCascadedResponseRx(enables, j_omega2, converter_rate, hb1_coeff,
        hb2_coeff, hb3_coeff_data, hb3_coeff_size, dec_int3_coeff_data,
        dec_int3_coeff_size, b_tap_store_data, c_tap_store_size, rg);
      if (1.0 > Gpass + 1.0) {
        loop_ub = 0;
      } else {
        loop_ub = (int)(Gpass + 1.0);
      }

      i0 = i_omega2->size[0] * i_omega2->size[1];
      i_omega2->size[0] = 1;
      i_omega2->size[1] = loop_ub;
      emxEnsureCapacity((emxArray__common *)i_omega2, i0, (int)sizeof(double));
      for (i0 = 0; i0 < loop_ub; i0++) {
        i_omega2->data[i_omega2->size[0] * i0] = omega2->data[i0];
      }

      c_analogresp(i_omega2, converter_rate, b1_data, b1_size, a1, b2_data,
                   b2_size, a2, r0);
      i0 = r2->size[0] * r2->size[1];
      r2->size[0] = 1;
      r2->size[1] = r0->size[1];
      emxEnsureCapacity((emxArray__common *)r2, i0, (int)sizeof(creal_T));
      loop_ub = r0->size[0] * r0->size[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        absa = r0->data[i0].re;
        im = r0->data[i0].im;
        sigma = rg2->data[i0].re;
        sigmax = rg2->data[i0].im;
        r2->data[i0].re = absa * sigma - im * sigmax;
        r2->data[i0].im = absa * sigmax + im * sigma;
      }

      b_abs(r2, fg);
      if (Gpass + 2.0 > omega2->size[1]) {
        i0 = 0;
        i2 = 0;
      } else {
        i0 = (int)(Gpass + 2.0) - 1;
        i2 = omega2->size[1];
      }

      i3 = h_omega2->size[0] * h_omega2->size[1];
      h_omega2->size[0] = 1;
      h_omega2->size[1] = i2 - i0;
      emxEnsureCapacity((emxArray__common *)h_omega2, i3, (int)sizeof(double));
      loop_ub = i2 - i0;
      for (i2 = 0; i2 < loop_ub; i2++) {
        h_omega2->data[h_omega2->size[0] * i2] = omega2->data[i0 + i2];
      }

      c_analogresp(h_omega2, converter_rate, b1_data, b1_size, a1, b2_data,
                   b2_size, a2, r0);
      i0 = r1->size[0] * r1->size[1];
      r1->size[0] = 1;
      r1->size[1] = r0->size[1];
      emxEnsureCapacity((emxArray__common *)r1, i0, (int)sizeof(creal_T));
      loop_ub = r0->size[0] * r0->size[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        absa = r0->data[i0].re;
        im = r0->data[i0].im;
        sigma = rg->data[i0].re;
        sigmax = rg->data[i0].im;
        r1->data[i0].re = absa * sigma - im * sigmax;
        r1->data[i0].im = absa * sigmax + im * sigma;
      }

      b_abs(r1, omega);
    } else {
      /*  TX */
      if (1.0 > Gpass + 1.0) {
        loop_ub = 0;
      } else {
        loop_ub = (int)(Gpass + 1.0);
      }

      if (1 > ccoef->size[1]) {
        hb3 = 0;
      } else {
        hb3 = ccoef->size[1];
      }

      i0 = g_omega2->size[0] * g_omega2->size[1];
      g_omega2->size[0] = 1;
      g_omega2->size[1] = loop_ub;
      emxEnsureCapacity((emxArray__common *)g_omega2, i0, (int)sizeof(double));
      for (i0 = 0; i0 < loop_ub; i0++) {
        g_omega2->data[g_omega2->size[0] * i0] = omega2->data[i0];
      }

      b_tap_store_size[0] = 1;
      b_tap_store_size[1] = hb3;
      for (i0 = 0; i0 < hb3; i0++) {
        b_tap_store_data[i0] = tap_store_data[i + (int)clkFIR * i0];
      }

      c_generateCascadedResponseRx(enables, g_omega2, converter_rate, hb1_coeff,
        hb2_coeff, hb3_coeff_data, hb3_coeff_size, dec_int3_coeff_data,
        dec_int3_coeff_size, b_tap_store_data, b_tap_store_size, rg2);
      if (Gpass + 2.0 > omega2->size[1]) {
        i0 = 0;
        i2 = 0;
      } else {
        i0 = (int)(Gpass + 2.0) - 1;
        i2 = omega2->size[1];
      }

      if (1 > ccoef->size[1]) {
        loop_ub = 0;
      } else {
        loop_ub = ccoef->size[1];
      }

      i3 = f_omega2->size[0] * f_omega2->size[1];
      f_omega2->size[0] = 1;
      f_omega2->size[1] = i2 - i0;
      emxEnsureCapacity((emxArray__common *)f_omega2, i3, (int)sizeof(double));
      hb3 = i2 - i0;
      for (i2 = 0; i2 < hb3; i2++) {
        f_omega2->data[f_omega2->size[0] * i2] = omega2->data[i0 + i2];
      }

      d_tap_store_size[0] = 1;
      d_tap_store_size[1] = loop_ub;
      for (i0 = 0; i0 < loop_ub; i0++) {
        b_tap_store_data[i0] = tap_store_data[i + (int)clkFIR * i0];
      }

      c_generateCascadedResponseRx(enables, f_omega2, converter_rate, hb1_coeff,
        hb2_coeff, hb3_coeff_data, hb3_coeff_size, dec_int3_coeff_data,
        dec_int3_coeff_size, b_tap_store_data, d_tap_store_size, rg);
      if (1.0 > Gpass + 1.0) {
        loop_ub = 0;
      } else {
        loop_ub = (int)(Gpass + 1.0);
      }

      i0 = e_omega2->size[0] * e_omega2->size[1];
      e_omega2->size[0] = 1;
      e_omega2->size[1] = loop_ub;
      emxEnsureCapacity((emxArray__common *)e_omega2, i0, (int)sizeof(double));
      for (i0 = 0; i0 < loop_ub; i0++) {
        e_omega2->data[e_omega2->size[0] * i0] = omega2->data[i0];
      }

      d_analogresp(e_omega2, converter_rate, b1_data, b1_size, a1, b2_data,
                   b2_size, a2, r0);
      i0 = b_rg2->size[0] * b_rg2->size[1];
      b_rg2->size[0] = 1;
      b_rg2->size[1] = rg2->size[1];
      emxEnsureCapacity((emxArray__common *)b_rg2, i0, (int)sizeof(creal_T));
      loop_ub = rg2->size[0] * rg2->size[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        sigma = rg2->data[i0].re;
        sigmax = rg2->data[i0].im;
        absa = r0->data[i0].re;
        im = r0->data[i0].im;
        b_rg2->data[i0].re = sigma * absa - sigmax * im;
        b_rg2->data[i0].im = sigma * im + sigmax * absa;
      }

      b_abs(b_rg2, fg);
      if (Gpass + 2.0 > omega2->size[1]) {
        i0 = 0;
        i2 = 0;
      } else {
        i0 = (int)(Gpass + 2.0) - 1;
        i2 = omega2->size[1];
      }

      i3 = d_omega2->size[0] * d_omega2->size[1];
      d_omega2->size[0] = 1;
      d_omega2->size[1] = i2 - i0;
      emxEnsureCapacity((emxArray__common *)d_omega2, i3, (int)sizeof(double));
      loop_ub = i2 - i0;
      for (i2 = 0; i2 < loop_ub; i2++) {
        d_omega2->data[d_omega2->size[0] * i2] = omega2->data[i0 + i2];
      }

      d_analogresp(d_omega2, converter_rate, b1_data, b1_size, a1, b2_data,
                   b2_size, a2, r0);
      i0 = b_rg->size[0] * b_rg->size[1];
      b_rg->size[0] = 1;
      b_rg->size[1] = rg->size[1];
      emxEnsureCapacity((emxArray__common *)b_rg, i0, (int)sizeof(creal_T));
      loop_ub = rg->size[0] * rg->size[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        sigma = rg->data[i0].re;
        sigmax = rg->data[i0].im;
        absa = r0->data[i0].re;
        im = r0->data[i0].im;
        b_rg->data[i0].re = sigma * absa - sigmax * im;
        b_rg->data[i0].im = sigma * im + sigmax * absa;
      }

      b_abs(b_rg, omega);
    }

    /*  quantitative values about actual passband and stopband */
    dec_int3 = 1;
    n = fg->size[1];
    sigma = fg->data[0];
    if (fg->size[1] > 1) {
      if (rtIsNaN(fg->data[0])) {
        d_combinedResponse = 2;
        exitg5 = false;
        while ((!exitg5) && (d_combinedResponse <= n)) {
          dec_int3 = d_combinedResponse;
          if (!rtIsNaN(fg->data[d_combinedResponse - 1])) {
            sigma = fg->data[d_combinedResponse - 1];
            exitg5 = true;
          } else {
            d_combinedResponse++;
          }
        }
      }

      if (dec_int3 < fg->size[1]) {
        while (dec_int3 + 1 <= n) {
          if (fg->data[dec_int3] > sigma) {
            sigma = fg->data[dec_int3];
          }

          dec_int3++;
        }
      }
    }

    dec_int3 = 1;
    n = fg->size[1];
    sigmax = fg->data[0];
    if (fg->size[1] > 1) {
      if (rtIsNaN(fg->data[0])) {
        d_combinedResponse = 2;
        exitg4 = false;
        while ((!exitg4) && (d_combinedResponse <= n)) {
          dec_int3 = d_combinedResponse;
          if (!rtIsNaN(fg->data[d_combinedResponse - 1])) {
            sigmax = fg->data[d_combinedResponse - 1];
            exitg4 = true;
          } else {
            d_combinedResponse++;
          }
        }
      }

      if (dec_int3 < fg->size[1]) {
        while (dec_int3 + 1 <= n) {
          if (fg->data[dec_int3] < sigmax) {
            sigmax = fg->data[dec_int3];
          }

          dec_int3++;
        }
      }
    }

    Apass_actual_vector_data[i] = mag2db(sigma) - mag2db(sigmax);
    dec_int3 = 1;
    n = omega->size[1];
    sigma = omega->data[0];
    if (omega->size[1] > 1) {
      if (rtIsNaN(omega->data[0])) {
        d_combinedResponse = 2;
        exitg3 = false;
        while ((!exitg3) && (d_combinedResponse <= n)) {
          dec_int3 = d_combinedResponse;
          if (!rtIsNaN(omega->data[d_combinedResponse - 1])) {
            sigma = omega->data[d_combinedResponse - 1];
            exitg3 = true;
          } else {
            d_combinedResponse++;
          }
        }
      }

      if (dec_int3 < omega->size[1]) {
        while (dec_int3 + 1 <= n) {
          if (omega->data[dec_int3] > sigma) {
            sigma = omega->data[dec_int3];
          }

          dec_int3++;
        }
      }
    }

    Astop_actual_vector_data[i] = -mag2db(sigma);
    if (int_FIR == 0.0) {
      if (1 > ccoef->size[1]) {
        loop_ub = 0;
      } else {
        loop_ub = ccoef->size[1];
      }

      i0 = h->size[0] * h->size[1];
      h->size[0] = 1;
      h->size[1] = loop_ub;
      emxEnsureCapacity((emxArray__common *)h, i0, (int)sizeof(double));
      for (i0 = 0; i0 < loop_ub; i0++) {
        h->data[h->size[0] * i0] = tap_store_data[(int)clkFIR * i0];
      }

      /* Apass_actual = Apass_actual_vector(1); */
      /* Astop_actual = Astop_actual_vector(1); */
      exitg2 = 1;
    } else if ((Apass_actual_vector_data[0] > Apass) ||
               (Astop_actual_vector_data[0] < Astop)) {
      if (1.0 > apnd) {
        loop_ub = 0;
      } else {
        loop_ub = (int)apnd;
      }

      i0 = h->size[0] * h->size[1];
      h->size[0] = 1;
      h->size[1] = loop_ub;
      emxEnsureCapacity((emxArray__common *)h, i0, (int)sizeof(double));
      for (i0 = 0; i0 < loop_ub; i0++) {
        h->data[h->size[0] * i0] = tap_store_data[(int)clkFIR * i0];
      }

      /* Apass_actual = Apass_actual_vector(1); */
      /* Astop_actual = Astop_actual_vector(1); */
      exitg2 = 1;
    } else if ((Apass_actual_vector_data[i] > Apass) ||
               (Astop_actual_vector_data[i] < Astop)) {
      if (1.0 > apnd + 16.0) {
        loop_ub = 0;
      } else {
        loop_ub = (int)(apnd + 16.0);
      }

      i0 = h->size[0] * h->size[1];
      h->size[0] = 1;
      h->size[1] = loop_ub;
      emxEnsureCapacity((emxArray__common *)h, i0, (int)sizeof(double));
      for (i0 = 0; i0 < loop_ub; i0++) {
        h->data[h->size[0] * i0] = tap_store_data[(i + (int)clkFIR * i0) - 1];
      }

      /* Apass_actual = Apass_actual_vector(i-1); */
      /* Astop_actual = Astop_actual_vector(i-1); */
      exitg2 = 1;
    } else {
      apnd -= 16.0;
      i++;
    }
  } while (exitg2 == 0);

  emxFree_real_T(&c_W1);
  emxFree_real_T(&c_F1);
  emxFree_real_T(&c_A1);
  emxFree_real_T(&b_W1);
  emxFree_real_T(&b_F1);
  emxFree_real_T(&b_A1);
  emxFree_real_T(&b_sw);
  emxFree_real_T(&b_F3);
  emxFree_real_T(&b_omega);
  emxFree_real_T(&k_omega2);
  emxFree_real_T(&j_omega2);
  emxFree_real_T(&i_omega2);
  emxFree_creal_T(&r2);
  emxFree_real_T(&h_omega2);
  emxFree_creal_T(&r1);
  emxFree_real_T(&g_omega2);
  emxFree_real_T(&f_omega2);
  emxFree_real_T(&e_omega2);
  emxFree_creal_T(&b_rg2);
  emxFree_real_T(&d_omega2);
  emxFree_creal_T(&b_rg);
  emxFree_creal_T(&r0);
  emxFree_real_T(&F4);
  emxFree_real_T(&F3);
  emxFree_real_T(&sw);
  emxFree_real_T(&ccoef);
  emxFree_real_T(&W2);
  emxFree_real_T(&W1);
  emxFree_real_T(&A2);
  emxFree_real_T(&A1);
  emxFree_real_T(&F2);
  emxFree_real_T(&F1);
  emxFree_real_T(&weight);
  emxFree_creal_T(&rgN);
  emxFree_real_T(&omega2);
  emxFree_real_T(&fg2);
  emxFree_creal_T(&rg);
  emxFree_creal_T(&rg2);
  emxFree_real_T(&omega);
  emxFree_real_T(&fg);
  emxFree_creal_T(&a2);
  emxFree_creal_T(&a1);
  if (c_strcmp(RxTx)) {
    emxInit_real_T(&r3, 2);
    emxInit_real_T(&r4, 2);
    if ((int_FIR == 1.0) && (FIR == 2.0)) {
      if (rt_remd_snf(h->size[1], 32.0) != 0.0) {
        i0 = r4->size[0] * r4->size[1];
        r4->size[0] = 1;
        r4->size[1] = 16 + h->size[1];
        emxEnsureCapacity((emxArray__common *)r4, i0, (int)sizeof(double));
        for (i0 = 0; i0 < 8; i0++) {
          r4->data[r4->size[0] * i0] = 0.0;
        }

        loop_ub = h->size[1];
        for (i0 = 0; i0 < loop_ub; i0++) {
          r4->data[r4->size[0] * (i0 + 8)] = h->data[h->size[0] * i0];
        }

        for (i0 = 0; i0 < 8; i0++) {
          r4->data[r4->size[0] * ((i0 + h->size[1]) + 8)] = 0.0;
        }

        i0 = h->size[0] * h->size[1];
        h->size[0] = 1;
        h->size[1] = r4->size[1];
        emxEnsureCapacity((emxArray__common *)h, i0, (int)sizeof(double));
        loop_ub = r4->size[1];
        for (i0 = 0; i0 < loop_ub; i0++) {
          h->data[h->size[0] * i0] = r4->data[r4->size[0] * i0];
        }
      }
    } else {
      if ((int_FIR == 1.0) && (FIR == 4.0) && (rt_remd_snf(h->size[1], 64.0) !=
           0.0)) {
        sigma = (ceil((double)h->size[1] / 64.0) * 64.0 - (double)h->size[1]) /
          2.0;
        i0 = r3->size[0] * r3->size[1];
        r3->size[0] = 1;
        r3->size[1] = ((int)sigma + h->size[1]) + (int)sigma;
        emxEnsureCapacity((emxArray__common *)r3, i0, (int)sizeof(double));
        loop_ub = (int)sigma;
        for (i0 = 0; i0 < loop_ub; i0++) {
          r3->data[r3->size[0] * i0] = 0.0;
        }

        loop_ub = h->size[1];
        for (i0 = 0; i0 < loop_ub; i0++) {
          r3->data[r3->size[0] * (i0 + (int)sigma)] = h->data[h->size[0] * i0];
        }

        loop_ub = (int)sigma;
        for (i0 = 0; i0 < loop_ub; i0++) {
          r3->data[r3->size[0] * ((i0 + (int)sigma) + h->size[1])] = 0.0;
        }

        i0 = h->size[0] * h->size[1];
        h->size[0] = 1;
        h->size[1] = r3->size[1];
        emxEnsureCapacity((emxArray__common *)h, i0, (int)sizeof(double));
        loop_ub = r3->size[1];
        for (i0 = 0; i0 < loop_ub; i0++) {
          h->data[h->size[0] * i0] = r3->data[r3->size[0] * i0];
        }
      }
    }

    emxFree_real_T(&r4);
    emxFree_real_T(&r3);
  }

  /*  There will always be 128 taps output */
  memset(&firTapsPreScale[0], 0, sizeof(double) << 7);
  loop_ub = h->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    firTapsPreScale[i0] = h->data[h->size[0] * i0];
  }

  emxFree_real_T(&h);

  /*  Calculate group delay */
  /*  Hmd = dec_int_func(input.FIR,h(1:128)); */
  /*  */
  /*  if ~isempty(ver('fixedpoint')) && license('test','fixed_point_toolbox') && license('checkout','fixed_point_toolbox') */
  /*      Hmd.Numerator = double(fi(Hmd.Numerator,true,16)); */
  /*  end */
  /*  if strcmp(input.RxTx, 'Rx') */
  /*      addStage(dfilter, Hmd); */
  /*  else */
  /*      addStage(dfilter, Hmd, 1); */
  /*  end */
  /* gd2c = grpdelay(Hmd,omega1,clkFIR).*(1/clkFIR); */
  /*  gd2 = grpdelay_cg(firTapsPreScale,1,omega1,clkFIR).'.*(1/clkFIR); */
  /*  */
  /*  if input.phEQ == -1 */
  /*      groupdelay = gd1 + gd2; */
  /*  else */
  /*      groupdelay = gd1 + gd2; */
  /*  end */
  /*  grpdelayvar = max(groupdelay)-min(groupdelay); */
  /*  Determine Gains */
  dec_int3 = 1;
  sigma = firTapsPreScale[0];
  if (rtIsNaN(firTapsPreScale[0])) {
    d_combinedResponse = 2;
    exitg1 = false;
    while ((!exitg1) && (d_combinedResponse < 129)) {
      dec_int3 = d_combinedResponse;
      if (!rtIsNaN(firTapsPreScale[d_combinedResponse - 1])) {
        sigma = firTapsPreScale[d_combinedResponse - 1];
        exitg1 = true;
      } else {
        d_combinedResponse++;
      }
    }
  }

  if (dec_int3 < 128) {
    while (dec_int3 + 1 < 129) {
      if (firTapsPreScale[dec_int3] > sigma) {
        sigma = firTapsPreScale[dec_int3];
      }

      dec_int3++;
    }
  }

  /*  switch aTFIR */
  /*      case 2 */
  /*          gain = 6; */
  /*      case 1 */
  /*          gain = 0; */
  /*      case 0 */
  /*          gain = -6; */
  /*      otherwise */
  /*          gain = -12; */
  /*  end */
  /*   */
  /*  if strcmp(input.RxTx, 'Rx') */
  /*      if aTFIR > 2 */
  /*          gain = 6; */
  /*      end */
  /*  else */
  /*      if input.FIR == 2 */
  /*          gain = gain+6; */
  /*      elseif input.FIR == 4 */
  /*          gain = gain+12; */
  /*      end */
  /*      if gain > 0 */
  /*          gain = 0; */
  /*      elseif gain < -6 */
  /*          gain = -6; */
  /*      end */
  /*  end */
  /*  Scale taps */
  sigmax = rt_powd_snf(2.0, 16.0 - (1.0 + ceil(b_log2(sigma))));

  /* output = input; */
  /*  %% Non-codegen outputs */
  /*  % externally accessible fields */
  /*  output.firtaps = firtaps; */
  /*  output.nfirtaps = length(h); */
  /*  %output.filter = dfilter; */
  /*  output.gain = gain; */
  /*  %output.Hm1 = Hm1; */
  /*  %output.Hm2 = Hm2; */
  /*  %output.Hm3 = Hm3; */
  /*  %output.Hm4 = Hm4; */
  /*  %output.Hmd = Hmd; */
  /*  output.enables = enables; */
  /*  */
  /*  % internal fields used by the GUI */
  /*  %output.Hanalog = Hanalog; */
  /*  output.Apass_actual = Apass_actual; */
  /*  output.Astop_actual = Astop_actual; */
  /*  %output.delay = delay; */
  /*  %output.grpdelayvar = grpdelayvar; */
  /*  %output.Hd1 = Hd1; */
  /*  %output.Hd2 = Hd2; */
  /*  %output.Hmiddle = Hmiddle; */
  /*  output.a1 = a1; */
  /*  output.b1 = b1; */
  /*  output.a2 = a2; */
  /*  output.b2 = b2; */
  /*  For codegen only output taps */
  for (i0 = 0; i0 < 128; i0++) {
    sigma = rt_roundd_snf(firTapsPreScale[i0] * sigmax);
    if (sigma < 32768.0) {
      if (sigma >= -32768.0) {
        i4 = (short)sigma;
      } else {
        i4 = MIN_int16_T;
      }
    } else if (sigma >= 32768.0) {
      i4 = MAX_int16_T;
    } else {
      i4 = 0;
    }

    outputTaps[i0] = i4;
  }
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void internal_design_filter_cg_initialize(void)
{
  rt_InitInfAndNaN(8U);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void internal_design_filter_cg_terminate(void)
{
  /* (no terminate code required) */
}

/*
 * File trailer for internal_design_filter_cg.c
 *
 * [EOF]
 */
