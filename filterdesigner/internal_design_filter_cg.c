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
#include <math.h>
#include <string.h>
#include "rt_defines.h"
#include <float.h>
#include "internal_design_filter_cg.h"
#include "internal_design_filter_cg_emxutil.h"
#include <stdio.h>

/* Type Definitions */
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
static bool anyNonFinite(const emxArray_creal_T *x);
static void b_abs(const emxArray_creal_T *x, emxArray_real_T *y);
static void b_acos(double *x);
static void b_analogresp(const char type[2], const emxArray_real_T *f, double
  Fconverter, const double b1_data[], const int b1_size[2], const
  emxArray_creal_T *a1, const double b2_data[], const int b2_size[2], const
  emxArray_creal_T *a2, emxArray_creal_T *abc);
static void b_butter_cg(double Wn, double num[4], emxArray_creal_T *den);
static void b_cos(emxArray_real_T *x);
static void b_determineBestFractionLength(const double tap_store[128], double
  taps[128]);
static void b_exp(creal_T x[2048]);
static void b_firfreqz(const double b[15], const struct_T *options, creal_T h
  [2048], double w[2048]);
static void b_firpm_cg(double order, const double ff[4], const emxArray_real_T
  *amplitudes, const emxArray_real_T *frequencies, const emxArray_real_T
  *weights, emxArray_real_T *h, bool *valid, double *err);
static void b_fix(emxArray_real_T *x);
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
static void b_polyval(const double p[7], const creal_T x[2048], creal_T y[2048]);
static void b_power(const double a[2048], double y[2048]);
static void b_rdivide_helper(const emxArray_real_T *x, const emxArray_real_T *y,
  emxArray_real_T *z);
static void b_sqrt(double *x);
static bool b_strcmp(const char a[2]);
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
static void c_freqz_cg(const double b[7], const double w[2048], double Fs,
  creal_T hh[2048]);
static void c_generateCascadedResponseRx(const char enables[4], const
  emxArray_real_T *w, double Fs, const double hb1_coeff[15], const double
  hb2_coeff[7], const double hb3_coeff_data[], const int hb3_coeff_size[2],
  const double dec_int3_coeff_data[], const int dec_int3_coeff_size[2], const
  emxArray_real_T *extraTaps, emxArray_creal_T *combinedResponse);
static void c_polyval(const double p_data[], const int p_size[2], const creal_T
                      x[2048], creal_T y[2048]);
static void c_power(const emxArray_real_T *a, emxArray_real_T *y);
static void c_sqrt(creal_T *x);
static bool c_strcmp(const char a[2]);
static double c_sum(const double x[128]);
static void c_us(const double o_data[], const int o_size[2], double u_data[],
                 int u_size[2]);
static double dBinv(double dBinput);
static void d_abs(const double x[128], double y[128]);
static void d_analogresp(const emxArray_real_T *f, double Fconverter, const
  double b1_data[], const int b1_size[2], const emxArray_creal_T *a1, const
  double b2_data[], const int b2_size[2], const emxArray_creal_T *a2,
  emxArray_creal_T *abc);
static void d_firfreqz(double b_data[], int b_size[2], const struct_T *options,
  creal_T h[2048], double w[2048]);
static void d_freqz_cg(const double b_data[], const int b_size[2], const double
  w[2048], double Fs, creal_T hh[2048]);
static void d_polyval(const double p[29], const creal_T x[2048], creal_T y[2048]);
static void d_us(const double o[15], double u[29]);
static double db2mag(double ydb);
static void determineBestFractionLength(const emxArray_real_T *tap_store, double
  i, double M, emxArray_real_T *taps);
static int div_s32_floor(int numerator, int denominator);
static void e_firfreqz(const double b[29], const struct_T *options, creal_T h
  [2048], double w[2048]);
static void e_freqz_cg(const double b[29], const double w[2048], double Fs,
  creal_T hh[2048]);
static void e_polyval(const double p[13], const creal_T x[2048], creal_T y[2048]);
static void e_us(const double o[7], double u[13]);
static int eml_zlahqr(emxArray_creal_T *h);
static void f_firfreqz(const double b[13], const struct_T *options, creal_T h
  [2048], double w[2048]);
static void f_freqz_cg(const double b[13], const double w[2048], double Fs,
  creal_T hh[2048]);
static void f_polyval(const double p[57], const creal_T x[2048], creal_T y[2048]);
static void f_us(const double o[15], double u[57]);
static void firfreqz(const struct_T *options, creal_T h[2048], double w[2048]);
static void firpm_cg(double order, const double ff[4], const emxArray_real_T
                     *amplitudes, const emxArray_real_T *frequencies, const
                     emxArray_real_T *weights, emxArray_real_T *h);
static void firpmgrid_cg(double nfilt, const double ff[4], emxArray_real_T
  *gridactual);
static void freqs_cg(const double b_data[], const int b_size[2], const
                     emxArray_creal_T *a, const double w[2048], creal_T h[2048]);
static void freqz_cg(const double w[2048], double Fs, creal_T hh[2048]);
static void g_firfreqz(const double b[57], const struct_T *options, creal_T h
  [2048], double w[2048]);
static void g_freqz_cg(const double b[57], const double w[2048], double Fs,
  creal_T hh[2048]);
static void g_polyval(const double p[43], const creal_T x[2048], creal_T y[2048]);
static void g_us(const double o[15], double u[43]);
static void generateCascadedResponseRx(const char enables[4], const double w
  [2048], double Fs, const double hb1_coeff[15], const double hb2_coeff[7],
  const double hb3_coeff_data[], const int hb3_coeff_size[2], const double
  dec_int3_coeff_data[], const int dec_int3_coeff_size[2], creal_T
  combinedResponse[2048]);
static void h_firfreqz(const double b[43], const struct_T *options, creal_T h
  [2048], double w[2048]);
static void h_freqz_cg(const double b[43], const double w[2048], double Fs,
  creal_T hh[2048]);
static void h_polyval(const double p[19], const creal_T x[2048], creal_T y[2048]);
static void h_us(const double o[7], double u[19]);
static void i_firfreqz(const double b[19], const struct_T *options, creal_T h
  [2048], double w[2048]);
static void i_freqz_cg(const double b[19], const double w[2048], double Fs,
  creal_T hh[2048]);
static void i_polyval(const double p[85], const creal_T x[2048], creal_T y[2048]);
static void i_us(const double o[15], double u[85]);
static void interp1(const emxArray_real_T *varargin_1, const emxArray_real_T
                    *varargin_2, const emxArray_real_T *varargin_3,
                    emxArray_real_T *Vq);
static void j_firfreqz(const double b[85], const struct_T *options, creal_T h
  [2048], double w[2048]);
static void j_freqz_cg(const double b[85], const double w[2048], double Fs,
  creal_T hh[2048]);
static void j_polyval(const emxArray_real_T *p, const emxArray_creal_T *x,
                      emxArray_creal_T *y);
static void k_freqz_cg(const emxArray_real_T *w, double Fs, emxArray_creal_T *hh);
static void l_freqz_cg(const double b[15], const emxArray_real_T *w, double Fs,
  emxArray_creal_T *hh);
static void lp2lp_cg(const emxArray_creal_T *a, const emxArray_real_T *b, double
                     wo, emxArray_creal_T *at, emxArray_real_T *bt, double *dt);
static void m_freqz_cg(const double b[7], const emxArray_real_T *w, double Fs,
  emxArray_creal_T *hh);
static double mag2db(double y);
static double mpower(double a, double b);
static void n_freqz_cg(const emxArray_real_T *b, const emxArray_real_T *w,
  double Fs, emxArray_creal_T *hh);
static void o_freqz_cg(const double b[29], const emxArray_real_T *w, double Fs,
  emxArray_creal_T *hh);
static void p_freqz_cg(const double b[13], const emxArray_real_T *w, double Fs,
  emxArray_creal_T *hh);
static void poly(const emxArray_creal_T *x, emxArray_creal_T *c);
static void polyval(const double p[15], const creal_T x[2048], creal_T y[2048]);
static void power(const double a[2048], double y[2048]);
static void q_freqz_cg(const double b[57], const emxArray_real_T *w, double Fs,
  emxArray_creal_T *hh);
static void r_freqz_cg(const double b[43], const emxArray_real_T *w, double Fs,
  emxArray_creal_T *hh);
static void rdivide_helper(const emxArray_creal_T *x, const emxArray_creal_T *y,
  emxArray_creal_T *z);
static creal_T recip(const creal_T y);
static double remezdd(double k, double n, double m, const emxArray_real_T *x);
static void remezm(double nfilt, const double edge[4], const emxArray_real_T
                   *grid, emxArray_real_T *des, emxArray_real_T *wt,
                   emxArray_real_T *h, double *dev, bool *valid);
static void removeTrailingZero(const double b_data[], const int b_size[2], const
  emxArray_creal_T *a, double bR_data[], int bR_size[2], emxArray_creal_T *aR);
static double rt_atan2d_snf(double u0, double u1);
static double rt_hypotd_snf(double u0, double u1);
static double rt_powd_snf(double u0, double u1);
static double rt_remd_snf(double u0, double u1);
static double rt_roundd_snf(double u);
static void s_freqz_cg(const double b[19], const emxArray_real_T *w, double Fs,
  emxArray_creal_T *hh);
static void sinc(double x[2048]);
static double sum(const double x[2048]);
static void t_freqz_cg(const double b[85], const emxArray_real_T *w, double Fs,
  emxArray_creal_T *hh);
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
static void zp2ss_cg(emxArray_creal_T *a, emxArray_real_T *b, emxArray_real_T *c,
                     double *d);

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
  bool b_bool;
  int kstr;
  int exitg1;
  static const char cv33[2] = { 'T', 'x' };

  static const char cv34[2] = { 'R', 'x' };

  double dv14[2048];
  double dv15[2048];
  static creal_T dcv1[2048];
  double abc_im;
  double im;
  b_bool = false;
  kstr = 0;
  do {
    exitg1 = 0;
    if (kstr < 2) {
      if (type[kstr] != cv33[kstr]) {
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
    kstr = 0;
  } else {
    b_bool = false;
    kstr = 0;
    do {
      exitg1 = 0;
      if (kstr < 2) {
        if (type[kstr] != cv34[kstr]) {
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
      dv14[kstr] = f[kstr] / Fconverter;
    }

    sinc(dv14);
    for (kstr = 0; kstr < 2048; kstr++) {
      dv15[kstr] = 6.2831853071795862 * f[kstr];
    }

    freqs_cg(b1_data, b1_size, a1, dv15, abc);
    for (kstr = 0; kstr < 2048; kstr++) {
      dv15[kstr] = 6.2831853071795862 * f[kstr];
    }

    freqs_cg(b2_data, b2_size, a2, dv15, dcv1);
    for (kstr = 0; kstr < 2048; kstr++) {
      abc_im = dv14[kstr] * abc[kstr].re;
      im = dv14[kstr] * abc[kstr].im;
      abc[kstr].re = abc_im * dcv1[kstr].re - im * dcv1[kstr].im;
      abc[kstr].im = abc_im * dcv1[kstr].im + im * dcv1[kstr].re;
    }
    break;

   case 1:
    for (kstr = 0; kstr < 2048; kstr++) {
      dv14[kstr] = 6.2831853071795862 * f[kstr];
    }

    freqs_cg(b1_data, b1_size, a1, dv14, abc);
    for (kstr = 0; kstr < 2048; kstr++) {
      dv14[kstr] = 6.2831853071795862 * f[kstr];
    }

    freqs_cg(b2_data, b2_size, a2, dv14, dcv1);
    for (kstr = 0; kstr < 2048; kstr++) {
      dv14[kstr] = f[kstr] / Fconverter;
    }

    sinc(dv14);
    power(dv14, dv15);
    for (kstr = 0; kstr < 2048; kstr++) {
      abc_im = abc[kstr].re * dcv1[kstr].im + abc[kstr].im * dcv1[kstr].re;
      abc[kstr].re = dv15[kstr] * (abc[kstr].re * dcv1[kstr].re - abc[kstr].im *
        dcv1[kstr].im);
      abc[kstr].im = dv15[kstr] * abc_im;
    }
    break;

   default:
    /*  Default to Rx */
    for (kstr = 0; kstr < 2048; kstr++) {
      dv14[kstr] = 6.2831853071795862 * f[kstr];
    }

    freqs_cg(b1_data, b1_size, a1, dv14, abc);
    for (kstr = 0; kstr < 2048; kstr++) {
      dv14[kstr] = 6.2831853071795862 * f[kstr];
    }

    freqs_cg(b2_data, b2_size, a2, dv14, dcv1);
    for (kstr = 0; kstr < 2048; kstr++) {
      dv14[kstr] = f[kstr] / Fconverter;
    }

    sinc(dv14);
    power(dv14, dv15);
    for (kstr = 0; kstr < 2048; kstr++) {
      abc_im = abc[kstr].re * dcv1[kstr].im + abc[kstr].im * dcv1[kstr].re;
      abc[kstr].re = dv15[kstr] * (abc[kstr].re * dcv1[kstr].re - abc[kstr].im *
        dcv1[kstr].im);
      abc[kstr].im = dv15[kstr] * abc_im;
    }
    break;
  }
}

/*
 * Arguments    : const emxArray_creal_T *x
 * Return Type  : bool
 */
static bool anyNonFinite(const emxArray_creal_T *x)
{
  bool p;
  int nx;
  int k;
  nx = x->size[0] * x->size[1];
  p = true;
  for (k = 0; k < nx; k++) {
    if (p && ((!rtIsInf(x->data[k].re)) && (!rtIsInf(x->data[k].im)) &&
              ((!rtIsNaN(x->data[k].re)) && (!rtIsNaN(x->data[k].im))))) {
      p = true;
    } else {
      p = false;
    }
  }

  return !p;
}

/*
 * Arguments    : const emxArray_creal_T *x
 *                emxArray_real_T *y
 * Return Type  : void
 */
static void b_abs(const emxArray_creal_T *x, emxArray_real_T *y)
{
  int nx;
  int k;
  nx = x->size[1];
  k = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = x->size[1];
  emxEnsureCapacity_real_T(y, k);
  for (k = 0; k < nx; k++) {
    y->data[k] = rt_hypotd_snf(x->data[k].re, x->data[k].im);
  }
}

/*
 * Arguments    : double *x
 * Return Type  : void
 */
static void b_acos(double *x)
{
  *x = acos(*x);
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
  bool b_bool;
  int kstr;
  int exitg1;
  static const char cv46[2] = { 'T', 'x' };

  emxArray_creal_T *r6;
  emxArray_real_T *x_tmp;
  emxArray_real_T *r7;
  static const char cv47[2] = { 'R', 'x' };

  int i42;
  double abc_re;
  double abc_im;
  double re;
  double im;
  b_bool = false;
  kstr = 0;
  do {
    exitg1 = 0;
    if (kstr < 2) {
      if (type[kstr] != cv46[kstr]) {
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
    kstr = 0;
  } else {
    b_bool = false;
    kstr = 0;
    do {
      exitg1 = 0;
      if (kstr < 2) {
        if (type[kstr] != cv47[kstr]) {
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

  emxInit_creal_T(&r6, 2);
  emxInit_real_T(&x_tmp, 2);
  emxInit_real_T(&r7, 2);
  switch (kstr) {
   case 0:
    i42 = x_tmp->size[0] * x_tmp->size[1];
    x_tmp->size[0] = 1;
    x_tmp->size[1] = f->size[1];
    emxEnsureCapacity_real_T(x_tmp, i42);
    kstr = f->size[0] * f->size[1];
    for (i42 = 0; i42 < kstr; i42++) {
      x_tmp->data[i42] = f->data[i42] / Fconverter;
    }

    i42 = x_tmp->size[1];
    for (kstr = 0; kstr < i42; kstr++) {
      if (fabs(x_tmp->data[kstr]) < 1.0020841800044864E-292) {
        x_tmp->data[kstr] = 1.0;
      } else {
        x_tmp->data[kstr] *= 3.1415926535897931;
        x_tmp->data[kstr] = sin(x_tmp->data[kstr]) / x_tmp->data[kstr];
      }
    }

    i42 = r7->size[0] * r7->size[1];
    r7->size[0] = 1;
    r7->size[1] = f->size[1];
    emxEnsureCapacity_real_T(r7, i42);
    kstr = f->size[0] * f->size[1];
    for (i42 = 0; i42 < kstr; i42++) {
      r7->data[i42] = 6.2831853071795862 * f->data[i42];
    }

    b_freqs_cg(b1_data, b1_size, a1, r7, abc);
    i42 = r7->size[0] * r7->size[1];
    r7->size[0] = 1;
    r7->size[1] = f->size[1];
    emxEnsureCapacity_real_T(r7, i42);
    kstr = f->size[0] * f->size[1];
    for (i42 = 0; i42 < kstr; i42++) {
      r7->data[i42] = 6.2831853071795862 * f->data[i42];
    }

    b_freqs_cg(b2_data, b2_size, a2, r7, r6);
    i42 = x_tmp->size[0] * x_tmp->size[1];
    kstr = abc->size[0] * abc->size[1];
    abc->size[0] = 1;
    abc->size[1] = x_tmp->size[1];
    emxEnsureCapacity_creal_T(abc, kstr);
    kstr = i42 - 1;
    for (i42 = 0; i42 <= kstr; i42++) {
      abc_re = x_tmp->data[i42] * abc->data[i42].re;
      abc_im = x_tmp->data[i42] * abc->data[i42].im;
      re = r6->data[i42].re;
      im = r6->data[i42].im;
      abc->data[i42].re = abc_re * re - abc_im * im;
      abc->data[i42].im = abc_re * im + abc_im * re;
    }
    break;

   case 1:
    i42 = r7->size[0] * r7->size[1];
    r7->size[0] = 1;
    r7->size[1] = f->size[1];
    emxEnsureCapacity_real_T(r7, i42);
    kstr = f->size[0] * f->size[1];
    for (i42 = 0; i42 < kstr; i42++) {
      r7->data[i42] = 6.2831853071795862 * f->data[i42];
    }

    b_freqs_cg(b1_data, b1_size, a1, r7, abc);
    i42 = r7->size[0] * r7->size[1];
    r7->size[0] = 1;
    r7->size[1] = f->size[1];
    emxEnsureCapacity_real_T(r7, i42);
    kstr = f->size[0] * f->size[1];
    for (i42 = 0; i42 < kstr; i42++) {
      r7->data[i42] = 6.2831853071795862 * f->data[i42];
    }

    b_freqs_cg(b2_data, b2_size, a2, r7, r6);
    i42 = x_tmp->size[0] * x_tmp->size[1];
    x_tmp->size[0] = 1;
    x_tmp->size[1] = f->size[1];
    emxEnsureCapacity_real_T(x_tmp, i42);
    kstr = f->size[0] * f->size[1];
    for (i42 = 0; i42 < kstr; i42++) {
      x_tmp->data[i42] = f->data[i42] / Fconverter;
    }

    i42 = x_tmp->size[1];
    for (kstr = 0; kstr < i42; kstr++) {
      if (fabs(x_tmp->data[kstr]) < 1.0020841800044864E-292) {
        x_tmp->data[kstr] = 1.0;
      } else {
        x_tmp->data[kstr] *= 3.1415926535897931;
        x_tmp->data[kstr] = sin(x_tmp->data[kstr]) / x_tmp->data[kstr];
      }
    }

    c_power(x_tmp, r7);
    i42 = abc->size[0] * abc->size[1];
    kstr = abc->size[0] * abc->size[1];
    abc->size[0] = 1;
    emxEnsureCapacity_creal_T(abc, kstr);
    kstr = i42 - 1;
    for (i42 = 0; i42 <= kstr; i42++) {
      abc_re = abc->data[i42].re * r6->data[i42].re - abc->data[i42].im *
        r6->data[i42].im;
      abc_im = abc->data[i42].re * r6->data[i42].im + abc->data[i42].im *
        r6->data[i42].re;
      abc->data[i42].re = r7->data[i42] * abc_re;
      abc->data[i42].im = r7->data[i42] * abc_im;
    }
    break;

   default:
    /*  Default to Rx */
    i42 = r7->size[0] * r7->size[1];
    r7->size[0] = 1;
    r7->size[1] = f->size[1];
    emxEnsureCapacity_real_T(r7, i42);
    kstr = f->size[0] * f->size[1];
    for (i42 = 0; i42 < kstr; i42++) {
      r7->data[i42] = 6.2831853071795862 * f->data[i42];
    }

    b_freqs_cg(b1_data, b1_size, a1, r7, abc);
    i42 = r7->size[0] * r7->size[1];
    r7->size[0] = 1;
    r7->size[1] = f->size[1];
    emxEnsureCapacity_real_T(r7, i42);
    kstr = f->size[0] * f->size[1];
    for (i42 = 0; i42 < kstr; i42++) {
      r7->data[i42] = 6.2831853071795862 * f->data[i42];
    }

    b_freqs_cg(b2_data, b2_size, a2, r7, r6);
    i42 = x_tmp->size[0] * x_tmp->size[1];
    x_tmp->size[0] = 1;
    x_tmp->size[1] = f->size[1];
    emxEnsureCapacity_real_T(x_tmp, i42);
    kstr = f->size[0] * f->size[1];
    for (i42 = 0; i42 < kstr; i42++) {
      x_tmp->data[i42] = f->data[i42] / Fconverter;
    }

    i42 = x_tmp->size[1];
    for (kstr = 0; kstr < i42; kstr++) {
      if (fabs(x_tmp->data[kstr]) < 1.0020841800044864E-292) {
        x_tmp->data[kstr] = 1.0;
      } else {
        x_tmp->data[kstr] *= 3.1415926535897931;
        x_tmp->data[kstr] = sin(x_tmp->data[kstr]) / x_tmp->data[kstr];
      }
    }

    c_power(x_tmp, r7);
    i42 = abc->size[0] * abc->size[1];
    kstr = abc->size[0] * abc->size[1];
    abc->size[0] = 1;
    emxEnsureCapacity_creal_T(abc, kstr);
    kstr = i42 - 1;
    for (i42 = 0; i42 <= kstr; i42++) {
      abc_re = abc->data[i42].re * r6->data[i42].re - abc->data[i42].im *
        r6->data[i42].im;
      abc_im = abc->data[i42].re * r6->data[i42].im + abc->data[i42].im *
        r6->data[i42].re;
      abc->data[i42].re = r7->data[i42] * abc_re;
      abc->data[i42].im = r7->data[i42] * abc_im;
    }
    break;
  }

  emxFree_real_T(&r7);
  emxFree_real_T(&x_tmp);
  emxFree_creal_T(&r6);
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
  emxArray_real_T *c;
  emxArray_creal_T *b_a;
  emxArray_real_T *b_b;
  double d;
  double y_re;
  double y_im;
  int i5;
  int k;
  emxInit_creal_T(&a, 2);
  emxInit_real_T(&b, 1);
  emxInit_real_T(&c, 2);
  emxInit_creal_T(&b_a, 2);
  emxInit_real_T(&b_b, 1);

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
  lp2lp_cg(a, b, Wn, b_a, b_b, &d);

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
  i5 = den->size[1];
  emxFree_real_T(&b_b);
  emxFree_creal_T(&b_a);
  emxFree_real_T(&c);
  emxFree_real_T(&b);
  emxFree_creal_T(&a);
  for (k = 0; k <= i5 - 2; k++) {
    d = 0.0 * y_im + 0.0 * y_re;
    y_re = (0.0 * y_re - 0.0 * y_im) + den->data[k + 1].re;
    y_im = d + den->data[k + 1].im;
  }

  d = 0.0 * y_re;
  if (0.0 * y_im == 0.0) {
    d /= 0.037037037037037035;
  } else if (d == 0.0) {
    d = 0.0;
  } else {
    d = rtNaN;
  }

  num[0] = d;
  d = 0.0 * y_re;
  if (0.0 * y_im == 0.0) {
    d /= 0.037037037037037035;
  } else if (d == 0.0) {
    d = 0.0;
  } else {
    d = rtNaN;
  }

  num[1] = d;
  d = 0.0 * y_re;
  if (0.0 * y_im == 0.0) {
    d /= 0.037037037037037035;
  } else if (d == 0.0) {
    d = 0.0;
  } else {
    d = rtNaN;
  }

  num[2] = d;
  d = 0.037037037037037035 * y_re;
  if (0.037037037037037035 * y_im == 0.0) {
    d /= 0.037037037037037035;
  } else if (d == 0.0) {
    d = 0.0;
  } else {
    d /= 0.037037037037037035;
  }

  num[3] = d;

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
  for (k = 0; k < nx; k++) {
    x->data[k] = cos(x->data[k]);
  }
}

/*
 * Codegen workaround for fixed fi call requirements
 * Arguments    : const double tap_store[128]
 *                double taps[128]
 * Return Type  : void
 */
static void b_determineBestFractionLength(const double tap_store[128], double
  taps[128])
{
  double r[2048];
  int i64;
  double b_r[128];
  double dv16[128];
  double u;
  double e[16];
  double v;
  double dv17[128];
  int r_tmp;
  short i65;
  double dv18[128];
  double dv19[128];
  double dv20[128];
  double dv21[128];
  double dv22[128];
  double dv23[128];
  double dv24[128];
  double dv25[128];
  double dv26[128];
  double dv27[128];
  double dv28[128];
  double dv29[128];
  double dv30[128];
  double dv31[128];
  int k;
  bool exitg1;
  memset(&r[0], 0, sizeof(double) << 11);
  for (i64 = 0; i64 < 128; i64++) {
    u = tap_store[i64] * 2.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i65 = (short)u;
      } else {
        i65 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i65 = MAX_int16_T;
    } else {
      i65 = 0;
    }

    u = (double)i65 * 0.5;
    r_tmp = i64 << 4;
    r[r_tmp] = u;
    b_r[i64] = u - tap_store[i64];
    u = tap_store[i64] * 4.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i65 = (short)u;
      } else {
        i65 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i65 = MAX_int16_T;
    } else {
      i65 = 0;
    }

    r[1 + r_tmp] = (double)i65 * 0.25;
  }

  d_abs(b_r, dv16);
  e[0] = c_sum(dv16);
  for (i64 = 0; i64 < 128; i64++) {
    r_tmp = i64 << 4;
    b_r[i64] = r[1 + r_tmp] - tap_store[i64];
    u = tap_store[i64] * 8.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i65 = (short)u;
      } else {
        i65 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i65 = MAX_int16_T;
    } else {
      i65 = 0;
    }

    r[2 + r_tmp] = (double)i65 * 0.125;
  }

  d_abs(b_r, dv17);
  e[1] = c_sum(dv17);
  for (i64 = 0; i64 < 128; i64++) {
    r_tmp = i64 << 4;
    b_r[i64] = r[2 + r_tmp] - tap_store[i64];
    u = tap_store[i64] * 16.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i65 = (short)u;
      } else {
        i65 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i65 = MAX_int16_T;
    } else {
      i65 = 0;
    }

    r[3 + r_tmp] = (double)i65 * 0.0625;
  }

  d_abs(b_r, dv18);
  e[2] = c_sum(dv18);
  for (i64 = 0; i64 < 128; i64++) {
    r_tmp = i64 << 4;
    b_r[i64] = r[3 + r_tmp] - tap_store[i64];
    u = tap_store[i64] * 32.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i65 = (short)u;
      } else {
        i65 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i65 = MAX_int16_T;
    } else {
      i65 = 0;
    }

    r[4 + r_tmp] = (double)i65 * 0.03125;
  }

  d_abs(b_r, dv19);
  e[3] = c_sum(dv19);
  for (i64 = 0; i64 < 128; i64++) {
    r_tmp = i64 << 4;
    b_r[i64] = r[4 + r_tmp] - tap_store[i64];
    u = tap_store[i64] * 64.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i65 = (short)u;
      } else {
        i65 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i65 = MAX_int16_T;
    } else {
      i65 = 0;
    }

    r[5 + r_tmp] = (double)i65 * 0.015625;
  }

  d_abs(b_r, dv20);
  e[4] = c_sum(dv20);
  for (i64 = 0; i64 < 128; i64++) {
    r_tmp = i64 << 4;
    b_r[i64] = r[5 + r_tmp] - tap_store[i64];
    u = tap_store[i64] * 128.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i65 = (short)u;
      } else {
        i65 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i65 = MAX_int16_T;
    } else {
      i65 = 0;
    }

    r[6 + r_tmp] = (double)i65 * 0.0078125;
  }

  d_abs(b_r, dv21);
  e[5] = c_sum(dv21);
  for (i64 = 0; i64 < 128; i64++) {
    r_tmp = i64 << 4;
    b_r[i64] = r[6 + r_tmp] - tap_store[i64];
    u = tap_store[i64] * 256.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i65 = (short)u;
      } else {
        i65 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i65 = MAX_int16_T;
    } else {
      i65 = 0;
    }

    r[7 + r_tmp] = (double)i65 * 0.00390625;
  }

  d_abs(b_r, dv22);
  e[6] = c_sum(dv22);
  for (i64 = 0; i64 < 128; i64++) {
    r_tmp = i64 << 4;
    b_r[i64] = r[7 + r_tmp] - tap_store[i64];
    u = tap_store[i64] * 512.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i65 = (short)u;
      } else {
        i65 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i65 = MAX_int16_T;
    } else {
      i65 = 0;
    }

    r[8 + r_tmp] = (double)i65 * 0.001953125;
  }

  d_abs(b_r, dv23);
  e[7] = c_sum(dv23);
  for (i64 = 0; i64 < 128; i64++) {
    r_tmp = i64 << 4;
    b_r[i64] = r[8 + r_tmp] - tap_store[i64];
    u = tap_store[i64] * 1024.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i65 = (short)u;
      } else {
        i65 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i65 = MAX_int16_T;
    } else {
      i65 = 0;
    }

    r[9 + r_tmp] = (double)i65 * 0.0009765625;
  }

  d_abs(b_r, dv24);
  e[8] = c_sum(dv24);
  for (i64 = 0; i64 < 128; i64++) {
    r_tmp = i64 << 4;
    b_r[i64] = r[9 + r_tmp] - tap_store[i64];
    u = tap_store[i64] * 2048.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i65 = (short)u;
      } else {
        i65 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i65 = MAX_int16_T;
    } else {
      i65 = 0;
    }

    r[10 + r_tmp] = (double)i65 * 0.00048828125;
  }

  d_abs(b_r, dv25);
  e[9] = c_sum(dv25);
  for (i64 = 0; i64 < 128; i64++) {
    r_tmp = i64 << 4;
    b_r[i64] = r[10 + r_tmp] - tap_store[i64];
    u = tap_store[i64] * 4096.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i65 = (short)u;
      } else {
        i65 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i65 = MAX_int16_T;
    } else {
      i65 = 0;
    }

    r[11 + r_tmp] = (double)i65 * 0.000244140625;
  }

  d_abs(b_r, dv26);
  e[10] = c_sum(dv26);
  for (i64 = 0; i64 < 128; i64++) {
    r_tmp = i64 << 4;
    b_r[i64] = r[11 + r_tmp] - tap_store[i64];
    u = tap_store[i64] * 8192.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i65 = (short)u;
      } else {
        i65 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i65 = MAX_int16_T;
    } else {
      i65 = 0;
    }

    r[12 + r_tmp] = (double)i65 * 0.0001220703125;
  }

  d_abs(b_r, dv27);
  e[11] = c_sum(dv27);
  for (i64 = 0; i64 < 128; i64++) {
    r_tmp = i64 << 4;
    b_r[i64] = r[12 + r_tmp] - tap_store[i64];
    u = tap_store[i64] * 16384.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i65 = (short)u;
      } else {
        i65 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i65 = MAX_int16_T;
    } else {
      i65 = 0;
    }

    r[13 + r_tmp] = (double)i65 * 6.103515625E-5;
  }

  d_abs(b_r, dv28);
  e[12] = c_sum(dv28);
  for (i64 = 0; i64 < 128; i64++) {
    r_tmp = i64 << 4;
    b_r[i64] = r[13 + r_tmp] - tap_store[i64];
    u = tap_store[i64] * 32768.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i65 = (short)u;
      } else {
        i65 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i65 = MAX_int16_T;
    } else {
      i65 = 0;
    }

    r[14 + r_tmp] = (double)i65 * 3.0517578125E-5;
  }

  d_abs(b_r, dv29);
  e[13] = c_sum(dv29);
  for (i64 = 0; i64 < 128; i64++) {
    r_tmp = i64 << 4;
    b_r[i64] = r[14 + r_tmp] - tap_store[i64];
    u = tap_store[i64] * 65536.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i65 = (short)u;
      } else {
        i65 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i65 = MAX_int16_T;
    } else {
      i65 = 0;
    }

    r[15 + r_tmp] = (double)i65 * 1.52587890625E-5;
  }

  d_abs(b_r, dv30);
  e[14] = c_sum(dv30);
  for (i64 = 0; i64 < 128; i64++) {
    b_r[i64] = r[15 + (i64 << 4)] - tap_store[i64];
  }

  d_abs(b_r, dv31);
  e[15] = c_sum(dv31);
  if (!rtIsNaN(e[0])) {
    r_tmp = 1;
  } else {
    r_tmp = 0;
    k = 2;
    exitg1 = false;
    while ((!exitg1) && (k < 17)) {
      if (!rtIsNaN(e[k - 1])) {
        r_tmp = k;
        exitg1 = true;
      } else {
        k++;
      }
    }
  }

  if (r_tmp == 0) {
    r_tmp = 1;
  } else {
    v = e[r_tmp - 1];
    i64 = r_tmp + 1;
    for (k = i64; k < 17; k++) {
      u = e[k - 1];
      if (v > u) {
        v = u;
        r_tmp = k;
      }
    }
  }

  for (i64 = 0; i64 < 128; i64++) {
    taps[i64] = r[(r_tmp + (i64 << 4)) - 1];
  }
}

/*
 * Arguments    : creal_T x[2048]
 * Return Type  : void
 */
static void b_exp(creal_T x[2048])
{
  int k;
  double r;
  double x_im;
  for (k = 0; k < 2048; k++) {
    if (x[k].im == 0.0) {
      x[k].re = exp(x[k].re);
      x[k].im = 0.0;
    } else if (rtIsInf(x[k].im) && rtIsInf(x[k].re) && (x[k].re < 0.0)) {
      x[k].re = 0.0;
      x[k].im = 0.0;
    } else {
      r = exp(x[k].re / 2.0);
      x_im = x[k].im;
      x[k].re = r * (r * cos(x[k].im));
      x[k].im = r * (r * sin(x_im));
    }
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
  int i73;
  creal_T dcv2[2048];
  double re_tmp;
  double re;
  creal_T s_tmp[2048];
  double h_re;
  double brm;

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
  for (i73 = 0; i73 < 2048; i73++) {
    w[i73] = options->w[i73];
    re_tmp = 6.2831853071795862 * options->w[i73] / options->Fs;
    re = re_tmp * 0.0;
    s_tmp[i73].re = re;
    s_tmp[i73].im = re_tmp;
    dcv2[i73].re = re;
    dcv2[i73].im = re_tmp;
  }

  b_exp(dcv2);
  polyval(b, dcv2, h);
  for (i73 = 0; i73 < 2048; i73++) {
    s_tmp[i73].re *= 14.0;
    s_tmp[i73].im *= 14.0;
  }

  b_exp(s_tmp);
  for (i73 = 0; i73 < 2048; i73++) {
    h_re = h[i73].re;
    if (s_tmp[i73].im == 0.0) {
      if (h[i73].im == 0.0) {
        h[i73].re /= s_tmp[i73].re;
        h[i73].im = 0.0;
      } else if (h[i73].re == 0.0) {
        h[i73].re = 0.0;
        h[i73].im /= s_tmp[i73].re;
      } else {
        h[i73].re /= s_tmp[i73].re;
        h[i73].im /= s_tmp[i73].re;
      }
    } else if (s_tmp[i73].re == 0.0) {
      if (h[i73].re == 0.0) {
        h[i73].re = h[i73].im / s_tmp[i73].im;
        h[i73].im = 0.0;
      } else if (h[i73].im == 0.0) {
        h[i73].re = 0.0;
        h[i73].im = -(h_re / s_tmp[i73].im);
      } else {
        h[i73].re = h[i73].im / s_tmp[i73].im;
        h[i73].im = -(h_re / s_tmp[i73].im);
      }
    } else {
      brm = fabs(s_tmp[i73].re);
      re = fabs(s_tmp[i73].im);
      if (brm > re) {
        re = s_tmp[i73].im / s_tmp[i73].re;
        re_tmp = s_tmp[i73].re + re * s_tmp[i73].im;
        h[i73].re = (h[i73].re + re * h[i73].im) / re_tmp;
        h[i73].im = (h[i73].im - re * h_re) / re_tmp;
      } else if (re == brm) {
        if (s_tmp[i73].re > 0.0) {
          re = 0.5;
        } else {
          re = -0.5;
        }

        if (s_tmp[i73].im > 0.0) {
          re_tmp = 0.5;
        } else {
          re_tmp = -0.5;
        }

        h[i73].re = (h[i73].re * re + h[i73].im * re_tmp) / brm;
        h[i73].im = (h[i73].im * re - h_re * re_tmp) / brm;
      } else {
        re = s_tmp[i73].re / s_tmp[i73].im;
        re_tmp = s_tmp[i73].im + re * s_tmp[i73].re;
        h[i73].re = (re * h[i73].re + h[i73].im) / re_tmp;
        h[i73].im = (re * h[i73].im - h_re) / re_tmp;
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
 *                bool *valid
 *                double *err
 * Return Type  : void
 */
static void b_firpm_cg(double order, const double ff[4], const emxArray_real_T
  *amplitudes, const emxArray_real_T *frequencies, const emxArray_real_T
  *weights, emxArray_real_T *h, bool *valid, double *err)
{
  emxArray_real_T *grid;
  emxArray_real_T *r13;
  emxArray_real_T *r14;
  double b_ff[4];
  emxArray_real_T *b_grid;
  int i57;
  int loop_ub;
  emxArray_real_T *b_h;
  double x;
  int h_idx_0;
  double d2;
  int i58;
  int i59;
  emxArray_real_T *c_h;
  emxInit_real_T(&grid, 2);
  emxInit_real_T(&r13, 2);
  emxInit_real_T(&r14, 2);

  /*  */
  firpmgrid_cg(order + 1.0, ff, grid);

  /*      orgFreqIndx = frequencies; */
  /*   */
  /*      positionsOfNewFreqIndx = zeros(size(grid)); */
  /*      for ind = 1:length(positionsOfNewFreqIndx) */
  /*   */
  /*          [~,indx] = min( abs(orgFreqIndx-grid(ind)) ); */
  /*   */
  /*          positionsOfNewFreqIndx(ind) = indx; */
  /*      end */
  /*   */
  /*      wt = weights(positionsOfNewFreqIndx); */
  /*      des = amplitudes(positionsOfNewFreqIndx); */
  /*  Workaround */
  /* ftype = 2; */
  /* sign_val = 1; */
  /*  Always bandpass designs */
  /*  cast to enforce precision rules */
  /*  Call actual design algorithm */
  interp1(frequencies, amplitudes, grid, r13);
  interp1(frequencies, weights, grid, r14);
  b_ff[0] = ff[0] / 2.0;
  b_ff[1] = ff[1] / 2.0;
  b_ff[2] = ff[2] / 2.0;
  b_ff[3] = ff[3] / 2.0;
  emxInit_real_T(&b_grid, 2);
  i57 = b_grid->size[0] * b_grid->size[1];
  b_grid->size[0] = 1;
  b_grid->size[1] = grid->size[1];
  emxEnsureCapacity_real_T(b_grid, i57);
  loop_ub = grid->size[0] * grid->size[1];
  for (i57 = 0; i57 < loop_ub; i57++) {
    b_grid->data[i57] = grid->data[i57] / 2.0;
  }

  emxFree_real_T(&grid);
  emxInit_real_T(&b_h, 2);
  remezm(order + 1.0, b_ff, b_grid, r13, r14, b_h, &x, valid);
  h_idx_0 = b_h->size[0] * b_h->size[1];
  i57 = h->size[0] * h->size[1];
  h->size[0] = 1;
  h->size[1] = h_idx_0;
  emxEnsureCapacity_real_T(h, i57);
  emxFree_real_T(&b_grid);
  emxFree_real_T(&r14);
  emxFree_real_T(&r13);
  for (i57 = 0; i57 < h_idx_0; i57++) {
    h->data[i57] = b_h->data[i57];
  }

  emxFree_real_T(&b_h);

  /*  make it a row */
  d2 = (double)h->size[1] - rt_remd_snf(order + 1.0, 2.0);
  if (1.0 > d2) {
    i57 = 1;
    i58 = 1;
    i59 = 0;
  } else {
    i57 = (int)d2;
    i58 = -1;
    i59 = 1;
  }

  emxInit_real_T(&c_h, 2);
  h_idx_0 = c_h->size[0] * c_h->size[1];
  c_h->size[0] = 1;
  loop_ub = div_s32_floor(i59 - i57, i58);
  c_h->size[1] = (h->size[1] + loop_ub) + 1;
  emxEnsureCapacity_real_T(c_h, h_idx_0);
  h_idx_0 = h->size[1];
  for (i59 = 0; i59 < h_idx_0; i59++) {
    c_h->data[i59] = h->data[i59];
  }

  for (i59 = 0; i59 <= loop_ub; i59++) {
    c_h->data[i59 + h->size[1]] = h->data[(i57 + i58 * i59) - 1];
  }

  i57 = h->size[0] * h->size[1];
  h->size[0] = 1;
  h->size[1] = c_h->size[1];
  emxEnsureCapacity_real_T(h, i57);
  loop_ub = c_h->size[0] * c_h->size[1];
  for (i57 = 0; i57 < loop_ub; i57++) {
    h->data[i57] = c_h->data[i57];
  }

  i57 = h->size[1];
  i58 = c_h->size[0] * c_h->size[1];
  c_h->size[0] = 1;
  loop_ub = div_s32_floor(1 - i57, -1);
  c_h->size[1] = loop_ub + 1;
  emxEnsureCapacity_real_T(c_h, i58);
  for (i58 = 0; i58 <= loop_ub; i58++) {
    c_h->data[i58] = h->data[(i57 - i58) - 1];
  }

  i57 = h->size[0] * h->size[1];
  h->size[0] = 1;
  h->size[1] = c_h->size[1];
  emxEnsureCapacity_real_T(h, i57);
  loop_ub = c_h->size[0] * c_h->size[1];
  for (i57 = 0; i57 < loop_ub; i57++) {
    h->data[i57] = c_h->data[i57];
  }

  emxFree_real_T(&c_h);
  *err = fabs(x);
}

/*
 * Arguments    : emxArray_real_T *x
 * Return Type  : void
 */
static void b_fix(emxArray_real_T *x)
{
  int nx;
  int k;
  nx = x->size[1];
  for (k = 0; k < nx; k++) {
    x->data[k] = trunc(x->data[k]);
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
  int i43;
  int loop_ub;
  emxArray_creal_T *y;
  emxArray_real_T c_b_data;
  int k;
  int i44;
  double a_re;
  double a_im;
  double s_re;
  double s_im;
  emxInit_creal_T(&s, 2);
  emxInit_creal_T(&b_a, 2);
  removeTrailingZero(b_data, b_size, a, b_b_data, b_b_size, b_a);
  i43 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = w->size[1];
  emxEnsureCapacity_creal_T(s, i43);
  loop_ub = w->size[0] * w->size[1];
  for (i43 = 0; i43 < loop_ub; i43++) {
    s->data[i43].re = w->data[i43] * 0.0;
    s->data[i43].im = w->data[i43];
  }

  emxInit_creal_T(&y, 2);
  i43 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = s->size[1];
  emxEnsureCapacity_creal_T(y, i43);
  if ((y->size[1] == 0) || (b_a->size[1] == 0)) {
  } else {
    i43 = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity_creal_T(y, i43);
    loop_ub = y->size[1];
    for (i43 = 0; i43 < loop_ub; i43++) {
      y->data[i43] = b_a->data[0];
    }

    i43 = b_a->size[1];
    for (k = 0; k <= i43 - 2; k++) {
      i44 = s->size[0] * s->size[1];
      loop_ub = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = s->size[1];
      emxEnsureCapacity_creal_T(y, loop_ub);
      a_re = b_a->data[k + 1].re;
      a_im = b_a->data[k + 1].im;
      loop_ub = i44 - 1;
      for (i44 = 0; i44 <= loop_ub; i44++) {
        s_re = s->data[i44].re * y->data[i44].re - s->data[i44].im * y->data[i44]
          .im;
        s_im = s->data[i44].re * y->data[i44].im + s->data[i44].im * y->data[i44]
          .re;
        y->data[i44].re = s_re + a_re;
        y->data[i44].im = s_im + a_im;
      }
    }
  }

  c_b_data.data = &b_b_data[0];
  c_b_data.size = &b_b_size[0];
  c_b_data.allocatedSize = 4;
  c_b_data.numDimensions = 2;
  c_b_data.canFreeData = false;
  j_polyval(&c_b_data, s, b_a);
  rdivide_helper(b_a, y, h);
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
  int i10;
  static const char cv14[8] = { 'o', 'n', 'e', 's', 'i', 'd', 'e', 'd' };

  double b_w[2048];
  static const char cv15[7] = { 'o', 'm', 'i', 't', 't', 'e', 'd' };

  /*  Cast to enforce precision rules */
  options.Fs = Fs;
  memcpy(&options.nfft[0], &w[0], sizeof(double) << 11);
  memcpy(&options.w[0], &w[0], sizeof(double) << 11);

  /*  Remaining are default or for advanced use */
  options.fvflag = 1.0;
  for (i10 = 0; i10 < 8; i10++) {
    options.range[i10] = cv14[i10];
  }

  options.centerdc = 0.0;
  for (i10 = 0; i10 < 7; i10++) {
    options.configlevel[i10] = cv15[i10];
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
  bool b_bool;
  int ix;
  int exitg1;
  emxArray_creal_T *d2;
  emxArray_creal_T *d3;
  static const char cv35[4] = { '2', '1', '1', '1' };

  double u[15];
  int iy;
  double tmp_data[29];
  int tmp_size[2];
  double b_u[30];
  double c_u[14];
  double d_u[60];
  double e_u[45];
  double f_u[21];
  double g_u[90];
  emxArray_real_T b_tmp_data;
  double h_u[7];
  int k;
  static const char cv36[4] = { '1', '2', '1', '1' };

  static const char cv37[4] = { '1', '1', '2', '1' };

  double combinedResponse_re;
  double combinedResponse_im;
  double d2_re;
  double d2_im;
  static const char cv38[4] = { '2', '2', '1', '1' };

  static const char cv39[4] = { '2', '1', '2', '1' };

  static const char cv40[4] = { '1', '2', '2', '1' };

  static const char cv41[4] = { '2', '2', '2', '1' };

  static const char cv42[4] = { '1', '1', '1', '3' };

  static const char cv43[4] = { '2', '1', '1', '3' };

  static const char cv44[4] = { '1', '2', '1', '3' };

  static const char cv45[4] = { '2', '2', '1', '3' };

  /*  Cast */
  b_bool = false;
  ix = 0;
  do {
    exitg1 = 0;
    if (ix < 4) {
      if (enables[ix] != '1') {
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
    ix = 0;
  } else {
    b_bool = false;
    ix = 0;
    do {
      exitg1 = 0;
      if (ix < 4) {
        if (enables[ix] != cv35[ix]) {
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
      ix = 1;
    } else {
      b_bool = false;
      ix = 0;
      do {
        exitg1 = 0;
        if (ix < 4) {
          if (enables[ix] != cv36[ix]) {
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
        ix = 2;
      } else {
        b_bool = false;
        ix = 0;
        do {
          exitg1 = 0;
          if (ix < 4) {
            if (enables[ix] != cv37[ix]) {
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
          ix = 3;
        } else {
          b_bool = false;
          ix = 0;
          do {
            exitg1 = 0;
            if (ix < 4) {
              if (enables[ix] != cv38[ix]) {
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
            ix = 4;
          } else {
            b_bool = false;
            ix = 0;
            do {
              exitg1 = 0;
              if (ix < 4) {
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
              ix = 5;
            } else {
              b_bool = false;
              ix = 0;
              do {
                exitg1 = 0;
                if (ix < 4) {
                  if (enables[ix] != cv40[ix]) {
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
                ix = 6;
              } else {
                b_bool = false;
                ix = 0;
                do {
                  exitg1 = 0;
                  if (ix < 4) {
                    if (enables[ix] != cv41[ix]) {
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
                  ix = 7;
                } else {
                  b_bool = false;
                  ix = 0;
                  do {
                    exitg1 = 0;
                    if (ix < 4) {
                      if (enables[ix] != cv42[ix]) {
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
                    ix = 8;
                  } else {
                    b_bool = false;
                    ix = 0;
                    do {
                      exitg1 = 0;
                      if (ix < 4) {
                        if (enables[ix] != cv43[ix]) {
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
                      ix = 9;
                    } else {
                      b_bool = false;
                      ix = 0;
                      do {
                        exitg1 = 0;
                        if (ix < 4) {
                          if (enables[ix] != cv44[ix]) {
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
                        ix = 10;
                      } else {
                        b_bool = false;
                        ix = 0;
                        do {
                          exitg1 = 0;
                          if (ix < 4) {
                            if (enables[ix] != cv45[ix]) {
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
    k_freqz_cg(w, Fs, combinedResponse);
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

    l_freqz_cg(u, w, Fs, combinedResponse);
    break;

   case 2:
    /*  Hb2 */
    for (iy = 0; iy < 7; iy++) {
      h_u[iy] = 0.0;
    }

    ix = 0;
    iy = 0;
    for (k = 0; k < 7; k++) {
      h_u[iy] = hb2_coeff[ix];
      ix++;
      iy++;
    }

    m_freqz_cg(h_u, w, Fs, combinedResponse);
    break;

   case 3:
    /*  Hb3 */
    c_us(hb3_coeff_data, hb3_coeff_size, tmp_data, tmp_size);
    b_tmp_data.data = &tmp_data[0];
    b_tmp_data.size = &tmp_size[0];
    b_tmp_data.allocatedSize = 29;
    b_tmp_data.numDimensions = 2;
    b_tmp_data.canFreeData = false;
    n_freqz_cg(&b_tmp_data, w, Fs, combinedResponse);
    break;

   case 4:
    /*  Hb2,Hb1 */
    memset(&b_u[0], 0, 30U * sizeof(double));
    ix = 0;
    iy = 0;
    for (k = 0; k < 15; k++) {
      b_u[iy] = hb1_coeff[ix];
      ix++;
      iy += 2;
    }

    for (iy = 0; iy < 7; iy++) {
      h_u[iy] = 0.0;
    }

    ix = 0;
    iy = 0;
    for (k = 0; k < 7; k++) {
      h_u[iy] = hb2_coeff[ix];
      ix++;
      iy++;
    }

    o_freqz_cg(*(double (*)[29])&b_u[0], w, Fs, combinedResponse);
    m_freqz_cg(h_u, w, Fs, d2);
    iy = combinedResponse->size[0] * combinedResponse->size[1];
    ix = combinedResponse->size[0] * combinedResponse->size[1];
    combinedResponse->size[0] = 1;
    emxEnsureCapacity_creal_T(combinedResponse, ix);
    ix = iy - 1;
    for (iy = 0; iy <= ix; iy++) {
      combinedResponse_re = combinedResponse->data[iy].re;
      combinedResponse_im = combinedResponse->data[iy].im;
      d2_re = d2->data[iy].re;
      d2_im = d2->data[iy].im;
      combinedResponse->data[iy].re = combinedResponse_re * d2_re -
        combinedResponse_im * d2_im;
      combinedResponse->data[iy].im = combinedResponse_re * d2_im +
        combinedResponse_im * d2_re;
    }
    break;

   case 5:
    /*  Hb3,Hb1 */
    memset(&b_u[0], 0, 30U * sizeof(double));
    ix = 0;
    iy = 0;
    for (k = 0; k < 15; k++) {
      b_u[iy] = hb1_coeff[ix];
      ix++;
      iy += 2;
    }

    c_us(hb3_coeff_data, hb3_coeff_size, tmp_data, tmp_size);
    b_tmp_data.data = &tmp_data[0];
    b_tmp_data.size = &tmp_size[0];
    b_tmp_data.allocatedSize = 29;
    b_tmp_data.numDimensions = 2;
    b_tmp_data.canFreeData = false;
    n_freqz_cg(&b_tmp_data, w, Fs, d2);
    o_freqz_cg(*(double (*)[29])&b_u[0], w, Fs, combinedResponse);
    iy = combinedResponse->size[0] * combinedResponse->size[1];
    ix = combinedResponse->size[0] * combinedResponse->size[1];
    combinedResponse->size[0] = 1;
    emxEnsureCapacity_creal_T(combinedResponse, ix);
    ix = iy - 1;
    for (iy = 0; iy <= ix; iy++) {
      combinedResponse_re = combinedResponse->data[iy].re;
      combinedResponse_im = combinedResponse->data[iy].im;
      d2_re = d2->data[iy].re;
      d2_im = d2->data[iy].im;
      combinedResponse->data[iy].re = combinedResponse_re * d2_re -
        combinedResponse_im * d2_im;
      combinedResponse->data[iy].im = combinedResponse_re * d2_im +
        combinedResponse_im * d2_re;
    }
    break;

   case 6:
    /*  Hb3,Hb2 */
    memset(&c_u[0], 0, 14U * sizeof(double));
    ix = 0;
    iy = 0;
    for (k = 0; k < 7; k++) {
      c_u[iy] = hb2_coeff[ix];
      ix++;
      iy += 2;
    }

    c_us(hb3_coeff_data, hb3_coeff_size, tmp_data, tmp_size);
    b_tmp_data.data = &tmp_data[0];
    b_tmp_data.size = &tmp_size[0];
    b_tmp_data.allocatedSize = 29;
    b_tmp_data.numDimensions = 2;
    b_tmp_data.canFreeData = false;
    n_freqz_cg(&b_tmp_data, w, Fs, d2);
    p_freqz_cg(*(double (*)[13])&c_u[0], w, Fs, combinedResponse);
    iy = combinedResponse->size[0] * combinedResponse->size[1];
    ix = combinedResponse->size[0] * combinedResponse->size[1];
    combinedResponse->size[0] = 1;
    emxEnsureCapacity_creal_T(combinedResponse, ix);
    ix = iy - 1;
    for (iy = 0; iy <= ix; iy++) {
      combinedResponse_re = combinedResponse->data[iy].re;
      combinedResponse_im = combinedResponse->data[iy].im;
      d2_re = d2->data[iy].re;
      d2_im = d2->data[iy].im;
      combinedResponse->data[iy].re = combinedResponse_re * d2_re -
        combinedResponse_im * d2_im;
      combinedResponse->data[iy].im = combinedResponse_re * d2_im +
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

    memset(&c_u[0], 0, 14U * sizeof(double));
    ix = 0;
    iy = 0;
    for (k = 0; k < 7; k++) {
      c_u[iy] = hb2_coeff[ix];
      ix++;
      iy += 2;
    }

    c_us(hb3_coeff_data, hb3_coeff_size, tmp_data, tmp_size);
    b_tmp_data.data = &tmp_data[0];
    b_tmp_data.size = &tmp_size[0];
    b_tmp_data.allocatedSize = 29;
    b_tmp_data.numDimensions = 2;
    b_tmp_data.canFreeData = false;
    n_freqz_cg(&b_tmp_data, w, Fs, d3);
    q_freqz_cg(*(double (*)[57])&d_u[0], w, Fs, combinedResponse);
    p_freqz_cg(*(double (*)[13])&c_u[0], w, Fs, d2);
    iy = combinedResponse->size[0] * combinedResponse->size[1];
    ix = combinedResponse->size[0] * combinedResponse->size[1];
    combinedResponse->size[0] = 1;
    emxEnsureCapacity_creal_T(combinedResponse, ix);
    ix = iy - 1;
    for (iy = 0; iy <= ix; iy++) {
      combinedResponse_re = combinedResponse->data[iy].re * d2->data[iy].re -
        combinedResponse->data[iy].im * d2->data[iy].im;
      combinedResponse_im = combinedResponse->data[iy].re * d2->data[iy].im +
        combinedResponse->data[iy].im * d2->data[iy].re;
      d2_re = d3->data[iy].re;
      d2_im = d3->data[iy].im;
      combinedResponse->data[iy].re = combinedResponse_re * d2_re -
        combinedResponse_im * d2_im;
      combinedResponse->data[iy].im = combinedResponse_re * d2_im +
        combinedResponse_im * d2_re;
    }
    break;

   case 8:
    /*  Dec/Int3 */
    c_us(dec_int3_coeff_data, dec_int3_coeff_size, tmp_data, tmp_size);
    b_tmp_data.data = &tmp_data[0];
    b_tmp_data.size = &tmp_size[0];
    b_tmp_data.allocatedSize = 29;
    b_tmp_data.numDimensions = 2;
    b_tmp_data.canFreeData = false;
    n_freqz_cg(&b_tmp_data, w, Fs, combinedResponse);

    /*  RECHECK ALL DEC BY 3     */
    break;

   case 9:
    /*  Dec/Int3,Hb1 */
    memset(&e_u[0], 0, 45U * sizeof(double));
    ix = 0;
    iy = 0;
    for (k = 0; k < 15; k++) {
      e_u[iy] = hb1_coeff[ix];
      ix++;
      iy += 3;
    }

    c_us(dec_int3_coeff_data, dec_int3_coeff_size, tmp_data, tmp_size);
    b_tmp_data.data = &tmp_data[0];
    b_tmp_data.size = &tmp_size[0];
    b_tmp_data.allocatedSize = 29;
    b_tmp_data.numDimensions = 2;
    b_tmp_data.canFreeData = false;
    n_freqz_cg(&b_tmp_data, w, Fs, d2);
    r_freqz_cg(*(double (*)[43])&e_u[0], w, Fs, combinedResponse);
    iy = combinedResponse->size[0] * combinedResponse->size[1];
    ix = combinedResponse->size[0] * combinedResponse->size[1];
    combinedResponse->size[0] = 1;
    emxEnsureCapacity_creal_T(combinedResponse, ix);
    ix = iy - 1;
    for (iy = 0; iy <= ix; iy++) {
      combinedResponse_re = combinedResponse->data[iy].re;
      combinedResponse_im = combinedResponse->data[iy].im;
      d2_re = d2->data[iy].re;
      d2_im = d2->data[iy].im;
      combinedResponse->data[iy].re = combinedResponse_re * d2_re -
        combinedResponse_im * d2_im;
      combinedResponse->data[iy].im = combinedResponse_re * d2_im +
        combinedResponse_im * d2_re;
    }
    break;

   case 10:
    /*  Dec/Int3,Hb2 */
    memset(&f_u[0], 0, 21U * sizeof(double));
    ix = 0;
    iy = 0;
    for (k = 0; k < 7; k++) {
      f_u[iy] = hb2_coeff[ix];
      ix++;
      iy += 3;
    }

    c_us(dec_int3_coeff_data, dec_int3_coeff_size, tmp_data, tmp_size);
    b_tmp_data.data = &tmp_data[0];
    b_tmp_data.size = &tmp_size[0];
    b_tmp_data.allocatedSize = 29;
    b_tmp_data.numDimensions = 2;
    b_tmp_data.canFreeData = false;
    n_freqz_cg(&b_tmp_data, w, Fs, d2);
    s_freqz_cg(*(double (*)[19])&f_u[0], w, Fs, combinedResponse);
    iy = combinedResponse->size[0] * combinedResponse->size[1];
    ix = combinedResponse->size[0] * combinedResponse->size[1];
    combinedResponse->size[0] = 1;
    emxEnsureCapacity_creal_T(combinedResponse, ix);
    ix = iy - 1;
    for (iy = 0; iy <= ix; iy++) {
      combinedResponse_re = combinedResponse->data[iy].re;
      combinedResponse_im = combinedResponse->data[iy].im;
      d2_re = d2->data[iy].re;
      d2_im = d2->data[iy].im;
      combinedResponse->data[iy].re = combinedResponse_re * d2_re -
        combinedResponse_im * d2_im;
      combinedResponse->data[iy].im = combinedResponse_re * d2_im +
        combinedResponse_im * d2_re;
    }
    break;

   case 11:
    /*  Dec/Int3,Hb2,Hb1 {Hm4,Hm2c34,Hm1} */
    memset(&g_u[0], 0, 90U * sizeof(double));
    ix = 0;
    iy = 0;
    for (k = 0; k < 15; k++) {
      g_u[iy] = hb1_coeff[ix];
      ix++;
      iy += 6;
    }

    memset(&f_u[0], 0, 21U * sizeof(double));
    ix = 0;
    iy = 0;
    for (k = 0; k < 7; k++) {
      f_u[iy] = hb2_coeff[ix];
      ix++;
      iy += 3;
    }

    c_us(dec_int3_coeff_data, dec_int3_coeff_size, tmp_data, tmp_size);
    b_tmp_data.data = &tmp_data[0];
    b_tmp_data.size = &tmp_size[0];
    b_tmp_data.allocatedSize = 29;
    b_tmp_data.numDimensions = 2;
    b_tmp_data.canFreeData = false;
    n_freqz_cg(&b_tmp_data, w, Fs, d3);
    t_freqz_cg(*(double (*)[85])&g_u[0], w, Fs, combinedResponse);
    s_freqz_cg(*(double (*)[19])&f_u[0], w, Fs, d2);
    iy = combinedResponse->size[0] * combinedResponse->size[1];
    ix = combinedResponse->size[0] * combinedResponse->size[1];
    combinedResponse->size[0] = 1;
    emxEnsureCapacity_creal_T(combinedResponse, ix);
    ix = iy - 1;
    for (iy = 0; iy <= ix; iy++) {
      combinedResponse_re = combinedResponse->data[iy].re * d2->data[iy].re -
        combinedResponse->data[iy].im * d2->data[iy].im;
      combinedResponse_im = combinedResponse->data[iy].re * d2->data[iy].im +
        combinedResponse->data[iy].im * d2->data[iy].re;
      d2_re = d3->data[iy].re;
      d2_im = d3->data[iy].im;
      combinedResponse->data[iy].re = combinedResponse_re * d2_re -
        combinedResponse_im * d2_im;
      combinedResponse->data[iy].im = combinedResponse_re * d2_im +
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
  bool b10;
  bool b11;
  double t;
  int eint;
  if (x == 0.0) {
    f = rtMinusInf;
  } else if (x < 0.0) {
    f = rtNaN;
  } else {
    b10 = !rtIsInf(x);
    b11 = !rtIsNaN(x);
    if (b10 && b11) {
      t = frexp(x, &eint);
      if (t == 0.5) {
        f = (double)eint - 1.0;
      } else if ((eint == 1) && (t < 0.75)) {
        f = log(2.0 * t) / 0.69314718055994529;
      } else {
        f = log(t) / 0.69314718055994529 + (double)eint;
      }
    } else {
      f = x;
    }
  }

  return f;
}

/*
 * Arguments    : const double p[7]
 *                const creal_T x[2048]
 *                creal_T y[2048]
 * Return Type  : void
 */
static void b_polyval(const double p[7], const creal_T x[2048], creal_T y[2048])
{
  int i11;
  int k;
  double p_re;
  double x_im;
  for (i11 = 0; i11 < 2048; i11++) {
    y[i11].re = p[0];
    y[i11].im = 0.0;
  }

  for (k = 0; k < 6; k++) {
    p_re = p[k + 1];
    for (i11 = 0; i11 < 2048; i11++) {
      x_im = x[i11].re * y[i11].im + x[i11].im * y[i11].re;
      y[i11].re = (x[i11].re * y[i11].re - x[i11].im * y[i11].im) + p_re;
      y[i11].im = x_im;
    }
  }
}

/*
 * Arguments    : const double a[2048]
 *                double y[2048]
 * Return Type  : void
 */
static void b_power(const double a[2048], double y[2048])
{
  int k;
  for (k = 0; k < 2048; k++) {
    y[k] = a[k] * a[k];
  }
}

/*
 * Arguments    : const emxArray_real_T *x
 *                const emxArray_real_T *y
 *                emxArray_real_T *z
 * Return Type  : void
 */
static void b_rdivide_helper(const emxArray_real_T *x, const emxArray_real_T *y,
  emxArray_real_T *z)
{
  int i54;
  int loop_ub;
  i54 = z->size[0] * z->size[1];
  z->size[0] = 1;
  z->size[1] = x->size[1];
  emxEnsureCapacity_real_T(z, i54);
  loop_ub = x->size[0] * x->size[1];
  for (i54 = 0; i54 < loop_ub; i54++) {
    z->data[i54] = x->data[i54] / y->data[i54];
  }
}

/*
 * Arguments    : double *x
 * Return Type  : void
 */
static void b_sqrt(double *x)
{
  *x = sqrt(*x);
}

/*
 * Arguments    : const char a[2]
 * Return Type  : bool
 */
static bool b_strcmp(const char a[2])
{
  bool b_bool;
  int kstr;
  int exitg1;
  static const char cv0[2] = { 'R', 'x' };

  b_bool = false;
  kstr = 0;
  do {
    exitg1 = 0;
    if (kstr < 2) {
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
  int vlen;
  int k;
  vlen = x->size[1];
  if (x->size[1] == 0) {
    y = 0.0;
  } else {
    y = x->data[0];
    for (k = 2; k <= vlen; k++) {
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
  int i71;
  int k;
  double x_re;
  double x_im;
  if (incx >= 1) {
    i71 = ix0 + incx * (n - 1);
    for (k = ix0; incx < 0 ? k >= i71 : k <= i71; k += incx) {
      x_re = x->data[k - 1].re;
      x_im = x->data[k - 1].im;
      x->data[k - 1].re = a.re * x_re - a.im * x_im;
      x->data[k - 1].im = a.re * x_im + a.im * x_re;
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
  double y_tmp;
  double scale;
  double b_y_tmp;
  double f2s;
  double f2;
  double fs_re;
  double fs_im;
  double gs_re;
  double gs_im;
  bool guard1 = false;
  double g2s;
  y_tmp = fabs(f.re);
  scale = y_tmp;
  b_y_tmp = fabs(f.im);
  if (b_y_tmp > y_tmp) {
    scale = b_y_tmp;
  }

  f2s = fabs(g.re);
  f2 = fabs(g.im);
  if (f2 > f2s) {
    f2s = f2;
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
    f2 = fs_re * fs_re + fs_im * fs_im;
    scale = gs_re * gs_re + gs_im * gs_im;
    f2s = scale;
    if (1.0 > scale) {
      f2s = 1.0;
    }

    if (f2 <= f2s * 2.0041683600089728E-292) {
      if ((f.re == 0.0) && (f.im == 0.0)) {
        *cs = 0.0;
        scale = rt_hypotd_snf(gs_re, gs_im);
        sn->re = gs_re / scale;
        sn->im = -gs_im / scale;
      } else {
        g2s = sqrt(scale);
        *cs = rt_hypotd_snf(fs_re, fs_im) / g2s;
        if (b_y_tmp > y_tmp) {
          y_tmp = b_y_tmp;
        }

        if (y_tmp > 1.0) {
          scale = rt_hypotd_snf(f.re, f.im);
          fs_re = f.re / scale;
          fs_im = f.im / scale;
        } else {
          f2 = 7.4428285367870146E+137 * f.re;
          f2s = 7.4428285367870146E+137 * f.im;
          scale = rt_hypotd_snf(f2, f2s);
          fs_re = f2 / scale;
          fs_im = f2s / scale;
        }

        gs_re /= g2s;
        gs_im = -gs_im / g2s;
        sn->re = fs_re * gs_re - fs_im * gs_im;
        sn->im = fs_re * gs_im + fs_im * gs_re;
      }
    } else {
      f2s = sqrt(1.0 + scale / f2);
      fs_re *= f2s;
      fs_im *= f2s;
      *cs = 1.0 / f2s;
      scale += f2;
      fs_re /= scale;
      fs_im /= scale;
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
  emxArray_real_T *c;
  int i4;
  int tmp_size[2];
  creal_T tmp_data[1];
  int b_tmp_size[1];
  double b_tmp_data[1];
  emxArray_creal_T *a;
  emxArray_real_T *b;
  emxArray_creal_T *r1;
  emxArray_creal_T c_tmp_data;
  emxArray_real_T d_tmp_data;
  double d;
  int loop_ub;
  double y_im;
  double ar;
  emxInit_real_T(&c, 2);

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
  i4 = c->size[0] * c->size[1];
  c->size[0] = 1;
  c->size[1] = 1;
  emxEnsureCapacity_real_T(c, i4);
  c->data[0] = 1.0;
  tmp_size[0] = 1;
  tmp_size[1] = 1;
  tmp_data[0].re = -1.0;
  tmp_data[0].im = 0.0;
  b_tmp_size[0] = 1;
  b_tmp_data[0] = 1.0;
  emxInit_creal_T(&a, 2);
  emxInit_real_T(&b, 1);
  emxInit_creal_T(&r1, 2);
  c_tmp_data.data = &tmp_data[0];
  c_tmp_data.size = &tmp_size[0];
  c_tmp_data.allocatedSize = 1;
  c_tmp_data.numDimensions = 2;
  c_tmp_data.canFreeData = false;
  d_tmp_data.data = &b_tmp_data[0];
  d_tmp_data.size = &b_tmp_size[0];
  d_tmp_data.allocatedSize = 1;
  d_tmp_data.numDimensions = 1;
  d_tmp_data.canFreeData = false;
  lp2lp_cg(&c_tmp_data, &d_tmp_data, Wn, a, b, &d);

  /*  step 5: Use Bilinear transformation to find discrete equivalent: */
  /*  nargout <= 3 */
  /*  Transform to zero-pole-gain and polynomial forms: */
  /*  nargout <= 2 */
  poly(a, r1);
  den_size[0] = 1;
  den_size[1] = r1->size[1];
  loop_ub = r1->size[0] * r1->size[1];
  emxFree_real_T(&c);
  emxFree_real_T(&b);
  emxFree_creal_T(&a);
  for (i4 = 0; i4 < loop_ub; i4++) {
    den_data[i4] = r1->data[i4];
  }

  emxFree_creal_T(&r1);

  /*  This internal function returns more exact numerator vectors */
  /*  for the num/den case. */
  /*  Wn input is two element band edge vector */
  /* --------------------------------- */
  /*  lowpass */
  d = (0.0 * den_data[0].re - 0.0 * den_data[0].im) + den_data[1].re;
  y_im = (0.0 * den_data[0].im + 0.0 * den_data[0].re) + den_data[1].im;
  ar = 0.0 * d;
  if (!(0.0 * y_im == 0.0)) {
    if (ar == 0.0) {
      ar = 0.0;
    } else {
      ar = rtNaN;
    }
  }

  num[0] = ar;
  if (y_im == 0.0) {
    ar = d;
  } else if (d == 0.0) {
    ar = 0.0;
  } else {
    ar = d;
  }

  num[1] = ar;

  /*  num = poly(a-b*c)+(d-1)*den; */
}

/*
 * Arguments    : const emxArray_real_T *x
 *                emxArray_real_T *y
 * Return Type  : void
 */
static void c_abs(const emxArray_real_T *x, emxArray_real_T *y)
{
  int nx;
  int k;
  nx = x->size[1];
  k = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = x->size[1];
  emxEnsureCapacity_real_T(y, k);
  for (k = 0; k < nx; k++) {
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
  emxArray_real_T *r16;
  int i62;
  int loop_ub;
  emxArray_creal_T *r17;
  emxArray_real_T *x;
  double abc_re;
  double abc_im;
  emxInit_real_T(&r16, 2);
  i62 = r16->size[0] * r16->size[1];
  r16->size[0] = 1;
  r16->size[1] = f->size[1];
  emxEnsureCapacity_real_T(r16, i62);
  loop_ub = f->size[0] * f->size[1];
  for (i62 = 0; i62 < loop_ub; i62++) {
    r16->data[i62] = 6.2831853071795862 * f->data[i62];
  }

  b_freqs_cg(b1_data, b1_size, a1, r16, abc);
  i62 = r16->size[0] * r16->size[1];
  r16->size[0] = 1;
  r16->size[1] = f->size[1];
  emxEnsureCapacity_real_T(r16, i62);
  loop_ub = f->size[0] * f->size[1];
  for (i62 = 0; i62 < loop_ub; i62++) {
    r16->data[i62] = 6.2831853071795862 * f->data[i62];
  }

  emxInit_creal_T(&r17, 2);
  emxInit_real_T(&x, 2);
  b_freqs_cg(b2_data, b2_size, a2, r16, r17);
  i62 = x->size[0] * x->size[1];
  x->size[0] = 1;
  x->size[1] = f->size[1];
  emxEnsureCapacity_real_T(x, i62);
  loop_ub = f->size[0] * f->size[1];
  for (i62 = 0; i62 < loop_ub; i62++) {
    x->data[i62] = f->data[i62] / Fconverter;
  }

  i62 = x->size[1];
  for (loop_ub = 0; loop_ub < i62; loop_ub++) {
    if (fabs(x->data[loop_ub]) < 1.0020841800044864E-292) {
      x->data[loop_ub] = 1.0;
    } else {
      x->data[loop_ub] *= 3.1415926535897931;
      x->data[loop_ub] = sin(x->data[loop_ub]) / x->data[loop_ub];
    }
  }

  c_power(x, r16);
  i62 = abc->size[0] * abc->size[1];
  loop_ub = abc->size[0] * abc->size[1];
  abc->size[0] = 1;
  emxEnsureCapacity_creal_T(abc, loop_ub);
  loop_ub = i62 - 1;
  emxFree_real_T(&x);
  for (i62 = 0; i62 <= loop_ub; i62++) {
    abc_re = abc->data[i62].re * r17->data[i62].re - abc->data[i62].im *
      r17->data[i62].im;
    abc_im = abc->data[i62].re * r17->data[i62].im + abc->data[i62].im *
      r17->data[i62].re;
    abc->data[i62].re = r16->data[i62] * abc_re;
    abc->data[i62].im = r16->data[i62] * abc_im;
  }

  emxFree_real_T(&r16);
  emxFree_creal_T(&r17);
}

/*
 * Arguments    : emxArray_creal_T *x
 * Return Type  : void
 */
static void c_exp(emxArray_creal_T *x)
{
  int nx;
  int k;
  double r;
  double x_im;
  double b_x_im;
  nx = x->size[1];
  for (k = 0; k < nx; k++) {
    if (x->data[k].im == 0.0) {
      r = x->data[k].re;
      x->data[k].re = exp(r);
      x->data[k].im = 0.0;
    } else if (rtIsInf(x->data[k].im) && rtIsInf(x->data[k].re) && (x->data[k].
                re < 0.0)) {
      x->data[k].re = 0.0;
      x->data[k].im = 0.0;
    } else {
      r = exp(x->data[k].re / 2.0);
      x_im = x->data[k].im;
      b_x_im = x->data[k].im;
      x->data[k].re = r * (r * cos(x_im));
      x->data[k].im = r * (r * sin(b_x_im));
    }
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
  int i74;
  creal_T dcv3[2048];
  double re_tmp;
  double re;
  creal_T s_tmp[2048];
  double h_re;
  double brm;

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
  for (i74 = 0; i74 < 2048; i74++) {
    w[i74] = options->w[i74];
    re_tmp = 6.2831853071795862 * options->w[i74] / options->Fs;
    re = re_tmp * 0.0;
    s_tmp[i74].re = re;
    s_tmp[i74].im = re_tmp;
    dcv3[i74].re = re;
    dcv3[i74].im = re_tmp;
  }

  b_exp(dcv3);
  b_polyval(b, dcv3, h);
  for (i74 = 0; i74 < 2048; i74++) {
    s_tmp[i74].re *= 6.0;
    s_tmp[i74].im *= 6.0;
  }

  b_exp(s_tmp);
  for (i74 = 0; i74 < 2048; i74++) {
    h_re = h[i74].re;
    if (s_tmp[i74].im == 0.0) {
      if (h[i74].im == 0.0) {
        h[i74].re /= s_tmp[i74].re;
        h[i74].im = 0.0;
      } else if (h[i74].re == 0.0) {
        h[i74].re = 0.0;
        h[i74].im /= s_tmp[i74].re;
      } else {
        h[i74].re /= s_tmp[i74].re;
        h[i74].im /= s_tmp[i74].re;
      }
    } else if (s_tmp[i74].re == 0.0) {
      if (h[i74].re == 0.0) {
        h[i74].re = h[i74].im / s_tmp[i74].im;
        h[i74].im = 0.0;
      } else if (h[i74].im == 0.0) {
        h[i74].re = 0.0;
        h[i74].im = -(h_re / s_tmp[i74].im);
      } else {
        h[i74].re = h[i74].im / s_tmp[i74].im;
        h[i74].im = -(h_re / s_tmp[i74].im);
      }
    } else {
      brm = fabs(s_tmp[i74].re);
      re = fabs(s_tmp[i74].im);
      if (brm > re) {
        re = s_tmp[i74].im / s_tmp[i74].re;
        re_tmp = s_tmp[i74].re + re * s_tmp[i74].im;
        h[i74].re = (h[i74].re + re * h[i74].im) / re_tmp;
        h[i74].im = (h[i74].im - re * h_re) / re_tmp;
      } else if (re == brm) {
        if (s_tmp[i74].re > 0.0) {
          re = 0.5;
        } else {
          re = -0.5;
        }

        if (s_tmp[i74].im > 0.0) {
          re_tmp = 0.5;
        } else {
          re_tmp = -0.5;
        }

        h[i74].re = (h[i74].re * re + h[i74].im * re_tmp) / brm;
        h[i74].im = (h[i74].im * re - h_re * re_tmp) / brm;
      } else {
        re = s_tmp[i74].re / s_tmp[i74].im;
        re_tmp = s_tmp[i74].im + re * s_tmp[i74].re;
        h[i74].re = (re * h[i74].re + h[i74].im) / re_tmp;
        h[i74].im = (re * h[i74].im - h_re) / re_tmp;
      }
    }
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
  int i12;
  static const char cv16[8] = { 'o', 'n', 'e', 's', 'i', 'd', 'e', 'd' };

  double b_w[2048];
  static const char cv17[7] = { 'o', 'm', 'i', 't', 't', 'e', 'd' };

  /*  Cast to enforce precision rules */
  options.Fs = Fs;
  memcpy(&options.nfft[0], &w[0], sizeof(double) << 11);
  memcpy(&options.w[0], &w[0], sizeof(double) << 11);

  /*  Remaining are default or for advanced use */
  options.fvflag = 1.0;
  for (i12 = 0; i12 < 8; i12++) {
    options.range[i12] = cv16[i12];
  }

  options.centerdc = 0.0;
  for (i12 = 0; i12 < 7; i12++) {
    options.configlevel[i12] = cv17[i12];
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
 *                const emxArray_real_T *extraTaps
 *                emxArray_creal_T *combinedResponse
 * Return Type  : void
 */
static void c_generateCascadedResponseRx(const char enables[4], const
  emxArray_real_T *w, double Fs, const double hb1_coeff[15], const double
  hb2_coeff[7], const double hb3_coeff_data[], const int hb3_coeff_size[2],
  const double dec_int3_coeff_data[], const int dec_int3_coeff_size[2], const
  emxArray_real_T *extraTaps, emxArray_creal_T *combinedResponse)
{
  bool b_bool;
  int ix;
  int exitg1;
  emxArray_creal_T *d2;
  emxArray_creal_T *d3;
  static const char cv48[4] = { '2', '1', '1', '1' };

  double u[15];
  int iy;
  double tmp_data[29];
  int tmp_size[2];
  double b_u[30];
  double c_u[14];
  double d_u[60];
  double e_u[45];
  double f_u[21];
  double g_u[90];
  int pre_FIR_us;
  emxArray_real_T b_tmp_data;
  double h_u[7];
  emxArray_real_T *i_u;
  int k;
  int vlenx;
  static const char cv49[4] = { '1', '2', '1', '1' };

  static const char cv50[4] = { '1', '1', '2', '1' };

  double combinedResponse_re;
  double combinedResponse_im;
  double d2_re;
  double d2_im;
  static const char cv51[4] = { '2', '2', '1', '1' };

  static const char cv52[4] = { '2', '1', '2', '1' };

  static const char cv53[4] = { '1', '2', '2', '1' };

  static const char cv54[4] = { '2', '2', '2', '1' };

  static const char cv55[4] = { '1', '1', '1', '3' };

  static const char cv56[4] = { '2', '1', '1', '3' };

  static const char cv57[4] = { '1', '2', '1', '3' };

  static const char cv58[4] = { '2', '2', '1', '3' };

  /*  Cast */
  b_bool = false;
  ix = 0;
  do {
    exitg1 = 0;
    if (ix < 4) {
      if (enables[ix] != '1') {
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
    ix = 0;
  } else {
    b_bool = false;
    ix = 0;
    do {
      exitg1 = 0;
      if (ix < 4) {
        if (enables[ix] != cv48[ix]) {
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
      ix = 1;
    } else {
      b_bool = false;
      ix = 0;
      do {
        exitg1 = 0;
        if (ix < 4) {
          if (enables[ix] != cv49[ix]) {
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
        ix = 2;
      } else {
        b_bool = false;
        ix = 0;
        do {
          exitg1 = 0;
          if (ix < 4) {
            if (enables[ix] != cv50[ix]) {
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
          ix = 3;
        } else {
          b_bool = false;
          ix = 0;
          do {
            exitg1 = 0;
            if (ix < 4) {
              if (enables[ix] != cv51[ix]) {
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
            ix = 4;
          } else {
            b_bool = false;
            ix = 0;
            do {
              exitg1 = 0;
              if (ix < 4) {
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
              ix = 5;
            } else {
              b_bool = false;
              ix = 0;
              do {
                exitg1 = 0;
                if (ix < 4) {
                  if (enables[ix] != cv53[ix]) {
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
                ix = 6;
              } else {
                b_bool = false;
                ix = 0;
                do {
                  exitg1 = 0;
                  if (ix < 4) {
                    if (enables[ix] != cv54[ix]) {
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
                  ix = 7;
                } else {
                  b_bool = false;
                  ix = 0;
                  do {
                    exitg1 = 0;
                    if (ix < 4) {
                      if (enables[ix] != cv55[ix]) {
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
                    ix = 8;
                  } else {
                    b_bool = false;
                    ix = 0;
                    do {
                      exitg1 = 0;
                      if (ix < 4) {
                        if (enables[ix] != cv56[ix]) {
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
                      ix = 9;
                    } else {
                      b_bool = false;
                      ix = 0;
                      do {
                        exitg1 = 0;
                        if (ix < 4) {
                          if (enables[ix] != cv57[ix]) {
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
                        ix = 10;
                      } else {
                        b_bool = false;
                        ix = 0;
                        do {
                          exitg1 = 0;
                          if (ix < 4) {
                            if (enables[ix] != cv58[ix]) {
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
    k_freqz_cg(w, Fs, combinedResponse);
    pre_FIR_us = 1;
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

    l_freqz_cg(u, w, Fs, combinedResponse);
    pre_FIR_us = 2;
    break;

   case 2:
    /*  Hb2 */
    for (iy = 0; iy < 7; iy++) {
      h_u[iy] = 0.0;
    }

    ix = 0;
    iy = 0;
    for (k = 0; k < 7; k++) {
      h_u[iy] = hb2_coeff[ix];
      ix++;
      iy++;
    }

    m_freqz_cg(h_u, w, Fs, combinedResponse);
    pre_FIR_us = 2;
    break;

   case 3:
    /*  Hb3 */
    c_us(hb3_coeff_data, hb3_coeff_size, tmp_data, tmp_size);
    b_tmp_data.data = &tmp_data[0];
    b_tmp_data.size = &tmp_size[0];
    b_tmp_data.allocatedSize = 29;
    b_tmp_data.numDimensions = 2;
    b_tmp_data.canFreeData = false;
    n_freqz_cg(&b_tmp_data, w, Fs, combinedResponse);
    pre_FIR_us = 2;
    break;

   case 4:
    /*  Hb2,Hb1 */
    memset(&b_u[0], 0, 30U * sizeof(double));
    ix = 0;
    iy = 0;
    for (k = 0; k < 15; k++) {
      b_u[iy] = hb1_coeff[ix];
      ix++;
      iy += 2;
    }

    for (iy = 0; iy < 7; iy++) {
      h_u[iy] = 0.0;
    }

    ix = 0;
    iy = 0;
    for (k = 0; k < 7; k++) {
      h_u[iy] = hb2_coeff[ix];
      ix++;
      iy++;
    }

    o_freqz_cg(*(double (*)[29])&b_u[0], w, Fs, combinedResponse);
    m_freqz_cg(h_u, w, Fs, d2);
    iy = combinedResponse->size[0] * combinedResponse->size[1];
    ix = combinedResponse->size[0] * combinedResponse->size[1];
    combinedResponse->size[0] = 1;
    emxEnsureCapacity_creal_T(combinedResponse, ix);
    ix = iy - 1;
    for (iy = 0; iy <= ix; iy++) {
      combinedResponse_re = combinedResponse->data[iy].re;
      combinedResponse_im = combinedResponse->data[iy].im;
      d2_re = d2->data[iy].re;
      d2_im = d2->data[iy].im;
      combinedResponse->data[iy].re = combinedResponse_re * d2_re -
        combinedResponse_im * d2_im;
      combinedResponse->data[iy].im = combinedResponse_re * d2_im +
        combinedResponse_im * d2_re;
    }

    pre_FIR_us = 4;
    break;

   case 5:
    /*  Hb3,Hb1 */
    memset(&b_u[0], 0, 30U * sizeof(double));
    ix = 0;
    iy = 0;
    for (k = 0; k < 15; k++) {
      b_u[iy] = hb1_coeff[ix];
      ix++;
      iy += 2;
    }

    c_us(hb3_coeff_data, hb3_coeff_size, tmp_data, tmp_size);
    b_tmp_data.data = &tmp_data[0];
    b_tmp_data.size = &tmp_size[0];
    b_tmp_data.allocatedSize = 29;
    b_tmp_data.numDimensions = 2;
    b_tmp_data.canFreeData = false;
    n_freqz_cg(&b_tmp_data, w, Fs, d2);
    o_freqz_cg(*(double (*)[29])&b_u[0], w, Fs, combinedResponse);
    iy = combinedResponse->size[0] * combinedResponse->size[1];
    ix = combinedResponse->size[0] * combinedResponse->size[1];
    combinedResponse->size[0] = 1;
    emxEnsureCapacity_creal_T(combinedResponse, ix);
    ix = iy - 1;
    for (iy = 0; iy <= ix; iy++) {
      combinedResponse_re = combinedResponse->data[iy].re;
      combinedResponse_im = combinedResponse->data[iy].im;
      d2_re = d2->data[iy].re;
      d2_im = d2->data[iy].im;
      combinedResponse->data[iy].re = combinedResponse_re * d2_re -
        combinedResponse_im * d2_im;
      combinedResponse->data[iy].im = combinedResponse_re * d2_im +
        combinedResponse_im * d2_re;
    }

    pre_FIR_us = 4;
    break;

   case 6:
    /*  Hb3,Hb2 */
    memset(&c_u[0], 0, 14U * sizeof(double));
    ix = 0;
    iy = 0;
    for (k = 0; k < 7; k++) {
      c_u[iy] = hb2_coeff[ix];
      ix++;
      iy += 2;
    }

    c_us(hb3_coeff_data, hb3_coeff_size, tmp_data, tmp_size);
    b_tmp_data.data = &tmp_data[0];
    b_tmp_data.size = &tmp_size[0];
    b_tmp_data.allocatedSize = 29;
    b_tmp_data.numDimensions = 2;
    b_tmp_data.canFreeData = false;
    n_freqz_cg(&b_tmp_data, w, Fs, d2);
    p_freqz_cg(*(double (*)[13])&c_u[0], w, Fs, combinedResponse);
    iy = combinedResponse->size[0] * combinedResponse->size[1];
    ix = combinedResponse->size[0] * combinedResponse->size[1];
    combinedResponse->size[0] = 1;
    emxEnsureCapacity_creal_T(combinedResponse, ix);
    ix = iy - 1;
    for (iy = 0; iy <= ix; iy++) {
      combinedResponse_re = combinedResponse->data[iy].re;
      combinedResponse_im = combinedResponse->data[iy].im;
      d2_re = d2->data[iy].re;
      d2_im = d2->data[iy].im;
      combinedResponse->data[iy].re = combinedResponse_re * d2_re -
        combinedResponse_im * d2_im;
      combinedResponse->data[iy].im = combinedResponse_re * d2_im +
        combinedResponse_im * d2_re;
    }

    pre_FIR_us = 4;
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

    memset(&c_u[0], 0, 14U * sizeof(double));
    ix = 0;
    iy = 0;
    for (k = 0; k < 7; k++) {
      c_u[iy] = hb2_coeff[ix];
      ix++;
      iy += 2;
    }

    c_us(hb3_coeff_data, hb3_coeff_size, tmp_data, tmp_size);
    b_tmp_data.data = &tmp_data[0];
    b_tmp_data.size = &tmp_size[0];
    b_tmp_data.allocatedSize = 29;
    b_tmp_data.numDimensions = 2;
    b_tmp_data.canFreeData = false;
    n_freqz_cg(&b_tmp_data, w, Fs, d3);
    q_freqz_cg(*(double (*)[57])&d_u[0], w, Fs, combinedResponse);
    p_freqz_cg(*(double (*)[13])&c_u[0], w, Fs, d2);
    iy = combinedResponse->size[0] * combinedResponse->size[1];
    ix = combinedResponse->size[0] * combinedResponse->size[1];
    combinedResponse->size[0] = 1;
    emxEnsureCapacity_creal_T(combinedResponse, ix);
    ix = iy - 1;
    for (iy = 0; iy <= ix; iy++) {
      combinedResponse_re = combinedResponse->data[iy].re * d2->data[iy].re -
        combinedResponse->data[iy].im * d2->data[iy].im;
      combinedResponse_im = combinedResponse->data[iy].re * d2->data[iy].im +
        combinedResponse->data[iy].im * d2->data[iy].re;
      d2_re = d3->data[iy].re;
      d2_im = d3->data[iy].im;
      combinedResponse->data[iy].re = combinedResponse_re * d2_re -
        combinedResponse_im * d2_im;
      combinedResponse->data[iy].im = combinedResponse_re * d2_im +
        combinedResponse_im * d2_re;
    }

    pre_FIR_us = 8;
    break;

   case 8:
    /*  Dec/Int3 */
    c_us(dec_int3_coeff_data, dec_int3_coeff_size, tmp_data, tmp_size);
    b_tmp_data.data = &tmp_data[0];
    b_tmp_data.size = &tmp_size[0];
    b_tmp_data.allocatedSize = 29;
    b_tmp_data.numDimensions = 2;
    b_tmp_data.canFreeData = false;
    n_freqz_cg(&b_tmp_data, w, Fs, combinedResponse);
    pre_FIR_us = 3;

    /*  RECHECK ALL DEC BY 3     */
    break;

   case 9:
    /*  Dec/Int3,Hb1 */
    memset(&e_u[0], 0, 45U * sizeof(double));
    ix = 0;
    iy = 0;
    for (k = 0; k < 15; k++) {
      e_u[iy] = hb1_coeff[ix];
      ix++;
      iy += 3;
    }

    c_us(dec_int3_coeff_data, dec_int3_coeff_size, tmp_data, tmp_size);
    b_tmp_data.data = &tmp_data[0];
    b_tmp_data.size = &tmp_size[0];
    b_tmp_data.allocatedSize = 29;
    b_tmp_data.numDimensions = 2;
    b_tmp_data.canFreeData = false;
    n_freqz_cg(&b_tmp_data, w, Fs, d2);
    r_freqz_cg(*(double (*)[43])&e_u[0], w, Fs, combinedResponse);
    iy = combinedResponse->size[0] * combinedResponse->size[1];
    ix = combinedResponse->size[0] * combinedResponse->size[1];
    combinedResponse->size[0] = 1;
    emxEnsureCapacity_creal_T(combinedResponse, ix);
    ix = iy - 1;
    for (iy = 0; iy <= ix; iy++) {
      combinedResponse_re = combinedResponse->data[iy].re;
      combinedResponse_im = combinedResponse->data[iy].im;
      d2_re = d2->data[iy].re;
      d2_im = d2->data[iy].im;
      combinedResponse->data[iy].re = combinedResponse_re * d2_re -
        combinedResponse_im * d2_im;
      combinedResponse->data[iy].im = combinedResponse_re * d2_im +
        combinedResponse_im * d2_re;
    }

    pre_FIR_us = 6;
    break;

   case 10:
    /*  Dec/Int3,Hb2 */
    memset(&f_u[0], 0, 21U * sizeof(double));
    ix = 0;
    iy = 0;
    for (k = 0; k < 7; k++) {
      f_u[iy] = hb2_coeff[ix];
      ix++;
      iy += 3;
    }

    c_us(dec_int3_coeff_data, dec_int3_coeff_size, tmp_data, tmp_size);
    b_tmp_data.data = &tmp_data[0];
    b_tmp_data.size = &tmp_size[0];
    b_tmp_data.allocatedSize = 29;
    b_tmp_data.numDimensions = 2;
    b_tmp_data.canFreeData = false;
    n_freqz_cg(&b_tmp_data, w, Fs, d2);
    s_freqz_cg(*(double (*)[19])&f_u[0], w, Fs, combinedResponse);
    iy = combinedResponse->size[0] * combinedResponse->size[1];
    ix = combinedResponse->size[0] * combinedResponse->size[1];
    combinedResponse->size[0] = 1;
    emxEnsureCapacity_creal_T(combinedResponse, ix);
    ix = iy - 1;
    for (iy = 0; iy <= ix; iy++) {
      combinedResponse_re = combinedResponse->data[iy].re;
      combinedResponse_im = combinedResponse->data[iy].im;
      d2_re = d2->data[iy].re;
      d2_im = d2->data[iy].im;
      combinedResponse->data[iy].re = combinedResponse_re * d2_re -
        combinedResponse_im * d2_im;
      combinedResponse->data[iy].im = combinedResponse_re * d2_im +
        combinedResponse_im * d2_re;
    }

    pre_FIR_us = 6;
    break;

   case 11:
    /*  Dec/Int3,Hb2,Hb1 {Hm4,Hm2c34,Hm1} */
    memset(&g_u[0], 0, 90U * sizeof(double));
    ix = 0;
    iy = 0;
    for (k = 0; k < 15; k++) {
      g_u[iy] = hb1_coeff[ix];
      ix++;
      iy += 6;
    }

    memset(&f_u[0], 0, 21U * sizeof(double));
    ix = 0;
    iy = 0;
    for (k = 0; k < 7; k++) {
      f_u[iy] = hb2_coeff[ix];
      ix++;
      iy += 3;
    }

    c_us(dec_int3_coeff_data, dec_int3_coeff_size, tmp_data, tmp_size);
    b_tmp_data.data = &tmp_data[0];
    b_tmp_data.size = &tmp_size[0];
    b_tmp_data.allocatedSize = 29;
    b_tmp_data.numDimensions = 2;
    b_tmp_data.canFreeData = false;
    n_freqz_cg(&b_tmp_data, w, Fs, d3);
    t_freqz_cg(*(double (*)[85])&g_u[0], w, Fs, combinedResponse);
    s_freqz_cg(*(double (*)[19])&f_u[0], w, Fs, d2);
    iy = combinedResponse->size[0] * combinedResponse->size[1];
    ix = combinedResponse->size[0] * combinedResponse->size[1];
    combinedResponse->size[0] = 1;
    emxEnsureCapacity_creal_T(combinedResponse, ix);
    ix = iy - 1;
    for (iy = 0; iy <= ix; iy++) {
      combinedResponse_re = combinedResponse->data[iy].re * d2->data[iy].re -
        combinedResponse->data[iy].im * d2->data[iy].im;
      combinedResponse_im = combinedResponse->data[iy].re * d2->data[iy].im +
        combinedResponse->data[iy].im * d2->data[iy].re;
      d2_re = d3->data[iy].re;
      d2_im = d3->data[iy].im;
      combinedResponse->data[iy].re = combinedResponse_re * d2_re -
        combinedResponse_im * d2_im;
      combinedResponse->data[iy].im = combinedResponse_re * d2_im +
        combinedResponse_im * d2_re;
    }

    pre_FIR_us = 12;
    break;
  }

  emxFree_creal_T(&d3);

  /*  Add filter extra filter to end of cascade */
  if (extraTaps->size[1] != 0) {
    emxInit_real_T(&i_u, 2);
    vlenx = extraTaps->size[1];
    ix = pre_FIR_us * extraTaps->size[1];
    iy = i_u->size[0] * i_u->size[1];
    i_u->size[0] = 1;
    i_u->size[1] = ix;
    emxEnsureCapacity_real_T(i_u, iy);
    for (iy = 0; iy < ix; iy++) {
      i_u->data[iy] = 0.0;
    }

    ix = 0;
    iy = 0;
    for (k = 0; k < vlenx; k++) {
      i_u->data[iy] = extraTaps->data[ix];
      ix++;
      iy += pre_FIR_us;
    }

    iy = (i_u->size[1] - pre_FIR_us) + 1;
    ix = i_u->size[0] * i_u->size[1];
    if (1 > iy) {
      i_u->size[1] = 0;
    } else {
      i_u->size[1] = iy;
    }

    emxEnsureCapacity_real_T(i_u, ix);
    n_freqz_cg(i_u, w, Fs, d2);
    iy = combinedResponse->size[0] * combinedResponse->size[1];
    ix = combinedResponse->size[0] * combinedResponse->size[1];
    combinedResponse->size[0] = 1;
    emxEnsureCapacity_creal_T(combinedResponse, ix);
    ix = iy - 1;
    emxFree_real_T(&i_u);
    for (iy = 0; iy <= ix; iy++) {
      combinedResponse_re = combinedResponse->data[iy].re;
      combinedResponse_im = combinedResponse->data[iy].im;
      d2_re = d2->data[iy].re;
      d2_im = d2->data[iy].im;
      combinedResponse->data[iy].re = combinedResponse_re * d2_re -
        combinedResponse_im * d2_im;
      combinedResponse->data[iy].im = combinedResponse_re * d2_im +
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
  int i13;
  int k;
  double p_re;
  int i14;
  double x_im;
  if (p_size[1] != 0) {
    for (i13 = 0; i13 < 2048; i13++) {
      y[i13].re = p_data[0];
      y[i13].im = 0.0;
    }

    i13 = p_size[1];
    for (k = 0; k <= i13 - 2; k++) {
      p_re = p_data[k + 1];
      for (i14 = 0; i14 < 2048; i14++) {
        x_im = x[i14].re * y[i14].im + x[i14].im * y[i14].re;
        y[i14].re = (x[i14].re * y[i14].re - x[i14].im * y[i14].im) + p_re;
        y[i14].im = x_im;
      }
    }
  }
}

/*
 * Arguments    : const emxArray_real_T *a
 *                emxArray_real_T *y
 * Return Type  : void
 */
static void c_power(const emxArray_real_T *a, emxArray_real_T *y)
{
  int nx;
  int k;
  nx = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = a->size[1];
  emxEnsureCapacity_real_T(y, nx);
  nx = a->size[1];
  for (k = 0; k < nx; k++) {
    y->data[k] = rt_powd_snf(a->data[k], 3.0);
  }
}

/*
 * Arguments    : creal_T *x
 * Return Type  : void
 */
static void c_sqrt(creal_T *x)
{
  double xr;
  double xi;
  double yr;
  double absxr;
  xr = x->re;
  xi = x->im;
  if (xi == 0.0) {
    if (xr < 0.0) {
      yr = 0.0;
      xr = sqrt(-xr);
    } else {
      yr = sqrt(xr);
      xr = 0.0;
    }
  } else if (xr == 0.0) {
    if (xi < 0.0) {
      yr = sqrt(-xi / 2.0);
      xr = -yr;
    } else {
      yr = sqrt(xi / 2.0);
      xr = yr;
    }
  } else if (rtIsNaN(xr)) {
    yr = xr;
  } else if (rtIsNaN(xi)) {
    yr = xi;
    xr = xi;
  } else if (rtIsInf(xi)) {
    yr = fabs(xi);
    xr = xi;
  } else if (rtIsInf(xr)) {
    if (xr < 0.0) {
      yr = 0.0;
      xr = xi * -xr;
    } else {
      yr = xr;
      xr = 0.0;
    }
  } else {
    absxr = fabs(xr);
    yr = fabs(xi);
    if ((absxr > 4.4942328371557893E+307) || (yr > 4.4942328371557893E+307)) {
      absxr *= 0.5;
      yr = rt_hypotd_snf(absxr, yr * 0.5);
      if (yr > absxr) {
        yr = sqrt(yr) * sqrt(1.0 + absxr / yr);
      } else {
        yr = sqrt(yr) * 1.4142135623730951;
      }
    } else {
      yr = sqrt((rt_hypotd_snf(absxr, yr) + absxr) * 0.5);
    }

    if (xr > 0.0) {
      xr = 0.5 * (xi / yr);
    } else {
      if (xi < 0.0) {
        xr = -yr;
      } else {
        xr = yr;
      }

      yr = 0.5 * (xi / xr);
    }
  }

  x->re = yr;
  x->im = xr;
}

/*
 * Arguments    : const char a[2]
 * Return Type  : bool
 */
static bool c_strcmp(const char a[2])
{
  bool b_bool;
  int kstr;
  int exitg1;
  static const char cv32[2] = { 'T', 'x' };

  b_bool = false;
  kstr = 0;
  do {
    exitg1 = 0;
    if (kstr < 2) {
      if (a[kstr] != cv32[kstr]) {
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
 * Arguments    : const double x[128]
 * Return Type  : double
 */
static double c_sum(const double x[128])
{
  double y;
  int k;
  y = x[0];
  for (k = 0; k < 127; k++) {
    y += x[k + 1];
  }

  return y;
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
  int vlenx;
  int ix;
  int iy;
  int k;
  vlenx = o_size[1];
  u_size[0] = 1;
  u_size[1] = (signed char)o_size[1];
  ix = (signed char)o_size[1];
  if (0 <= ix - 1) {
    memset(&u_data[0], 0, (unsigned int)(ix * (int)sizeof(double)));
  }

  ix = 0;
  iy = 0;
  for (k = 0; k < vlenx; k++) {
    u_data[iy] = o_data[ix];
    ix++;
    iy++;
  }
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
 * Arguments    : const double x[128]
 *                double y[128]
 * Return Type  : void
 */
static void d_abs(const double x[128], double y[128])
{
  int k;
  for (k = 0; k < 128; k++) {
    y[k] = fabs(x[k]);
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
static void d_analogresp(const emxArray_real_T *f, double Fconverter, const
  double b1_data[], const int b1_size[2], const emxArray_creal_T *a1, const
  double b2_data[], const int b2_size[2], const emxArray_creal_T *a2,
  emxArray_creal_T *abc)
{
  emxArray_real_T *x;
  int i63;
  int loop_ub;
  emxArray_real_T *r18;
  emxArray_creal_T *r19;
  double x_re;
  double x_im;
  double re;
  double im;
  emxInit_real_T(&x, 2);
  i63 = x->size[0] * x->size[1];
  x->size[0] = 1;
  x->size[1] = f->size[1];
  emxEnsureCapacity_real_T(x, i63);
  loop_ub = f->size[0] * f->size[1];
  for (i63 = 0; i63 < loop_ub; i63++) {
    x->data[i63] = f->data[i63] / Fconverter;
  }

  i63 = x->size[1];
  for (loop_ub = 0; loop_ub < i63; loop_ub++) {
    if (fabs(x->data[loop_ub]) < 1.0020841800044864E-292) {
      x->data[loop_ub] = 1.0;
    } else {
      x->data[loop_ub] *= 3.1415926535897931;
      x->data[loop_ub] = sin(x->data[loop_ub]) / x->data[loop_ub];
    }
  }

  emxInit_real_T(&r18, 2);
  i63 = r18->size[0] * r18->size[1];
  r18->size[0] = 1;
  r18->size[1] = f->size[1];
  emxEnsureCapacity_real_T(r18, i63);
  loop_ub = f->size[0] * f->size[1];
  for (i63 = 0; i63 < loop_ub; i63++) {
    r18->data[i63] = 6.2831853071795862 * f->data[i63];
  }

  b_freqs_cg(b1_data, b1_size, a1, r18, abc);
  i63 = r18->size[0] * r18->size[1];
  r18->size[0] = 1;
  r18->size[1] = f->size[1];
  emxEnsureCapacity_real_T(r18, i63);
  loop_ub = f->size[0] * f->size[1];
  for (i63 = 0; i63 < loop_ub; i63++) {
    r18->data[i63] = 6.2831853071795862 * f->data[i63];
  }

  emxInit_creal_T(&r19, 2);
  b_freqs_cg(b2_data, b2_size, a2, r18, r19);
  i63 = x->size[0] * x->size[1];
  loop_ub = abc->size[0] * abc->size[1];
  abc->size[0] = 1;
  abc->size[1] = x->size[1];
  emxEnsureCapacity_creal_T(abc, loop_ub);
  loop_ub = i63 - 1;
  emxFree_real_T(&r18);
  for (i63 = 0; i63 <= loop_ub; i63++) {
    x_re = x->data[i63] * abc->data[i63].re;
    x_im = x->data[i63] * abc->data[i63].im;
    re = r19->data[i63].re;
    im = r19->data[i63].im;
    abc->data[i63].re = x_re * re - x_im * im;
    abc->data[i63].im = x_re * im + x_im * re;
  }

  emxFree_real_T(&x);
  emxFree_creal_T(&r19);
}

/*
 * Make b a row
 * Arguments    : double b_data[]
 *                int b_size[2]
 *                const struct_T *options
 *                creal_T h[2048]
 *                double w[2048]
 * Return Type  : void
 */
static void d_firfreqz(double b_data[], int b_size[2], const struct_T *options,
  creal_T h[2048], double w[2048])
{
  double b_b_data[29];
  int i75;
  creal_T dcv4[2048];
  double re_tmp;
  double re;
  creal_T s_tmp[2048];
  double h_re;
  double brm;

  /* -------------------------------------------------------------------------- */
  if (0 <= b_size[1] - 1) {
    memcpy(&b_b_data[0], &b_data[0], (unsigned int)(b_size[1] * (int)sizeof
            (double)));
  }

  b_size[0] = 1;
  if (0 <= b_size[1] - 1) {
    memcpy(&b_data[0], &b_b_data[0], (unsigned int)(b_size[1] * (int)sizeof
            (double)));
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
  for (i75 = 0; i75 < 2048; i75++) {
    w[i75] = options->w[i75];
    re_tmp = 6.2831853071795862 * options->w[i75] / options->Fs;
    re = re_tmp * 0.0;
    s_tmp[i75].re = re;
    s_tmp[i75].im = re_tmp;
    dcv4[i75].re = re;
    dcv4[i75].im = re_tmp;
  }

  b_exp(dcv4);
  c_polyval(b_data, b_size, dcv4, h);
  for (i75 = 0; i75 < 2048; i75++) {
    s_tmp[i75].re *= (double)b_size[1] - 1.0;
    s_tmp[i75].im *= (double)b_size[1] - 1.0;
  }

  b_exp(s_tmp);
  for (i75 = 0; i75 < 2048; i75++) {
    h_re = h[i75].re;
    if (s_tmp[i75].im == 0.0) {
      if (h[i75].im == 0.0) {
        h[i75].re /= s_tmp[i75].re;
        h[i75].im = 0.0;
      } else if (h[i75].re == 0.0) {
        h[i75].re = 0.0;
        h[i75].im /= s_tmp[i75].re;
      } else {
        h[i75].re /= s_tmp[i75].re;
        h[i75].im /= s_tmp[i75].re;
      }
    } else if (s_tmp[i75].re == 0.0) {
      if (h[i75].re == 0.0) {
        h[i75].re = h[i75].im / s_tmp[i75].im;
        h[i75].im = 0.0;
      } else if (h[i75].im == 0.0) {
        h[i75].re = 0.0;
        h[i75].im = -(h_re / s_tmp[i75].im);
      } else {
        h[i75].re = h[i75].im / s_tmp[i75].im;
        h[i75].im = -(h_re / s_tmp[i75].im);
      }
    } else {
      brm = fabs(s_tmp[i75].re);
      re = fabs(s_tmp[i75].im);
      if (brm > re) {
        re = s_tmp[i75].im / s_tmp[i75].re;
        re_tmp = s_tmp[i75].re + re * s_tmp[i75].im;
        h[i75].re = (h[i75].re + re * h[i75].im) / re_tmp;
        h[i75].im = (h[i75].im - re * h_re) / re_tmp;
      } else if (re == brm) {
        if (s_tmp[i75].re > 0.0) {
          re = 0.5;
        } else {
          re = -0.5;
        }

        if (s_tmp[i75].im > 0.0) {
          re_tmp = 0.5;
        } else {
          re_tmp = -0.5;
        }

        h[i75].re = (h[i75].re * re + h[i75].im * re_tmp) / brm;
        h[i75].im = (h[i75].im * re - h_re * re_tmp) / brm;
      } else {
        re = s_tmp[i75].re / s_tmp[i75].im;
        re_tmp = s_tmp[i75].im + re * s_tmp[i75].re;
        h[i75].re = (re * h[i75].re + h[i75].im) / re_tmp;
        h[i75].im = (re * h[i75].im - h_re) / re_tmp;
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
  int loop_ub;
  static const char cv18[8] = { 'o', 'n', 'e', 's', 'i', 'd', 'e', 'd' };

  int b_b_size[2];
  static const char cv19[7] = { 'o', 'm', 'i', 't', 't', 'e', 'd' };

  double b_b_data[29];
  double b_w[2048];

  /*  Cast to enforce precision rules */
  options.Fs = Fs;
  memcpy(&options.nfft[0], &w[0], sizeof(double) << 11);
  memcpy(&options.w[0], &w[0], sizeof(double) << 11);

  /*  Remaining are default or for advanced use */
  options.fvflag = 1.0;
  for (loop_ub = 0; loop_ub < 8; loop_ub++) {
    options.range[loop_ub] = cv18[loop_ub];
  }

  options.centerdc = 0.0;
  for (loop_ub = 0; loop_ub < 7; loop_ub++) {
    options.configlevel[loop_ub] = cv19[loop_ub];
  }

  b_b_size[0] = 1;
  b_b_size[1] = b_size[1];
  loop_ub = b_size[0] * b_size[1];
  if (0 <= loop_ub - 1) {
    memcpy(&b_b_data[0], &b_data[0], (unsigned int)(loop_ub * (int)sizeof(double)));
  }

  d_firfreqz(b_b_data, b_b_size, &options, hh, b_w);

  /*  Generate the default structure to pass to freqzplot */
  /*  If rad/sample, Fs is empty */
}

/*
 * Arguments    : const double p[29]
 *                const creal_T x[2048]
 *                creal_T y[2048]
 * Return Type  : void
 */
static void d_polyval(const double p[29], const creal_T x[2048], creal_T y[2048])
{
  int i15;
  int k;
  double p_re;
  double x_im;
  for (i15 = 0; i15 < 2048; i15++) {
    y[i15].re = p[0];
    y[i15].im = 0.0;
  }

  for (k = 0; k < 28; k++) {
    p_re = p[k + 1];
    for (i15 = 0; i15 < 2048; i15++) {
      x_im = x[i15].re * y[i15].im + x[i15].im * y[i15].re;
      y[i15].re = (x[i15].re * y[i15].re - x[i15].im * y[i15].im) + p_re;
      y[i15].im = x_im;
    }
  }
}

/*
 * Arguments    : const double o[15]
 *                double u[29]
 * Return Type  : void
 */
static void d_us(const double o[15], double u[29])
{
  double b_u[30];
  int ix;
  int iy;
  int k;
  memset(&b_u[0], 0, 30U * sizeof(double));
  ix = 0;
  iy = 0;
  for (k = 0; k < 15; k++) {
    b_u[iy] = o[ix];
    ix++;
    iy += 2;
  }

  memcpy(&u[0], &b_u[0], 29U * sizeof(double));
}

/*
 * Arguments    : double ydb
 * Return Type  : double
 */
static double db2mag(double ydb)
{
  return rt_powd_snf(10.0, ydb / 20.0);
}

/*
 * Codegen workaround for fixed fi call requirements
 * Arguments    : const emxArray_real_T *tap_store
 *                double i
 *                double M
 *                emxArray_real_T *taps
 * Return Type  : void
 */
static void determineBestFractionLength(const emxArray_real_T *tap_store, double
  i, double M, emxArray_real_T *taps)
{
  int loop_ub;
  emxArray_real_T *org;
  int i60;
  emxArray_real_T *r;
  double u;
  double v;
  emxArray_real_T *b_r;
  short i61;
  emxArray_real_T *r15;
  double e[16];
  int idx;
  bool exitg1;
  if (1.0 > M) {
    loop_ub = 0;
  } else {
    loop_ub = (int)M;
  }

  emxInit_real_T(&org, 2);
  i60 = org->size[0] * org->size[1];
  org->size[0] = 1;
  org->size[1] = loop_ub;
  emxEnsureCapacity_real_T(org, i60);
  for (i60 = 0; i60 < loop_ub; i60++) {
    org->data[i60] = tap_store->data[((int)i + tap_store->size[0] * i60) - 1];
  }

  emxInit_real_T(&r, 2);
  i60 = r->size[0] * r->size[1];
  r->size[0] = 16;
  loop_ub = (int)M;
  r->size[1] = loop_ub;
  emxEnsureCapacity_real_T(r, i60);
  loop_ub <<= 4;
  for (i60 = 0; i60 < loop_ub; i60++) {
    r->data[i60] = 0.0;
  }

  if (1.0 > M) {
    loop_ub = -1;
  } else {
    loop_ub = (int)M - 1;
  }

  for (i60 = 0; i60 <= loop_ub; i60++) {
    u = tap_store->data[((int)i + tap_store->size[0] * i60) - 1] * 2.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i61 = (short)u;
      } else {
        i61 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i61 = MAX_int16_T;
    } else {
      i61 = 0;
    }

    r->data[i60 << 4] = (double)i61 * 0.5;
  }

  if (1.0 > M) {
    loop_ub = 0;
  } else {
    loop_ub = (int)M;
  }

  emxInit_real_T(&b_r, 2);
  i60 = b_r->size[0] * b_r->size[1];
  b_r->size[0] = 1;
  b_r->size[1] = loop_ub;
  emxEnsureCapacity_real_T(b_r, i60);
  for (i60 = 0; i60 < loop_ub; i60++) {
    b_r->data[i60] = r->data[i60 << 4] - org->data[i60];
  }

  emxInit_real_T(&r15, 2);
  c_abs(b_r, r15);
  e[0] = b_sum(r15);
  if (1.0 > M) {
    loop_ub = -1;
  } else {
    loop_ub = (int)M - 1;
  }

  for (i60 = 0; i60 <= loop_ub; i60++) {
    u = tap_store->data[((int)i + tap_store->size[0] * i60) - 1] * 4.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i61 = (short)u;
      } else {
        i61 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i61 = MAX_int16_T;
    } else {
      i61 = 0;
    }

    r->data[1 + (i60 << 4)] = (double)i61 * 0.25;
  }

  if (1.0 > M) {
    loop_ub = 0;
  } else {
    loop_ub = (int)M;
  }

  i60 = b_r->size[0] * b_r->size[1];
  b_r->size[0] = 1;
  b_r->size[1] = loop_ub;
  emxEnsureCapacity_real_T(b_r, i60);
  for (i60 = 0; i60 < loop_ub; i60++) {
    b_r->data[i60] = r->data[1 + (i60 << 4)] - org->data[i60];
  }

  c_abs(b_r, r15);
  e[1] = b_sum(r15);
  if (1.0 > M) {
    loop_ub = -1;
  } else {
    loop_ub = (int)M - 1;
  }

  for (i60 = 0; i60 <= loop_ub; i60++) {
    u = tap_store->data[((int)i + tap_store->size[0] * i60) - 1] * 8.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i61 = (short)u;
      } else {
        i61 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i61 = MAX_int16_T;
    } else {
      i61 = 0;
    }

    r->data[2 + (i60 << 4)] = (double)i61 * 0.125;
  }

  if (1.0 > M) {
    loop_ub = 0;
  } else {
    loop_ub = (int)M;
  }

  i60 = b_r->size[0] * b_r->size[1];
  b_r->size[0] = 1;
  b_r->size[1] = loop_ub;
  emxEnsureCapacity_real_T(b_r, i60);
  for (i60 = 0; i60 < loop_ub; i60++) {
    b_r->data[i60] = r->data[2 + (i60 << 4)] - org->data[i60];
  }

  c_abs(b_r, r15);
  e[2] = b_sum(r15);
  if (1.0 > M) {
    loop_ub = -1;
  } else {
    loop_ub = (int)M - 1;
  }

  for (i60 = 0; i60 <= loop_ub; i60++) {
    u = tap_store->data[((int)i + tap_store->size[0] * i60) - 1] * 16.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i61 = (short)u;
      } else {
        i61 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i61 = MAX_int16_T;
    } else {
      i61 = 0;
    }

    r->data[3 + (i60 << 4)] = (double)i61 * 0.0625;
  }

  if (1.0 > M) {
    loop_ub = 0;
  } else {
    loop_ub = (int)M;
  }

  i60 = b_r->size[0] * b_r->size[1];
  b_r->size[0] = 1;
  b_r->size[1] = loop_ub;
  emxEnsureCapacity_real_T(b_r, i60);
  for (i60 = 0; i60 < loop_ub; i60++) {
    b_r->data[i60] = r->data[3 + (i60 << 4)] - org->data[i60];
  }

  c_abs(b_r, r15);
  e[3] = b_sum(r15);
  if (1.0 > M) {
    loop_ub = -1;
  } else {
    loop_ub = (int)M - 1;
  }

  for (i60 = 0; i60 <= loop_ub; i60++) {
    u = tap_store->data[((int)i + tap_store->size[0] * i60) - 1] * 32.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i61 = (short)u;
      } else {
        i61 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i61 = MAX_int16_T;
    } else {
      i61 = 0;
    }

    r->data[4 + (i60 << 4)] = (double)i61 * 0.03125;
  }

  if (1.0 > M) {
    loop_ub = 0;
  } else {
    loop_ub = (int)M;
  }

  i60 = b_r->size[0] * b_r->size[1];
  b_r->size[0] = 1;
  b_r->size[1] = loop_ub;
  emxEnsureCapacity_real_T(b_r, i60);
  for (i60 = 0; i60 < loop_ub; i60++) {
    b_r->data[i60] = r->data[4 + (i60 << 4)] - org->data[i60];
  }

  c_abs(b_r, r15);
  e[4] = b_sum(r15);
  if (1.0 > M) {
    loop_ub = -1;
  } else {
    loop_ub = (int)M - 1;
  }

  for (i60 = 0; i60 <= loop_ub; i60++) {
    u = tap_store->data[((int)i + tap_store->size[0] * i60) - 1] * 64.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i61 = (short)u;
      } else {
        i61 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i61 = MAX_int16_T;
    } else {
      i61 = 0;
    }

    r->data[5 + (i60 << 4)] = (double)i61 * 0.015625;
  }

  if (1.0 > M) {
    loop_ub = 0;
  } else {
    loop_ub = (int)M;
  }

  i60 = b_r->size[0] * b_r->size[1];
  b_r->size[0] = 1;
  b_r->size[1] = loop_ub;
  emxEnsureCapacity_real_T(b_r, i60);
  for (i60 = 0; i60 < loop_ub; i60++) {
    b_r->data[i60] = r->data[5 + (i60 << 4)] - org->data[i60];
  }

  c_abs(b_r, r15);
  e[5] = b_sum(r15);
  if (1.0 > M) {
    loop_ub = -1;
  } else {
    loop_ub = (int)M - 1;
  }

  for (i60 = 0; i60 <= loop_ub; i60++) {
    u = tap_store->data[((int)i + tap_store->size[0] * i60) - 1] * 128.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i61 = (short)u;
      } else {
        i61 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i61 = MAX_int16_T;
    } else {
      i61 = 0;
    }

    r->data[6 + (i60 << 4)] = (double)i61 * 0.0078125;
  }

  if (1.0 > M) {
    loop_ub = 0;
  } else {
    loop_ub = (int)M;
  }

  i60 = b_r->size[0] * b_r->size[1];
  b_r->size[0] = 1;
  b_r->size[1] = loop_ub;
  emxEnsureCapacity_real_T(b_r, i60);
  for (i60 = 0; i60 < loop_ub; i60++) {
    b_r->data[i60] = r->data[6 + (i60 << 4)] - org->data[i60];
  }

  c_abs(b_r, r15);
  e[6] = b_sum(r15);
  if (1.0 > M) {
    loop_ub = -1;
  } else {
    loop_ub = (int)M - 1;
  }

  for (i60 = 0; i60 <= loop_ub; i60++) {
    u = tap_store->data[((int)i + tap_store->size[0] * i60) - 1] * 256.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i61 = (short)u;
      } else {
        i61 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i61 = MAX_int16_T;
    } else {
      i61 = 0;
    }

    r->data[7 + (i60 << 4)] = (double)i61 * 0.00390625;
  }

  if (1.0 > M) {
    loop_ub = 0;
  } else {
    loop_ub = (int)M;
  }

  i60 = b_r->size[0] * b_r->size[1];
  b_r->size[0] = 1;
  b_r->size[1] = loop_ub;
  emxEnsureCapacity_real_T(b_r, i60);
  for (i60 = 0; i60 < loop_ub; i60++) {
    b_r->data[i60] = r->data[7 + (i60 << 4)] - org->data[i60];
  }

  c_abs(b_r, r15);
  e[7] = b_sum(r15);
  if (1.0 > M) {
    loop_ub = -1;
  } else {
    loop_ub = (int)M - 1;
  }

  for (i60 = 0; i60 <= loop_ub; i60++) {
    u = tap_store->data[((int)i + tap_store->size[0] * i60) - 1] * 512.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i61 = (short)u;
      } else {
        i61 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i61 = MAX_int16_T;
    } else {
      i61 = 0;
    }

    r->data[8 + (i60 << 4)] = (double)i61 * 0.001953125;
  }

  if (1.0 > M) {
    loop_ub = 0;
  } else {
    loop_ub = (int)M;
  }

  i60 = b_r->size[0] * b_r->size[1];
  b_r->size[0] = 1;
  b_r->size[1] = loop_ub;
  emxEnsureCapacity_real_T(b_r, i60);
  for (i60 = 0; i60 < loop_ub; i60++) {
    b_r->data[i60] = r->data[8 + (i60 << 4)] - org->data[i60];
  }

  c_abs(b_r, r15);
  e[8] = b_sum(r15);
  if (1.0 > M) {
    loop_ub = -1;
  } else {
    loop_ub = (int)M - 1;
  }

  for (i60 = 0; i60 <= loop_ub; i60++) {
    u = tap_store->data[((int)i + tap_store->size[0] * i60) - 1] * 1024.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i61 = (short)u;
      } else {
        i61 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i61 = MAX_int16_T;
    } else {
      i61 = 0;
    }

    r->data[9 + (i60 << 4)] = (double)i61 * 0.0009765625;
  }

  if (1.0 > M) {
    loop_ub = 0;
  } else {
    loop_ub = (int)M;
  }

  i60 = b_r->size[0] * b_r->size[1];
  b_r->size[0] = 1;
  b_r->size[1] = loop_ub;
  emxEnsureCapacity_real_T(b_r, i60);
  for (i60 = 0; i60 < loop_ub; i60++) {
    b_r->data[i60] = r->data[9 + (i60 << 4)] - org->data[i60];
  }

  c_abs(b_r, r15);
  e[9] = b_sum(r15);
  if (1.0 > M) {
    loop_ub = -1;
  } else {
    loop_ub = (int)M - 1;
  }

  for (i60 = 0; i60 <= loop_ub; i60++) {
    u = tap_store->data[((int)i + tap_store->size[0] * i60) - 1] * 2048.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i61 = (short)u;
      } else {
        i61 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i61 = MAX_int16_T;
    } else {
      i61 = 0;
    }

    r->data[10 + (i60 << 4)] = (double)i61 * 0.00048828125;
  }

  if (1.0 > M) {
    loop_ub = 0;
  } else {
    loop_ub = (int)M;
  }

  i60 = b_r->size[0] * b_r->size[1];
  b_r->size[0] = 1;
  b_r->size[1] = loop_ub;
  emxEnsureCapacity_real_T(b_r, i60);
  for (i60 = 0; i60 < loop_ub; i60++) {
    b_r->data[i60] = r->data[10 + (i60 << 4)] - org->data[i60];
  }

  c_abs(b_r, r15);
  e[10] = b_sum(r15);
  if (1.0 > M) {
    loop_ub = -1;
  } else {
    loop_ub = (int)M - 1;
  }

  for (i60 = 0; i60 <= loop_ub; i60++) {
    u = tap_store->data[((int)i + tap_store->size[0] * i60) - 1] * 4096.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i61 = (short)u;
      } else {
        i61 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i61 = MAX_int16_T;
    } else {
      i61 = 0;
    }

    r->data[11 + (i60 << 4)] = (double)i61 * 0.000244140625;
  }

  if (1.0 > M) {
    loop_ub = 0;
  } else {
    loop_ub = (int)M;
  }

  i60 = b_r->size[0] * b_r->size[1];
  b_r->size[0] = 1;
  b_r->size[1] = loop_ub;
  emxEnsureCapacity_real_T(b_r, i60);
  for (i60 = 0; i60 < loop_ub; i60++) {
    b_r->data[i60] = r->data[11 + (i60 << 4)] - org->data[i60];
  }

  c_abs(b_r, r15);
  e[11] = b_sum(r15);
  if (1.0 > M) {
    loop_ub = -1;
  } else {
    loop_ub = (int)M - 1;
  }

  for (i60 = 0; i60 <= loop_ub; i60++) {
    u = tap_store->data[((int)i + tap_store->size[0] * i60) - 1] * 8192.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i61 = (short)u;
      } else {
        i61 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i61 = MAX_int16_T;
    } else {
      i61 = 0;
    }

    r->data[12 + (i60 << 4)] = (double)i61 * 0.0001220703125;
  }

  if (1.0 > M) {
    loop_ub = 0;
  } else {
    loop_ub = (int)M;
  }

  i60 = b_r->size[0] * b_r->size[1];
  b_r->size[0] = 1;
  b_r->size[1] = loop_ub;
  emxEnsureCapacity_real_T(b_r, i60);
  for (i60 = 0; i60 < loop_ub; i60++) {
    b_r->data[i60] = r->data[12 + (i60 << 4)] - org->data[i60];
  }

  c_abs(b_r, r15);
  e[12] = b_sum(r15);
  if (1.0 > M) {
    loop_ub = -1;
  } else {
    loop_ub = (int)M - 1;
  }

  for (i60 = 0; i60 <= loop_ub; i60++) {
    u = tap_store->data[((int)i + tap_store->size[0] * i60) - 1] * 16384.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i61 = (short)u;
      } else {
        i61 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i61 = MAX_int16_T;
    } else {
      i61 = 0;
    }

    r->data[13 + (i60 << 4)] = (double)i61 * 6.103515625E-5;
  }

  if (1.0 > M) {
    loop_ub = 0;
  } else {
    loop_ub = (int)M;
  }

  i60 = b_r->size[0] * b_r->size[1];
  b_r->size[0] = 1;
  b_r->size[1] = loop_ub;
  emxEnsureCapacity_real_T(b_r, i60);
  for (i60 = 0; i60 < loop_ub; i60++) {
    b_r->data[i60] = r->data[13 + (i60 << 4)] - org->data[i60];
  }

  c_abs(b_r, r15);
  e[13] = b_sum(r15);
  if (1.0 > M) {
    loop_ub = -1;
  } else {
    loop_ub = (int)M - 1;
  }

  for (i60 = 0; i60 <= loop_ub; i60++) {
    u = tap_store->data[((int)i + tap_store->size[0] * i60) - 1] * 32768.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i61 = (short)u;
      } else {
        i61 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i61 = MAX_int16_T;
    } else {
      i61 = 0;
    }

    r->data[14 + (i60 << 4)] = (double)i61 * 3.0517578125E-5;
  }

  if (1.0 > M) {
    loop_ub = 0;
  } else {
    loop_ub = (int)M;
  }

  i60 = b_r->size[0] * b_r->size[1];
  b_r->size[0] = 1;
  b_r->size[1] = loop_ub;
  emxEnsureCapacity_real_T(b_r, i60);
  for (i60 = 0; i60 < loop_ub; i60++) {
    b_r->data[i60] = r->data[14 + (i60 << 4)] - org->data[i60];
  }

  c_abs(b_r, r15);
  e[14] = b_sum(r15);
  if (1.0 > M) {
    loop_ub = -1;
  } else {
    loop_ub = (int)M - 1;
  }

  for (i60 = 0; i60 <= loop_ub; i60++) {
    u = tap_store->data[((int)i + tap_store->size[0] * i60) - 1] * 65536.0;
    v = fabs(u);
    if (v < 4.503599627370496E+15) {
      if (v >= 0.5) {
        u = floor(u + 0.5);
      } else {
        u *= 0.0;
      }
    }

    if (u < 32768.0) {
      if (u >= -32768.0) {
        i61 = (short)u;
      } else {
        i61 = MIN_int16_T;
      }
    } else if (u >= 32768.0) {
      i61 = MAX_int16_T;
    } else {
      i61 = 0;
    }

    r->data[15 + (i60 << 4)] = (double)i61 * 1.52587890625E-5;
  }

  if (1.0 > M) {
    loop_ub = 0;
  } else {
    loop_ub = (int)M;
  }

  i60 = b_r->size[0] * b_r->size[1];
  b_r->size[0] = 1;
  b_r->size[1] = loop_ub;
  emxEnsureCapacity_real_T(b_r, i60);
  for (i60 = 0; i60 < loop_ub; i60++) {
    b_r->data[i60] = r->data[15 + (i60 << 4)] - org->data[i60];
  }

  emxFree_real_T(&org);
  c_abs(b_r, r15);
  e[15] = b_sum(r15);
  emxFree_real_T(&b_r);
  emxFree_real_T(&r15);
  if (!rtIsNaN(e[0])) {
    idx = 1;
  } else {
    idx = 0;
    loop_ub = 2;
    exitg1 = false;
    while ((!exitg1) && (loop_ub < 17)) {
      if (!rtIsNaN(e[loop_ub - 1])) {
        idx = loop_ub;
        exitg1 = true;
      } else {
        loop_ub++;
      }
    }
  }

  if (idx == 0) {
    idx = 1;
  } else {
    u = e[idx - 1];
    i60 = idx + 1;
    for (loop_ub = i60; loop_ub < 17; loop_ub++) {
      v = e[loop_ub - 1];
      if (u > v) {
        u = v;
        idx = loop_ub;
      }
    }
  }

  if (1.0 > M) {
    loop_ub = 0;
  } else {
    loop_ub = (int)M;
  }

  i60 = taps->size[0] * taps->size[1];
  taps->size[0] = 1;
  taps->size[1] = loop_ub;
  emxEnsureCapacity_real_T(taps, i60);
  for (i60 = 0; i60 < loop_ub; i60++) {
    taps->data[i60] = r->data[(idx + (i60 << 4)) - 1];
  }

  emxFree_real_T(&r);
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
  bool quotientNeedsNegation;
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
 * Arguments    : const double b[29]
 *                const struct_T *options
 *                creal_T h[2048]
 *                double w[2048]
 * Return Type  : void
 */
static void e_firfreqz(const double b[29], const struct_T *options, creal_T h
  [2048], double w[2048])
{
  int i76;
  creal_T dcv5[2048];
  double re_tmp;
  double re;
  creal_T s_tmp[2048];
  double h_re;
  double brm;

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
  for (i76 = 0; i76 < 2048; i76++) {
    w[i76] = options->w[i76];
    re_tmp = 6.2831853071795862 * options->w[i76] / options->Fs;
    re = re_tmp * 0.0;
    s_tmp[i76].re = re;
    s_tmp[i76].im = re_tmp;
    dcv5[i76].re = re;
    dcv5[i76].im = re_tmp;
  }

  b_exp(dcv5);
  d_polyval(b, dcv5, h);
  for (i76 = 0; i76 < 2048; i76++) {
    s_tmp[i76].re *= 28.0;
    s_tmp[i76].im *= 28.0;
  }

  b_exp(s_tmp);
  for (i76 = 0; i76 < 2048; i76++) {
    h_re = h[i76].re;
    if (s_tmp[i76].im == 0.0) {
      if (h[i76].im == 0.0) {
        h[i76].re /= s_tmp[i76].re;
        h[i76].im = 0.0;
      } else if (h[i76].re == 0.0) {
        h[i76].re = 0.0;
        h[i76].im /= s_tmp[i76].re;
      } else {
        h[i76].re /= s_tmp[i76].re;
        h[i76].im /= s_tmp[i76].re;
      }
    } else if (s_tmp[i76].re == 0.0) {
      if (h[i76].re == 0.0) {
        h[i76].re = h[i76].im / s_tmp[i76].im;
        h[i76].im = 0.0;
      } else if (h[i76].im == 0.0) {
        h[i76].re = 0.0;
        h[i76].im = -(h_re / s_tmp[i76].im);
      } else {
        h[i76].re = h[i76].im / s_tmp[i76].im;
        h[i76].im = -(h_re / s_tmp[i76].im);
      }
    } else {
      brm = fabs(s_tmp[i76].re);
      re = fabs(s_tmp[i76].im);
      if (brm > re) {
        re = s_tmp[i76].im / s_tmp[i76].re;
        re_tmp = s_tmp[i76].re + re * s_tmp[i76].im;
        h[i76].re = (h[i76].re + re * h[i76].im) / re_tmp;
        h[i76].im = (h[i76].im - re * h_re) / re_tmp;
      } else if (re == brm) {
        if (s_tmp[i76].re > 0.0) {
          re = 0.5;
        } else {
          re = -0.5;
        }

        if (s_tmp[i76].im > 0.0) {
          re_tmp = 0.5;
        } else {
          re_tmp = -0.5;
        }

        h[i76].re = (h[i76].re * re + h[i76].im * re_tmp) / brm;
        h[i76].im = (h[i76].im * re - h_re * re_tmp) / brm;
      } else {
        re = s_tmp[i76].re / s_tmp[i76].im;
        re_tmp = s_tmp[i76].im + re * s_tmp[i76].re;
        h[i76].re = (re * h[i76].re + h[i76].im) / re_tmp;
        h[i76].im = (re * h[i76].im - h_re) / re_tmp;
      }
    }
  }
}

/*
 * FREQZ_CG Frequency response of digital filter with codegen support
 *
 *  This function is based on 'freqz' by The MathWorks Inc.
 * Arguments    : const double b[29]
 *                const double w[2048]
 *                double Fs
 *                creal_T hh[2048]
 * Return Type  : void
 */
static void e_freqz_cg(const double b[29], const double w[2048], double Fs,
  creal_T hh[2048])
{
  static struct_T options;
  int i16;
  static const char cv20[8] = { 'o', 'n', 'e', 's', 'i', 'd', 'e', 'd' };

  double b_w[2048];
  static const char cv21[7] = { 'o', 'm', 'i', 't', 't', 'e', 'd' };

  /*  Cast to enforce precision rules */
  options.Fs = Fs;
  memcpy(&options.nfft[0], &w[0], sizeof(double) << 11);
  memcpy(&options.w[0], &w[0], sizeof(double) << 11);

  /*  Remaining are default or for advanced use */
  options.fvflag = 1.0;
  for (i16 = 0; i16 < 8; i16++) {
    options.range[i16] = cv20[i16];
  }

  options.centerdc = 0.0;
  for (i16 = 0; i16 < 7; i16++) {
    options.configlevel[i16] = cv21[i16];
  }

  e_firfreqz(b, &options, hh, b_w);

  /*  Generate the default structure to pass to freqzplot */
  /*  If rad/sample, Fs is empty */
}

/*
 * Arguments    : const double p[13]
 *                const creal_T x[2048]
 *                creal_T y[2048]
 * Return Type  : void
 */
static void e_polyval(const double p[13], const creal_T x[2048], creal_T y[2048])
{
  int i17;
  int k;
  double p_re;
  double x_im;
  for (i17 = 0; i17 < 2048; i17++) {
    y[i17].re = p[0];
    y[i17].im = 0.0;
  }

  for (k = 0; k < 12; k++) {
    p_re = p[k + 1];
    for (i17 = 0; i17 < 2048; i17++) {
      x_im = x[i17].re * y[i17].im + x[i17].im * y[i17].re;
      y[i17].re = (x[i17].re * y[i17].re - x[i17].im * y[i17].im) + p_re;
      y[i17].im = x_im;
    }
  }
}

/*
 * Arguments    : const double o[7]
 *                double u[13]
 * Return Type  : void
 */
static void e_us(const double o[7], double u[13])
{
  double b_u[14];
  int ix;
  int iy;
  int k;
  memset(&b_u[0], 0, 14U * sizeof(double));
  ix = 0;
  iy = 0;
  for (k = 0; k < 7; k++) {
    b_u[iy] = o[ix];
    ix++;
    iy += 2;
  }

  memcpy(&u[0], &b_u[0], 13U * sizeof(double));
}

/*
 * Arguments    : emxArray_creal_T *h
 * Return Type  : int
 */
static int eml_zlahqr(emxArray_creal_T *h)
{
  int info;
  int n;
  double itmax;
  int ldh;
  int knt;
  int j;
  int i;
  double SMLNUM;
  double aa;
  bool exitg1;
  double tst;
  double br;
  int L;
  creal_T sc;
  bool goto140;
  int its;
  bool exitg2;
  int k;
  bool exitg3;
  creal_T b_sc;
  double htmp1;
  double t_re;
  double ab;
  bool goto70;
  double ba;
  int m;
  double u_re;
  double x_im;
  double u_im;
  double s;
  int b_k;
  creal_T v[2];
  n = h->size[0];
  itmax = 30.0 * fmax(10.0, h->size[0]);
  ldh = h->size[0];
  info = 0;
  if ((h->size[0] != 0) && (1 != h->size[0])) {
    knt = h->size[0];
    for (j = 0; j <= knt - 4; j++) {
      h->data[(j + h->size[0] * j) + 2].re = 0.0;
      h->data[(j + h->size[0] * j) + 2].im = 0.0;
      h->data[(j + h->size[0] * j) + 3].re = 0.0;
      h->data[(j + h->size[0] * j) + 3].im = 0.0;
    }

    if (1 <= n - 2) {
      h->data[(n + h->size[0] * (n - 3)) - 1].re = 0.0;
      h->data[(n + h->size[0] * (n - 3)) - 1].im = 0.0;
    }

    for (i = 2; i <= n; i++) {
      if (h->data[(i + h->size[0] * (i - 2)) - 1].im != 0.0) {
        aa = h->data[(i + h->size[0] * (i - 2)) - 1].re;
        tst = h->data[(i + h->size[0] * (i - 2)) - 1].im;
        br = fabs(h->data[(i + h->size[0] * (i - 2)) - 1].re) + fabs(h->data[(i
          + h->size[0] * (i - 2)) - 1].im);
        if (tst == 0.0) {
          sc.re = aa / br;
          sc.im = 0.0;
        } else if (aa == 0.0) {
          sc.re = 0.0;
          sc.im = tst / br;
        } else {
          sc.re = aa / br;
          sc.im = tst / br;
        }

        br = rt_hypotd_snf(sc.re, sc.im);
        if (-sc.im == 0.0) {
          sc.re /= br;
          sc.im = 0.0;
        } else if (sc.re == 0.0) {
          sc.re = 0.0;
          sc.im = -sc.im / br;
        } else {
          sc.re /= br;
          sc.im = -sc.im / br;
        }

        tst = h->data[(i + h->size[0] * (i - 2)) - 1].re;
        aa = h->data[(i + h->size[0] * (i - 2)) - 1].im;
        h->data[(i + h->size[0] * (i - 2)) - 1].re = rt_hypotd_snf(tst, aa);
        h->data[(i + h->size[0] * (i - 2)) - 1].im = 0.0;
        b_xscal((n - i) + 1, sc, h, i + (i - 1) * ldh, ldh);
        b_sc.re = sc.re;
        b_sc.im = -sc.im;
        knt = i + 1;
        if (n < knt) {
          knt = n;
        }

        xscal(knt, b_sc, h, 1 + (i - 1) * ldh);
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
      while ((!exitg2) && (its <= (int)itmax)) {
        k = i;
        exitg3 = false;
        while ((!exitg3) && ((k + 1 > L + 2) && (!(fabs(h->data[k + h->size[0] *
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
            br = h->data[(k + h->size[0] * (k - 1)) - 1].re - h->data[k +
              h->size[0] * k].re;
            x_im = h->data[(k + h->size[0] * (k - 1)) - 1].im - h->data[k +
              h->size[0] * k].im;
            tst = fabs(br) + fabs(x_im);
            if (htmp1 > tst) {
              aa = htmp1;
              htmp1 = tst;
            } else {
              aa = tst;
            }

            s = aa + ab;
            if (ba * (ab / s) <= fmax(SMLNUM, 2.2204460492503131E-16 * (htmp1 *
                  (aa / s)))) {
              exitg3 = true;
            } else {
              k--;
            }
          } else {
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
            ab = h->data[k + h->size[0] * k].im;
          } else if (its == 20) {
            t_re = 0.75 * fabs(h->data[i + h->size[0] * (i - 1)].re) + h->data[i
              + h->size[0] * i].re;
            ab = h->data[i + h->size[0] * i].im;
          } else {
            t_re = h->data[i + h->size[0] * i].re;
            ab = h->data[i + h->size[0] * i].im;
            b_sc = h->data[(i + h->size[0] * i) - 1];
            c_sqrt(&b_sc);
            sc = h->data[i + h->size[0] * (i - 1)];
            c_sqrt(&sc);
            u_re = b_sc.re * sc.re - b_sc.im * sc.im;
            u_im = b_sc.re * sc.im + b_sc.im * sc.re;
            s = fabs(u_re) + fabs(u_im);
            if (s != 0.0) {
              tst = h->data[(i + h->size[0] * (i - 1)) - 1].re - h->data[i +
                h->size[0] * i].re;
              aa = h->data[(i + h->size[0] * (i - 1)) - 1].im - h->data[i +
                h->size[0] * i].im;
              br = 0.5 * tst;
              x_im = 0.5 * aa;
              tst = fabs(br) + fabs(x_im);
              s = fmax(s, tst);
              if (x_im == 0.0) {
                t_re = br / s;
                ab = 0.0;
              } else if (br == 0.0) {
                t_re = 0.0;
                ab = x_im / s;
              } else {
                t_re = br / s;
                ab = x_im / s;
              }

              htmp1 = t_re;
              t_re = t_re * t_re - ab * ab;
              ab = htmp1 * ab + ab * htmp1;
              if (u_im == 0.0) {
                sc.re = u_re / s;
                sc.im = 0.0;
              } else if (u_re == 0.0) {
                sc.re = 0.0;
                sc.im = u_im / s;
              } else {
                sc.re = u_re / s;
                sc.im = u_im / s;
              }

              b_sc.re = t_re + (sc.re * sc.re - sc.im * sc.im);
              b_sc.im = ab + (sc.re * sc.im + sc.im * sc.re);
              c_sqrt(&b_sc);
              sc.re = s * b_sc.re;
              sc.im = s * b_sc.im;
              if (tst > 0.0) {
                if (x_im == 0.0) {
                  t_re = br / tst;
                  ab = 0.0;
                } else if (br == 0.0) {
                  t_re = 0.0;
                  ab = x_im / tst;
                } else {
                  t_re = br / tst;
                  ab = x_im / tst;
                }

                if (t_re * sc.re + ab * sc.im < 0.0) {
                  sc.re = -sc.re;
                  sc.im = -sc.im;
                }
              }

              br += sc.re;
              htmp1 = x_im + sc.im;
              if (htmp1 == 0.0) {
                if (u_im == 0.0) {
                  ba = u_re / br;
                  tst = 0.0;
                } else if (u_re == 0.0) {
                  ba = 0.0;
                  tst = u_im / br;
                } else {
                  ba = u_re / br;
                  tst = u_im / br;
                }
              } else if (br == 0.0) {
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
                ab = fabs(br);
                tst = fabs(htmp1);
                if (ab > tst) {
                  s = htmp1 / br;
                  tst = br + s * htmp1;
                  ba = (u_re + s * u_im) / tst;
                  tst = (u_im - s * u_re) / tst;
                } else if (tst == ab) {
                  if (br > 0.0) {
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
                  s = br / htmp1;
                  tst = htmp1 + s * br;
                  ba = (s * u_re + u_im) / tst;
                  tst = (s * u_im - u_re) / tst;
                }
              }

              t_re = h->data[i + h->size[0] * i].re - (u_re * ba - u_im * tst);
              ab = h->data[i + h->size[0] * i].im - (u_re * tst + u_im * ba);
            }
          }

          goto70 = false;
          m = i;
          exitg3 = false;
          while ((!exitg3) && (m > k + 1)) {
            sc.re = h->data[(m + h->size[0] * (m - 1)) - 1].re - t_re;
            sc.im = h->data[(m + h->size[0] * (m - 1)) - 1].im - ab;
            tst = h->data[m + h->size[0] * (m - 1)].re;
            s = (fabs(sc.re) + fabs(sc.im)) + fabs(tst);
            if (sc.im == 0.0) {
              sc.re /= s;
              sc.im = 0.0;
            } else if (sc.re == 0.0) {
              sc.re = 0.0;
              sc.im /= s;
            } else {
              sc.re /= s;
              sc.im /= s;
            }

            tst /= s;
            v[0] = sc;
            v[1].re = tst;
            v[1].im = 0.0;
            if (fabs(h->data[(m + h->size[0] * (m - 2)) - 1].re) * fabs(tst) <=
                2.2204460492503131E-16 * ((fabs(sc.re) + fabs(sc.im)) * ((fabs
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
            sc.re = h->data[k + h->size[0] * k].re - t_re;
            sc.im = h->data[k + h->size[0] * k].im - ab;
            tst = h->data[(k + h->size[0] * k) + 1].re;
            s = (fabs(sc.re) + fabs(sc.im)) + fabs(tst);
            if (sc.im == 0.0) {
              v[0].re = sc.re / s;
              v[0].im = 0.0;
            } else if (sc.re == 0.0) {
              v[0].re = 0.0;
              v[0].im = sc.im / s;
            } else {
              v[0].re = sc.re / s;
              v[0].im = sc.im / s;
            }

            tst /= s;
            v[1].re = tst;
            v[1].im = 0.0;
          }

          for (b_k = m; b_k <= i; b_k++) {
            if (b_k > m) {
              v[0] = h->data[(b_k + h->size[0] * (b_k - 2)) - 1];
              v[1] = h->data[b_k + h->size[0] * (b_k - 2)];
            }

            sc = v[0];
            br = v[1].re;
            x_im = v[1].im;
            t_re = 0.0;
            ab = 0.0;
            tst = rt_hypotd_snf(v[1].re, v[1].im);
            if ((tst != 0.0) || (v[0].im != 0.0)) {
              htmp1 = xdlapy3(v[0].re, v[0].im, tst);
              if (v[0].re >= 0.0) {
                htmp1 = -htmp1;
              }

              if (fabs(htmp1) < 1.0020841800044864E-292) {
                knt = -1;
                do {
                  knt++;
                  br *= 9.9792015476736E+291;
                  x_im *= 9.9792015476736E+291;
                  htmp1 *= 9.9792015476736E+291;
                  sc.re *= 9.9792015476736E+291;
                  sc.im *= 9.9792015476736E+291;
                } while (!(fabs(htmp1) >= 1.0020841800044864E-292));

                htmp1 = xdlapy3(sc.re, sc.im, rt_hypotd_snf(br, x_im));
                if (sc.re >= 0.0) {
                  htmp1 = -htmp1;
                }

                aa = htmp1 - sc.re;
                if (0.0 - sc.im == 0.0) {
                  t_re = aa / htmp1;
                  ab = 0.0;
                } else if (aa == 0.0) {
                  t_re = 0.0;
                  ab = (0.0 - sc.im) / htmp1;
                } else {
                  t_re = aa / htmp1;
                  ab = (0.0 - sc.im) / htmp1;
                }

                b_sc.re = sc.re - htmp1;
                b_sc.im = sc.im;
                sc = recip(b_sc);
                tst = br;
                br = sc.re * br - sc.im * x_im;
                x_im = sc.re * x_im + sc.im * tst;
                for (j = 0; j <= knt; j++) {
                  htmp1 *= 1.0020841800044864E-292;
                }

                sc.re = htmp1;
                sc.im = 0.0;
              } else {
                aa = htmp1 - v[0].re;
                if (0.0 - v[0].im == 0.0) {
                  t_re = aa / htmp1;
                  ab = 0.0;
                } else if (aa == 0.0) {
                  t_re = 0.0;
                  ab = (0.0 - v[0].im) / htmp1;
                } else {
                  t_re = aa / htmp1;
                  ab = (0.0 - v[0].im) / htmp1;
                }

                b_sc.re = v[0].re - htmp1;
                b_sc.im = v[0].im;
                b_sc = recip(b_sc);
                br = b_sc.re * v[1].re - b_sc.im * v[1].im;
                x_im = b_sc.re * v[1].im + b_sc.im * v[1].re;
                sc.re = htmp1;
                sc.im = 0.0;
              }
            }

            v[0] = sc;
            v[1].re = br;
            v[1].im = x_im;
            if (b_k > m) {
              h->data[(b_k + h->size[0] * (b_k - 2)) - 1] = sc;
              h->data[b_k + h->size[0] * (b_k - 2)].re = 0.0;
              h->data[b_k + h->size[0] * (b_k - 2)].im = 0.0;
            }

            htmp1 = t_re * br - ab * x_im;
            for (j = b_k; j <= n; j++) {
              tst = t_re * h->data[(b_k + h->size[0] * (j - 1)) - 1].re - -ab *
                h->data[(b_k + h->size[0] * (j - 1)) - 1].im;
              aa = t_re * h->data[(b_k + h->size[0] * (j - 1)) - 1].im + -ab *
                h->data[(b_k + h->size[0] * (j - 1)) - 1].re;
              sc.re = tst + htmp1 * h->data[b_k + h->size[0] * (j - 1)].re;
              sc.im = aa + htmp1 * h->data[b_k + h->size[0] * (j - 1)].im;
              h->data[(b_k + h->size[0] * (j - 1)) - 1].re -= sc.re;
              h->data[(b_k + h->size[0] * (j - 1)) - 1].im -= sc.im;
              h->data[b_k + h->size[0] * (j - 1)].re -= sc.re * br - sc.im *
                x_im;
              h->data[b_k + h->size[0] * (j - 1)].im -= sc.re * x_im + sc.im *
                br;
            }

            if (b_k + 2 < i + 1) {
              knt = b_k + 1;
            } else {
              knt = i;
            }

            for (j = 0; j <= knt; j++) {
              tst = t_re * h->data[j + h->size[0] * (b_k - 1)].re - ab * h->
                data[j + h->size[0] * (b_k - 1)].im;
              aa = t_re * h->data[j + h->size[0] * (b_k - 1)].im + ab * h->
                data[j + h->size[0] * (b_k - 1)].re;
              sc.re = tst + htmp1 * h->data[j + h->size[0] * b_k].re;
              sc.im = aa + htmp1 * h->data[j + h->size[0] * b_k].im;
              h->data[j + h->size[0] * (b_k - 1)].re -= sc.re;
              h->data[j + h->size[0] * (b_k - 1)].im -= sc.im;
              h->data[j + h->size[0] * b_k].re -= sc.re * br - sc.im * -x_im;
              h->data[j + h->size[0] * b_k].im -= sc.re * -x_im + sc.im * br;
            }

            if ((b_k == m) && (m > k + 1)) {
              br = rt_hypotd_snf(1.0 - t_re, 0.0 - ab);
              if (0.0 - ab == 0.0) {
                sc.re = (1.0 - t_re) / br;
                sc.im = 0.0;
              } else if (1.0 - t_re == 0.0) {
                sc.re = 0.0;
                sc.im = (0.0 - ab) / br;
              } else {
                sc.re = (1.0 - t_re) / br;
                sc.im = (0.0 - ab) / br;
              }

              tst = h->data[m + h->size[0] * (m - 1)].re;
              aa = h->data[m + h->size[0] * (m - 1)].im;
              h->data[m + h->size[0] * (m - 1)].re = tst * sc.re - aa * -sc.im;
              h->data[m + h->size[0] * (m - 1)].im = tst * -sc.im + aa * sc.re;
              if (m + 2 <= i + 1) {
                tst = h->data[(m + h->size[0] * m) + 1].re;
                aa = h->data[(m + h->size[0] * m) + 1].im;
                h->data[(m + h->size[0] * m) + 1].re = tst * sc.re - aa * sc.im;
                h->data[(m + h->size[0] * m) + 1].im = tst * sc.im + aa * sc.re;
              }

              for (j = m; j <= i + 1; j++) {
                if (j != m + 1) {
                  if (n > j) {
                    b_xscal(n - j, sc, h, j + j * ldh, ldh);
                  }

                  b_sc.re = sc.re;
                  b_sc.im = -sc.im;
                  xscal(j - 1, b_sc, h, 1 + (j - 1) * ldh);
                }
              }
            }
          }

          sc = h->data[i + h->size[0] * (i - 1)];
          if (h->data[i + h->size[0] * (i - 1)].im != 0.0) {
            tst = rt_hypotd_snf(h->data[i + h->size[0] * (i - 1)].re, h->data[i
                                + h->size[0] * (i - 1)].im);
            h->data[i + h->size[0] * (i - 1)].re = tst;
            h->data[i + h->size[0] * (i - 1)].im = 0.0;
            if (sc.im == 0.0) {
              sc.re /= tst;
              sc.im = 0.0;
            } else if (sc.re == 0.0) {
              sc.re = 0.0;
              sc.im /= tst;
            } else {
              sc.re /= tst;
              sc.im /= tst;
            }

            if (n > i + 1) {
              b_sc.re = sc.re;
              b_sc.im = -sc.im;
              b_xscal((n - i) - 1, b_sc, h, (i + (i + 1) * ldh) + 1, ldh);
            }

            xscal(i, sc, h, 1 + i * ldh);
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
 * Make b a row
 * Arguments    : const double b[13]
 *                const struct_T *options
 *                creal_T h[2048]
 *                double w[2048]
 * Return Type  : void
 */
static void f_firfreqz(const double b[13], const struct_T *options, creal_T h
  [2048], double w[2048])
{
  int i77;
  creal_T dcv6[2048];
  double re_tmp;
  double re;
  creal_T s_tmp[2048];
  double h_re;
  double brm;

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
  for (i77 = 0; i77 < 2048; i77++) {
    w[i77] = options->w[i77];
    re_tmp = 6.2831853071795862 * options->w[i77] / options->Fs;
    re = re_tmp * 0.0;
    s_tmp[i77].re = re;
    s_tmp[i77].im = re_tmp;
    dcv6[i77].re = re;
    dcv6[i77].im = re_tmp;
  }

  b_exp(dcv6);
  e_polyval(b, dcv6, h);
  for (i77 = 0; i77 < 2048; i77++) {
    s_tmp[i77].re *= 12.0;
    s_tmp[i77].im *= 12.0;
  }

  b_exp(s_tmp);
  for (i77 = 0; i77 < 2048; i77++) {
    h_re = h[i77].re;
    if (s_tmp[i77].im == 0.0) {
      if (h[i77].im == 0.0) {
        h[i77].re /= s_tmp[i77].re;
        h[i77].im = 0.0;
      } else if (h[i77].re == 0.0) {
        h[i77].re = 0.0;
        h[i77].im /= s_tmp[i77].re;
      } else {
        h[i77].re /= s_tmp[i77].re;
        h[i77].im /= s_tmp[i77].re;
      }
    } else if (s_tmp[i77].re == 0.0) {
      if (h[i77].re == 0.0) {
        h[i77].re = h[i77].im / s_tmp[i77].im;
        h[i77].im = 0.0;
      } else if (h[i77].im == 0.0) {
        h[i77].re = 0.0;
        h[i77].im = -(h_re / s_tmp[i77].im);
      } else {
        h[i77].re = h[i77].im / s_tmp[i77].im;
        h[i77].im = -(h_re / s_tmp[i77].im);
      }
    } else {
      brm = fabs(s_tmp[i77].re);
      re = fabs(s_tmp[i77].im);
      if (brm > re) {
        re = s_tmp[i77].im / s_tmp[i77].re;
        re_tmp = s_tmp[i77].re + re * s_tmp[i77].im;
        h[i77].re = (h[i77].re + re * h[i77].im) / re_tmp;
        h[i77].im = (h[i77].im - re * h_re) / re_tmp;
      } else if (re == brm) {
        if (s_tmp[i77].re > 0.0) {
          re = 0.5;
        } else {
          re = -0.5;
        }

        if (s_tmp[i77].im > 0.0) {
          re_tmp = 0.5;
        } else {
          re_tmp = -0.5;
        }

        h[i77].re = (h[i77].re * re + h[i77].im * re_tmp) / brm;
        h[i77].im = (h[i77].im * re - h_re * re_tmp) / brm;
      } else {
        re = s_tmp[i77].re / s_tmp[i77].im;
        re_tmp = s_tmp[i77].im + re * s_tmp[i77].re;
        h[i77].re = (re * h[i77].re + h[i77].im) / re_tmp;
        h[i77].im = (re * h[i77].im - h_re) / re_tmp;
      }
    }
  }
}

/*
 * FREQZ_CG Frequency response of digital filter with codegen support
 *
 *  This function is based on 'freqz' by The MathWorks Inc.
 * Arguments    : const double b[13]
 *                const double w[2048]
 *                double Fs
 *                creal_T hh[2048]
 * Return Type  : void
 */
static void f_freqz_cg(const double b[13], const double w[2048], double Fs,
  creal_T hh[2048])
{
  static struct_T options;
  int i18;
  static const char cv22[8] = { 'o', 'n', 'e', 's', 'i', 'd', 'e', 'd' };

  double b_w[2048];
  static const char cv23[7] = { 'o', 'm', 'i', 't', 't', 'e', 'd' };

  /*  Cast to enforce precision rules */
  options.Fs = Fs;
  memcpy(&options.nfft[0], &w[0], sizeof(double) << 11);
  memcpy(&options.w[0], &w[0], sizeof(double) << 11);

  /*  Remaining are default or for advanced use */
  options.fvflag = 1.0;
  for (i18 = 0; i18 < 8; i18++) {
    options.range[i18] = cv22[i18];
  }

  options.centerdc = 0.0;
  for (i18 = 0; i18 < 7; i18++) {
    options.configlevel[i18] = cv23[i18];
  }

  f_firfreqz(b, &options, hh, b_w);

  /*  Generate the default structure to pass to freqzplot */
  /*  If rad/sample, Fs is empty */
}

/*
 * Arguments    : const double p[57]
 *                const creal_T x[2048]
 *                creal_T y[2048]
 * Return Type  : void
 */
static void f_polyval(const double p[57], const creal_T x[2048], creal_T y[2048])
{
  int i19;
  int k;
  double p_re;
  double x_im;
  for (i19 = 0; i19 < 2048; i19++) {
    y[i19].re = p[0];
    y[i19].im = 0.0;
  }

  for (k = 0; k < 56; k++) {
    p_re = p[k + 1];
    for (i19 = 0; i19 < 2048; i19++) {
      x_im = x[i19].re * y[i19].im + x[i19].im * y[i19].re;
      y[i19].re = (x[i19].re * y[i19].re - x[i19].im * y[i19].im) + p_re;
      y[i19].im = x_im;
    }
  }
}

/*
 * Arguments    : const double o[15]
 *                double u[57]
 * Return Type  : void
 */
static void f_us(const double o[15], double u[57])
{
  double b_u[60];
  int ix;
  int iy;
  int k;
  memset(&b_u[0], 0, 60U * sizeof(double));
  ix = 0;
  iy = 0;
  for (k = 0; k < 15; k++) {
    b_u[iy] = o[ix];
    ix++;
    iy += 4;
  }

  memcpy(&u[0], &b_u[0], 57U * sizeof(double));
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
  int i72;
  double re_tmp;
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
  for (i72 = 0; i72 < 2048; i72++) {
    w[i72] = options->w[i72];
    re_tmp = 6.2831853071795862 * options->w[i72] / options->Fs;
    h[i72].re = 0.0 * (re_tmp * 0.0);
    h[i72].im = 0.0 * re_tmp;
  }

  b_exp(h);
  for (i72 = 0; i72 < 2048; i72++) {
    if (h[i72].im == 0.0) {
      h[i72].re = 1.0 / h[i72].re;
      h[i72].im = 0.0;
    } else if (h[i72].re == 0.0) {
      h[i72].re = 0.0;
      h[i72].im = -(1.0 / h[i72].im);
    } else {
      brm = fabs(h[i72].re);
      re_tmp = fabs(h[i72].im);
      if (brm > re_tmp) {
        re_tmp = h[i72].im / h[i72].re;
        d = h[i72].re + re_tmp * h[i72].im;
        h[i72].re = (1.0 + re_tmp * 0.0) / d;
        h[i72].im = (0.0 - re_tmp) / d;
      } else if (re_tmp == brm) {
        if (h[i72].re > 0.0) {
          re_tmp = 0.5;
        } else {
          re_tmp = -0.5;
        }

        if (h[i72].im > 0.0) {
          d = 0.5;
        } else {
          d = -0.5;
        }

        h[i72].re = (re_tmp + 0.0 * d) / brm;
        h[i72].im = (0.0 * re_tmp - d) / brm;
      } else {
        re_tmp = h[i72].re / h[i72].im;
        d = h[i72].im + re_tmp * h[i72].re;
        h[i72].re = re_tmp / d;
        h[i72].im = (re_tmp * 0.0 - 1.0) / d;
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
  emxArray_real_T *r8;
  emxArray_real_T *r9;
  double b_ff[4];
  emxArray_real_T *b_grid;
  int i45;
  int loop_ub;
  emxArray_real_T *b_h;
  double err;
  bool valid;
  int h_idx_0;
  int i46;
  int i47;
  emxArray_real_T *c_h;
  emxInit_real_T(&grid, 2);
  emxInit_real_T(&r8, 2);
  emxInit_real_T(&r9, 2);

  /*  */
  firpmgrid_cg(order + 1.0, ff, grid);

  /*      orgFreqIndx = frequencies; */
  /*   */
  /*      positionsOfNewFreqIndx = zeros(size(grid)); */
  /*      for ind = 1:length(positionsOfNewFreqIndx) */
  /*   */
  /*          [~,indx] = min( abs(orgFreqIndx-grid(ind)) ); */
  /*   */
  /*          positionsOfNewFreqIndx(ind) = indx; */
  /*      end */
  /*   */
  /*      wt = weights(positionsOfNewFreqIndx); */
  /*      des = amplitudes(positionsOfNewFreqIndx); */
  /*  Workaround */
  /* ftype = 2; */
  /* sign_val = 1; */
  /*  Always bandpass designs */
  /*  cast to enforce precision rules */
  /*  Call actual design algorithm */
  interp1(frequencies, amplitudes, grid, r8);
  interp1(frequencies, weights, grid, r9);
  b_ff[0] = ff[0] / 2.0;
  b_ff[1] = ff[1] / 2.0;
  b_ff[2] = ff[2] / 2.0;
  b_ff[3] = ff[3] / 2.0;
  emxInit_real_T(&b_grid, 2);
  i45 = b_grid->size[0] * b_grid->size[1];
  b_grid->size[0] = 1;
  b_grid->size[1] = grid->size[1];
  emxEnsureCapacity_real_T(b_grid, i45);
  loop_ub = grid->size[0] * grid->size[1];
  for (i45 = 0; i45 < loop_ub; i45++) {
    b_grid->data[i45] = grid->data[i45] / 2.0;
  }

  emxFree_real_T(&grid);
  emxInit_real_T(&b_h, 2);
  remezm(order + 1.0, b_ff, b_grid, r8, r9, b_h, &err, &valid);
  h_idx_0 = b_h->size[0] * b_h->size[1];
  i45 = h->size[0] * h->size[1];
  h->size[0] = 1;
  h->size[1] = h_idx_0;
  emxEnsureCapacity_real_T(h, i45);
  emxFree_real_T(&b_grid);
  emxFree_real_T(&r9);
  emxFree_real_T(&r8);
  for (i45 = 0; i45 < h_idx_0; i45++) {
    h->data[i45] = b_h->data[i45];
  }

  emxFree_real_T(&b_h);

  /*  make it a row */
  err = (double)h->size[1] - rt_remd_snf(order + 1.0, 2.0);
  if (1.0 > err) {
    i45 = 1;
    i46 = 1;
    i47 = 0;
  } else {
    i45 = (int)err;
    i46 = -1;
    i47 = 1;
  }

  emxInit_real_T(&c_h, 2);
  h_idx_0 = c_h->size[0] * c_h->size[1];
  c_h->size[0] = 1;
  loop_ub = div_s32_floor(i47 - i45, i46);
  c_h->size[1] = (h->size[1] + loop_ub) + 1;
  emxEnsureCapacity_real_T(c_h, h_idx_0);
  h_idx_0 = h->size[1];
  for (i47 = 0; i47 < h_idx_0; i47++) {
    c_h->data[i47] = h->data[i47];
  }

  for (i47 = 0; i47 <= loop_ub; i47++) {
    c_h->data[i47 + h->size[1]] = h->data[(i45 + i46 * i47) - 1];
  }

  i45 = h->size[0] * h->size[1];
  h->size[0] = 1;
  h->size[1] = c_h->size[1];
  emxEnsureCapacity_real_T(h, i45);
  loop_ub = c_h->size[0] * c_h->size[1];
  for (i45 = 0; i45 < loop_ub; i45++) {
    h->data[i45] = c_h->data[i45];
  }

  if (1 > h->size[1]) {
    i45 = 1;
    i46 = 1;
    i47 = 0;
  } else {
    i45 = h->size[1];
    i46 = -1;
    i47 = 1;
  }

  h_idx_0 = c_h->size[0] * c_h->size[1];
  c_h->size[0] = 1;
  loop_ub = div_s32_floor(i47 - i45, i46);
  c_h->size[1] = loop_ub + 1;
  emxEnsureCapacity_real_T(c_h, h_idx_0);
  for (i47 = 0; i47 <= loop_ub; i47++) {
    c_h->data[i47] = h->data[(i45 + i46 * i47) - 1];
  }

  i45 = h->size[0] * h->size[1];
  h->size[0] = 1;
  h->size[1] = c_h->size[1];
  emxEnsureCapacity_real_T(h, i45);
  loop_ub = c_h->size[0] * c_h->size[1];
  for (i45 = 0; i45 < loop_ub; i45++) {
    h->data[i45] = c_h->data[i45];
  }

  emxFree_real_T(&c_h);
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
  double delf;
  double j;
  double l;
  double gridSize;
  emxArray_real_T *newgrid;
  emxArray_int32_T *r10;
  double a_tmp;
  double a;
  double ngrid;
  double b_tmp;
  int i48;
  int nm1d2;
  double ndbl;
  double apnd;
  double cdiff;
  int k;
  int n;

  /* -------------------------------------------------------------------------- */
  memset(&grid[0], 0, 100000U * sizeof(double));

  /*  Make large initial memory */
  /*     Generate frequency grid */
  grid[0] = ff[0];
  delf = 1.0 / (16.0 * trunc(nfilt / 2.0));

  /*  If value at frequency 0 is constrained, make sure first grid point */
  /*  is not too small: */
  j = 1.0;
  l = 1.0;
  gridSize = 1.0;
  emxInit_real_T(&newgrid, 2);
  emxInit_int32_T(&r10, 2);
  while (l + 1.0 <= 4.0) {
    a_tmp = grid[(int)j - 1];
    a = a_tmp + delf;
    b_tmp = ff[(int)(l + 1.0) - 1];
    ngrid = b_tmp + delf;
    if (rtIsNaN(a) || rtIsNaN(delf) || rtIsNaN(ngrid)) {
      i48 = newgrid->size[0] * newgrid->size[1];
      newgrid->size[0] = 1;
      newgrid->size[1] = 1;
      emxEnsureCapacity_real_T(newgrid, i48);
      newgrid->data[0] = rtNaN;
    } else if ((delf == 0.0) || ((a < ngrid) && (delf < 0.0)) || ((ngrid < a) &&
                (delf > 0.0))) {
      newgrid->size[0] = 1;
      newgrid->size[1] = 0;
    } else if ((rtIsInf(a) || rtIsInf(ngrid)) && (rtIsInf(delf) || (a == ngrid)))
    {
      i48 = newgrid->size[0] * newgrid->size[1];
      newgrid->size[0] = 1;
      newgrid->size[1] = 1;
      emxEnsureCapacity_real_T(newgrid, i48);
      newgrid->data[0] = rtNaN;
    } else if (rtIsInf(delf)) {
      i48 = newgrid->size[0] * newgrid->size[1];
      newgrid->size[0] = 1;
      newgrid->size[1] = 1;
      emxEnsureCapacity_real_T(newgrid, i48);
      newgrid->data[0] = a;
    } else if ((floor(a) == a) && (floor(delf) == delf)) {
      i48 = newgrid->size[0] * newgrid->size[1];
      newgrid->size[0] = 1;
      nm1d2 = (int)floor((ngrid - a) / delf);
      newgrid->size[1] = nm1d2 + 1;
      emxEnsureCapacity_real_T(newgrid, i48);
      for (i48 = 0; i48 <= nm1d2; i48++) {
        newgrid->data[i48] = a + delf * (double)i48;
      }
    } else {
      ndbl = floor((ngrid - a) / delf + 0.5);
      apnd = a + ndbl * delf;
      if (delf > 0.0) {
        cdiff = apnd - ngrid;
      } else {
        cdiff = ngrid - apnd;
      }

      if (fabs(cdiff) < 4.4408920985006262E-16 * fmax(fabs(a), fabs(ngrid))) {
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

      i48 = newgrid->size[0] * newgrid->size[1];
      newgrid->size[0] = 1;
      newgrid->size[1] = n;
      emxEnsureCapacity_real_T(newgrid, i48);
      if (n > 0) {
        newgrid->data[0] = a;
        if (n > 1) {
          newgrid->data[n - 1] = apnd;
          nm1d2 = (n - 1) / 2;
          for (k = 0; k <= nm1d2 - 2; k++) {
            ngrid = (1.0 + (double)k) * delf;
            newgrid->data[1 + k] = a + ngrid;
            newgrid->data[(n - k) - 2] = apnd - ngrid;
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
      j = ((ff[(int)(l + 1.0) - 1] + delf) - (grid[(int)j - 1] + delf)) / 10.0;
      a = a_tmp + j;
      ngrid = b_tmp + j;
      if (rtIsNaN(a) || rtIsNaN(j) || rtIsNaN(ngrid)) {
        i48 = newgrid->size[0] * newgrid->size[1];
        newgrid->size[0] = 1;
        newgrid->size[1] = 1;
        emxEnsureCapacity_real_T(newgrid, i48);
        newgrid->data[0] = rtNaN;
      } else if ((j == 0.0) || ((a < ngrid) && (j < 0.0)) || ((ngrid < a) && (j >
        0.0))) {
        newgrid->size[0] = 1;
        newgrid->size[1] = 0;
      } else if ((rtIsInf(a) || rtIsInf(ngrid)) && (rtIsInf(j) || (a == ngrid)))
      {
        i48 = newgrid->size[0] * newgrid->size[1];
        newgrid->size[0] = 1;
        newgrid->size[1] = 1;
        emxEnsureCapacity_real_T(newgrid, i48);
        newgrid->data[0] = rtNaN;
      } else if (rtIsInf(j)) {
        i48 = newgrid->size[0] * newgrid->size[1];
        newgrid->size[0] = 1;
        newgrid->size[1] = 1;
        emxEnsureCapacity_real_T(newgrid, i48);
        newgrid->data[0] = a;
      } else if ((floor(a) == a) && (floor(j) == j)) {
        i48 = newgrid->size[0] * newgrid->size[1];
        newgrid->size[0] = 1;
        nm1d2 = (int)floor((ngrid - a) / j);
        newgrid->size[1] = nm1d2 + 1;
        emxEnsureCapacity_real_T(newgrid, i48);
        for (i48 = 0; i48 <= nm1d2; i48++) {
          newgrid->data[i48] = a + j * (double)i48;
        }
      } else {
        ndbl = floor((ngrid - a) / j + 0.5);
        apnd = a + ndbl * j;
        if (j > 0.0) {
          cdiff = apnd - ngrid;
        } else {
          cdiff = ngrid - apnd;
        }

        if (fabs(cdiff) < 4.4408920985006262E-16 * fmax(fabs(a), fabs(ngrid))) {
          ndbl++;
          apnd = ngrid;
        } else if (cdiff > 0.0) {
          apnd = a + (ndbl - 1.0) * j;
        } else {
          ndbl++;
        }

        if (ndbl >= 0.0) {
          n = (int)ndbl;
        } else {
          n = 0;
        }

        i48 = newgrid->size[0] * newgrid->size[1];
        newgrid->size[0] = 1;
        newgrid->size[1] = n;
        emxEnsureCapacity_real_T(newgrid, i48);
        if (n > 0) {
          newgrid->data[0] = a;
          if (n > 1) {
            newgrid->data[n - 1] = apnd;
            nm1d2 = (n - 1) / 2;
            for (k = 0; k <= nm1d2 - 2; k++) {
              ngrid = (1.0 + (double)k) * j;
              newgrid->data[1 + k] = a + ngrid;
              newgrid->data[(n - k) - 2] = apnd - ngrid;
            }

            if (nm1d2 << 1 == n - 1) {
              newgrid->data[nm1d2] = (a + apnd) / 2.0;
            } else {
              ngrid = (double)nm1d2 * j;
              newgrid->data[nm1d2] = a + ngrid;
              newgrid->data[nm1d2 + 1] = apnd - ngrid;
            }
          }
        }
      }
    }

    /* grid = [grid newgrid]; */
    i48 = newgrid->size[1];
    k = r10->size[0] * r10->size[1];
    r10->size[0] = 1;
    nm1d2 = (int)((double)i48 - 1.0);
    r10->size[1] = nm1d2 + 1;
    emxEnsureCapacity_int32_T(r10, k);
    for (i48 = 0; i48 <= nm1d2; i48++) {
      r10->data[i48] = (int)(gridSize + (1.0 + (double)i48));
    }

    nm1d2 = newgrid->size[0] * newgrid->size[1];
    for (i48 = 0; i48 < nm1d2; i48++) {
      grid[r10->data[i48] - 1] = newgrid->data[i48];
    }

    gridSize += (double)newgrid->size[1];

    /* jend = length(grid); */
    /* length(grid); */
    if (gridSize > 1.0) {
      grid[(int)(gridSize - 1.0) - 1] = b_tmp;
      j = gridSize;
    } else {
      j = 2.0;
    }

    l += 2.0;
    if (l + 1.0 <= 4.0) {
      grid[(int)j - 1] = ff[(int)l - 1];
    }
  }

  emxFree_int32_T(&r10);
  emxFree_real_T(&newgrid);
  ngrid = j - 1.0;

  /*  If value at frequency 1 is constrained, remove that grid point: */
  i48 = (int)(j - 1.0) - 1;
  if (grid[i48] > 1.0 - delf) {
    if (ff[2] < 1.0 - delf) {
      ngrid = (j - 1.0) - 1.0;
    } else {
      grid[i48] = ff[2];
    }
  }

  if (1.0 > ngrid) {
    nm1d2 = 0;
  } else {
    nm1d2 = (int)ngrid;
  }

  i48 = gridactual->size[0] * gridactual->size[1];
  gridactual->size[0] = 1;
  gridactual->size[1] = nm1d2;
  emxEnsureCapacity_real_T(gridactual, i48);
  for (i48 = 0; i48 < nm1d2; i48++) {
    gridactual->data[i48] = grid[i48];
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
  int i27;
  bool b0;
  static creal_T s[2048];
  creal_T y[2048];
  int k;
  double h_re;
  double a_re;
  double a_im;
  int i28;
  double s_im;
  emxInit_creal_T(&b_a, 2);
  removeTrailingZero(b_data, b_size, a, b_b_data, b_b_size, b_a);
  for (i27 = 0; i27 < 2048; i27++) {
    s[i27].re = w[i27] * 0.0;
    s[i27].im = w[i27];
  }

  b0 = (b_a->size[1] == 0);
  if (!b0) {
    for (i27 = 0; i27 < 2048; i27++) {
      y[i27] = b_a->data[0];
    }

    i27 = b_a->size[1];
    for (k = 0; k <= i27 - 2; k++) {
      a_re = b_a->data[k + 1].re;
      a_im = b_a->data[k + 1].im;
      for (i28 = 0; i28 < 2048; i28++) {
        s_im = s[i28].re * y[i28].im + s[i28].im * y[i28].re;
        y[i28].re = (s[i28].re * y[i28].re - s[i28].im * y[i28].im) + a_re;
        y[i28].im = s_im + a_im;
      }
    }
  }

  emxFree_creal_T(&b_a);
  c_polyval(b_b_data, b_b_size, s, h);
  for (i27 = 0; i27 < 2048; i27++) {
    h_re = h[i27].re;
    if (y[i27].im == 0.0) {
      if (h[i27].im == 0.0) {
        h[i27].re /= y[i27].re;
        h[i27].im = 0.0;
      } else if (h[i27].re == 0.0) {
        h[i27].re = 0.0;
        h[i27].im /= y[i27].re;
      } else {
        h[i27].re /= y[i27].re;
        h[i27].im /= y[i27].re;
      }
    } else if (y[i27].re == 0.0) {
      if (h[i27].re == 0.0) {
        h[i27].re = h[i27].im / y[i27].im;
        h[i27].im = 0.0;
      } else if (h[i27].im == 0.0) {
        h[i27].re = 0.0;
        h[i27].im = -(h_re / y[i27].im);
      } else {
        h[i27].re = h[i27].im / y[i27].im;
        h[i27].im = -(h_re / y[i27].im);
      }
    } else {
      s_im = fabs(y[i27].re);
      a_re = fabs(y[i27].im);
      if (s_im > a_re) {
        a_re = y[i27].im / y[i27].re;
        a_im = y[i27].re + a_re * y[i27].im;
        h[i27].re = (h[i27].re + a_re * h[i27].im) / a_im;
        h[i27].im = (h[i27].im - a_re * h_re) / a_im;
      } else if (a_re == s_im) {
        if (y[i27].re > 0.0) {
          a_re = 0.5;
        } else {
          a_re = -0.5;
        }

        if (y[i27].im > 0.0) {
          a_im = 0.5;
        } else {
          a_im = -0.5;
        }

        h[i27].re = (h[i27].re * a_re + h[i27].im * a_im) / s_im;
        h[i27].im = (h[i27].im * a_re - h_re * a_im) / s_im;
      } else {
        a_re = y[i27].re / y[i27].im;
        a_im = y[i27].im + a_re * y[i27].re;
        h[i27].re = (a_re * h[i27].re + h[i27].im) / a_im;
        h[i27].im = (a_re * h[i27].im - h_re) / a_im;
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
  int i8;
  static const char cv12[8] = { 'o', 'n', 'e', 's', 'i', 'd', 'e', 'd' };

  double b_w[2048];
  static const char cv13[7] = { 'o', 'm', 'i', 't', 't', 'e', 'd' };

  /*  Cast to enforce precision rules */
  options.Fs = Fs;
  memcpy(&options.nfft[0], &w[0], sizeof(double) << 11);
  memcpy(&options.w[0], &w[0], sizeof(double) << 11);

  /*  Remaining are default or for advanced use */
  options.fvflag = 1.0;
  for (i8 = 0; i8 < 8; i8++) {
    options.range[i8] = cv12[i8];
  }

  options.centerdc = 0.0;
  for (i8 = 0; i8 < 7; i8++) {
    options.configlevel[i8] = cv13[i8];
  }

  firfreqz(&options, hh, b_w);

  /*  Generate the default structure to pass to freqzplot */
  /*  If rad/sample, Fs is empty */
}

/*
 * Make b a row
 * Arguments    : const double b[57]
 *                const struct_T *options
 *                creal_T h[2048]
 *                double w[2048]
 * Return Type  : void
 */
static void g_firfreqz(const double b[57], const struct_T *options, creal_T h
  [2048], double w[2048])
{
  int i78;
  creal_T dcv7[2048];
  double re_tmp;
  double re;
  creal_T s_tmp[2048];
  double h_re;
  double brm;

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
  for (i78 = 0; i78 < 2048; i78++) {
    w[i78] = options->w[i78];
    re_tmp = 6.2831853071795862 * options->w[i78] / options->Fs;
    re = re_tmp * 0.0;
    s_tmp[i78].re = re;
    s_tmp[i78].im = re_tmp;
    dcv7[i78].re = re;
    dcv7[i78].im = re_tmp;
  }

  b_exp(dcv7);
  f_polyval(b, dcv7, h);
  for (i78 = 0; i78 < 2048; i78++) {
    s_tmp[i78].re *= 56.0;
    s_tmp[i78].im *= 56.0;
  }

  b_exp(s_tmp);
  for (i78 = 0; i78 < 2048; i78++) {
    h_re = h[i78].re;
    if (s_tmp[i78].im == 0.0) {
      if (h[i78].im == 0.0) {
        h[i78].re /= s_tmp[i78].re;
        h[i78].im = 0.0;
      } else if (h[i78].re == 0.0) {
        h[i78].re = 0.0;
        h[i78].im /= s_tmp[i78].re;
      } else {
        h[i78].re /= s_tmp[i78].re;
        h[i78].im /= s_tmp[i78].re;
      }
    } else if (s_tmp[i78].re == 0.0) {
      if (h[i78].re == 0.0) {
        h[i78].re = h[i78].im / s_tmp[i78].im;
        h[i78].im = 0.0;
      } else if (h[i78].im == 0.0) {
        h[i78].re = 0.0;
        h[i78].im = -(h_re / s_tmp[i78].im);
      } else {
        h[i78].re = h[i78].im / s_tmp[i78].im;
        h[i78].im = -(h_re / s_tmp[i78].im);
      }
    } else {
      brm = fabs(s_tmp[i78].re);
      re = fabs(s_tmp[i78].im);
      if (brm > re) {
        re = s_tmp[i78].im / s_tmp[i78].re;
        re_tmp = s_tmp[i78].re + re * s_tmp[i78].im;
        h[i78].re = (h[i78].re + re * h[i78].im) / re_tmp;
        h[i78].im = (h[i78].im - re * h_re) / re_tmp;
      } else if (re == brm) {
        if (s_tmp[i78].re > 0.0) {
          re = 0.5;
        } else {
          re = -0.5;
        }

        if (s_tmp[i78].im > 0.0) {
          re_tmp = 0.5;
        } else {
          re_tmp = -0.5;
        }

        h[i78].re = (h[i78].re * re + h[i78].im * re_tmp) / brm;
        h[i78].im = (h[i78].im * re - h_re * re_tmp) / brm;
      } else {
        re = s_tmp[i78].re / s_tmp[i78].im;
        re_tmp = s_tmp[i78].im + re * s_tmp[i78].re;
        h[i78].re = (re * h[i78].re + h[i78].im) / re_tmp;
        h[i78].im = (re * h[i78].im - h_re) / re_tmp;
      }
    }
  }
}

/*
 * FREQZ_CG Frequency response of digital filter with codegen support
 *
 *  This function is based on 'freqz' by The MathWorks Inc.
 * Arguments    : const double b[57]
 *                const double w[2048]
 *                double Fs
 *                creal_T hh[2048]
 * Return Type  : void
 */
static void g_freqz_cg(const double b[57], const double w[2048], double Fs,
  creal_T hh[2048])
{
  static struct_T options;
  int i20;
  static const char cv24[8] = { 'o', 'n', 'e', 's', 'i', 'd', 'e', 'd' };

  double b_w[2048];
  static const char cv25[7] = { 'o', 'm', 'i', 't', 't', 'e', 'd' };

  /*  Cast to enforce precision rules */
  options.Fs = Fs;
  memcpy(&options.nfft[0], &w[0], sizeof(double) << 11);
  memcpy(&options.w[0], &w[0], sizeof(double) << 11);

  /*  Remaining are default or for advanced use */
  options.fvflag = 1.0;
  for (i20 = 0; i20 < 8; i20++) {
    options.range[i20] = cv24[i20];
  }

  options.centerdc = 0.0;
  for (i20 = 0; i20 < 7; i20++) {
    options.configlevel[i20] = cv25[i20];
  }

  g_firfreqz(b, &options, hh, b_w);

  /*  Generate the default structure to pass to freqzplot */
  /*  If rad/sample, Fs is empty */
}

/*
 * Arguments    : const double p[43]
 *                const creal_T x[2048]
 *                creal_T y[2048]
 * Return Type  : void
 */
static void g_polyval(const double p[43], const creal_T x[2048], creal_T y[2048])
{
  int i21;
  int k;
  double p_re;
  double x_im;
  for (i21 = 0; i21 < 2048; i21++) {
    y[i21].re = p[0];
    y[i21].im = 0.0;
  }

  for (k = 0; k < 42; k++) {
    p_re = p[k + 1];
    for (i21 = 0; i21 < 2048; i21++) {
      x_im = x[i21].re * y[i21].im + x[i21].im * y[i21].re;
      y[i21].re = (x[i21].re * y[i21].re - x[i21].im * y[i21].im) + p_re;
      y[i21].im = x_im;
    }
  }
}

/*
 * Arguments    : const double o[15]
 *                double u[43]
 * Return Type  : void
 */
static void g_us(const double o[15], double u[43])
{
  double b_u[45];
  int ix;
  int iy;
  int k;
  memset(&b_u[0], 0, 45U * sizeof(double));
  ix = 0;
  iy = 0;
  for (k = 0; k < 15; k++) {
    b_u[iy] = o[ix];
    ix++;
    iy += 3;
  }

  memcpy(&u[0], &b_u[0], 43U * sizeof(double));
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
  bool b_bool;
  int kstr;
  int exitg1;
  double dv2[15];
  double dv3[7];
  double tmp_data[29];
  int tmp_size[2];
  double dv4[29];
  double dv5[29];
  double dv6[13];
  double dv7[57];
  double dv8[43];
  double dv9[19];
  double dv10[85];
  static const char cv1[4] = { '2', '1', '1', '1' };

  double dv11[7];
  double dv12[13];
  double dv13[19];
  static creal_T d2[2048];
  static creal_T d3[2048];
  double combinedResponse_re;
  double combinedResponse_im;
  static const char cv2[4] = { '1', '2', '1', '1' };

  static const char cv3[4] = { '1', '1', '2', '1' };

  static const char cv4[4] = { '2', '2', '1', '1' };

  static const char cv5[4] = { '2', '1', '2', '1' };

  static const char cv6[4] = { '1', '2', '2', '1' };

  static const char cv7[4] = { '2', '2', '2', '1' };

  static const char cv8[4] = { '1', '1', '1', '3' };

  static const char cv9[4] = { '2', '1', '1', '3' };

  static const char cv10[4] = { '1', '2', '1', '3' };

  static const char cv11[4] = { '2', '2', '1', '3' };

  /*  Cast */
  b_bool = false;
  kstr = 0;
  do {
    exitg1 = 0;
    if (kstr < 4) {
      if (enables[kstr] != '1') {
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
    kstr = 0;
  } else {
    b_bool = false;
    kstr = 0;
    do {
      exitg1 = 0;
      if (kstr < 4) {
        if (enables[kstr] != cv1[kstr]) {
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
      b_bool = false;
      kstr = 0;
      do {
        exitg1 = 0;
        if (kstr < 4) {
          if (enables[kstr] != cv2[kstr]) {
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
        kstr = 2;
      } else {
        b_bool = false;
        kstr = 0;
        do {
          exitg1 = 0;
          if (kstr < 4) {
            if (enables[kstr] != cv3[kstr]) {
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
          kstr = 3;
        } else {
          b_bool = false;
          kstr = 0;
          do {
            exitg1 = 0;
            if (kstr < 4) {
              if (enables[kstr] != cv4[kstr]) {
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
            kstr = 4;
          } else {
            b_bool = false;
            kstr = 0;
            do {
              exitg1 = 0;
              if (kstr < 4) {
                if (enables[kstr] != cv5[kstr]) {
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
              kstr = 5;
            } else {
              b_bool = false;
              kstr = 0;
              do {
                exitg1 = 0;
                if (kstr < 4) {
                  if (enables[kstr] != cv6[kstr]) {
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
                kstr = 6;
              } else {
                b_bool = false;
                kstr = 0;
                do {
                  exitg1 = 0;
                  if (kstr < 4) {
                    if (enables[kstr] != cv7[kstr]) {
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
                  kstr = 7;
                } else {
                  b_bool = false;
                  kstr = 0;
                  do {
                    exitg1 = 0;
                    if (kstr < 4) {
                      if (enables[kstr] != cv8[kstr]) {
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
                    kstr = 8;
                  } else {
                    b_bool = false;
                    kstr = 0;
                    do {
                      exitg1 = 0;
                      if (kstr < 4) {
                        if (enables[kstr] != cv9[kstr]) {
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
                      kstr = 9;
                    } else {
                      b_bool = false;
                      kstr = 0;
                      do {
                        exitg1 = 0;
                        if (kstr < 4) {
                          if (enables[kstr] != cv10[kstr]) {
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
                        kstr = 10;
                      } else {
                        b_bool = false;
                        kstr = 0;
                        do {
                          exitg1 = 0;
                          if (kstr < 4) {
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

   case 6:
    /*  Hb3,Hb2 */
    e_us(hb2_coeff, dv6);
    f_freqz_cg(dv6, w, Fs, combinedResponse);
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
    f_us(hb1_coeff, dv7);
    g_freqz_cg(dv7, w, Fs, combinedResponse);
    e_us(hb2_coeff, dv12);
    f_freqz_cg(dv12, w, Fs, d2);
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

    /*  RECHECK ALL DEC BY 3     */
    break;

   case 9:
    /*  Dec/Int3,Hb1 */
    g_us(hb1_coeff, dv8);
    h_freqz_cg(dv8, w, Fs, combinedResponse);
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
    h_us(hb2_coeff, dv9);
    i_freqz_cg(dv9, w, Fs, combinedResponse);
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
    /*  Dec/Int3,Hb2,Hb1 {Hm4,Hm2c34,Hm1} */
    i_us(hb1_coeff, dv10);
    j_freqz_cg(dv10, w, Fs, combinedResponse);
    h_us(hb2_coeff, dv13);
    i_freqz_cg(dv13, w, Fs, d2);
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
 * Make b a row
 * Arguments    : const double b[43]
 *                const struct_T *options
 *                creal_T h[2048]
 *                double w[2048]
 * Return Type  : void
 */
static void h_firfreqz(const double b[43], const struct_T *options, creal_T h
  [2048], double w[2048])
{
  int i79;
  creal_T dcv8[2048];
  double re_tmp;
  double re;
  creal_T s_tmp[2048];
  double h_re;
  double brm;

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
  for (i79 = 0; i79 < 2048; i79++) {
    w[i79] = options->w[i79];
    re_tmp = 6.2831853071795862 * options->w[i79] / options->Fs;
    re = re_tmp * 0.0;
    s_tmp[i79].re = re;
    s_tmp[i79].im = re_tmp;
    dcv8[i79].re = re;
    dcv8[i79].im = re_tmp;
  }

  b_exp(dcv8);
  g_polyval(b, dcv8, h);
  for (i79 = 0; i79 < 2048; i79++) {
    s_tmp[i79].re *= 42.0;
    s_tmp[i79].im *= 42.0;
  }

  b_exp(s_tmp);
  for (i79 = 0; i79 < 2048; i79++) {
    h_re = h[i79].re;
    if (s_tmp[i79].im == 0.0) {
      if (h[i79].im == 0.0) {
        h[i79].re /= s_tmp[i79].re;
        h[i79].im = 0.0;
      } else if (h[i79].re == 0.0) {
        h[i79].re = 0.0;
        h[i79].im /= s_tmp[i79].re;
      } else {
        h[i79].re /= s_tmp[i79].re;
        h[i79].im /= s_tmp[i79].re;
      }
    } else if (s_tmp[i79].re == 0.0) {
      if (h[i79].re == 0.0) {
        h[i79].re = h[i79].im / s_tmp[i79].im;
        h[i79].im = 0.0;
      } else if (h[i79].im == 0.0) {
        h[i79].re = 0.0;
        h[i79].im = -(h_re / s_tmp[i79].im);
      } else {
        h[i79].re = h[i79].im / s_tmp[i79].im;
        h[i79].im = -(h_re / s_tmp[i79].im);
      }
    } else {
      brm = fabs(s_tmp[i79].re);
      re = fabs(s_tmp[i79].im);
      if (brm > re) {
        re = s_tmp[i79].im / s_tmp[i79].re;
        re_tmp = s_tmp[i79].re + re * s_tmp[i79].im;
        h[i79].re = (h[i79].re + re * h[i79].im) / re_tmp;
        h[i79].im = (h[i79].im - re * h_re) / re_tmp;
      } else if (re == brm) {
        if (s_tmp[i79].re > 0.0) {
          re = 0.5;
        } else {
          re = -0.5;
        }

        if (s_tmp[i79].im > 0.0) {
          re_tmp = 0.5;
        } else {
          re_tmp = -0.5;
        }

        h[i79].re = (h[i79].re * re + h[i79].im * re_tmp) / brm;
        h[i79].im = (h[i79].im * re - h_re * re_tmp) / brm;
      } else {
        re = s_tmp[i79].re / s_tmp[i79].im;
        re_tmp = s_tmp[i79].im + re * s_tmp[i79].re;
        h[i79].re = (re * h[i79].re + h[i79].im) / re_tmp;
        h[i79].im = (re * h[i79].im - h_re) / re_tmp;
      }
    }
  }
}

/*
 * FREQZ_CG Frequency response of digital filter with codegen support
 *
 *  This function is based on 'freqz' by The MathWorks Inc.
 * Arguments    : const double b[43]
 *                const double w[2048]
 *                double Fs
 *                creal_T hh[2048]
 * Return Type  : void
 */
static void h_freqz_cg(const double b[43], const double w[2048], double Fs,
  creal_T hh[2048])
{
  static struct_T options;
  int i22;
  static const char cv26[8] = { 'o', 'n', 'e', 's', 'i', 'd', 'e', 'd' };

  double b_w[2048];
  static const char cv27[7] = { 'o', 'm', 'i', 't', 't', 'e', 'd' };

  /*  Cast to enforce precision rules */
  options.Fs = Fs;
  memcpy(&options.nfft[0], &w[0], sizeof(double) << 11);
  memcpy(&options.w[0], &w[0], sizeof(double) << 11);

  /*  Remaining are default or for advanced use */
  options.fvflag = 1.0;
  for (i22 = 0; i22 < 8; i22++) {
    options.range[i22] = cv26[i22];
  }

  options.centerdc = 0.0;
  for (i22 = 0; i22 < 7; i22++) {
    options.configlevel[i22] = cv27[i22];
  }

  h_firfreqz(b, &options, hh, b_w);

  /*  Generate the default structure to pass to freqzplot */
  /*  If rad/sample, Fs is empty */
}

/*
 * Arguments    : const double p[19]
 *                const creal_T x[2048]
 *                creal_T y[2048]
 * Return Type  : void
 */
static void h_polyval(const double p[19], const creal_T x[2048], creal_T y[2048])
{
  int i23;
  int k;
  double p_re;
  double x_im;
  for (i23 = 0; i23 < 2048; i23++) {
    y[i23].re = p[0];
    y[i23].im = 0.0;
  }

  for (k = 0; k < 18; k++) {
    p_re = p[k + 1];
    for (i23 = 0; i23 < 2048; i23++) {
      x_im = x[i23].re * y[i23].im + x[i23].im * y[i23].re;
      y[i23].re = (x[i23].re * y[i23].re - x[i23].im * y[i23].im) + p_re;
      y[i23].im = x_im;
    }
  }
}

/*
 * Arguments    : const double o[7]
 *                double u[19]
 * Return Type  : void
 */
static void h_us(const double o[7], double u[19])
{
  double b_u[21];
  int ix;
  int iy;
  int k;
  memset(&b_u[0], 0, 21U * sizeof(double));
  ix = 0;
  iy = 0;
  for (k = 0; k < 7; k++) {
    b_u[iy] = o[ix];
    ix++;
    iy += 3;
  }

  memcpy(&u[0], &b_u[0], 19U * sizeof(double));
}

/*
 * Make b a row
 * Arguments    : const double b[19]
 *                const struct_T *options
 *                creal_T h[2048]
 *                double w[2048]
 * Return Type  : void
 */
static void i_firfreqz(const double b[19], const struct_T *options, creal_T h
  [2048], double w[2048])
{
  int i80;
  creal_T dcv9[2048];
  double re_tmp;
  double re;
  creal_T s_tmp[2048];
  double h_re;
  double brm;

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
  for (i80 = 0; i80 < 2048; i80++) {
    w[i80] = options->w[i80];
    re_tmp = 6.2831853071795862 * options->w[i80] / options->Fs;
    re = re_tmp * 0.0;
    s_tmp[i80].re = re;
    s_tmp[i80].im = re_tmp;
    dcv9[i80].re = re;
    dcv9[i80].im = re_tmp;
  }

  b_exp(dcv9);
  h_polyval(b, dcv9, h);
  for (i80 = 0; i80 < 2048; i80++) {
    s_tmp[i80].re *= 18.0;
    s_tmp[i80].im *= 18.0;
  }

  b_exp(s_tmp);
  for (i80 = 0; i80 < 2048; i80++) {
    h_re = h[i80].re;
    if (s_tmp[i80].im == 0.0) {
      if (h[i80].im == 0.0) {
        h[i80].re /= s_tmp[i80].re;
        h[i80].im = 0.0;
      } else if (h[i80].re == 0.0) {
        h[i80].re = 0.0;
        h[i80].im /= s_tmp[i80].re;
      } else {
        h[i80].re /= s_tmp[i80].re;
        h[i80].im /= s_tmp[i80].re;
      }
    } else if (s_tmp[i80].re == 0.0) {
      if (h[i80].re == 0.0) {
        h[i80].re = h[i80].im / s_tmp[i80].im;
        h[i80].im = 0.0;
      } else if (h[i80].im == 0.0) {
        h[i80].re = 0.0;
        h[i80].im = -(h_re / s_tmp[i80].im);
      } else {
        h[i80].re = h[i80].im / s_tmp[i80].im;
        h[i80].im = -(h_re / s_tmp[i80].im);
      }
    } else {
      brm = fabs(s_tmp[i80].re);
      re = fabs(s_tmp[i80].im);
      if (brm > re) {
        re = s_tmp[i80].im / s_tmp[i80].re;
        re_tmp = s_tmp[i80].re + re * s_tmp[i80].im;
        h[i80].re = (h[i80].re + re * h[i80].im) / re_tmp;
        h[i80].im = (h[i80].im - re * h_re) / re_tmp;
      } else if (re == brm) {
        if (s_tmp[i80].re > 0.0) {
          re = 0.5;
        } else {
          re = -0.5;
        }

        if (s_tmp[i80].im > 0.0) {
          re_tmp = 0.5;
        } else {
          re_tmp = -0.5;
        }

        h[i80].re = (h[i80].re * re + h[i80].im * re_tmp) / brm;
        h[i80].im = (h[i80].im * re - h_re * re_tmp) / brm;
      } else {
        re = s_tmp[i80].re / s_tmp[i80].im;
        re_tmp = s_tmp[i80].im + re * s_tmp[i80].re;
        h[i80].re = (re * h[i80].re + h[i80].im) / re_tmp;
        h[i80].im = (re * h[i80].im - h_re) / re_tmp;
      }
    }
  }
}

/*
 * FREQZ_CG Frequency response of digital filter with codegen support
 *
 *  This function is based on 'freqz' by The MathWorks Inc.
 * Arguments    : const double b[19]
 *                const double w[2048]
 *                double Fs
 *                creal_T hh[2048]
 * Return Type  : void
 */
static void i_freqz_cg(const double b[19], const double w[2048], double Fs,
  creal_T hh[2048])
{
  static struct_T options;
  int i24;
  static const char cv28[8] = { 'o', 'n', 'e', 's', 'i', 'd', 'e', 'd' };

  double b_w[2048];
  static const char cv29[7] = { 'o', 'm', 'i', 't', 't', 'e', 'd' };

  /*  Cast to enforce precision rules */
  options.Fs = Fs;
  memcpy(&options.nfft[0], &w[0], sizeof(double) << 11);
  memcpy(&options.w[0], &w[0], sizeof(double) << 11);

  /*  Remaining are default or for advanced use */
  options.fvflag = 1.0;
  for (i24 = 0; i24 < 8; i24++) {
    options.range[i24] = cv28[i24];
  }

  options.centerdc = 0.0;
  for (i24 = 0; i24 < 7; i24++) {
    options.configlevel[i24] = cv29[i24];
  }

  i_firfreqz(b, &options, hh, b_w);

  /*  Generate the default structure to pass to freqzplot */
  /*  If rad/sample, Fs is empty */
}

/*
 * Arguments    : const double p[85]
 *                const creal_T x[2048]
 *                creal_T y[2048]
 * Return Type  : void
 */
static void i_polyval(const double p[85], const creal_T x[2048], creal_T y[2048])
{
  int i25;
  int k;
  double p_re;
  double x_im;
  for (i25 = 0; i25 < 2048; i25++) {
    y[i25].re = p[0];
    y[i25].im = 0.0;
  }

  for (k = 0; k < 84; k++) {
    p_re = p[k + 1];
    for (i25 = 0; i25 < 2048; i25++) {
      x_im = x[i25].re * y[i25].im + x[i25].im * y[i25].re;
      y[i25].re = (x[i25].re * y[i25].re - x[i25].im * y[i25].im) + p_re;
      y[i25].im = x_im;
    }
  }
}

/*
 * Arguments    : const double o[15]
 *                double u[85]
 * Return Type  : void
 */
static void i_us(const double o[15], double u[85])
{
  double b_u[90];
  int ix;
  int iy;
  int k;
  memset(&b_u[0], 0, 90U * sizeof(double));
  ix = 0;
  iy = 0;
  for (k = 0; k < 15; k++) {
    b_u[iy] = o[ix];
    ix++;
    iy += 6;
  }

  memcpy(&u[0], &b_u[0], 85U * sizeof(double));
}

/*
 * Arguments    : const emxArray_real_T *varargin_1
 *                const emxArray_real_T *varargin_2
 *                const emxArray_real_T *varargin_3
 *                emxArray_real_T *Vq
 * Return Type  : void
 */
static void interp1(const emxArray_real_T *varargin_1, const emxArray_real_T
                    *varargin_2, const emxArray_real_T *varargin_3,
                    emxArray_real_T *Vq)
{
  emxArray_real_T *y;
  int j2;
  int outsize_idx_1;
  emxArray_real_T *x;
  int nx;
  int k;
  int exitg1;
  int low_ip1;
  double xtmp;
  int mid_i;
  emxInit_real_T(&y, 2);
  j2 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = varargin_2->size[1];
  emxEnsureCapacity_real_T(y, j2);
  outsize_idx_1 = varargin_2->size[0] * varargin_2->size[1];
  for (j2 = 0; j2 < outsize_idx_1; j2++) {
    y->data[j2] = varargin_2->data[j2];
  }

  emxInit_real_T(&x, 2);
  j2 = x->size[0] * x->size[1];
  x->size[0] = 1;
  x->size[1] = varargin_1->size[1];
  emxEnsureCapacity_real_T(x, j2);
  outsize_idx_1 = varargin_1->size[0] * varargin_1->size[1];
  for (j2 = 0; j2 < outsize_idx_1; j2++) {
    x->data[j2] = varargin_1->data[j2];
  }

  nx = varargin_1->size[1] - 1;
  outsize_idx_1 = varargin_3->size[1];
  j2 = Vq->size[0] * Vq->size[1];
  Vq->size[0] = 1;
  Vq->size[1] = outsize_idx_1;
  emxEnsureCapacity_real_T(Vq, j2);
  for (j2 = 0; j2 < outsize_idx_1; j2++) {
    Vq->data[j2] = rtNaN;
  }

  if (varargin_3->size[1] != 0) {
    k = 0;
    do {
      exitg1 = 0;
      if (k <= nx) {
        if (rtIsNaN(varargin_1->data[k])) {
          exitg1 = 1;
        } else {
          k++;
        }
      } else {
        if (varargin_1->data[1] < varargin_1->data[0]) {
          j2 = (nx + 1) >> 1;
          for (low_ip1 = 0; low_ip1 < j2; low_ip1++) {
            xtmp = x->data[low_ip1];
            outsize_idx_1 = nx - low_ip1;
            x->data[low_ip1] = x->data[outsize_idx_1];
            x->data[outsize_idx_1] = xtmp;
          }

          outsize_idx_1 = varargin_2->size[1] >> 1;
          for (low_ip1 = 0; low_ip1 < outsize_idx_1; low_ip1++) {
            j2 = (varargin_2->size[1] - low_ip1) - 1;
            xtmp = y->data[low_ip1];
            y->data[low_ip1] = y->data[j2];
            y->data[j2] = xtmp;
          }
        }

        outsize_idx_1 = varargin_3->size[1];
        for (k = 0; k < outsize_idx_1; k++) {
          if (rtIsNaN(varargin_3->data[k])) {
            Vq->data[k] = rtNaN;
          } else {
            if ((!(varargin_3->data[k] > x->data[x->size[1] - 1])) &&
                (!(varargin_3->data[k] < x->data[0]))) {
              j2 = x->size[1];
              nx = 1;
              low_ip1 = 2;
              while (j2 > low_ip1) {
                mid_i = (nx >> 1) + (j2 >> 1);
                if (((nx & 1) == 1) && ((j2 & 1) == 1)) {
                  mid_i++;
                }

                if (varargin_3->data[k] >= x->data[mid_i - 1]) {
                  nx = mid_i;
                  low_ip1 = mid_i + 1;
                } else {
                  j2 = mid_i;
                }
              }

              xtmp = (varargin_3->data[k] - x->data[nx - 1]) / (x->data[nx] -
                x->data[nx - 1]);
              if (xtmp == 0.0) {
                Vq->data[k] = y->data[nx - 1];
              } else if (xtmp == 1.0) {
                Vq->data[k] = y->data[nx];
              } else if (y->data[nx - 1] == y->data[nx]) {
                Vq->data[k] = y->data[nx - 1];
              } else {
                Vq->data[k] = (1.0 - xtmp) * y->data[nx - 1] + xtmp * y->data[nx];
              }
            }
          }
        }

        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  emxFree_real_T(&x);
  emxFree_real_T(&y);
}

/*
 * Make b a row
 * Arguments    : const double b[85]
 *                const struct_T *options
 *                creal_T h[2048]
 *                double w[2048]
 * Return Type  : void
 */
static void j_firfreqz(const double b[85], const struct_T *options, creal_T h
  [2048], double w[2048])
{
  int i81;
  creal_T dcv10[2048];
  double re_tmp;
  double re;
  creal_T s_tmp[2048];
  double h_re;
  double brm;

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
  for (i81 = 0; i81 < 2048; i81++) {
    w[i81] = options->w[i81];
    re_tmp = 6.2831853071795862 * options->w[i81] / options->Fs;
    re = re_tmp * 0.0;
    s_tmp[i81].re = re;
    s_tmp[i81].im = re_tmp;
    dcv10[i81].re = re;
    dcv10[i81].im = re_tmp;
  }

  b_exp(dcv10);
  i_polyval(b, dcv10, h);
  for (i81 = 0; i81 < 2048; i81++) {
    s_tmp[i81].re *= 84.0;
    s_tmp[i81].im *= 84.0;
  }

  b_exp(s_tmp);
  for (i81 = 0; i81 < 2048; i81++) {
    h_re = h[i81].re;
    if (s_tmp[i81].im == 0.0) {
      if (h[i81].im == 0.0) {
        h[i81].re /= s_tmp[i81].re;
        h[i81].im = 0.0;
      } else if (h[i81].re == 0.0) {
        h[i81].re = 0.0;
        h[i81].im /= s_tmp[i81].re;
      } else {
        h[i81].re /= s_tmp[i81].re;
        h[i81].im /= s_tmp[i81].re;
      }
    } else if (s_tmp[i81].re == 0.0) {
      if (h[i81].re == 0.0) {
        h[i81].re = h[i81].im / s_tmp[i81].im;
        h[i81].im = 0.0;
      } else if (h[i81].im == 0.0) {
        h[i81].re = 0.0;
        h[i81].im = -(h_re / s_tmp[i81].im);
      } else {
        h[i81].re = h[i81].im / s_tmp[i81].im;
        h[i81].im = -(h_re / s_tmp[i81].im);
      }
    } else {
      brm = fabs(s_tmp[i81].re);
      re = fabs(s_tmp[i81].im);
      if (brm > re) {
        re = s_tmp[i81].im / s_tmp[i81].re;
        re_tmp = s_tmp[i81].re + re * s_tmp[i81].im;
        h[i81].re = (h[i81].re + re * h[i81].im) / re_tmp;
        h[i81].im = (h[i81].im - re * h_re) / re_tmp;
      } else if (re == brm) {
        if (s_tmp[i81].re > 0.0) {
          re = 0.5;
        } else {
          re = -0.5;
        }

        if (s_tmp[i81].im > 0.0) {
          re_tmp = 0.5;
        } else {
          re_tmp = -0.5;
        }

        h[i81].re = (h[i81].re * re + h[i81].im * re_tmp) / brm;
        h[i81].im = (h[i81].im * re - h_re * re_tmp) / brm;
      } else {
        re = s_tmp[i81].re / s_tmp[i81].im;
        re_tmp = s_tmp[i81].im + re * s_tmp[i81].re;
        h[i81].re = (re * h[i81].re + h[i81].im) / re_tmp;
        h[i81].im = (re * h[i81].im - h_re) / re_tmp;
      }
    }
  }
}

/*
 * FREQZ_CG Frequency response of digital filter with codegen support
 *
 *  This function is based on 'freqz' by The MathWorks Inc.
 * Arguments    : const double b[85]
 *                const double w[2048]
 *                double Fs
 *                creal_T hh[2048]
 * Return Type  : void
 */
static void j_freqz_cg(const double b[85], const double w[2048], double Fs,
  creal_T hh[2048])
{
  static struct_T options;
  int i26;
  static const char cv30[8] = { 'o', 'n', 'e', 's', 'i', 'd', 'e', 'd' };

  double b_w[2048];
  static const char cv31[7] = { 'o', 'm', 'i', 't', 't', 'e', 'd' };

  /*  Cast to enforce precision rules */
  options.Fs = Fs;
  memcpy(&options.nfft[0], &w[0], sizeof(double) << 11);
  memcpy(&options.w[0], &w[0], sizeof(double) << 11);

  /*  Remaining are default or for advanced use */
  options.fvflag = 1.0;
  for (i26 = 0; i26 < 8; i26++) {
    options.range[i26] = cv30[i26];
  }

  options.centerdc = 0.0;
  for (i26 = 0; i26 < 7; i26++) {
    options.configlevel[i26] = cv31[i26];
  }

  j_firfreqz(b, &options, hh, b_w);

  /*  Generate the default structure to pass to freqzplot */
  /*  If rad/sample, Fs is empty */
}

/*
 * Arguments    : const emxArray_real_T *p
 *                const emxArray_creal_T *x
 *                emxArray_creal_T *y
 * Return Type  : void
 */
static void j_polyval(const emxArray_real_T *p, const emxArray_creal_T *x,
                      emxArray_creal_T *y)
{
  int i34;
  int loop_ub;
  int k;
  int i35;
  double p_re;
  double x_re;
  double x_im;
  i34 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = x->size[1];
  emxEnsureCapacity_creal_T(y, i34);
  if ((y->size[1] == 0) || (p->size[1] == 0)) {
  } else {
    i34 = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity_creal_T(y, i34);
    loop_ub = y->size[1];
    for (i34 = 0; i34 < loop_ub; i34++) {
      y->data[i34].re = p->data[0];
      y->data[i34].im = 0.0;
    }

    i34 = p->size[1];
    for (k = 0; k <= i34 - 2; k++) {
      i35 = x->size[0] * x->size[1];
      loop_ub = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = x->size[1];
      emxEnsureCapacity_creal_T(y, loop_ub);
      p_re = p->data[k + 1];
      loop_ub = i35 - 1;
      for (i35 = 0; i35 <= loop_ub; i35++) {
        x_re = x->data[i35].re * y->data[i35].re - x->data[i35].im * y->data[i35]
          .im;
        x_im = x->data[i35].re * y->data[i35].im + x->data[i35].im * y->data[i35]
          .re;
        y->data[i35].re = x_re + p_re;
        y->data[i35].im = x_im;
      }
    }
  }
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
static void k_freqz_cg(const emxArray_real_T *w, double Fs, emxArray_creal_T *hh)
{
  emxArray_real_T *digw;
  int i29;
  int loop_ub;
  emxArray_creal_T *y;
  emxArray_creal_T *r3;
  bool b1;
  double re;
  double im;
  emxInit_real_T(&digw, 2);

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
  i29 = digw->size[0] * digw->size[1];
  digw->size[0] = 1;
  digw->size[1] = w->size[1];
  emxEnsureCapacity_real_T(digw, i29);
  loop_ub = w->size[0] * w->size[1];
  for (i29 = 0; i29 < loop_ub; i29++) {
    digw->data[i29] = 6.2831853071795862 * w->data[i29] / Fs;
  }

  emxInit_creal_T(&y, 2);
  emxInit_creal_T(&r3, 2);

  /*  Convert from Hz to rad/sample for computational purposes */
  /*  Digital frequency must be used for this calculation */
  i29 = r3->size[0] * r3->size[1];
  r3->size[0] = 1;
  r3->size[1] = digw->size[1];
  emxEnsureCapacity_creal_T(r3, i29);
  loop_ub = digw->size[0] * digw->size[1];
  for (i29 = 0; i29 < loop_ub; i29++) {
    r3->data[i29].re = digw->data[i29] * 0.0;
    r3->data[i29].im = digw->data[i29];
  }

  c_exp(r3);
  i29 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = r3->size[1];
  emxEnsureCapacity_creal_T(y, i29);
  b1 = (y->size[1] == 0);
  if (!b1) {
    i29 = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity_creal_T(y, i29);
    loop_ub = y->size[1];
    for (i29 = 0; i29 < loop_ub; i29++) {
      y->data[i29].re = 1.0;
      y->data[i29].im = 0.0;
    }
  }

  i29 = r3->size[0] * r3->size[1];
  r3->size[0] = 1;
  r3->size[1] = digw->size[1];
  emxEnsureCapacity_creal_T(r3, i29);
  loop_ub = digw->size[0] * digw->size[1];
  for (i29 = 0; i29 < loop_ub; i29++) {
    re = digw->data[i29] * 0.0;
    im = digw->data[i29];
    r3->data[i29].re = 0.0 * re;
    r3->data[i29].im = 0.0 * im;
  }

  emxFree_real_T(&digw);
  c_exp(r3);
  rdivide_helper(y, r3, hh);

  /*  Generate the default structure to pass to freqzplot */
  /*  If rad/sample, Fs is empty */
  emxFree_creal_T(&r3);
  emxFree_creal_T(&y);
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
static void l_freqz_cg(const double b[15], const emxArray_real_T *w, double Fs,
  emxArray_creal_T *hh)
{
  emxArray_real_T *digw;
  int i31;
  int loop_ub;
  emxArray_creal_T *s;
  emxArray_creal_T *y;
  bool b2;
  int k;
  double b_re;
  double s_re;
  double s_im;
  emxInit_real_T(&digw, 2);

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
  i31 = digw->size[0] * digw->size[1];
  digw->size[0] = 1;
  digw->size[1] = w->size[1];
  emxEnsureCapacity_real_T(digw, i31);
  loop_ub = w->size[0] * w->size[1];
  for (i31 = 0; i31 < loop_ub; i31++) {
    digw->data[i31] = 6.2831853071795862 * w->data[i31] / Fs;
  }

  emxInit_creal_T(&s, 2);

  /*  Convert from Hz to rad/sample for computational purposes */
  i31 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = digw->size[1];
  emxEnsureCapacity_creal_T(s, i31);
  loop_ub = digw->size[0] * digw->size[1];
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
  emxEnsureCapacity_creal_T(y, i31);
  b2 = (y->size[1] == 0);
  if (!b2) {
    i31 = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity_creal_T(y, i31);
    loop_ub = y->size[1];
    for (i31 = 0; i31 < loop_ub; i31++) {
      y->data[i31].re = b[0];
      y->data[i31].im = 0.0;
    }

    for (k = 0; k < 14; k++) {
      i31 = s->size[0] * s->size[1];
      loop_ub = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = s->size[1];
      emxEnsureCapacity_creal_T(y, loop_ub);
      b_re = b[k + 1];
      loop_ub = i31 - 1;
      for (i31 = 0; i31 <= loop_ub; i31++) {
        s_re = s->data[i31].re * y->data[i31].re - s->data[i31].im * y->data[i31]
          .im;
        s_im = s->data[i31].re * y->data[i31].im + s->data[i31].im * y->data[i31]
          .re;
        y->data[i31].re = s_re + b_re;
        y->data[i31].im = s_im;
      }
    }
  }

  i31 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = digw->size[1];
  emxEnsureCapacity_creal_T(s, i31);
  loop_ub = digw->size[0] * digw->size[1];
  for (i31 = 0; i31 < loop_ub; i31++) {
    b_re = digw->data[i31] * 0.0;
    s_re = digw->data[i31];
    s->data[i31].re = 14.0 * b_re;
    s->data[i31].im = 14.0 * s_re;
  }

  emxFree_real_T(&digw);
  c_exp(s);
  rdivide_helper(y, s, hh);

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
 *                double wo
 *                emxArray_creal_T *at
 *                emxArray_real_T *bt
 *                double *dt
 * Return Type  : void
 */
static void lp2lp_cg(const emxArray_creal_T *a, const emxArray_real_T *b, double
                     wo, emxArray_creal_T *at, emxArray_real_T *bt, double *dt)
{
  int i66;
  int loop_ub;

  /*  Transform lowpass to lowpass */
  i66 = at->size[0] * at->size[1];
  at->size[0] = a->size[0];
  at->size[1] = a->size[1];
  emxEnsureCapacity_creal_T(at, i66);
  loop_ub = a->size[0] * a->size[1];
  for (i66 = 0; i66 < loop_ub; i66++) {
    at->data[i66].re = wo * a->data[i66].re;
    at->data[i66].im = wo * a->data[i66].im;
  }

  i66 = bt->size[0];
  bt->size[0] = b->size[0];
  emxEnsureCapacity_real_T(bt, i66);
  loop_ub = b->size[0];
  for (i66 = 0; i66 < loop_ub; i66++) {
    bt->data[i66] = wo * b->data[i66];
  }

  *dt = 0.0;
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
static void m_freqz_cg(const double b[7], const emxArray_real_T *w, double Fs,
  emxArray_creal_T *hh)
{
  emxArray_real_T *digw;
  int i32;
  int loop_ub;
  emxArray_creal_T *s;
  emxArray_creal_T *y;
  bool b3;
  int k;
  double b_re;
  double s_re;
  double s_im;
  emxInit_real_T(&digw, 2);

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
  i32 = digw->size[0] * digw->size[1];
  digw->size[0] = 1;
  digw->size[1] = w->size[1];
  emxEnsureCapacity_real_T(digw, i32);
  loop_ub = w->size[0] * w->size[1];
  for (i32 = 0; i32 < loop_ub; i32++) {
    digw->data[i32] = 6.2831853071795862 * w->data[i32] / Fs;
  }

  emxInit_creal_T(&s, 2);

  /*  Convert from Hz to rad/sample for computational purposes */
  i32 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = digw->size[1];
  emxEnsureCapacity_creal_T(s, i32);
  loop_ub = digw->size[0] * digw->size[1];
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
  emxEnsureCapacity_creal_T(y, i32);
  b3 = (y->size[1] == 0);
  if (!b3) {
    i32 = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity_creal_T(y, i32);
    loop_ub = y->size[1];
    for (i32 = 0; i32 < loop_ub; i32++) {
      y->data[i32].re = b[0];
      y->data[i32].im = 0.0;
    }

    for (k = 0; k < 6; k++) {
      i32 = s->size[0] * s->size[1];
      loop_ub = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = s->size[1];
      emxEnsureCapacity_creal_T(y, loop_ub);
      b_re = b[k + 1];
      loop_ub = i32 - 1;
      for (i32 = 0; i32 <= loop_ub; i32++) {
        s_re = s->data[i32].re * y->data[i32].re - s->data[i32].im * y->data[i32]
          .im;
        s_im = s->data[i32].re * y->data[i32].im + s->data[i32].im * y->data[i32]
          .re;
        y->data[i32].re = s_re + b_re;
        y->data[i32].im = s_im;
      }
    }
  }

  i32 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = digw->size[1];
  emxEnsureCapacity_creal_T(s, i32);
  loop_ub = digw->size[0] * digw->size[1];
  for (i32 = 0; i32 < loop_ub; i32++) {
    b_re = digw->data[i32] * 0.0;
    s_re = digw->data[i32];
    s->data[i32].re = 6.0 * b_re;
    s->data[i32].im = 6.0 * s_re;
  }

  emxFree_real_T(&digw);
  c_exp(s);
  rdivide_helper(y, s, hh);

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
 * Arguments    : double a
 *                double b
 * Return Type  : double
 */
static double mpower(double a, double b)
{
  return rt_powd_snf(a, b);
}

/*
 * FREQZ_CG Frequency response of digital filter with codegen support
 *
 *  This function is based on 'freqz' by The MathWorks Inc.
 * Arguments    : const emxArray_real_T *b
 *                const emxArray_real_T *w
 *                double Fs
 *                emxArray_creal_T *hh
 * Return Type  : void
 */
static void n_freqz_cg(const emxArray_real_T *b, const emxArray_real_T *w,
  double Fs, emxArray_creal_T *hh)
{
  emxArray_real_T *b_b;
  int b_idx_0;
  int i33;
  emxArray_real_T *digw;
  emxArray_creal_T *r4;
  emxArray_creal_T *r5;
  double c_b;
  double re;
  double im;
  emxInit_real_T(&b_b, 2);

  /*  Cast to enforce precision rules */
  /*  Remaining are default or for advanced use */
  /*  Make b a row */
  /* -------------------------------------------------------------------------- */
  b_idx_0 = b->size[1];
  i33 = b_b->size[0] * b_b->size[1];
  b_b->size[0] = 1;
  b_b->size[1] = b_idx_0;
  emxEnsureCapacity_real_T(b_b, i33);
  for (i33 = 0; i33 < b_idx_0; i33++) {
    b_b->data[i33] = b->data[i33];
  }

  emxInit_real_T(&digw, 2);

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
  i33 = digw->size[0] * digw->size[1];
  digw->size[0] = 1;
  digw->size[1] = w->size[1];
  emxEnsureCapacity_real_T(digw, i33);
  b_idx_0 = w->size[0] * w->size[1];
  for (i33 = 0; i33 < b_idx_0; i33++) {
    digw->data[i33] = 6.2831853071795862 * w->data[i33] / Fs;
  }

  emxInit_creal_T(&r4, 2);

  /*  Convert from Hz to rad/sample for computational purposes */
  /*  Digital frequency must be used for this calculation */
  i33 = r4->size[0] * r4->size[1];
  r4->size[0] = 1;
  r4->size[1] = digw->size[1];
  emxEnsureCapacity_creal_T(r4, i33);
  b_idx_0 = digw->size[0] * digw->size[1];
  for (i33 = 0; i33 < b_idx_0; i33++) {
    r4->data[i33].re = digw->data[i33] * 0.0;
    r4->data[i33].im = digw->data[i33];
  }

  emxInit_creal_T(&r5, 2);
  c_exp(r4);
  j_polyval(b_b, r4, r5);
  i33 = r4->size[0] * r4->size[1];
  r4->size[0] = 1;
  r4->size[1] = digw->size[1];
  emxEnsureCapacity_creal_T(r4, i33);
  c_b = (double)b_b->size[1] - 1.0;
  b_idx_0 = digw->size[0] * digw->size[1];
  emxFree_real_T(&b_b);
  for (i33 = 0; i33 < b_idx_0; i33++) {
    re = digw->data[i33] * 0.0;
    im = digw->data[i33];
    r4->data[i33].re = c_b * re;
    r4->data[i33].im = c_b * im;
  }

  emxFree_real_T(&digw);
  c_exp(r4);
  rdivide_helper(r5, r4, hh);

  /*  Generate the default structure to pass to freqzplot */
  /*  If rad/sample, Fs is empty */
  emxFree_creal_T(&r5);
  emxFree_creal_T(&r4);
}

/*
 * FREQZ_CG Frequency response of digital filter with codegen support
 *
 *  This function is based on 'freqz' by The MathWorks Inc.
 * Arguments    : const double b[29]
 *                const emxArray_real_T *w
 *                double Fs
 *                emxArray_creal_T *hh
 * Return Type  : void
 */
static void o_freqz_cg(const double b[29], const emxArray_real_T *w, double Fs,
  emxArray_creal_T *hh)
{
  emxArray_real_T *digw;
  int i36;
  int loop_ub;
  emxArray_creal_T *s;
  emxArray_creal_T *y;
  bool b4;
  int k;
  double b_re;
  double s_re;
  double s_im;
  emxInit_real_T(&digw, 2);

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
  i36 = digw->size[0] * digw->size[1];
  digw->size[0] = 1;
  digw->size[1] = w->size[1];
  emxEnsureCapacity_real_T(digw, i36);
  loop_ub = w->size[0] * w->size[1];
  for (i36 = 0; i36 < loop_ub; i36++) {
    digw->data[i36] = 6.2831853071795862 * w->data[i36] / Fs;
  }

  emxInit_creal_T(&s, 2);

  /*  Convert from Hz to rad/sample for computational purposes */
  i36 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = digw->size[1];
  emxEnsureCapacity_creal_T(s, i36);
  loop_ub = digw->size[0] * digw->size[1];
  for (i36 = 0; i36 < loop_ub; i36++) {
    s->data[i36].re = digw->data[i36] * 0.0;
    s->data[i36].im = digw->data[i36];
  }

  emxInit_creal_T(&y, 2);
  c_exp(s);

  /*  Digital frequency must be used for this calculation */
  i36 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = s->size[1];
  emxEnsureCapacity_creal_T(y, i36);
  b4 = (y->size[1] == 0);
  if (!b4) {
    i36 = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity_creal_T(y, i36);
    loop_ub = y->size[1];
    for (i36 = 0; i36 < loop_ub; i36++) {
      y->data[i36].re = b[0];
      y->data[i36].im = 0.0;
    }

    for (k = 0; k < 28; k++) {
      i36 = s->size[0] * s->size[1];
      loop_ub = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = s->size[1];
      emxEnsureCapacity_creal_T(y, loop_ub);
      b_re = b[k + 1];
      loop_ub = i36 - 1;
      for (i36 = 0; i36 <= loop_ub; i36++) {
        s_re = s->data[i36].re * y->data[i36].re - s->data[i36].im * y->data[i36]
          .im;
        s_im = s->data[i36].re * y->data[i36].im + s->data[i36].im * y->data[i36]
          .re;
        y->data[i36].re = s_re + b_re;
        y->data[i36].im = s_im;
      }
    }
  }

  i36 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = digw->size[1];
  emxEnsureCapacity_creal_T(s, i36);
  loop_ub = digw->size[0] * digw->size[1];
  for (i36 = 0; i36 < loop_ub; i36++) {
    b_re = digw->data[i36] * 0.0;
    s_re = digw->data[i36];
    s->data[i36].re = 28.0 * b_re;
    s->data[i36].im = 28.0 * s_re;
  }

  emxFree_real_T(&digw);
  c_exp(s);
  rdivide_helper(y, s, hh);

  /*  Generate the default structure to pass to freqzplot */
  /*  If rad/sample, Fs is empty */
  emxFree_creal_T(&y);
  emxFree_creal_T(&s);
}

/*
 * FREQZ_CG Frequency response of digital filter with codegen support
 *
 *  This function is based on 'freqz' by The MathWorks Inc.
 * Arguments    : const double b[13]
 *                const emxArray_real_T *w
 *                double Fs
 *                emxArray_creal_T *hh
 * Return Type  : void
 */
static void p_freqz_cg(const double b[13], const emxArray_real_T *w, double Fs,
  emxArray_creal_T *hh)
{
  emxArray_real_T *digw;
  int i37;
  int loop_ub;
  emxArray_creal_T *s;
  emxArray_creal_T *y;
  bool b5;
  int k;
  double b_re;
  double s_re;
  double s_im;
  emxInit_real_T(&digw, 2);

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
  i37 = digw->size[0] * digw->size[1];
  digw->size[0] = 1;
  digw->size[1] = w->size[1];
  emxEnsureCapacity_real_T(digw, i37);
  loop_ub = w->size[0] * w->size[1];
  for (i37 = 0; i37 < loop_ub; i37++) {
    digw->data[i37] = 6.2831853071795862 * w->data[i37] / Fs;
  }

  emxInit_creal_T(&s, 2);

  /*  Convert from Hz to rad/sample for computational purposes */
  i37 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = digw->size[1];
  emxEnsureCapacity_creal_T(s, i37);
  loop_ub = digw->size[0] * digw->size[1];
  for (i37 = 0; i37 < loop_ub; i37++) {
    s->data[i37].re = digw->data[i37] * 0.0;
    s->data[i37].im = digw->data[i37];
  }

  emxInit_creal_T(&y, 2);
  c_exp(s);

  /*  Digital frequency must be used for this calculation */
  i37 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = s->size[1];
  emxEnsureCapacity_creal_T(y, i37);
  b5 = (y->size[1] == 0);
  if (!b5) {
    i37 = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity_creal_T(y, i37);
    loop_ub = y->size[1];
    for (i37 = 0; i37 < loop_ub; i37++) {
      y->data[i37].re = b[0];
      y->data[i37].im = 0.0;
    }

    for (k = 0; k < 12; k++) {
      i37 = s->size[0] * s->size[1];
      loop_ub = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = s->size[1];
      emxEnsureCapacity_creal_T(y, loop_ub);
      b_re = b[k + 1];
      loop_ub = i37 - 1;
      for (i37 = 0; i37 <= loop_ub; i37++) {
        s_re = s->data[i37].re * y->data[i37].re - s->data[i37].im * y->data[i37]
          .im;
        s_im = s->data[i37].re * y->data[i37].im + s->data[i37].im * y->data[i37]
          .re;
        y->data[i37].re = s_re + b_re;
        y->data[i37].im = s_im;
      }
    }
  }

  i37 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = digw->size[1];
  emxEnsureCapacity_creal_T(s, i37);
  loop_ub = digw->size[0] * digw->size[1];
  for (i37 = 0; i37 < loop_ub; i37++) {
    b_re = digw->data[i37] * 0.0;
    s_re = digw->data[i37];
    s->data[i37].re = 12.0 * b_re;
    s->data[i37].im = 12.0 * s_re;
  }

  emxFree_real_T(&digw);
  c_exp(s);
  rdivide_helper(y, s, hh);

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
  int i;
  bool p;
  int info;
  bool exitg2;
  emxArray_creal_T *beta1;
  emxArray_creal_T *T;
  unsigned int unnamed_idx_0;
  unsigned int unnamed_idx_1;
  int exitg1;
  double x_re;
  double alpha1_re;
  double x_im;
  double alpha1_im;
  bool b_x;
  double beta1_re;
  double beta1_im;
  int m;
  double brm;
  int istart;
  int jend;
  emxInit_creal_T(&alpha1, 1);
  if ((x->size[0] == 0) || (x->size[1] == 0)) {
    i = alpha1->size[0];
    alpha1->size[0] = x->size[0];
    emxEnsureCapacity_creal_T(alpha1, i);
    info = x->size[0];
    for (i = 0; i < info; i++) {
      alpha1->data[i].re = 0.0;
      alpha1->data[i].im = 0.0;
    }
  } else if (anyNonFinite(x)) {
    if ((x->size[0] == 1) && (x->size[1] == 1)) {
      i = alpha1->size[0];
      alpha1->size[0] = 1;
      emxEnsureCapacity_creal_T(alpha1, i);
      alpha1->data[0].re = rtNaN;
      alpha1->data[0].im = 0.0;
    } else {
      i = alpha1->size[0];
      alpha1->size[0] = x->size[0];
      emxEnsureCapacity_creal_T(alpha1, i);
      info = x->size[0];
      for (i = 0; i < info; i++) {
        alpha1->data[i].re = rtNaN;
        alpha1->data[i].im = 0.0;
      }
    }
  } else if ((x->size[0] == 1) && (x->size[1] == 1)) {
    i = alpha1->size[0];
    alpha1->size[0] = 1;
    emxEnsureCapacity_creal_T(alpha1, i);
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
      emxInit_creal_T(&T, 2);
      if (anyNonFinite(x)) {
        unnamed_idx_0 = (unsigned int)x->size[0];
        unnamed_idx_1 = (unsigned int)x->size[1];
        i = T->size[0] * T->size[1];
        T->size[0] = (int)unnamed_idx_0;
        T->size[1] = (int)unnamed_idx_1;
        emxEnsureCapacity_creal_T(T, i);
        info = (int)unnamed_idx_0 * (int)unnamed_idx_1;
        for (i = 0; i < info; i++) {
          T->data[i].re = rtNaN;
          T->data[i].im = 0.0;
        }

        m = T->size[0];
        if (1 < T->size[0]) {
          istart = 2;
          if (T->size[0] - 2 < T->size[1] - 1) {
            jend = T->size[0] - 1;
          } else {
            jend = T->size[1];
          }

          for (info = 0; info < jend; info++) {
            for (i = istart; i <= m; i++) {
              T->data[(i + T->size[0] * info) - 1].re = 0.0;
              T->data[(i + T->size[0] * info) - 1].im = 0.0;
            }

            istart++;
          }
        }
      } else {
        i = T->size[0] * T->size[1];
        T->size[0] = x->size[0];
        T->size[1] = x->size[1];
        emxEnsureCapacity_creal_T(T, i);
        info = x->size[0] * x->size[1];
        for (i = 0; i < info; i++) {
          T->data[i] = x->data[i];
        }

        xgehrd(T);
        eml_zlahqr(T);
        m = T->size[0];
        if ((T->size[0] == 0) || (T->size[1] == 0) || (3 >= T->size[0])) {
        } else {
          istart = 4;
          if (T->size[0] - 4 < T->size[1] - 1) {
            jend = T->size[0] - 3;
          } else {
            jend = T->size[1];
          }

          for (info = 0; info < jend; info++) {
            for (i = istart; i <= m; i++) {
              T->data[(i + T->size[0] * info) - 1].re = 0.0;
              T->data[(i + T->size[0] * info) - 1].im = 0.0;
            }

            istart++;
          }
        }
      }

      m = T->size[0];
      i = alpha1->size[0];
      alpha1->size[0] = T->size[0];
      emxEnsureCapacity_creal_T(alpha1, i);
      for (info = 0; info < m; info++) {
        alpha1->data[info] = T->data[info + T->size[0] * info];
      }

      emxFree_creal_T(&T);
    } else {
      emxInit_creal_T(&beta1, 1);
      xzgeev(x, &info, alpha1, beta1);
      i = alpha1->size[0];
      emxEnsureCapacity_creal_T(alpha1, i);
      info = alpha1->size[0];
      for (i = 0; i < info; i++) {
        alpha1_re = alpha1->data[i].re;
        alpha1_im = alpha1->data[i].im;
        beta1_re = beta1->data[i].re;
        beta1_im = beta1->data[i].im;
        if (beta1_im == 0.0) {
          if (alpha1_im == 0.0) {
            alpha1->data[i].re = alpha1_re / beta1_re;
            alpha1->data[i].im = 0.0;
          } else if (alpha1_re == 0.0) {
            alpha1->data[i].re = 0.0;
            alpha1->data[i].im = alpha1_im / beta1_re;
          } else {
            alpha1->data[i].re = alpha1_re / beta1_re;
            alpha1->data[i].im = alpha1_im / beta1_re;
          }
        } else if (beta1_re == 0.0) {
          if (alpha1_re == 0.0) {
            alpha1->data[i].re = alpha1_im / beta1_im;
            alpha1->data[i].im = 0.0;
          } else if (alpha1_im == 0.0) {
            alpha1->data[i].re = 0.0;
            alpha1->data[i].im = -(alpha1_re / beta1_im);
          } else {
            alpha1->data[i].re = alpha1_im / beta1_im;
            alpha1->data[i].im = -(alpha1_re / beta1_im);
          }
        } else {
          brm = fabs(beta1_re);
          x_re = fabs(beta1_im);
          if (brm > x_re) {
            x_im = beta1_im / beta1_re;
            x_re = beta1_re + x_im * beta1_im;
            alpha1->data[i].re = (alpha1_re + x_im * alpha1_im) / x_re;
            alpha1->data[i].im = (alpha1_im - x_im * alpha1_re) / x_re;
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

            alpha1->data[i].re = (alpha1_re * x_im + alpha1_im * x_re) / brm;
            alpha1->data[i].im = (alpha1_im * x_im - alpha1_re * x_re) / brm;
          } else {
            x_im = beta1_re / beta1_im;
            x_re = beta1_im + x_im * beta1_re;
            alpha1->data[i].re = (x_im * alpha1_re + alpha1_im) / x_re;
            alpha1->data[i].im = (x_im * alpha1_im - alpha1_re) / x_re;
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
  int i9;
  int k;
  double p_re;
  double x_im;
  for (i9 = 0; i9 < 2048; i9++) {
    y[i9].re = p[0];
    y[i9].im = 0.0;
  }

  for (k = 0; k < 14; k++) {
    p_re = p[k + 1];
    for (i9 = 0; i9 < 2048; i9++) {
      x_im = x[i9].re * y[i9].im + x[i9].im * y[i9].re;
      y[i9].re = (x[i9].re * y[i9].re - x[i9].im * y[i9].im) + p_re;
      y[i9].im = x_im;
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
    y[k] = rt_powd_snf(a[k], 3.0);
  }
}

/*
 * FREQZ_CG Frequency response of digital filter with codegen support
 *
 *  This function is based on 'freqz' by The MathWorks Inc.
 * Arguments    : const double b[57]
 *                const emxArray_real_T *w
 *                double Fs
 *                emxArray_creal_T *hh
 * Return Type  : void
 */
static void q_freqz_cg(const double b[57], const emxArray_real_T *w, double Fs,
  emxArray_creal_T *hh)
{
  emxArray_real_T *digw;
  int i38;
  int loop_ub;
  emxArray_creal_T *s;
  emxArray_creal_T *y;
  bool b6;
  int k;
  double b_re;
  double s_re;
  double s_im;
  emxInit_real_T(&digw, 2);

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
  i38 = digw->size[0] * digw->size[1];
  digw->size[0] = 1;
  digw->size[1] = w->size[1];
  emxEnsureCapacity_real_T(digw, i38);
  loop_ub = w->size[0] * w->size[1];
  for (i38 = 0; i38 < loop_ub; i38++) {
    digw->data[i38] = 6.2831853071795862 * w->data[i38] / Fs;
  }

  emxInit_creal_T(&s, 2);

  /*  Convert from Hz to rad/sample for computational purposes */
  i38 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = digw->size[1];
  emxEnsureCapacity_creal_T(s, i38);
  loop_ub = digw->size[0] * digw->size[1];
  for (i38 = 0; i38 < loop_ub; i38++) {
    s->data[i38].re = digw->data[i38] * 0.0;
    s->data[i38].im = digw->data[i38];
  }

  emxInit_creal_T(&y, 2);
  c_exp(s);

  /*  Digital frequency must be used for this calculation */
  i38 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = s->size[1];
  emxEnsureCapacity_creal_T(y, i38);
  b6 = (y->size[1] == 0);
  if (!b6) {
    i38 = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity_creal_T(y, i38);
    loop_ub = y->size[1];
    for (i38 = 0; i38 < loop_ub; i38++) {
      y->data[i38].re = b[0];
      y->data[i38].im = 0.0;
    }

    for (k = 0; k < 56; k++) {
      i38 = s->size[0] * s->size[1];
      loop_ub = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = s->size[1];
      emxEnsureCapacity_creal_T(y, loop_ub);
      b_re = b[k + 1];
      loop_ub = i38 - 1;
      for (i38 = 0; i38 <= loop_ub; i38++) {
        s_re = s->data[i38].re * y->data[i38].re - s->data[i38].im * y->data[i38]
          .im;
        s_im = s->data[i38].re * y->data[i38].im + s->data[i38].im * y->data[i38]
          .re;
        y->data[i38].re = s_re + b_re;
        y->data[i38].im = s_im;
      }
    }
  }

  i38 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = digw->size[1];
  emxEnsureCapacity_creal_T(s, i38);
  loop_ub = digw->size[0] * digw->size[1];
  for (i38 = 0; i38 < loop_ub; i38++) {
    b_re = digw->data[i38] * 0.0;
    s_re = digw->data[i38];
    s->data[i38].re = 56.0 * b_re;
    s->data[i38].im = 56.0 * s_re;
  }

  emxFree_real_T(&digw);
  c_exp(s);
  rdivide_helper(y, s, hh);

  /*  Generate the default structure to pass to freqzplot */
  /*  If rad/sample, Fs is empty */
  emxFree_creal_T(&y);
  emxFree_creal_T(&s);
}

/*
 * FREQZ_CG Frequency response of digital filter with codegen support
 *
 *  This function is based on 'freqz' by The MathWorks Inc.
 * Arguments    : const double b[43]
 *                const emxArray_real_T *w
 *                double Fs
 *                emxArray_creal_T *hh
 * Return Type  : void
 */
static void r_freqz_cg(const double b[43], const emxArray_real_T *w, double Fs,
  emxArray_creal_T *hh)
{
  emxArray_real_T *digw;
  int i39;
  int loop_ub;
  emxArray_creal_T *s;
  emxArray_creal_T *y;
  bool b7;
  int k;
  double b_re;
  double s_re;
  double s_im;
  emxInit_real_T(&digw, 2);

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
  i39 = digw->size[0] * digw->size[1];
  digw->size[0] = 1;
  digw->size[1] = w->size[1];
  emxEnsureCapacity_real_T(digw, i39);
  loop_ub = w->size[0] * w->size[1];
  for (i39 = 0; i39 < loop_ub; i39++) {
    digw->data[i39] = 6.2831853071795862 * w->data[i39] / Fs;
  }

  emxInit_creal_T(&s, 2);

  /*  Convert from Hz to rad/sample for computational purposes */
  i39 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = digw->size[1];
  emxEnsureCapacity_creal_T(s, i39);
  loop_ub = digw->size[0] * digw->size[1];
  for (i39 = 0; i39 < loop_ub; i39++) {
    s->data[i39].re = digw->data[i39] * 0.0;
    s->data[i39].im = digw->data[i39];
  }

  emxInit_creal_T(&y, 2);
  c_exp(s);

  /*  Digital frequency must be used for this calculation */
  i39 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = s->size[1];
  emxEnsureCapacity_creal_T(y, i39);
  b7 = (y->size[1] == 0);
  if (!b7) {
    i39 = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity_creal_T(y, i39);
    loop_ub = y->size[1];
    for (i39 = 0; i39 < loop_ub; i39++) {
      y->data[i39].re = b[0];
      y->data[i39].im = 0.0;
    }

    for (k = 0; k < 42; k++) {
      i39 = s->size[0] * s->size[1];
      loop_ub = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = s->size[1];
      emxEnsureCapacity_creal_T(y, loop_ub);
      b_re = b[k + 1];
      loop_ub = i39 - 1;
      for (i39 = 0; i39 <= loop_ub; i39++) {
        s_re = s->data[i39].re * y->data[i39].re - s->data[i39].im * y->data[i39]
          .im;
        s_im = s->data[i39].re * y->data[i39].im + s->data[i39].im * y->data[i39]
          .re;
        y->data[i39].re = s_re + b_re;
        y->data[i39].im = s_im;
      }
    }
  }

  i39 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = digw->size[1];
  emxEnsureCapacity_creal_T(s, i39);
  loop_ub = digw->size[0] * digw->size[1];
  for (i39 = 0; i39 < loop_ub; i39++) {
    b_re = digw->data[i39] * 0.0;
    s_re = digw->data[i39];
    s->data[i39].re = 42.0 * b_re;
    s->data[i39].im = 42.0 * s_re;
  }

  emxFree_real_T(&digw);
  c_exp(s);
  rdivide_helper(y, s, hh);

  /*  Generate the default structure to pass to freqzplot */
  /*  If rad/sample, Fs is empty */
  emxFree_creal_T(&y);
  emxFree_creal_T(&s);
}

/*
 * Arguments    : const emxArray_creal_T *x
 *                const emxArray_creal_T *y
 *                emxArray_creal_T *z
 * Return Type  : void
 */
static void rdivide_helper(const emxArray_creal_T *x, const emxArray_creal_T *y,
  emxArray_creal_T *z)
{
  int i30;
  int loop_ub;
  double x_re;
  double x_im;
  double y_re;
  double y_im;
  double brm;
  double bim;
  double s;
  i30 = z->size[0] * z->size[1];
  z->size[0] = 1;
  z->size[1] = x->size[1];
  emxEnsureCapacity_creal_T(z, i30);
  loop_ub = x->size[0] * x->size[1];
  for (i30 = 0; i30 < loop_ub; i30++) {
    x_re = x->data[i30].re;
    x_im = x->data[i30].im;
    y_re = y->data[i30].re;
    y_im = y->data[i30].im;
    if (y_im == 0.0) {
      if (x_im == 0.0) {
        z->data[i30].re = x_re / y_re;
        z->data[i30].im = 0.0;
      } else if (x_re == 0.0) {
        z->data[i30].re = 0.0;
        z->data[i30].im = x_im / y_re;
      } else {
        z->data[i30].re = x_re / y_re;
        z->data[i30].im = x_im / y_re;
      }
    } else if (y_re == 0.0) {
      if (x_re == 0.0) {
        z->data[i30].re = x_im / y_im;
        z->data[i30].im = 0.0;
      } else if (x_im == 0.0) {
        z->data[i30].re = 0.0;
        z->data[i30].im = -(x_re / y_im);
      } else {
        z->data[i30].re = x_im / y_im;
        z->data[i30].im = -(x_re / y_im);
      }
    } else {
      brm = fabs(y_re);
      bim = fabs(y_im);
      if (brm > bim) {
        s = y_im / y_re;
        bim = y_re + s * y_im;
        z->data[i30].re = (x_re + s * x_im) / bim;
        z->data[i30].im = (x_im - s * x_re) / bim;
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

        z->data[i30].re = (x_re * s + x_im * bim) / brm;
        z->data[i30].im = (x_im * s - x_re * bim) / brm;
      } else {
        s = y_re / y_im;
        bim = y_im + s * y_re;
        z->data[i30].re = (s * x_re + x_im) / bim;
        z->data[i30].im = (s * x_im - x_re) / bim;
      }
    }
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
  int i55;
  emxArray_real_T *xx;
  emxArray_int32_T *r12;
  int l;
  int i56;
  int vlen;
  int i;
  double b_x;
  int end;
  int loop_ub;

  /*  */
  /*    Author: T. Krauss 1993 */
  /*        Was Revision: 1.4, Date: 1994/01/25 17:59:44 */
  y = 1.0;
  i55 = (int)m;
  emxInit_real_T(&xx, 2);
  emxInit_int32_T(&r12, 2);
  for (l = 0; l < i55; l++) {
    if ((m == 0.0) || (((m > 0.0) && (1.0 + (double)l > n)) || ((0.0 > m) && (n >
           1.0 + (double)l)))) {
      i56 = 1;
      vlen = 1;
      i = 0;
    } else {
      i56 = (int)(1U + l);
      vlen = i55;
      i = (int)n;
    }

    b_x = x->data[(int)k - 1];
    end = xx->size[0] * xx->size[1];
    xx->size[0] = 1;
    loop_ub = div_s32_floor(i - i56, vlen);
    xx->size[1] = loop_ub + 1;
    emxEnsureCapacity_real_T(xx, end);
    for (i = 0; i <= loop_ub; i++) {
      xx->data[i] = 2.0 * (b_x - x->data[(i56 + vlen * i) - 1]);
    }

    end = xx->size[1] - 1;
    vlen = 0;
    for (i = 0; i <= end; i++) {
      if (xx->data[i] != 0.0) {
        vlen++;
      }
    }

    i56 = r12->size[0] * r12->size[1];
    r12->size[0] = 1;
    r12->size[1] = vlen;
    emxEnsureCapacity_int32_T(r12, i56);
    vlen = 0;
    for (i = 0; i <= end; i++) {
      if (xx->data[i] != 0.0) {
        r12->data[vlen] = i + 1;
        vlen++;
      }
    }

    vlen = r12->size[1];
    if (r12->size[1] == 0) {
      b_x = 1.0;
    } else {
      b_x = xx->data[r12->data[0] - 1];
      for (i = 2; i <= vlen; i++) {
        b_x *= xx->data[r12->data[i - 1] - 1];
      }
    }

    y *= b_x;
  }

  emxFree_int32_T(&r12);
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
 *                bool *valid
 * Return Type  : void
 */
static void remezm(double nfilt, const double edge[4], const emxArray_real_T
                   *grid, emxArray_real_T *des, emxArray_real_T *wt,
                   emxArray_real_T *h, double *dev, bool *valid)
{
  double nodd;
  double nfcns;
  int varargin_2;
  int ngrid;
  emxArray_real_T *j;
  emxArray_real_T *sel;
  int i49;
  double ft;
  int loop_ub;
  emxArray_real_T *iext;
  emxArray_real_T *y;
  emxArray_real_T *x;
  int i50;
  double comp;
  double dtemp;
  double b_y1;
  int luck;
  int nut1;
  double err;
  double devl;
  emxArray_real_T *ad;
  int niter;
  double jchnge;
  double b_x;
  emxArray_real_T *l;
  emxArray_int32_T *r11;
  emxArray_real_T *b;
  emxArray_real_T *a;
  emxArray_real_T *b_wt;
  emxArray_int32_T *b_iext;
  bool guard1 = false;
  int exitg1;
  int b_loop_ub;
  int nut;
  double k1;
  double knz;
  double b_l;
  int kkk;
  int nu;
  double dnum;
  double dden;
  int i51;
  int i52;
  int l_tmp;
  int i53;
  bool guard2 = false;
  double b_j;
  int flag34;
  int exitg3;
  int flag;
  bool exitg2;
  int x_tmp;
  bool exitg4;

  /*  */
  *valid = true;
  nodd = rt_remd_snf(nfilt, 2.0);

  /*  nodd == 1 ==> filter length is odd */
  /*  nodd == 0 ==> filter length is even */
  nfcns = trunc(nfilt / 2.0);
  if (nodd == 1.0) {
    nfcns++;
  }

  varargin_2 = grid->size[1];
  ngrid = grid->size[1];
  emxInit_real_T(&j, 2);
  emxInit_real_T(&sel, 2);
  if (nodd != 1.0) {
    i49 = sel->size[0] * sel->size[1];
    sel->size[0] = 1;
    sel->size[1] = grid->size[1];
    emxEnsureCapacity_real_T(sel, i49);
    loop_ub = grid->size[0] * grid->size[1];
    for (i49 = 0; i49 < loop_ub; i49++) {
      sel->data[i49] = 3.1415926535897931 * grid->data[i49];
    }

    b_cos(sel);
    b_rdivide_helper(des, sel, j);
    i49 = des->size[0] * des->size[1];
    des->size[0] = 1;
    des->size[1] = j->size[1];
    emxEnsureCapacity_real_T(des, i49);
    loop_ub = j->size[0] * j->size[1];
    for (i49 = 0; i49 < loop_ub; i49++) {
      des->data[i49] = j->data[i49];
    }

    i49 = sel->size[0] * sel->size[1];
    sel->size[0] = 1;
    sel->size[1] = grid->size[1];
    emxEnsureCapacity_real_T(sel, i49);
    loop_ub = grid->size[0] * grid->size[1];
    for (i49 = 0; i49 < loop_ub; i49++) {
      sel->data[i49] = 3.1415926535897931 * grid->data[i49];
    }

    b_cos(sel);
    i49 = wt->size[0] * wt->size[1];
    i50 = wt->size[0] * wt->size[1];
    wt->size[0] = 1;
    emxEnsureCapacity_real_T(wt, i50);
    loop_ub = i49 - 1;
    for (i49 = 0; i49 <= loop_ub; i49++) {
      wt->data[i49] *= sel->data[i49];
    }
  }

  ft = ((double)grid->size[1] - 1.0) / nfcns;
  if (rtIsNaN(nfcns)) {
    i49 = j->size[0] * j->size[1];
    j->size[0] = 1;
    j->size[1] = 1;
    emxEnsureCapacity_real_T(j, i49);
    j->data[0] = rtNaN;
  } else if (nfcns < 1.0) {
    j->size[0] = 1;
    j->size[1] = 0;
  } else if (rtIsInf(nfcns) && (1.0 == nfcns)) {
    i49 = j->size[0] * j->size[1];
    j->size[0] = 1;
    j->size[1] = 1;
    emxEnsureCapacity_real_T(j, i49);
    j->data[0] = rtNaN;
  } else {
    i49 = j->size[0] * j->size[1];
    j->size[0] = 1;
    loop_ub = (int)(nfcns - 1.0);
    j->size[1] = loop_ub + 1;
    emxEnsureCapacity_real_T(j, i49);
    for (i49 = 0; i49 <= loop_ub; i49++) {
      j->data[i49] = 1.0 + (double)i49;
    }
  }

  i49 = sel->size[0] * sel->size[1];
  sel->size[0] = 1;
  sel->size[1] = j->size[1] + 1;
  emxEnsureCapacity_real_T(sel, i49);
  loop_ub = j->size[1];
  for (i49 = 0; i49 < loop_ub; i49++) {
    sel->data[i49] = ft * (j->data[i49] - 1.0) + 1.0;
  }

  emxInit_real_T(&iext, 1);
  sel->data[j->size[1]] = grid->size[1];
  b_fix(sel);
  i49 = iext->size[0];
  iext->size[0] = sel->size[1] + 1;
  emxEnsureCapacity_real_T(iext, i49);
  loop_ub = sel->size[1];
  for (i49 = 0; i49 < loop_ub; i49++) {
    iext->data[i49] = sel->data[i49];
  }

  emxInit_real_T(&y, 2);
  emxInit_real_T(&x, 2);
  iext->data[sel->size[1]] = 0.0;

  /*  Remez exchange loop */
  comp = -1.0;
  dtemp = -1.0;
  b_y1 = -1.0;
  luck = -1;
  nut1 = -1;
  err = -1.0;
  i49 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = 1;
  emxEnsureCapacity_real_T(y, i49);
  y->data[0] = -1.0;
  *dev = -1.0;
  devl = -1.0;
  i49 = x->size[0] * x->size[1];
  x->size[0] = 1;
  loop_ub = (int)(nfcns + 1.0);
  x->size[1] = loop_ub;
  emxEnsureCapacity_real_T(x, i49);
  for (i49 = 0; i49 < loop_ub; i49++) {
    x->data[i49] = 0.0;
  }

  emxInit_real_T(&ad, 2);
  niter = 0;
  jchnge = 1.0;
  b_x = trunc((nfcns - 1.0) / 15.0);
  i49 = ad->size[0] * ad->size[1];
  ad->size[0] = 1;
  ad->size[1] = loop_ub;
  emxEnsureCapacity_real_T(ad, i49);
  for (i49 = 0; i49 < loop_ub; i49++) {
    ad->data[i49] = 0.0;
  }

  /*  index manager(s) */
  emxInit_real_T(&l, 2);
  emxInit_int32_T(&r11, 2);
  emxInit_real_T(&b, 1);
  emxInit_real_T(&a, 2);
  emxInit_real_T(&b_wt, 2);
  emxInit_int32_T(&b_iext, 1);
  guard1 = false;
  do {
    exitg1 = 0;
    if (jchnge > 0.0) {
      i49 = (int)((nfcns + 1.0) + 1.0) - 1;
      iext->data[i49] = (double)varargin_2 + 1.0;
      niter++;
      if (niter > 250) {
        guard1 = true;
        exitg1 = 1;
      } else {
        if (1.0 > nfcns + 1.0) {
          b_loop_ub = 0;
        } else {
          b_loop_ub = (int)(nfcns + 1.0);
        }

        i50 = l->size[0] * l->size[1];
        l->size[0] = 1;
        l->size[1] = b_loop_ub;
        emxEnsureCapacity_real_T(l, i50);
        for (i50 = 0; i50 < b_loop_ub; i50++) {
          l->data[i50] = iext->data[i50];
        }

        i50 = x->size[0] * x->size[1];
        x->size[0] = 1;
        x->size[1] = l->size[1];
        emxEnsureCapacity_real_T(x, i50);
        nut = l->size[0] * l->size[1];
        for (i50 = 0; i50 < nut; i50++) {
          x->data[i50] = 6.2831853071795862 * grid->data[(int)l->data[i50] - 1];
        }

        b_cos(x);
        for (nu = 0; nu < loop_ub; nu++) {
          ad->data[nu] = remezdd(1.0 + (double)nu, nfcns + 1.0, b_x + 1.0, x);
        }

        ft = ad->size[1];
        i50 = j->size[0] * j->size[1];
        j->size[0] = 1;
        j->size[1] = (int)ft;
        emxEnsureCapacity_real_T(j, i50);
        nut = (int)ft;
        for (i50 = 0; i50 < nut; i50++) {
          j->data[i50] = 1.0;
        }

        if (2.0 > nfcns + 1.0) {
          i50 = 1;
          i51 = 1;
          i52 = 0;
          l_tmp = 0;
          i53 = 1;
        } else {
          i50 = 2;
          i51 = 2;
          i52 = loop_ub;
          l_tmp = 1;
          i53 = 2;
        }

        nut = div_s32_floor(i52 - i50, i51);
        for (i50 = 0; i50 <= nut; i50++) {
          j->data[l_tmp + i53 * i50] = -1.0;
        }

        i50 = b->size[0];
        b->size[0] = l->size[1];
        emxEnsureCapacity_real_T(b, i50);
        nut = l->size[1];
        for (i50 = 0; i50 < nut; i50++) {
          b->data[i50] = des->data[(int)l->data[i50] - 1];
        }

        guard2 = false;
        if (ad->size[1] == 1) {
          guard2 = true;
        } else {
          i50 = b_iext->size[0];
          b_iext->size[0] = b_loop_ub;
          emxEnsureCapacity_int32_T(b_iext, i50);
          for (i50 = 0; i50 < b_loop_ub; i50++) {
            b_iext->data[i50] = (int)iext->data[i50];
          }

          if (b_iext->size[0] == 1) {
            guard2 = true;
          } else {
            dnum = 0.0;
            b_loop_ub = ad->size[1];
            for (i50 = 0; i50 < b_loop_ub; i50++) {
              dnum += ad->data[i50] * b->data[i50];
            }
          }
        }

        if (guard2) {
          dnum = 0.0;
          b_loop_ub = ad->size[1];
          for (i50 = 0; i50 < b_loop_ub; i50++) {
            dnum += ad->data[i50] * b->data[i50];
          }
        }

        i50 = b_wt->size[0] * b_wt->size[1];
        b_wt->size[0] = 1;
        b_wt->size[1] = l->size[1];
        emxEnsureCapacity_real_T(b_wt, i50);
        b_loop_ub = l->size[0] * l->size[1];
        for (i50 = 0; i50 < b_loop_ub; i50++) {
          b_wt->data[i50] = wt->data[(int)l->data[i50] - 1];
        }

        b_rdivide_helper(ad, b_wt, sel);
        i50 = b->size[0];
        b->size[0] = sel->size[1];
        emxEnsureCapacity_real_T(b, i50);
        b_loop_ub = sel->size[1];
        for (i50 = 0; i50 < b_loop_ub; i50++) {
          b->data[i50] = sel->data[i50];
        }

        if ((j->size[1] == 1) || (b->size[0] == 1)) {
          dden = 0.0;
          b_loop_ub = j->size[1];
          for (i50 = 0; i50 < b_loop_ub; i50++) {
            dden += j->data[i50] * b->data[i50];
          }
        } else {
          dden = 0.0;
          b_loop_ub = j->size[1];
          for (i50 = 0; i50 < b_loop_ub; i50++) {
            dden += j->data[i50] * b->data[i50];
          }
        }

        *dev = dnum / dden;
        nu = 1;
        if (*dev > 0.0) {
          nu = -1;
        }

        *dev *= -(double)nu;
        ft = (double)nu * *dev;
        i50 = a->size[0] * a->size[1];
        a->size[0] = 1;
        a->size[1] = j->size[1];
        emxEnsureCapacity_real_T(a, i50);
        b_loop_ub = j->size[0] * j->size[1];
        for (i50 = 0; i50 < b_loop_ub; i50++) {
          a->data[i50] = ft * j->data[i50];
        }

        i50 = b_wt->size[0] * b_wt->size[1];
        b_wt->size[0] = 1;
        b_wt->size[1] = l->size[1];
        emxEnsureCapacity_real_T(b_wt, i50);
        b_loop_ub = l->size[0] * l->size[1];
        for (i50 = 0; i50 < b_loop_ub; i50++) {
          b_wt->data[i50] = wt->data[(int)l->data[i50] - 1];
        }

        b_rdivide_helper(a, b_wt, sel);
        i50 = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = l->size[1];
        emxEnsureCapacity_real_T(y, i50);
        b_loop_ub = l->size[0] * l->size[1];
        for (i50 = 0; i50 < b_loop_ub; i50++) {
          y->data[i50] = des->data[(int)l->data[i50] - 1] + sel->data[i50];
        }

        if (*dev <= devl) {
          /* warning(message('signal:firpm:DidNotConverge',niter)) */
          printf("%s\n", "DidNotConverge");
          fflush(stdout);
          i49 = h->size[0] * h->size[1];
          loop_ub = (int)nfilt;
          h->size[0] = loop_ub;
          h->size[1] = 1;
          emxEnsureCapacity_real_T(h, i49);
          for (i49 = 0; i49 < loop_ub; i49++) {
            h->data[i49] = 0.0;
          }

          *dev = -1.0;

          /* iext */
          *valid = false;
          exitg1 = 1;
        } else {
          devl = *dev;
          jchnge = 0.0;
          k1 = iext->data[0];
          kkk = loop_ub - 1;
          knz = iext->data[kkk];
          dden = 0.0;
          nut = -nu;
          b_j = 1.0;
          flag34 = 1;
          while (b_j < (nfcns + 1.0) + 1.0) {
            dnum = iext->data[(int)(unsigned int)b_j];
            l_tmp = (int)b_j - 1;
            b_l = iext->data[l_tmp] + 1.0;
            nut = -nut;
            if (b_j == 2.0) {
              b_y1 = comp;
            }

            comp = *dev;
            flag = 1;
            if (iext->data[l_tmp] + 1.0 < iext->data[(int)b_j]) {
              /*  gee */
              ft = cos(6.2831853071795862 * grid->data[(int)(iext->data[l_tmp] +
                        1.0) - 1]);
              i50 = b_wt->size[0] * b_wt->size[1];
              b_wt->size[0] = 1;
              b_wt->size[1] = x->size[1];
              emxEnsureCapacity_real_T(b_wt, i50);
              b_loop_ub = x->size[0] * x->size[1];
              for (i50 = 0; i50 < b_loop_ub; i50++) {
                b_wt->data[i50] = ft - x->data[i50];
              }

              b_rdivide_helper(ad, b_wt, j);
              i50 = b->size[0];
              b->size[0] = y->size[1];
              emxEnsureCapacity_real_T(b, i50);
              b_loop_ub = y->size[1];
              for (i50 = 0; i50 < b_loop_ub; i50++) {
                b->data[i50] = y->data[i50];
              }

              if ((j->size[1] == 1) || (b->size[0] == 1)) {
                ft = 0.0;
                b_loop_ub = j->size[1];
                for (i50 = 0; i50 < b_loop_ub; i50++) {
                  ft += j->data[i50] * b->data[i50];
                }
              } else {
                ft = 0.0;
                b_loop_ub = j->size[1];
                for (i50 = 0; i50 < b_loop_ub; i50++) {
                  ft += j->data[i50] * b->data[i50];
                }
              }

              err = (ft / b_sum(j) - des->data[(int)(iext->data[l_tmp] + 1.0) -
                     1]) * wt->data[(int)(iext->data[l_tmp] + 1.0) - 1];
              ft = (double)nut * err;
              dtemp = ft - *dev;
              if (dtemp > 0.0) {
                comp = ft;
                b_l = (iext->data[l_tmp] + 1.0) + 1.0;
                exitg2 = false;
                while ((!exitg2) && (b_l < dnum)) {
                  /*  gee */
                  x_tmp = (int)b_l - 1;
                  ft = cos(6.2831853071795862 * grid->data[x_tmp]);
                  i50 = b_wt->size[0] * b_wt->size[1];
                  b_wt->size[0] = 1;
                  b_wt->size[1] = x->size[1];
                  emxEnsureCapacity_real_T(b_wt, i50);
                  b_loop_ub = x->size[0] * x->size[1];
                  for (i50 = 0; i50 < b_loop_ub; i50++) {
                    b_wt->data[i50] = ft - x->data[i50];
                  }

                  b_rdivide_helper(ad, b_wt, j);
                  i50 = b->size[0];
                  b->size[0] = y->size[1];
                  emxEnsureCapacity_real_T(b, i50);
                  b_loop_ub = y->size[1];
                  for (i50 = 0; i50 < b_loop_ub; i50++) {
                    b->data[i50] = y->data[i50];
                  }

                  if ((j->size[1] == 1) || (b->size[0] == 1)) {
                    ft = 0.0;
                    b_loop_ub = j->size[1];
                    for (i50 = 0; i50 < b_loop_ub; i50++) {
                      ft += j->data[i50] * b->data[i50];
                    }
                  } else {
                    ft = 0.0;
                    b_loop_ub = j->size[1];
                    for (i50 = 0; i50 < b_loop_ub; i50++) {
                      ft += j->data[i50] * b->data[i50];
                    }
                  }

                  err = (ft / b_sum(j) - des->data[x_tmp]) * wt->data[x_tmp];
                  dtemp = (double)nut * err - comp;
                  if (dtemp > 0.0) {
                    comp = (double)nut * err;
                    b_l++;
                  } else {
                    exitg2 = true;
                  }
                }

                iext->data[l_tmp] = b_l - 1.0;
                b_j++;
                dden = b_l - 1.0;
                jchnge++;
                flag = 0;
              }
            }

            if (flag != 0) {
              b_l -= 2.0;
              exitg2 = false;
              while ((!exitg2) && (b_l > dden)) {
                /*  gee */
                x_tmp = (int)b_l - 1;
                ft = cos(6.2831853071795862 * grid->data[x_tmp]);
                i50 = b_wt->size[0] * b_wt->size[1];
                b_wt->size[0] = 1;
                b_wt->size[1] = x->size[1];
                emxEnsureCapacity_real_T(b_wt, i50);
                b_loop_ub = x->size[0] * x->size[1];
                for (i50 = 0; i50 < b_loop_ub; i50++) {
                  b_wt->data[i50] = ft - x->data[i50];
                }

                b_rdivide_helper(ad, b_wt, j);
                i50 = b->size[0];
                b->size[0] = y->size[1];
                emxEnsureCapacity_real_T(b, i50);
                b_loop_ub = y->size[1];
                for (i50 = 0; i50 < b_loop_ub; i50++) {
                  b->data[i50] = y->data[i50];
                }

                if ((j->size[1] == 1) || (b->size[0] == 1)) {
                  ft = 0.0;
                  b_loop_ub = j->size[1];
                  for (i50 = 0; i50 < b_loop_ub; i50++) {
                    ft += j->data[i50] * b->data[i50];
                  }
                } else {
                  ft = 0.0;
                  b_loop_ub = j->size[1];
                  for (i50 = 0; i50 < b_loop_ub; i50++) {
                    ft += j->data[i50] * b->data[i50];
                  }
                }

                err = (ft / b_sum(j) - des->data[x_tmp]) * wt->data[x_tmp];
                dtemp = (double)nut * err - comp;
                if ((dtemp > 0.0) || (jchnge > 0.0)) {
                  exitg2 = true;
                } else {
                  b_l--;
                }
              }

              if (b_l <= dden) {
                b_l = iext->data[(int)b_j - 1] + 1.0;
                if (jchnge > 0.0) {
                  iext->data[(int)b_j - 1] = (iext->data[(int)b_j - 1] + 1.0) -
                    1.0;
                  b_j++;
                  dden = b_l - 1.0;
                  jchnge++;
                } else {
                  b_l = (iext->data[(int)b_j - 1] + 1.0) + 1.0;
                  exitg2 = false;
                  while ((!exitg2) && (b_l < dnum)) {
                    /*  gee */
                    x_tmp = (int)b_l - 1;
                    ft = cos(6.2831853071795862 * grid->data[x_tmp]);
                    i50 = b_wt->size[0] * b_wt->size[1];
                    b_wt->size[0] = 1;
                    b_wt->size[1] = x->size[1];
                    emxEnsureCapacity_real_T(b_wt, i50);
                    b_loop_ub = x->size[0] * x->size[1];
                    for (i50 = 0; i50 < b_loop_ub; i50++) {
                      b_wt->data[i50] = ft - x->data[i50];
                    }

                    b_rdivide_helper(ad, b_wt, j);
                    i50 = b->size[0];
                    b->size[0] = y->size[1];
                    emxEnsureCapacity_real_T(b, i50);
                    b_loop_ub = y->size[1];
                    for (i50 = 0; i50 < b_loop_ub; i50++) {
                      b->data[i50] = y->data[i50];
                    }

                    if ((j->size[1] == 1) || (b->size[0] == 1)) {
                      ft = 0.0;
                      b_loop_ub = j->size[1];
                      for (i50 = 0; i50 < b_loop_ub; i50++) {
                        ft += j->data[i50] * b->data[i50];
                      }
                    } else {
                      ft = 0.0;
                      b_loop_ub = j->size[1];
                      for (i50 = 0; i50 < b_loop_ub; i50++) {
                        ft += j->data[i50] * b->data[i50];
                      }
                    }

                    err = (ft / b_sum(j) - des->data[x_tmp]) * wt->data[x_tmp];
                    dtemp = (double)nut * err - comp;
                    if (dtemp > 0.0) {
                      exitg2 = true;
                    } else {
                      b_l++;
                    }
                  }

                  if ((b_l < dnum) && (dtemp > 0.0)) {
                    comp = (double)nut * err;
                    b_l++;
                    exitg2 = false;
                    while ((!exitg2) && (b_l < dnum)) {
                      /*  gee */
                      x_tmp = (int)b_l - 1;
                      ft = cos(6.2831853071795862 * grid->data[x_tmp]);
                      i50 = b_wt->size[0] * b_wt->size[1];
                      b_wt->size[0] = 1;
                      b_wt->size[1] = x->size[1];
                      emxEnsureCapacity_real_T(b_wt, i50);
                      b_loop_ub = x->size[0] * x->size[1];
                      for (i50 = 0; i50 < b_loop_ub; i50++) {
                        b_wt->data[i50] = ft - x->data[i50];
                      }

                      b_rdivide_helper(ad, b_wt, j);
                      i50 = b->size[0];
                      b->size[0] = y->size[1];
                      emxEnsureCapacity_real_T(b, i50);
                      b_loop_ub = y->size[1];
                      for (i50 = 0; i50 < b_loop_ub; i50++) {
                        b->data[i50] = y->data[i50];
                      }

                      if ((j->size[1] == 1) || (b->size[0] == 1)) {
                        ft = 0.0;
                        b_loop_ub = j->size[1];
                        for (i50 = 0; i50 < b_loop_ub; i50++) {
                          ft += j->data[i50] * b->data[i50];
                        }
                      } else {
                        ft = 0.0;
                        b_loop_ub = j->size[1];
                        for (i50 = 0; i50 < b_loop_ub; i50++) {
                          ft += j->data[i50] * b->data[i50];
                        }
                      }

                      err = (ft / b_sum(j) - des->data[x_tmp]) * wt->data[x_tmp];
                      ft = (double)nut * err;
                      dtemp = ft - comp;
                      if (dtemp > 0.0) {
                        comp = ft;
                        b_l++;
                      } else {
                        exitg2 = true;
                      }
                    }

                    iext->data[(int)b_j - 1] = b_l - 1.0;
                    b_j++;
                    dden = b_l - 1.0;
                    jchnge = 1.0;
                  } else {
                    dden = iext->data[(int)b_j - 1];
                    b_j++;
                  }
                }
              } else if (dtemp > 0.0) {
                comp = (double)nut * err;
                b_l--;
                exitg2 = false;
                while ((!exitg2) && (b_l > dden)) {
                  /*  gee */
                  x_tmp = (int)b_l - 1;
                  ft = cos(6.2831853071795862 * grid->data[x_tmp]);
                  i50 = b_wt->size[0] * b_wt->size[1];
                  b_wt->size[0] = 1;
                  b_wt->size[1] = x->size[1];
                  emxEnsureCapacity_real_T(b_wt, i50);
                  b_loop_ub = x->size[0] * x->size[1];
                  for (i50 = 0; i50 < b_loop_ub; i50++) {
                    b_wt->data[i50] = ft - x->data[i50];
                  }

                  b_rdivide_helper(ad, b_wt, j);
                  i50 = b->size[0];
                  b->size[0] = y->size[1];
                  emxEnsureCapacity_real_T(b, i50);
                  b_loop_ub = y->size[1];
                  for (i50 = 0; i50 < b_loop_ub; i50++) {
                    b->data[i50] = y->data[i50];
                  }

                  if ((j->size[1] == 1) || (b->size[0] == 1)) {
                    ft = 0.0;
                    b_loop_ub = j->size[1];
                    for (i50 = 0; i50 < b_loop_ub; i50++) {
                      ft += j->data[i50] * b->data[i50];
                    }
                  } else {
                    ft = 0.0;
                    b_loop_ub = j->size[1];
                    for (i50 = 0; i50 < b_loop_ub; i50++) {
                      ft += j->data[i50] * b->data[i50];
                    }
                  }

                  err = (ft / b_sum(j) - des->data[x_tmp]) * wt->data[x_tmp];
                  ft = (double)nut * err;
                  dtemp = ft - comp;
                  if (dtemp > 0.0) {
                    comp = ft;
                    b_l--;
                  } else {
                    exitg2 = true;
                  }
                }

                dden = iext->data[(int)b_j - 1];
                iext->data[(int)b_j - 1] = b_l + 1.0;
                b_j++;
                jchnge++;
              } else {
                dden = iext->data[(int)b_j - 1];
                b_j++;
              }
            }
          }

          do {
            exitg3 = 0;
            if (b_j == (nfcns + 1.0) + 1.0) {
              ft = iext->data[0];
              if ((k1 > ft) || (rtIsNaN(k1) && (!rtIsNaN(ft)))) {
                k1 = ft;
              }

              ft = iext->data[kkk];
              if ((knz < ft) || (rtIsNaN(knz) && (!rtIsNaN(ft)))) {
                knz = ft;
              }

              nut1 = nut;
              nut = -nu;
              comp *= 1.00001;
              luck = 1;
              flag = 1;
              b_l = 1.0;
              exitg2 = false;
              while ((!exitg2) && (b_l < k1)) {
                /*  gee */
                x_tmp = (int)b_l - 1;
                ft = cos(6.2831853071795862 * grid->data[x_tmp]);
                i50 = b_wt->size[0] * b_wt->size[1];
                b_wt->size[0] = 1;
                b_wt->size[1] = x->size[1];
                emxEnsureCapacity_real_T(b_wt, i50);
                b_loop_ub = x->size[0] * x->size[1];
                for (i50 = 0; i50 < b_loop_ub; i50++) {
                  b_wt->data[i50] = ft - x->data[i50];
                }

                b_rdivide_helper(ad, b_wt, j);
                i50 = b->size[0];
                b->size[0] = y->size[1];
                emxEnsureCapacity_real_T(b, i50);
                b_loop_ub = y->size[1];
                for (i50 = 0; i50 < b_loop_ub; i50++) {
                  b->data[i50] = y->data[i50];
                }

                if ((j->size[1] == 1) || (b->size[0] == 1)) {
                  ft = 0.0;
                  b_loop_ub = j->size[1];
                  for (i50 = 0; i50 < b_loop_ub; i50++) {
                    ft += j->data[i50] * b->data[i50];
                  }
                } else {
                  ft = 0.0;
                  b_loop_ub = j->size[1];
                  for (i50 = 0; i50 < b_loop_ub; i50++) {
                    ft += j->data[i50] * b->data[i50];
                  }
                }

                err = (ft / b_sum(j) - des->data[x_tmp]) * wt->data[x_tmp];
                dtemp = err * -(double)nu - comp;
                if (dtemp > 0.0) {
                  comp = -(double)nu * err;
                  b_l++;
                  exitg4 = false;
                  while ((!exitg4) && (b_l < k1)) {
                    /*  gee */
                    x_tmp = (int)b_l - 1;
                    ft = cos(6.2831853071795862 * grid->data[x_tmp]);
                    i50 = b_wt->size[0] * b_wt->size[1];
                    b_wt->size[0] = 1;
                    b_wt->size[1] = x->size[1];
                    emxEnsureCapacity_real_T(b_wt, i50);
                    b_loop_ub = x->size[0] * x->size[1];
                    for (i50 = 0; i50 < b_loop_ub; i50++) {
                      b_wt->data[i50] = ft - x->data[i50];
                    }

                    b_rdivide_helper(ad, b_wt, j);
                    i50 = b->size[0];
                    b->size[0] = y->size[1];
                    emxEnsureCapacity_real_T(b, i50);
                    b_loop_ub = y->size[1];
                    for (i50 = 0; i50 < b_loop_ub; i50++) {
                      b->data[i50] = y->data[i50];
                    }

                    if ((j->size[1] == 1) || (b->size[0] == 1)) {
                      b_j = 0.0;
                      b_loop_ub = j->size[1];
                      for (i50 = 0; i50 < b_loop_ub; i50++) {
                        b_j += j->data[i50] * b->data[i50];
                      }
                    } else {
                      b_j = 0.0;
                      b_loop_ub = j->size[1];
                      for (i50 = 0; i50 < b_loop_ub; i50++) {
                        b_j += j->data[i50] * b->data[i50];
                      }
                    }

                    err = (b_j / b_sum(j) - des->data[x_tmp]) * wt->data[x_tmp];
                    dtemp = -(double)nu * err - comp;
                    if (dtemp > 0.0) {
                      comp = -(double)nu * err;
                      b_l++;
                    } else {
                      exitg4 = true;
                    }
                  }

                  iext->data[i49] = b_l - 1.0;
                  b_j = ((nfcns + 1.0) + 1.0) + 1.0;
                  jchnge++;
                  flag = 0;
                  exitg2 = true;
                } else {
                  b_l++;
                }
              }

              if (flag != 0) {
                luck = 6;
                nut = -nut1;
                comp = b_y1 * 1.00001;
                b_l = ((double)ngrid + 1.0) - 1.0;
                exitg2 = false;
                while ((!exitg2) && (b_l > knz)) {
                  /*  gee */
                  x_tmp = (int)b_l - 1;
                  ft = cos(6.2831853071795862 * grid->data[x_tmp]);
                  i50 = b_wt->size[0] * b_wt->size[1];
                  b_wt->size[0] = 1;
                  b_wt->size[1] = x->size[1];
                  emxEnsureCapacity_real_T(b_wt, i50);
                  b_loop_ub = x->size[0] * x->size[1];
                  for (i50 = 0; i50 < b_loop_ub; i50++) {
                    b_wt->data[i50] = ft - x->data[i50];
                  }

                  b_rdivide_helper(ad, b_wt, j);
                  i50 = b->size[0];
                  b->size[0] = y->size[1];
                  emxEnsureCapacity_real_T(b, i50);
                  b_loop_ub = y->size[1];
                  for (i50 = 0; i50 < b_loop_ub; i50++) {
                    b->data[i50] = y->data[i50];
                  }

                  if ((j->size[1] == 1) || (b->size[0] == 1)) {
                    ft = 0.0;
                    b_loop_ub = j->size[1];
                    for (i50 = 0; i50 < b_loop_ub; i50++) {
                      ft += j->data[i50] * b->data[i50];
                    }
                  } else {
                    ft = 0.0;
                    b_loop_ub = j->size[1];
                    for (i50 = 0; i50 < b_loop_ub; i50++) {
                      ft += j->data[i50] * b->data[i50];
                    }
                  }

                  err = (ft / b_sum(j) - des->data[x_tmp]) * wt->data[x_tmp];
                  dtemp = err * -(double)nut1 - comp;
                  if (dtemp > 0.0) {
                    comp = -(double)nut1 * err;
                    luck = 16;
                    b_l--;
                    exitg4 = false;
                    while ((!exitg4) && (b_l > knz)) {
                      /*  gee */
                      x_tmp = (int)b_l - 1;
                      ft = cos(6.2831853071795862 * grid->data[x_tmp]);
                      i50 = b_wt->size[0] * b_wt->size[1];
                      b_wt->size[0] = 1;
                      b_wt->size[1] = x->size[1];
                      emxEnsureCapacity_real_T(b_wt, i50);
                      b_loop_ub = x->size[0] * x->size[1];
                      for (i50 = 0; i50 < b_loop_ub; i50++) {
                        b_wt->data[i50] = ft - x->data[i50];
                      }

                      b_rdivide_helper(ad, b_wt, j);
                      i50 = b->size[0];
                      b->size[0] = y->size[1];
                      emxEnsureCapacity_real_T(b, i50);
                      b_loop_ub = y->size[1];
                      for (i50 = 0; i50 < b_loop_ub; i50++) {
                        b->data[i50] = y->data[i50];
                      }

                      if ((j->size[1] == 1) || (b->size[0] == 1)) {
                        b_j = 0.0;
                        b_loop_ub = j->size[1];
                        for (i50 = 0; i50 < b_loop_ub; i50++) {
                          b_j += j->data[i50] * b->data[i50];
                        }
                      } else {
                        b_j = 0.0;
                        b_loop_ub = j->size[1];
                        for (i50 = 0; i50 < b_loop_ub; i50++) {
                          b_j += j->data[i50] * b->data[i50];
                        }
                      }

                      err = (b_j / b_sum(j) - des->data[x_tmp]) * wt->data[x_tmp];
                      dtemp = -(double)nut1 * err - comp;
                      if (dtemp > 0.0) {
                        comp = -(double)nut1 * err;
                        b_l--;
                      } else {
                        exitg4 = true;
                      }
                    }

                    iext->data[i49] = b_l + 1.0;
                    b_j = ((nfcns + 1.0) + 1.0) + 1.0;
                    jchnge++;
                    flag = 0;
                    exitg2 = true;
                  } else {
                    b_l--;
                  }
                }

                if (flag != 0) {
                  flag34 = 0;
                  if (luck != 6) {
                    dden = (nfcns + 1.0) - nfcns;
                    if (2.0 > dden) {
                      i50 = -2;
                      i51 = 0;
                    } else {
                      i50 = -1;
                      i51 = (int)dden;
                    }

                    if (dden > (nfcns + 1.0) - 1.0) {
                      i52 = 1;
                      l_tmp = 0;
                    } else {
                      i52 = (int)dden;
                      l_tmp = (int)((nfcns + 1.0) - 1.0);
                    }

                    /*  Update index */
                    if (2.0 > dden) {
                      i53 = -2;
                      nu = 0;
                    } else {
                      i53 = -1;
                      nu = (int)dden;
                    }

                    if (dden > (nfcns + 1.0) - 1.0) {
                      nut = 1;
                      flag = 0;
                    } else {
                      nut = (int)dden;
                      flag = (int)((nfcns + 1.0) - 1.0);
                    }

                    x_tmp = sel->size[0] * sel->size[1];
                    sel->size[0] = 1;
                    b_loop_ub = i51 - i50;
                    sel->size[1] = (b_loop_ub + l_tmp) - i52;
                    emxEnsureCapacity_real_T(sel, x_tmp);
                    sel->data[0] = k1;
                    for (x_tmp = 0; x_tmp <= b_loop_ub - 3; x_tmp++) {
                      sel->data[x_tmp + 1] = iext->data[(i50 + x_tmp) + 2];
                    }

                    b_loop_ub = l_tmp - i52;
                    for (l_tmp = 0; l_tmp <= b_loop_ub; l_tmp++) {
                      sel->data[((l_tmp + i51) - i50) - 1] = iext->data[(i52 +
                        l_tmp) - 1];
                    }

                    b_loop_ub = sel->size[1];
                    i50 = r11->size[0] * r11->size[1];
                    r11->size[0] = 1;
                    r11->size[1] = b_loop_ub;
                    emxEnsureCapacity_int32_T(r11, i50);
                    for (i50 = 0; i50 < b_loop_ub; i50++) {
                      r11->data[i50] = i50;
                    }

                    i50 = b->size[0];
                    b_loop_ub = nu - i53;
                    b->size[0] = (b_loop_ub + flag) - nut;
                    emxEnsureCapacity_real_T(b, i50);
                    b->data[0] = k1;
                    for (i50 = 0; i50 <= b_loop_ub - 3; i50++) {
                      b->data[i50 + 1] = iext->data[(i53 + i50) + 2];
                    }

                    b_loop_ub = flag - nut;
                    for (i50 = 0; i50 <= b_loop_ub; i50++) {
                      b->data[((i50 + nu) - i53) - 1] = iext->data[(nut + i50) -
                        1];
                    }

                    b_loop_ub = r11->size[0] * r11->size[1];
                    for (i50 = 0; i50 < b_loop_ub; i50++) {
                      iext->data[r11->data[i50]] = b->data[i50];
                    }

                    jchnge++;
                  }

                  exitg3 = 1;
                }
              }
            } else {
              exitg3 = 1;
            }
          } while (exitg3 == 0);

          if ((flag34 != 0) && (b_j > (nfcns + 1.0) + 1.0)) {
            if (luck > 9) {
              if (2.0 > nfcns + 1.0) {
                i50 = 0;
                i51 = 0;
              } else {
                i50 = 1;
                i51 = loop_ub;
              }

              if (nfcns + 1.0 > (nfcns + 1.0) - 1.0) {
                i52 = 1;
                l_tmp = 0;
              } else {
                i52 = loop_ub;
                l_tmp = (int)((nfcns + 1.0) - 1.0);
              }

              /*  Update index */
              if (2.0 > nfcns + 1.0) {
                i53 = 0;
                nu = 0;
              } else {
                i53 = 1;
                nu = loop_ub;
              }

              if (nfcns + 1.0 > (nfcns + 1.0) - 1.0) {
                nut = 1;
                flag = 0;
              } else {
                nut = loop_ub;
                flag = (int)((nfcns + 1.0) - 1.0);
              }

              x_tmp = sel->size[0] * sel->size[1];
              sel->size[0] = 1;
              b_loop_ub = i51 - i50;
              kkk = (b_loop_ub + l_tmp) - i52;
              sel->size[1] = kkk + 3;
              emxEnsureCapacity_real_T(sel, x_tmp);
              for (x_tmp = 0; x_tmp < b_loop_ub; x_tmp++) {
                sel->data[x_tmp] = iext->data[i50 + x_tmp];
              }

              b_loop_ub = l_tmp - i52;
              for (l_tmp = 0; l_tmp <= b_loop_ub; l_tmp++) {
                sel->data[(l_tmp + i51) - i50] = iext->data[(i52 + l_tmp) - 1];
              }

              sel->data[kkk + 1] = iext->data[i49];
              sel->data[kkk + 2] = iext->data[i49];
              b_loop_ub = sel->size[1];
              i50 = r11->size[0] * r11->size[1];
              r11->size[0] = 1;
              r11->size[1] = b_loop_ub;
              emxEnsureCapacity_int32_T(r11, i50);
              for (i50 = 0; i50 < b_loop_ub; i50++) {
                r11->data[i50] = i50;
              }

              dden = iext->data[i49];
              ft = iext->data[i49];
              i49 = b->size[0];
              b_loop_ub = nu - i53;
              i50 = (b_loop_ub + flag) - nut;
              b->size[0] = i50 + 3;
              emxEnsureCapacity_real_T(b, i49);
              for (i49 = 0; i49 < b_loop_ub; i49++) {
                b->data[i49] = iext->data[i53 + i49];
              }

              b_loop_ub = flag - nut;
              for (i49 = 0; i49 <= b_loop_ub; i49++) {
                b->data[(i49 + nu) - i53] = iext->data[(nut + i49) - 1];
              }

              b->data[i50 + 1] = dden;
              b->data[i50 + 2] = ft;
              b_loop_ub = r11->size[0] * r11->size[1];
              for (i49 = 0; i49 < b_loop_ub; i49++) {
                iext->data[r11->data[i49]] = b->data[i49];
              }

              jchnge++;
            } else {
              if ((b_y1 < comp) || (rtIsNaN(b_y1) && (!rtIsNaN(comp)))) {
                b_y1 = comp;
              }

              k1 = iext->data[i49];
              comp = b_y1 * 1.00001;
              b_l = ((double)ngrid + 1.0) - 1.0;
              exitg2 = false;
              while ((!exitg2) && (b_l > knz)) {
                /*  gee */
                x_tmp = (int)b_l - 1;
                ft = cos(6.2831853071795862 * grid->data[x_tmp]);
                i50 = b_wt->size[0] * b_wt->size[1];
                b_wt->size[0] = 1;
                b_wt->size[1] = x->size[1];
                emxEnsureCapacity_real_T(b_wt, i50);
                b_loop_ub = x->size[0] * x->size[1];
                for (i50 = 0; i50 < b_loop_ub; i50++) {
                  b_wt->data[i50] = ft - x->data[i50];
                }

                b_rdivide_helper(ad, b_wt, j);
                i50 = b->size[0];
                b->size[0] = y->size[1];
                emxEnsureCapacity_real_T(b, i50);
                b_loop_ub = y->size[1];
                for (i50 = 0; i50 < b_loop_ub; i50++) {
                  b->data[i50] = y->data[i50];
                }

                if ((j->size[1] == 1) || (b->size[0] == 1)) {
                  b_j = 0.0;
                  b_loop_ub = j->size[1];
                  for (i50 = 0; i50 < b_loop_ub; i50++) {
                    b_j += j->data[i50] * b->data[i50];
                  }
                } else {
                  b_j = 0.0;
                  b_loop_ub = j->size[1];
                  for (i50 = 0; i50 < b_loop_ub; i50++) {
                    b_j += j->data[i50] * b->data[i50];
                  }
                }

                err = (b_j / b_sum(j) - des->data[x_tmp]) * wt->data[x_tmp];
                dtemp = err * -(double)nut1 - comp;
                if (dtemp > 0.0) {
                  comp = -(double)nut1 * err;
                  luck += 10;
                  b_l--;
                  exitg4 = false;
                  while ((!exitg4) && (b_l > knz)) {
                    /*  gee */
                    x_tmp = (int)b_l - 1;
                    ft = cos(6.2831853071795862 * grid->data[x_tmp]);
                    i50 = b_wt->size[0] * b_wt->size[1];
                    b_wt->size[0] = 1;
                    b_wt->size[1] = x->size[1];
                    emxEnsureCapacity_real_T(b_wt, i50);
                    b_loop_ub = x->size[0] * x->size[1];
                    for (i50 = 0; i50 < b_loop_ub; i50++) {
                      b_wt->data[i50] = ft - x->data[i50];
                    }

                    b_rdivide_helper(ad, b_wt, j);
                    i50 = b->size[0];
                    b->size[0] = y->size[1];
                    emxEnsureCapacity_real_T(b, i50);
                    b_loop_ub = y->size[1];
                    for (i50 = 0; i50 < b_loop_ub; i50++) {
                      b->data[i50] = y->data[i50];
                    }

                    if ((j->size[1] == 1) || (b->size[0] == 1)) {
                      b_j = 0.0;
                      b_loop_ub = j->size[1];
                      for (i50 = 0; i50 < b_loop_ub; i50++) {
                        b_j += j->data[i50] * b->data[i50];
                      }
                    } else {
                      b_j = 0.0;
                      b_loop_ub = j->size[1];
                      for (i50 = 0; i50 < b_loop_ub; i50++) {
                        b_j += j->data[i50] * b->data[i50];
                      }
                    }

                    err = (b_j / b_sum(j) - des->data[x_tmp]) * wt->data[x_tmp];
                    ft = -(double)nut1 * err;
                    dtemp = ft - comp;
                    if (dtemp > 0.0) {
                      comp = ft;
                      b_l--;
                    } else {
                      exitg4 = true;
                    }
                  }

                  iext->data[i49] = b_l + 1.0;
                  jchnge++;
                  if (2.0 > nfcns + 1.0) {
                    i50 = 0;
                    i51 = 0;
                  } else {
                    i50 = 1;
                    i51 = loop_ub;
                  }

                  if (nfcns + 1.0 > (nfcns + 1.0) - 1.0) {
                    i52 = 1;
                    l_tmp = 0;
                  } else {
                    i52 = loop_ub;
                    l_tmp = (int)((nfcns + 1.0) - 1.0);
                  }

                  /*  Update index */
                  if (2.0 > nfcns + 1.0) {
                    i53 = 0;
                    nu = 0;
                  } else {
                    i53 = 1;
                    nu = loop_ub;
                  }

                  if (nfcns + 1.0 > (nfcns + 1.0) - 1.0) {
                    nut = 1;
                    flag = 0;
                  } else {
                    nut = loop_ub;
                    flag = (int)((nfcns + 1.0) - 1.0);
                  }

                  x_tmp = sel->size[0] * sel->size[1];
                  sel->size[0] = 1;
                  b_loop_ub = i51 - i50;
                  kkk = (b_loop_ub + l_tmp) - i52;
                  sel->size[1] = kkk + 3;
                  emxEnsureCapacity_real_T(sel, x_tmp);
                  for (x_tmp = 0; x_tmp < b_loop_ub; x_tmp++) {
                    sel->data[x_tmp] = iext->data[i50 + x_tmp];
                  }

                  b_loop_ub = l_tmp - i52;
                  for (l_tmp = 0; l_tmp <= b_loop_ub; l_tmp++) {
                    sel->data[(l_tmp + i51) - i50] = iext->data[(i52 + l_tmp) -
                      1];
                  }

                  sel->data[kkk + 1] = iext->data[i49];
                  sel->data[kkk + 2] = iext->data[i49];
                  b_loop_ub = sel->size[1];
                  i50 = r11->size[0] * r11->size[1];
                  r11->size[0] = 1;
                  r11->size[1] = b_loop_ub;
                  emxEnsureCapacity_int32_T(r11, i50);
                  for (i50 = 0; i50 < b_loop_ub; i50++) {
                    r11->data[i50] = i50;
                  }

                  dden = iext->data[i49];
                  ft = iext->data[i49];
                  i49 = b->size[0];
                  b_loop_ub = nu - i53;
                  i50 = (b_loop_ub + flag) - nut;
                  b->size[0] = i50 + 3;
                  emxEnsureCapacity_real_T(b, i49);
                  for (i49 = 0; i49 < b_loop_ub; i49++) {
                    b->data[i49] = iext->data[i53 + i49];
                  }

                  b_loop_ub = flag - nut;
                  for (i49 = 0; i49 <= b_loop_ub; i49++) {
                    b->data[(i49 + nu) - i53] = iext->data[(nut + i49) - 1];
                  }

                  b->data[i50 + 1] = dden;
                  b->data[i50 + 2] = ft;
                  b_loop_ub = r11->size[0] * r11->size[1];
                  for (i49 = 0; i49 < b_loop_ub; i49++) {
                    iext->data[r11->data[i49]] = b->data[i49];
                  }

                  exitg2 = true;
                } else {
                  b_l--;
                }
              }

              if (luck != 6) {
                dden = (nfcns + 1.0) - nfcns;
                if (2.0 > dden) {
                  i49 = -2;
                  i50 = 0;
                } else {
                  i49 = -1;
                  i50 = (int)dden;
                }

                dden = (nfcns + 1.0) - nfcns;
                if (dden > (nfcns + 1.0) - 1.0) {
                  i51 = 1;
                  i52 = 0;
                } else {
                  i51 = (int)dden;
                  i52 = (int)((nfcns + 1.0) - 1.0);
                }

                /*  Update index */
                dden = (nfcns + 1.0) - nfcns;
                if (2.0 > dden) {
                  l_tmp = -2;
                  i53 = 0;
                } else {
                  l_tmp = -1;
                  i53 = (int)dden;
                }

                dden = (nfcns + 1.0) - nfcns;
                if (dden > (nfcns + 1.0) - 1.0) {
                  nu = 1;
                  nut = 0;
                } else {
                  nu = (int)dden;
                  nut = (int)((nfcns + 1.0) - 1.0);
                }

                flag = sel->size[0] * sel->size[1];
                sel->size[0] = 1;
                b_loop_ub = i50 - i49;
                sel->size[1] = (b_loop_ub + i52) - i51;
                emxEnsureCapacity_real_T(sel, flag);
                sel->data[0] = k1;
                for (flag = 0; flag <= b_loop_ub - 3; flag++) {
                  sel->data[flag + 1] = iext->data[(i49 + flag) + 2];
                }

                b_loop_ub = i52 - i51;
                for (i52 = 0; i52 <= b_loop_ub; i52++) {
                  sel->data[((i52 + i50) - i49) - 1] = iext->data[(i51 + i52) -
                    1];
                }

                b_loop_ub = sel->size[1];
                i49 = r11->size[0] * r11->size[1];
                r11->size[0] = 1;
                r11->size[1] = b_loop_ub;
                emxEnsureCapacity_int32_T(r11, i49);
                for (i49 = 0; i49 < b_loop_ub; i49++) {
                  r11->data[i49] = i49;
                }

                i49 = b->size[0];
                b_loop_ub = i53 - l_tmp;
                b->size[0] = (b_loop_ub + nut) - nu;
                emxEnsureCapacity_real_T(b, i49);
                b->data[0] = k1;
                for (i49 = 0; i49 <= b_loop_ub - 3; i49++) {
                  b->data[i49 + 1] = iext->data[(l_tmp + i49) + 2];
                }

                b_loop_ub = nut - nu;
                for (i49 = 0; i49 <= b_loop_ub; i49++) {
                  b->data[((i49 + i53) - l_tmp) - 1] = iext->data[(nu + i49) - 1];
                }

                b_loop_ub = r11->size[0] * r11->size[1];
                for (i49 = 0; i49 < b_loop_ub; i49++) {
                  iext->data[r11->data[i49]] = b->data[i49];
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
    jchnge = -1.0;
    err = -1.0;

    /*  initialize memory */
    /* x(nzz) = -2; */
    i49 = l->size[0] * l->size[1];
    l->size[0] = 1;
    l->size[1] = x->size[1] + 1;
    emxEnsureCapacity_real_T(l, i49);
    b_loop_ub = x->size[1];
    for (i49 = 0; i49 < b_loop_ub; i49++) {
      l->data[i49] = x->data[i49];
    }

    l->data[x->size[1]] = -2.0;
    k1 = 2.0 * nfcns - 1.0;
    knz = 1.0 / k1;
    b_l = 1.0;
    kkk = 0;
    if (((edge[0] == 0.0) && (edge[3] == 0.5)) || (nfcns <= 3.0)) {
      kkk = 1;
    }

    if (kkk != 1) {
      dtemp = cos(6.2831853071795862 * grid->data[0]);
      dnum = cos(6.2831853071795862 * grid->data[varargin_2 - 1]);
      ft = dtemp - dnum;
      err = 2.0 / ft;
      jchnge = -(dtemp + dnum) / ft;
    }

    b_loop_ub = (int)nfcns;
    i49 = a->size[0] * a->size[1];
    a->size[0] = 1;
    a->size[1] = b_loop_ub;
    emxEnsureCapacity_real_T(a, i49);
    for (nu = 0; nu < b_loop_ub; nu++) {
      ft = ((1.0 + (double)nu) - 1.0) * knz;
      dnum = cos(6.2831853071795862 * ft);
      if (kkk != 1) {
        dnum = (dnum - jchnge) / err;
        dden = dnum;
        b_acos(&dden);
        ft = dden / 6.2831853071795862;
      }

      dden = l->data[(int)b_l - 1];
      while ((dnum <= dden) && (dden - dnum >= 1.0E-6)) {
        b_l++;
        dden = l->data[(int)b_l - 1];
      }

      if (fabs(dnum - dden) < 1.0E-6) {
        a->data[nu] = y->data[(int)b_l - 1];
      } else {
        /*  gee */
        b_x = cos(6.2831853071795862 * ft);
        i49 = b_wt->size[0] * b_wt->size[1];
        b_wt->size[0] = 1;
        b_wt->size[1] = loop_ub;
        emxEnsureCapacity_real_T(b_wt, i49);
        for (i49 = 0; i49 < loop_ub; i49++) {
          b_wt->data[i49] = b_x - l->data[i49];
        }

        b_rdivide_helper(ad, b_wt, j);
        i49 = b->size[0];
        b->size[0] = y->size[1];
        emxEnsureCapacity_real_T(b, i49);
        nut = y->size[1];
        for (i49 = 0; i49 < nut; i49++) {
          b->data[i49] = y->data[i49];
        }

        if ((j->size[1] == 1) || (b->size[0] == 1)) {
          b_j = 0.0;
          nut = j->size[1];
          for (i49 = 0; i49 < nut; i49++) {
            b_j += j->data[i49] * b->data[i49];
          }
        } else {
          b_j = 0.0;
          nut = j->size[1];
          for (i49 = 0; i49 < nut; i49++) {
            b_j += j->data[i49] * b->data[i49];
          }
        }

        a->data[nu] = b_j / b_sum(j);
      }

      if (1 < (int)(b_l - 1.0)) {
        b_l = (int)(b_l - 1.0);
      } else {
        b_l = 1.0;
      }
    }

    dden = 6.2831853071795862 / k1;
    i49 = j->size[0] * j->size[1];
    j->size[0] = 1;
    j->size[1] = b_loop_ub;
    emxEnsureCapacity_real_T(j, i49);
    for (nu = 0; nu < b_loop_ub; nu++) {
      dnum = ((1.0 + (double)nu) - 1.0) * dden;
      if (nfcns - 1.0 < 1.0) {
        j->data[nu] = a->data[0];
      } else {
        if (2.0 > nfcns) {
          i49 = 0;
          i50 = 0;
          y->size[0] = 1;
          y->size[1] = 0;
        } else {
          i49 = 1;
          i50 = b_loop_ub;
          i51 = y->size[0] * y->size[1];
          y->size[0] = 1;
          loop_ub = (int)((nfcns - 1.0) - 1.0);
          y->size[1] = loop_ub + 1;
          emxEnsureCapacity_real_T(y, i51);
          for (i51 = 0; i51 <= loop_ub; i51++) {
            y->data[i51] = 1.0 + (double)i51;
          }
        }

        i51 = y->size[0] * y->size[1];
        i52 = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity_real_T(y, i52);
        loop_ub = i51 - 1;
        for (i51 = 0; i51 <= loop_ub; i51++) {
          y->data[i51] *= dnum;
        }

        b_cos(y);
        i51 = y->size[0] * y->size[1];
        i52 = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity_real_T(y, i52);
        loop_ub = i51 - 1;
        for (i51 = 0; i51 <= loop_ub; i51++) {
          y->data[i51] *= 2.0;
        }

        i51 = b->size[0];
        loop_ub = i50 - i49;
        b->size[0] = loop_ub;
        emxEnsureCapacity_real_T(b, i51);
        for (i50 = 0; i50 < loop_ub; i50++) {
          b->data[i50] = a->data[i49 + i50];
        }

        if ((y->size[1] == 1) || (loop_ub == 1)) {
          ft = 0.0;
          loop_ub = y->size[1];
          for (i49 = 0; i49 < loop_ub; i49++) {
            ft += y->data[i49] * b->data[i49];
          }
        } else {
          ft = 0.0;
          loop_ub = y->size[1];
          for (i49 = 0; i49 < loop_ub; i49++) {
            ft += y->data[i49] * b->data[i49];
          }
        }

        j->data[nu] = a->data[0] + ft;
      }
    }

    if (2.0 > nfcns) {
      i49 = -1;
      i50 = 0;
    } else {
      i49 = 0;
      i50 = b_loop_ub;
    }

    i51 = iext->size[0];
    loop_ub = i50 - i49;
    iext->size[0] = loop_ub;
    emxEnsureCapacity_real_T(iext, i51);
    iext->data[0] = j->data[0] / k1;
    for (i50 = 0; i50 <= loop_ub - 2; i50++) {
      iext->data[i50 + 1] = 2.0 * j->data[(i49 + i50) + 1] / k1;
    }

    i49 = j->size[0] * j->size[1];
    j->size[0] = 1;
    j->size[1] = b_loop_ub;
    emxEnsureCapacity_real_T(j, i49);
    for (i49 = 0; i49 < b_loop_ub; i49++) {
      j->data[i49] = 0.0;
    }

    i49 = l->size[0] * l->size[1];
    l->size[0] = 1;
    loop_ub = b_loop_ub - 1;
    l->size[1] = loop_ub;
    emxEnsureCapacity_real_T(l, i49);
    for (i49 = 0; i49 < loop_ub; i49++) {
      l->data[i49] = 0.0;
    }

    if (kkk != 1) {
      j->data[0] = 2.0 * iext->data[loop_ub] * jchnge + iext->data[b_loop_ub - 2];
      j->data[1] = 2.0 * err * iext->data[loop_ub];
      l->data[0] = iext->data[b_loop_ub - 3] - iext->data[loop_ub];
      for (nu = 0; nu <= b_loop_ub - 3; nu++) {
        if (2 + nu == loop_ub) {
          err /= 2.0;
          jchnge /= 2.0;
        }

        j->data[nu + 2] = 0.0;
        i49 = x->size[0] * x->size[1];
        x->size[0] = 1;
        nut = (int)((2.0 + (double)nu) - 1.0);
        x->size[1] = nut + 1;
        emxEnsureCapacity_real_T(x, i49);
        for (i49 = 0; i49 <= nut; i49++) {
          x->data[i49] = 1.0 + (double)i49;
        }

        nut = x->size[0] * x->size[1] - 1;
        for (i49 = 0; i49 <= nut; i49++) {
          a->data[(int)x->data[i49] - 1] = j->data[(int)x->data[i49] - 1];
        }

        ft = 2.0 * jchnge;
        nut = x->size[0] * x->size[1] - 1;
        for (i49 = 0; i49 <= nut; i49++) {
          j->data[(int)x->data[i49] - 1] = ft * a->data[(int)x->data[i49] - 1];
        }

        j->data[1] += 2.0 * a->data[0] * err;
        i49 = sel->size[0] * sel->size[1];
        sel->size[0] = 1;
        sel->size[1] = nu + 1;
        emxEnsureCapacity_real_T(sel, i49);
        for (i49 = 0; i49 <= nu; i49++) {
          sel->data[i49] = 1.0 + (double)i49;
        }

        i49 = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = sel->size[1];
        emxEnsureCapacity_real_T(y, i49);
        nut = sel->size[0] * sel->size[1];
        for (i49 = 0; i49 < nut; i49++) {
          y->data[i49] = a->data[(int)sel->data[i49]];
        }

        nut = sel->size[0] * sel->size[1];
        i49 = b->size[0];
        b->size[0] = nut;
        emxEnsureCapacity_real_T(b, i49);
        for (i49 = 0; i49 < nut; i49++) {
          b->data[i49] = (j->data[(int)sel->data[i49] - 1] + l->data[(int)
                          sel->data[i49] - 1]) + err * y->data[i49];
        }

        nut = b->size[0];
        for (i49 = 0; i49 < nut; i49++) {
          j->data[(int)sel->data[i49] - 1] = b->data[i49];
        }

        i49 = sel->size[0] * sel->size[1];
        sel->size[0] = 1;
        sel->size[1] = nu + 1;
        emxEnsureCapacity_real_T(sel, i49);
        for (i49 = 0; i49 <= nu; i49++) {
          sel->data[i49] = 3.0 + (double)i49;
        }

        i49 = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = sel->size[1];
        emxEnsureCapacity_real_T(y, i49);
        nut = sel->size[0] * sel->size[1];
        for (i49 = 0; i49 < nut; i49++) {
          y->data[i49] = a->data[(int)sel->data[i49] - 2];
        }

        nut = sel->size[0] * sel->size[1];
        i49 = b->size[0];
        b->size[0] = nut;
        emxEnsureCapacity_real_T(b, i49);
        for (i49 = 0; i49 < nut; i49++) {
          b->data[i49] = j->data[(int)sel->data[i49] - 1] + err * y->data[i49];
        }

        nut = b->size[0];
        for (i49 = 0; i49 < nut; i49++) {
          j->data[(int)sel->data[i49] - 1] = b->data[i49];
        }

        if (2 + nu != loop_ub) {
          nut = x->size[0] * x->size[1] - 1;
          for (i49 = 0; i49 <= nut; i49++) {
            l->data[(int)x->data[i49] - 1] = -a->data[(int)x->data[i49] - 1];
          }

          l->data[0] += iext->data[(b_loop_ub - nu) - 4];
        }
      }

      for (i49 = 0; i49 < b_loop_ub; i49++) {
        iext->data[i49] = j->data[i49];
      }
    }

    /*  alpha must be at lease >=3 */
    if (nfcns <= 3.0) {
      /* alpha(nfcns + 1) = 0; */
      /* alpha(nfcns + 2) = 0; */
      i49 = j->size[0] * j->size[1];
      j->size[0] = 1;
      j->size[1] = iext->size[0] + 2;
      emxEnsureCapacity_real_T(j, i49);
      b_loop_ub = iext->size[0];
      for (i49 = 0; i49 < b_loop_ub; i49++) {
        j->data[i49] = iext->data[i49];
      }

      j->data[iext->size[0]] = 0.0;
      j->data[iext->size[0] + 1] = 0.0;
    } else {
      i49 = j->size[0] * j->size[1];
      j->size[0] = 1;
      j->size[1] = iext->size[0];
      emxEnsureCapacity_real_T(j, i49);
      b_loop_ub = iext->size[0];
      for (i49 = 0; i49 < b_loop_ub; i49++) {
        j->data[i49] = iext->data[i49];
      }
    }

    /* alpha=alpha'; */
    /*  now that's done! */
    if (nodd != 0.0) {
      i49 = y->size[0] * y->size[1];
      y->size[0] = 1;
      loop_ub = (int)-((0.0 - (nfcns - 1.0)) - -1.0);
      y->size[1] = loop_ub + 1;
      emxEnsureCapacity_real_T(y, i49);
      for (i49 = 0; i49 <= loop_ub; i49++) {
        y->data[i49] = j->data[(int)((nfcns + 1.0) + (-1.0 - (double)i49)) - 1];
      }

      i49 = h->size[0] * h->size[1];
      h->size[0] = 1;
      h->size[1] = y->size[1] + 1;
      emxEnsureCapacity_real_T(h, i49);
      loop_ub = y->size[1];
      for (i49 = 0; i49 < loop_ub; i49++) {
        h->data[h->size[0] * i49] = 0.5 * y->data[i49];
      }

      h->data[h->size[0] * y->size[1]] = j->data[0];
    } else {
      if ((nfcns - (nfcns - 1.0)) + 2.0 > nfcns) {
        i49 = 0;
        i50 = 1;
      } else {
        i49 = loop_ub;
        i50 = -1;
      }

      i51 = sel->size[0] * sel->size[1];
      sel->size[0] = 1;
      b_loop_ub = (int)-((0.0 - (nfcns - 1.0)) - -2.0);
      sel->size[1] = b_loop_ub + 1;
      emxEnsureCapacity_real_T(sel, i51);
      for (i51 = 0; i51 <= b_loop_ub; i51++) {
        sel->data[i51] = j->data[(int)((nfcns + 1.0) + (double)(-2 - i51)) - 1];
      }

      i51 = h->size[0] * h->size[1];
      h->size[0] = 1;
      h->size[1] = 2 + sel->size[1];
      emxEnsureCapacity_real_T(h, i51);
      h->data[0] = 0.25 * j->data[loop_ub];
      loop_ub = sel->size[1];
      for (i51 = 0; i51 < loop_ub; i51++) {
        h->data[h->size[0] * (i51 + 1)] = 0.25 * (sel->data[i51] + j->data[i49 +
          i50 * i51]);
      }

      h->data[h->size[0] * (1 + sel->size[1])] = 0.25 * (2.0 * j->data[0] +
        j->data[1]);
    }
  }

  emxFree_int32_T(&b_iext);
  emxFree_real_T(&b_wt);
  emxFree_real_T(&a);
  emxFree_real_T(&b);
  emxFree_int32_T(&r11);
  emxFree_real_T(&sel);
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
  int ii;
  int idx;
  int ii_size_idx_1;
  bool exitg1;
  signed char ii_data[1];
  int b_ii_size_idx_1;
  bool b_a;
  int b_ii_data[1];
  signed char varargin_1_data[1];
  int varargin_2_data[1];
  int c_ii_size_idx_1;
  double d1;

  /*  end freqs */
  ii = b_size[1];
  idx = 0;
  ii_size_idx_1 = 1;
  exitg1 = false;
  while ((!exitg1) && (ii > 0)) {
    if (b_data[ii - 1] != 0.0) {
      idx = 1;
      ii_data[0] = (signed char)ii;
      exitg1 = true;
    } else {
      ii--;
    }
  }

  if (idx == 0) {
    ii_size_idx_1 = 0;
  }

  ii = a->size[1];
  idx = 0;
  b_ii_size_idx_1 = 1;
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
    b_ii_size_idx_1 = 0;
  }

  for (idx = 0; idx < ii_size_idx_1; idx++) {
    varargin_1_data[0] = (signed char)((signed char)b_size[1] - ii_data[0]);
  }

  ii = a->size[1];
  for (idx = 0; idx < b_ii_size_idx_1; idx++) {
    varargin_2_data[0] = ii - b_ii_data[0];
  }

  if (ii_size_idx_1 <= b_ii_size_idx_1) {
    c_ii_size_idx_1 = ii_size_idx_1;
  } else {
    c_ii_size_idx_1 = 0;
  }

  if (0 <= c_ii_size_idx_1 - 1) {
    b_ii_data[0] = (int)fmin(varargin_1_data[0], varargin_2_data[0]);
  }

  if (b_ii_data[0] > 0) {
    d1 = (double)b_size[1] - (double)b_ii_data[0];
    if (1.0 > d1) {
      ii = 0;
    } else {
      ii = (int)d1;
    }

    bR_size[0] = 1;
    bR_size[1] = ii;
    if (0 <= ii - 1) {
      memcpy(&bR_data[0], &b_data[0], (unsigned int)(ii * (int)sizeof(double)));
    }

    d1 = (double)a->size[1] - (double)b_ii_data[0];
    if (1.0 > d1) {
      ii = 0;
    } else {
      ii = (int)d1;
    }

    idx = aR->size[0] * aR->size[1];
    aR->size[0] = 1;
    aR->size[1] = ii;
    emxEnsureCapacity_creal_T(aR, idx);
    for (idx = 0; idx < ii; idx++) {
      aR->data[idx] = a->data[idx];
    }
  } else {
    bR_size[0] = 1;
    bR_size[1] = b_size[1];
    ii = b_size[0] * b_size[1];
    if (0 <= ii - 1) {
      memcpy(&bR_data[0], &b_data[0], (unsigned int)(ii * (int)sizeof(double)));
    }

    idx = aR->size[0] * aR->size[1];
    aR->size[0] = 1;
    aR->size[1] = a->size[1];
    emxEnsureCapacity_creal_T(aR, idx);
    ii = a->size[0] * a->size[1];
    for (idx = 0; idx < ii; idx++) {
      aR->data[idx] = a->data[idx];
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
  double d3;
  double d4;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    d3 = fabs(u0);
    d4 = fabs(u1);
    if (rtIsInf(u1)) {
      if (d3 == 1.0) {
        y = 1.0;
      } else if (d3 > 1.0) {
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
    } else if (d4 == 0.0) {
      y = 1.0;
    } else if (d4 == 1.0) {
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
  double q;
  if (rtIsNaN(u0) || rtIsInf(u0) || (rtIsNaN(u1) || rtIsInf(u1))) {
    y = rtNaN;
  } else if ((u1 != 0.0) && (u1 != trunc(u1))) {
    q = fabs(u0 / u1);
    if (fabs(q - floor(q + 0.5)) <= DBL_EPSILON * q) {
      y = 0.0 * u0;
    } else {
      y = fmod(u0, u1);
    }
  } else {
    y = fmod(u0, u1);
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
 * FREQZ_CG Frequency response of digital filter with codegen support
 *
 *  This function is based on 'freqz' by The MathWorks Inc.
 * Arguments    : const double b[19]
 *                const emxArray_real_T *w
 *                double Fs
 *                emxArray_creal_T *hh
 * Return Type  : void
 */
static void s_freqz_cg(const double b[19], const emxArray_real_T *w, double Fs,
  emxArray_creal_T *hh)
{
  emxArray_real_T *digw;
  int i40;
  int loop_ub;
  emxArray_creal_T *s;
  emxArray_creal_T *y;
  bool b8;
  int k;
  double b_re;
  double s_re;
  double s_im;
  emxInit_real_T(&digw, 2);

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
  i40 = digw->size[0] * digw->size[1];
  digw->size[0] = 1;
  digw->size[1] = w->size[1];
  emxEnsureCapacity_real_T(digw, i40);
  loop_ub = w->size[0] * w->size[1];
  for (i40 = 0; i40 < loop_ub; i40++) {
    digw->data[i40] = 6.2831853071795862 * w->data[i40] / Fs;
  }

  emxInit_creal_T(&s, 2);

  /*  Convert from Hz to rad/sample for computational purposes */
  i40 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = digw->size[1];
  emxEnsureCapacity_creal_T(s, i40);
  loop_ub = digw->size[0] * digw->size[1];
  for (i40 = 0; i40 < loop_ub; i40++) {
    s->data[i40].re = digw->data[i40] * 0.0;
    s->data[i40].im = digw->data[i40];
  }

  emxInit_creal_T(&y, 2);
  c_exp(s);

  /*  Digital frequency must be used for this calculation */
  i40 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = s->size[1];
  emxEnsureCapacity_creal_T(y, i40);
  b8 = (y->size[1] == 0);
  if (!b8) {
    i40 = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity_creal_T(y, i40);
    loop_ub = y->size[1];
    for (i40 = 0; i40 < loop_ub; i40++) {
      y->data[i40].re = b[0];
      y->data[i40].im = 0.0;
    }

    for (k = 0; k < 18; k++) {
      i40 = s->size[0] * s->size[1];
      loop_ub = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = s->size[1];
      emxEnsureCapacity_creal_T(y, loop_ub);
      b_re = b[k + 1];
      loop_ub = i40 - 1;
      for (i40 = 0; i40 <= loop_ub; i40++) {
        s_re = s->data[i40].re * y->data[i40].re - s->data[i40].im * y->data[i40]
          .im;
        s_im = s->data[i40].re * y->data[i40].im + s->data[i40].im * y->data[i40]
          .re;
        y->data[i40].re = s_re + b_re;
        y->data[i40].im = s_im;
      }
    }
  }

  i40 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = digw->size[1];
  emxEnsureCapacity_creal_T(s, i40);
  loop_ub = digw->size[0] * digw->size[1];
  for (i40 = 0; i40 < loop_ub; i40++) {
    b_re = digw->data[i40] * 0.0;
    s_re = digw->data[i40];
    s->data[i40].re = 18.0 * b_re;
    s->data[i40].im = 18.0 * s_re;
  }

  emxFree_real_T(&digw);
  c_exp(s);
  rdivide_helper(y, s, hh);

  /*  Generate the default structure to pass to freqzplot */
  /*  If rad/sample, Fs is empty */
  emxFree_creal_T(&y);
  emxFree_creal_T(&s);
}

/*
 * Arguments    : double x[2048]
 * Return Type  : void
 */
static void sinc(double x[2048])
{
  int k;
  for (k = 0; k < 2048; k++) {
    if (fabs(x[k]) < 1.0020841800044864E-292) {
      x[k] = 1.0;
    } else {
      x[k] *= 3.1415926535897931;
      x[k] = sin(x[k]) / x[k];
    }
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
 * FREQZ_CG Frequency response of digital filter with codegen support
 *
 *  This function is based on 'freqz' by The MathWorks Inc.
 * Arguments    : const double b[85]
 *                const emxArray_real_T *w
 *                double Fs
 *                emxArray_creal_T *hh
 * Return Type  : void
 */
static void t_freqz_cg(const double b[85], const emxArray_real_T *w, double Fs,
  emxArray_creal_T *hh)
{
  emxArray_real_T *digw;
  int i41;
  int loop_ub;
  emxArray_creal_T *s;
  emxArray_creal_T *y;
  bool b9;
  int k;
  double b_re;
  double s_re;
  double s_im;
  emxInit_real_T(&digw, 2);

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
  i41 = digw->size[0] * digw->size[1];
  digw->size[0] = 1;
  digw->size[1] = w->size[1];
  emxEnsureCapacity_real_T(digw, i41);
  loop_ub = w->size[0] * w->size[1];
  for (i41 = 0; i41 < loop_ub; i41++) {
    digw->data[i41] = 6.2831853071795862 * w->data[i41] / Fs;
  }

  emxInit_creal_T(&s, 2);

  /*  Convert from Hz to rad/sample for computational purposes */
  i41 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = digw->size[1];
  emxEnsureCapacity_creal_T(s, i41);
  loop_ub = digw->size[0] * digw->size[1];
  for (i41 = 0; i41 < loop_ub; i41++) {
    s->data[i41].re = digw->data[i41] * 0.0;
    s->data[i41].im = digw->data[i41];
  }

  emxInit_creal_T(&y, 2);
  c_exp(s);

  /*  Digital frequency must be used for this calculation */
  i41 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = s->size[1];
  emxEnsureCapacity_creal_T(y, i41);
  b9 = (y->size[1] == 0);
  if (!b9) {
    i41 = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity_creal_T(y, i41);
    loop_ub = y->size[1];
    for (i41 = 0; i41 < loop_ub; i41++) {
      y->data[i41].re = b[0];
      y->data[i41].im = 0.0;
    }

    for (k = 0; k < 84; k++) {
      i41 = s->size[0] * s->size[1];
      loop_ub = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = s->size[1];
      emxEnsureCapacity_creal_T(y, loop_ub);
      b_re = b[k + 1];
      loop_ub = i41 - 1;
      for (i41 = 0; i41 <= loop_ub; i41++) {
        s_re = s->data[i41].re * y->data[i41].re - s->data[i41].im * y->data[i41]
          .im;
        s_im = s->data[i41].re * y->data[i41].im + s->data[i41].im * y->data[i41]
          .re;
        y->data[i41].re = s_re + b_re;
        y->data[i41].im = s_im;
      }
    }
  }

  i41 = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = digw->size[1];
  emxEnsureCapacity_creal_T(s, i41);
  loop_ub = digw->size[0] * digw->size[1];
  for (i41 = 0; i41 < loop_ub; i41++) {
    b_re = digw->data[i41] * 0.0;
    s_re = digw->data[i41];
    s->data[i41].re = 84.0 * b_re;
    s->data[i41].im = 84.0 * s_re;
  }

  emxFree_real_T(&digw);
  c_exp(s);
  rdivide_helper(y, s, hh);

  /*  Generate the default structure to pass to freqzplot */
  /*  If rad/sample, Fs is empty */
  emxFree_creal_T(&y);
  emxFree_creal_T(&s);
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
  emxEnsureCapacity_creal_T(c, k);
  c->data[0].re = 1.0;
  c->data[0].im = 0.0;
  for (unnamed_idx_1 = 0; unnamed_idx_1 < n; unnamed_idx_1++) {
    x_re = -x->data[unnamed_idx_1].re;
    x_im = -x->data[unnamed_idx_1].im;
    c_re = c->data[unnamed_idx_1].re;
    c_im = c->data[unnamed_idx_1].im;
    c->data[unnamed_idx_1 + 1].re = x_re * c_re - x_im * c_im;
    c->data[unnamed_idx_1 + 1].im = x_re * c_im + x_im * c_re;
    for (k = unnamed_idx_1 + 1; k >= 2; k--) {
      x_re = x->data[unnamed_idx_1].re * c->data[k - 2].re - x->
        data[unnamed_idx_1].im * c->data[k - 2].im;
      x_im = x->data[unnamed_idx_1].re * c->data[k - 2].im + x->
        data[unnamed_idx_1].im * c->data[k - 2].re;
      c->data[k - 1].re -= x_re;
      c->data[k - 1].im -= x_im;
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
  y = fmax(a, b);
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
  int i67;
  int i;
  int im1n_tmp;
  int in;
  creal_T alpha1;
  int n_tmp;
  int lastc;
  double xnorm;
  double beta1;
  int iv0_tmp;
  bool b_tau;
  int knt;
  double beta1_im;
  int i68;
  int lastv;
  int b_lastc;
  bool exitg1;
  int k;
  creal_T c;
  double temp_im;
  int ix;
  int exitg5;
  int exitg2;
  int i69;
  int exitg4;
  n = a->size[0];
  if (a->size[0] < 1) {
    ntau = 0;
  } else {
    ntau = a->size[0] - 1;
  }

  emxInit_creal_T(&tau, 1);
  emxInit_creal_T(&work, 1);
  i67 = tau->size[0];
  tau->size[0] = ntau;
  emxEnsureCapacity_creal_T(tau, i67);
  ntau = a->size[0];
  i67 = work->size[0];
  work->size[0] = ntau;
  emxEnsureCapacity_creal_T(work, i67);
  for (i67 = 0; i67 < ntau; i67++) {
    work->data[i67].re = 0.0;
    work->data[i67].im = 0.0;
  }

  i67 = a->size[0];
  for (i = 0; i <= i67 - 2; i++) {
    im1n_tmp = i * n;
    in = (i + 1) * n;
    alpha1 = a->data[(i + a->size[0] * i) + 1];
    ntau = i + 3;
    if (ntau >= n) {
      ntau = n;
    }

    ntau += im1n_tmp;
    n_tmp = n - i;
    lastc = n_tmp - 2;
    tau->data[i].re = 0.0;
    tau->data[i].im = 0.0;
    if (lastc + 1 > 0) {
      xnorm = xnrm2(lastc, a, ntau);
      if ((xnorm != 0.0) || (a->data[(i + a->size[0] * i) + 1].im != 0.0)) {
        beta1 = xdlapy3(a->data[(i + a->size[0] * i) + 1].re, a->data[(i +
          a->size[0] * i) + 1].im, xnorm);
        if (a->data[(i + a->size[0] * i) + 1].re >= 0.0) {
          beta1 = -beta1;
        }

        if (fabs(beta1) < 1.0020841800044864E-292) {
          knt = -1;
          i68 = (ntau + lastc) - 1;
          do {
            knt++;
            for (k = ntau; k <= i68; k++) {
              xnorm = a->data[k - 1].re;
              beta1_im = a->data[k - 1].im;
              a->data[k - 1].re = 9.9792015476736E+291 * xnorm - 0.0 * beta1_im;
              a->data[k - 1].im = 9.9792015476736E+291 * beta1_im + 0.0 * xnorm;
            }

            beta1 *= 9.9792015476736E+291;
            alpha1.re *= 9.9792015476736E+291;
            alpha1.im *= 9.9792015476736E+291;
          } while (!(fabs(beta1) >= 1.0020841800044864E-292));

          beta1 = xdlapy3(alpha1.re, alpha1.im, xnrm2(lastc, a, ntau));
          if (alpha1.re >= 0.0) {
            beta1 = -beta1;
          }

          xnorm = beta1 - alpha1.re;
          if (0.0 - alpha1.im == 0.0) {
            tau->data[i].re = xnorm / beta1;
            tau->data[i].im = 0.0;
          } else if (xnorm == 0.0) {
            tau->data[i].re = 0.0;
            tau->data[i].im = (0.0 - alpha1.im) / beta1;
          } else {
            tau->data[i].re = xnorm / beta1;
            tau->data[i].im = (0.0 - alpha1.im) / beta1;
          }

          c.re = alpha1.re - beta1;
          c.im = alpha1.im;
          xscal(lastc, recip(c), a, ntau);
          for (k = 0; k <= knt; k++) {
            beta1 *= 1.0020841800044864E-292;
          }

          alpha1.re = beta1;
          alpha1.im = 0.0;
        } else {
          xnorm = beta1 - a->data[(i + a->size[0] * i) + 1].re;
          beta1_im = 0.0 - a->data[(i + a->size[0] * i) + 1].im;
          if (beta1_im == 0.0) {
            tau->data[i].re = xnorm / beta1;
            tau->data[i].im = 0.0;
          } else if (xnorm == 0.0) {
            tau->data[i].re = 0.0;
            tau->data[i].im = beta1_im / beta1;
          } else {
            tau->data[i].re = xnorm / beta1;
            tau->data[i].im = beta1_im / beta1;
          }

          c.re = a->data[(i + a->size[0] * i) + 1].re - beta1;
          c.im = a->data[(i + a->size[0] * i) + 1].im;
          xscal(lastc, recip(c), a, ntau);
          alpha1.re = beta1;
          alpha1.im = 0.0;
        }
      }
    }

    a->data[(i + a->size[0] * i) + 1].re = 1.0;
    a->data[(i + a->size[0] * i) + 1].im = 0.0;
    n_tmp -= 3;
    iv0_tmp = (i + im1n_tmp) + 2;
    im1n_tmp = in + 1;
    b_tau = ((tau->data[i].re != 0.0) || (tau->data[i].im != 0.0));
    if (b_tau) {
      lastv = n_tmp + 2;
      ntau = iv0_tmp + n_tmp;
      exitg1 = false;
      while ((!exitg1) && (lastv > 0)) {
        b_tau = ((a->data[ntau].re == 0.0) && (a->data[ntau].im == 0.0));
        if (b_tau) {
          lastv--;
          ntau--;
        } else {
          exitg1 = true;
        }
      }

      b_lastc = n;
      exitg1 = false;
      while ((!exitg1) && (b_lastc > 0)) {
        ntau = in + b_lastc;
        knt = ntau + (lastv - 1) * n;
        do {
          exitg2 = 0;
          if (((n > 0) && (ntau <= knt)) || ((n < 0) && (ntau >= knt))) {
            b_tau = ((a->data[ntau - 1].re != 0.0) || (a->data[ntau - 1].im !=
                      0.0));
            if (b_tau) {
              exitg2 = 1;
            } else {
              ntau += n;
            }
          } else {
            b_lastc--;
            exitg2 = 2;
          }
        } while (exitg2 == 0);

        if (exitg2 == 1) {
          exitg1 = true;
        }
      }
    } else {
      lastv = 0;
      b_lastc = 0;
    }

    if (lastv > 0) {
      if (b_lastc != 0) {
        for (ntau = 0; ntau < b_lastc; ntau++) {
          work->data[ntau].re = 0.0;
          work->data[ntau].im = 0.0;
        }

        ix = iv0_tmp;
        i68 = (in + n * (lastv - 1)) + 1;
        for (knt = im1n_tmp; n < 0 ? knt >= i68 : knt <= i68; knt += n) {
          c.re = a->data[ix - 1].re - 0.0 * a->data[ix - 1].im;
          c.im = a->data[ix - 1].im + 0.0 * a->data[ix - 1].re;
          ntau = 0;
          i69 = (knt + b_lastc) - 1;
          for (k = knt; k <= i69; k++) {
            xnorm = a->data[k - 1].re * c.re - a->data[k - 1].im * c.im;
            beta1_im = a->data[k - 1].re * c.im + a->data[k - 1].im * c.re;
            work->data[ntau].re += xnorm;
            work->data[ntau].im += beta1_im;
            ntau++;
          }

          ix++;
        }
      }

      c.re = -tau->data[i].re;
      c.im = -tau->data[i].im;
      if ((!(c.re == 0.0)) || (!(c.im == 0.0))) {
        ntau = in;
        knt = iv0_tmp - 1;
        for (k = 0; k < lastv; k++) {
          b_tau = ((a->data[knt].re != 0.0) || (a->data[knt].im != 0.0));
          if (b_tau) {
            beta1 = a->data[knt].re * c.re + a->data[knt].im * c.im;
            temp_im = a->data[knt].re * c.im - a->data[knt].im * c.re;
            ix = 0;
            i68 = ntau + 1;
            i69 = b_lastc + ntau;
            for (im1n_tmp = i68; im1n_tmp <= i69; im1n_tmp++) {
              xnorm = work->data[ix].re * beta1 - work->data[ix].im * temp_im;
              beta1_im = work->data[ix].re * temp_im + work->data[ix].im * beta1;
              a->data[im1n_tmp - 1].re += xnorm;
              a->data[im1n_tmp - 1].im += beta1_im;
              ix++;
            }
          }

          knt++;
          ntau += n;
        }
      }
    }

    im1n_tmp = (i + in) + 2;
    beta1 = tau->data[i].re;
    temp_im = -tau->data[i].im;
    if ((beta1 != 0.0) || (temp_im != 0.0)) {
      lastv = n_tmp + 2;
      ntau = iv0_tmp + n_tmp;
      do {
        exitg5 = 0;
        if (lastv > 0) {
          b_tau = ((a->data[ntau].re == 0.0) && (a->data[ntau].im == 0.0));
          if (b_tau) {
            lastv--;
            ntau--;
          } else {
            exitg5 = 1;
          }
        } else {
          exitg5 = 2;
        }
      } while (exitg5 == 0);

      do {
        exitg4 = 0;
        if (lastc + 1 > 0) {
          ntau = im1n_tmp + lastc * n;
          k = ntau;
          do {
            exitg2 = 0;
            if (k <= (ntau + lastv) - 1) {
              b_tau = ((a->data[k - 1].re != 0.0) || (a->data[k - 1].im != 0.0));
              if (b_tau) {
                exitg2 = 1;
              } else {
                k++;
              }
            } else {
              lastc--;
              exitg2 = 2;
            }
          } while (exitg2 == 0);

          if (exitg2 == 1) {
            exitg4 = 1;
          }
        } else {
          exitg4 = 1;
        }
      } while (exitg4 == 0);
    } else {
      lastv = 0;
      lastc = -1;
    }

    if (lastv > 0) {
      if (lastc + 1 != 0) {
        for (ntau = 0; ntau <= lastc; ntau++) {
          work->data[ntau].re = 0.0;
          work->data[ntau].im = 0.0;
        }

        ntau = 0;
        i68 = im1n_tmp + n * lastc;
        for (knt = im1n_tmp; n < 0 ? knt >= i68 : knt <= i68; knt += n) {
          ix = iv0_tmp - 1;
          c.re = 0.0;
          c.im = 0.0;
          i69 = (knt + lastv) - 1;
          for (k = knt; k <= i69; k++) {
            c.re += a->data[k - 1].re * a->data[ix].re + a->data[k - 1].im *
              a->data[ix].im;
            c.im += a->data[k - 1].re * a->data[ix].im - a->data[k - 1].im *
              a->data[ix].re;
            ix++;
          }

          work->data[ntau].re += c.re - 0.0 * c.im;
          work->data[ntau].im += c.im + 0.0 * c.re;
          ntau++;
        }
      }

      c.re = -beta1;
      c.im = -temp_im;
      if ((!(-beta1 == 0.0)) || (!(-temp_im == 0.0))) {
        ntau = im1n_tmp - 1;
        knt = 0;
        for (k = 0; k <= lastc; k++) {
          b_tau = ((work->data[knt].re != 0.0) || (work->data[knt].im != 0.0));
          if (b_tau) {
            beta1 = work->data[knt].re * c.re + work->data[knt].im * c.im;
            temp_im = work->data[knt].re * c.im - work->data[knt].im * c.re;
            ix = iv0_tmp;
            i68 = ntau + 1;
            i69 = lastv + ntau;
            for (im1n_tmp = i68; im1n_tmp <= i69; im1n_tmp++) {
              xnorm = a->data[ix - 1].re * beta1 - a->data[ix - 1].im * temp_im;
              beta1_im = a->data[ix - 1].re * temp_im + a->data[ix - 1].im *
                beta1;
              a->data[im1n_tmp - 1].re += xnorm;
              a->data[im1n_tmp - 1].im += beta1_im;
              ix++;
            }
          }

          knt++;
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
  if (n >= 1) {
    if (n == 1) {
      y = rt_hypotd_snf(x->data[ix0 - 1].re, x->data[ix0 - 1].im);
    } else {
      scale = 3.3121686421112381E-170;
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
  int i70;
  int k;
  double x_re;
  double x_im;
  i70 = (ix0 + n) - 1;
  for (k = ix0; k <= i70; k++) {
    x_re = x->data[k - 1].re;
    x_im = x->data[k - 1].im;
    x->data[k - 1].re = a.re * x_re - a.im * x_im;
    x->data[k - 1].im = a.re * x_im + a.im * x_re;
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
  int nzcount;
  int ii;
  double anrm;
  bool exitg1;
  double absxk;
  bool ilascl;
  double anrmto;
  int ilo;
  double ctoc;
  int ihi;
  bool notdone;
  int exitg3;
  double cfrom1;
  int i;
  int n;
  double cto1;
  int j;
  double At_re;
  int jrow;
  creal_T b_At;
  creal_T c_At;
  bool exitg4;
  double c;
  creal_T atmp;
  int exitg2;
  bool d_At;
  double stemp_re;
  emxInit_creal_T(&At, 2);
  nzcount = At->size[0] * At->size[1];
  At->size[0] = A->size[0];
  At->size[1] = A->size[1];
  emxEnsureCapacity_creal_T(At, nzcount);
  ii = A->size[0] * A->size[1];
  for (nzcount = 0; nzcount < ii; nzcount++) {
    At->data[nzcount] = A->data[nzcount];
  }

  *info = 0;
  nzcount = alpha1->size[0];
  alpha1->size[0] = At->size[0];
  emxEnsureCapacity_creal_T(alpha1, nzcount);
  ii = At->size[0];
  for (nzcount = 0; nzcount < ii; nzcount++) {
    alpha1->data[nzcount].re = 0.0;
    alpha1->data[nzcount].im = 0.0;
  }

  nzcount = beta1->size[0];
  beta1->size[0] = At->size[0];
  emxEnsureCapacity_creal_T(beta1, nzcount);
  ii = At->size[0];
  for (nzcount = 0; nzcount < ii; nzcount++) {
    beta1->data[nzcount].re = 0.0;
    beta1->data[nzcount].im = 0.0;
  }

  if ((At->size[0] != 0) && (At->size[1] != 0)) {
    anrm = 0.0;
    ii = 0;
    exitg1 = false;
    while ((!exitg1) && (ii <= At->size[0] * At->size[1] - 1)) {
      absxk = rt_hypotd_snf(At->data[ii].re, At->data[ii].im);
      if (rtIsNaN(absxk)) {
        anrm = rtNaN;
        exitg1 = true;
      } else {
        if (absxk > anrm) {
          anrm = absxk;
        }

        ii++;
      }
    }

    if (rtIsInf(anrm) || rtIsNaN(anrm)) {
      nzcount = alpha1->size[0];
      alpha1->size[0] = At->size[0];
      emxEnsureCapacity_creal_T(alpha1, nzcount);
      ii = At->size[0];
      for (nzcount = 0; nzcount < ii; nzcount++) {
        alpha1->data[nzcount].re = rtNaN;
        alpha1->data[nzcount].im = 0.0;
      }

      nzcount = beta1->size[0];
      beta1->size[0] = At->size[0];
      emxEnsureCapacity_creal_T(beta1, nzcount);
      ii = At->size[0];
      for (nzcount = 0; nzcount < ii; nzcount++) {
        beta1->data[nzcount].re = rtNaN;
        beta1->data[nzcount].im = 0.0;
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
            At_re = 2.0041683600089728E-292;
            absxk = cfrom1;
          } else if (cto1 > absxk) {
            At_re = 4.9896007738368E+291;
            ctoc = cto1;
          } else {
            At_re = ctoc / absxk;
            notdone = false;
          }

          nzcount = At->size[0] * At->size[1];
          ii = At->size[0] * At->size[1];
          emxEnsureCapacity_creal_T(At, ii);
          ii = nzcount - 1;
          for (nzcount = 0; nzcount <= ii; nzcount++) {
            At->data[nzcount].re *= At_re;
            At->data[nzcount].im *= At_re;
          }
        }
      }

      ilo = 1;
      ihi = At->size[0];
      if (At->size[0] <= 1) {
        ihi = 1;
      } else {
        do {
          exitg3 = 0;
          i = 0;
          j = 0;
          notdone = false;
          ii = ihi;
          exitg1 = false;
          while ((!exitg1) && (ii > 0)) {
            nzcount = 0;
            i = ii;
            j = ihi;
            jrow = 0;
            exitg4 = false;
            while ((!exitg4) && (jrow <= ihi - 1)) {
              d_At = ((At->data[(ii + At->size[0] * jrow) - 1].re != 0.0) ||
                      (At->data[(ii + At->size[0] * jrow) - 1].im != 0.0));
              if (d_At || (ii == jrow + 1)) {
                if (nzcount == 0) {
                  j = jrow + 1;
                  nzcount = 1;
                  jrow++;
                } else {
                  nzcount = 2;
                  exitg4 = true;
                }
              } else {
                jrow++;
              }
            }

            if (nzcount < 2) {
              notdone = true;
              exitg1 = true;
            } else {
              ii--;
            }
          }

          if (!notdone) {
            exitg3 = 2;
          } else {
            n = At->size[0];
            if (i != ihi) {
              for (ii = 1; ii <= n; ii++) {
                atmp = At->data[(i + At->size[0] * (ii - 1)) - 1];
                At->data[(i + At->size[0] * (ii - 1)) - 1] = At->data[(ihi +
                  At->size[0] * (ii - 1)) - 1];
                At->data[(ihi + At->size[0] * (ii - 1)) - 1] = atmp;
              }
            }

            if (j != ihi) {
              for (ii = 0; ii < ihi; ii++) {
                atmp = At->data[ii + At->size[0] * (j - 1)];
                At->data[ii + At->size[0] * (j - 1)] = At->data[ii + At->size[0]
                  * (ihi - 1)];
                At->data[ii + At->size[0] * (ihi - 1)] = atmp;
              }
            }

            ihi--;
            if (ihi == 1) {
              exitg3 = 1;
            }
          }
        } while (exitg3 == 0);

        if (exitg3 == 1) {
        } else {
          do {
            exitg2 = 0;
            i = 0;
            j = 0;
            notdone = false;
            jrow = ilo;
            exitg1 = false;
            while ((!exitg1) && (jrow <= ihi)) {
              nzcount = 0;
              i = ihi;
              j = jrow;
              ii = ilo;
              exitg4 = false;
              while ((!exitg4) && (ii <= ihi)) {
                d_At = ((At->data[(ii + At->size[0] * (jrow - 1)) - 1].re != 0.0)
                        || (At->data[(ii + At->size[0] * (jrow - 1)) - 1].im !=
                            0.0));
                if (d_At || (ii == jrow)) {
                  if (nzcount == 0) {
                    i = ii;
                    nzcount = 1;
                    ii++;
                  } else {
                    nzcount = 2;
                    exitg4 = true;
                  }
                } else {
                  ii++;
                }
              }

              if (nzcount < 2) {
                notdone = true;
                exitg1 = true;
              } else {
                jrow++;
              }
            }

            if (!notdone) {
              exitg2 = 1;
            } else {
              n = At->size[0];
              if (i != ilo) {
                for (ii = ilo; ii <= n; ii++) {
                  atmp = At->data[(i + At->size[0] * (ii - 1)) - 1];
                  At->data[(i + At->size[0] * (ii - 1)) - 1] = At->data[(ilo +
                    At->size[0] * (ii - 1)) - 1];
                  At->data[(ilo + At->size[0] * (ii - 1)) - 1] = atmp;
                }
              }

              if (j != ilo) {
                for (ii = 0; ii < ihi; ii++) {
                  atmp = At->data[ii + At->size[0] * (j - 1)];
                  At->data[ii + At->size[0] * (j - 1)] = At->data[ii + At->size
                    [0] * (ilo - 1)];
                  At->data[ii + At->size[0] * (ilo - 1)] = atmp;
                }
              }

              ilo++;
              if (ilo == ihi) {
                exitg2 = 1;
              }
            }
          } while (exitg2 == 0);
        }
      }

      n = At->size[0];
      if ((At->size[0] > 1) && (ihi >= ilo + 2)) {
        for (ii = ilo - 1; ii + 1 < ihi - 1; ii++) {
          nzcount = ii + 2;
          for (jrow = ihi - 1; jrow + 1 > ii + 2; jrow--) {
            b_At = At->data[(jrow + At->size[0] * ii) - 1];
            c_At = At->data[jrow + At->size[0] * ii];
            xzlartg(b_At, c_At, &c, &atmp, &At->data[(jrow + At->size[0] * ii) -
                    1]);
            At->data[jrow + At->size[0] * ii].re = 0.0;
            At->data[jrow + At->size[0] * ii].im = 0.0;
            for (j = nzcount; j <= n; j++) {
              absxk = atmp.re * At->data[jrow + At->size[0] * (j - 1)].re -
                atmp.im * At->data[jrow + At->size[0] * (j - 1)].im;
              ctoc = atmp.re * At->data[jrow + At->size[0] * (j - 1)].im +
                atmp.im * At->data[jrow + At->size[0] * (j - 1)].re;
              stemp_re = c * At->data[(jrow + At->size[0] * (j - 1)) - 1].re +
                absxk;
              absxk = c * At->data[(jrow + At->size[0] * (j - 1)) - 1].im + ctoc;
              ctoc = At->data[(jrow + At->size[0] * (j - 1)) - 1].re;
              cfrom1 = At->data[(jrow + At->size[0] * (j - 1)) - 1].im;
              cto1 = At->data[(jrow + At->size[0] * (j - 1)) - 1].im;
              At_re = At->data[(jrow + At->size[0] * (j - 1)) - 1].re;
              At->data[jrow + At->size[0] * (j - 1)].re = c * At->data[jrow +
                At->size[0] * (j - 1)].re - (atmp.re * ctoc + atmp.im * cfrom1);
              At->data[jrow + At->size[0] * (j - 1)].im = c * At->data[jrow +
                At->size[0] * (j - 1)].im - (atmp.re * cto1 - atmp.im * At_re);
              At->data[(jrow + At->size[0] * (j - 1)) - 1].re = stemp_re;
              At->data[(jrow + At->size[0] * (j - 1)) - 1].im = absxk;
            }

            atmp.re = -atmp.re;
            atmp.im = -atmp.im;
            for (i = 1; i <= ihi; i++) {
              absxk = atmp.re * At->data[(i + At->size[0] * (jrow - 1)) - 1].re
                - atmp.im * At->data[(i + At->size[0] * (jrow - 1)) - 1].im;
              ctoc = atmp.re * At->data[(i + At->size[0] * (jrow - 1)) - 1].im +
                atmp.im * At->data[(i + At->size[0] * (jrow - 1)) - 1].re;
              stemp_re = c * At->data[(i + At->size[0] * jrow) - 1].re + absxk;
              absxk = c * At->data[(i + At->size[0] * jrow) - 1].im + ctoc;
              ctoc = At->data[(i + At->size[0] * jrow) - 1].re;
              cfrom1 = At->data[(i + At->size[0] * jrow) - 1].im;
              cto1 = At->data[(i + At->size[0] * jrow) - 1].im;
              At_re = At->data[(i + At->size[0] * jrow) - 1].re;
              At->data[(i + At->size[0] * (jrow - 1)) - 1].re = c * At->data[(i
                + At->size[0] * (jrow - 1)) - 1].re - (atmp.re * ctoc + atmp.im *
                cfrom1);
              At->data[(i + At->size[0] * (jrow - 1)) - 1].im = c * At->data[(i
                + At->size[0] * (jrow - 1)) - 1].im - (atmp.re * cto1 - atmp.im *
                At_re);
              At->data[(i + At->size[0] * jrow) - 1].re = stemp_re;
              At->data[(i + At->size[0] * jrow) - 1].im = absxk;
            }
          }
        }
      }

      xzhgeqz(At, ilo, ihi, info, alpha1, beta1);
      if ((*info == 0) && ilascl) {
        notdone = true;
        while (notdone) {
          cfrom1 = anrmto * 2.0041683600089728E-292;
          cto1 = anrm / 4.9896007738368E+291;
          if ((cfrom1 > anrm) && (anrm != 0.0)) {
            At_re = 2.0041683600089728E-292;
            anrmto = cfrom1;
          } else if (cto1 > anrmto) {
            At_re = 4.9896007738368E+291;
            anrm = cto1;
          } else {
            At_re = anrm / anrmto;
            notdone = false;
          }

          nzcount = alpha1->size[0];
          emxEnsureCapacity_creal_T(alpha1, nzcount);
          ii = alpha1->size[0];
          for (nzcount = 0; nzcount < ii; nzcount++) {
            alpha1->data[nzcount].re *= At_re;
            alpha1->data[nzcount].im *= At_re;
          }
        }
      }
    }
  }

  emxFree_creal_T(&At);
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
  int n;
  double eshift_re;
  double eshift_im;
  creal_T ctemp;
  double anorm;
  double scale;
  double reAij;
  double sumsq;
  double b_atol;
  bool firstNonZero;
  int j;
  double ascale;
  double bscale;
  double imAij;
  bool guard1 = false;
  bool guard2 = false;
  double temp2;
  int ifirst;
  int istart;
  int ilast;
  int ilastm1;
  int ifrstm;
  int ilastm;
  int iiter;
  bool goto60;
  bool goto70;
  bool goto90;
  int jiter;
  int exitg1;
  bool b_guard1 = false;
  bool guard3 = false;
  bool exitg2;
  creal_T b_ascale;
  creal_T shift;
  creal_T c_A;
  double ad22_re;
  double ad22_im;
  double t1_im;
  emxInit_creal_T(&b_A, 2);
  jm1 = b_A->size[0] * b_A->size[1];
  b_A->size[0] = A->size[0];
  b_A->size[1] = A->size[1];
  emxEnsureCapacity_creal_T(b_A, jm1);
  jp1 = A->size[0] * A->size[1];
  for (jm1 = 0; jm1 < jp1; jm1++) {
    b_A->data[jm1] = A->data[jm1];
  }

  *info = -1;
  if ((A->size[0] == 1) && (A->size[1] == 1)) {
    ihi = 1;
  }

  n = A->size[0];
  jm1 = alpha1->size[0];
  alpha1->size[0] = A->size[0];
  emxEnsureCapacity_creal_T(alpha1, jm1);
  jp1 = A->size[0];
  for (jm1 = 0; jm1 < jp1; jm1++) {
    alpha1->data[jm1].re = 0.0;
    alpha1->data[jm1].im = 0.0;
  }

  jm1 = beta1->size[0];
  beta1->size[0] = A->size[0];
  emxEnsureCapacity_creal_T(beta1, jm1);
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
  if (ilo <= ihi) {
    scale = 0.0;
    sumsq = 0.0;
    firstNonZero = true;
    for (j = ilo; j <= ihi; j++) {
      jm1 = j + 1;
      if (ihi < j + 1) {
        jm1 = ihi;
      }

      for (jp1 = ilo; jp1 <= jm1; jp1++) {
        reAij = A->data[(jp1 + A->size[0] * (j - 1)) - 1].re;
        imAij = A->data[(jp1 + A->size[0] * (j - 1)) - 1].im;
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
  firstNonZero = true;
  jm1 = ihi + 1;
  for (j = jm1; j <= n; j++) {
    alpha1->data[j - 1] = A->data[(j + A->size[0] * (j - 1)) - 1];
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
    jiter = 0;
    do {
      exitg1 = 0;
      if (jiter <= 30 * ((ihi - ilo) + 1) - 1) {
        b_guard1 = false;
        if (ilast + 1 == ilo) {
          goto60 = true;
          b_guard1 = true;
        } else if (fabs(b_A->data[ilast + b_A->size[0] * ilastm1].re) + fabs
                   (b_A->data[ilast + b_A->size[0] * ilastm1].im) <= b_atol) {
          b_A->data[ilast + b_A->size[0] * ilastm1].re = 0.0;
          b_A->data[ilast + b_A->size[0] * ilastm1].im = 0.0;
          goto60 = true;
          b_guard1 = true;
        } else {
          j = ilastm1;
          guard3 = false;
          exitg2 = false;
          while ((!exitg2) && (j + 1 >= ilo)) {
            if (j + 1 == ilo) {
              guard3 = true;
              exitg2 = true;
            } else if (fabs(b_A->data[j + b_A->size[0] * (j - 1)].re) + fabs
                       (b_A->data[j + b_A->size[0] * (j - 1)].im) <= b_atol) {
              b_A->data[j + b_A->size[0] * (j - 1)].re = 0.0;
              b_A->data[j + b_A->size[0] * (j - 1)].im = 0.0;
              guard3 = true;
              exitg2 = true;
            } else {
              j--;
              guard3 = false;
            }
          }

          if (guard3) {
            ifirst = j + 1;
            goto70 = true;
          }

          if (goto70) {
            b_guard1 = true;
          } else {
            jp1 = alpha1->size[0];
            jm1 = alpha1->size[0];
            alpha1->size[0] = jp1;
            emxEnsureCapacity_creal_T(alpha1, jm1);
            for (jm1 = 0; jm1 < jp1; jm1++) {
              alpha1->data[jm1].re = rtNaN;
              alpha1->data[jm1].im = 0.0;
            }

            jp1 = beta1->size[0];
            jm1 = beta1->size[0];
            beta1->size[0] = jp1;
            emxEnsureCapacity_creal_T(beta1, jm1);
            for (jm1 = 0; jm1 < jp1; jm1++) {
              beta1->data[jm1].re = rtNaN;
              beta1->data[jm1].im = 0.0;
            }

            *info = 0;
            exitg1 = 1;
          }
        }

        if (b_guard1) {
          if (goto60) {
            goto60 = false;
            alpha1->data[ilast] = b_A->data[ilast + b_A->size[0] * ilast];
            ilast = ilastm1;
            ilastm1--;
            if (ilast + 1 < ilo) {
              firstNonZero = false;
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

              jiter++;
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
                c_sqrt(&shift);
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
                  b_ascale = b_A->data[(j + b_A->size[0] * jm1) - 1];
                  c_A = b_A->data[j + b_A->size[0] * jm1];
                  xzlartg(b_ascale, c_A, &imAij, &shift, &b_A->data[(j +
                           b_A->size[0] * jm1) - 1]);
                  b_A->data[j + b_A->size[0] * jm1].re = 0.0;
                  b_A->data[j + b_A->size[0] * jm1].im = 0.0;
                }

                for (n = j; n <= ilastm; n++) {
                  anorm = shift.re * b_A->data[j + b_A->size[0] * (n - 1)].re -
                    shift.im * b_A->data[j + b_A->size[0] * (n - 1)].im;
                  reAij = shift.re * b_A->data[j + b_A->size[0] * (n - 1)].im +
                    shift.im * b_A->data[j + b_A->size[0] * (n - 1)].re;
                  ad22_re = imAij * b_A->data[(j + b_A->size[0] * (n - 1)) - 1].
                    re + anorm;
                  ad22_im = imAij * b_A->data[(j + b_A->size[0] * (n - 1)) - 1].
                    im + reAij;
                  anorm = b_A->data[(j + b_A->size[0] * (n - 1)) - 1].re;
                  reAij = b_A->data[(j + b_A->size[0] * (n - 1)) - 1].im;
                  scale = b_A->data[(j + b_A->size[0] * (n - 1)) - 1].im;
                  sumsq = b_A->data[(j + b_A->size[0] * (n - 1)) - 1].re;
                  b_A->data[j + b_A->size[0] * (n - 1)].re = imAij * b_A->data[j
                    + b_A->size[0] * (n - 1)].re - (shift.re * anorm + shift.im *
                    reAij);
                  b_A->data[j + b_A->size[0] * (n - 1)].im = imAij * b_A->data[j
                    + b_A->size[0] * (n - 1)].im - (shift.re * scale - shift.im *
                    sumsq);
                  b_A->data[(j + b_A->size[0] * (n - 1)) - 1].re = ad22_re;
                  b_A->data[(j + b_A->size[0] * (n - 1)) - 1].im = ad22_im;
                }

                shift.re = -shift.re;
                shift.im = -shift.im;
                n = j;
                if (ilast + 1 < j + 2) {
                  n = ilast - 1;
                }

                for (jp1 = ifrstm; jp1 <= n + 2; jp1++) {
                  anorm = shift.re * b_A->data[(jp1 + b_A->size[0] * (j - 1)) -
                    1].re - shift.im * b_A->data[(jp1 + b_A->size[0] * (j - 1))
                    - 1].im;
                  reAij = shift.re * b_A->data[(jp1 + b_A->size[0] * (j - 1)) -
                    1].im + shift.im * b_A->data[(jp1 + b_A->size[0] * (j - 1))
                    - 1].re;
                  ad22_re = imAij * b_A->data[(jp1 + b_A->size[0] * j) - 1].re +
                    anorm;
                  ad22_im = imAij * b_A->data[(jp1 + b_A->size[0] * j) - 1].im +
                    reAij;
                  anorm = b_A->data[(jp1 + b_A->size[0] * j) - 1].re;
                  reAij = b_A->data[(jp1 + b_A->size[0] * j) - 1].im;
                  scale = b_A->data[(jp1 + b_A->size[0] * j) - 1].im;
                  sumsq = b_A->data[(jp1 + b_A->size[0] * j) - 1].re;
                  b_A->data[(jp1 + b_A->size[0] * (j - 1)) - 1].re = imAij *
                    b_A->data[(jp1 + b_A->size[0] * (j - 1)) - 1].re - (shift.re
                    * anorm + shift.im * reAij);
                  b_A->data[(jp1 + b_A->size[0] * (j - 1)) - 1].im = imAij *
                    b_A->data[(jp1 + b_A->size[0] * (j - 1)) - 1].im - (shift.re
                    * scale - shift.im * sumsq);
                  b_A->data[(jp1 + b_A->size[0] * j) - 1].re = ad22_re;
                  b_A->data[(jp1 + b_A->size[0] * j) - 1].im = ad22_im;
                }

                jm1 = j - 1;
                j++;
              }
            }

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
    if (firstNonZero) {
      *info = ilast;
      for (jp1 = 0; jp1 <= ilast; jp1++) {
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
    for (j = 0; j <= ilo - 2; j++) {
      alpha1->data[j] = b_A->data[j + b_A->size[0] * j];
    }
  }

  emxFree_creal_T(&b_A);
  (*info)++;
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
  double y_tmp;
  double scale;
  double b_y_tmp;
  double f2s;
  double f2;
  double fs_re;
  double fs_im;
  double gs_re;
  double gs_im;
  int count;
  int rescaledir;
  bool guard1 = false;
  double g2s;
  y_tmp = fabs(f.re);
  scale = y_tmp;
  b_y_tmp = fabs(f.im);
  if (b_y_tmp > y_tmp) {
    scale = b_y_tmp;
  }

  f2s = fabs(g.re);
  f2 = fabs(g.im);
  if (f2 > f2s) {
    f2s = f2;
  }

  if (f2s > scale) {
    scale = f2s;
  }

  fs_re = f.re;
  fs_im = f.im;
  gs_re = g.re;
  gs_im = g.im;
  count = -1;
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
    f2 = fs_re * fs_re + fs_im * fs_im;
    scale = gs_re * gs_re + gs_im * gs_im;
    f2s = scale;
    if (1.0 > scale) {
      f2s = 1.0;
    }

    if (f2 <= f2s * 2.0041683600089728E-292) {
      if ((f.re == 0.0) && (f.im == 0.0)) {
        *cs = 0.0;
        r->re = rt_hypotd_snf(g.re, g.im);
        r->im = 0.0;
        scale = rt_hypotd_snf(gs_re, gs_im);
        sn->re = gs_re / scale;
        sn->im = -gs_im / scale;
      } else {
        g2s = sqrt(scale);
        *cs = rt_hypotd_snf(fs_re, fs_im) / g2s;
        if (b_y_tmp > y_tmp) {
          y_tmp = b_y_tmp;
        }

        if (y_tmp > 1.0) {
          scale = rt_hypotd_snf(f.re, f.im);
          fs_re = f.re / scale;
          fs_im = f.im / scale;
        } else {
          f2 = 7.4428285367870146E+137 * f.re;
          f2s = 7.4428285367870146E+137 * f.im;
          scale = rt_hypotd_snf(f2, f2s);
          fs_re = f2 / scale;
          fs_im = f2s / scale;
        }

        gs_re /= g2s;
        gs_im = -gs_im / g2s;
        sn->re = fs_re * gs_re - fs_im * gs_im;
        sn->im = fs_re * gs_im + fs_im * gs_re;
        r->re = *cs * f.re + (sn->re * g.re - sn->im * g.im);
        r->im = *cs * f.im + (sn->re * g.im + sn->im * g.re);
      }
    } else {
      f2s = sqrt(1.0 + scale / f2);
      r->re = f2s * fs_re;
      r->im = f2s * fs_im;
      *cs = 1.0 / f2s;
      scale += f2;
      f2 = r->re / scale;
      f2s = r->im / scale;
      sn->re = f2 * gs_re - f2s * -gs_im;
      sn->im = f2 * -gs_im + f2s * gs_re;
      if (rescaledir > 0) {
        for (rescaledir = 0; rescaledir <= count; rescaledir++) {
          r->re *= 7.4428285367870146E+137;
          r->im *= 7.4428285367870146E+137;
        }
      } else {
        if (rescaledir < 0) {
          for (rescaledir = 0; rescaledir <= count; rescaledir++) {
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
 *                emxArray_real_T *c
 *                double *d
 * Return Type  : void
 */
static void zp2ss_cg(emxArray_creal_T *a, emxArray_real_T *b, emxArray_real_T *c,
                     double *d)
{
  int i6;
  creal_T b_c[3];
  int k;
  double d0;
  double Y[4];
  emxArray_int8_T *reshapes_f2;
  emxArray_cint8_T *result;
  int loop_ub;
  int i7;
  signed char sizes_idx_1;
  signed char input_sizes_idx_0;
  emxArray_real_T *b1;
  emxArray_real_T *varargin_2;
  double B[4];
  int result_re;
  int result_im;
  emxArray_real_T *r2;

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
  i6 = a->size[0] * a->size[1];
  a->size[0] = 1;
  a->size[1] = 1;
  emxEnsureCapacity_creal_T(a, i6);
  a->data[0].re = -1.0;
  a->data[0].im = 0.0;
  i6 = b->size[0];
  b->size[0] = 1;
  emxEnsureCapacity_real_T(b, i6);
  b->data[0] = 1.0;
  i6 = c->size[0] * c->size[1];
  c->size[0] = 1;
  c->size[1] = 1;
  emxEnsureCapacity_real_T(c, i6);
  c->data[0] = 1.0;

  /*  If odd number of zeros only, convert the zero at the */
  /*  end, along with a pole-pair into state-space. */
  /*    H(s) = (s+num(2))/(s^2+den(2)s+den(3)) */
  /*  Now we have an even number of poles and zeros, although not */
  /*  necessarily the same number - there may be more poles. */
  /*    H(s) = (s^2+num(2)s+num(3))/(s^2+den(2)s+den(3)) */
  /*  Loop through rest of pairs, connecting in series to build the model. */
  /*  Take care of any left over unmatched pole pairs. */
  /*    H(s) = 1/(s^2+den(2)s+den(3)) */
  b_c[1].re = 0.49999999999999978;
  for (k = 2; k >= 2; k--) {
    b_c[1].re -= -0.49999999999999978;
  }

  /*  Balancing transformation */
  d0 = (1.0 - -b_c[1].re * 0.0) / 1.0000000000000002;
  Y[1] = d0;
  Y[0] = -b_c[1].re - d0 * 0.0;
  Y[3] = 0.0;
  Y[2] = -0.99999999999999989;
  emxInit_int8_T(&reshapes_f2, 2);

  /*  [a,b,c,d] = series(a,b,c,d,a1,b1,c1,d1); */
  /*  Next lines perform series connection */
  i6 = reshapes_f2->size[0] * reshapes_f2->size[1];
  reshapes_f2->size[0] = 1;
  reshapes_f2->size[1] = 2;
  emxEnsureCapacity_int8_T(reshapes_f2, i6);
  for (i6 = 0; i6 < 2; i6++) {
    reshapes_f2->data[i6] = 0;
  }

  emxInit_cint8_T(&result, 2);
  i6 = result->size[0] * result->size[1];
  result->size[0] = 1;
  result->size[1] = 1 + reshapes_f2->size[1];
  emxEnsureCapacity_cint8_T(result, i6);
  for (i6 = 0; i6 < 1; i6++) {
    for (i7 = 0; i7 < 1; i7++) {
      sizes_idx_1 = (signed char)a->data[0].re;
      input_sizes_idx_0 = (signed char)a->data[0].im;
      result->data[0].re = sizes_idx_1;
      result->data[0].im = input_sizes_idx_0;
    }
  }

  loop_ub = reshapes_f2->size[1];
  for (i6 = 0; i6 < loop_ub; i6++) {
    k = reshapes_f2->size[0];
    for (i7 = 0; i7 < k; i7++) {
      result->data[i7 + result->size[0] * (i6 + 1)].re = reshapes_f2->data[i7 +
        reshapes_f2->size[0] * i6];
      result->data[i7 + result->size[0] * (i6 + 1)].im = 0;
    }
  }

  emxFree_int8_T(&reshapes_f2);
  emxInit_real_T(&b1, 2);
  i6 = b1->size[0] * b1->size[1];
  b1->size[0] = 2;
  b1->size[1] = c->size[1];
  emxEnsureCapacity_real_T(b1, i6);
  loop_ub = c->size[1];
  for (i6 = 0; i6 < loop_ub; i6++) {
    b1->data[i6 << 1] = c->data[i6];
  }

  loop_ub = c->size[1];
  for (i6 = 0; i6 < loop_ub; i6++) {
    b1->data[1 + (i6 << 1)] = 0.0 * c->data[i6];
  }

  for (i6 = 0; i6 < 2; i6++) {
    B[i6] = 0.0;
    d0 = Y[i6 + 2];
    B[i6] = Y[i6] + d0 * 0.0;
    B[i6 + 2] = 0.0;
    B[i6 + 2] = Y[i6] * 0.0 + d0 * 1.0000000000000002;
  }

  emxInit_real_T(&varargin_2, 2);
  i6 = varargin_2->size[0] * varargin_2->size[1];
  varargin_2->size[0] = 2;
  varargin_2->size[1] = b1->size[1] + 2;
  emxEnsureCapacity_real_T(varargin_2, i6);
  loop_ub = b1->size[1];
  for (i6 = 0; i6 < loop_ub; i6++) {
    i7 = i6 << 1;
    varargin_2->data[i7] = b1->data[i7];
    i7++;
    varargin_2->data[i7] = b1->data[i7];
  }

  varargin_2->data[b1->size[1] << 1] = B[0];
  varargin_2->data[1 + (b1->size[1] << 1)] = B[1];
  varargin_2->data[(1 + b1->size[1]) << 1] = B[2];
  varargin_2->data[1 + ((1 + b1->size[1]) << 1)] = B[3];
  emxFree_real_T(&b1);
  if (result->size[1] != 0) {
    sizes_idx_1 = (signed char)result->size[1];
  } else {
    sizes_idx_1 = 3;
  }

  input_sizes_idx_0 = (signed char)(result->size[1] != 0);
  k = input_sizes_idx_0;
  i6 = a->size[0] * a->size[1];
  a->size[0] = input_sizes_idx_0 + 2;
  a->size[1] = sizes_idx_1;
  emxEnsureCapacity_creal_T(a, i6);
  loop_ub = sizes_idx_1;
  for (i6 = 0; i6 < loop_ub; i6++) {
    for (i7 = 0; i7 < k; i7++) {
      result_re = result->data[input_sizes_idx_0 * i6].re;
      result_im = result->data[input_sizes_idx_0 * i6].im;
      a->data[a->size[0] * i6].re = result_re;
      a->data[a->size[0] * i6].im = result_im;
    }
  }

  emxFree_cint8_T(&result);
  loop_ub = sizes_idx_1;
  for (i6 = 0; i6 < loop_ub; i6++) {
    for (i7 = 0; i7 < 2; i7++) {
      a->data[(i7 + input_sizes_idx_0) + a->size[0] * i6].re = varargin_2->
        data[i7 + (i6 << 1)];
      a->data[(i7 + input_sizes_idx_0) + a->size[0] * i6].im = 0.0;
    }
  }

  emxFree_real_T(&varargin_2);
  i6 = b->size[0];
  i7 = b->size[0];
  b->size[0] = i6 + 2;
  emxEnsureCapacity_real_T(b, i7);
  b->data[i6] = 0.0;
  b->data[i6 + 1] = 0.0;
  emxInit_real_T(&r2, 2);
  i6 = r2->size[0] * r2->size[1];
  r2->size[0] = 1;
  r2->size[1] = c->size[1] + 2;
  emxEnsureCapacity_real_T(r2, i6);
  loop_ub = c->size[1];
  for (i6 = 0; i6 < loop_ub; i6++) {
    r2->data[i6] = 0.0 * c->data[i6];
  }

  r2->data[c->size[1]] = 0.0;
  r2->data[1 + c->size[1]] = 1.0000000000000002;
  i6 = c->size[0] * c->size[1];
  c->size[0] = 1;
  c->size[1] = r2->size[1];
  emxEnsureCapacity_real_T(c, i6);
  loop_ub = r2->size[0] * r2->size[1];
  for (i6 = 0; i6 < loop_ub; i6++) {
    c->data[i6] = r2->data[i6];
  }

  emxFree_real_T(&r2);

  /*  Apply gain k: */
  i6 = c->size[0] * c->size[1];
  i7 = c->size[0] * c->size[1];
  c->size[0] = 1;
  emxEnsureCapacity_real_T(c, i7);
  loop_ub = i6 - 1;
  for (i6 = 0; i6 <= loop_ub; i6++) {
    c->data[i6] *= 0.99999999999999989;
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
 *                double maxTaps
 *                short outputTaps[128]
 *                double *numOutputTaps
 *                double *filterGain
 *                double *Apass_actual
 *                double *Astop_actual
 * Return Type  : void
 */
void internal_design_filter_cg(double Rdata, double Fpass, double Fstop, double
  caldiv, double FIR, double HB1, double PLL_mult, double Apass, double Astop,
  double phEQ, double HB2, double HB3, const char Type[7], const char RxTx[2],
  double RFbw, double DAC_div, double converter_rate, double PLL_rate, double
  Fcenter, double wnom, double FIRdBmin, double int_FIR, double maxTaps, short
  outputTaps[128], double *numOutputTaps, double *filterGain, double
  *Apass_actual, double *Astop_actual)
{
  emxArray_creal_T *a1;
  emxArray_creal_T *a2;
  double b1[4];
  double b2[2];
  creal_T a2_data[2];
  int a2_size[2];
  int b1_size[2];
  double b1_data[4];
  int i0;
  int b2_size[2];
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
  static const double y[5] = { 0.0625, 0.25, 0.375, 0.25, 0.0625 };

  double dec_int3_coeff_data[29];
  static const double b_y[17] = { 0.00335693359375, 0.00506591796875, 0.0,
    -0.02398681640625, -0.035400390625, 0.0, 0.1168212890625, 0.24664306640625,
    0.3125, 0.24664306640625, 0.1168212890625, 0.0, -0.035400390625,
    -0.02398681640625, 0.0, 0.00506591796875, 0.00335693359375 };

  static const double c_y[29] = { 0.00146484375, -0.00077311197916666663, 0.0,
    -0.00634765625, -0.00048828125, 0.0, 0.019490559895833332,
    0.0090738932291666661, 0.0, -0.0494384765625, -0.0404052734375, 0.0,
    0.14522298177083331, 0.25541178385416663, 0.33333333333333331,
    0.25541178385416663, 0.14522298177083331, 0.0, -0.0404052734375,
    -0.0494384765625, 0.0, 0.0090738932291666661, 0.019490559895833332, 0.0,
    -0.00048828125, -0.00634765625, 0.0, -0.00077311197916666663, 0.00146484375
  };

  int nm1d2;
  int idx;
  double sigmax;
  signed char i1;
  char enables[4];
  double w[2048];
  double phi[2048];
  static const double hb2_coeff[7] = { -0.03515625, 0.0, 0.28515625, 0.5,
    0.28515625, 0.0, -0.03515625 };

  static creal_T combinedResponse[2048];
  static creal_T dcv0[2048];
  double b_combinedResponse[2048];
  double invariance[2048];
  double apnd;
  double b_phi[2048];
  double sigma;
  emxArray_real_T *fg;
  double b_w[2048];
  double clkFIR;
  double Gstop;
  double Gpass;
  int loop_ub;
  emxArray_real_T *omega;
  emxArray_creal_T *rg1;
  emxArray_creal_T *c_combinedResponse;
  int i2;
  emxArray_creal_T *r0;
  double rg1_re;
  double rg1_im;
  double re;
  double im;
  emxArray_real_T *F4;
  emxArray_real_T *fg2;
  int n;
  emxArray_real_T *sw;
  int k;
  emxArray_real_T *omega2;
  emxArray_real_T *F3;
  emxArray_creal_T *rgN;
  emxArray_real_T *b_omega2;
  emxArray_creal_T *d_combinedResponse;
  emxArray_real_T *weight;
  bool exitg1;
  emxArray_real_T *F1;
  emxArray_real_T *F2;
  emxArray_real_T *A1;
  emxArray_real_T *A2;
  emxArray_real_T *W1;
  emxArray_real_T *W2;
  emxArray_real_T *tap_store;
  emxArray_real_T *Apass_actual_vector;
  emxArray_real_T *Astop_actual_vector;
  unsigned int i;
  emxArray_real_T *ccoef;
  emxArray_real_T *b_ccoef;
  emxArray_real_T *b_F1;
  emxArray_real_T *b_W1;
  int exitg2;
  bool valid;
  double firTapsPreScale[128];
  double b_firTapsPreScale[128];
  short i3;
  (void)caldiv;
  (void)PLL_mult;
  (void)Type;
  (void)RFbw;
  (void)DAC_div;
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
    b1_data[0] = b2[0];
    b1_data[1] = b2[1];
    i0 = a1->size[0] * a1->size[1];
    a1->size[0] = 1;
    a1->size[1] = 2;
    emxEnsureCapacity_creal_T(a1, i0);
    a1->data[0] = a2_data[0];
    a1->data[1] = a2_data[1];

    /*  1st order */
    b_butter_cg(6.2831853071795862 * wnom, b1, a2);
    b2_size[0] = 1;
    b2_size[1] = 4;
    b2_data[0] = b1[0];
    b2_data[1] = b1[1];
    b2_data[2] = b1[2];
    b2_data[3] = b1[3];

    /*  3rd order */
    /*  Define the digital filters with fixed coefficients */
    memcpy(&hb1_coeff[0], &dv1[0], 15U * sizeof(double));
    hb3_coeff_size[0] = 1;
    hb3_coeff_size[1] = 5;
    for (i0 = 0; i0 < 5; i0++) {
      hb3_coeff_data[i0] = y[i0];
    }

    dec_int3_coeff_size[0] = 1;
    dec_int3_coeff_size[1] = 17;
    memcpy(&dec_int3_coeff_data[0], &b_y[0], 17U * sizeof(double));
  } else {
    /*  Define the analog filters (for design purpose) */
    b_butter_cg(6.2831853071795862 * wnom, b1, a1);
    b1_size[0] = 1;
    b1_size[1] = 4;
    b1_data[0] = b1[0];
    b1_data[1] = b1[1];
    b1_data[2] = b1[2];
    b1_data[3] = b1[3];

    /*  3rd order */
    butter_cg(6.2831853071795862 * (wnom * 3.125), b2, a2_data, a2_size);
    b2_size[0] = 1;
    b2_size[1] = 2;
    b2_data[0] = b2[0];
    b2_data[1] = b2[1];
    i0 = a2->size[0] * a2->size[1];
    a2->size[0] = 1;
    a2->size[1] = 2;
    emxEnsureCapacity_creal_T(a2, i0);
    a2->data[0] = a2_data[0];
    a2->data[1] = a2_data[1];

    /*  1st order */
    /*  Define the digital filters with fixed coefficients */
    memcpy(&hb1_coeff[0], &dv0[0], 15U * sizeof(double));
    hb3_coeff_size[0] = 1;
    hb3_coeff_size[1] = 3;
    hb3_coeff_data[0] = 0.25;
    hb3_coeff_data[1] = 0.5;
    hb3_coeff_data[2] = 0.25;
    dec_int3_coeff_size[0] = 1;
    dec_int3_coeff_size[1] = 29;
    memcpy(&dec_int3_coeff_data[0], &c_y[0], 29U * sizeof(double));
  }

  /*  Configure staging of filters */
  if (HB3 == 2.0) {
    nm1d2 = 50;
    idx = 49;
  } else if (HB3 == 3.0) {
    nm1d2 = 49;
    idx = 51;
  } else {
    nm1d2 = 49;
    idx = 49;
  }

  /*  convert the enables into a string */
  sigmax = rt_roundd_snf(HB1);
  if (sigmax < 128.0) {
    if (sigmax >= -128.0) {
      i1 = (signed char)sigmax;
    } else {
      i1 = MIN_int8_T;
    }
  } else if (sigmax >= 128.0) {
    i1 = MAX_int8_T;
  } else {
    i1 = 0;
  }

  i0 = 48 + i1;
  if (i0 > 127) {
    i0 = 127;
  }

  enables[0] = (signed char)i0;
  sigmax = rt_roundd_snf(HB2);
  if (sigmax < 128.0) {
    if (sigmax >= -128.0) {
      i1 = (signed char)sigmax;
    } else {
      i1 = MIN_int8_T;
    }
  } else if (sigmax >= 128.0) {
    i1 = MAX_int8_T;
  } else {
    i1 = 0;
  }

  i0 = 48 + i1;
  if (i0 > 127) {
    i0 = 127;
  }

  enables[1] = (signed char)i0;
  enables[2] = (signed char)nm1d2;
  enables[3] = (signed char)idx;

  /*  Find out the best fit delay on passband */
  memset(&w[0], 0, sizeof(double) << 11);
  memset(&phi[0], 0, sizeof(double) << 11);
  w[0] = -Fpass;
  for (nm1d2 = 0; nm1d2 < 2047; nm1d2++) {
    w[nm1d2 + 1] = w[0] - 2.0 * w[0] * (2.0 + (double)nm1d2) / 2048.0;
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
    sigmax = combinedResponse[i0].re * dcv0[i0].re - combinedResponse[i0].im *
      dcv0[i0].im;
    apnd = combinedResponse[i0].re * dcv0[i0].im + combinedResponse[i0].im *
      dcv0[i0].re;
    b_combinedResponse[i0] = sigmax;
    combinedResponse[i0].re = sigmax;
    combinedResponse[i0].im = apnd;
  }

  b_power(b_combinedResponse, invariance);
  for (i0 = 0; i0 < 2048; i0++) {
    b_combinedResponse[i0] = combinedResponse[i0].im;
  }

  b_power(b_combinedResponse, b_phi);
  for (i0 = 0; i0 < 2048; i0++) {
    invariance[i0] += b_phi[i0];
  }

  phi[0] = rt_atan2d_snf(combinedResponse[0].im, combinedResponse[0].re);
  for (nm1d2 = 0; nm1d2 < 2047; nm1d2++) {
    sigma = rt_atan2d_snf(combinedResponse[nm1d2 + 1].im, combinedResponse[nm1d2
                          + 1].re) - phi[nm1d2];
    phi[nm1d2 + 1] = phi[nm1d2] + (sigma - 6.2831853071795862 * floor(sigma /
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

  emxInit_real_T(&fg, 2);

  /*  Design the FIR */
  clkFIR = Rdata * FIR;
  Gstop = ceil(16384.0 * Fstop / clkFIR);
  Gpass = fmin(floor(16384.0 * Fpass / clkFIR), Gstop - 1.0);
  i0 = fg->size[0] * fg->size[1];
  fg->size[0] = 1;
  loop_ub = (int)(Gpass + 1.0);
  fg->size[1] = loop_ub;
  emxEnsureCapacity_real_T(fg, i0);
  for (i0 = 0; i0 < loop_ub; i0++) {
    fg->data[i0] = 0.0;
  }

  emxInit_real_T(&omega, 2);
  i0 = omega->size[0] * omega->size[1];
  omega->size[0] = 1;
  omega->size[1] = loop_ub;
  emxEnsureCapacity_real_T(omega, i0);
  for (i0 = 0; i0 < loop_ub; i0++) {
    omega->data[i0] = 0.0;
  }

  /*  passband */
  for (nm1d2 = 0; nm1d2 < loop_ub; nm1d2++) {
    fg->data[nm1d2] = ((1.0 + (double)nm1d2) - 1.0) / 16384.0;
    omega->data[nm1d2] = fg->data[nm1d2] * clkFIR;
  }

  emxInit_creal_T(&rg1, 2);
  emxInit_creal_T(&c_combinedResponse, 2);

  /*  Generate responses then convolve */
  b_generateCascadedResponseRx(enables, omega, converter_rate, hb1_coeff,
    hb2_coeff, hb3_coeff_data, hb3_coeff_size, dec_int3_coeff_data,
    dec_int3_coeff_size, c_combinedResponse);

  /*  Determine overall response with analog filters inline */
  b_analogresp(RxTx, omega, converter_rate, b1_data, b1_size, a1, b2_data,
               b2_size, a2, rg1);
  i0 = c_combinedResponse->size[0] * c_combinedResponse->size[1];
  i2 = rg1->size[0] * rg1->size[1];
  rg1->size[0] = 1;
  rg1->size[1] = c_combinedResponse->size[1];
  emxEnsureCapacity_creal_T(rg1, i2);
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    sigmax = c_combinedResponse->data[i0].re;
    apnd = c_combinedResponse->data[i0].im;
    rg1_re = rg1->data[i0].re;
    rg1_im = rg1->data[i0].im;
    rg1->data[i0].re = sigmax * rg1_re - apnd * rg1_im;
    rg1->data[i0].im = sigmax * rg1_im + apnd * rg1_re;
  }

  emxInit_creal_T(&r0, 2);
  i0 = r0->size[0] * r0->size[1];
  r0->size[0] = 1;
  r0->size[1] = omega->size[1];
  emxEnsureCapacity_creal_T(r0, i0);
  loop_ub = omega->size[0] * omega->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    re = omega->data[i0] * -0.0;
    im = omega->data[i0] * -6.2831853071795862;
    r0->data[i0].re = sigma * re;
    r0->data[i0].im = sigma * im;
  }

  c_exp(r0);
  rdivide_helper(r0, rg1, c_combinedResponse);
  sigma = Gpass + 1.0;

  /*  Expand memory correctly */
  emxInit_real_T(&F4, 2);
  if (rtIsNaN(Gstop)) {
    i0 = F4->size[0] * F4->size[1];
    F4->size[1] = 1;
    emxEnsureCapacity_real_T(F4, i0);
  } else if (8192.0 < Gstop) {
    i0 = F4->size[0] * F4->size[1];
    F4->size[1] = 0;
    emxEnsureCapacity_real_T(F4, i0);
  } else if (rtIsInf(Gstop) && (Gstop == 8192.0)) {
    i0 = F4->size[0] * F4->size[1];
    F4->size[1] = 1;
    emxEnsureCapacity_real_T(F4, i0);
  } else if (Gstop == Gstop) {
    i0 = F4->size[0] * F4->size[1];
    F4->size[1] = (int)(8192.0 - Gstop) + 1;
    emxEnsureCapacity_real_T(F4, i0);
  } else {
    sigmax = floor((8192.0 - Gstop) + 0.5);
    apnd = Gstop + sigmax;
    if (fabs(apnd - 8192.0) < 4.4408920985006262E-16 * fmax(fabs(Gstop), 8192.0))
    {
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

    i0 = F4->size[0] * F4->size[1];
    F4->size[0] = 1;
    F4->size[1] = n;
    emxEnsureCapacity_real_T(F4, i0);
    if ((n > 0) && (n > 1)) {
      F4->data[n - 1] = apnd;
      nm1d2 = (n - 1) / 2;
      for (k = 0; k <= nm1d2 - 2; k++) {
        F4->data[1 + k] = Gstop + (1.0 + (double)k);
        F4->data[(n - k) - 2] = apnd - (1.0 + (double)k);
      }

      if (nm1d2 << 1 == n - 1) {
        F4->data[nm1d2] = (Gstop + apnd) / 2.0;
      } else {
        F4->data[nm1d2] = Gstop + (double)nm1d2;
        F4->data[nm1d2 + 1] = apnd - (double)nm1d2;
      }
    }
  }

  emxInit_real_T(&fg2, 2);
  i0 = fg2->size[0] * fg2->size[1];
  fg2->size[0] = 1;
  fg2->size[1] = (int)((unsigned int)F4->size[1] + fg->size[1]);
  emxEnsureCapacity_real_T(fg2, i0);
  loop_ub = (int)((unsigned int)F4->size[1] + fg->size[1]);
  for (i0 = 0; i0 < loop_ub; i0++) {
    fg2->data[i0] = 0.0;
  }

  loop_ub = fg->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    fg2->data[i0] = fg->data[i0];
  }

  emxInit_real_T(&sw, 2);
  if (rtIsNaN(Gstop)) {
    i0 = sw->size[0] * sw->size[1];
    sw->size[1] = 1;
    emxEnsureCapacity_real_T(sw, i0);
  } else if (8192.0 < Gstop) {
    i0 = sw->size[0] * sw->size[1];
    sw->size[1] = 0;
    emxEnsureCapacity_real_T(sw, i0);
  } else if (rtIsInf(Gstop) && (Gstop == 8192.0)) {
    i0 = sw->size[0] * sw->size[1];
    sw->size[1] = 1;
    emxEnsureCapacity_real_T(sw, i0);
  } else if (Gstop == Gstop) {
    i0 = sw->size[0] * sw->size[1];
    sw->size[1] = (int)(8192.0 - Gstop) + 1;
    emxEnsureCapacity_real_T(sw, i0);
  } else {
    sigmax = floor((8192.0 - Gstop) + 0.5);
    apnd = Gstop + sigmax;
    if (fabs(apnd - 8192.0) < 4.4408920985006262E-16 * fmax(fabs(Gstop), 8192.0))
    {
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
    emxEnsureCapacity_real_T(sw, i0);
    if (n > 0) {
      sw->data[0] = Gstop;
      if (n > 1) {
        sw->data[n - 1] = apnd;
        nm1d2 = (n - 1) / 2;
        for (k = 0; k <= nm1d2 - 2; k++) {
          sw->data[1 + k] = Gstop + (1.0 + (double)k);
          sw->data[(n - k) - 2] = apnd - (1.0 + (double)k);
        }

        if (nm1d2 << 1 == n - 1) {
          sw->data[nm1d2] = (Gstop + apnd) / 2.0;
        } else {
          sw->data[nm1d2] = Gstop + (double)nm1d2;
          sw->data[nm1d2 + 1] = apnd - (double)nm1d2;
        }
      }
    }
  }

  emxInit_real_T(&omega2, 2);
  i0 = omega2->size[0] * omega2->size[1];
  omega2->size[0] = 1;
  omega2->size[1] = (int)((unsigned int)sw->size[1] + omega->size[1]);
  emxEnsureCapacity_real_T(omega2, i0);
  loop_ub = (int)((unsigned int)sw->size[1] + omega->size[1]);
  for (i0 = 0; i0 < loop_ub; i0++) {
    omega2->data[i0] = 0.0;
  }

  loop_ub = omega->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    omega2->data[i0] = omega->data[i0];
  }

  emxInit_real_T(&F3, 2);
  if (rtIsNaN(Gstop)) {
    i0 = F3->size[0] * F3->size[1];
    F3->size[1] = 1;
    emxEnsureCapacity_real_T(F3, i0);
  } else if (8192.0 < Gstop) {
    i0 = F3->size[0] * F3->size[1];
    F3->size[1] = 0;
    emxEnsureCapacity_real_T(F3, i0);
  } else if (rtIsInf(Gstop) && (Gstop == 8192.0)) {
    i0 = F3->size[0] * F3->size[1];
    F3->size[1] = 1;
    emxEnsureCapacity_real_T(F3, i0);
  } else if (Gstop == Gstop) {
    i0 = F3->size[0] * F3->size[1];
    F3->size[1] = (int)(8192.0 - Gstop) + 1;
    emxEnsureCapacity_real_T(F3, i0);
  } else {
    sigmax = floor((8192.0 - Gstop) + 0.5);
    apnd = Gstop + sigmax;
    if (fabs(apnd - 8192.0) < 4.4408920985006262E-16 * fmax(fabs(Gstop), 8192.0))
    {
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

    i0 = F3->size[0] * F3->size[1];
    F3->size[0] = 1;
    F3->size[1] = n;
    emxEnsureCapacity_real_T(F3, i0);
    if (n > 0) {
      F3->data[0] = Gstop;
      if (n > 1) {
        F3->data[n - 1] = apnd;
        nm1d2 = (n - 1) / 2;
        for (k = 0; k <= nm1d2 - 2; k++) {
          F3->data[1 + k] = Gstop + (1.0 + (double)k);
          F3->data[(n - k) - 2] = apnd - (1.0 + (double)k);
        }

        if (nm1d2 << 1 == n - 1) {
          F3->data[nm1d2] = (Gstop + apnd) / 2.0;
        } else {
          F3->data[nm1d2] = Gstop + (double)nm1d2;
          F3->data[nm1d2 + 1] = apnd - (double)nm1d2;
        }
      }
    }
  }

  emxInit_creal_T(&rgN, 2);
  i0 = rgN->size[0] * rgN->size[1];
  rgN->size[0] = 1;
  rgN->size[1] = (int)((unsigned int)F3->size[1] + c_combinedResponse->size[1]);
  emxEnsureCapacity_creal_T(rgN, i0);
  loop_ub = (int)((unsigned int)F3->size[1] + c_combinedResponse->size[1]);
  for (i0 = 0; i0 < loop_ub; i0++) {
    rgN->data[i0].re = 0.0;
    rgN->data[i0].im = 0.0;
  }

  loop_ub = c_combinedResponse->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    rgN->data[i0] = c_combinedResponse->data[i0];
  }

  /*  stop band */
  i0 = (int)(8192.0 + (1.0 - Gstop));
  for (idx = 0; idx < i0; idx++) {
    sigma++;
    i2 = (int)sigma - 1;
    fg2->data[i2] = (Gstop + (double)idx) / 16384.0;
    omega2->data[i2] = fg2->data[i2] * clkFIR;
    rgN->data[i2].re = 0.0;
    rgN->data[i2].im = 0.0;
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
  idx = b_omega2->size[0] * b_omega2->size[1];
  b_omega2->size[0] = 1;
  loop_ub = i2 - i0;
  b_omega2->size[1] = loop_ub;
  emxEnsureCapacity_real_T(b_omega2, idx);
  for (i2 = 0; i2 < loop_ub; i2++) {
    b_omega2->data[i2] = omega2->data[i0 + i2];
  }

  b_generateCascadedResponseRx(enables, b_omega2, converter_rate, hb1_coeff,
    hb2_coeff, hb3_coeff_data, hb3_coeff_size, dec_int3_coeff_data,
    dec_int3_coeff_size, c_combinedResponse);
  if (Gpass + 2.0 > omega2->size[1]) {
    i0 = 0;
    i2 = 0;
  } else {
    i0 = (int)(Gpass + 2.0) - 1;
    i2 = omega2->size[1];
  }

  idx = b_omega2->size[0] * b_omega2->size[1];
  b_omega2->size[0] = 1;
  loop_ub = i2 - i0;
  b_omega2->size[1] = loop_ub;
  emxEnsureCapacity_real_T(b_omega2, idx);
  for (i2 = 0; i2 < loop_ub; i2++) {
    b_omega2->data[i2] = omega2->data[i0 + i2];
  }

  emxInit_creal_T(&d_combinedResponse, 2);
  b_analogresp(RxTx, b_omega2, converter_rate, b1_data, b1_size, a1, b2_data,
               b2_size, a2, r0);
  i0 = d_combinedResponse->size[0] * d_combinedResponse->size[1];
  d_combinedResponse->size[0] = 1;
  d_combinedResponse->size[1] = c_combinedResponse->size[1];
  emxEnsureCapacity_creal_T(d_combinedResponse, i0);
  loop_ub = c_combinedResponse->size[0] * c_combinedResponse->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    sigmax = c_combinedResponse->data[i0].re;
    apnd = c_combinedResponse->data[i0].im;
    re = r0->data[i0].re;
    im = r0->data[i0].im;
    d_combinedResponse->data[i0].re = sigmax * re - apnd * im;
    d_combinedResponse->data[i0].im = sigmax * im + apnd * re;
  }

  b_abs(d_combinedResponse, fg);
  if (b_strcmp(RxTx)) {
    i0 = fg->size[0] * fg->size[1];
    i2 = fg->size[0] * fg->size[1];
    fg->size[0] = 1;
    emxEnsureCapacity_real_T(fg, i2);
    sigmax = dBinv(-Astop);
    loop_ub = i0 - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      fg->data[i0] /= sigmax;
    }
  } else {
    sigma = FIR;
    b_sqrt(&sigma);
    i0 = fg->size[0] * fg->size[1];
    i2 = fg->size[0] * fg->size[1];
    fg->size[0] = 1;
    emxEnsureCapacity_real_T(fg, i2);
    sigmax = dBinv(-Astop);
    loop_ub = i0 - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      fg->data[i0] = sigma * fg->data[i0] / sigmax;
    }
  }

  sigma = dBinv(FIRdBmin);
  i0 = omega->size[0] * omega->size[1];
  omega->size[0] = 1;
  omega->size[1] = fg->size[1];
  emxEnsureCapacity_real_T(omega, i0);
  nm1d2 = fg->size[1];
  for (k = 0; k < nm1d2; k++) {
    omega->data[k] = fmax(fg->data[k], sigma);
  }

  if (phEQ == -1.0) {
    b_abs(rgN, F4);
    i0 = rgN->size[0] * rgN->size[1];
    rgN->size[0] = 1;
    rgN->size[1] = F4->size[1];
    emxEnsureCapacity_creal_T(rgN, i0);
    loop_ub = F4->size[0] * F4->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      rgN->data[i0].re = F4->data[i0];
      rgN->data[i0].im = 0.0;
    }
  }

  emxInit_real_T(&weight, 2);
  b_abs(rg1, F4);
  sigmax = dBinv(Apass / 2.0) - 1.0;
  i0 = weight->size[0] * weight->size[1];
  weight->size[0] = 1;
  weight->size[1] = F4->size[1] + omega->size[1];
  emxEnsureCapacity_real_T(weight, i0);
  loop_ub = F4->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    weight->data[i0] = F4->data[i0] / sigmax;
  }

  loop_ub = omega->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    weight->data[i0 + F4->size[1]] = omega->data[i0];
  }

  n = weight->size[1];
  if (weight->size[1] <= 2) {
    if (weight->size[1] == 1) {
      apnd = weight->data[0];
    } else if ((weight->data[0] < weight->data[1]) || (rtIsNaN(weight->data[0]) &&
                (!rtIsNaN(weight->data[1])))) {
      apnd = weight->data[1];
    } else {
      apnd = weight->data[0];
    }
  } else {
    if (!rtIsNaN(weight->data[0])) {
      idx = 1;
    } else {
      idx = 0;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= weight->size[1])) {
        if (!rtIsNaN(weight->data[k - 1])) {
          idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }

    if (idx == 0) {
      apnd = weight->data[0];
    } else {
      apnd = weight->data[idx - 1];
      i0 = idx + 1;
      for (k = i0; k <= n; k++) {
        if (apnd < weight->data[k - 1]) {
          apnd = weight->data[k - 1];
        }
      }
    }
  }

  i0 = weight->size[0] * weight->size[1];
  i2 = weight->size[0] * weight->size[1];
  weight->size[0] = 1;
  emxEnsureCapacity_real_T(weight, i2);
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    weight->data[i0] /= apnd;
  }

  /*  Set up design for FIR filter */
  i0 = fg->size[0] * fg->size[1];
  fg->size[0] = 1;
  fg->size[1] = rgN->size[1];
  emxEnsureCapacity_real_T(fg, i0);
  loop_ub = rgN->size[0] * rgN->size[1];
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
  emxEnsureCapacity_real_T(F1, i0);
  for (i0 = 0; i0 < loop_ub; i0++) {
    F1->data[i0] = fg2->data[i0] * 2.0;
  }

  if (Gpass + 2.0 > fg2->size[1]) {
    i0 = 0;
    i2 = 0;
  } else {
    i0 = (int)(Gpass + 2.0) - 1;
    i2 = fg2->size[1];
  }

  emxInit_real_T(&F2, 2);
  idx = F2->size[0] * F2->size[1];
  F2->size[0] = 1;
  loop_ub = i2 - i0;
  F2->size[1] = loop_ub;
  emxEnsureCapacity_real_T(F2, idx);
  for (i2 = 0; i2 < loop_ub; i2++) {
    F2->data[i2] = fg2->data[i0 + i2] * 2.0;
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
  emxEnsureCapacity_real_T(A1, i0);
  for (i0 = 0; i0 < loop_ub; i0++) {
    A1->data[i0] = fg->data[i0];
  }

  if (Gpass + 2.0 > fg->size[1]) {
    i0 = 0;
    i2 = 0;
  } else {
    i0 = (int)(Gpass + 2.0) - 1;
    i2 = fg->size[1];
  }

  emxInit_real_T(&A2, 2);
  idx = A2->size[0] * A2->size[1];
  A2->size[0] = 1;
  loop_ub = i2 - i0;
  A2->size[1] = loop_ub;
  emxEnsureCapacity_real_T(A2, idx);
  for (i2 = 0; i2 < loop_ub; i2++) {
    A2->data[i2] = fg->data[i0 + i2];
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
  emxEnsureCapacity_real_T(W1, i0);
  for (i0 = 0; i0 < loop_ub; i0++) {
    W1->data[i0] = weight->data[i0];
  }

  if (Gpass + 2.0 > weight->size[1]) {
    i0 = 0;
    i2 = 0;
  } else {
    i0 = (int)(Gpass + 2.0) - 1;
    i2 = weight->size[1];
  }

  emxInit_real_T(&W2, 2);
  idx = W2->size[0] * W2->size[1];
  W2->size[0] = 1;
  loop_ub = i2 - i0;
  W2->size[1] = loop_ub;
  emxEnsureCapacity_real_T(W2, idx);
  for (i2 = 0; i2 < loop_ub; i2++) {
    W2->data[i2] = weight->data[i0 + i2];
  }

  emxInit_real_T(&tap_store, 2);

  /*  Determine the number of taps for FIR */
  /*  if strcmp(input.RxTx, 'Rx') */
  /*      if hb3 == 1 */
  /*          N = min(16*floor(input.converter_rate/(input.Rdata)),128); */
  /*      else */
  /*          N = min(16*floor(input.converter_rate/(2*input.Rdata)),128); */
  /*      end */
  /*  else */
  /*      switch input.FIR */
  /*          case 1 */
  /*              Nmax = 64; */
  /*          case 2 */
  /*              Nmax = 128; */
  /*          case 4 */
  /*              Nmax = 128; */
  /*          otherwise */
  /*              error('Wrong FIR Type'); */
  /*      end */
  /*      N = min(16*floor(input.converter_rate*input.DAC_div/(2*input.Rdata)),Nmax); */
  /*  end */
  clkFIR = maxTaps;
  i0 = tap_store->size[0] * tap_store->size[1];
  loop_ub = (int)(maxTaps / 16.0);
  tap_store->size[0] = loop_ub;
  i2 = (int)maxTaps;
  tap_store->size[1] = i2;
  emxEnsureCapacity_real_T(tap_store, i0);
  n = loop_ub * i2;
  for (i0 = 0; i0 < n; i0++) {
    tap_store->data[i0] = 0.0;
  }

  emxInit_real_T(&Apass_actual_vector, 1);
  i0 = Apass_actual_vector->size[0];
  Apass_actual_vector->size[0] = loop_ub;
  emxEnsureCapacity_real_T(Apass_actual_vector, i0);
  for (i0 = 0; i0 < loop_ub; i0++) {
    Apass_actual_vector->data[i0] = 0.0;
  }

  emxInit_real_T(&Astop_actual_vector, 1);
  i0 = Astop_actual_vector->size[0];
  Astop_actual_vector->size[0] = loop_ub;
  emxEnsureCapacity_real_T(Astop_actual_vector, i0);
  for (i0 = 0; i0 < loop_ub; i0++) {
    Astop_actual_vector->data[i0] = 0.0;
  }

  i = 1U;

  /*  Design filter */
  emxInit_real_T(&ccoef, 2);
  emxInit_real_T(&b_ccoef, 2);
  emxInit_real_T(&b_F1, 2);
  emxInit_real_T(&b_W1, 2);
  do {
    exitg2 = 0;
    if (int_FIR != 0.0) {
      b1[0] = F1->data[0];
      b1[1] = F1->data[F1->size[1] - 1];
      b1[2] = F2->data[0];
      b1[3] = F2->data[F2->size[1] - 1];
      i0 = b_omega2->size[0] * b_omega2->size[1];
      b_omega2->size[0] = 1;
      b_omega2->size[1] = A1->size[1] + A2->size[1];
      emxEnsureCapacity_real_T(b_omega2, i0);
      loop_ub = A1->size[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        b_omega2->data[i0] = A1->data[i0];
      }

      loop_ub = A2->size[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        b_omega2->data[i0 + A1->size[1]] = A2->data[i0];
      }

      i0 = b_F1->size[0] * b_F1->size[1];
      b_F1->size[0] = 1;
      b_F1->size[1] = F1->size[1] + F2->size[1];
      emxEnsureCapacity_real_T(b_F1, i0);
      loop_ub = F1->size[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        b_F1->data[i0] = F1->data[i0];
      }

      loop_ub = F2->size[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        b_F1->data[i0 + F1->size[1]] = F2->data[i0];
      }

      i0 = b_W1->size[0] * b_W1->size[1];
      b_W1->size[0] = 1;
      b_W1->size[1] = W1->size[1] + W2->size[1];
      emxEnsureCapacity_real_T(b_W1, i0);
      loop_ub = W1->size[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        b_W1->data[i0] = W1->data[i0];
      }

      loop_ub = W2->size[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        b_W1->data[i0 + W1->size[1]] = W2->data[i0];
      }

      firpm_cg(clkFIR - 1.0, b1, b_omega2, b_F1, b_W1, ccoef);
    } else {
      /*  Check different designs until we reach required ripple condition */
      sigma = db2mag(-Astop);

      /*  Peak Ripple */
      i0 = ccoef->size[0] * ccoef->size[1];
      ccoef->size[0] = 1;
      ccoef->size[1] = 1;
      emxEnsureCapacity_real_T(ccoef, i0);
      ccoef->data[0] = 0.0;

      /*  Predef type */
      k = 0;
      exitg1 = false;
      while ((!exitg1) && (k < 126)) {
        b1[0] = F1->data[0];
        b1[1] = F1->data[F1->size[1] - 1];
        b1[2] = F2->data[0];
        b1[3] = F2->data[F2->size[1] - 1];
        i0 = b_omega2->size[0] * b_omega2->size[1];
        b_omega2->size[0] = 1;
        b_omega2->size[1] = A1->size[1] + A2->size[1];
        emxEnsureCapacity_real_T(b_omega2, i0);
        loop_ub = A1->size[1];
        for (i0 = 0; i0 < loop_ub; i0++) {
          b_omega2->data[i0] = A1->data[i0];
        }

        loop_ub = A2->size[1];
        for (i0 = 0; i0 < loop_ub; i0++) {
          b_omega2->data[i0 + A1->size[1]] = A2->data[i0];
        }

        i0 = b_F1->size[0] * b_F1->size[1];
        b_F1->size[0] = 1;
        b_F1->size[1] = F1->size[1] + F2->size[1];
        emxEnsureCapacity_real_T(b_F1, i0);
        loop_ub = F1->size[1];
        for (i0 = 0; i0 < loop_ub; i0++) {
          b_F1->data[i0] = F1->data[i0];
        }

        loop_ub = F2->size[1];
        for (i0 = 0; i0 < loop_ub; i0++) {
          b_F1->data[i0 + F1->size[1]] = F2->data[i0];
        }

        i0 = b_W1->size[0] * b_W1->size[1];
        b_W1->size[0] = 1;
        b_W1->size[1] = W1->size[1] + W2->size[1];
        emxEnsureCapacity_real_T(b_W1, i0);
        loop_ub = W1->size[1];
        for (i0 = 0; i0 < loop_ub; i0++) {
          b_W1->data[i0] = W1->data[i0];
        }

        loop_ub = W2->size[1];
        for (i0 = 0; i0 < loop_ub; i0++) {
          b_W1->data[i0 + W1->size[1]] = W2->data[i0];
        }

        b_firpm_cg(3.0 + (double)k, b1, b_omega2, b_F1, b_W1, b_ccoef, &valid,
                   &sigmax);
        i0 = ccoef->size[0] * ccoef->size[1];
        ccoef->size[0] = 1;
        ccoef->size[1] = b_ccoef->size[1];
        emxEnsureCapacity_real_T(ccoef, i0);
        loop_ub = b_ccoef->size[0] * b_ccoef->size[1];
        for (i0 = 0; i0 < loop_ub; i0++) {
          ccoef->data[i0] = b_ccoef->data[i0];
        }

        /*  Check if design meets specs */
        if ((sigmax < sigma) && valid) {
          exitg1 = true;
        } else {
          k++;
        }
      }
    }

    /*  Enable phase equalization and apply update to taps */
    if (phEQ != -1.0) {
      if (1 > fg2->size[1]) {
        i0 = 1;
        i2 = 1;
        idx = 0;
      } else {
        i0 = fg2->size[1];
        i2 = -1;
        idx = 1;
      }

      nm1d2 = fg->size[0] * fg->size[1];
      fg->size[0] = 1;
      loop_ub = div_s32_floor(idx - i0, i2);
      fg->size[1] = loop_ub + 1;
      emxEnsureCapacity_real_T(fg, nm1d2);
      for (idx = 0; idx <= loop_ub; idx++) {
        fg->data[idx] = 0.5 - fg2->data[(i0 + i2 * idx) - 1];
      }

      if (1 > rgN->size[1]) {
        i0 = 1;
        i2 = 1;
        idx = 0;
      } else {
        i0 = rgN->size[1];
        i2 = -1;
        idx = 1;
      }

      nm1d2 = omega->size[0] * omega->size[1];
      omega->size[0] = 1;
      loop_ub = div_s32_floor(idx - i0, i2);
      omega->size[1] = loop_ub + 1;
      emxEnsureCapacity_real_T(omega, nm1d2);
      for (idx = 0; idx <= loop_ub; idx++) {
        omega->data[idx] = rgN->data[(i0 + i2 * idx) - 1].im;
      }

      if (1 > weight->size[1]) {
        i0 = 1;
        i2 = 1;
        idx = 0;
      } else {
        i0 = weight->size[1];
        i2 = -1;
        idx = 1;
      }

      nm1d2 = sw->size[0] * sw->size[1];
      sw->size[0] = 1;
      loop_ub = div_s32_floor(idx - i0, i2);
      idx = loop_ub + 1;
      sw->size[1] = idx;
      emxEnsureCapacity_real_T(sw, nm1d2);
      for (nm1d2 = 0; nm1d2 <= loop_ub; nm1d2++) {
        sw->data[nm1d2] = weight->data[(i0 + i2 * nm1d2) - 1];
      }

      if (1.0 > (8192.0 - Gstop) + 1.0) {
        n = 0;
      } else {
        n = (int)((8192.0 - Gstop) + 1.0);
      }

      i0 = F3->size[0] * F3->size[1];
      F3->size[0] = 1;
      F3->size[1] = n;
      emxEnsureCapacity_real_T(F3, i0);
      for (i0 = 0; i0 < n; i0++) {
        F3->data[i0] = fg->data[i0] * 2.0;
      }

      if ((8192.0 - Gstop) + 2.0 > fg->size[1]) {
        i0 = 0;
        i2 = 0;
      } else {
        i0 = (int)((8192.0 - Gstop) + 2.0) - 1;
        i2 = fg->size[1];
      }

      nm1d2 = F4->size[0] * F4->size[1];
      F4->size[0] = 1;
      n = i2 - i0;
      F4->size[1] = n;
      emxEnsureCapacity_real_T(F4, nm1d2);
      for (i2 = 0; i2 < n; i2++) {
        F4->data[i2] = fg->data[i0 + i2] * 2.0;
      }

      if (1.0 > (8192.0 - Gstop) + 1.0) {
        n = 0;
      } else {
        n = (int)((8192.0 - Gstop) + 1.0);
      }

      if ((8192.0 - Gstop) + 2.0 > omega->size[1]) {
        i0 = 0;
        i2 = 0;
      } else {
        i0 = (int)((8192.0 - Gstop) + 2.0) - 1;
        i2 = omega->size[1];
      }

      if (1.0 > (8192.0 - Gstop) + 1.0) {
        k = 0;
      } else {
        k = (int)((8192.0 - Gstop) + 1.0);
      }

      if ((8192.0 - Gstop) + 2.0 > idx) {
        idx = 0;
        loop_ub = -1;
      } else {
        idx = (int)((8192.0 - Gstop) + 2.0) - 1;
      }

      if (int_FIR != 0.0) {
        sigma = clkFIR - 1.0;
      } else {
        sigma = (double)ccoef->size[1] - 1.0;
      }

      b1[0] = F3->data[0];
      b1[1] = F3->data[F3->size[1] - 1];
      b1[2] = F4->data[0];
      b1[3] = F4->data[F4->size[1] - 1];
      nm1d2 = b_W1->size[0] * b_W1->size[1];
      b_W1->size[0] = 1;
      b_W1->size[1] = (n + i2) - i0;
      emxEnsureCapacity_real_T(b_W1, nm1d2);
      for (nm1d2 = 0; nm1d2 < n; nm1d2++) {
        b_W1->data[nm1d2] = omega->data[nm1d2];
      }

      nm1d2 = i2 - i0;
      for (i2 = 0; i2 < nm1d2; i2++) {
        b_W1->data[i2 + n] = omega->data[i0 + i2];
      }

      i0 = b_F1->size[0] * b_F1->size[1];
      b_F1->size[0] = 1;
      b_F1->size[1] = F3->size[1] + F4->size[1];
      emxEnsureCapacity_real_T(b_F1, i0);
      n = F3->size[1];
      for (i0 = 0; i0 < n; i0++) {
        b_F1->data[i0] = F3->data[i0];
      }

      n = F4->size[1];
      for (i0 = 0; i0 < n; i0++) {
        b_F1->data[i0 + F3->size[1]] = F4->data[i0];
      }

      i0 = b_omega2->size[0] * b_omega2->size[1];
      b_omega2->size[0] = 1;
      b_omega2->size[1] = ((k + loop_ub) - idx) + 1;
      emxEnsureCapacity_real_T(b_omega2, i0);
      for (i0 = 0; i0 < k; i0++) {
        b_omega2->data[i0] = sw->data[i0];
      }

      loop_ub -= idx;
      for (i0 = 0; i0 <= loop_ub; i0++) {
        b_omega2->data[i0 + k] = sw->data[idx + i0];
      }

      firpm_cg(sigma, b1, b_W1, b_F1, b_omega2, fg);
      i0 = fg->size[1];
      for (k = 0; k < i0; k++) {
        fg->data[k] = -fg->data[k] * mpower(-1.0, (1.0 + (double)k) - 1.0);
      }
    } else {
      b2[1] = ccoef->size[1];
      i0 = fg->size[0] * fg->size[1];
      fg->size[0] = 1;
      fg->size[1] = (int)b2[1];
      emxEnsureCapacity_real_T(fg, i0);
      loop_ub = (int)b2[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        fg->data[i0] = 0.0;
      }
    }

    loop_ub = ccoef->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      tap_store->data[((int)i + tap_store->size[0] * i0) - 1] = ccoef->data[i0]
        + fg->data[i0];
    }

    /*  scoef ==0 when no EQ */
    determineBestFractionLength(tap_store, i, ccoef->size[1], b_F1);
    loop_ub = b_F1->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      tap_store->data[((int)i + tap_store->size[0] * i0) - 1] = b_F1->data[i0];
    }

    if (b_strcmp(RxTx)) {
      if (1.0 > Gpass + 1.0) {
        loop_ub = 0;
      } else {
        loop_ub = (int)(Gpass + 1.0);
      }

      if (1 > ccoef->size[1]) {
        n = 0;
      } else {
        n = ccoef->size[1];
      }

      i0 = b_omega2->size[0] * b_omega2->size[1];
      b_omega2->size[0] = 1;
      b_omega2->size[1] = loop_ub;
      emxEnsureCapacity_real_T(b_omega2, i0);
      for (i0 = 0; i0 < loop_ub; i0++) {
        b_omega2->data[i0] = omega2->data[i0];
      }

      i0 = b_F1->size[0] * b_F1->size[1];
      b_F1->size[0] = 1;
      b_F1->size[1] = n;
      emxEnsureCapacity_real_T(b_F1, i0);
      for (i0 = 0; i0 < n; i0++) {
        b_F1->data[i0] = tap_store->data[((int)i + tap_store->size[0] * i0) - 1];
      }

      c_generateCascadedResponseRx(enables, b_omega2, converter_rate, hb1_coeff,
        hb2_coeff, hb3_coeff_data, hb3_coeff_size, dec_int3_coeff_data,
        dec_int3_coeff_size, b_F1, c_combinedResponse);
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

      idx = b_omega2->size[0] * b_omega2->size[1];
      b_omega2->size[0] = 1;
      n = i2 - i0;
      b_omega2->size[1] = n;
      emxEnsureCapacity_real_T(b_omega2, idx);
      for (i2 = 0; i2 < n; i2++) {
        b_omega2->data[i2] = omega2->data[i0 + i2];
      }

      i0 = b_F1->size[0] * b_F1->size[1];
      b_F1->size[0] = 1;
      b_F1->size[1] = loop_ub;
      emxEnsureCapacity_real_T(b_F1, i0);
      for (i0 = 0; i0 < loop_ub; i0++) {
        b_F1->data[i0] = tap_store->data[((int)i + tap_store->size[0] * i0) - 1];
      }

      c_generateCascadedResponseRx(enables, b_omega2, converter_rate, hb1_coeff,
        hb2_coeff, hb3_coeff_data, hb3_coeff_size, dec_int3_coeff_data,
        dec_int3_coeff_size, b_F1, rg1);
      if (1.0 > Gpass + 1.0) {
        loop_ub = 0;
      } else {
        loop_ub = (int)(Gpass + 1.0);
      }

      i0 = b_omega2->size[0] * b_omega2->size[1];
      b_omega2->size[0] = 1;
      b_omega2->size[1] = loop_ub;
      emxEnsureCapacity_real_T(b_omega2, i0);
      for (i0 = 0; i0 < loop_ub; i0++) {
        b_omega2->data[i0] = omega2->data[i0];
      }

      c_analogresp(b_omega2, converter_rate, b1_data, b1_size, a1, b2_data,
                   b2_size, a2, r0);
      i0 = d_combinedResponse->size[0] * d_combinedResponse->size[1];
      d_combinedResponse->size[0] = 1;
      d_combinedResponse->size[1] = r0->size[1];
      emxEnsureCapacity_creal_T(d_combinedResponse, i0);
      loop_ub = r0->size[0] * r0->size[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        re = r0->data[i0].re;
        im = r0->data[i0].im;
        sigmax = c_combinedResponse->data[i0].re;
        apnd = c_combinedResponse->data[i0].im;
        d_combinedResponse->data[i0].re = re * sigmax - im * apnd;
        d_combinedResponse->data[i0].im = re * apnd + im * sigmax;
      }

      b_abs(d_combinedResponse, fg);
      if (Gpass + 2.0 > omega2->size[1]) {
        i0 = 0;
        i2 = 0;
      } else {
        i0 = (int)(Gpass + 2.0) - 1;
        i2 = omega2->size[1];
      }

      idx = b_omega2->size[0] * b_omega2->size[1];
      b_omega2->size[0] = 1;
      loop_ub = i2 - i0;
      b_omega2->size[1] = loop_ub;
      emxEnsureCapacity_real_T(b_omega2, idx);
      for (i2 = 0; i2 < loop_ub; i2++) {
        b_omega2->data[i2] = omega2->data[i0 + i2];
      }

      c_analogresp(b_omega2, converter_rate, b1_data, b1_size, a1, b2_data,
                   b2_size, a2, r0);
      i0 = d_combinedResponse->size[0] * d_combinedResponse->size[1];
      d_combinedResponse->size[0] = 1;
      d_combinedResponse->size[1] = r0->size[1];
      emxEnsureCapacity_creal_T(d_combinedResponse, i0);
      loop_ub = r0->size[0] * r0->size[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        re = r0->data[i0].re;
        im = r0->data[i0].im;
        rg1_re = rg1->data[i0].re;
        rg1_im = rg1->data[i0].im;
        d_combinedResponse->data[i0].re = re * rg1_re - im * rg1_im;
        d_combinedResponse->data[i0].im = re * rg1_im + im * rg1_re;
      }

      b_abs(d_combinedResponse, omega);
    } else {
      /*  TX */
      if (1.0 > Gpass + 1.0) {
        loop_ub = 0;
      } else {
        loop_ub = (int)(Gpass + 1.0);
      }

      if (1 > ccoef->size[1]) {
        n = 0;
      } else {
        n = ccoef->size[1];
      }

      i0 = b_omega2->size[0] * b_omega2->size[1];
      b_omega2->size[0] = 1;
      b_omega2->size[1] = loop_ub;
      emxEnsureCapacity_real_T(b_omega2, i0);
      for (i0 = 0; i0 < loop_ub; i0++) {
        b_omega2->data[i0] = omega2->data[i0];
      }

      i0 = b_F1->size[0] * b_F1->size[1];
      b_F1->size[0] = 1;
      b_F1->size[1] = n;
      emxEnsureCapacity_real_T(b_F1, i0);
      for (i0 = 0; i0 < n; i0++) {
        b_F1->data[i0] = tap_store->data[((int)i + tap_store->size[0] * i0) - 1];
      }

      c_generateCascadedResponseRx(enables, b_omega2, converter_rate, hb1_coeff,
        hb2_coeff, hb3_coeff_data, hb3_coeff_size, dec_int3_coeff_data,
        dec_int3_coeff_size, b_F1, c_combinedResponse);
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

      idx = b_omega2->size[0] * b_omega2->size[1];
      b_omega2->size[0] = 1;
      n = i2 - i0;
      b_omega2->size[1] = n;
      emxEnsureCapacity_real_T(b_omega2, idx);
      for (i2 = 0; i2 < n; i2++) {
        b_omega2->data[i2] = omega2->data[i0 + i2];
      }

      i0 = b_F1->size[0] * b_F1->size[1];
      b_F1->size[0] = 1;
      b_F1->size[1] = loop_ub;
      emxEnsureCapacity_real_T(b_F1, i0);
      for (i0 = 0; i0 < loop_ub; i0++) {
        b_F1->data[i0] = tap_store->data[((int)i + tap_store->size[0] * i0) - 1];
      }

      c_generateCascadedResponseRx(enables, b_omega2, converter_rate, hb1_coeff,
        hb2_coeff, hb3_coeff_data, hb3_coeff_size, dec_int3_coeff_data,
        dec_int3_coeff_size, b_F1, rg1);
      if (1.0 > Gpass + 1.0) {
        loop_ub = 0;
      } else {
        loop_ub = (int)(Gpass + 1.0);
      }

      i0 = b_omega2->size[0] * b_omega2->size[1];
      b_omega2->size[0] = 1;
      b_omega2->size[1] = loop_ub;
      emxEnsureCapacity_real_T(b_omega2, i0);
      for (i0 = 0; i0 < loop_ub; i0++) {
        b_omega2->data[i0] = omega2->data[i0];
      }

      d_analogresp(b_omega2, converter_rate, b1_data, b1_size, a1, b2_data,
                   b2_size, a2, r0);
      i0 = d_combinedResponse->size[0] * d_combinedResponse->size[1];
      d_combinedResponse->size[0] = 1;
      d_combinedResponse->size[1] = c_combinedResponse->size[1];
      emxEnsureCapacity_creal_T(d_combinedResponse, i0);
      loop_ub = c_combinedResponse->size[0] * c_combinedResponse->size[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        sigmax = c_combinedResponse->data[i0].re;
        apnd = c_combinedResponse->data[i0].im;
        re = r0->data[i0].re;
        im = r0->data[i0].im;
        d_combinedResponse->data[i0].re = sigmax * re - apnd * im;
        d_combinedResponse->data[i0].im = sigmax * im + apnd * re;
      }

      b_abs(d_combinedResponse, fg);
      if (Gpass + 2.0 > omega2->size[1]) {
        i0 = 0;
        i2 = 0;
      } else {
        i0 = (int)(Gpass + 2.0) - 1;
        i2 = omega2->size[1];
      }

      idx = b_omega2->size[0] * b_omega2->size[1];
      b_omega2->size[0] = 1;
      loop_ub = i2 - i0;
      b_omega2->size[1] = loop_ub;
      emxEnsureCapacity_real_T(b_omega2, idx);
      for (i2 = 0; i2 < loop_ub; i2++) {
        b_omega2->data[i2] = omega2->data[i0 + i2];
      }

      d_analogresp(b_omega2, converter_rate, b1_data, b1_size, a1, b2_data,
                   b2_size, a2, r0);
      i0 = d_combinedResponse->size[0] * d_combinedResponse->size[1];
      d_combinedResponse->size[0] = 1;
      d_combinedResponse->size[1] = rg1->size[1];
      emxEnsureCapacity_creal_T(d_combinedResponse, i0);
      loop_ub = rg1->size[0] * rg1->size[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        rg1_re = rg1->data[i0].re;
        rg1_im = rg1->data[i0].im;
        re = r0->data[i0].re;
        im = r0->data[i0].im;
        d_combinedResponse->data[i0].re = rg1_re * re - rg1_im * im;
        d_combinedResponse->data[i0].im = rg1_re * im + rg1_im * re;
      }

      b_abs(d_combinedResponse, omega);
    }

    /*  quantitative values about actual passband and stopband */
    n = fg->size[1];
    if (fg->size[1] <= 2) {
      if (fg->size[1] == 1) {
        apnd = fg->data[0];
      } else if ((fg->data[0] < fg->data[1]) || (rtIsNaN(fg->data[0]) &&
                  (!rtIsNaN(fg->data[1])))) {
        apnd = fg->data[1];
      } else {
        apnd = fg->data[0];
      }
    } else {
      if (!rtIsNaN(fg->data[0])) {
        idx = 1;
      } else {
        idx = 0;
        k = 2;
        exitg1 = false;
        while ((!exitg1) && (k <= fg->size[1])) {
          if (!rtIsNaN(fg->data[k - 1])) {
            idx = k;
            exitg1 = true;
          } else {
            k++;
          }
        }
      }

      if (idx == 0) {
        apnd = fg->data[0];
      } else {
        apnd = fg->data[idx - 1];
        i0 = idx + 1;
        for (k = i0; k <= n; k++) {
          if (apnd < fg->data[k - 1]) {
            apnd = fg->data[k - 1];
          }
        }
      }
    }

    n = fg->size[1];
    if (fg->size[1] <= 2) {
      if (fg->size[1] == 1) {
        sigma = fg->data[0];
      } else if ((fg->data[0] > fg->data[1]) || (rtIsNaN(fg->data[0]) &&
                  (!rtIsNaN(fg->data[1])))) {
        sigma = fg->data[1];
      } else {
        sigma = fg->data[0];
      }
    } else {
      if (!rtIsNaN(fg->data[0])) {
        idx = 1;
      } else {
        idx = 0;
        k = 2;
        exitg1 = false;
        while ((!exitg1) && (k <= fg->size[1])) {
          if (!rtIsNaN(fg->data[k - 1])) {
            idx = k;
            exitg1 = true;
          } else {
            k++;
          }
        }
      }

      if (idx == 0) {
        sigma = fg->data[0];
      } else {
        sigma = fg->data[idx - 1];
        i0 = idx + 1;
        for (k = i0; k <= n; k++) {
          if (sigma > fg->data[k - 1]) {
            sigma = fg->data[k - 1];
          }
        }
      }
    }

    i0 = (int)i - 1;
    Apass_actual_vector->data[i0] = mag2db(apnd) - mag2db(sigma);
    n = omega->size[1];
    if (omega->size[1] <= 2) {
      if (omega->size[1] == 1) {
        apnd = omega->data[0];
      } else if ((omega->data[0] < omega->data[1]) || (rtIsNaN(omega->data[0]) &&
                  (!rtIsNaN(omega->data[1])))) {
        apnd = omega->data[1];
      } else {
        apnd = omega->data[0];
      }
    } else {
      if (!rtIsNaN(omega->data[0])) {
        idx = 1;
      } else {
        idx = 0;
        k = 2;
        exitg1 = false;
        while ((!exitg1) && (k <= omega->size[1])) {
          if (!rtIsNaN(omega->data[k - 1])) {
            idx = k;
            exitg1 = true;
          } else {
            k++;
          }
        }
      }

      if (idx == 0) {
        apnd = omega->data[0];
      } else {
        apnd = omega->data[idx - 1];
        i2 = idx + 1;
        for (k = i2; k <= n; k++) {
          if (apnd < omega->data[k - 1]) {
            apnd = omega->data[k - 1];
          }
        }
      }
    }

    Astop_actual_vector->data[i0] = -mag2db(apnd);
    if (int_FIR == 0.0) {
      if (1 > ccoef->size[1]) {
        loop_ub = 0;
      } else {
        loop_ub = ccoef->size[1];
      }

      i0 = fg->size[0] * fg->size[1];
      fg->size[0] = 1;
      fg->size[1] = loop_ub;
      emxEnsureCapacity_real_T(fg, i0);
      for (i0 = 0; i0 < loop_ub; i0++) {
        fg->data[i0] = tap_store->data[tap_store->size[0] * i0];
      }

      *Apass_actual = Apass_actual_vector->data[0];
      *Astop_actual = Astop_actual_vector->data[0];
      exitg2 = 1;
    } else if ((Apass_actual_vector->data[0] > Apass) ||
               (Astop_actual_vector->data[0] < Astop)) {
      if (1.0 > clkFIR) {
        loop_ub = 0;
      } else {
        loop_ub = (int)clkFIR;
      }

      i0 = fg->size[0] * fg->size[1];
      fg->size[0] = 1;
      fg->size[1] = loop_ub;
      emxEnsureCapacity_real_T(fg, i0);
      for (i0 = 0; i0 < loop_ub; i0++) {
        fg->data[i0] = tap_store->data[tap_store->size[0] * i0];
      }

      *Apass_actual = Apass_actual_vector->data[0];
      *Astop_actual = Astop_actual_vector->data[0];
      exitg2 = 1;
    } else if ((Apass_actual_vector->data[i0] > Apass) ||
               (Astop_actual_vector->data[i0] < Astop)) {
      if (1.0 > clkFIR + 16.0) {
        loop_ub = 0;
      } else {
        loop_ub = (int)(clkFIR + 16.0);
      }

      i0 = fg->size[0] * fg->size[1];
      fg->size[0] = 1;
      fg->size[1] = loop_ub;
      emxEnsureCapacity_real_T(fg, i0);
      for (i0 = 0; i0 < loop_ub; i0++) {
        fg->data[i0] = tap_store->data[((int)i + tap_store->size[0] * i0) - 2];
      }

      nm1d2 = (int)i - 2;
      *Apass_actual = Apass_actual_vector->data[nm1d2];
      *Astop_actual = Astop_actual_vector->data[nm1d2];
      exitg2 = 1;
    } else {
      clkFIR -= 16.0;
      i++;
    }
  } while (exitg2 == 0);

  emxFree_real_T(&b_omega2);
  emxFree_creal_T(&d_combinedResponse);
  emxFree_real_T(&b_W1);
  emxFree_real_T(&b_F1);
  emxFree_creal_T(&r0);
  emxFree_creal_T(&c_combinedResponse);
  emxFree_real_T(&b_ccoef);
  emxFree_real_T(&F4);
  emxFree_real_T(&F3);
  emxFree_real_T(&sw);
  emxFree_real_T(&ccoef);
  emxFree_real_T(&Astop_actual_vector);
  emxFree_real_T(&Apass_actual_vector);
  emxFree_real_T(&tap_store);
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
  emxFree_creal_T(&rg1);
  emxFree_creal_T(&a2);
  emxFree_creal_T(&a1);
  if (c_strcmp(RxTx)) {
    if ((int_FIR == 1.0) && (FIR == 2.0)) {
      if (rt_remd_snf(fg->size[1], 32.0) != 0.0) {
        i0 = omega->size[0] * omega->size[1];
        omega->size[0] = 1;
        omega->size[1] = 16 + fg->size[1];
        emxEnsureCapacity_real_T(omega, i0);
        for (i0 = 0; i0 < 8; i0++) {
          omega->data[i0] = 0.0;
        }

        loop_ub = fg->size[1];
        for (i0 = 0; i0 < loop_ub; i0++) {
          omega->data[i0 + 8] = fg->data[i0];
        }

        for (i0 = 0; i0 < 8; i0++) {
          omega->data[(i0 + fg->size[1]) + 8] = 0.0;
        }

        i0 = fg->size[0] * fg->size[1];
        fg->size[0] = 1;
        fg->size[1] = omega->size[1];
        emxEnsureCapacity_real_T(fg, i0);
        loop_ub = omega->size[0] * omega->size[1];
        for (i0 = 0; i0 < loop_ub; i0++) {
          fg->data[i0] = omega->data[i0];
        }
      }
    } else {
      if ((int_FIR == 1.0) && (FIR == 4.0) && (rt_remd_snf(fg->size[1], 64.0) !=
           0.0)) {
        sigma = (ceil((double)fg->size[1] / 64.0) * 64.0 - (double)fg->size[1]) /
          2.0;
        i0 = omega->size[0] * omega->size[1];
        omega->size[0] = 1;
        loop_ub = (int)sigma;
        omega->size[1] = (loop_ub + fg->size[1]) + loop_ub;
        emxEnsureCapacity_real_T(omega, i0);
        for (i0 = 0; i0 < loop_ub; i0++) {
          omega->data[i0] = 0.0;
        }

        n = fg->size[1];
        for (i0 = 0; i0 < n; i0++) {
          omega->data[i0 + loop_ub] = fg->data[i0];
        }

        for (i0 = 0; i0 < loop_ub; i0++) {
          omega->data[(i0 + loop_ub) + fg->size[1]] = 0.0;
        }

        i0 = fg->size[0] * fg->size[1];
        fg->size[0] = 1;
        fg->size[1] = omega->size[1];
        emxEnsureCapacity_real_T(fg, i0);
        loop_ub = omega->size[0] * omega->size[1];
        for (i0 = 0; i0 < loop_ub; i0++) {
          fg->data[i0] = omega->data[i0];
        }
      }
    }
  }

  emxFree_real_T(&omega);

  /*  There will always be 128 taps output */
  memset(&firTapsPreScale[0], 0, sizeof(double) << 7);
  loop_ub = fg->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    firTapsPreScale[i0] = fg->data[i0];
  }

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
  if (!rtIsNaN(firTapsPreScale[0])) {
    idx = 1;
  } else {
    idx = 0;
    k = 2;
    exitg1 = false;
    while ((!exitg1) && (k < 129)) {
      if (!rtIsNaN(firTapsPreScale[k - 1])) {
        idx = k;
        exitg1 = true;
      } else {
        k++;
      }
    }
  }

  if (idx == 0) {
    apnd = firTapsPreScale[0];
  } else {
    apnd = firTapsPreScale[idx - 1];
    i0 = idx + 1;
    for (k = i0; k < 129; k++) {
      sigmax = firTapsPreScale[k - 1];
      if (apnd < sigmax) {
        apnd = sigmax;
      }
    }
  }

  sigma = b_log2(apnd);
  sigma = ceil(sigma);
  switch ((int)(1.0 + sigma)) {
   case 2:
    nm1d2 = 6;
    break;

   case 1:
    nm1d2 = 0;
    break;

   case 0:
    nm1d2 = -6;
    break;

   default:
    nm1d2 = -12;
    break;
  }

  if (b_strcmp(RxTx)) {
    if (1.0 + sigma > 2.0) {
      nm1d2 = 6;
    }
  } else {
    if (FIR == 2.0) {
      nm1d2 += 6;
    } else {
      if (FIR == 4.0) {
        nm1d2 += 12;
      }
    }

    if (nm1d2 > 0) {
      nm1d2 = 0;
    } else {
      if (nm1d2 < -6) {
        nm1d2 = -6;
      }
    }
  }

  /*  Scale taps */
  memcpy(&b_firTapsPreScale[0], &firTapsPreScale[0], sizeof(double) << 7);
  b_determineBestFractionLength(b_firTapsPreScale, firTapsPreScale);
  sigmax = mpower(2.0, 16.0 - (1.0 + sigma));
  for (i0 = 0; i0 < 128; i0++) {
    sigma = rt_roundd_snf(firTapsPreScale[i0] * sigmax);
    if (sigma < 32768.0) {
      if (sigma >= -32768.0) {
        i3 = (short)sigma;
      } else {
        i3 = MIN_int16_T;
      }
    } else if (sigma >= 32768.0) {
      i3 = MAX_int16_T;
    } else {
      i3 = 0;
    }

    outputTaps[i0] = i3;
  }

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
  *numOutputTaps = fg->size[1];
  *filterGain = nm1d2;
  emxFree_real_T(&fg);
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
