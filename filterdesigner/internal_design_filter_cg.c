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
static boolean_T anyNonFinite(const emxArray_creal_T *x);
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
static void b_polyval(const double p[7], const creal_T x[2048],
                      creal_T y[2048]);
static void b_power(const double a[2048], double y[2048]);
static void b_rdivide(const emxArray_creal_T *x, const emxArray_creal_T *y,
                      emxArray_creal_T *z);
static void b_sinc(emxArray_real_T *x);
static void b_sqrt(double *x);
static boolean_T b_strcmp(const char a[2]);
static double b_sum(const emxArray_real_T *x);
static void b_us(const double o[7], double u[7]);
static void b_xscal(int n, const creal_T a, emxArray_creal_T *x, int ix0, int
                    incx);
static void b_xzlartg(const creal_T f, const creal_T g, double *cs,
                      creal_T *sn);
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
        emxArray_real_T *extraTaps, emxArray_creal_T *combinedResponse);
static void c_polyval(const double p_data[], const int p_size[2], const creal_T
                      x[2048], creal_T y[2048]);
static void c_power(const emxArray_real_T *a, emxArray_real_T *y);
static void c_rdivide(const emxArray_real_T *x, const emxArray_real_T *y,
                      emxArray_real_T *z);
static void c_sqrt(creal_T *x);
static boolean_T c_strcmp(const char a[2]);
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
static void d_polyval(const double p[29], const creal_T x[2048],
                      creal_T y[2048]);
static void d_us(const double o[15], double u[29]);
static double db2mag(double ydb);
static void determineBestFractionLength(const emxArray_real_T *tap_store, double
                                        i, double M, emxArray_real_T *taps);
static int div_s32_floor(int numerator, int denominator);
static void e_firfreqz(const double b[29], const struct_T *options, creal_T h
                       [2048], double w[2048]);
static void e_freqz_cg(const double b[29], const double w[2048], double Fs,
                       creal_T hh[2048]);
static void e_polyval(const double p[13], const creal_T x[2048],
                      creal_T y[2048]);
static void e_us(const double o[7], double u[13]);
static void eig(const emxArray_creal_T *A, emxArray_creal_T *V);
static int eml_zlahqr(emxArray_creal_T *h);
static void f_firfreqz(const double b[13], const struct_T *options, creal_T h
                       [2048], double w[2048]);
static void f_freqz_cg(const double b[13], const double w[2048], double Fs,
                       creal_T hh[2048]);
static void f_polyval(const double p[57], const creal_T x[2048],
                      creal_T y[2048]);
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
static void g_polyval(const double p[43], const creal_T x[2048],
                      creal_T y[2048]);
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
static void h_polyval(const double p[19], const creal_T x[2048],
                      creal_T y[2048]);
static void h_us(const double o[7], double u[19]);
static void i_firfreqz(const double b[19], const struct_T *options, creal_T h
                       [2048], double w[2048]);
static void i_freqz_cg(const double b[19], const double w[2048], double Fs,
                       creal_T hh[2048]);
static void i_polyval(const double p[85], const creal_T x[2048],
                      creal_T y[2048]);
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
static void k_freqz_cg(const emxArray_real_T *w, double Fs,
                       emxArray_creal_T *hh);
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
static creal_T xzlarfg(creal_T *alpha1, creal_T *x);
static void xzlartg(const creal_T f, const creal_T g, double *cs, creal_T *sn,
                    creal_T *r);
static void zp2ss_cg(emxArray_creal_T *a, emxArray_real_T *b,
                     emxArray_real_T *c,
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
    boolean_T b_bool;
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
        if (kstr + 1 < 3) {
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
            if (kstr + 1 < 3) {
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
 * Return Type  : boolean_T
 */
static boolean_T anyNonFinite(const emxArray_creal_T *x)
{
    boolean_T p;
    int nx;
    int k;
    nx = x->size[0] * x->size[1];
    p = true;
    for (k = 0; k + 1 <= nx; k++) {
        if (p && ((!(rtIsInf(x->data[k].re) || rtIsInf(x->data[k].im))) &&
                  (!(rtIsNaN(x->data[k].re) || rtIsNaN(x->data[k].im))))) {
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
    int k;
    k = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = x->size[1];
    emxEnsureCapacity_real_T(y, k);
    for (k = 0; k + 1 <= x->size[1]; k++) {
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
    boolean_T b_bool;
    int kstr;
    int exitg1;
    static const char cv46[2] = { 'T', 'x' };

    emxArray_creal_T *r17;
    emxArray_real_T *r18;
    emxArray_real_T *r19;
    static const char cv47[2] = { 'R', 'x' };

    int loop_ub;
    double abc_re;
    double abc_im;
    double re;
    double im;
    b_bool = false;
    kstr = 0;
    do {
        exitg1 = 0;
        if (kstr + 1 < 3) {
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
            if (kstr + 1 < 3) {
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

    emxInit_creal_T(&r17, 2);
    emxInit_real_T(&r18, 2);
    emxInit_real_T(&r19, 2);
    switch (kstr) {
    case 0:
        rdivide(f, Fconverter, r19);
        b_sinc(r19);
        kstr = r18->size[0] * r18->size[1];
        r18->size[0] = 1;
        r18->size[1] = f->size[1];
        emxEnsureCapacity_real_T(r18, kstr);
        loop_ub = f->size[0] * f->size[1];
        for (kstr = 0; kstr < loop_ub; kstr++) {
            r18->data[kstr] = 6.2831853071795862 * f->data[kstr];
        }

        b_freqs_cg(b1_data, b1_size, a1, r18, abc);
        kstr = r18->size[0] * r18->size[1];
        r18->size[0] = 1;
        r18->size[1] = f->size[1];
        emxEnsureCapacity_real_T(r18, kstr);
        loop_ub = f->size[0] * f->size[1];
        for (kstr = 0; kstr < loop_ub; kstr++) {
            r18->data[kstr] = 6.2831853071795862 * f->data[kstr];
        }

        b_freqs_cg(b2_data, b2_size, a2, r18, r17);
        kstr = abc->size[0] * abc->size[1];
        abc->size[0] = 1;
        abc->size[1] = r19->size[1];
        emxEnsureCapacity_creal_T(abc, kstr);
        loop_ub = r19->size[0] * r19->size[1];
        for (kstr = 0; kstr < loop_ub; kstr++) {
            abc_re = r19->data[kstr] * abc->data[kstr].re;
            abc_im = r19->data[kstr] * abc->data[kstr].im;
            re = r17->data[kstr].re;
            im = r17->data[kstr].im;
            abc->data[kstr].re = abc_re * re - abc_im * im;
            abc->data[kstr].im = abc_re * im + abc_im * re;
        }
        break;

    case 1:
        kstr = r19->size[0] * r19->size[1];
        r19->size[0] = 1;
        r19->size[1] = f->size[1];
        emxEnsureCapacity_real_T(r19, kstr);
        loop_ub = f->size[0] * f->size[1];
        for (kstr = 0; kstr < loop_ub; kstr++) {
            r19->data[kstr] = 6.2831853071795862 * f->data[kstr];
        }

        b_freqs_cg(b1_data, b1_size, a1, r19, abc);
        kstr = r19->size[0] * r19->size[1];
        r19->size[0] = 1;
        r19->size[1] = f->size[1];
        emxEnsureCapacity_real_T(r19, kstr);
        loop_ub = f->size[0] * f->size[1];
        for (kstr = 0; kstr < loop_ub; kstr++) {
            r19->data[kstr] = 6.2831853071795862 * f->data[kstr];
        }

        b_freqs_cg(b2_data, b2_size, a2, r19, r17);
        rdivide(f, Fconverter, r19);
        b_sinc(r19);
        c_power(r19, r18);
        kstr = abc->size[0] * abc->size[1];
        abc->size[0] = 1;
        emxEnsureCapacity_creal_T(abc, kstr);
        kstr = abc->size[0];
        loop_ub = abc->size[1];
        loop_ub *= kstr;
        for (kstr = 0; kstr < loop_ub; kstr++) {
            abc_re = abc->data[kstr].re * r17->data[kstr].re - abc->data[kstr].im *
                     r17->data[kstr].im;
            abc_im = abc->data[kstr].re * r17->data[kstr].im + abc->data[kstr].im *
                     r17->data[kstr].re;
            abc->data[kstr].re = r18->data[kstr] * abc_re;
            abc->data[kstr].im = r18->data[kstr] * abc_im;
        }
        break;

    default:
        /*  Default to Rx */
        kstr = r19->size[0] * r19->size[1];
        r19->size[0] = 1;
        r19->size[1] = f->size[1];
        emxEnsureCapacity_real_T(r19, kstr);
        loop_ub = f->size[0] * f->size[1];
        for (kstr = 0; kstr < loop_ub; kstr++) {
            r19->data[kstr] = 6.2831853071795862 * f->data[kstr];
        }

        b_freqs_cg(b1_data, b1_size, a1, r19, abc);
        kstr = r19->size[0] * r19->size[1];
        r19->size[0] = 1;
        r19->size[1] = f->size[1];
        emxEnsureCapacity_real_T(r19, kstr);
        loop_ub = f->size[0] * f->size[1];
        for (kstr = 0; kstr < loop_ub; kstr++) {
            r19->data[kstr] = 6.2831853071795862 * f->data[kstr];
        }

        b_freqs_cg(b2_data, b2_size, a2, r19, r17);
        rdivide(f, Fconverter, r19);
        b_sinc(r19);
        c_power(r19, r18);
        kstr = abc->size[0] * abc->size[1];
        abc->size[0] = 1;
        emxEnsureCapacity_creal_T(abc, kstr);
        kstr = abc->size[0];
        loop_ub = abc->size[1];
        loop_ub *= kstr;
        for (kstr = 0; kstr < loop_ub; kstr++) {
            abc_re = abc->data[kstr].re * r17->data[kstr].re - abc->data[kstr].im *
                     r17->data[kstr].im;
            abc_im = abc->data[kstr].re * r17->data[kstr].im + abc->data[kstr].im *
                     r17->data[kstr].re;
            abc->data[kstr].re = r18->data[kstr] * abc_re;
            abc->data[kstr].im = r18->data[kstr] * abc_im;
        }
        break;
    }

    emxFree_real_T(&r19);
    emxFree_real_T(&r18);
    emxFree_creal_T(&r17);
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
    int k;
    static const double c_a[4] = { 0.0, 0.0, 0.0, 0.037037037037037035 };

    double ai;
    emxInit_creal_T(&a, 2);
    emxInit_real_T1(&b, 1);
    emxInit_real_T(&c, 2);
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
    k = 0;
    emxFree_real_T(&b_b);
    emxFree_creal_T(&b_a);
    emxFree_real_T(&c);
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
        ai = c_a[k] * y_im;
        if (ai == 0.0) {
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
 * Codegen workaround for fixed fi call requirements
 * Arguments    : const double tap_store[128]
 *                double taps[128]
 * Return Type  : void
 */
static void b_determineBestFractionLength(const double tap_store[128], double
        taps[128])
{
    double r[2048];
    int ixstart;
    double b_r[128];
    double dv16[128];
    double u;
    double e[16];
    double v;
    double dv17[128];
    short i57;
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
    int itmp;
    int ix;
    boolean_T exitg1;
    memset(&r[0], 0, sizeof(double) << 11);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        u = tap_store[ixstart] * 2.0;
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
                i57 = (short)u;
            } else {
                i57 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i57 = MAX_int16_T;
        } else {
            i57 = 0;
        }

        r[ixstart << 4] = (double)i57 * 0.5;
        b_r[ixstart] = r[ixstart << 4] - tap_store[ixstart];
        u = tap_store[ixstart] * 4.0;
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
                i57 = (short)u;
            } else {
                i57 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i57 = MAX_int16_T;
        } else {
            i57 = 0;
        }

        r[1 + (ixstart << 4)] = (double)i57 * 0.25;
    }

    d_abs(b_r, dv16);
    e[0] = c_sum(dv16);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        b_r[ixstart] = r[1 + (ixstart << 4)] - tap_store[ixstart];
        u = tap_store[ixstart] * 8.0;
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
                i57 = (short)u;
            } else {
                i57 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i57 = MAX_int16_T;
        } else {
            i57 = 0;
        }

        r[2 + (ixstart << 4)] = (double)i57 * 0.125;
    }

    d_abs(b_r, dv17);
    e[1] = c_sum(dv17);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        b_r[ixstart] = r[2 + (ixstart << 4)] - tap_store[ixstart];
        u = tap_store[ixstart] * 16.0;
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
                i57 = (short)u;
            } else {
                i57 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i57 = MAX_int16_T;
        } else {
            i57 = 0;
        }

        r[3 + (ixstart << 4)] = (double)i57 * 0.0625;
    }

    d_abs(b_r, dv18);
    e[2] = c_sum(dv18);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        b_r[ixstart] = r[3 + (ixstart << 4)] - tap_store[ixstart];
        u = tap_store[ixstart] * 32.0;
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
                i57 = (short)u;
            } else {
                i57 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i57 = MAX_int16_T;
        } else {
            i57 = 0;
        }

        r[4 + (ixstart << 4)] = (double)i57 * 0.03125;
    }

    d_abs(b_r, dv19);
    e[3] = c_sum(dv19);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        b_r[ixstart] = r[4 + (ixstart << 4)] - tap_store[ixstart];
        u = tap_store[ixstart] * 64.0;
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
                i57 = (short)u;
            } else {
                i57 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i57 = MAX_int16_T;
        } else {
            i57 = 0;
        }

        r[5 + (ixstart << 4)] = (double)i57 * 0.015625;
    }

    d_abs(b_r, dv20);
    e[4] = c_sum(dv20);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        b_r[ixstart] = r[5 + (ixstart << 4)] - tap_store[ixstart];
        u = tap_store[ixstart] * 128.0;
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
                i57 = (short)u;
            } else {
                i57 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i57 = MAX_int16_T;
        } else {
            i57 = 0;
        }

        r[6 + (ixstart << 4)] = (double)i57 * 0.0078125;
    }

    d_abs(b_r, dv21);
    e[5] = c_sum(dv21);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        b_r[ixstart] = r[6 + (ixstart << 4)] - tap_store[ixstart];
        u = tap_store[ixstart] * 256.0;
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
                i57 = (short)u;
            } else {
                i57 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i57 = MAX_int16_T;
        } else {
            i57 = 0;
        }

        r[7 + (ixstart << 4)] = (double)i57 * 0.00390625;
    }

    d_abs(b_r, dv22);
    e[6] = c_sum(dv22);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        b_r[ixstart] = r[7 + (ixstart << 4)] - tap_store[ixstart];
        u = tap_store[ixstart] * 512.0;
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
                i57 = (short)u;
            } else {
                i57 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i57 = MAX_int16_T;
        } else {
            i57 = 0;
        }

        r[8 + (ixstart << 4)] = (double)i57 * 0.001953125;
    }

    d_abs(b_r, dv23);
    e[7] = c_sum(dv23);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        b_r[ixstart] = r[8 + (ixstart << 4)] - tap_store[ixstart];
        u = tap_store[ixstart] * 1024.0;
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
                i57 = (short)u;
            } else {
                i57 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i57 = MAX_int16_T;
        } else {
            i57 = 0;
        }

        r[9 + (ixstart << 4)] = (double)i57 * 0.0009765625;
    }

    d_abs(b_r, dv24);
    e[8] = c_sum(dv24);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        b_r[ixstart] = r[9 + (ixstart << 4)] - tap_store[ixstart];
        u = tap_store[ixstart] * 2048.0;
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
                i57 = (short)u;
            } else {
                i57 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i57 = MAX_int16_T;
        } else {
            i57 = 0;
        }

        r[10 + (ixstart << 4)] = (double)i57 * 0.00048828125;
    }

    d_abs(b_r, dv25);
    e[9] = c_sum(dv25);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        b_r[ixstart] = r[10 + (ixstart << 4)] - tap_store[ixstart];
        u = tap_store[ixstart] * 4096.0;
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
                i57 = (short)u;
            } else {
                i57 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i57 = MAX_int16_T;
        } else {
            i57 = 0;
        }

        r[11 + (ixstart << 4)] = (double)i57 * 0.000244140625;
    }

    d_abs(b_r, dv26);
    e[10] = c_sum(dv26);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        b_r[ixstart] = r[11 + (ixstart << 4)] - tap_store[ixstart];
        u = tap_store[ixstart] * 8192.0;
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
                i57 = (short)u;
            } else {
                i57 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i57 = MAX_int16_T;
        } else {
            i57 = 0;
        }

        r[12 + (ixstart << 4)] = (double)i57 * 0.0001220703125;
    }

    d_abs(b_r, dv27);
    e[11] = c_sum(dv27);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        b_r[ixstart] = r[12 + (ixstart << 4)] - tap_store[ixstart];
        u = tap_store[ixstart] * 16384.0;
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
                i57 = (short)u;
            } else {
                i57 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i57 = MAX_int16_T;
        } else {
            i57 = 0;
        }

        r[13 + (ixstart << 4)] = (double)i57 * 6.103515625E-5;
    }

    d_abs(b_r, dv28);
    e[12] = c_sum(dv28);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        b_r[ixstart] = r[13 + (ixstart << 4)] - tap_store[ixstart];
        u = tap_store[ixstart] * 32768.0;
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
                i57 = (short)u;
            } else {
                i57 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i57 = MAX_int16_T;
        } else {
            i57 = 0;
        }

        r[14 + (ixstart << 4)] = (double)i57 * 3.0517578125E-5;
    }

    d_abs(b_r, dv29);
    e[13] = c_sum(dv29);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        b_r[ixstart] = r[14 + (ixstart << 4)] - tap_store[ixstart];
        u = tap_store[ixstart] * 65536.0;
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
                i57 = (short)u;
            } else {
                i57 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i57 = MAX_int16_T;
        } else {
            i57 = 0;
        }

        r[15 + (ixstart << 4)] = (double)i57 * 1.52587890625E-5;
    }

    d_abs(b_r, dv30);
    e[14] = c_sum(dv30);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        b_r[ixstart] = r[15 + (ixstart << 4)] - tap_store[ixstart];
    }

    d_abs(b_r, dv31);
    e[15] = c_sum(dv31);
    ixstart = 1;
    u = e[0];
    itmp = 0;
    if (rtIsNaN(e[0])) {
        ix = 2;
        exitg1 = false;
        while ((!exitg1) && (ix < 17)) {
            ixstart = ix;
            if (!rtIsNaN(e[ix - 1])) {
                u = e[ix - 1];
                itmp = ix - 1;
                exitg1 = true;
            } else {
                ix++;
            }
        }
    }

    if (ixstart < 16) {
        while (ixstart + 1 < 17) {
            if (e[ixstart] < u) {
                u = e[ixstart];
                itmp = ixstart;
            }

            ixstart++;
        }
    }

    for (ixstart = 0; ixstart < 128; ixstart++) {
        taps[ixstart] = r[itmp + (ixstart << 4)];
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
    int i64;
    creal_T dcv2[2048];
    double digw;
    double b_digw[2048];
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
    for (i64 = 0; i64 < 2048; i64++) {
        w[i64] = options->w[i64];
        digw = 6.2831853071795862 * options->w[i64] / options->Fs;
        dcv2[i64].re = digw * 0.0;
        dcv2[i64].im = digw;
        b_digw[i64] = digw;
    }

    b_exp(dcv2);
    polyval(b, dcv2, h);
    for (i64 = 0; i64 < 2048; i64++) {
        dcv2[i64].re = 14.0 * (b_digw[i64] * 0.0);
        dcv2[i64].im = 14.0 * b_digw[i64];
    }

    b_exp(dcv2);
    for (i64 = 0; i64 < 2048; i64++) {
        h_re = h[i64].re;
        if (dcv2[i64].im == 0.0) {
            if (h[i64].im == 0.0) {
                h[i64].re /= dcv2[i64].re;
                h[i64].im = 0.0;
            } else if (h[i64].re == 0.0) {
                h[i64].re = 0.0;
                h[i64].im /= dcv2[i64].re;
            } else {
                h[i64].re /= dcv2[i64].re;
                h[i64].im /= dcv2[i64].re;
            }
        } else if (dcv2[i64].re == 0.0) {
            if (h[i64].re == 0.0) {
                h[i64].re = h[i64].im / dcv2[i64].im;
                h[i64].im = 0.0;
            } else if (h[i64].im == 0.0) {
                h[i64].re = 0.0;
                h[i64].im = -(h_re / dcv2[i64].im);
            } else {
                h[i64].re = h[i64].im / dcv2[i64].im;
                h[i64].im = -(h_re / dcv2[i64].im);
            }
        } else {
            brm = fabs(dcv2[i64].re);
            digw = fabs(dcv2[i64].im);
            if (brm > digw) {
                digw = dcv2[i64].im / dcv2[i64].re;
                d = dcv2[i64].re + digw * dcv2[i64].im;
                h[i64].re = (h[i64].re + digw * h[i64].im) / d;
                h[i64].im = (h[i64].im - digw * h_re) / d;
            } else if (digw == brm) {
                if (dcv2[i64].re > 0.0) {
                    digw = 0.5;
                } else {
                    digw = -0.5;
                }

                if (dcv2[i64].im > 0.0) {
                    d = 0.5;
                } else {
                    d = -0.5;
                }

                h[i64].re = (h[i64].re * digw + h[i64].im * d) / brm;
                h[i64].im = (h[i64].im * digw - h_re * d) / brm;
            } else {
                digw = dcv2[i64].re / dcv2[i64].im;
                d = dcv2[i64].im + digw * dcv2[i64].re;
                h[i64].re = (digw * h[i64].re + h[i64].im) / d;
                h[i64].im = (digw * h[i64].im - h_re) / d;
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
    emxArray_real_T *b_h;
    int i52;
    emxArray_real_T *c_h;
    double b_ff[4];
    double x;
    boolean_T b_valid;
    int h_idx_0;
    double d2;
    int i53;
    int i54;
    int loop_ub;
    emxInit_real_T(&grid, 2);
    emxInit_real_T(&des, 2);
    emxInit_real_T(&wt, 2);
    emxInit_real_T(&b_h, 2);

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
    interp1(frequencies, amplitudes, grid, des);
    interp1(frequencies, weights, grid, wt);

    /*  Workaround */
    /* ftype = 2; */
    /* sign_val = 1; */
    /*  Always bandpass designs */
    /*  cast to enforce precision rules */
    /*  Call actual design algorithm */
    rdivide(grid, 2.0, b_h);
    emxFree_real_T(&grid);
    for (i52 = 0; i52 < 4; i52++) {
        b_ff[i52] = ff[i52] / 2.0;
    }

    emxInit_real_T(&c_h, 2);
    remezm(order + 1.0, b_ff, b_h, des, wt, c_h, &x, &b_valid);
    h_idx_0 = c_h->size[0] * c_h->size[1];
    i52 = h->size[0] * h->size[1];
    h->size[0] = 1;
    h->size[1] = h_idx_0;
    emxEnsureCapacity_real_T(h, i52);
    emxFree_real_T(&wt);
    emxFree_real_T(&des);
    for (i52 = 0; i52 < h_idx_0; i52++) {
        h->data[h->size[0] * i52] = c_h->data[i52];
    }

    emxFree_real_T(&c_h);

    /*  make it a row */
    d2 = (double)h->size[1] - rt_remd_snf(order + 1.0, 2.0);
    if (1.0 > d2) {
        i52 = 1;
        h_idx_0 = 1;
        i53 = 0;
    } else {
        i52 = (int)d2;
        h_idx_0 = -1;
        i53 = 1;
    }

    i54 = b_h->size[0] * b_h->size[1];
    b_h->size[0] = 1;
    b_h->size[1] = (h->size[1] + div_s32_floor(i53 - i52, h_idx_0)) + 1;
    emxEnsureCapacity_real_T(b_h, i54);
    loop_ub = h->size[1];
    for (i54 = 0; i54 < loop_ub; i54++) {
        b_h->data[b_h->size[0] * i54] = h->data[h->size[0] * i54];
    }

    loop_ub = div_s32_floor(i53 - i52, h_idx_0);
    for (i53 = 0; i53 <= loop_ub; i53++) {
        b_h->data[b_h->size[0] * (i53 + h->size[1])] = h->data[(i52 + h_idx_0 * i53)
                - 1];
    }

    i52 = h->size[0] * h->size[1];
    h->size[0] = 1;
    h->size[1] = b_h->size[1];
    emxEnsureCapacity_real_T(h, i52);
    loop_ub = b_h->size[1];
    for (i52 = 0; i52 < loop_ub; i52++) {
        h->data[h->size[0] * i52] = b_h->data[b_h->size[0] * i52];
    }

    if (1 > h->size[1]) {
        i52 = 1;
        h_idx_0 = 1;
        i53 = 0;
    } else {
        i52 = h->size[1];
        h_idx_0 = -1;
        i53 = 1;
    }

    i54 = b_h->size[0] * b_h->size[1];
    b_h->size[0] = 1;
    b_h->size[1] = div_s32_floor(i53 - i52, h_idx_0) + 1;
    emxEnsureCapacity_real_T(b_h, i54);
    loop_ub = div_s32_floor(i53 - i52, h_idx_0);
    for (i53 = 0; i53 <= loop_ub; i53++) {
        b_h->data[b_h->size[0] * i53] = h->data[(i52 + h_idx_0 * i53) - 1];
    }

    i52 = h->size[0] * h->size[1];
    h->size[0] = 1;
    h->size[1] = b_h->size[1];
    emxEnsureCapacity_real_T(h, i52);
    loop_ub = b_h->size[1];
    for (i52 = 0; i52 < loop_ub; i52++) {
        h->data[h->size[0] * i52] = b_h->data[b_h->size[0] * i52];
    }

    emxFree_real_T(&b_h);
    *valid = b_valid;
    *err = fabs(x);
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
    int i40;
    int loop_ub;
    emxArray_creal_T *y;
    boolean_T b11;
    emxArray_real_T c_b_data;
    int k;
    double a_re;
    double a_im;
    double s_re;
    double s_im;
    emxInit_creal_T(&s, 2);
    emxInit_creal_T(&b_a, 2);
    removeTrailingZero(b_data, b_size, a, b_b_data, b_b_size, b_a);
    i40 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = w->size[1];
    emxEnsureCapacity_creal_T(s, i40);
    loop_ub = w->size[0] * w->size[1];
    for (i40 = 0; i40 < loop_ub; i40++) {
        s->data[i40].re = w->data[i40] * 0.0;
        s->data[i40].im = w->data[i40];
    }

    emxInit_creal_T(&y, 2);
    i40 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = s->size[1];
    emxEnsureCapacity_creal_T(y, i40);
    if ((y->size[1] == 0) || (b_a->size[1] == 0)) {
        b11 = true;
    } else {
        b11 = false;
    }

    if (!b11) {
        i40 = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity_creal_T(y, i40);
        loop_ub = y->size[1];
        for (i40 = 0; i40 < loop_ub; i40++) {
            y->data[y->size[0] * i40] = b_a->data[0];
        }

        for (k = 0; k <= b_a->size[1] - 2; k++) {
            i40 = y->size[0] * y->size[1];
            y->size[0] = 1;
            y->size[1] = s->size[1];
            emxEnsureCapacity_creal_T(y, i40);
            a_re = b_a->data[k + 1].re;
            a_im = b_a->data[k + 1].im;
            loop_ub = s->size[0] * s->size[1];
            for (i40 = 0; i40 < loop_ub; i40++) {
                s_re = s->data[i40].re * y->data[i40].re - s->data[i40].im * y->data[i40]
                       .im;
                s_im = s->data[i40].re * y->data[i40].im + s->data[i40].im * y->data[i40]
                       .re;
                y->data[i40].re = s_re + a_re;
                y->data[i40].im = s_im + a_im;
            }
        }
    }

    c_b_data.data = (double *)&b_b_data;
    c_b_data.size = (int *)&b_b_size;
    c_b_data.allocatedSize = 4;
    c_b_data.numDimensions = 2;
    c_b_data.canFreeData = false;
    j_polyval(&c_b_data, s, b_a);
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
    boolean_T b_bool;
    int ix;
    int exitg1;
    emxArray_creal_T *d2;
    emxArray_creal_T *d3;
    static const char cv35[4] = { '2', '1', '1', '1' };

    double u[15];
    double tmp_data[29];
    int tmp_size[2];
    double b_u[30];
    double c_u[14];
    double d_u[60];
    double e_u[45];
    double f_u[21];
    double g_u[90];
    emxArray_real_T b_tmp_data;
    int iy;
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
    ix = 1;
    do {
        exitg1 = 0;
        if (ix < 5) {
            if (enables[ix - 1] != '1') {
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
            if (ix + 1 < 5) {
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
                if (ix + 1 < 5) {
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
                    if (ix + 1 < 5) {
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
                        if (ix + 1 < 5) {
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
                            ix = 5;
                        } else {
                            b_bool = false;
                            ix = 0;
                            do {
                                exitg1 = 0;
                                if (ix + 1 < 5) {
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
                                    if (ix + 1 < 5) {
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
                                        if (ix + 1 < 5) {
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
                                            if (ix + 1 < 5) {
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
                                                if (ix + 1 < 5) {
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
                                                    if (ix + 1 < 5) {
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
        for (ix = 0; ix < 7; ix++) {
            h_u[ix] = 0.0;
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
        b_tmp_data.data = (double *)&tmp_data;
        b_tmp_data.size = (int *)&tmp_size;
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

        o_freqz_cg(*(double (*)[29])&b_u[0], w, Fs, combinedResponse);
        for (ix = 0; ix < 7; ix++) {
            h_u[ix] = 0.0;
        }

        ix = 0;
        iy = 0;
        for (k = 0; k < 7; k++) {
            h_u[iy] = hb2_coeff[ix];
            ix++;
            iy++;
        }

        m_freqz_cg(h_u, w, Fs, d2);
        ix = combinedResponse->size[0] * combinedResponse->size[1];
        combinedResponse->size[0] = 1;
        emxEnsureCapacity_creal_T(combinedResponse, ix);
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
        memset(&b_u[0], 0, 30U * sizeof(double));
        ix = 0;
        iy = 0;
        for (k = 0; k < 15; k++) {
            b_u[iy] = hb1_coeff[ix];
            ix++;
            iy += 2;
        }

        o_freqz_cg(*(double (*)[29])&b_u[0], w, Fs, combinedResponse);
        c_us(hb3_coeff_data, hb3_coeff_size, tmp_data, tmp_size);
        b_tmp_data.data = (double *)&tmp_data;
        b_tmp_data.size = (int *)&tmp_size;
        b_tmp_data.allocatedSize = 29;
        b_tmp_data.numDimensions = 2;
        b_tmp_data.canFreeData = false;
        n_freqz_cg(&b_tmp_data, w, Fs, d2);
        ix = combinedResponse->size[0] * combinedResponse->size[1];
        combinedResponse->size[0] = 1;
        emxEnsureCapacity_creal_T(combinedResponse, ix);
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
        memset(&c_u[0], 0, 14U * sizeof(double));
        ix = 0;
        iy = 0;
        for (k = 0; k < 7; k++) {
            c_u[iy] = hb2_coeff[ix];
            ix++;
            iy += 2;
        }

        p_freqz_cg(*(double (*)[13])&c_u[0], w, Fs, combinedResponse);
        c_us(hb3_coeff_data, hb3_coeff_size, tmp_data, tmp_size);
        b_tmp_data.data = (double *)&tmp_data;
        b_tmp_data.size = (int *)&tmp_size;
        b_tmp_data.allocatedSize = 29;
        b_tmp_data.numDimensions = 2;
        b_tmp_data.canFreeData = false;
        n_freqz_cg(&b_tmp_data, w, Fs, d2);
        ix = combinedResponse->size[0] * combinedResponse->size[1];
        combinedResponse->size[0] = 1;
        emxEnsureCapacity_creal_T(combinedResponse, ix);
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

        q_freqz_cg(*(double (*)[57])&d_u[0], w, Fs, combinedResponse);
        memset(&c_u[0], 0, 14U * sizeof(double));
        ix = 0;
        iy = 0;
        for (k = 0; k < 7; k++) {
            c_u[iy] = hb2_coeff[ix];
            ix++;
            iy += 2;
        }

        p_freqz_cg(*(double (*)[13])&c_u[0], w, Fs, d2);
        c_us(hb3_coeff_data, hb3_coeff_size, tmp_data, tmp_size);
        b_tmp_data.data = (double *)&tmp_data;
        b_tmp_data.size = (int *)&tmp_size;
        b_tmp_data.allocatedSize = 29;
        b_tmp_data.numDimensions = 2;
        b_tmp_data.canFreeData = false;
        n_freqz_cg(&b_tmp_data, w, Fs, d3);
        ix = combinedResponse->size[0] * combinedResponse->size[1];
        combinedResponse->size[0] = 1;
        emxEnsureCapacity_creal_T(combinedResponse, ix);
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
        b_tmp_data.data = (double *)&tmp_data;
        b_tmp_data.size = (int *)&tmp_size;
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

        r_freqz_cg(*(double (*)[43])&e_u[0], w, Fs, combinedResponse);
        c_us(dec_int3_coeff_data, dec_int3_coeff_size, tmp_data, tmp_size);
        b_tmp_data.data = (double *)&tmp_data;
        b_tmp_data.size = (int *)&tmp_size;
        b_tmp_data.allocatedSize = 29;
        b_tmp_data.numDimensions = 2;
        b_tmp_data.canFreeData = false;
        n_freqz_cg(&b_tmp_data, w, Fs, d2);
        ix = combinedResponse->size[0] * combinedResponse->size[1];
        combinedResponse->size[0] = 1;
        emxEnsureCapacity_creal_T(combinedResponse, ix);
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
        memset(&f_u[0], 0, 21U * sizeof(double));
        ix = 0;
        iy = 0;
        for (k = 0; k < 7; k++) {
            f_u[iy] = hb2_coeff[ix];
            ix++;
            iy += 3;
        }

        s_freqz_cg(*(double (*)[19])&f_u[0], w, Fs, combinedResponse);
        c_us(dec_int3_coeff_data, dec_int3_coeff_size, tmp_data, tmp_size);
        b_tmp_data.data = (double *)&tmp_data;
        b_tmp_data.size = (int *)&tmp_size;
        b_tmp_data.allocatedSize = 29;
        b_tmp_data.numDimensions = 2;
        b_tmp_data.canFreeData = false;
        n_freqz_cg(&b_tmp_data, w, Fs, d2);
        ix = combinedResponse->size[0] * combinedResponse->size[1];
        combinedResponse->size[0] = 1;
        emxEnsureCapacity_creal_T(combinedResponse, ix);
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
        /*  Dec/Int3,Hb2,Hb1 {Hm4,Hm2c34,Hm1} */
        memset(&g_u[0], 0, 90U * sizeof(double));
        ix = 0;
        iy = 0;
        for (k = 0; k < 15; k++) {
            g_u[iy] = hb1_coeff[ix];
            ix++;
            iy += 6;
        }

        t_freqz_cg(*(double (*)[85])&g_u[0], w, Fs, combinedResponse);
        memset(&f_u[0], 0, 21U * sizeof(double));
        ix = 0;
        iy = 0;
        for (k = 0; k < 7; k++) {
            f_u[iy] = hb2_coeff[ix];
            ix++;
            iy += 3;
        }

        s_freqz_cg(*(double (*)[19])&f_u[0], w, Fs, d2);
        c_us(dec_int3_coeff_data, dec_int3_coeff_size, tmp_data, tmp_size);
        b_tmp_data.data = (double *)&tmp_data;
        b_tmp_data.size = (int *)&tmp_size;
        b_tmp_data.allocatedSize = 29;
        b_tmp_data.numDimensions = 2;
        b_tmp_data.canFreeData = false;
        n_freqz_cg(&b_tmp_data, w, Fs, d3);
        ix = combinedResponse->size[0] * combinedResponse->size[1];
        combinedResponse->size[0] = 1;
        emxEnsureCapacity_creal_T(combinedResponse, ix);
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
        } else if ((eint == 1) && (t < 0.75)) {
            f = log(2.0 * t) / 0.69314718055994529;
        } else {
            f = log(t) / 0.69314718055994529 + (double)eint;
        }
    } else {
        f = x;
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
    double x_im;
    for (i11 = 0; i11 < 2048; i11++) {
        y[i11].re = p[0];
        y[i11].im = 0.0;
    }

    for (k = 0; k < 6; k++) {
        for (i11 = 0; i11 < 2048; i11++) {
            x_im = x[i11].re * y[i11].im + x[i11].im * y[i11].re;
            y[i11].re = (x[i11].re * y[i11].re - x[i11].im * y[i11].im) + p[k + 1];
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
 * Arguments    : const emxArray_creal_T *x
 *                const emxArray_creal_T *y
 *                emxArray_creal_T *z
 * Return Type  : void
 */
static void b_rdivide(const emxArray_creal_T *x, const emxArray_creal_T *y,
                      emxArray_creal_T *z)
{
    int i29;
    int loop_ub;
    double x_re;
    double x_im;
    double y_re;
    double y_im;
    double brm;
    double bim;
    double s;
    i29 = z->size[0] * z->size[1];
    z->size[0] = 1;
    z->size[1] = x->size[1];
    emxEnsureCapacity_creal_T(z, i29);
    loop_ub = x->size[0] * x->size[1];
    for (i29 = 0; i29 < loop_ub; i29++) {
        x_re = x->data[i29].re;
        x_im = x->data[i29].im;
        y_re = y->data[i29].re;
        y_im = y->data[i29].im;
        if (y_im == 0.0) {
            if (x_im == 0.0) {
                z->data[i29].re = x_re / y_re;
                z->data[i29].im = 0.0;
            } else if (x_re == 0.0) {
                z->data[i29].re = 0.0;
                z->data[i29].im = x_im / y_re;
            } else {
                z->data[i29].re = x_re / y_re;
                z->data[i29].im = x_im / y_re;
            }
        } else if (y_re == 0.0) {
            if (x_re == 0.0) {
                z->data[i29].re = x_im / y_im;
                z->data[i29].im = 0.0;
            } else if (x_im == 0.0) {
                z->data[i29].re = 0.0;
                z->data[i29].im = -(x_re / y_im);
            } else {
                z->data[i29].re = x_im / y_im;
                z->data[i29].im = -(x_re / y_im);
            }
        } else {
            brm = fabs(y_re);
            bim = fabs(y_im);
            if (brm > bim) {
                s = y_im / y_re;
                bim = y_re + s * y_im;
                z->data[i29].re = (x_re + s * x_im) / bim;
                z->data[i29].im = (x_im - s * x_re) / bim;
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

                z->data[i29].re = (x_re * s + x_im * bim) / brm;
                z->data[i29].im = (x_im * s - x_re * bim) / brm;
            } else {
                s = y_re / y_im;
                bim = y_im + s * y_re;
                z->data[i29].re = (s * x_re + x_im) / bim;
                z->data[i29].im = (s * x_im - x_re) / bim;
            }
        }
    }
}

/*
 * Arguments    : emxArray_real_T *x
 * Return Type  : void
 */
static void b_sinc(emxArray_real_T *x)
{
    int i72;
    int k;
    i72 = x->size[1];
    for (k = 0; k < i72; k++) {
        if (fabs(x->data[k]) < 1.0020841800044864E-292) {
            x->data[k] = 1.0;
        } else {
            x->data[k] *= 3.1415926535897931;
            x->data[k] = sin(x->data[k]) / x->data[k];
        }
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
    int i62;
    int k;
    double x_re;
    double x_im;
    if (!(incx < 1)) {
        i62 = ix0 + incx * (n - 1);
        for (k = ix0; k <= i62; k += incx) {
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
    double scale;
    double f2s;
    double x;
    double fs_re;
    double fs_im;
    double gs_re;
    double gs_im;
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
        x = g2;
        if (1.0 > g2) {
            x = 1.0;
        }

        if (scale <= x * 2.0041683600089728E-292) {
            if ((f.re == 0.0) && (f.im == 0.0)) {
                *cs = 0.0;
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
            }
        } else {
            f2s = sqrt(1.0 + g2 / scale);
            fs_re *= f2s;
            fs_im *= f2s;
            *cs = 1.0 / f2s;
            g2 += scale;
            fs_re /= g2;
            fs_im /= g2;
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
    int i5;
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
    double ai;
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
    i5 = c->size[0] * c->size[1];
    c->size[0] = 1;
    c->size[1] = 1;
    emxEnsureCapacity_real_T(c, i5);
    c->data[0] = 1.0;
    tmp_size[0] = 1;
    tmp_size[1] = 1;
    tmp_data[0].re = -1.0;
    tmp_data[0].im = 0.0;
    b_tmp_size[0] = 1;
    b_tmp_data[0] = 1.0;
    emxInit_creal_T(&a, 2);
    emxInit_real_T1(&b, 1);
    emxInit_creal_T(&r1, 2);
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
    for (i5 = 0; i5 < loop_ub; i5++) {
        den_data[i5] = r1->data[i5];
    }

    emxFree_creal_T(&r1);

    /*  This internal function returns more exact numerator vectors */
    /*  for the num/den case. */
    /*  Wn input is two element band edge vector */
    /* --------------------------------- */
    /*  lowpass */
    d = (0.0 * den_data[0].re - 0.0 * den_data[0].im) + den_data[1].re;
    y_im = (0.0 * den_data[0].im + 0.0 * den_data[0].re) + den_data[1].im;
    for (i5 = 0; i5 < 2; i5++) {
        ar = (double)i5 * d;
        ai = (double)i5 * y_im;
        if ((!(ai == 0.0)) && (ar == 0.0)) {
            ar = 0.0;
        }

        num[i5] = ar;
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
    emxEnsureCapacity_real_T(y, k);
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
    emxArray_real_T *r25;
    int b_abc;
    int loop_ub;
    emxArray_creal_T *r26;
    emxArray_real_T *r27;
    double abc_re;
    double abc_im;
    emxInit_real_T(&r25, 2);
    b_abc = r25->size[0] * r25->size[1];
    r25->size[0] = 1;
    r25->size[1] = f->size[1];
    emxEnsureCapacity_real_T(r25, b_abc);
    loop_ub = f->size[0] * f->size[1];
    for (b_abc = 0; b_abc < loop_ub; b_abc++) {
        r25->data[b_abc] = 6.2831853071795862 * f->data[b_abc];
    }

    b_freqs_cg(b1_data, b1_size, a1, r25, abc);
    b_abc = r25->size[0] * r25->size[1];
    r25->size[0] = 1;
    r25->size[1] = f->size[1];
    emxEnsureCapacity_real_T(r25, b_abc);
    loop_ub = f->size[0] * f->size[1];
    for (b_abc = 0; b_abc < loop_ub; b_abc++) {
        r25->data[b_abc] = 6.2831853071795862 * f->data[b_abc];
    }

    emxInit_creal_T(&r26, 2);
    emxInit_real_T(&r27, 2);
    b_freqs_cg(b2_data, b2_size, a2, r25, r26);
    rdivide(f, Fconverter, r25);
    b_sinc(r25);
    c_power(r25, r27);
    b_abc = abc->size[0] * abc->size[1];
    abc->size[0] = 1;
    emxEnsureCapacity_creal_T(abc, b_abc);
    b_abc = abc->size[0];
    loop_ub = abc->size[1];
    loop_ub *= b_abc;
    emxFree_real_T(&r25);
    for (b_abc = 0; b_abc < loop_ub; b_abc++) {
        abc_re = abc->data[b_abc].re * r26->data[b_abc].re - abc->data[b_abc].im *
                 r26->data[b_abc].im;
        abc_im = abc->data[b_abc].re * r26->data[b_abc].im + abc->data[b_abc].im *
                 r26->data[b_abc].re;
        abc->data[b_abc].re = r27->data[b_abc] * abc_re;
        abc->data[b_abc].im = r27->data[b_abc] * abc_im;
    }

    emxFree_real_T(&r27);
    emxFree_creal_T(&r26);
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
    int i65;
    creal_T dcv3[2048];
    double digw;
    double b_digw[2048];
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
    for (i65 = 0; i65 < 2048; i65++) {
        w[i65] = options->w[i65];
        digw = 6.2831853071795862 * options->w[i65] / options->Fs;
        dcv3[i65].re = digw * 0.0;
        dcv3[i65].im = digw;
        b_digw[i65] = digw;
    }

    b_exp(dcv3);
    b_polyval(b, dcv3, h);
    for (i65 = 0; i65 < 2048; i65++) {
        dcv3[i65].re = 6.0 * (b_digw[i65] * 0.0);
        dcv3[i65].im = 6.0 * b_digw[i65];
    }

    b_exp(dcv3);
    for (i65 = 0; i65 < 2048; i65++) {
        h_re = h[i65].re;
        if (dcv3[i65].im == 0.0) {
            if (h[i65].im == 0.0) {
                h[i65].re /= dcv3[i65].re;
                h[i65].im = 0.0;
            } else if (h[i65].re == 0.0) {
                h[i65].re = 0.0;
                h[i65].im /= dcv3[i65].re;
            } else {
                h[i65].re /= dcv3[i65].re;
                h[i65].im /= dcv3[i65].re;
            }
        } else if (dcv3[i65].re == 0.0) {
            if (h[i65].re == 0.0) {
                h[i65].re = h[i65].im / dcv3[i65].im;
                h[i65].im = 0.0;
            } else if (h[i65].im == 0.0) {
                h[i65].re = 0.0;
                h[i65].im = -(h_re / dcv3[i65].im);
            } else {
                h[i65].re = h[i65].im / dcv3[i65].im;
                h[i65].im = -(h_re / dcv3[i65].im);
            }
        } else {
            brm = fabs(dcv3[i65].re);
            digw = fabs(dcv3[i65].im);
            if (brm > digw) {
                digw = dcv3[i65].im / dcv3[i65].re;
                d = dcv3[i65].re + digw * dcv3[i65].im;
                h[i65].re = (h[i65].re + digw * h[i65].im) / d;
                h[i65].im = (h[i65].im - digw * h_re) / d;
            } else if (digw == brm) {
                if (dcv3[i65].re > 0.0) {
                    digw = 0.5;
                } else {
                    digw = -0.5;
                }

                if (dcv3[i65].im > 0.0) {
                    d = 0.5;
                } else {
                    d = -0.5;
                }

                h[i65].re = (h[i65].re * digw + h[i65].im * d) / brm;
                h[i65].im = (h[i65].im * digw - h_re * d) / brm;
            } else {
                digw = dcv3[i65].re / dcv3[i65].im;
                d = dcv3[i65].im + digw * dcv3[i65].re;
                h[i65].re = (digw * h[i65].re + h[i65].im) / d;
                h[i65].im = (digw * h[i65].im - h_re) / d;
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
    boolean_T b_bool;
    int iy;
    int exitg1;
    emxArray_creal_T *d2;
    emxArray_creal_T *d3;
    static const char cv48[4] = { '2', '1', '1', '1' };

    double u[15];
    int stages;
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
    emxArray_real_T *i_u;
    int k;
    int c;
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
    iy = 1;
    do {
        exitg1 = 0;
        if (iy < 5) {
            if (enables[iy - 1] != '1') {
                exitg1 = 1;
            } else {
                iy++;
            }
        } else {
            b_bool = true;
            exitg1 = 1;
        }
    } while (exitg1 == 0);

    if (b_bool) {
        iy = 0;
    } else {
        b_bool = false;
        iy = 0;
        do {
            exitg1 = 0;
            if (iy + 1 < 5) {
                if (enables[iy] != cv48[iy]) {
                    exitg1 = 1;
                } else {
                    iy++;
                }
            } else {
                b_bool = true;
                exitg1 = 1;
            }
        } while (exitg1 == 0);

        if (b_bool) {
            iy = 1;
        } else {
            b_bool = false;
            iy = 0;
            do {
                exitg1 = 0;
                if (iy + 1 < 5) {
                    if (enables[iy] != cv49[iy]) {
                        exitg1 = 1;
                    } else {
                        iy++;
                    }
                } else {
                    b_bool = true;
                    exitg1 = 1;
                }
            } while (exitg1 == 0);

            if (b_bool) {
                iy = 2;
            } else {
                b_bool = false;
                iy = 0;
                do {
                    exitg1 = 0;
                    if (iy + 1 < 5) {
                        if (enables[iy] != cv50[iy]) {
                            exitg1 = 1;
                        } else {
                            iy++;
                        }
                    } else {
                        b_bool = true;
                        exitg1 = 1;
                    }
                } while (exitg1 == 0);

                if (b_bool) {
                    iy = 3;
                } else {
                    b_bool = false;
                    iy = 0;
                    do {
                        exitg1 = 0;
                        if (iy + 1 < 5) {
                            if (enables[iy] != cv51[iy]) {
                                exitg1 = 1;
                            } else {
                                iy++;
                            }
                        } else {
                            b_bool = true;
                            exitg1 = 1;
                        }
                    } while (exitg1 == 0);

                    if (b_bool) {
                        iy = 4;
                    } else {
                        b_bool = false;
                        iy = 0;
                        do {
                            exitg1 = 0;
                            if (iy + 1 < 5) {
                                if (enables[iy] != cv52[iy]) {
                                    exitg1 = 1;
                                } else {
                                    iy++;
                                }
                            } else {
                                b_bool = true;
                                exitg1 = 1;
                            }
                        } while (exitg1 == 0);

                        if (b_bool) {
                            iy = 5;
                        } else {
                            b_bool = false;
                            iy = 0;
                            do {
                                exitg1 = 0;
                                if (iy + 1 < 5) {
                                    if (enables[iy] != cv53[iy]) {
                                        exitg1 = 1;
                                    } else {
                                        iy++;
                                    }
                                } else {
                                    b_bool = true;
                                    exitg1 = 1;
                                }
                            } while (exitg1 == 0);

                            if (b_bool) {
                                iy = 6;
                            } else {
                                b_bool = false;
                                iy = 0;
                                do {
                                    exitg1 = 0;
                                    if (iy + 1 < 5) {
                                        if (enables[iy] != cv54[iy]) {
                                            exitg1 = 1;
                                        } else {
                                            iy++;
                                        }
                                    } else {
                                        b_bool = true;
                                        exitg1 = 1;
                                    }
                                } while (exitg1 == 0);

                                if (b_bool) {
                                    iy = 7;
                                } else {
                                    b_bool = false;
                                    iy = 0;
                                    do {
                                        exitg1 = 0;
                                        if (iy + 1 < 5) {
                                            if (enables[iy] != cv55[iy]) {
                                                exitg1 = 1;
                                            } else {
                                                iy++;
                                            }
                                        } else {
                                            b_bool = true;
                                            exitg1 = 1;
                                        }
                                    } while (exitg1 == 0);

                                    if (b_bool) {
                                        iy = 8;
                                    } else {
                                        b_bool = false;
                                        iy = 0;
                                        do {
                                            exitg1 = 0;
                                            if (iy + 1 < 5) {
                                                if (enables[iy] != cv56[iy]) {
                                                    exitg1 = 1;
                                                } else {
                                                    iy++;
                                                }
                                            } else {
                                                b_bool = true;
                                                exitg1 = 1;
                                            }
                                        } while (exitg1 == 0);

                                        if (b_bool) {
                                            iy = 9;
                                        } else {
                                            b_bool = false;
                                            iy = 0;
                                            do {
                                                exitg1 = 0;
                                                if (iy + 1 < 5) {
                                                    if (enables[iy] != cv57[iy]) {
                                                        exitg1 = 1;
                                                    } else {
                                                        iy++;
                                                    }
                                                } else {
                                                    b_bool = true;
                                                    exitg1 = 1;
                                                }
                                            } while (exitg1 == 0);

                                            if (b_bool) {
                                                iy = 10;
                                            } else {
                                                b_bool = false;
                                                iy = 0;
                                                do {
                                                    exitg1 = 0;
                                                    if (iy + 1 < 5) {
                                                        if (enables[iy] != cv58[iy]) {
                                                            exitg1 = 1;
                                                        } else {
                                                            iy++;
                                                        }
                                                    } else {
                                                        b_bool = true;
                                                        exitg1 = 1;
                                                    }
                                                } while (exitg1 == 0);

                                                if (b_bool) {
                                                    iy = 11;
                                                } else {
                                                    iy = -1;
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
    switch (iy) {
    case 0:
        /*  only FIR */
        k_freqz_cg(w, Fs, combinedResponse);
        stages = 1;
        break;

    case 1:
        /*  Hb1 */
        memset(&u[0], 0, 15U * sizeof(double));
        stages = 0;
        iy = 0;
        for (k = 0; k < 15; k++) {
            u[iy] = hb1_coeff[stages];
            stages++;
            iy++;
        }

        l_freqz_cg(u, w, Fs, combinedResponse);
        stages = 1;
        break;

    case 2:
        /*  Hb2 */
        for (stages = 0; stages < 7; stages++) {
            h_u[stages] = 0.0;
        }

        stages = 0;
        iy = 0;
        for (k = 0; k < 7; k++) {
            h_u[iy] = hb2_coeff[stages];
            stages++;
            iy++;
        }

        m_freqz_cg(h_u, w, Fs, combinedResponse);
        stages = 1;
        break;

    case 3:
        /*  Hb3 */
        c_us(hb3_coeff_data, hb3_coeff_size, tmp_data, tmp_size);
        b_tmp_data.data = (double *)&tmp_data;
        b_tmp_data.size = (int *)&tmp_size;
        b_tmp_data.allocatedSize = 29;
        b_tmp_data.numDimensions = 2;
        b_tmp_data.canFreeData = false;
        n_freqz_cg(&b_tmp_data, w, Fs, combinedResponse);
        stages = 1;
        break;

    case 4:
        /*  Hb2,Hb1 */
        memset(&b_u[0], 0, 30U * sizeof(double));
        stages = 0;
        iy = 0;
        for (k = 0; k < 15; k++) {
            b_u[iy] = hb1_coeff[stages];
            stages++;
            iy += 2;
        }

        o_freqz_cg(*(double (*)[29])&b_u[0], w, Fs, combinedResponse);
        for (stages = 0; stages < 7; stages++) {
            h_u[stages] = 0.0;
        }

        stages = 0;
        iy = 0;
        for (k = 0; k < 7; k++) {
            h_u[iy] = hb2_coeff[stages];
            stages++;
            iy++;
        }

        m_freqz_cg(h_u, w, Fs, d2);
        stages = combinedResponse->size[0] * combinedResponse->size[1];
        combinedResponse->size[0] = 1;
        emxEnsureCapacity_creal_T(combinedResponse, stages);
        iy = combinedResponse->size[0];
        stages = combinedResponse->size[1];
        iy *= stages;
        for (stages = 0; stages < iy; stages++) {
            combinedResponse_re = combinedResponse->data[stages].re;
            combinedResponse_im = combinedResponse->data[stages].im;
            d2_re = d2->data[stages].re;
            d2_im = d2->data[stages].im;
            combinedResponse->data[stages].re = combinedResponse_re * d2_re -
                                                combinedResponse_im * d2_im;
            combinedResponse->data[stages].im = combinedResponse_re * d2_im +
                                                combinedResponse_im * d2_re;
        }

        stages = 2;
        break;

    case 5:
        /*  Hb3,Hb1 */
        memset(&b_u[0], 0, 30U * sizeof(double));
        stages = 0;
        iy = 0;
        for (k = 0; k < 15; k++) {
            b_u[iy] = hb1_coeff[stages];
            stages++;
            iy += 2;
        }

        o_freqz_cg(*(double (*)[29])&b_u[0], w, Fs, combinedResponse);
        c_us(hb3_coeff_data, hb3_coeff_size, tmp_data, tmp_size);
        b_tmp_data.data = (double *)&tmp_data;
        b_tmp_data.size = (int *)&tmp_size;
        b_tmp_data.allocatedSize = 29;
        b_tmp_data.numDimensions = 2;
        b_tmp_data.canFreeData = false;
        n_freqz_cg(&b_tmp_data, w, Fs, d2);
        stages = combinedResponse->size[0] * combinedResponse->size[1];
        combinedResponse->size[0] = 1;
        emxEnsureCapacity_creal_T(combinedResponse, stages);
        iy = combinedResponse->size[0];
        stages = combinedResponse->size[1];
        iy *= stages;
        for (stages = 0; stages < iy; stages++) {
            combinedResponse_re = combinedResponse->data[stages].re;
            combinedResponse_im = combinedResponse->data[stages].im;
            d2_re = d2->data[stages].re;
            d2_im = d2->data[stages].im;
            combinedResponse->data[stages].re = combinedResponse_re * d2_re -
                                                combinedResponse_im * d2_im;
            combinedResponse->data[stages].im = combinedResponse_re * d2_im +
                                                combinedResponse_im * d2_re;
        }

        stages = 2;
        break;

    case 6:
        /*  Hb3,Hb2 */
        memset(&c_u[0], 0, 14U * sizeof(double));
        stages = 0;
        iy = 0;
        for (k = 0; k < 7; k++) {
            c_u[iy] = hb2_coeff[stages];
            stages++;
            iy += 2;
        }

        p_freqz_cg(*(double (*)[13])&c_u[0], w, Fs, combinedResponse);
        c_us(hb3_coeff_data, hb3_coeff_size, tmp_data, tmp_size);
        b_tmp_data.data = (double *)&tmp_data;
        b_tmp_data.size = (int *)&tmp_size;
        b_tmp_data.allocatedSize = 29;
        b_tmp_data.numDimensions = 2;
        b_tmp_data.canFreeData = false;
        n_freqz_cg(&b_tmp_data, w, Fs, d2);
        stages = combinedResponse->size[0] * combinedResponse->size[1];
        combinedResponse->size[0] = 1;
        emxEnsureCapacity_creal_T(combinedResponse, stages);
        iy = combinedResponse->size[0];
        stages = combinedResponse->size[1];
        iy *= stages;
        for (stages = 0; stages < iy; stages++) {
            combinedResponse_re = combinedResponse->data[stages].re;
            combinedResponse_im = combinedResponse->data[stages].im;
            d2_re = d2->data[stages].re;
            d2_im = d2->data[stages].im;
            combinedResponse->data[stages].re = combinedResponse_re * d2_re -
                                                combinedResponse_im * d2_im;
            combinedResponse->data[stages].im = combinedResponse_re * d2_im +
                                                combinedResponse_im * d2_re;
        }

        stages = 2;
        break;

    case 7:
        /*  Hb3,Hb2,Hb1 */
        memset(&d_u[0], 0, 60U * sizeof(double));
        stages = 0;
        iy = 0;
        for (k = 0; k < 15; k++) {
            d_u[iy] = hb1_coeff[stages];
            stages++;
            iy += 4;
        }

        q_freqz_cg(*(double (*)[57])&d_u[0], w, Fs, combinedResponse);
        memset(&c_u[0], 0, 14U * sizeof(double));
        stages = 0;
        iy = 0;
        for (k = 0; k < 7; k++) {
            c_u[iy] = hb2_coeff[stages];
            stages++;
            iy += 2;
        }

        p_freqz_cg(*(double (*)[13])&c_u[0], w, Fs, d2);
        c_us(hb3_coeff_data, hb3_coeff_size, tmp_data, tmp_size);
        b_tmp_data.data = (double *)&tmp_data;
        b_tmp_data.size = (int *)&tmp_size;
        b_tmp_data.allocatedSize = 29;
        b_tmp_data.numDimensions = 2;
        b_tmp_data.canFreeData = false;
        n_freqz_cg(&b_tmp_data, w, Fs, d3);
        stages = combinedResponse->size[0] * combinedResponse->size[1];
        combinedResponse->size[0] = 1;
        emxEnsureCapacity_creal_T(combinedResponse, stages);
        iy = combinedResponse->size[0];
        stages = combinedResponse->size[1];
        iy *= stages;
        for (stages = 0; stages < iy; stages++) {
            combinedResponse_re = combinedResponse->data[stages].re * d2->data[stages]
                                  .re - combinedResponse->data[stages].im * d2->data[stages].im;
            combinedResponse_im = combinedResponse->data[stages].re * d2->data[stages]
                                  .im + combinedResponse->data[stages].im * d2->data[stages].re;
            d2_re = d3->data[stages].re;
            d2_im = d3->data[stages].im;
            combinedResponse->data[stages].re = combinedResponse_re * d2_re -
                                                combinedResponse_im * d2_im;
            combinedResponse->data[stages].im = combinedResponse_re * d2_im +
                                                combinedResponse_im * d2_re;
        }

        stages = 3;
        break;

    case 8:
        /*  Dec/Int3 */
        c_us(dec_int3_coeff_data, dec_int3_coeff_size, tmp_data, tmp_size);
        b_tmp_data.data = (double *)&tmp_data;
        b_tmp_data.size = (int *)&tmp_size;
        b_tmp_data.allocatedSize = 29;
        b_tmp_data.numDimensions = 2;
        b_tmp_data.canFreeData = false;
        n_freqz_cg(&b_tmp_data, w, Fs, combinedResponse);
        stages = 1;

        /*  RECHECK ALL DEC BY 3     */
        break;

    case 9:
        /*  Dec/Int3,Hb1 */
        memset(&e_u[0], 0, 45U * sizeof(double));
        stages = 0;
        iy = 0;
        for (k = 0; k < 15; k++) {
            e_u[iy] = hb1_coeff[stages];
            stages++;
            iy += 3;
        }

        r_freqz_cg(*(double (*)[43])&e_u[0], w, Fs, combinedResponse);
        c_us(dec_int3_coeff_data, dec_int3_coeff_size, tmp_data, tmp_size);
        b_tmp_data.data = (double *)&tmp_data;
        b_tmp_data.size = (int *)&tmp_size;
        b_tmp_data.allocatedSize = 29;
        b_tmp_data.numDimensions = 2;
        b_tmp_data.canFreeData = false;
        n_freqz_cg(&b_tmp_data, w, Fs, d2);
        stages = combinedResponse->size[0] * combinedResponse->size[1];
        combinedResponse->size[0] = 1;
        emxEnsureCapacity_creal_T(combinedResponse, stages);
        iy = combinedResponse->size[0];
        stages = combinedResponse->size[1];
        iy *= stages;
        for (stages = 0; stages < iy; stages++) {
            combinedResponse_re = combinedResponse->data[stages].re;
            combinedResponse_im = combinedResponse->data[stages].im;
            d2_re = d2->data[stages].re;
            d2_im = d2->data[stages].im;
            combinedResponse->data[stages].re = combinedResponse_re * d2_re -
                                                combinedResponse_im * d2_im;
            combinedResponse->data[stages].im = combinedResponse_re * d2_im +
                                                combinedResponse_im * d2_re;
        }

        stages = 2;
        break;

    case 10:
        /*  Dec/Int3,Hb2 */
        memset(&f_u[0], 0, 21U * sizeof(double));
        stages = 0;
        iy = 0;
        for (k = 0; k < 7; k++) {
            f_u[iy] = hb2_coeff[stages];
            stages++;
            iy += 3;
        }

        s_freqz_cg(*(double (*)[19])&f_u[0], w, Fs, combinedResponse);
        c_us(dec_int3_coeff_data, dec_int3_coeff_size, tmp_data, tmp_size);
        b_tmp_data.data = (double *)&tmp_data;
        b_tmp_data.size = (int *)&tmp_size;
        b_tmp_data.allocatedSize = 29;
        b_tmp_data.numDimensions = 2;
        b_tmp_data.canFreeData = false;
        n_freqz_cg(&b_tmp_data, w, Fs, d2);
        stages = combinedResponse->size[0] * combinedResponse->size[1];
        combinedResponse->size[0] = 1;
        emxEnsureCapacity_creal_T(combinedResponse, stages);
        iy = combinedResponse->size[0];
        stages = combinedResponse->size[1];
        iy *= stages;
        for (stages = 0; stages < iy; stages++) {
            combinedResponse_re = combinedResponse->data[stages].re;
            combinedResponse_im = combinedResponse->data[stages].im;
            d2_re = d2->data[stages].re;
            d2_im = d2->data[stages].im;
            combinedResponse->data[stages].re = combinedResponse_re * d2_re -
                                                combinedResponse_im * d2_im;
            combinedResponse->data[stages].im = combinedResponse_re * d2_im +
                                                combinedResponse_im * d2_re;
        }

        stages = 2;
        break;

    case 11:
        /*  Dec/Int3,Hb2,Hb1 {Hm4,Hm2c34,Hm1} */
        memset(&g_u[0], 0, 90U * sizeof(double));
        stages = 0;
        iy = 0;
        for (k = 0; k < 15; k++) {
            g_u[iy] = hb1_coeff[stages];
            stages++;
            iy += 6;
        }

        t_freqz_cg(*(double (*)[85])&g_u[0], w, Fs, combinedResponse);
        memset(&f_u[0], 0, 21U * sizeof(double));
        stages = 0;
        iy = 0;
        for (k = 0; k < 7; k++) {
            f_u[iy] = hb2_coeff[stages];
            stages++;
            iy += 3;
        }

        s_freqz_cg(*(double (*)[19])&f_u[0], w, Fs, d2);
        c_us(dec_int3_coeff_data, dec_int3_coeff_size, tmp_data, tmp_size);
        b_tmp_data.data = (double *)&tmp_data;
        b_tmp_data.size = (int *)&tmp_size;
        b_tmp_data.allocatedSize = 29;
        b_tmp_data.numDimensions = 2;
        b_tmp_data.canFreeData = false;
        n_freqz_cg(&b_tmp_data, w, Fs, d3);
        stages = combinedResponse->size[0] * combinedResponse->size[1];
        combinedResponse->size[0] = 1;
        emxEnsureCapacity_creal_T(combinedResponse, stages);
        iy = combinedResponse->size[0];
        stages = combinedResponse->size[1];
        iy *= stages;
        for (stages = 0; stages < iy; stages++) {
            combinedResponse_re = combinedResponse->data[stages].re * d2->data[stages]
                                  .re - combinedResponse->data[stages].im * d2->data[stages].im;
            combinedResponse_im = combinedResponse->data[stages].re * d2->data[stages]
                                  .im + combinedResponse->data[stages].im * d2->data[stages].re;
            d2_re = d3->data[stages].re;
            d2_im = d3->data[stages].im;
            combinedResponse->data[stages].re = combinedResponse_re * d2_re -
                                                combinedResponse_im * d2_im;
            combinedResponse->data[stages].im = combinedResponse_re * d2_im +
                                                combinedResponse_im * d2_re;
        }

        stages = 3;
        break;
    }

    emxFree_creal_T(&d3);

    /*  Add filter extra filter to end of cascade */
    if (!(extraTaps->size[1] == 0)) {
        emxInit_real_T(&i_u, 2);
        c = (int)rt_powd_snf(2.0, stages);
        iy = c * extraTaps->size[1];
        stages = i_u->size[0] * i_u->size[1];
        i_u->size[0] = 1;
        i_u->size[1] = iy;
        emxEnsureCapacity_real_T(i_u, stages);
        for (stages = 0; stages < iy; stages++) {
            i_u->data[stages] = 0.0;
        }

        stages = 0;
        iy = 0;
        for (k = 1; k <= extraTaps->size[1]; k++) {
            i_u->data[iy] = extraTaps->data[stages];
            stages++;
            iy += c;
        }

        stages = (i_u->size[1] - c) + 1;
        iy = i_u->size[0] * i_u->size[1];
        if (1 > stages) {
            i_u->size[1] = 0;
        } else {
            i_u->size[1] = stages;
        }

        emxEnsureCapacity_real_T(i_u, iy);
        n_freqz_cg(i_u, w, Fs, d2);
        stages = combinedResponse->size[0] * combinedResponse->size[1];
        combinedResponse->size[0] = 1;
        emxEnsureCapacity_creal_T(combinedResponse, stages);
        iy = combinedResponse->size[0];
        stages = combinedResponse->size[1];
        iy *= stages;
        emxFree_real_T(&i_u);
        for (stages = 0; stages < iy; stages++) {
            combinedResponse_re = combinedResponse->data[stages].re;
            combinedResponse_im = combinedResponse->data[stages].im;
            d2_re = d2->data[stages].re;
            d2_im = d2->data[stages].im;
            combinedResponse->data[stages].re = combinedResponse_re * d2_re -
                                                combinedResponse_im * d2_im;
            combinedResponse->data[stages].im = combinedResponse_re * d2_im +
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
    double x_im;
    if (!(p_size[1] == 0)) {
        for (i13 = 0; i13 < 2048; i13++) {
            y[i13].re = p_data[0];
            y[i13].im = 0.0;
        }

        for (k = 0; k <= p_size[1] - 2; k++) {
            for (i13 = 0; i13 < 2048; i13++) {
                x_im = x[i13].re * y[i13].im + x[i13].im * y[i13].re;
                y[i13].re = (x[i13].re * y[i13].re - x[i13].im * y[i13].im) + p_data[k +
                            1];
                y[i13].im = x_im;
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
    int unnamed_idx_1;
    int k;
    unnamed_idx_1 = a->size[1];
    k = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = a->size[1];
    emxEnsureCapacity_real_T(y, k);
    for (k = 0; k + 1 <= unnamed_idx_1; k++) {
        y->data[k] = rt_powd_snf(a->data[k], 3.0);
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
    int i50;
    int loop_ub;
    i50 = z->size[0] * z->size[1];
    z->size[0] = 1;
    z->size[1] = x->size[1];
    emxEnsureCapacity_real_T(z, i50);
    loop_ub = x->size[0] * x->size[1];
    for (i50 = 0; i50 < loop_ub; i50++) {
        z->data[i50] = x->data[i50] / y->data[i50];
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
    double absxi;
    double absxr;
    xr = x->re;
    xi = x->im;
    if (xi == 0.0) {
        if (xr < 0.0) {
            absxi = 0.0;
            xr = sqrt(-xr);
        } else {
            absxi = sqrt(xr);
            xr = 0.0;
        }
    } else if (xr == 0.0) {
        if (xi < 0.0) {
            absxi = sqrt(-xi / 2.0);
            xr = -absxi;
        } else {
            absxi = sqrt(xi / 2.0);
            xr = absxi;
        }
    } else if (rtIsNaN(xr)) {
        absxi = xr;
    } else if (rtIsNaN(xi)) {
        absxi = xi;
        xr = xi;
    } else if (rtIsInf(xi)) {
        absxi = fabs(xi);
        xr = xi;
    } else if (rtIsInf(xr)) {
        if (xr < 0.0) {
            absxi = 0.0;
            xr = xi * -xr;
        } else {
            absxi = xr;
            xr = 0.0;
        }
    } else {
        absxr = fabs(xr);
        absxi = fabs(xi);
        if ((absxr > 4.4942328371557893E+307) || (absxi > 4.4942328371557893E+307)) {
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

        if (xr > 0.0) {
            xr = 0.5 * (xi / absxi);
        } else {
            if (xi < 0.0) {
                xr = -absxi;
            } else {
                xr = absxi;
            }

            absxi = 0.5 * (xi / xr);
        }
    }

    x->re = absxi;
    x->im = xr;
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
    static const char cv32[2] = { 'T', 'x' };

    b_bool = false;
    kstr = 0;
    do {
        exitg1 = 0;
        if (kstr + 1 < 3) {
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
    int ix;
    int iy;
    int k;
    u_size[0] = 1;
    u_size[1] = (signed char)o_size[1];
    ix = (signed char)o_size[1];
    if (0 <= ix - 1) {
        memset(&u_data[0], 0, (unsigned int)(ix * (int)sizeof(double)));
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
    emxArray_real_T *r28;
    emxArray_real_T *r29;
    int i56;
    int loop_ub;
    emxArray_creal_T *r30;
    double re;
    double im;
    double b_re;
    double b_im;
    emxInit_real_T(&r28, 2);
    emxInit_real_T(&r29, 2);
    rdivide(f, Fconverter, r28);
    b_sinc(r28);
    i56 = r29->size[0] * r29->size[1];
    r29->size[0] = 1;
    r29->size[1] = f->size[1];
    emxEnsureCapacity_real_T(r29, i56);
    loop_ub = f->size[0] * f->size[1];
    for (i56 = 0; i56 < loop_ub; i56++) {
        r29->data[i56] = 6.2831853071795862 * f->data[i56];
    }

    b_freqs_cg(b1_data, b1_size, a1, r29, abc);
    i56 = r29->size[0] * r29->size[1];
    r29->size[0] = 1;
    r29->size[1] = f->size[1];
    emxEnsureCapacity_real_T(r29, i56);
    loop_ub = f->size[0] * f->size[1];
    for (i56 = 0; i56 < loop_ub; i56++) {
        r29->data[i56] = 6.2831853071795862 * f->data[i56];
    }

    emxInit_creal_T(&r30, 2);
    b_freqs_cg(b2_data, b2_size, a2, r29, r30);
    i56 = abc->size[0] * abc->size[1];
    abc->size[0] = 1;
    abc->size[1] = r28->size[1];
    emxEnsureCapacity_creal_T(abc, i56);
    loop_ub = r28->size[0] * r28->size[1];
    emxFree_real_T(&r29);
    for (i56 = 0; i56 < loop_ub; i56++) {
        re = r28->data[i56] * abc->data[i56].re;
        im = r28->data[i56] * abc->data[i56].im;
        b_re = r30->data[i56].re;
        b_im = r30->data[i56].im;
        abc->data[i56].re = re * b_re - im * b_im;
        abc->data[i56].im = re * b_im + im * b_re;
    }

    emxFree_real_T(&r28);
    emxFree_creal_T(&r30);
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
    int b_idx_0;
    int loop_ub;
    double b_b_data[29];
    creal_T dcv4[2048];
    double digw;
    double b_digw[2048];
    double h_re;
    double brm;
    double d;

    /* -------------------------------------------------------------------------- */
    b_idx_0 = b_size[1];
    loop_ub = b_size[1];
    if (0 <= loop_ub - 1) {
        memcpy(&b_b_data[0], &b_data[0], (unsigned int)(loop_ub * (int)sizeof(double)));
    }

    b_size[0] = 1;
    if (0 <= b_idx_0 - 1) {
        memcpy(&b_data[0], &b_b_data[0], (unsigned int)(b_idx_0 * (int)sizeof(double)));
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
    for (b_idx_0 = 0; b_idx_0 < 2048; b_idx_0++) {
        w[b_idx_0] = options->w[b_idx_0];
        digw = 6.2831853071795862 * options->w[b_idx_0] / options->Fs;
        dcv4[b_idx_0].re = digw * 0.0;
        dcv4[b_idx_0].im = digw;
        b_digw[b_idx_0] = digw;
    }

    b_exp(dcv4);
    c_polyval(b_data, b_size, dcv4, h);
    for (b_idx_0 = 0; b_idx_0 < 2048; b_idx_0++) {
        dcv4[b_idx_0].re = ((double)b_size[1] - 1.0) * (b_digw[b_idx_0] * 0.0);
        dcv4[b_idx_0].im = ((double)b_size[1] - 1.0) * b_digw[b_idx_0];
    }

    b_exp(dcv4);
    for (b_idx_0 = 0; b_idx_0 < 2048; b_idx_0++) {
        h_re = h[b_idx_0].re;
        if (dcv4[b_idx_0].im == 0.0) {
            if (h[b_idx_0].im == 0.0) {
                h[b_idx_0].re /= dcv4[b_idx_0].re;
                h[b_idx_0].im = 0.0;
            } else if (h[b_idx_0].re == 0.0) {
                h[b_idx_0].re = 0.0;
                h[b_idx_0].im /= dcv4[b_idx_0].re;
            } else {
                h[b_idx_0].re /= dcv4[b_idx_0].re;
                h[b_idx_0].im /= dcv4[b_idx_0].re;
            }
        } else if (dcv4[b_idx_0].re == 0.0) {
            if (h[b_idx_0].re == 0.0) {
                h[b_idx_0].re = h[b_idx_0].im / dcv4[b_idx_0].im;
                h[b_idx_0].im = 0.0;
            } else if (h[b_idx_0].im == 0.0) {
                h[b_idx_0].re = 0.0;
                h[b_idx_0].im = -(h_re / dcv4[b_idx_0].im);
            } else {
                h[b_idx_0].re = h[b_idx_0].im / dcv4[b_idx_0].im;
                h[b_idx_0].im = -(h_re / dcv4[b_idx_0].im);
            }
        } else {
            brm = fabs(dcv4[b_idx_0].re);
            digw = fabs(dcv4[b_idx_0].im);
            if (brm > digw) {
                digw = dcv4[b_idx_0].im / dcv4[b_idx_0].re;
                d = dcv4[b_idx_0].re + digw * dcv4[b_idx_0].im;
                h[b_idx_0].re = (h[b_idx_0].re + digw * h[b_idx_0].im) / d;
                h[b_idx_0].im = (h[b_idx_0].im - digw * h_re) / d;
            } else if (digw == brm) {
                if (dcv4[b_idx_0].re > 0.0) {
                    digw = 0.5;
                } else {
                    digw = -0.5;
                }

                if (dcv4[b_idx_0].im > 0.0) {
                    d = 0.5;
                } else {
                    d = -0.5;
                }

                h[b_idx_0].re = (h[b_idx_0].re * digw + h[b_idx_0].im * d) / brm;
                h[b_idx_0].im = (h[b_idx_0].im * digw - h_re * d) / brm;
            } else {
                digw = dcv4[b_idx_0].re / dcv4[b_idx_0].im;
                d = dcv4[b_idx_0].im + digw * dcv4[b_idx_0].re;
                h[b_idx_0].re = (digw * h[b_idx_0].re + h[b_idx_0].im) / d;
                h[b_idx_0].im = (digw * h[b_idx_0].im - h_re) / d;
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
static void d_polyval(const double p[29], const creal_T x[2048],
                      creal_T y[2048])
{
    int i14;
    int k;
    double x_im;
    for (i14 = 0; i14 < 2048; i14++) {
        y[i14].re = p[0];
        y[i14].im = 0.0;
    }

    for (k = 0; k < 28; k++) {
        for (i14 = 0; i14 < 2048; i14++) {
            x_im = x[i14].re * y[i14].im + x[i14].im * y[i14].re;
            y[i14].re = (x[i14].re * y[i14].re - x[i14].im * y[i14].im) + p[k + 1];
            y[i14].im = x_im;
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
    int ixstart;
    emxArray_real_T *org;
    int ix;
    emxArray_real_T *r;
    double u;
    double v;
    emxArray_real_T *b_r;
    short i55;
    emxArray_real_T *r24;
    double e[16];
    int itmp;
    boolean_T exitg1;
    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    emxInit_real_T(&org, 2);
    ix = org->size[0] * org->size[1];
    org->size[0] = 1;
    org->size[1] = ixstart;
    emxEnsureCapacity_real_T(org, ix);
    for (ix = 0; ix < ixstart; ix++) {
        org->data[org->size[0] * ix] = tap_store->data[((int)i + tap_store->size[0] *
                                       ix) - 1];
    }

    emxInit_real_T(&r, 2);
    ix = r->size[0] * r->size[1];
    r->size[0] = 16;
    r->size[1] = (int)M;
    emxEnsureCapacity_real_T(r, ix);
    ixstart = (int)M << 4;
    for (ix = 0; ix < ixstart; ix++) {
        r->data[ix] = 0.0;
    }

    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        u = tap_store->data[((int)i + tap_store->size[0] * ix) - 1] * 2.0;
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
                i55 = (short)u;
            } else {
                i55 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i55 = MAX_int16_T;
        } else {
            i55 = 0;
        }

        r->data[r->size[0] * ix] = (double)i55 * 0.5;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    emxInit_real_T(&b_r, 2);
    ix = b_r->size[0] * b_r->size[1];
    b_r->size[0] = 1;
    b_r->size[1] = ixstart;
    emxEnsureCapacity_real_T(b_r, ix);
    for (ix = 0; ix < ixstart; ix++) {
        b_r->data[b_r->size[0] * ix] = r->data[r->size[0] * ix] - org->data
                                       [org->size[0] * ix];
    }

    emxInit_real_T(&r24, 2);
    c_abs(b_r, r24);
    e[0] = b_sum(r24);
    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        u = tap_store->data[((int)i + tap_store->size[0] * ix) - 1] * 4.0;
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
                i55 = (short)u;
            } else {
                i55 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i55 = MAX_int16_T;
        } else {
            i55 = 0;
        }

        r->data[1 + r->size[0] * ix] = (double)i55 * 0.25;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    ix = b_r->size[0] * b_r->size[1];
    b_r->size[0] = 1;
    b_r->size[1] = ixstart;
    emxEnsureCapacity_real_T(b_r, ix);
    for (ix = 0; ix < ixstart; ix++) {
        b_r->data[b_r->size[0] * ix] = r->data[1 + r->size[0] * ix] - org->data
                                       [org->size[0] * ix];
    }

    c_abs(b_r, r24);
    e[1] = b_sum(r24);
    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        u = tap_store->data[((int)i + tap_store->size[0] * ix) - 1] * 8.0;
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
                i55 = (short)u;
            } else {
                i55 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i55 = MAX_int16_T;
        } else {
            i55 = 0;
        }

        r->data[2 + r->size[0] * ix] = (double)i55 * 0.125;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    ix = b_r->size[0] * b_r->size[1];
    b_r->size[0] = 1;
    b_r->size[1] = ixstart;
    emxEnsureCapacity_real_T(b_r, ix);
    for (ix = 0; ix < ixstart; ix++) {
        b_r->data[b_r->size[0] * ix] = r->data[2 + r->size[0] * ix] - org->data
                                       [org->size[0] * ix];
    }

    c_abs(b_r, r24);
    e[2] = b_sum(r24);
    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        u = tap_store->data[((int)i + tap_store->size[0] * ix) - 1] * 16.0;
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
                i55 = (short)u;
            } else {
                i55 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i55 = MAX_int16_T;
        } else {
            i55 = 0;
        }

        r->data[3 + r->size[0] * ix] = (double)i55 * 0.0625;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    ix = b_r->size[0] * b_r->size[1];
    b_r->size[0] = 1;
    b_r->size[1] = ixstart;
    emxEnsureCapacity_real_T(b_r, ix);
    for (ix = 0; ix < ixstart; ix++) {
        b_r->data[b_r->size[0] * ix] = r->data[3 + r->size[0] * ix] - org->data
                                       [org->size[0] * ix];
    }

    c_abs(b_r, r24);
    e[3] = b_sum(r24);
    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        u = tap_store->data[((int)i + tap_store->size[0] * ix) - 1] * 32.0;
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
                i55 = (short)u;
            } else {
                i55 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i55 = MAX_int16_T;
        } else {
            i55 = 0;
        }

        r->data[4 + r->size[0] * ix] = (double)i55 * 0.03125;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    ix = b_r->size[0] * b_r->size[1];
    b_r->size[0] = 1;
    b_r->size[1] = ixstart;
    emxEnsureCapacity_real_T(b_r, ix);
    for (ix = 0; ix < ixstart; ix++) {
        b_r->data[b_r->size[0] * ix] = r->data[4 + r->size[0] * ix] - org->data
                                       [org->size[0] * ix];
    }

    c_abs(b_r, r24);
    e[4] = b_sum(r24);
    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        u = tap_store->data[((int)i + tap_store->size[0] * ix) - 1] * 64.0;
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
                i55 = (short)u;
            } else {
                i55 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i55 = MAX_int16_T;
        } else {
            i55 = 0;
        }

        r->data[5 + r->size[0] * ix] = (double)i55 * 0.015625;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    ix = b_r->size[0] * b_r->size[1];
    b_r->size[0] = 1;
    b_r->size[1] = ixstart;
    emxEnsureCapacity_real_T(b_r, ix);
    for (ix = 0; ix < ixstart; ix++) {
        b_r->data[b_r->size[0] * ix] = r->data[5 + r->size[0] * ix] - org->data
                                       [org->size[0] * ix];
    }

    c_abs(b_r, r24);
    e[5] = b_sum(r24);
    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        u = tap_store->data[((int)i + tap_store->size[0] * ix) - 1] * 128.0;
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
                i55 = (short)u;
            } else {
                i55 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i55 = MAX_int16_T;
        } else {
            i55 = 0;
        }

        r->data[6 + r->size[0] * ix] = (double)i55 * 0.0078125;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    ix = b_r->size[0] * b_r->size[1];
    b_r->size[0] = 1;
    b_r->size[1] = ixstart;
    emxEnsureCapacity_real_T(b_r, ix);
    for (ix = 0; ix < ixstart; ix++) {
        b_r->data[b_r->size[0] * ix] = r->data[6 + r->size[0] * ix] - org->data
                                       [org->size[0] * ix];
    }

    c_abs(b_r, r24);
    e[6] = b_sum(r24);
    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        u = tap_store->data[((int)i + tap_store->size[0] * ix) - 1] * 256.0;
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
                i55 = (short)u;
            } else {
                i55 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i55 = MAX_int16_T;
        } else {
            i55 = 0;
        }

        r->data[7 + r->size[0] * ix] = (double)i55 * 0.00390625;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    ix = b_r->size[0] * b_r->size[1];
    b_r->size[0] = 1;
    b_r->size[1] = ixstart;
    emxEnsureCapacity_real_T(b_r, ix);
    for (ix = 0; ix < ixstart; ix++) {
        b_r->data[b_r->size[0] * ix] = r->data[7 + r->size[0] * ix] - org->data
                                       [org->size[0] * ix];
    }

    c_abs(b_r, r24);
    e[7] = b_sum(r24);
    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        u = tap_store->data[((int)i + tap_store->size[0] * ix) - 1] * 512.0;
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
                i55 = (short)u;
            } else {
                i55 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i55 = MAX_int16_T;
        } else {
            i55 = 0;
        }

        r->data[8 + r->size[0] * ix] = (double)i55 * 0.001953125;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    ix = b_r->size[0] * b_r->size[1];
    b_r->size[0] = 1;
    b_r->size[1] = ixstart;
    emxEnsureCapacity_real_T(b_r, ix);
    for (ix = 0; ix < ixstart; ix++) {
        b_r->data[b_r->size[0] * ix] = r->data[8 + r->size[0] * ix] - org->data
                                       [org->size[0] * ix];
    }

    c_abs(b_r, r24);
    e[8] = b_sum(r24);
    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        u = tap_store->data[((int)i + tap_store->size[0] * ix) - 1] * 1024.0;
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
                i55 = (short)u;
            } else {
                i55 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i55 = MAX_int16_T;
        } else {
            i55 = 0;
        }

        r->data[9 + r->size[0] * ix] = (double)i55 * 0.0009765625;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    ix = b_r->size[0] * b_r->size[1];
    b_r->size[0] = 1;
    b_r->size[1] = ixstart;
    emxEnsureCapacity_real_T(b_r, ix);
    for (ix = 0; ix < ixstart; ix++) {
        b_r->data[b_r->size[0] * ix] = r->data[9 + r->size[0] * ix] - org->data
                                       [org->size[0] * ix];
    }

    c_abs(b_r, r24);
    e[9] = b_sum(r24);
    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        u = tap_store->data[((int)i + tap_store->size[0] * ix) - 1] * 2048.0;
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
                i55 = (short)u;
            } else {
                i55 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i55 = MAX_int16_T;
        } else {
            i55 = 0;
        }

        r->data[10 + r->size[0] * ix] = (double)i55 * 0.00048828125;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    ix = b_r->size[0] * b_r->size[1];
    b_r->size[0] = 1;
    b_r->size[1] = ixstart;
    emxEnsureCapacity_real_T(b_r, ix);
    for (ix = 0; ix < ixstart; ix++) {
        b_r->data[b_r->size[0] * ix] = r->data[10 + r->size[0] * ix] - org->data
                                       [org->size[0] * ix];
    }

    c_abs(b_r, r24);
    e[10] = b_sum(r24);
    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        u = tap_store->data[((int)i + tap_store->size[0] * ix) - 1] * 4096.0;
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
                i55 = (short)u;
            } else {
                i55 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i55 = MAX_int16_T;
        } else {
            i55 = 0;
        }

        r->data[11 + r->size[0] * ix] = (double)i55 * 0.000244140625;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    ix = b_r->size[0] * b_r->size[1];
    b_r->size[0] = 1;
    b_r->size[1] = ixstart;
    emxEnsureCapacity_real_T(b_r, ix);
    for (ix = 0; ix < ixstart; ix++) {
        b_r->data[b_r->size[0] * ix] = r->data[11 + r->size[0] * ix] - org->data
                                       [org->size[0] * ix];
    }

    c_abs(b_r, r24);
    e[11] = b_sum(r24);
    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        u = tap_store->data[((int)i + tap_store->size[0] * ix) - 1] * 8192.0;
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
                i55 = (short)u;
            } else {
                i55 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i55 = MAX_int16_T;
        } else {
            i55 = 0;
        }

        r->data[12 + r->size[0] * ix] = (double)i55 * 0.0001220703125;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    ix = b_r->size[0] * b_r->size[1];
    b_r->size[0] = 1;
    b_r->size[1] = ixstart;
    emxEnsureCapacity_real_T(b_r, ix);
    for (ix = 0; ix < ixstart; ix++) {
        b_r->data[b_r->size[0] * ix] = r->data[12 + r->size[0] * ix] - org->data
                                       [org->size[0] * ix];
    }

    c_abs(b_r, r24);
    e[12] = b_sum(r24);
    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        u = tap_store->data[((int)i + tap_store->size[0] * ix) - 1] * 16384.0;
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
                i55 = (short)u;
            } else {
                i55 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i55 = MAX_int16_T;
        } else {
            i55 = 0;
        }

        r->data[13 + r->size[0] * ix] = (double)i55 * 6.103515625E-5;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    ix = b_r->size[0] * b_r->size[1];
    b_r->size[0] = 1;
    b_r->size[1] = ixstart;
    emxEnsureCapacity_real_T(b_r, ix);
    for (ix = 0; ix < ixstart; ix++) {
        b_r->data[b_r->size[0] * ix] = r->data[13 + r->size[0] * ix] - org->data
                                       [org->size[0] * ix];
    }

    c_abs(b_r, r24);
    e[13] = b_sum(r24);
    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        u = tap_store->data[((int)i + tap_store->size[0] * ix) - 1] * 32768.0;
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
                i55 = (short)u;
            } else {
                i55 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i55 = MAX_int16_T;
        } else {
            i55 = 0;
        }

        r->data[14 + r->size[0] * ix] = (double)i55 * 3.0517578125E-5;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    ix = b_r->size[0] * b_r->size[1];
    b_r->size[0] = 1;
    b_r->size[1] = ixstart;
    emxEnsureCapacity_real_T(b_r, ix);
    for (ix = 0; ix < ixstart; ix++) {
        b_r->data[b_r->size[0] * ix] = r->data[14 + r->size[0] * ix] - org->data
                                       [org->size[0] * ix];
    }

    c_abs(b_r, r24);
    e[14] = b_sum(r24);
    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        u = tap_store->data[((int)i + tap_store->size[0] * ix) - 1] * 65536.0;
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
                i55 = (short)u;
            } else {
                i55 = MIN_int16_T;
            }
        } else if (u >= 32768.0) {
            i55 = MAX_int16_T;
        } else {
            i55 = 0;
        }

        r->data[15 + r->size[0] * ix] = (double)i55 * 1.52587890625E-5;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    ix = b_r->size[0] * b_r->size[1];
    b_r->size[0] = 1;
    b_r->size[1] = ixstart;
    emxEnsureCapacity_real_T(b_r, ix);
    for (ix = 0; ix < ixstart; ix++) {
        b_r->data[b_r->size[0] * ix] = r->data[15 + r->size[0] * ix] - org->data
                                       [org->size[0] * ix];
    }

    emxFree_real_T(&org);
    c_abs(b_r, r24);
    e[15] = b_sum(r24);
    ixstart = 1;
    u = e[0];
    itmp = 0;
    emxFree_real_T(&b_r);
    emxFree_real_T(&r24);
    if (rtIsNaN(e[0])) {
        ix = 2;
        exitg1 = false;
        while ((!exitg1) && (ix < 17)) {
            ixstart = ix;
            if (!rtIsNaN(e[ix - 1])) {
                u = e[ix - 1];
                itmp = ix - 1;
                exitg1 = true;
            } else {
                ix++;
            }
        }
    }

    if (ixstart < 16) {
        while (ixstart + 1 < 17) {
            if (e[ixstart] < u) {
                u = e[ixstart];
                itmp = ixstart;
            }

            ixstart++;
        }
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    ix = taps->size[0] * taps->size[1];
    taps->size[0] = 1;
    taps->size[1] = ixstart;
    emxEnsureCapacity_real_T(taps, ix);
    for (ix = 0; ix < ixstart; ix++) {
        taps->data[taps->size[0] * ix] = r->data[itmp + r->size[0] * ix];
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
 * Arguments    : const double b[29]
 *                const struct_T *options
 *                creal_T h[2048]
 *                double w[2048]
 * Return Type  : void
 */
static void e_firfreqz(const double b[29], const struct_T *options, creal_T h
                       [2048], double w[2048])
{
    int i66;
    creal_T dcv5[2048];
    double digw;
    double b_digw[2048];
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
    for (i66 = 0; i66 < 2048; i66++) {
        w[i66] = options->w[i66];
        digw = 6.2831853071795862 * options->w[i66] / options->Fs;
        dcv5[i66].re = digw * 0.0;
        dcv5[i66].im = digw;
        b_digw[i66] = digw;
    }

    b_exp(dcv5);
    d_polyval(b, dcv5, h);
    for (i66 = 0; i66 < 2048; i66++) {
        dcv5[i66].re = 28.0 * (b_digw[i66] * 0.0);
        dcv5[i66].im = 28.0 * b_digw[i66];
    }

    b_exp(dcv5);
    for (i66 = 0; i66 < 2048; i66++) {
        h_re = h[i66].re;
        if (dcv5[i66].im == 0.0) {
            if (h[i66].im == 0.0) {
                h[i66].re /= dcv5[i66].re;
                h[i66].im = 0.0;
            } else if (h[i66].re == 0.0) {
                h[i66].re = 0.0;
                h[i66].im /= dcv5[i66].re;
            } else {
                h[i66].re /= dcv5[i66].re;
                h[i66].im /= dcv5[i66].re;
            }
        } else if (dcv5[i66].re == 0.0) {
            if (h[i66].re == 0.0) {
                h[i66].re = h[i66].im / dcv5[i66].im;
                h[i66].im = 0.0;
            } else if (h[i66].im == 0.0) {
                h[i66].re = 0.0;
                h[i66].im = -(h_re / dcv5[i66].im);
            } else {
                h[i66].re = h[i66].im / dcv5[i66].im;
                h[i66].im = -(h_re / dcv5[i66].im);
            }
        } else {
            brm = fabs(dcv5[i66].re);
            digw = fabs(dcv5[i66].im);
            if (brm > digw) {
                digw = dcv5[i66].im / dcv5[i66].re;
                d = dcv5[i66].re + digw * dcv5[i66].im;
                h[i66].re = (h[i66].re + digw * h[i66].im) / d;
                h[i66].im = (h[i66].im - digw * h_re) / d;
            } else if (digw == brm) {
                if (dcv5[i66].re > 0.0) {
                    digw = 0.5;
                } else {
                    digw = -0.5;
                }

                if (dcv5[i66].im > 0.0) {
                    d = 0.5;
                } else {
                    d = -0.5;
                }

                h[i66].re = (h[i66].re * digw + h[i66].im * d) / brm;
                h[i66].im = (h[i66].im * digw - h_re * d) / brm;
            } else {
                digw = dcv5[i66].re / dcv5[i66].im;
                d = dcv5[i66].im + digw * dcv5[i66].re;
                h[i66].re = (digw * h[i66].re + h[i66].im) / d;
                h[i66].im = (digw * h[i66].im - h_re) / d;
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
 * Arguments    : const double p[13]
 *                const creal_T x[2048]
 *                creal_T y[2048]
 * Return Type  : void
 */
static void e_polyval(const double p[13], const creal_T x[2048],
                      creal_T y[2048])
{
    int i16;
    int k;
    double x_im;
    for (i16 = 0; i16 < 2048; i16++) {
        y[i16].re = p[0];
        y[i16].im = 0.0;
    }

    for (k = 0; k < 12; k++) {
        for (i16 = 0; i16 < 2048; i16++) {
            x_im = x[i16].re * y[i16].im + x[i16].im * y[i16].re;
            y[i16].re = (x[i16].re * y[i16].re - x[i16].im * y[i16].im) + p[k + 1];
            y[i16].im = x_im;
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
 * Arguments    : const emxArray_creal_T *A
 *                emxArray_creal_T *V
 * Return Type  : void
 */
static void eig(const emxArray_creal_T *A, emxArray_creal_T *V)
{
    int info;
    boolean_T p;
    int i;
    boolean_T exitg2;
    emxArray_creal_T *T;
    unsigned int uv0[2];
    int exitg1;
    double A_re;
    double V_re;
    double A_im;
    double V_im;
    boolean_T b_A;
    double beta1_re;
    double beta1_im;
    int m;
    double brm;
    int istart;
    int jend;
    if ((A->size[0] == 0) || (A->size[1] == 0)) {
        info = V->size[0];
        V->size[0] = A->size[0];
        emxEnsureCapacity_creal_T1(V, info);
        i = A->size[0];
        for (info = 0; info < i; info++) {
            V->data[info].re = 0.0;
            V->data[info].im = 0.0;
        }
    } else if (anyNonFinite(A)) {
        if ((A->size[0] == 1) && (A->size[1] == 1)) {
            info = V->size[0];
            V->size[0] = 1;
            emxEnsureCapacity_creal_T1(V, info);
            V->data[0].re = rtNaN;
            V->data[0].im = 0.0;
        } else {
            info = V->size[0];
            V->size[0] = A->size[0];
            emxEnsureCapacity_creal_T1(V, info);
            i = A->size[0];
            for (info = 0; info < i; info++) {
                V->data[info].re = rtNaN;
                V->data[info].im = 0.0;
            }
        }
    } else if ((A->size[0] == 1) && (A->size[1] == 1)) {
        info = V->size[0];
        V->size[0] = 1;
        emxEnsureCapacity_creal_T1(V, info);
        V->data[0] = A->data[0];
    } else {
        p = (A->size[0] == A->size[1]);
        if (p) {
            info = 0;
            exitg2 = false;
            while ((!exitg2) && (info <= A->size[1] - 1)) {
                i = 0;
                do {
                    exitg1 = 0;
                    if (i <= info) {
                        A_re = A->data[info + A->size[0] * i].re;
                        A_im = -A->data[info + A->size[0] * i].im;
                        b_A = ((A->data[i + A->size[0] * info].re == A_re) && (A->data[i +
                                A->size[0] * info].im == A_im));
                        if (!b_A) {
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
            if (anyNonFinite(A)) {
                for (info = 0; info < 2; info++) {
                    uv0[info] = (unsigned int)A->size[info];
                }

                info = T->size[0] * T->size[1];
                T->size[0] = (int)uv0[0];
                T->size[1] = (int)uv0[1];
                emxEnsureCapacity_creal_T(T, info);
                i = (int)uv0[0] * (int)uv0[1];
                for (info = 0; info < i; info++) {
                    T->data[info].re = rtNaN;
                    T->data[info].im = 0.0;
                }

                m = T->size[0];
                if (!(1 >= T->size[0])) {
                    istart = 2;
                    if (T->size[0] - 2 < T->size[1] - 1) {
                        jend = T->size[0] - 1;
                    } else {
                        jend = T->size[1];
                    }

                    for (info = 1; info <= jend; info++) {
                        for (i = istart; i <= m; i++) {
                            T->data[(i + T->size[0] * (info - 1)) - 1].re = 0.0;
                            T->data[(i + T->size[0] * (info - 1)) - 1].im = 0.0;
                        }

                        istart++;
                    }
                }
            } else {
                info = T->size[0] * T->size[1];
                T->size[0] = A->size[0];
                T->size[1] = A->size[1];
                emxEnsureCapacity_creal_T(T, info);
                i = A->size[0] * A->size[1];
                for (info = 0; info < i; info++) {
                    T->data[info] = A->data[info];
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

                    for (info = 1; info <= jend; info++) {
                        for (i = istart; i <= m; i++) {
                            T->data[(i + T->size[0] * (info - 1)) - 1].re = 0.0;
                            T->data[(i + T->size[0] * (info - 1)) - 1].im = 0.0;
                        }

                        istart++;
                    }
                }
            }

            info = V->size[0];
            V->size[0] = T->size[0];
            emxEnsureCapacity_creal_T1(V, info);
            for (info = 0; info + 1 <= T->size[0]; info++) {
                V->data[info] = T->data[info + T->size[0] * info];
            }

            emxFree_creal_T(&T);
        } else {
            emxInit_creal_T1(&T, 1);
            xzgeev(A, &info, V, T);
            info = V->size[0];
            emxEnsureCapacity_creal_T1(V, info);
            i = V->size[0];
            for (info = 0; info < i; info++) {
                V_re = V->data[info].re;
                V_im = V->data[info].im;
                beta1_re = T->data[info].re;
                beta1_im = T->data[info].im;
                if (beta1_im == 0.0) {
                    if (V_im == 0.0) {
                        V->data[info].re = V_re / beta1_re;
                        V->data[info].im = 0.0;
                    } else if (V_re == 0.0) {
                        V->data[info].re = 0.0;
                        V->data[info].im = V_im / beta1_re;
                    } else {
                        V->data[info].re = V_re / beta1_re;
                        V->data[info].im = V_im / beta1_re;
                    }
                } else if (beta1_re == 0.0) {
                    if (V_re == 0.0) {
                        V->data[info].re = V_im / beta1_im;
                        V->data[info].im = 0.0;
                    } else if (V_im == 0.0) {
                        V->data[info].re = 0.0;
                        V->data[info].im = -(V_re / beta1_im);
                    } else {
                        V->data[info].re = V_im / beta1_im;
                        V->data[info].im = -(V_re / beta1_im);
                    }
                } else {
                    brm = fabs(beta1_re);
                    A_re = fabs(beta1_im);
                    if (brm > A_re) {
                        A_im = beta1_im / beta1_re;
                        A_re = beta1_re + A_im * beta1_im;
                        V->data[info].re = (V_re + A_im * V_im) / A_re;
                        V->data[info].im = (V_im - A_im * V_re) / A_re;
                    } else if (A_re == brm) {
                        if (beta1_re > 0.0) {
                            A_im = 0.5;
                        } else {
                            A_im = -0.5;
                        }

                        if (beta1_im > 0.0) {
                            A_re = 0.5;
                        } else {
                            A_re = -0.5;
                        }

                        V->data[info].re = (V_re * A_im + V_im * A_re) / brm;
                        V->data[info].im = (V_im * A_im - V_re * A_re) / brm;
                    } else {
                        A_im = beta1_re / beta1_im;
                        A_re = beta1_im + A_im * beta1_re;
                        V->data[info].re = (A_im * V_re + V_im) / A_re;
                        V->data[info].im = (A_im * V_im - V_re) / A_re;
                    }
                }
            }

            emxFree_creal_T(&T);
        }
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
    int maxval;
    int itmax;
    int ldh;
    int i;
    double SMLNUM;
    double tst;
    boolean_T exitg1;
    double aa;
    double ba;
    int L;
    creal_T u2;
    boolean_T goto140;
    int its;
    boolean_T exitg2;
    int k;
    boolean_T exitg3;
    double htmp1;
    creal_T y;
    double ab;
    boolean_T goto70;
    int m;
    double x_re;
    double u_re;
    double x_im;
    double u_im;
    double s;
    int b_k;
    creal_T v[2];
    double b_SMLNUM;
    int i61;
    n = h->size[0];
    if (h->size[0] < 10) {
        maxval = 10;
    } else {
        maxval = h->size[0];
    }

    itmax = 30 * maxval;
    ldh = h->size[0];
    info = 0;
    if ((h->size[0] != 0) && (1 != h->size[0])) {
        for (maxval = 0; maxval + 1 <= n - 3; maxval++) {
            h->data[(maxval + h->size[0] * maxval) + 2].re = 0.0;
            h->data[(maxval + h->size[0] * maxval) + 2].im = 0.0;
            h->data[(maxval + h->size[0] * maxval) + 3].re = 0.0;
            h->data[(maxval + h->size[0] * maxval) + 3].im = 0.0;
        }

        if (1 <= n - 2) {
            h->data[(n + h->size[0] * (n - 3)) - 1].re = 0.0;
            h->data[(n + h->size[0] * (n - 3)) - 1].im = 0.0;
        }

        for (i = 1; i + 1 <= n; i++) {
            if (h->data[i + h->size[0] * (i - 1)].im != 0.0) {
                tst = h->data[i + h->size[0] * (i - 1)].re;
                aa = h->data[i + h->size[0] * (i - 1)].im;
                ba = fabs(h->data[i + h->size[0] * (i - 1)].re) + fabs(h->data[i +
                        h->size[0] * (i - 1)].im);
                if (aa == 0.0) {
                    u2.re = tst / ba;
                    u2.im = 0.0;
                } else if (tst == 0.0) {
                    u2.re = 0.0;
                    u2.im = aa / ba;
                } else {
                    u2.re = tst / ba;
                    u2.im = aa / ba;
                }

                ba = rt_hypotd_snf(u2.re, u2.im);
                if (-u2.im == 0.0) {
                    u2.re /= ba;
                    u2.im = 0.0;
                } else if (u2.re == 0.0) {
                    u2.re = 0.0;
                    u2.im = -u2.im / ba;
                } else {
                    u2.re /= ba;
                    u2.im = -u2.im / ba;
                }

                tst = h->data[i + h->size[0] * (i - 1)].re;
                htmp1 = h->data[i + h->size[0] * (i - 1)].im;
                h->data[i + h->size[0] * (i - 1)].re = rt_hypotd_snf(tst, htmp1);
                h->data[i + h->size[0] * (i - 1)].im = 0.0;
                b_xscal(n - i, u2, h, (i + i * ldh) + 1, ldh);
                y.re = u2.re;
                y.im = -u2.im;
                maxval = i + 2;
                if (n < maxval) {
                    maxval = n;
                }

                xscal(maxval, y, h, 1 + i * ldh);
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
            while ((!exitg2) && (its <= itmax)) {
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
                        x_re = h->data[(k + h->size[0] * (k - 1)) - 1].re - h->data[k +
                                h->size[0] * k].re;
                        x_im = h->data[(k + h->size[0] * (k - 1)) - 1].im - h->data[k +
                                h->size[0] * k].im;
                        tst = fabs(x_re) + fabs(x_im);
                        if (htmp1 > tst) {
                            aa = htmp1;
                            htmp1 = tst;
                        } else {
                            aa = tst;
                        }

                        s = aa + ab;
                        tst = 2.2204460492503131E-16 * (htmp1 * (aa / s));
                        if ((SMLNUM > tst) || rtIsNaN(tst)) {
                            b_SMLNUM = SMLNUM;
                        } else {
                            b_SMLNUM = tst;
                        }

                        if (ba * (ab / s) <= b_SMLNUM) {
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
                        ba = 0.75 * fabs(h->data[(k + h->size[0] * k) + 1].re) + h->data[k +
                                h->size[0] * k].re;
                        ab = h->data[k + h->size[0] * k].im;
                    } else if (its == 20) {
                        ba = 0.75 * fabs(h->data[i + h->size[0] * (i - 1)].re) + h->data[i +
                                h->size[0] * i].re;
                        ab = h->data[i + h->size[0] * i].im;
                    } else {
                        ba = h->data[i + h->size[0] * i].re;
                        ab = h->data[i + h->size[0] * i].im;
                        y = h->data[(i + h->size[0] * i) - 1];
                        c_sqrt(&y);
                        u2 = h->data[i + h->size[0] * (i - 1)];
                        c_sqrt(&u2);
                        u_re = y.re * u2.re - y.im * u2.im;
                        u_im = y.re * u2.im + y.im * u2.re;
                        s = fabs(u_re) + fabs(u_im);
                        if (s != 0.0) {
                            tst = h->data[(i + h->size[0] * (i - 1)) - 1].re - h->data[i +
                                    h->size[0] * i].re;
                            htmp1 = h->data[(i + h->size[0] * (i - 1)) - 1].im - h->data[i +
                                    h->size[0] * i].im;
                            x_re = 0.5 * tst;
                            x_im = 0.5 * htmp1;
                            aa = fabs(x_re) + fabs(x_im);
                            tst = fabs(x_re) + fabs(x_im);
                            if (!((s > tst) || rtIsNaN(tst))) {
                                s = tst;
                            }

                            if (x_im == 0.0) {
                                ba = x_re / s;
                                ab = 0.0;
                            } else if (x_re == 0.0) {
                                ba = 0.0;
                                ab = x_im / s;
                            } else {
                                ba = x_re / s;
                                ab = x_im / s;
                            }

                            tst = ba;
                            ba = ba * ba - ab * ab;
                            ab = tst * ab + ab * tst;
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

                            y.re = ba + (u2.re * u2.re - u2.im * u2.im);
                            y.im = ab + (u2.re * u2.im + u2.im * u2.re);
                            c_sqrt(&y);
                            y.re *= s;
                            y.im *= s;
                            if (aa > 0.0) {
                                if (x_im == 0.0) {
                                    ba = x_re / aa;
                                    ab = 0.0;
                                } else if (x_re == 0.0) {
                                    ba = 0.0;
                                    ab = x_im / aa;
                                } else {
                                    ba = x_re / aa;
                                    ab = x_im / aa;
                                }

                                if (ba * y.re + ab * y.im < 0.0) {
                                    y.re = -y.re;
                                    y.im = -y.im;
                                }
                            }

                            ba = x_re + y.re;
                            htmp1 = x_im + y.im;
                            if (htmp1 == 0.0) {
                                if (u_im == 0.0) {
                                    x_re = u_re / ba;
                                    tst = 0.0;
                                } else if (u_re == 0.0) {
                                    x_re = 0.0;
                                    tst = u_im / ba;
                                } else {
                                    x_re = u_re / ba;
                                    tst = u_im / ba;
                                }
                            } else if (ba == 0.0) {
                                if (u_re == 0.0) {
                                    x_re = u_im / htmp1;
                                    tst = 0.0;
                                } else if (u_im == 0.0) {
                                    x_re = 0.0;
                                    tst = -(u_re / htmp1);
                                } else {
                                    x_re = u_im / htmp1;
                                    tst = -(u_re / htmp1);
                                }
                            } else {
                                ab = fabs(ba);
                                tst = fabs(htmp1);
                                if (ab > tst) {
                                    s = htmp1 / ba;
                                    tst = ba + s * htmp1;
                                    x_re = (u_re + s * u_im) / tst;
                                    tst = (u_im - s * u_re) / tst;
                                } else if (tst == ab) {
                                    if (ba > 0.0) {
                                        aa = 0.5;
                                    } else {
                                        aa = -0.5;
                                    }

                                    if (htmp1 > 0.0) {
                                        tst = 0.5;
                                    } else {
                                        tst = -0.5;
                                    }

                                    x_re = (u_re * aa + u_im * tst) / ab;
                                    tst = (u_im * aa - u_re * tst) / ab;
                                } else {
                                    s = ba / htmp1;
                                    tst = htmp1 + s * ba;
                                    x_re = (s * u_re + u_im) / tst;
                                    tst = (s * u_im - u_re) / tst;
                                }
                            }

                            ba = h->data[i + h->size[0] * i].re - (u_re * x_re - u_im * tst);
                            ab = h->data[i + h->size[0] * i].im - (u_re * tst + u_im * x_re);
                        }
                    }

                    goto70 = false;
                    m = i;
                    exitg3 = false;
                    while ((!exitg3) && (m > k + 1)) {
                        u2.re = h->data[(m + h->size[0] * (m - 1)) - 1].re - ba;
                        u2.im = h->data[(m + h->size[0] * (m - 1)) - 1].im - ab;
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
                                                                      h->size[0] * m].re) + fabs(h->data[m + h->size[0] * m].im))))) {
                            goto70 = true;
                            exitg3 = true;
                        } else {
                            m--;
                        }
                    }

                    if (!goto70) {
                        u2.re = h->data[k + h->size[0] * k].re - ba;
                        u2.im = h->data[k + h->size[0] * k].im - ab;
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

                        u2 = xzlarfg(&v[0], &v[1]);
                        if (b_k > m) {
                            h->data[(b_k + h->size[0] * (b_k - 2)) - 1] = v[0];
                            h->data[b_k + h->size[0] * (b_k - 2)].re = 0.0;
                            h->data[b_k + h->size[0] * (b_k - 2)].im = 0.0;
                        }

                        htmp1 = u2.re * v[1].re - u2.im * v[1].im;
                        for (maxval = b_k - 1; maxval + 1 <= n; maxval++) {
                            tst = u2.re * h->data[(b_k + h->size[0] * maxval) - 1].re - -u2.im
                                  * h->data[(b_k + h->size[0] * maxval) - 1].im;
                            aa = u2.re * h->data[(b_k + h->size[0] * maxval) - 1].im + -u2.im *
                                 h->data[(b_k + h->size[0] * maxval) - 1].re;
                            ba = tst + htmp1 * h->data[b_k + h->size[0] * maxval].re;
                            ab = aa + htmp1 * h->data[b_k + h->size[0] * maxval].im;
                            h->data[(b_k + h->size[0] * maxval) - 1].re -= ba;
                            h->data[(b_k + h->size[0] * maxval) - 1].im -= ab;
                            h->data[b_k + h->size[0] * maxval].re -= ba * v[1].re - ab * v[1].
                                    im;
                            h->data[b_k + h->size[0] * maxval].im -= ba * v[1].im + ab * v[1].
                                    re;
                        }

                        if (b_k + 2 < i + 1) {
                            i61 = b_k;
                        } else {
                            i61 = i - 1;
                        }

                        for (maxval = 0; maxval + 1 <= i61 + 2; maxval++) {
                            tst = u2.re * h->data[maxval + h->size[0] * (b_k - 1)].re - u2.im *
                                  h->data[maxval + h->size[0] * (b_k - 1)].im;
                            aa = u2.re * h->data[maxval + h->size[0] * (b_k - 1)].im + u2.im *
                                 h->data[maxval + h->size[0] * (b_k - 1)].re;
                            ba = tst + htmp1 * h->data[maxval + h->size[0] * b_k].re;
                            ab = aa + htmp1 * h->data[maxval + h->size[0] * b_k].im;
                            h->data[maxval + h->size[0] * (b_k - 1)].re -= ba;
                            h->data[maxval + h->size[0] * (b_k - 1)].im -= ab;
                            h->data[maxval + h->size[0] * b_k].re -= ba * v[1].re - ab * -v[1]
                                    .im;
                            h->data[maxval + h->size[0] * b_k].im -= ba * -v[1].im + ab * v[1]
                                    .re;
                        }

                        if ((b_k == m) && (m > k + 1)) {
                            u2.re = 1.0 - u2.re;
                            u2.im = 0.0 - u2.im;
                            ba = rt_hypotd_snf(u2.re, u2.im);
                            if (u2.im == 0.0) {
                                u2.re /= ba;
                                u2.im = 0.0;
                            } else if (u2.re == 0.0) {
                                u2.re = 0.0;
                                u2.im /= ba;
                            } else {
                                u2.re /= ba;
                                u2.im /= ba;
                            }

                            tst = h->data[m + h->size[0] * (m - 1)].re;
                            htmp1 = h->data[m + h->size[0] * (m - 1)].im;
                            h->data[m + h->size[0] * (m - 1)].re = tst * u2.re - htmp1 *
                                                                   -u2.im;
                            h->data[m + h->size[0] * (m - 1)].im = tst * -u2.im + htmp1 *
                                                                   u2.re;
                            if (m + 2 <= i + 1) {
                                tst = h->data[(m + h->size[0] * m) + 1].re;
                                htmp1 = h->data[(m + h->size[0] * m) + 1].im;
                                h->data[(m + h->size[0] * m) + 1].re = tst * u2.re - htmp1 *
                                                                       u2.im;
                                h->data[(m + h->size[0] * m) + 1].im = tst * u2.im + htmp1 *
                                                                       u2.re;
                            }

                            for (maxval = m; maxval <= i + 1; maxval++) {
                                if (maxval != m + 1) {
                                    if (n > maxval) {
                                        b_xscal(n - maxval, u2, h, maxval + maxval * ldh, ldh);
                                    }

                                    y.re = u2.re;
                                    y.im = -u2.im;
                                    xscal(maxval - 1, y, h, 1 + (maxval - 1) * ldh);
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
                            y.re = u2.re;
                            y.im = -u2.im;
                            b_xscal((n - i) - 1, y, h, (i + (i + 1) * ldh) + 1, ldh);
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
    int i67;
    creal_T dcv6[2048];
    double digw;
    double b_digw[2048];
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
    for (i67 = 0; i67 < 2048; i67++) {
        w[i67] = options->w[i67];
        digw = 6.2831853071795862 * options->w[i67] / options->Fs;
        dcv6[i67].re = digw * 0.0;
        dcv6[i67].im = digw;
        b_digw[i67] = digw;
    }

    b_exp(dcv6);
    e_polyval(b, dcv6, h);
    for (i67 = 0; i67 < 2048; i67++) {
        dcv6[i67].re = 12.0 * (b_digw[i67] * 0.0);
        dcv6[i67].im = 12.0 * b_digw[i67];
    }

    b_exp(dcv6);
    for (i67 = 0; i67 < 2048; i67++) {
        h_re = h[i67].re;
        if (dcv6[i67].im == 0.0) {
            if (h[i67].im == 0.0) {
                h[i67].re /= dcv6[i67].re;
                h[i67].im = 0.0;
            } else if (h[i67].re == 0.0) {
                h[i67].re = 0.0;
                h[i67].im /= dcv6[i67].re;
            } else {
                h[i67].re /= dcv6[i67].re;
                h[i67].im /= dcv6[i67].re;
            }
        } else if (dcv6[i67].re == 0.0) {
            if (h[i67].re == 0.0) {
                h[i67].re = h[i67].im / dcv6[i67].im;
                h[i67].im = 0.0;
            } else if (h[i67].im == 0.0) {
                h[i67].re = 0.0;
                h[i67].im = -(h_re / dcv6[i67].im);
            } else {
                h[i67].re = h[i67].im / dcv6[i67].im;
                h[i67].im = -(h_re / dcv6[i67].im);
            }
        } else {
            brm = fabs(dcv6[i67].re);
            digw = fabs(dcv6[i67].im);
            if (brm > digw) {
                digw = dcv6[i67].im / dcv6[i67].re;
                d = dcv6[i67].re + digw * dcv6[i67].im;
                h[i67].re = (h[i67].re + digw * h[i67].im) / d;
                h[i67].im = (h[i67].im - digw * h_re) / d;
            } else if (digw == brm) {
                if (dcv6[i67].re > 0.0) {
                    digw = 0.5;
                } else {
                    digw = -0.5;
                }

                if (dcv6[i67].im > 0.0) {
                    d = 0.5;
                } else {
                    d = -0.5;
                }

                h[i67].re = (h[i67].re * digw + h[i67].im * d) / brm;
                h[i67].im = (h[i67].im * digw - h_re * d) / brm;
            } else {
                digw = dcv6[i67].re / dcv6[i67].im;
                d = dcv6[i67].im + digw * dcv6[i67].re;
                h[i67].re = (digw * h[i67].re + h[i67].im) / d;
                h[i67].im = (digw * h[i67].im - h_re) / d;
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
 * Arguments    : const double p[57]
 *                const creal_T x[2048]
 *                creal_T y[2048]
 * Return Type  : void
 */
static void f_polyval(const double p[57], const creal_T x[2048],
                      creal_T y[2048])
{
    int i18;
    int k;
    double x_im;
    for (i18 = 0; i18 < 2048; i18++) {
        y[i18].re = p[0];
        y[i18].im = 0.0;
    }

    for (k = 0; k < 56; k++) {
        for (i18 = 0; i18 < 2048; i18++) {
            x_im = x[i18].re * y[i18].im + x[i18].im * y[i18].re;
            y[i18].re = (x[i18].re * y[i18].re - x[i18].im * y[i18].im) + p[k + 1];
            y[i18].im = x_im;
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
    int i63;
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
    for (i63 = 0; i63 < 2048; i63++) {
        w[i63] = options->w[i63];
        h[i63].re = 0.0 * (6.2831853071795862 * options->w[i63] / options->Fs * 0.0);
        h[i63].im = 0.0 * (6.2831853071795862 * options->w[i63] / options->Fs);
    }

    b_exp(h);
    for (i63 = 0; i63 < 2048; i63++) {
        if (h[i63].im == 0.0) {
            h[i63].re = 1.0 / h[i63].re;
            h[i63].im = 0.0;
        } else if (h[i63].re == 0.0) {
            h[i63].re = 0.0;
            h[i63].im = -(1.0 / h[i63].im);
        } else {
            brm = fabs(h[i63].re);
            bim = fabs(h[i63].im);
            if (brm > bim) {
                bim = h[i63].im / h[i63].re;
                d = h[i63].re + bim * h[i63].im;
                h[i63].re = (1.0 + bim * 0.0) / d;
                h[i63].im = (0.0 - bim) / d;
            } else if (bim == brm) {
                if (h[i63].re > 0.0) {
                    bim = 0.5;
                } else {
                    bim = -0.5;
                }

                if (h[i63].im > 0.0) {
                    d = 0.5;
                } else {
                    d = -0.5;
                }

                h[i63].re = (bim + 0.0 * d) / brm;
                h[i63].im = (0.0 * bim - d) / brm;
            } else {
                bim = h[i63].re / h[i63].im;
                d = h[i63].im + bim * h[i63].re;
                h[i63].re = bim / d;
                h[i63].im = (bim * 0.0 - 1.0) / d;
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
    emxArray_real_T *b_h;
    int i41;
    emxArray_real_T *c_h;
    double b_ff[4];
    double err;
    boolean_T valid;
    int h_idx_0;
    int i42;
    int i43;
    int loop_ub;
    emxInit_real_T(&grid, 2);
    emxInit_real_T(&des, 2);
    emxInit_real_T(&wt, 2);
    emxInit_real_T(&b_h, 2);

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
    interp1(frequencies, amplitudes, grid, des);
    interp1(frequencies, weights, grid, wt);

    /*  Workaround */
    /* ftype = 2; */
    /* sign_val = 1; */
    /*  Always bandpass designs */
    /*  cast to enforce precision rules */
    /*  Call actual design algorithm */
    rdivide(grid, 2.0, b_h);
    emxFree_real_T(&grid);
    for (i41 = 0; i41 < 4; i41++) {
        b_ff[i41] = ff[i41] / 2.0;
    }

    emxInit_real_T(&c_h, 2);
    remezm(order + 1.0, b_ff, b_h, des, wt, c_h, &err, &valid);
    h_idx_0 = c_h->size[0] * c_h->size[1];
    i41 = h->size[0] * h->size[1];
    h->size[0] = 1;
    h->size[1] = h_idx_0;
    emxEnsureCapacity_real_T(h, i41);
    emxFree_real_T(&wt);
    emxFree_real_T(&des);
    for (i41 = 0; i41 < h_idx_0; i41++) {
        h->data[h->size[0] * i41] = c_h->data[i41];
    }

    emxFree_real_T(&c_h);

    /*  make it a row */
    err = (double)h->size[1] - rt_remd_snf(order + 1.0, 2.0);
    if (1.0 > err) {
        i41 = 1;
        h_idx_0 = 1;
        i42 = 0;
    } else {
        i41 = (int)err;
        h_idx_0 = -1;
        i42 = 1;
    }

    i43 = b_h->size[0] * b_h->size[1];
    b_h->size[0] = 1;
    b_h->size[1] = (h->size[1] + div_s32_floor(i42 - i41, h_idx_0)) + 1;
    emxEnsureCapacity_real_T(b_h, i43);
    loop_ub = h->size[1];
    for (i43 = 0; i43 < loop_ub; i43++) {
        b_h->data[b_h->size[0] * i43] = h->data[h->size[0] * i43];
    }

    loop_ub = div_s32_floor(i42 - i41, h_idx_0);
    for (i42 = 0; i42 <= loop_ub; i42++) {
        b_h->data[b_h->size[0] * (i42 + h->size[1])] = h->data[(i41 + h_idx_0 * i42)
                - 1];
    }

    i41 = h->size[0] * h->size[1];
    h->size[0] = 1;
    h->size[1] = b_h->size[1];
    emxEnsureCapacity_real_T(h, i41);
    loop_ub = b_h->size[1];
    for (i41 = 0; i41 < loop_ub; i41++) {
        h->data[h->size[0] * i41] = b_h->data[b_h->size[0] * i41];
    }

    if (1 > h->size[1]) {
        i41 = 1;
        h_idx_0 = 1;
        i42 = 0;
    } else {
        i41 = h->size[1];
        h_idx_0 = -1;
        i42 = 1;
    }

    i43 = b_h->size[0] * b_h->size[1];
    b_h->size[0] = 1;
    b_h->size[1] = div_s32_floor(i42 - i41, h_idx_0) + 1;
    emxEnsureCapacity_real_T(b_h, i43);
    loop_ub = div_s32_floor(i42 - i41, h_idx_0);
    for (i42 = 0; i42 <= loop_ub; i42++) {
        b_h->data[b_h->size[0] * i42] = h->data[(i41 + h_idx_0 * i42) - 1];
    }

    i41 = h->size[0] * h->size[1];
    h->size[0] = 1;
    h->size[1] = b_h->size[1];
    emxEnsureCapacity_real_T(h, i41);
    loop_ub = b_h->size[1];
    for (i41 = 0; i41 < loop_ub; i41++) {
        h->data[h->size[0] * i41] = b_h->data[b_h->size[0] * i41];
    }

    emxFree_real_T(&b_h);
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
    emxArray_int32_T *r20;
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
    emxInit_int32_T(&r20, 2);
    while (l + 1.0 <= 4.0) {
        a = grid[(int)j - 1] + delf;
        ngrid = ff[(int)(l + 1.0) - 1] + delf;
        if (rtIsNaN(a) || rtIsNaN(delf) || rtIsNaN(ngrid)) {
            k = newgrid->size[0] * newgrid->size[1];
            newgrid->size[0] = 1;
            newgrid->size[1] = 1;
            emxEnsureCapacity_real_T(newgrid, k);
            newgrid->data[0] = rtNaN;
        } else if ((delf == 0.0) || ((a < ngrid) && (delf < 0.0)) || ((ngrid < a) &&
                   (delf > 0.0))) {
            k = newgrid->size[0] * newgrid->size[1];
            newgrid->size[0] = 1;
            newgrid->size[1] = 0;
            emxEnsureCapacity_real_T(newgrid, k);
        } else if ((rtIsInf(a) || rtIsInf(ngrid)) && (rtIsInf(delf) || (a == ngrid))) {
            k = newgrid->size[0] * newgrid->size[1];
            newgrid->size[0] = 1;
            newgrid->size[1] = 1;
            emxEnsureCapacity_real_T(newgrid, k);
            newgrid->data[0] = rtNaN;
        } else if (rtIsInf(delf)) {
            k = newgrid->size[0] * newgrid->size[1];
            newgrid->size[0] = 1;
            newgrid->size[1] = 1;
            emxEnsureCapacity_real_T(newgrid, k);
            newgrid->data[0] = a;
        } else if ((floor(a) == a) && (floor(delf) == delf)) {
            k = newgrid->size[0] * newgrid->size[1];
            newgrid->size[0] = 1;
            newgrid->size[1] = (int)floor((ngrid - a) / delf) + 1;
            emxEnsureCapacity_real_T(newgrid, k);
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
            if ((absa > absb) || rtIsNaN(absb)) {
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
            emxEnsureCapacity_real_T(newgrid, k);
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
                emxEnsureCapacity_real_T(newgrid, k);
                newgrid->data[0] = rtNaN;
            } else if ((delf1 == 0.0) || ((a < ngrid) && (delf1 < 0.0)) || ((ngrid < a)
                       && (delf1 > 0.0))) {
                k = newgrid->size[0] * newgrid->size[1];
                newgrid->size[0] = 1;
                newgrid->size[1] = 0;
                emxEnsureCapacity_real_T(newgrid, k);
            } else if ((rtIsInf(a) || rtIsInf(ngrid)) && (rtIsInf(delf1) || (a ==
                       ngrid))) {
                k = newgrid->size[0] * newgrid->size[1];
                newgrid->size[0] = 1;
                newgrid->size[1] = 1;
                emxEnsureCapacity_real_T(newgrid, k);
                newgrid->data[0] = rtNaN;
            } else if (rtIsInf(delf1)) {
                k = newgrid->size[0] * newgrid->size[1];
                newgrid->size[0] = 1;
                newgrid->size[1] = 1;
                emxEnsureCapacity_real_T(newgrid, k);
                newgrid->data[0] = a;
            } else if ((floor(a) == a) && (floor(delf1) == delf1)) {
                k = newgrid->size[0] * newgrid->size[1];
                newgrid->size[0] = 1;
                newgrid->size[1] = (int)floor((ngrid - a) / delf1) + 1;
                emxEnsureCapacity_real_T(newgrid, k);
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
                if ((absa > absb) || rtIsNaN(absb)) {
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
                emxEnsureCapacity_real_T(newgrid, k);
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
        nm1d2 = r20->size[0] * r20->size[1];
        r20->size[0] = 1;
        r20->size[1] = (int)((double)k - 1.0) + 1;
        emxEnsureCapacity_int32_T(r20, nm1d2);
        nm1d2 = (int)((double)k - 1.0);
        for (k = 0; k <= nm1d2; k++) {
            r20->data[r20->size[0] * k] = (int)(gridSize + (1.0 + (double)k));
        }

        nm1d2 = newgrid->size[0] * newgrid->size[1];
        for (k = 0; k < nm1d2; k++) {
            grid[r20->data[k] - 1] = newgrid->data[k];
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

    emxFree_int32_T(&r20);
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
    emxEnsureCapacity_real_T(gridactual, k);
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
    int i26;
    boolean_T b0;
    static creal_T s[2048];
    int k;
    creal_T y[2048];
    double h_re;
    double bim;
    double d;
    double brm;
    emxInit_creal_T(&b_a, 2);
    removeTrailingZero(b_data, b_size, a, b_b_data, b_b_size, b_a);
    for (i26 = 0; i26 < 2048; i26++) {
        s[i26].re = w[i26] * 0.0;
        s[i26].im = w[i26];
    }

    b0 = (b_a->size[1] == 0);
    if (!b0) {
        for (i26 = 0; i26 < 2048; i26++) {
            y[i26] = b_a->data[0];
        }

        for (k = 0; k <= b_a->size[1] - 2; k++) {
            bim = b_a->data[k + 1].re;
            d = b_a->data[k + 1].im;
            for (i26 = 0; i26 < 2048; i26++) {
                brm = s[i26].re * y[i26].im + s[i26].im * y[i26].re;
                y[i26].re = (s[i26].re * y[i26].re - s[i26].im * y[i26].im) + bim;
                y[i26].im = brm + d;
            }
        }
    }

    emxFree_creal_T(&b_a);
    c_polyval(b_b_data, b_b_size, s, h);
    for (i26 = 0; i26 < 2048; i26++) {
        h_re = h[i26].re;
        if (y[i26].im == 0.0) {
            if (h[i26].im == 0.0) {
                h[i26].re /= y[i26].re;
                h[i26].im = 0.0;
            } else if (h[i26].re == 0.0) {
                h[i26].re = 0.0;
                h[i26].im /= y[i26].re;
            } else {
                h[i26].re /= y[i26].re;
                h[i26].im /= y[i26].re;
            }
        } else if (y[i26].re == 0.0) {
            if (h[i26].re == 0.0) {
                h[i26].re = h[i26].im / y[i26].im;
                h[i26].im = 0.0;
            } else if (h[i26].im == 0.0) {
                h[i26].re = 0.0;
                h[i26].im = -(h_re / y[i26].im);
            } else {
                h[i26].re = h[i26].im / y[i26].im;
                h[i26].im = -(h_re / y[i26].im);
            }
        } else {
            brm = fabs(y[i26].re);
            bim = fabs(y[i26].im);
            if (brm > bim) {
                bim = y[i26].im / y[i26].re;
                d = y[i26].re + bim * y[i26].im;
                h[i26].re = (h[i26].re + bim * h[i26].im) / d;
                h[i26].im = (h[i26].im - bim * h_re) / d;
            } else if (bim == brm) {
                if (y[i26].re > 0.0) {
                    bim = 0.5;
                } else {
                    bim = -0.5;
                }

                if (y[i26].im > 0.0) {
                    d = 0.5;
                } else {
                    d = -0.5;
                }

                h[i26].re = (h[i26].re * bim + h[i26].im * d) / brm;
                h[i26].im = (h[i26].im * bim - h_re * d) / brm;
            } else {
                bim = y[i26].re / y[i26].im;
                d = y[i26].im + bim * y[i26].re;
                h[i26].re = (bim * h[i26].re + h[i26].im) / d;
                h[i26].im = (bim * h[i26].im - h_re) / d;
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
    int i68;
    creal_T dcv7[2048];
    double digw;
    double b_digw[2048];
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
    for (i68 = 0; i68 < 2048; i68++) {
        w[i68] = options->w[i68];
        digw = 6.2831853071795862 * options->w[i68] / options->Fs;
        dcv7[i68].re = digw * 0.0;
        dcv7[i68].im = digw;
        b_digw[i68] = digw;
    }

    b_exp(dcv7);
    f_polyval(b, dcv7, h);
    for (i68 = 0; i68 < 2048; i68++) {
        dcv7[i68].re = 56.0 * (b_digw[i68] * 0.0);
        dcv7[i68].im = 56.0 * b_digw[i68];
    }

    b_exp(dcv7);
    for (i68 = 0; i68 < 2048; i68++) {
        h_re = h[i68].re;
        if (dcv7[i68].im == 0.0) {
            if (h[i68].im == 0.0) {
                h[i68].re /= dcv7[i68].re;
                h[i68].im = 0.0;
            } else if (h[i68].re == 0.0) {
                h[i68].re = 0.0;
                h[i68].im /= dcv7[i68].re;
            } else {
                h[i68].re /= dcv7[i68].re;
                h[i68].im /= dcv7[i68].re;
            }
        } else if (dcv7[i68].re == 0.0) {
            if (h[i68].re == 0.0) {
                h[i68].re = h[i68].im / dcv7[i68].im;
                h[i68].im = 0.0;
            } else if (h[i68].im == 0.0) {
                h[i68].re = 0.0;
                h[i68].im = -(h_re / dcv7[i68].im);
            } else {
                h[i68].re = h[i68].im / dcv7[i68].im;
                h[i68].im = -(h_re / dcv7[i68].im);
            }
        } else {
            brm = fabs(dcv7[i68].re);
            digw = fabs(dcv7[i68].im);
            if (brm > digw) {
                digw = dcv7[i68].im / dcv7[i68].re;
                d = dcv7[i68].re + digw * dcv7[i68].im;
                h[i68].re = (h[i68].re + digw * h[i68].im) / d;
                h[i68].im = (h[i68].im - digw * h_re) / d;
            } else if (digw == brm) {
                if (dcv7[i68].re > 0.0) {
                    digw = 0.5;
                } else {
                    digw = -0.5;
                }

                if (dcv7[i68].im > 0.0) {
                    d = 0.5;
                } else {
                    d = -0.5;
                }

                h[i68].re = (h[i68].re * digw + h[i68].im * d) / brm;
                h[i68].im = (h[i68].im * digw - h_re * d) / brm;
            } else {
                digw = dcv7[i68].re / dcv7[i68].im;
                d = dcv7[i68].im + digw * dcv7[i68].re;
                h[i68].re = (digw * h[i68].re + h[i68].im) / d;
                h[i68].im = (digw * h[i68].im - h_re) / d;
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
 * Arguments    : const double p[43]
 *                const creal_T x[2048]
 *                creal_T y[2048]
 * Return Type  : void
 */
static void g_polyval(const double p[43], const creal_T x[2048],
                      creal_T y[2048])
{
    int i20;
    int k;
    double x_im;
    for (i20 = 0; i20 < 2048; i20++) {
        y[i20].re = p[0];
        y[i20].im = 0.0;
    }

    for (k = 0; k < 42; k++) {
        for (i20 = 0; i20 < 2048; i20++) {
            x_im = x[i20].re * y[i20].im + x[i20].im * y[i20].re;
            y[i20].re = (x[i20].re * y[i20].re - x[i20].im * y[i20].im) + p[k + 1];
            y[i20].im = x_im;
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
    boolean_T b_bool;
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
    kstr = 1;
    do {
        exitg1 = 0;
        if (kstr < 5) {
            if (enables[kstr - 1] != '1') {
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
            if (kstr + 1 < 5) {
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
                if (kstr + 1 < 5) {
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
                    if (kstr + 1 < 5) {
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
                        if (kstr + 1 < 5) {
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
                            if (kstr + 1 < 5) {
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
                                if (kstr + 1 < 5) {
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
                                    if (kstr + 1 < 5) {
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
                                        if (kstr + 1 < 5) {
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
                                            if (kstr + 1 < 5) {
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
                                                if (kstr + 1 < 5) {
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
    int i69;
    creal_T dcv8[2048];
    double digw;
    double b_digw[2048];
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
    for (i69 = 0; i69 < 2048; i69++) {
        w[i69] = options->w[i69];
        digw = 6.2831853071795862 * options->w[i69] / options->Fs;
        dcv8[i69].re = digw * 0.0;
        dcv8[i69].im = digw;
        b_digw[i69] = digw;
    }

    b_exp(dcv8);
    g_polyval(b, dcv8, h);
    for (i69 = 0; i69 < 2048; i69++) {
        dcv8[i69].re = 42.0 * (b_digw[i69] * 0.0);
        dcv8[i69].im = 42.0 * b_digw[i69];
    }

    b_exp(dcv8);
    for (i69 = 0; i69 < 2048; i69++) {
        h_re = h[i69].re;
        if (dcv8[i69].im == 0.0) {
            if (h[i69].im == 0.0) {
                h[i69].re /= dcv8[i69].re;
                h[i69].im = 0.0;
            } else if (h[i69].re == 0.0) {
                h[i69].re = 0.0;
                h[i69].im /= dcv8[i69].re;
            } else {
                h[i69].re /= dcv8[i69].re;
                h[i69].im /= dcv8[i69].re;
            }
        } else if (dcv8[i69].re == 0.0) {
            if (h[i69].re == 0.0) {
                h[i69].re = h[i69].im / dcv8[i69].im;
                h[i69].im = 0.0;
            } else if (h[i69].im == 0.0) {
                h[i69].re = 0.0;
                h[i69].im = -(h_re / dcv8[i69].im);
            } else {
                h[i69].re = h[i69].im / dcv8[i69].im;
                h[i69].im = -(h_re / dcv8[i69].im);
            }
        } else {
            brm = fabs(dcv8[i69].re);
            digw = fabs(dcv8[i69].im);
            if (brm > digw) {
                digw = dcv8[i69].im / dcv8[i69].re;
                d = dcv8[i69].re + digw * dcv8[i69].im;
                h[i69].re = (h[i69].re + digw * h[i69].im) / d;
                h[i69].im = (h[i69].im - digw * h_re) / d;
            } else if (digw == brm) {
                if (dcv8[i69].re > 0.0) {
                    digw = 0.5;
                } else {
                    digw = -0.5;
                }

                if (dcv8[i69].im > 0.0) {
                    d = 0.5;
                } else {
                    d = -0.5;
                }

                h[i69].re = (h[i69].re * digw + h[i69].im * d) / brm;
                h[i69].im = (h[i69].im * digw - h_re * d) / brm;
            } else {
                digw = dcv8[i69].re / dcv8[i69].im;
                d = dcv8[i69].im + digw * dcv8[i69].re;
                h[i69].re = (digw * h[i69].re + h[i69].im) / d;
                h[i69].im = (digw * h[i69].im - h_re) / d;
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
    int i21;
    static const char cv26[8] = { 'o', 'n', 'e', 's', 'i', 'd', 'e', 'd' };

    double b_w[2048];
    static const char cv27[7] = { 'o', 'm', 'i', 't', 't', 'e', 'd' };

    /*  Cast to enforce precision rules */
    options.Fs = Fs;
    memcpy(&options.nfft[0], &w[0], sizeof(double) << 11);
    memcpy(&options.w[0], &w[0], sizeof(double) << 11);

    /*  Remaining are default or for advanced use */
    options.fvflag = 1.0;
    for (i21 = 0; i21 < 8; i21++) {
        options.range[i21] = cv26[i21];
    }

    options.centerdc = 0.0;
    for (i21 = 0; i21 < 7; i21++) {
        options.configlevel[i21] = cv27[i21];
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
static void h_polyval(const double p[19], const creal_T x[2048],
                      creal_T y[2048])
{
    int i22;
    int k;
    double x_im;
    for (i22 = 0; i22 < 2048; i22++) {
        y[i22].re = p[0];
        y[i22].im = 0.0;
    }

    for (k = 0; k < 18; k++) {
        for (i22 = 0; i22 < 2048; i22++) {
            x_im = x[i22].re * y[i22].im + x[i22].im * y[i22].re;
            y[i22].re = (x[i22].re * y[i22].re - x[i22].im * y[i22].im) + p[k + 1];
            y[i22].im = x_im;
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
    int i70;
    creal_T dcv9[2048];
    double digw;
    double b_digw[2048];
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
    for (i70 = 0; i70 < 2048; i70++) {
        w[i70] = options->w[i70];
        digw = 6.2831853071795862 * options->w[i70] / options->Fs;
        dcv9[i70].re = digw * 0.0;
        dcv9[i70].im = digw;
        b_digw[i70] = digw;
    }

    b_exp(dcv9);
    h_polyval(b, dcv9, h);
    for (i70 = 0; i70 < 2048; i70++) {
        dcv9[i70].re = 18.0 * (b_digw[i70] * 0.0);
        dcv9[i70].im = 18.0 * b_digw[i70];
    }

    b_exp(dcv9);
    for (i70 = 0; i70 < 2048; i70++) {
        h_re = h[i70].re;
        if (dcv9[i70].im == 0.0) {
            if (h[i70].im == 0.0) {
                h[i70].re /= dcv9[i70].re;
                h[i70].im = 0.0;
            } else if (h[i70].re == 0.0) {
                h[i70].re = 0.0;
                h[i70].im /= dcv9[i70].re;
            } else {
                h[i70].re /= dcv9[i70].re;
                h[i70].im /= dcv9[i70].re;
            }
        } else if (dcv9[i70].re == 0.0) {
            if (h[i70].re == 0.0) {
                h[i70].re = h[i70].im / dcv9[i70].im;
                h[i70].im = 0.0;
            } else if (h[i70].im == 0.0) {
                h[i70].re = 0.0;
                h[i70].im = -(h_re / dcv9[i70].im);
            } else {
                h[i70].re = h[i70].im / dcv9[i70].im;
                h[i70].im = -(h_re / dcv9[i70].im);
            }
        } else {
            brm = fabs(dcv9[i70].re);
            digw = fabs(dcv9[i70].im);
            if (brm > digw) {
                digw = dcv9[i70].im / dcv9[i70].re;
                d = dcv9[i70].re + digw * dcv9[i70].im;
                h[i70].re = (h[i70].re + digw * h[i70].im) / d;
                h[i70].im = (h[i70].im - digw * h_re) / d;
            } else if (digw == brm) {
                if (dcv9[i70].re > 0.0) {
                    digw = 0.5;
                } else {
                    digw = -0.5;
                }

                if (dcv9[i70].im > 0.0) {
                    d = 0.5;
                } else {
                    d = -0.5;
                }

                h[i70].re = (h[i70].re * digw + h[i70].im * d) / brm;
                h[i70].im = (h[i70].im * digw - h_re * d) / brm;
            } else {
                digw = dcv9[i70].re / dcv9[i70].im;
                d = dcv9[i70].im + digw * dcv9[i70].re;
                h[i70].re = (digw * h[i70].re + h[i70].im) / d;
                h[i70].im = (digw * h[i70].im - h_re) / d;
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
    int i23;
    static const char cv28[8] = { 'o', 'n', 'e', 's', 'i', 'd', 'e', 'd' };

    double b_w[2048];
    static const char cv29[7] = { 'o', 'm', 'i', 't', 't', 'e', 'd' };

    /*  Cast to enforce precision rules */
    options.Fs = Fs;
    memcpy(&options.nfft[0], &w[0], sizeof(double) << 11);
    memcpy(&options.w[0], &w[0], sizeof(double) << 11);

    /*  Remaining are default or for advanced use */
    options.fvflag = 1.0;
    for (i23 = 0; i23 < 8; i23++) {
        options.range[i23] = cv28[i23];
    }

    options.centerdc = 0.0;
    for (i23 = 0; i23 < 7; i23++) {
        options.configlevel[i23] = cv29[i23];
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
static void i_polyval(const double p[85], const creal_T x[2048],
                      creal_T y[2048])
{
    int i24;
    int k;
    double x_im;
    for (i24 = 0; i24 < 2048; i24++) {
        y[i24].re = p[0];
        y[i24].im = 0.0;
    }

    for (k = 0; k < 84; k++) {
        for (i24 = 0; i24 < 2048; i24++) {
            x_im = x[i24].re * y[i24].im + x[i24].im * y[i24].re;
            y[i24].re = (x[i24].re * y[i24].re - x[i24].im * y[i24].im) + p[k + 1];
            y[i24].im = x_im;
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
    int low_ip1;
    int nd2;
    emxArray_real_T *x;
    int nx;
    int outsize[2];
    int k;
    int exitg1;
    int mid_i;
    double r;
    emxInit_real_T(&y, 2);
    low_ip1 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = varargin_2->size[1];
    emxEnsureCapacity_real_T(y, low_ip1);
    nd2 = varargin_2->size[0] * varargin_2->size[1];
    for (low_ip1 = 0; low_ip1 < nd2; low_ip1++) {
        y->data[low_ip1] = varargin_2->data[low_ip1];
    }

    emxInit_real_T(&x, 2);
    low_ip1 = x->size[0] * x->size[1];
    x->size[0] = 1;
    x->size[1] = varargin_1->size[1];
    emxEnsureCapacity_real_T(x, low_ip1);
    nd2 = varargin_1->size[0] * varargin_1->size[1];
    for (low_ip1 = 0; low_ip1 < nd2; low_ip1++) {
        x->data[low_ip1] = varargin_1->data[low_ip1];
    }

    nx = varargin_1->size[1];
    for (low_ip1 = 0; low_ip1 < 2; low_ip1++) {
        outsize[low_ip1] = varargin_3->size[low_ip1];
    }

    low_ip1 = Vq->size[0] * Vq->size[1];
    Vq->size[0] = 1;
    Vq->size[1] = outsize[1];
    emxEnsureCapacity_real_T(Vq, low_ip1);
    nd2 = outsize[1];
    for (low_ip1 = 0; low_ip1 < nd2; low_ip1++) {
        Vq->data[low_ip1] = rtNaN;
    }

    if (varargin_3->size[1] != 0) {
        k = 1;
        do {
            exitg1 = 0;
            if (k <= nx) {
                if (rtIsNaN(varargin_1->data[k - 1])) {
                    exitg1 = 1;
                } else {
                    k++;
                }
            } else {
                if (varargin_1->data[1] < varargin_1->data[0]) {
                    low_ip1 = nx >> 1;
                    for (mid_i = 1; mid_i <= low_ip1; mid_i++) {
                        r = x->data[mid_i - 1];
                        x->data[mid_i - 1] = x->data[nx - mid_i];
                        x->data[nx - mid_i] = r;
                    }

                    nd2 = varargin_2->size[1] >> 1;
                    for (mid_i = 1; mid_i <= nd2; mid_i++) {
                        low_ip1 = varargin_2->size[1] - mid_i;
                        r = y->data[mid_i - 1];
                        y->data[mid_i - 1] = y->data[low_ip1];
                        y->data[low_ip1] = r;
                    }
                }

                for (k = 0; k + 1 <= varargin_3->size[1]; k++) {
                    if (rtIsNaN(varargin_3->data[k])) {
                        Vq->data[k] = rtNaN;
                    } else {
                        if ((!(varargin_3->data[k] > x->data[x->size[1] - 1])) &&
                            (!(varargin_3->data[k] < x->data[0]))) {
                            nd2 = 1;
                            low_ip1 = 2;
                            nx = x->size[1];
                            while (nx > low_ip1) {
                                mid_i = (nd2 >> 1) + (nx >> 1);
                                if (((nd2 & 1) == 1) && ((nx & 1) == 1)) {
                                    mid_i++;
                                }

                                if (varargin_3->data[k] >= x->data[mid_i - 1]) {
                                    nd2 = mid_i;
                                    low_ip1 = mid_i + 1;
                                } else {
                                    nx = mid_i;
                                }
                            }

                            r = (varargin_3->data[k] - x->data[nd2 - 1]) / (x->data[nd2] -
                                    x->data[nd2 - 1]);
                            if (r == 0.0) {
                                Vq->data[k] = y->data[nd2 - 1];
                            } else if (r == 1.0) {
                                Vq->data[k] = y->data[nd2];
                            } else if (y->data[nd2 - 1] == y->data[nd2]) {
                                Vq->data[k] = y->data[nd2 - 1];
                            } else {
                                Vq->data[k] = (1.0 - r) * y->data[nd2 - 1] + r * y->data[nd2];
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
    int i71;
    creal_T dcv10[2048];
    double digw;
    double b_digw[2048];
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
    for (i71 = 0; i71 < 2048; i71++) {
        w[i71] = options->w[i71];
        digw = 6.2831853071795862 * options->w[i71] / options->Fs;
        dcv10[i71].re = digw * 0.0;
        dcv10[i71].im = digw;
        b_digw[i71] = digw;
    }

    b_exp(dcv10);
    i_polyval(b, dcv10, h);
    for (i71 = 0; i71 < 2048; i71++) {
        dcv10[i71].re = 84.0 * (b_digw[i71] * 0.0);
        dcv10[i71].im = 84.0 * b_digw[i71];
    }

    b_exp(dcv10);
    for (i71 = 0; i71 < 2048; i71++) {
        h_re = h[i71].re;
        if (dcv10[i71].im == 0.0) {
            if (h[i71].im == 0.0) {
                h[i71].re /= dcv10[i71].re;
                h[i71].im = 0.0;
            } else if (h[i71].re == 0.0) {
                h[i71].re = 0.0;
                h[i71].im /= dcv10[i71].re;
            } else {
                h[i71].re /= dcv10[i71].re;
                h[i71].im /= dcv10[i71].re;
            }
        } else if (dcv10[i71].re == 0.0) {
            if (h[i71].re == 0.0) {
                h[i71].re = h[i71].im / dcv10[i71].im;
                h[i71].im = 0.0;
            } else if (h[i71].im == 0.0) {
                h[i71].re = 0.0;
                h[i71].im = -(h_re / dcv10[i71].im);
            } else {
                h[i71].re = h[i71].im / dcv10[i71].im;
                h[i71].im = -(h_re / dcv10[i71].im);
            }
        } else {
            brm = fabs(dcv10[i71].re);
            digw = fabs(dcv10[i71].im);
            if (brm > digw) {
                digw = dcv10[i71].im / dcv10[i71].re;
                d = dcv10[i71].re + digw * dcv10[i71].im;
                h[i71].re = (h[i71].re + digw * h[i71].im) / d;
                h[i71].im = (h[i71].im - digw * h_re) / d;
            } else if (digw == brm) {
                if (dcv10[i71].re > 0.0) {
                    digw = 0.5;
                } else {
                    digw = -0.5;
                }

                if (dcv10[i71].im > 0.0) {
                    d = 0.5;
                } else {
                    d = -0.5;
                }

                h[i71].re = (h[i71].re * digw + h[i71].im * d) / brm;
                h[i71].im = (h[i71].im * digw - h_re * d) / brm;
            } else {
                digw = dcv10[i71].re / dcv10[i71].im;
                d = dcv10[i71].im + digw * dcv10[i71].re;
                h[i71].re = (digw * h[i71].re + h[i71].im) / d;
                h[i71].im = (digw * h[i71].im - h_re) / d;
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
    int i25;
    static const char cv30[8] = { 'o', 'n', 'e', 's', 'i', 'd', 'e', 'd' };

    double b_w[2048];
    static const char cv31[7] = { 'o', 'm', 'i', 't', 't', 'e', 'd' };

    /*  Cast to enforce precision rules */
    options.Fs = Fs;
    memcpy(&options.nfft[0], &w[0], sizeof(double) << 11);
    memcpy(&options.w[0], &w[0], sizeof(double) << 11);

    /*  Remaining are default or for advanced use */
    options.fvflag = 1.0;
    for (i25 = 0; i25 < 8; i25++) {
        options.range[i25] = cv30[i25];
    }

    options.centerdc = 0.0;
    for (i25 = 0; i25 < 7; i25++) {
        options.configlevel[i25] = cv31[i25];
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
    int i33;
    boolean_T b4;
    int loop_ub;
    int k;
    double p_re;
    double x_re;
    double x_im;
    i33 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = x->size[1];
    emxEnsureCapacity_creal_T(y, i33);
    if ((y->size[1] == 0) || (p->size[1] == 0)) {
        b4 = true;
    } else {
        b4 = false;
    }

    if (!b4) {
        i33 = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity_creal_T(y, i33);
        loop_ub = y->size[1];
        for (i33 = 0; i33 < loop_ub; i33++) {
            y->data[y->size[0] * i33].re = p->data[0];
            y->data[y->size[0] * i33].im = 0.0;
        }

        for (k = 0; k <= p->size[1] - 2; k++) {
            i33 = y->size[0] * y->size[1];
            y->size[0] = 1;
            y->size[1] = x->size[1];
            emxEnsureCapacity_creal_T(y, i33);
            p_re = p->data[k + 1];
            loop_ub = x->size[0] * x->size[1];
            for (i33 = 0; i33 < loop_ub; i33++) {
                x_re = x->data[i33].re * y->data[i33].re - x->data[i33].im * y->data[i33]
                       .im;
                x_im = x->data[i33].re * y->data[i33].im + x->data[i33].im * y->data[i33]
                       .re;
                y->data[i33].re = x_re + p_re;
                y->data[i33].im = x_im;
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
static void k_freqz_cg(const emxArray_real_T *w, double Fs,
                       emxArray_creal_T *hh)
{
    emxArray_real_T *r4;
    int i28;
    int loop_ub;
    emxArray_real_T *digw;
    emxArray_creal_T *y;
    emxArray_creal_T *r5;
    boolean_T b1;
    double re;
    double im;
    emxInit_real_T(&r4, 2);

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
    i28 = r4->size[0] * r4->size[1];
    r4->size[0] = 1;
    r4->size[1] = w->size[1];
    emxEnsureCapacity_real_T(r4, i28);
    loop_ub = w->size[0] * w->size[1];
    for (i28 = 0; i28 < loop_ub; i28++) {
        r4->data[i28] = 6.2831853071795862 * w->data[i28];
    }

    emxInit_real_T(&digw, 2);
    emxInit_creal_T(&y, 2);
    emxInit_creal_T(&r5, 2);
    rdivide(r4, Fs, digw);

    /*  Convert from Hz to rad/sample for computational purposes */
    /*  Digital frequency must be used for this calculation */
    i28 = r5->size[0] * r5->size[1];
    r5->size[0] = 1;
    r5->size[1] = digw->size[1];
    emxEnsureCapacity_creal_T(r5, i28);
    loop_ub = digw->size[0] * digw->size[1];
    emxFree_real_T(&r4);
    for (i28 = 0; i28 < loop_ub; i28++) {
        r5->data[i28].re = digw->data[i28] * 0.0;
        r5->data[i28].im = digw->data[i28];
    }

    c_exp(r5);
    i28 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = r5->size[1];
    emxEnsureCapacity_creal_T(y, i28);
    b1 = (y->size[1] == 0);
    if (!b1) {
        i28 = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity_creal_T(y, i28);
        loop_ub = y->size[1];
        for (i28 = 0; i28 < loop_ub; i28++) {
            y->data[y->size[0] * i28].re = 1.0;
            y->data[y->size[0] * i28].im = 0.0;
        }
    }

    i28 = r5->size[0] * r5->size[1];
    r5->size[0] = 1;
    r5->size[1] = digw->size[1];
    emxEnsureCapacity_creal_T(r5, i28);
    loop_ub = digw->size[0] * digw->size[1];
    for (i28 = 0; i28 < loop_ub; i28++) {
        re = digw->data[i28] * 0.0;
        im = digw->data[i28];
        r5->data[i28].re = 0.0 * re;
        r5->data[i28].im = 0.0 * im;
    }

    emxFree_real_T(&digw);
    c_exp(r5);
    b_rdivide(y, r5, hh);

    /*  Generate the default structure to pass to freqzplot */
    /*  If rad/sample, Fs is empty */
    emxFree_creal_T(&r5);
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
    emxArray_real_T *r6;
    int i30;
    int loop_ub;
    emxArray_real_T *digw;
    emxArray_creal_T *s;
    emxArray_creal_T *y;
    boolean_T b2;
    int k;
    double s_re;
    double s_im;
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
    i30 = r6->size[0] * r6->size[1];
    r6->size[0] = 1;
    r6->size[1] = w->size[1];
    emxEnsureCapacity_real_T(r6, i30);
    loop_ub = w->size[0] * w->size[1];
    for (i30 = 0; i30 < loop_ub; i30++) {
        r6->data[i30] = 6.2831853071795862 * w->data[i30];
    }

    emxInit_real_T(&digw, 2);
    emxInit_creal_T(&s, 2);
    rdivide(r6, Fs, digw);

    /*  Convert from Hz to rad/sample for computational purposes */
    i30 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = digw->size[1];
    emxEnsureCapacity_creal_T(s, i30);
    loop_ub = digw->size[0] * digw->size[1];
    emxFree_real_T(&r6);
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
    emxEnsureCapacity_creal_T(y, i30);
    b2 = (y->size[1] == 0);
    if (!b2) {
        i30 = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity_creal_T(y, i30);
        loop_ub = y->size[1];
        for (i30 = 0; i30 < loop_ub; i30++) {
            y->data[y->size[0] * i30].re = b[0];
            y->data[y->size[0] * i30].im = 0.0;
        }

        for (k = 0; k < 14; k++) {
            i30 = y->size[0] * y->size[1];
            y->size[0] = 1;
            y->size[1] = s->size[1];
            emxEnsureCapacity_creal_T(y, i30);
            loop_ub = s->size[0] * s->size[1];
            for (i30 = 0; i30 < loop_ub; i30++) {
                s_re = s->data[i30].re * y->data[i30].re - s->data[i30].im * y->data[i30]
                       .im;
                s_im = s->data[i30].re * y->data[i30].im + s->data[i30].im * y->data[i30]
                       .re;
                y->data[i30].re = s_re + b[k + 1];
                y->data[i30].im = s_im;
            }
        }
    }

    i30 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = digw->size[1];
    emxEnsureCapacity_creal_T(s, i30);
    loop_ub = digw->size[0] * digw->size[1];
    for (i30 = 0; i30 < loop_ub; i30++) {
        s_re = digw->data[i30] * 0.0;
        s_im = digw->data[i30];
        s->data[i30].re = 14.0 * s_re;
        s->data[i30].im = 14.0 * s_im;
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
 *                double wo
 *                emxArray_creal_T *at
 *                emxArray_real_T *bt
 *                double *dt
 * Return Type  : void
 */
static void lp2lp_cg(const emxArray_creal_T *a, const emxArray_real_T *b, double
                     wo, emxArray_creal_T *at, emxArray_real_T *bt, double *dt)
{
    int i58;
    int loop_ub;

    /*  Transform lowpass to lowpass */
    i58 = at->size[0] * at->size[1];
    at->size[0] = a->size[0];
    at->size[1] = a->size[1];
    emxEnsureCapacity_creal_T(at, i58);
    loop_ub = a->size[0] * a->size[1];
    for (i58 = 0; i58 < loop_ub; i58++) {
        at->data[i58].re = wo * a->data[i58].re;
        at->data[i58].im = wo * a->data[i58].im;
    }

    i58 = bt->size[0];
    bt->size[0] = b->size[0];
    emxEnsureCapacity_real_T1(bt, i58);
    loop_ub = b->size[0];
    for (i58 = 0; i58 < loop_ub; i58++) {
        bt->data[i58] = wo * b->data[i58];
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
    emxArray_real_T *r7;
    int i31;
    int loop_ub;
    emxArray_real_T *digw;
    emxArray_creal_T *s;
    emxArray_creal_T *y;
    boolean_T b3;
    int k;
    double s_re;
    double s_im;
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
    i31 = r7->size[0] * r7->size[1];
    r7->size[0] = 1;
    r7->size[1] = w->size[1];
    emxEnsureCapacity_real_T(r7, i31);
    loop_ub = w->size[0] * w->size[1];
    for (i31 = 0; i31 < loop_ub; i31++) {
        r7->data[i31] = 6.2831853071795862 * w->data[i31];
    }

    emxInit_real_T(&digw, 2);
    emxInit_creal_T(&s, 2);
    rdivide(r7, Fs, digw);

    /*  Convert from Hz to rad/sample for computational purposes */
    i31 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = digw->size[1];
    emxEnsureCapacity_creal_T(s, i31);
    loop_ub = digw->size[0] * digw->size[1];
    emxFree_real_T(&r7);
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
    b3 = (y->size[1] == 0);
    if (!b3) {
        i31 = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity_creal_T(y, i31);
        loop_ub = y->size[1];
        for (i31 = 0; i31 < loop_ub; i31++) {
            y->data[y->size[0] * i31].re = b[0];
            y->data[y->size[0] * i31].im = 0.0;
        }

        for (k = 0; k < 6; k++) {
            i31 = y->size[0] * y->size[1];
            y->size[0] = 1;
            y->size[1] = s->size[1];
            emxEnsureCapacity_creal_T(y, i31);
            loop_ub = s->size[0] * s->size[1];
            for (i31 = 0; i31 < loop_ub; i31++) {
                s_re = s->data[i31].re * y->data[i31].re - s->data[i31].im * y->data[i31]
                       .im;
                s_im = s->data[i31].re * y->data[i31].im + s->data[i31].im * y->data[i31]
                       .re;
                y->data[i31].re = s_re + b[k + 1];
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
        s_re = digw->data[i31] * 0.0;
        s_im = digw->data[i31];
        s->data[i31].re = 6.0 * s_re;
        s->data[i31].im = 6.0 * s_im;
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
    emxArray_real_T *y;
    int b_idx_0;
    int i32;
    emxArray_real_T *r8;
    emxArray_real_T *digw;
    emxArray_creal_T *r9;
    emxArray_creal_T *r10;
    double b_y;
    double re;
    double im;
    emxInit_real_T(&y, 2);

    /*  Cast to enforce precision rules */
    /*  Remaining are default or for advanced use */
    /*  Make b a row */
    /* -------------------------------------------------------------------------- */
    b_idx_0 = b->size[1];
    i32 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = b_idx_0;
    emxEnsureCapacity_real_T(y, i32);
    for (i32 = 0; i32 < b_idx_0; i32++) {
        y->data[y->size[0] * i32] = b->data[i32];
    }

    emxInit_real_T(&r8, 2);

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
    i32 = r8->size[0] * r8->size[1];
    r8->size[0] = 1;
    r8->size[1] = w->size[1];
    emxEnsureCapacity_real_T(r8, i32);
    b_idx_0 = w->size[0] * w->size[1];
    for (i32 = 0; i32 < b_idx_0; i32++) {
        r8->data[i32] = 6.2831853071795862 * w->data[i32];
    }

    emxInit_real_T(&digw, 2);
    emxInit_creal_T(&r9, 2);
    rdivide(r8, Fs, digw);

    /*  Convert from Hz to rad/sample for computational purposes */
    /*  Digital frequency must be used for this calculation */
    i32 = r9->size[0] * r9->size[1];
    r9->size[0] = 1;
    r9->size[1] = digw->size[1];
    emxEnsureCapacity_creal_T(r9, i32);
    b_idx_0 = digw->size[0] * digw->size[1];
    emxFree_real_T(&r8);
    for (i32 = 0; i32 < b_idx_0; i32++) {
        r9->data[i32].re = digw->data[i32] * 0.0;
        r9->data[i32].im = digw->data[i32];
    }

    emxInit_creal_T(&r10, 2);
    c_exp(r9);
    j_polyval(y, r9, r10);
    i32 = r9->size[0] * r9->size[1];
    r9->size[0] = 1;
    r9->size[1] = digw->size[1];
    emxEnsureCapacity_creal_T(r9, i32);
    b_y = (double)y->size[1] - 1.0;
    b_idx_0 = digw->size[0] * digw->size[1];
    emxFree_real_T(&y);
    for (i32 = 0; i32 < b_idx_0; i32++) {
        re = digw->data[i32] * 0.0;
        im = digw->data[i32];
        r9->data[i32].re = b_y * re;
        r9->data[i32].im = b_y * im;
    }

    emxFree_real_T(&digw);
    c_exp(r9);
    b_rdivide(r10, r9, hh);

    /*  Generate the default structure to pass to freqzplot */
    /*  If rad/sample, Fs is empty */
    emxFree_creal_T(&r10);
    emxFree_creal_T(&r9);
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
    emxArray_real_T *r11;
    int i34;
    int loop_ub;
    emxArray_real_T *digw;
    emxArray_creal_T *s;
    emxArray_creal_T *y;
    boolean_T b5;
    int k;
    double s_re;
    double s_im;
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
    i34 = r11->size[0] * r11->size[1];
    r11->size[0] = 1;
    r11->size[1] = w->size[1];
    emxEnsureCapacity_real_T(r11, i34);
    loop_ub = w->size[0] * w->size[1];
    for (i34 = 0; i34 < loop_ub; i34++) {
        r11->data[i34] = 6.2831853071795862 * w->data[i34];
    }

    emxInit_real_T(&digw, 2);
    emxInit_creal_T(&s, 2);
    rdivide(r11, Fs, digw);

    /*  Convert from Hz to rad/sample for computational purposes */
    i34 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = digw->size[1];
    emxEnsureCapacity_creal_T(s, i34);
    loop_ub = digw->size[0] * digw->size[1];
    emxFree_real_T(&r11);
    for (i34 = 0; i34 < loop_ub; i34++) {
        s->data[i34].re = digw->data[i34] * 0.0;
        s->data[i34].im = digw->data[i34];
    }

    emxInit_creal_T(&y, 2);
    c_exp(s);

    /*  Digital frequency must be used for this calculation */
    i34 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = s->size[1];
    emxEnsureCapacity_creal_T(y, i34);
    b5 = (y->size[1] == 0);
    if (!b5) {
        i34 = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity_creal_T(y, i34);
        loop_ub = y->size[1];
        for (i34 = 0; i34 < loop_ub; i34++) {
            y->data[y->size[0] * i34].re = b[0];
            y->data[y->size[0] * i34].im = 0.0;
        }

        for (k = 0; k < 28; k++) {
            i34 = y->size[0] * y->size[1];
            y->size[0] = 1;
            y->size[1] = s->size[1];
            emxEnsureCapacity_creal_T(y, i34);
            loop_ub = s->size[0] * s->size[1];
            for (i34 = 0; i34 < loop_ub; i34++) {
                s_re = s->data[i34].re * y->data[i34].re - s->data[i34].im * y->data[i34]
                       .im;
                s_im = s->data[i34].re * y->data[i34].im + s->data[i34].im * y->data[i34]
                       .re;
                y->data[i34].re = s_re + b[k + 1];
                y->data[i34].im = s_im;
            }
        }
    }

    i34 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = digw->size[1];
    emxEnsureCapacity_creal_T(s, i34);
    loop_ub = digw->size[0] * digw->size[1];
    for (i34 = 0; i34 < loop_ub; i34++) {
        s_re = digw->data[i34] * 0.0;
        s_im = digw->data[i34];
        s->data[i34].re = 28.0 * s_re;
        s->data[i34].im = 28.0 * s_im;
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
 * Arguments    : const double b[13]
 *                const emxArray_real_T *w
 *                double Fs
 *                emxArray_creal_T *hh
 * Return Type  : void
 */
static void p_freqz_cg(const double b[13], const emxArray_real_T *w, double Fs,
                       emxArray_creal_T *hh)
{
    emxArray_real_T *r12;
    int i35;
    int loop_ub;
    emxArray_real_T *digw;
    emxArray_creal_T *s;
    emxArray_creal_T *y;
    boolean_T b6;
    int k;
    double s_re;
    double s_im;
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
    i35 = r12->size[0] * r12->size[1];
    r12->size[0] = 1;
    r12->size[1] = w->size[1];
    emxEnsureCapacity_real_T(r12, i35);
    loop_ub = w->size[0] * w->size[1];
    for (i35 = 0; i35 < loop_ub; i35++) {
        r12->data[i35] = 6.2831853071795862 * w->data[i35];
    }

    emxInit_real_T(&digw, 2);
    emxInit_creal_T(&s, 2);
    rdivide(r12, Fs, digw);

    /*  Convert from Hz to rad/sample for computational purposes */
    i35 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = digw->size[1];
    emxEnsureCapacity_creal_T(s, i35);
    loop_ub = digw->size[0] * digw->size[1];
    emxFree_real_T(&r12);
    for (i35 = 0; i35 < loop_ub; i35++) {
        s->data[i35].re = digw->data[i35] * 0.0;
        s->data[i35].im = digw->data[i35];
    }

    emxInit_creal_T(&y, 2);
    c_exp(s);

    /*  Digital frequency must be used for this calculation */
    i35 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = s->size[1];
    emxEnsureCapacity_creal_T(y, i35);
    b6 = (y->size[1] == 0);
    if (!b6) {
        i35 = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity_creal_T(y, i35);
        loop_ub = y->size[1];
        for (i35 = 0; i35 < loop_ub; i35++) {
            y->data[y->size[0] * i35].re = b[0];
            y->data[y->size[0] * i35].im = 0.0;
        }

        for (k = 0; k < 12; k++) {
            i35 = y->size[0] * y->size[1];
            y->size[0] = 1;
            y->size[1] = s->size[1];
            emxEnsureCapacity_creal_T(y, i35);
            loop_ub = s->size[0] * s->size[1];
            for (i35 = 0; i35 < loop_ub; i35++) {
                s_re = s->data[i35].re * y->data[i35].re - s->data[i35].im * y->data[i35]
                       .im;
                s_im = s->data[i35].re * y->data[i35].im + s->data[i35].im * y->data[i35]
                       .re;
                y->data[i35].re = s_re + b[k + 1];
                y->data[i35].im = s_im;
            }
        }
    }

    i35 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = digw->size[1];
    emxEnsureCapacity_creal_T(s, i35);
    loop_ub = digw->size[0] * digw->size[1];
    for (i35 = 0; i35 < loop_ub; i35++) {
        s_re = digw->data[i35] * 0.0;
        s_im = digw->data[i35];
        s->data[i35].re = 12.0 * s_re;
        s->data[i35].im = 12.0 * s_im;
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
    emxArray_creal_T *r2;
    emxInit_creal_T1(&r2, 1);
    eig(x, r2);
    vector_poly(r2, c);
    emxFree_creal_T(&r2);
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
    double x_im;
    for (i9 = 0; i9 < 2048; i9++) {
        y[i9].re = p[0];
        y[i9].im = 0.0;
    }

    for (k = 0; k < 14; k++) {
        for (i9 = 0; i9 < 2048; i9++) {
            x_im = x[i9].re * y[i9].im + x[i9].im * y[i9].re;
            y[i9].re = (x[i9].re * y[i9].re - x[i9].im * y[i9].im) + p[k + 1];
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
    emxArray_real_T *r13;
    int i36;
    int loop_ub;
    emxArray_real_T *digw;
    emxArray_creal_T *s;
    emxArray_creal_T *y;
    boolean_T b7;
    int k;
    double s_re;
    double s_im;
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
    i36 = r13->size[0] * r13->size[1];
    r13->size[0] = 1;
    r13->size[1] = w->size[1];
    emxEnsureCapacity_real_T(r13, i36);
    loop_ub = w->size[0] * w->size[1];
    for (i36 = 0; i36 < loop_ub; i36++) {
        r13->data[i36] = 6.2831853071795862 * w->data[i36];
    }

    emxInit_real_T(&digw, 2);
    emxInit_creal_T(&s, 2);
    rdivide(r13, Fs, digw);

    /*  Convert from Hz to rad/sample for computational purposes */
    i36 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = digw->size[1];
    emxEnsureCapacity_creal_T(s, i36);
    loop_ub = digw->size[0] * digw->size[1];
    emxFree_real_T(&r13);
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
    b7 = (y->size[1] == 0);
    if (!b7) {
        i36 = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity_creal_T(y, i36);
        loop_ub = y->size[1];
        for (i36 = 0; i36 < loop_ub; i36++) {
            y->data[y->size[0] * i36].re = b[0];
            y->data[y->size[0] * i36].im = 0.0;
        }

        for (k = 0; k < 56; k++) {
            i36 = y->size[0] * y->size[1];
            y->size[0] = 1;
            y->size[1] = s->size[1];
            emxEnsureCapacity_creal_T(y, i36);
            loop_ub = s->size[0] * s->size[1];
            for (i36 = 0; i36 < loop_ub; i36++) {
                s_re = s->data[i36].re * y->data[i36].re - s->data[i36].im * y->data[i36]
                       .im;
                s_im = s->data[i36].re * y->data[i36].im + s->data[i36].im * y->data[i36]
                       .re;
                y->data[i36].re = s_re + b[k + 1];
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
        s_re = digw->data[i36] * 0.0;
        s_im = digw->data[i36];
        s->data[i36].re = 56.0 * s_re;
        s->data[i36].im = 56.0 * s_im;
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
 * Arguments    : const double b[43]
 *                const emxArray_real_T *w
 *                double Fs
 *                emxArray_creal_T *hh
 * Return Type  : void
 */
static void r_freqz_cg(const double b[43], const emxArray_real_T *w, double Fs,
                       emxArray_creal_T *hh)
{
    emxArray_real_T *r14;
    int i37;
    int loop_ub;
    emxArray_real_T *digw;
    emxArray_creal_T *s;
    emxArray_creal_T *y;
    boolean_T b8;
    int k;
    double s_re;
    double s_im;
    emxInit_real_T(&r14, 2);

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
    i37 = r14->size[0] * r14->size[1];
    r14->size[0] = 1;
    r14->size[1] = w->size[1];
    emxEnsureCapacity_real_T(r14, i37);
    loop_ub = w->size[0] * w->size[1];
    for (i37 = 0; i37 < loop_ub; i37++) {
        r14->data[i37] = 6.2831853071795862 * w->data[i37];
    }

    emxInit_real_T(&digw, 2);
    emxInit_creal_T(&s, 2);
    rdivide(r14, Fs, digw);

    /*  Convert from Hz to rad/sample for computational purposes */
    i37 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = digw->size[1];
    emxEnsureCapacity_creal_T(s, i37);
    loop_ub = digw->size[0] * digw->size[1];
    emxFree_real_T(&r14);
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
    b8 = (y->size[1] == 0);
    if (!b8) {
        i37 = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity_creal_T(y, i37);
        loop_ub = y->size[1];
        for (i37 = 0; i37 < loop_ub; i37++) {
            y->data[y->size[0] * i37].re = b[0];
            y->data[y->size[0] * i37].im = 0.0;
        }

        for (k = 0; k < 42; k++) {
            i37 = y->size[0] * y->size[1];
            y->size[0] = 1;
            y->size[1] = s->size[1];
            emxEnsureCapacity_creal_T(y, i37);
            loop_ub = s->size[0] * s->size[1];
            for (i37 = 0; i37 < loop_ub; i37++) {
                s_re = s->data[i37].re * y->data[i37].re - s->data[i37].im * y->data[i37]
                       .im;
                s_im = s->data[i37].re * y->data[i37].im + s->data[i37].im * y->data[i37]
                       .re;
                y->data[i37].re = s_re + b[k + 1];
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
        s_re = digw->data[i37] * 0.0;
        s_im = digw->data[i37];
        s->data[i37].re = 42.0 * s_re;
        s->data[i37].im = 42.0 * s_im;
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
 * Arguments    : const emxArray_real_T *x
 *                double y
 *                emxArray_real_T *z
 * Return Type  : void
 */
static void rdivide(const emxArray_real_T *x, double y, emxArray_real_T *z)
{
    int i27;
    int loop_ub;
    i27 = z->size[0] * z->size[1];
    z->size[0] = 1;
    z->size[1] = x->size[1];
    emxEnsureCapacity_real_T(z, i27);
    loop_ub = x->size[0] * x->size[1];
    for (i27 = 0; i27 < loop_ub; i27++) {
        z->data[i27] = x->data[i27] / y;
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
    emxArray_int32_T *r23;
    int i51;
    int i;
    int end;
    double b_x;
    int vlen;

    /*  */
    /*    Author: T. Krauss 1993 */
    /*        Was Revision: 1.4, Date: 1994/01/25 17:59:44 */
    y = 1.0;
    l = 0;
    emxInit_real_T(&xx, 2);
    emxInit_int32_T(&r23, 2);
    while (l <= (int)m - 1) {
        if ((m == 0.0) || (((m > 0.0) && (1.0 + (double)l > n)) || ((0.0 > m) && (n >
                           1.0 + (double)l)))) {
            i51 = 1;
            i = 1;
            end = 0;
        } else {
            i51 = l + 1;
            i = (int)m;
            end = (int)n;
        }

        b_x = x->data[(int)k - 1];
        vlen = xx->size[0] * xx->size[1];
        xx->size[0] = 1;
        xx->size[1] = div_s32_floor(end - i51, i) + 1;
        emxEnsureCapacity_real_T(xx, vlen);
        vlen = div_s32_floor(end - i51, i);
        for (end = 0; end <= vlen; end++) {
            xx->data[xx->size[0] * end] = 2.0 * (b_x - x->data[(i51 + i * end) - 1]);
        }

        end = xx->size[1] - 1;
        vlen = 0;
        for (i = 0; i <= end; i++) {
            if (xx->data[i] != 0.0) {
                vlen++;
            }
        }

        i51 = r23->size[0] * r23->size[1];
        r23->size[0] = 1;
        r23->size[1] = vlen;
        emxEnsureCapacity_int32_T(r23, i51);
        vlen = 0;
        for (i = 0; i <= end; i++) {
            if (xx->data[i] != 0.0) {
                r23->data[vlen] = i + 1;
                vlen++;
            }
        }

        vlen = r23->size[1];
        if (r23->size[1] == 0) {
            b_x = 1.0;
        } else {
            b_x = xx->data[r23->data[0] - 1];
            for (i = 2; i <= vlen; i++) {
                b_x *= xx->data[r23->data[r23->size[0] * (i - 1)] - 1];
            }
        }

        y *= b_x;
        l++;
    }

    emxFree_int32_T(&r23);
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
    boolean_T b_valid;
    double nodd;
    double nfcns;
    int varargin_2;
    int ngrid;
    emxArray_real_T *j;
    emxArray_real_T *x2;
    int i44;
    double temp;
    int loop_ub;
    emxArray_real_T *iext;
    emxArray_real_T *y;
    emxArray_real_T *x;
    int ixstart;
    int flag;
    double comp;
    double dtemp;
    double b_y1;
    int luck;
    int nut1;
    double err;
    double b_dev;
    double devl;
    emxArray_real_T *ad;
    int niter;
    double jchnge;
    double d1;
    emxArray_real_T *l;
    emxArray_int32_T *r21;
    emxArray_real_T *b;
    emxArray_real_T *a;
    emxArray_real_T *b_wt;
    emxArray_int8_T *r22;
    emxArray_int32_T *b_iext;
    emxArray_real_T *mtmp;
    boolean_T guard1 = false;
    int exitg1;
    double k1;
    int nut;
    double b_j;
    double b_l;
    int nu;
    double varargin_1[2];
    double dnum;
    double kup;
    int i45;
    int i46;
    int i47;
    int i48;
    int i49;
    boolean_T guard2 = false;
    double b_y;
    int flag34;
    int exitg3;
    boolean_T exitg2;
    boolean_T exitg4;

    /*  */
    b_valid = true;
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
        i44 = x2->size[0] * x2->size[1];
        x2->size[0] = 1;
        x2->size[1] = grid->size[1];
        emxEnsureCapacity_real_T(x2, i44);
        loop_ub = grid->size[0] * grid->size[1];
        for (i44 = 0; i44 < loop_ub; i44++) {
            x2->data[i44] = 3.1415926535897931 * grid->data[i44];
        }

        b_cos(x2);
        c_rdivide(des, x2, j);
        i44 = des->size[0] * des->size[1];
        des->size[0] = 1;
        des->size[1] = j->size[1];
        emxEnsureCapacity_real_T(des, i44);
        loop_ub = j->size[0] * j->size[1];
        for (i44 = 0; i44 < loop_ub; i44++) {
            des->data[i44] = j->data[i44];
        }

        i44 = x2->size[0] * x2->size[1];
        x2->size[0] = 1;
        x2->size[1] = grid->size[1];
        emxEnsureCapacity_real_T(x2, i44);
        loop_ub = grid->size[0] * grid->size[1];
        for (i44 = 0; i44 < loop_ub; i44++) {
            x2->data[i44] = 3.1415926535897931 * grid->data[i44];
        }

        b_cos(x2);
        i44 = wt->size[0] * wt->size[1];
        wt->size[0] = 1;
        emxEnsureCapacity_real_T(wt, i44);
        ixstart = wt->size[0];
        flag = wt->size[1];
        loop_ub = ixstart * flag;
        for (i44 = 0; i44 < loop_ub; i44++) {
            wt->data[i44] *= x2->data[i44];
        }
    }

    temp = ((double)grid->size[1] - 1.0) / nfcns;
    if (rtIsNaN(nfcns)) {
        i44 = j->size[0] * j->size[1];
        j->size[0] = 1;
        j->size[1] = 1;
        emxEnsureCapacity_real_T(j, i44);
        j->data[0] = rtNaN;
    } else if (nfcns < 1.0) {
        i44 = j->size[0] * j->size[1];
        j->size[0] = 1;
        j->size[1] = 0;
        emxEnsureCapacity_real_T(j, i44);
    } else if (rtIsInf(nfcns) && (1.0 == nfcns)) {
        i44 = j->size[0] * j->size[1];
        j->size[0] = 1;
        j->size[1] = 1;
        emxEnsureCapacity_real_T(j, i44);
        j->data[0] = rtNaN;
    } else {
        i44 = j->size[0] * j->size[1];
        j->size[0] = 1;
        j->size[1] = (int)floor(nfcns - 1.0) + 1;
        emxEnsureCapacity_real_T(j, i44);
        loop_ub = (int)floor(nfcns - 1.0);
        for (i44 = 0; i44 <= loop_ub; i44++) {
            j->data[j->size[0] * i44] = 1.0 + (double)i44;
        }
    }

    i44 = x2->size[0] * x2->size[1];
    x2->size[0] = 1;
    x2->size[1] = j->size[1] + 1;
    emxEnsureCapacity_real_T(x2, i44);
    loop_ub = j->size[1];
    for (i44 = 0; i44 < loop_ub; i44++) {
        x2->data[x2->size[0] * i44] = temp * (j->data[j->size[0] * i44] - 1.0) + 1.0;
    }

    emxInit_real_T1(&iext, 1);
    x2->data[x2->size[0] * j->size[1]] = grid->size[1];
    c_fix(x2);
    i44 = iext->size[0];
    iext->size[0] = x2->size[1] + 1;
    emxEnsureCapacity_real_T1(iext, i44);
    loop_ub = x2->size[1];
    for (i44 = 0; i44 < loop_ub; i44++) {
        iext->data[i44] = x2->data[x2->size[0] * i44];
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
    i44 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = 1;
    emxEnsureCapacity_real_T(y, i44);
    y->data[0] = -1.0;
    b_dev = -1.0;
    devl = -1.0;
    i44 = x->size[0] * x->size[1];
    x->size[0] = 1;
    x->size[1] = (int)(nfcns + 1.0);
    emxEnsureCapacity_real_T(x, i44);
    loop_ub = (int)(nfcns + 1.0);
    for (i44 = 0; i44 < loop_ub; i44++) {
        x->data[i44] = 0.0;
    }

    emxInit_real_T(&ad, 2);
    niter = 0;
    jchnge = 1.0;
    d1 = (nfcns - 1.0) / 15.0;
    b_fix(&d1);
    i44 = ad->size[0] * ad->size[1];
    ad->size[0] = 1;
    ad->size[1] = (int)(nfcns + 1.0);
    emxEnsureCapacity_real_T(ad, i44);
    loop_ub = (int)(nfcns + 1.0);
    for (i44 = 0; i44 < loop_ub; i44++) {
        ad->data[i44] = 0.0;
    }

    /*  index manager(s) */
    emxInit_real_T(&l, 2);
    emxInit_int32_T(&r21, 2);
    emxInit_real_T1(&b, 1);
    emxInit_real_T(&a, 2);
    emxInit_real_T(&b_wt, 2);
    emxInit_int8_T(&r22, 2);
    emxInit_int32_T1(&b_iext, 1);
    emxInit_real_T1(&mtmp, 1);
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

                i44 = l->size[0] * l->size[1];
                l->size[0] = 1;
                l->size[1] = loop_ub;
                emxEnsureCapacity_real_T(l, i44);
                for (i44 = 0; i44 < loop_ub; i44++) {
                    l->data[l->size[0] * i44] = iext->data[i44];
                }

                i44 = x->size[0] * x->size[1];
                x->size[0] = 1;
                x->size[1] = l->size[1];
                emxEnsureCapacity_real_T(x, i44);
                nut = l->size[0] * l->size[1];
                for (i44 = 0; i44 < nut; i44++) {
                    x->data[i44] = 6.2831853071795862 * grid->data[(int)l->data[i44] - 1];
                }

                b_cos(x);
                for (ixstart = 0; ixstart < (int)(nfcns + 1.0); ixstart++) {
                    ad->data[ixstart] = remezdd(1.0 + (double)ixstart, nfcns + 1.0, d1 +
                                                1.0, x);
                }

                for (i44 = 0; i44 < 2; i44++) {
                    varargin_1[i44] = ad->size[i44];
                }

                i44 = j->size[0] * j->size[1];
                j->size[0] = 1;
                j->size[1] = (int)varargin_1[1];
                emxEnsureCapacity_real_T(j, i44);
                nut = (int)varargin_1[1];
                for (i44 = 0; i44 < nut; i44++) {
                    j->data[i44] = 1.0;
                }

                if (2.0 > nfcns + 1.0) {
                    i44 = 0;
                    i45 = 1;
                    i46 = 0;
                    i47 = 0;
                    i48 = 1;
                } else {
                    i44 = 1;
                    i45 = 2;
                    i46 = (int)(nfcns + 1.0);
                    i47 = 1;
                    i48 = 2;
                }

                i49 = r22->size[0] * r22->size[1];
                r22->size[0] = 1;
                r22->size[1] = (int)varargin_1[1];
                emxEnsureCapacity_int8_T(r22, i49);
                nut = (int)varargin_1[1];
                for (i49 = 0; i49 < nut; i49++) {
                    r22->data[r22->size[0] * i49] = 1;
                }

                nut = div_s32_floor((i46 - i44) - 1, i45);
                for (i46 = 0; i46 <= nut; i46++) {
                    j->data[i47 + i48 * i46] = -(double)r22->data[i44 + i45 * i46];
                }

                i44 = b->size[0];
                b->size[0] = l->size[1];
                emxEnsureCapacity_real_T1(b, i44);
                nut = l->size[1];
                for (i44 = 0; i44 < nut; i44++) {
                    b->data[i44] = des->data[(int)l->data[l->size[0] * i44] - 1];
                }

                guard2 = false;
                if (ad->size[1] == 1) {
                    guard2 = true;
                } else {
                    i44 = b_iext->size[0];
                    b_iext->size[0] = loop_ub;
                    emxEnsureCapacity_int32_T1(b_iext, i44);
                    for (i44 = 0; i44 < loop_ub; i44++) {
                        b_iext->data[i44] = (int)iext->data[i44];
                    }

                    if (b_iext->size[0] == 1) {
                        guard2 = true;
                    } else {
                        dnum = 0.0;
                        for (i44 = 0; i44 < ad->size[1]; i44++) {
                            dnum += ad->data[ad->size[0] * i44] * b->data[i44];
                        }
                    }
                }

                if (guard2) {
                    dnum = 0.0;
                    for (i44 = 0; i44 < ad->size[1]; i44++) {
                        dnum += ad->data[ad->size[0] * i44] * b->data[i44];
                    }
                }

                i44 = b_wt->size[0] * b_wt->size[1];
                b_wt->size[0] = 1;
                b_wt->size[1] = l->size[1];
                emxEnsureCapacity_real_T(b_wt, i44);
                loop_ub = l->size[0] * l->size[1];
                for (i44 = 0; i44 < loop_ub; i44++) {
                    b_wt->data[i44] = wt->data[(int)l->data[i44] - 1];
                }

                c_rdivide(ad, b_wt, x2);
                i44 = b->size[0];
                b->size[0] = x2->size[1];
                emxEnsureCapacity_real_T1(b, i44);
                loop_ub = x2->size[1];
                for (i44 = 0; i44 < loop_ub; i44++) {
                    b->data[i44] = x2->data[x2->size[0] * i44];
                }

                if ((j->size[1] == 1) || (b->size[0] == 1)) {
                    temp = 0.0;
                    for (i44 = 0; i44 < j->size[1]; i44++) {
                        temp += j->data[j->size[0] * i44] * b->data[i44];
                    }
                } else {
                    temp = 0.0;
                    for (i44 = 0; i44 < j->size[1]; i44++) {
                        temp += j->data[j->size[0] * i44] * b->data[i44];
                    }
                }

                b_dev = dnum / temp;
                nu = 1;
                if (b_dev > 0.0) {
                    nu = -1;
                }

                b_dev *= -(double)nu;
                temp = (double)nu * b_dev;
                i44 = a->size[0] * a->size[1];
                a->size[0] = 1;
                a->size[1] = j->size[1];
                emxEnsureCapacity_real_T(a, i44);
                loop_ub = j->size[0] * j->size[1];
                for (i44 = 0; i44 < loop_ub; i44++) {
                    a->data[i44] = temp * j->data[i44];
                }

                i44 = b_wt->size[0] * b_wt->size[1];
                b_wt->size[0] = 1;
                b_wt->size[1] = l->size[1];
                emxEnsureCapacity_real_T(b_wt, i44);
                loop_ub = l->size[0] * l->size[1];
                for (i44 = 0; i44 < loop_ub; i44++) {
                    b_wt->data[i44] = wt->data[(int)l->data[i44] - 1];
                }

                c_rdivide(a, b_wt, x2);
                i44 = y->size[0] * y->size[1];
                y->size[0] = 1;
                y->size[1] = l->size[1];
                emxEnsureCapacity_real_T(y, i44);
                loop_ub = l->size[0] * l->size[1];
                for (i44 = 0; i44 < loop_ub; i44++) {
                    y->data[i44] = des->data[(int)l->data[i44] - 1] + x2->data[i44];
                }

                if (b_dev <= devl) {
                    /* warning(message('signal:firpm:DidNotConverge',niter)) */
                    printf("%s\n", "DidNotConverge");
                    fflush(stdout);
                    i44 = h->size[0] * h->size[1];
                    h->size[0] = (int)nfilt;
                    h->size[1] = 1;
                    emxEnsureCapacity_real_T(h, i44);
                    loop_ub = (int)nfilt;
                    for (i44 = 0; i44 < loop_ub; i44++) {
                        h->data[i44] = 0.0;
                    }

                    b_dev = -1.0;

                    /* iext */
                    b_valid = false;
                    exitg1 = 1;
                } else {
                    devl = b_dev;
                    jchnge = 0.0;
                    k1 = iext->data[0];
                    dnum = iext->data[(int)(nfcns + 1.0) - 1];
                    temp = 0.0;
                    nut = -nu;
                    b_j = 1.0;
                    flag34 = 1;
                    while (b_j < (nfcns + 1.0) + 1.0) {
                        kup = iext->data[(int)(unsigned int)b_j];
                        b_l = iext->data[(int)b_j - 1] + 1.0;
                        nut = -nut;
                        if (b_j == 2.0) {
                            b_y1 = comp;
                        }

                        comp = b_dev;
                        flag = 1;
                        if (iext->data[(int)b_j - 1] + 1.0 < iext->data[(int)(b_j + 1.0) - 1]) {
                            /*  gee */
                            dtemp = cos(6.2831853071795862 * grid->data[(int)(iext->data[(int)
                                        b_j - 1] + 1.0) - 1]);
                            i44 = b_wt->size[0] * b_wt->size[1];
                            b_wt->size[0] = 1;
                            b_wt->size[1] = x->size[1];
                            emxEnsureCapacity_real_T(b_wt, i44);
                            loop_ub = x->size[0] * x->size[1];
                            for (i44 = 0; i44 < loop_ub; i44++) {
                                b_wt->data[i44] = dtemp - x->data[i44];
                            }

                            c_rdivide(ad, b_wt, j);
                            i44 = b->size[0];
                            b->size[0] = y->size[1];
                            emxEnsureCapacity_real_T1(b, i44);
                            loop_ub = y->size[1];
                            for (i44 = 0; i44 < loop_ub; i44++) {
                                b->data[i44] = y->data[y->size[0] * i44];
                            }

                            if ((j->size[1] == 1) || (b->size[0] == 1)) {
                                b_y = 0.0;
                                for (i44 = 0; i44 < j->size[1]; i44++) {
                                    b_y += j->data[j->size[0] * i44] * b->data[i44];
                                }
                            } else {
                                b_y = 0.0;
                                for (i44 = 0; i44 < j->size[1]; i44++) {
                                    b_y += j->data[j->size[0] * i44] * b->data[i44];
                                }
                            }

                            dtemp = b_sum(j);
                            err = (b_y / dtemp - des->data[(int)(iext->data[(int)b_j - 1] +
                                                                 1.0) - 1]) * wt->data[(int)(iext->data[(int)b_j - 1] + 1.0)
                                                                         - 1];
                            dtemp = (double)nut * err - b_dev;
                            if (dtemp > 0.0) {
                                comp = (double)nut * err;
                                b_l = (iext->data[(int)b_j - 1] + 1.0) + 1.0;
                                exitg2 = false;
                                while ((!exitg2) && (b_l < kup)) {
                                    /*  gee */
                                    dtemp = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                                    i44 = b_wt->size[0] * b_wt->size[1];
                                    b_wt->size[0] = 1;
                                    b_wt->size[1] = x->size[1];
                                    emxEnsureCapacity_real_T(b_wt, i44);
                                    loop_ub = x->size[0] * x->size[1];
                                    for (i44 = 0; i44 < loop_ub; i44++) {
                                        b_wt->data[i44] = dtemp - x->data[i44];
                                    }

                                    c_rdivide(ad, b_wt, j);
                                    i44 = b->size[0];
                                    b->size[0] = y->size[1];
                                    emxEnsureCapacity_real_T1(b, i44);
                                    loop_ub = y->size[1];
                                    for (i44 = 0; i44 < loop_ub; i44++) {
                                        b->data[i44] = y->data[y->size[0] * i44];
                                    }

                                    if ((j->size[1] == 1) || (b->size[0] == 1)) {
                                        b_y = 0.0;
                                        for (i44 = 0; i44 < j->size[1]; i44++) {
                                            b_y += j->data[j->size[0] * i44] * b->data[i44];
                                        }
                                    } else {
                                        b_y = 0.0;
                                        for (i44 = 0; i44 < j->size[1]; i44++) {
                                            b_y += j->data[j->size[0] * i44] * b->data[i44];
                                        }
                                    }

                                    dtemp = b_sum(j);
                                    err = (b_y / dtemp - des->data[(int)b_l - 1]) * wt->data[(int)
                                            b_l - 1];
                                    dtemp = (double)nut * err - comp;
                                    if (dtemp > 0.0) {
                                        comp = (double)nut * err;
                                        b_l++;
                                    } else {
                                        exitg2 = true;
                                    }
                                }

                                iext->data[(int)b_j - 1] = b_l - 1.0;
                                b_j++;
                                temp = b_l - 1.0;
                                jchnge++;
                                flag = 0;
                            }
                        }

                        if (flag != 0) {
                            b_l -= 2.0;
                            exitg2 = false;
                            while ((!exitg2) && (b_l > temp)) {
                                /*  gee */
                                dtemp = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                                i44 = b_wt->size[0] * b_wt->size[1];
                                b_wt->size[0] = 1;
                                b_wt->size[1] = x->size[1];
                                emxEnsureCapacity_real_T(b_wt, i44);
                                loop_ub = x->size[0] * x->size[1];
                                for (i44 = 0; i44 < loop_ub; i44++) {
                                    b_wt->data[i44] = dtemp - x->data[i44];
                                }

                                c_rdivide(ad, b_wt, j);
                                i44 = b->size[0];
                                b->size[0] = y->size[1];
                                emxEnsureCapacity_real_T1(b, i44);
                                loop_ub = y->size[1];
                                for (i44 = 0; i44 < loop_ub; i44++) {
                                    b->data[i44] = y->data[y->size[0] * i44];
                                }

                                if ((j->size[1] == 1) || (b->size[0] == 1)) {
                                    b_y = 0.0;
                                    for (i44 = 0; i44 < j->size[1]; i44++) {
                                        b_y += j->data[j->size[0] * i44] * b->data[i44];
                                    }
                                } else {
                                    b_y = 0.0;
                                    for (i44 = 0; i44 < j->size[1]; i44++) {
                                        b_y += j->data[j->size[0] * i44] * b->data[i44];
                                    }
                                }

                                dtemp = b_sum(j);
                                err = (b_y / dtemp - des->data[(int)b_l - 1]) * wt->data[(int)
                                        b_l - 1];
                                dtemp = (double)nut * err - comp;
                                if ((dtemp > 0.0) || (jchnge > 0.0)) {
                                    exitg2 = true;
                                } else {
                                    b_l--;
                                }
                            }

                            if (b_l <= temp) {
                                b_l = iext->data[(int)b_j - 1] + 1.0;
                                if (jchnge > 0.0) {
                                    iext->data[(int)b_j - 1] = (iext->data[(int)b_j - 1] + 1.0) -
                                                               1.0;
                                    b_j++;
                                    temp = b_l - 1.0;
                                    jchnge++;
                                } else {
                                    b_l = (iext->data[(int)b_j - 1] + 1.0) + 1.0;
                                    exitg2 = false;
                                    while ((!exitg2) && (b_l < kup)) {
                                        /*  gee */
                                        dtemp = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                                        i44 = b_wt->size[0] * b_wt->size[1];
                                        b_wt->size[0] = 1;
                                        b_wt->size[1] = x->size[1];
                                        emxEnsureCapacity_real_T(b_wt, i44);
                                        loop_ub = x->size[0] * x->size[1];
                                        for (i44 = 0; i44 < loop_ub; i44++) {
                                            b_wt->data[i44] = dtemp - x->data[i44];
                                        }

                                        c_rdivide(ad, b_wt, j);
                                        i44 = b->size[0];
                                        b->size[0] = y->size[1];
                                        emxEnsureCapacity_real_T1(b, i44);
                                        loop_ub = y->size[1];
                                        for (i44 = 0; i44 < loop_ub; i44++) {
                                            b->data[i44] = y->data[y->size[0] * i44];
                                        }

                                        if ((j->size[1] == 1) || (b->size[0] == 1)) {
                                            b_y = 0.0;
                                            for (i44 = 0; i44 < j->size[1]; i44++) {
                                                b_y += j->data[j->size[0] * i44] * b->data[i44];
                                            }
                                        } else {
                                            b_y = 0.0;
                                            for (i44 = 0; i44 < j->size[1]; i44++) {
                                                b_y += j->data[j->size[0] * i44] * b->data[i44];
                                            }
                                        }

                                        dtemp = b_sum(j);
                                        err = (b_y / dtemp - des->data[(int)b_l - 1]) * wt->data
                                              [(int)b_l - 1];
                                        dtemp = (double)nut * err - comp;
                                        if (dtemp > 0.0) {
                                            exitg2 = true;
                                        } else {
                                            b_l++;
                                        }
                                    }

                                    if ((b_l < kup) && (dtemp > 0.0)) {
                                        comp = (double)nut * err;
                                        b_l++;
                                        exitg2 = false;
                                        while ((!exitg2) && (b_l < kup)) {
                                            /*  gee */
                                            dtemp = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                                            i44 = b_wt->size[0] * b_wt->size[1];
                                            b_wt->size[0] = 1;
                                            b_wt->size[1] = x->size[1];
                                            emxEnsureCapacity_real_T(b_wt, i44);
                                            loop_ub = x->size[0] * x->size[1];
                                            for (i44 = 0; i44 < loop_ub; i44++) {
                                                b_wt->data[i44] = dtemp - x->data[i44];
                                            }

                                            c_rdivide(ad, b_wt, j);
                                            i44 = b->size[0];
                                            b->size[0] = y->size[1];
                                            emxEnsureCapacity_real_T1(b, i44);
                                            loop_ub = y->size[1];
                                            for (i44 = 0; i44 < loop_ub; i44++) {
                                                b->data[i44] = y->data[y->size[0] * i44];
                                            }

                                            if ((j->size[1] == 1) || (b->size[0] == 1)) {
                                                b_y = 0.0;
                                                for (i44 = 0; i44 < j->size[1]; i44++) {
                                                    b_y += j->data[j->size[0] * i44] * b->data[i44];
                                                }
                                            } else {
                                                b_y = 0.0;
                                                for (i44 = 0; i44 < j->size[1]; i44++) {
                                                    b_y += j->data[j->size[0] * i44] * b->data[i44];
                                                }
                                            }

                                            dtemp = b_sum(j);
                                            err = (b_y / dtemp - des->data[(int)b_l - 1]) * wt->data
                                                  [(int)b_l - 1];
                                            dtemp = (double)nut * err - comp;
                                            if (dtemp > 0.0) {
                                                comp = (double)nut * err;
                                                b_l++;
                                            } else {
                                                exitg2 = true;
                                            }
                                        }

                                        iext->data[(int)b_j - 1] = b_l - 1.0;
                                        b_j++;
                                        temp = b_l - 1.0;
                                        jchnge = 1.0;
                                    } else {
                                        temp = iext->data[(int)b_j - 1];
                                        b_j++;
                                    }
                                }
                            } else if (dtemp > 0.0) {
                                comp = (double)nut * err;
                                b_l--;
                                exitg2 = false;
                                while ((!exitg2) && (b_l > temp)) {
                                    /*  gee */
                                    dtemp = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                                    i44 = b_wt->size[0] * b_wt->size[1];
                                    b_wt->size[0] = 1;
                                    b_wt->size[1] = x->size[1];
                                    emxEnsureCapacity_real_T(b_wt, i44);
                                    loop_ub = x->size[0] * x->size[1];
                                    for (i44 = 0; i44 < loop_ub; i44++) {
                                        b_wt->data[i44] = dtemp - x->data[i44];
                                    }

                                    c_rdivide(ad, b_wt, j);
                                    i44 = b->size[0];
                                    b->size[0] = y->size[1];
                                    emxEnsureCapacity_real_T1(b, i44);
                                    loop_ub = y->size[1];
                                    for (i44 = 0; i44 < loop_ub; i44++) {
                                        b->data[i44] = y->data[y->size[0] * i44];
                                    }

                                    if ((j->size[1] == 1) || (b->size[0] == 1)) {
                                        b_y = 0.0;
                                        for (i44 = 0; i44 < j->size[1]; i44++) {
                                            b_y += j->data[j->size[0] * i44] * b->data[i44];
                                        }
                                    } else {
                                        b_y = 0.0;
                                        for (i44 = 0; i44 < j->size[1]; i44++) {
                                            b_y += j->data[j->size[0] * i44] * b->data[i44];
                                        }
                                    }

                                    dtemp = b_sum(j);
                                    err = (b_y / dtemp - des->data[(int)b_l - 1]) * wt->data[(int)
                                            b_l - 1];
                                    dtemp = (double)nut * err - comp;
                                    if (dtemp > 0.0) {
                                        comp = (double)nut * err;
                                        b_l--;
                                    } else {
                                        exitg2 = true;
                                    }
                                }

                                temp = iext->data[(int)b_j - 1];
                                iext->data[(int)b_j - 1] = b_l + 1.0;
                                b_j++;
                                jchnge++;
                            } else {
                                temp = iext->data[(int)b_j - 1];
                                b_j++;
                            }
                        }
                    }

                    do {
                        exitg3 = 0;
                        if (b_j == (nfcns + 1.0) + 1.0) {
                            varargin_1[1] = iext->data[0];
                            ixstart = 1;
                            if (rtIsNaN(k1)) {
                                flag = 2;
                                exitg2 = false;
                                while ((!exitg2) && (flag < 3)) {
                                    ixstart = 2;
                                    if (!rtIsNaN(varargin_1[1])) {
                                        k1 = varargin_1[1];
                                        exitg2 = true;
                                    } else {
                                        flag = 3;
                                    }
                                }
                            }

                            if ((ixstart < 2) && (varargin_1[1] < k1)) {
                                k1 = varargin_1[1];
                            }

                            varargin_1[1] = iext->data[(int)(nfcns + 1.0) - 1];
                            ixstart = 1;
                            if (rtIsNaN(dnum)) {
                                flag = 2;
                                exitg2 = false;
                                while ((!exitg2) && (flag < 3)) {
                                    ixstart = 2;
                                    if (!rtIsNaN(varargin_1[1])) {
                                        dnum = varargin_1[1];
                                        exitg2 = true;
                                    } else {
                                        flag = 3;
                                    }
                                }
                            }

                            if ((ixstart < 2) && (varargin_1[1] > dnum)) {
                                dnum = varargin_1[1];
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
                                dtemp = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                                i44 = b_wt->size[0] * b_wt->size[1];
                                b_wt->size[0] = 1;
                                b_wt->size[1] = x->size[1];
                                emxEnsureCapacity_real_T(b_wt, i44);
                                loop_ub = x->size[0] * x->size[1];
                                for (i44 = 0; i44 < loop_ub; i44++) {
                                    b_wt->data[i44] = dtemp - x->data[i44];
                                }

                                c_rdivide(ad, b_wt, j);
                                i44 = b->size[0];
                                b->size[0] = y->size[1];
                                emxEnsureCapacity_real_T1(b, i44);
                                loop_ub = y->size[1];
                                for (i44 = 0; i44 < loop_ub; i44++) {
                                    b->data[i44] = y->data[y->size[0] * i44];
                                }

                                if ((j->size[1] == 1) || (b->size[0] == 1)) {
                                    b_y = 0.0;
                                    for (i44 = 0; i44 < j->size[1]; i44++) {
                                        b_y += j->data[j->size[0] * i44] * b->data[i44];
                                    }
                                } else {
                                    b_y = 0.0;
                                    for (i44 = 0; i44 < j->size[1]; i44++) {
                                        b_y += j->data[j->size[0] * i44] * b->data[i44];
                                    }
                                }

                                dtemp = b_sum(j);
                                err = (b_y / dtemp - des->data[(int)b_l - 1]) * wt->data[(int)
                                        b_l - 1];
                                dtemp = err * -(double)nu - comp;
                                if (dtemp > 0.0) {
                                    comp = -(double)nu * err;
                                    b_l++;
                                    exitg4 = false;
                                    while ((!exitg4) && (b_l < k1)) {
                                        /*  gee */
                                        dtemp = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                                        i44 = b_wt->size[0] * b_wt->size[1];
                                        b_wt->size[0] = 1;
                                        b_wt->size[1] = x->size[1];
                                        emxEnsureCapacity_real_T(b_wt, i44);
                                        loop_ub = x->size[0] * x->size[1];
                                        for (i44 = 0; i44 < loop_ub; i44++) {
                                            b_wt->data[i44] = dtemp - x->data[i44];
                                        }

                                        c_rdivide(ad, b_wt, j);
                                        i44 = b->size[0];
                                        b->size[0] = y->size[1];
                                        emxEnsureCapacity_real_T1(b, i44);
                                        loop_ub = y->size[1];
                                        for (i44 = 0; i44 < loop_ub; i44++) {
                                            b->data[i44] = y->data[y->size[0] * i44];
                                        }

                                        if ((j->size[1] == 1) || (b->size[0] == 1)) {
                                            b_y = 0.0;
                                            for (i44 = 0; i44 < j->size[1]; i44++) {
                                                b_y += j->data[j->size[0] * i44] * b->data[i44];
                                            }
                                        } else {
                                            b_y = 0.0;
                                            for (i44 = 0; i44 < j->size[1]; i44++) {
                                                b_y += j->data[j->size[0] * i44] * b->data[i44];
                                            }
                                        }

                                        dtemp = b_sum(j);
                                        err = (b_y / dtemp - des->data[(int)b_l - 1]) * wt->data
                                              [(int)b_l - 1];
                                        dtemp = -(double)nu * err - comp;
                                        if (dtemp > 0.0) {
                                            comp = -(double)nu * err;
                                            b_l++;
                                        } else {
                                            exitg4 = true;
                                        }
                                    }

                                    iext->data[(int)((nfcns + 1.0) + 1.0) - 1] = b_l - 1.0;
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
                                while ((!exitg2) && (b_l > dnum)) {
                                    /*  gee */
                                    dtemp = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                                    i44 = b_wt->size[0] * b_wt->size[1];
                                    b_wt->size[0] = 1;
                                    b_wt->size[1] = x->size[1];
                                    emxEnsureCapacity_real_T(b_wt, i44);
                                    loop_ub = x->size[0] * x->size[1];
                                    for (i44 = 0; i44 < loop_ub; i44++) {
                                        b_wt->data[i44] = dtemp - x->data[i44];
                                    }

                                    c_rdivide(ad, b_wt, j);
                                    i44 = b->size[0];
                                    b->size[0] = y->size[1];
                                    emxEnsureCapacity_real_T1(b, i44);
                                    loop_ub = y->size[1];
                                    for (i44 = 0; i44 < loop_ub; i44++) {
                                        b->data[i44] = y->data[y->size[0] * i44];
                                    }

                                    if ((j->size[1] == 1) || (b->size[0] == 1)) {
                                        b_y = 0.0;
                                        for (i44 = 0; i44 < j->size[1]; i44++) {
                                            b_y += j->data[j->size[0] * i44] * b->data[i44];
                                        }
                                    } else {
                                        b_y = 0.0;
                                        for (i44 = 0; i44 < j->size[1]; i44++) {
                                            b_y += j->data[j->size[0] * i44] * b->data[i44];
                                        }
                                    }

                                    dtemp = b_sum(j);
                                    err = (b_y / dtemp - des->data[(int)b_l - 1]) * wt->data[(int)
                                            b_l - 1];
                                    dtemp = err * -(double)nut1 - comp;
                                    if (dtemp > 0.0) {
                                        comp = -(double)nut1 * err;
                                        luck = 16;
                                        b_l--;
                                        exitg4 = false;
                                        while ((!exitg4) && (b_l > dnum)) {
                                            /*  gee */
                                            dtemp = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                                            i44 = b_wt->size[0] * b_wt->size[1];
                                            b_wt->size[0] = 1;
                                            b_wt->size[1] = x->size[1];
                                            emxEnsureCapacity_real_T(b_wt, i44);
                                            loop_ub = x->size[0] * x->size[1];
                                            for (i44 = 0; i44 < loop_ub; i44++) {
                                                b_wt->data[i44] = dtemp - x->data[i44];
                                            }

                                            c_rdivide(ad, b_wt, j);
                                            i44 = b->size[0];
                                            b->size[0] = y->size[1];
                                            emxEnsureCapacity_real_T1(b, i44);
                                            loop_ub = y->size[1];
                                            for (i44 = 0; i44 < loop_ub; i44++) {
                                                b->data[i44] = y->data[y->size[0] * i44];
                                            }

                                            if ((j->size[1] == 1) || (b->size[0] == 1)) {
                                                b_y = 0.0;
                                                for (i44 = 0; i44 < j->size[1]; i44++) {
                                                    b_y += j->data[j->size[0] * i44] * b->data[i44];
                                                }
                                            } else {
                                                b_y = 0.0;
                                                for (i44 = 0; i44 < j->size[1]; i44++) {
                                                    b_y += j->data[j->size[0] * i44] * b->data[i44];
                                                }
                                            }

                                            dtemp = b_sum(j);
                                            err = (b_y / dtemp - des->data[(int)b_l - 1]) * wt->data
                                                  [(int)b_l - 1];
                                            dtemp = -(double)nut1 * err - comp;
                                            if (dtemp > 0.0) {
                                                comp = -(double)nut1 * err;
                                                b_l--;
                                            } else {
                                                exitg4 = true;
                                            }
                                        }

                                        iext->data[(int)((nfcns + 1.0) + 1.0) - 1] = b_l + 1.0;
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
                                        temp = (nfcns + 1.0) - nfcns;
                                        if (2.0 > temp) {
                                            i44 = -2;
                                            i45 = 0;
                                        } else {
                                            i44 = -1;
                                            i45 = (int)temp;
                                        }

                                        temp = (nfcns + 1.0) - nfcns;
                                        if (temp > (nfcns + 1.0) - 1.0) {
                                            i46 = 1;
                                            i47 = 0;
                                        } else {
                                            i46 = (int)temp;
                                            i47 = (int)((nfcns + 1.0) - 1.0);
                                        }

                                        /*  Update index */
                                        temp = (nfcns + 1.0) - nfcns;
                                        if (2.0 > temp) {
                                            i48 = -2;
                                            i49 = 0;
                                        } else {
                                            i48 = -1;
                                            i49 = (int)temp;
                                        }

                                        temp = (nfcns + 1.0) - nfcns;
                                        if (temp > (nfcns + 1.0) - 1.0) {
                                            ixstart = 1;
                                            flag = 0;
                                        } else {
                                            ixstart = (int)temp;
                                            flag = (int)((nfcns + 1.0) - 1.0);
                                        }

                                        nu = b_wt->size[0] * b_wt->size[1];
                                        b_wt->size[0] = 1;
                                        b_wt->size[1] = ((i45 - i44) + i47) - i46;
                                        emxEnsureCapacity_real_T(b_wt, nu);
                                        loop_ub = i45 - i44;
                                        for (nu = 0; nu <= loop_ub - 3; nu++) {
                                            b_wt->data[b_wt->size[0] * (nu + 1)] = iext->data[(i44 +
                                                                                   nu) + 2];
                                        }

                                        loop_ub = i47 - i46;
                                        for (i47 = 0; i47 <= loop_ub; i47++) {
                                            b_wt->data[b_wt->size[0] * (((i47 + i45) - i44) - 1)] =
                                                iext->data[(i46 + i47) - 1];
                                        }

                                        loop_ub = b_wt->size[1];
                                        i44 = r21->size[0] * r21->size[1];
                                        r21->size[0] = 1;
                                        r21->size[1] = loop_ub;
                                        emxEnsureCapacity_int32_T(r21, i44);
                                        for (i44 = 0; i44 < loop_ub; i44++) {
                                            r21->data[r21->size[0] * i44] = i44;
                                        }

                                        i44 = mtmp->size[0];
                                        mtmp->size[0] = ((i49 - i48) + flag) - ixstart;
                                        emxEnsureCapacity_real_T1(mtmp, i44);
                                        mtmp->data[0] = k1;
                                        loop_ub = i49 - i48;
                                        for (i44 = 0; i44 <= loop_ub - 3; i44++) {
                                            mtmp->data[i44 + 1] = iext->data[(i48 + i44) + 2];
                                        }

                                        loop_ub = flag - ixstart;
                                        for (i44 = 0; i44 <= loop_ub; i44++) {
                                            mtmp->data[((i44 + i49) - i48) - 1] = iext->data[(ixstart
                                                                                  + i44) - 1];
                                        }

                                        loop_ub = r21->size[1];
                                        for (i44 = 0; i44 < loop_ub; i44++) {
                                            iext->data[r21->data[r21->size[0] * i44]] = mtmp->data
                                                    [(*(int (*)[2])r21->size)[0] * i44];
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
                                i44 = 0;
                                i45 = 0;
                            } else {
                                i44 = 1;
                                i45 = (int)(nfcns + 1.0);
                            }

                            if (nfcns + 1.0 > (nfcns + 1.0) - 1.0) {
                                i46 = 0;
                                i47 = 0;
                            } else {
                                i46 = (int)(nfcns + 1.0) - 1;
                                i47 = (int)((nfcns + 1.0) - 1.0);
                            }

                            /*  Update index */
                            if (2.0 > nfcns + 1.0) {
                                i48 = 0;
                                i49 = 0;
                            } else {
                                i48 = 1;
                                i49 = (int)(nfcns + 1.0);
                            }

                            if (nfcns + 1.0 > (nfcns + 1.0) - 1.0) {
                                ixstart = 0;
                                flag = 0;
                            } else {
                                ixstart = (int)(nfcns + 1.0) - 1;
                                flag = (int)((nfcns + 1.0) - 1.0);
                            }

                            nu = b_wt->size[0] * b_wt->size[1];
                            b_wt->size[0] = 1;
                            b_wt->size[1] = (((i45 - i44) + i47) - i46) + 2;
                            emxEnsureCapacity_real_T(b_wt, nu);
                            loop_ub = i45 - i44;
                            for (nu = 0; nu < loop_ub; nu++) {
                                b_wt->data[b_wt->size[0] * nu] = iext->data[i44 + nu];
                            }

                            loop_ub = i47 - i46;
                            for (nu = 0; nu < loop_ub; nu++) {
                                b_wt->data[b_wt->size[0] * ((nu + i45) - i44)] = iext->data[i46
                                        + nu];
                            }

                            b_wt->data[b_wt->size[0] * (((i45 - i44) + i47) - i46)] =
                                iext->data[(int)((nfcns + 1.0) + 1.0) - 1];
                            b_wt->data[b_wt->size[0] * ((((i45 - i44) + i47) - i46) + 1)] =
                                iext->data[(int)((nfcns + 1.0) + 1.0) - 1];
                            loop_ub = b_wt->size[1];
                            i44 = r21->size[0] * r21->size[1];
                            r21->size[0] = 1;
                            r21->size[1] = loop_ub;
                            emxEnsureCapacity_int32_T(r21, i44);
                            for (i44 = 0; i44 < loop_ub; i44++) {
                                r21->data[r21->size[0] * i44] = i44;
                            }

                            temp = iext->data[(int)((nfcns + 1.0) + 1.0) - 1];
                            dnum = iext->data[(int)((nfcns + 1.0) + 1.0) - 1];
                            i44 = mtmp->size[0];
                            mtmp->size[0] = (((i49 - i48) + flag) - ixstart) + 2;
                            emxEnsureCapacity_real_T1(mtmp, i44);
                            loop_ub = i49 - i48;
                            for (i44 = 0; i44 < loop_ub; i44++) {
                                mtmp->data[i44] = iext->data[i48 + i44];
                            }

                            loop_ub = flag - ixstart;
                            for (i44 = 0; i44 < loop_ub; i44++) {
                                mtmp->data[(i44 + i49) - i48] = iext->data[ixstart + i44];
                            }

                            mtmp->data[((i49 - i48) + flag) - ixstart] = temp;
                            mtmp->data[(((i49 - i48) + flag) - ixstart) + 1] = dnum;
                            loop_ub = r21->size[1];
                            for (i44 = 0; i44 < loop_ub; i44++) {
                                iext->data[r21->data[r21->size[0] * i44]] = mtmp->data[(*(int (*)
                                        [2])r21->size)[0] * i44];
                            }

                            jchnge++;
                        } else {
                            ixstart = 1;
                            if (rtIsNaN(b_y1)) {
                                flag = 2;
                                exitg2 = false;
                                while ((!exitg2) && (flag < 3)) {
                                    ixstart = 2;
                                    if (!rtIsNaN(comp)) {
                                        b_y1 = comp;
                                        exitg2 = true;
                                    } else {
                                        flag = 3;
                                    }
                                }
                            }

                            if ((ixstart < 2) && (comp > b_y1)) {
                                b_y1 = comp;
                            }

                            k1 = iext->data[(int)((nfcns + 1.0) + 1.0) - 1];
                            comp = b_y1 * 1.00001;
                            b_l = ((double)ngrid + 1.0) - 1.0;
                            exitg2 = false;
                            while ((!exitg2) && (b_l > dnum)) {
                                /*  gee */
                                dtemp = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                                i44 = b_wt->size[0] * b_wt->size[1];
                                b_wt->size[0] = 1;
                                b_wt->size[1] = x->size[1];
                                emxEnsureCapacity_real_T(b_wt, i44);
                                loop_ub = x->size[0] * x->size[1];
                                for (i44 = 0; i44 < loop_ub; i44++) {
                                    b_wt->data[i44] = dtemp - x->data[i44];
                                }

                                c_rdivide(ad, b_wt, j);
                                i44 = b->size[0];
                                b->size[0] = y->size[1];
                                emxEnsureCapacity_real_T1(b, i44);
                                loop_ub = y->size[1];
                                for (i44 = 0; i44 < loop_ub; i44++) {
                                    b->data[i44] = y->data[y->size[0] * i44];
                                }

                                if ((j->size[1] == 1) || (b->size[0] == 1)) {
                                    b_y = 0.0;
                                    for (i44 = 0; i44 < j->size[1]; i44++) {
                                        b_y += j->data[j->size[0] * i44] * b->data[i44];
                                    }
                                } else {
                                    b_y = 0.0;
                                    for (i44 = 0; i44 < j->size[1]; i44++) {
                                        b_y += j->data[j->size[0] * i44] * b->data[i44];
                                    }
                                }

                                dtemp = b_sum(j);
                                err = (b_y / dtemp - des->data[(int)b_l - 1]) * wt->data[(int)
                                        b_l - 1];
                                dtemp = err * -(double)nut1 - comp;
                                if (dtemp > 0.0) {
                                    comp = -(double)nut1 * err;
                                    luck += 10;
                                    b_l--;
                                    exitg4 = false;
                                    while ((!exitg4) && (b_l > dnum)) {
                                        /*  gee */
                                        dtemp = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                                        i44 = b_wt->size[0] * b_wt->size[1];
                                        b_wt->size[0] = 1;
                                        b_wt->size[1] = x->size[1];
                                        emxEnsureCapacity_real_T(b_wt, i44);
                                        loop_ub = x->size[0] * x->size[1];
                                        for (i44 = 0; i44 < loop_ub; i44++) {
                                            b_wt->data[i44] = dtemp - x->data[i44];
                                        }

                                        c_rdivide(ad, b_wt, j);
                                        i44 = b->size[0];
                                        b->size[0] = y->size[1];
                                        emxEnsureCapacity_real_T1(b, i44);
                                        loop_ub = y->size[1];
                                        for (i44 = 0; i44 < loop_ub; i44++) {
                                            b->data[i44] = y->data[y->size[0] * i44];
                                        }

                                        if ((j->size[1] == 1) || (b->size[0] == 1)) {
                                            b_y = 0.0;
                                            for (i44 = 0; i44 < j->size[1]; i44++) {
                                                b_y += j->data[j->size[0] * i44] * b->data[i44];
                                            }
                                        } else {
                                            b_y = 0.0;
                                            for (i44 = 0; i44 < j->size[1]; i44++) {
                                                b_y += j->data[j->size[0] * i44] * b->data[i44];
                                            }
                                        }

                                        dtemp = b_sum(j);
                                        err = (b_y / dtemp - des->data[(int)b_l - 1]) * wt->data
                                              [(int)b_l - 1];
                                        dtemp = -(double)nut1 * err - comp;
                                        if (dtemp > 0.0) {
                                            comp = -(double)nut1 * err;
                                            b_l--;
                                        } else {
                                            exitg4 = true;
                                        }
                                    }

                                    iext->data[(int)((nfcns + 1.0) + 1.0) - 1] = b_l + 1.0;
                                    jchnge++;
                                    if (2.0 > nfcns + 1.0) {
                                        i44 = 0;
                                        i45 = 0;
                                    } else {
                                        i44 = 1;
                                        i45 = (int)(nfcns + 1.0);
                                    }

                                    if (nfcns + 1.0 > (nfcns + 1.0) - 1.0) {
                                        i46 = 0;
                                        i47 = 0;
                                    } else {
                                        i46 = (int)(nfcns + 1.0) - 1;
                                        i47 = (int)((nfcns + 1.0) - 1.0);
                                    }

                                    /*  Update index */
                                    if (2.0 > nfcns + 1.0) {
                                        i48 = 0;
                                        i49 = 0;
                                    } else {
                                        i48 = 1;
                                        i49 = (int)(nfcns + 1.0);
                                    }

                                    if (nfcns + 1.0 > (nfcns + 1.0) - 1.0) {
                                        ixstart = 0;
                                        flag = 0;
                                    } else {
                                        ixstart = (int)(nfcns + 1.0) - 1;
                                        flag = (int)((nfcns + 1.0) - 1.0);
                                    }

                                    nu = b_wt->size[0] * b_wt->size[1];
                                    b_wt->size[0] = 1;
                                    b_wt->size[1] = (((i45 - i44) + i47) - i46) + 2;
                                    emxEnsureCapacity_real_T(b_wt, nu);
                                    loop_ub = i45 - i44;
                                    for (nu = 0; nu < loop_ub; nu++) {
                                        b_wt->data[b_wt->size[0] * nu] = iext->data[i44 + nu];
                                    }

                                    loop_ub = i47 - i46;
                                    for (nu = 0; nu < loop_ub; nu++) {
                                        b_wt->data[b_wt->size[0] * ((nu + i45) - i44)] = iext->
                                                data[i46 + nu];
                                    }

                                    b_wt->data[b_wt->size[0] * (((i45 - i44) + i47) - i46)] =
                                        iext->data[(int)((nfcns + 1.0) + 1.0) - 1];
                                    b_wt->data[b_wt->size[0] * ((((i45 - i44) + i47) - i46) + 1)] =
                                        iext->data[(int)((nfcns + 1.0) + 1.0) - 1];
                                    loop_ub = b_wt->size[1];
                                    i44 = r21->size[0] * r21->size[1];
                                    r21->size[0] = 1;
                                    r21->size[1] = loop_ub;
                                    emxEnsureCapacity_int32_T(r21, i44);
                                    for (i44 = 0; i44 < loop_ub; i44++) {
                                        r21->data[r21->size[0] * i44] = i44;
                                    }

                                    temp = iext->data[(int)((nfcns + 1.0) + 1.0) - 1];
                                    dnum = iext->data[(int)((nfcns + 1.0) + 1.0) - 1];
                                    i44 = mtmp->size[0];
                                    mtmp->size[0] = (((i49 - i48) + flag) - ixstart) + 2;
                                    emxEnsureCapacity_real_T1(mtmp, i44);
                                    loop_ub = i49 - i48;
                                    for (i44 = 0; i44 < loop_ub; i44++) {
                                        mtmp->data[i44] = iext->data[i48 + i44];
                                    }

                                    loop_ub = flag - ixstart;
                                    for (i44 = 0; i44 < loop_ub; i44++) {
                                        mtmp->data[(i44 + i49) - i48] = iext->data[ixstart + i44];
                                    }

                                    mtmp->data[((i49 - i48) + flag) - ixstart] = temp;
                                    mtmp->data[(((i49 - i48) + flag) - ixstart) + 1] = dnum;
                                    loop_ub = r21->size[1];
                                    for (i44 = 0; i44 < loop_ub; i44++) {
                                        iext->data[r21->data[r21->size[0] * i44]] = mtmp->data
                                                [(*(int (*)[2])r21->size)[0] * i44];
                                    }

                                    exitg2 = true;
                                } else {
                                    b_l--;
                                }
                            }

                            if (luck != 6) {
                                temp = (nfcns + 1.0) - nfcns;
                                if (2.0 > temp) {
                                    i44 = -2;
                                    i45 = 0;
                                } else {
                                    i44 = -1;
                                    i45 = (int)temp;
                                }

                                temp = (nfcns + 1.0) - nfcns;
                                if (temp > (nfcns + 1.0) - 1.0) {
                                    i46 = 1;
                                    i47 = 0;
                                } else {
                                    i46 = (int)temp;
                                    i47 = (int)((nfcns + 1.0) - 1.0);
                                }

                                /*  Update index */
                                temp = (nfcns + 1.0) - nfcns;
                                if (2.0 > temp) {
                                    i48 = -2;
                                    i49 = 0;
                                } else {
                                    i48 = -1;
                                    i49 = (int)temp;
                                }

                                temp = (nfcns + 1.0) - nfcns;
                                if (temp > (nfcns + 1.0) - 1.0) {
                                    ixstart = 1;
                                    flag = 0;
                                } else {
                                    ixstart = (int)temp;
                                    flag = (int)((nfcns + 1.0) - 1.0);
                                }

                                nu = b_wt->size[0] * b_wt->size[1];
                                b_wt->size[0] = 1;
                                b_wt->size[1] = ((i45 - i44) + i47) - i46;
                                emxEnsureCapacity_real_T(b_wt, nu);
                                loop_ub = i45 - i44;
                                for (nu = 0; nu <= loop_ub - 3; nu++) {
                                    b_wt->data[b_wt->size[0] * (nu + 1)] = iext->data[(i44 + nu) +
                                                                           2];
                                }

                                loop_ub = i47 - i46;
                                for (i47 = 0; i47 <= loop_ub; i47++) {
                                    b_wt->data[b_wt->size[0] * (((i47 + i45) - i44) - 1)] =
                                        iext->data[(i46 + i47) - 1];
                                }

                                loop_ub = b_wt->size[1];
                                i44 = r21->size[0] * r21->size[1];
                                r21->size[0] = 1;
                                r21->size[1] = loop_ub;
                                emxEnsureCapacity_int32_T(r21, i44);
                                for (i44 = 0; i44 < loop_ub; i44++) {
                                    r21->data[r21->size[0] * i44] = i44;
                                }

                                i44 = mtmp->size[0];
                                mtmp->size[0] = ((i49 - i48) + flag) - ixstart;
                                emxEnsureCapacity_real_T1(mtmp, i44);
                                mtmp->data[0] = k1;
                                loop_ub = i49 - i48;
                                for (i44 = 0; i44 <= loop_ub - 3; i44++) {
                                    mtmp->data[i44 + 1] = iext->data[(i48 + i44) + 2];
                                }

                                loop_ub = flag - ixstart;
                                for (i44 = 0; i44 <= loop_ub; i44++) {
                                    mtmp->data[((i44 + i49) - i48) - 1] = iext->data[(ixstart +
                                                                          i44) - 1];
                                }

                                loop_ub = r21->size[1];
                                for (i44 = 0; i44 < loop_ub; i44++) {
                                    iext->data[r21->data[r21->size[0] * i44]] = mtmp->data[(*(int
                                            (*)[2])r21->size)[0] * i44];
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
        err = -1.0;
        k1 = -1.0;

        /*  initialize memory */
        /* x(nzz) = -2; */
        i44 = x2->size[0] * x2->size[1];
        x2->size[0] = 1;
        x2->size[1] = x->size[1] + 1;
        emxEnsureCapacity_real_T(x2, i44);
        loop_ub = x->size[1];
        for (i44 = 0; i44 < loop_ub; i44++) {
            x2->data[x2->size[0] * i44] = x->data[x->size[0] * i44];
        }

        x2->data[x2->size[0] * x->size[1]] = -2.0;
        jchnge = 2.0 * nfcns - 1.0;
        b_j = 1.0 / jchnge;
        b_l = 1.0;
        nu = 0;
        if (((edge[0] == 0.0) && (edge[3] == 0.5)) || (nfcns <= 3.0)) {
            nu = 1;
        }

        if (nu != 1) {
            dtemp = cos(6.2831853071795862 * grid->data[0]);
            dnum = cos(6.2831853071795862 * grid->data[varargin_2 - 1]);
            k1 = 2.0 / (dtemp - dnum);
            err = -(dtemp + dnum) / (dtemp - dnum);
        }

        i44 = a->size[0] * a->size[1];
        a->size[0] = 1;
        a->size[1] = (int)nfcns;
        emxEnsureCapacity_real_T(a, i44);
        for (nut = 0; nut < (int)nfcns; nut++) {
            temp = ((1.0 + (double)nut) - 1.0) * b_j;
            kup = cos(6.2831853071795862 * temp);
            if (nu != 1) {
                kup = (kup - err) / k1;
                dtemp = kup;
                b_acos(&dtemp);
                temp = dtemp / 6.2831853071795862;
            }

            dnum = x2->data[(int)b_l - 1];
            while ((kup <= dnum) && (dnum - kup >= 1.0E-6)) {
                b_l++;
                dnum = x2->data[(int)b_l - 1];
            }

            if (fabs(kup - dnum) < 1.0E-6) {
                a->data[nut] = y->data[(int)b_l - 1];
            } else {
                /*  gee */
                loop_ub = (int)(nfcns + 1.0);
                dtemp = cos(6.2831853071795862 * temp);
                i44 = b_wt->size[0] * b_wt->size[1];
                b_wt->size[0] = 1;
                b_wt->size[1] = (int)(nfcns + 1.0);
                emxEnsureCapacity_real_T(b_wt, i44);
                for (i44 = 0; i44 < loop_ub; i44++) {
                    b_wt->data[b_wt->size[0] * i44] = dtemp - x2->data[i44];
                }

                c_rdivide(ad, b_wt, j);
                i44 = b->size[0];
                b->size[0] = y->size[1];
                emxEnsureCapacity_real_T1(b, i44);
                loop_ub = y->size[1];
                for (i44 = 0; i44 < loop_ub; i44++) {
                    b->data[i44] = y->data[y->size[0] * i44];
                }

                if ((j->size[1] == 1) || (b->size[0] == 1)) {
                    b_y = 0.0;
                    for (i44 = 0; i44 < j->size[1]; i44++) {
                        b_y += j->data[j->size[0] * i44] * b->data[i44];
                    }
                } else {
                    b_y = 0.0;
                    for (i44 = 0; i44 < j->size[1]; i44++) {
                        b_y += j->data[j->size[0] * i44] * b->data[i44];
                    }
                }

                dtemp = b_sum(j);
                a->data[nut] = b_y / dtemp;
            }

            ixstart = 1;
            if ((int)(b_l - 1.0) > 1) {
                ixstart = (int)(b_l - 1.0);
            }

            b_l = ixstart;
        }

        temp = 6.2831853071795862 / jchnge;
        i44 = j->size[0] * j->size[1];
        j->size[0] = 1;
        j->size[1] = (int)nfcns;
        emxEnsureCapacity_real_T(j, i44);
        for (nut = 0; nut < (int)nfcns; nut++) {
            dnum = ((1.0 + (double)nut) - 1.0) * temp;
            if (nfcns - 1.0 < 1.0) {
                j->data[nut] = a->data[0];
            } else {
                if (2.0 > nfcns) {
                    i44 = 1;
                    i45 = 1;
                } else {
                    i44 = 2;
                    i45 = (int)nfcns + 1;
                }

                if (nfcns - 1.0 < 1.0) {
                    i46 = y->size[0] * y->size[1];
                    y->size[0] = 1;
                    y->size[1] = 0;
                    emxEnsureCapacity_real_T(y, i46);
                } else {
                    i46 = y->size[0] * y->size[1];
                    y->size[0] = 1;
                    y->size[1] = (int)floor((nfcns - 1.0) - 1.0) + 1;
                    emxEnsureCapacity_real_T(y, i46);
                    loop_ub = (int)floor((nfcns - 1.0) - 1.0);
                    for (i46 = 0; i46 <= loop_ub; i46++) {
                        y->data[y->size[0] * i46] = 1.0 + (double)i46;
                    }
                }

                i46 = y->size[0] * y->size[1];
                y->size[0] = 1;
                emxEnsureCapacity_real_T(y, i46);
                ixstart = y->size[0];
                flag = y->size[1];
                loop_ub = ixstart * flag;
                for (i46 = 0; i46 < loop_ub; i46++) {
                    y->data[i46] *= dnum;
                }

                b_cos(y);
                i46 = y->size[0] * y->size[1];
                y->size[0] = 1;
                emxEnsureCapacity_real_T(y, i46);
                ixstart = y->size[0];
                flag = y->size[1];
                loop_ub = ixstart * flag;
                for (i46 = 0; i46 < loop_ub; i46++) {
                    y->data[i46] *= 2.0;
                }

                i46 = b->size[0];
                b->size[0] = i45 - i44;
                emxEnsureCapacity_real_T1(b, i46);
                loop_ub = i45 - i44;
                for (i46 = 0; i46 < loop_ub; i46++) {
                    b->data[i46] = a->data[(i44 + i46) - 1];
                }

                if ((y->size[1] == 1) || (i45 - i44 == 1)) {
                    b_y = 0.0;
                    for (i44 = 0; i44 < y->size[1]; i44++) {
                        b_y += y->data[y->size[0] * i44] * b->data[i44];
                    }
                } else {
                    b_y = 0.0;
                    for (i44 = 0; i44 < y->size[1]; i44++) {
                        b_y += y->data[y->size[0] * i44] * b->data[i44];
                    }
                }

                j->data[nut] = a->data[0] + b_y;
            }
        }

        if (2.0 > nfcns) {
            i44 = -1;
            i45 = 0;
        } else {
            i44 = 0;
            i45 = (int)nfcns;
        }

        i46 = iext->size[0];
        iext->size[0] = i45 - i44;
        emxEnsureCapacity_real_T1(iext, i46);
        iext->data[0] = j->data[0] / jchnge;
        loop_ub = i45 - i44;
        for (i45 = 0; i45 <= loop_ub - 2; i45++) {
            iext->data[i45 + 1] = 2.0 * j->data[(i44 + i45) + 1] / jchnge;
        }

        i44 = j->size[0] * j->size[1];
        j->size[0] = 1;
        j->size[1] = (int)nfcns;
        emxEnsureCapacity_real_T(j, i44);
        loop_ub = (int)nfcns;
        for (i44 = 0; i44 < loop_ub; i44++) {
            j->data[i44] = 0.0;
        }

        i44 = l->size[0] * l->size[1];
        l->size[0] = 1;
        l->size[1] = (int)nfcns - 1;
        emxEnsureCapacity_real_T(l, i44);
        loop_ub = (int)nfcns - 1;
        for (i44 = 0; i44 < loop_ub; i44++) {
            l->data[i44] = 0.0;
        }

        if (nu != 1) {
            j->data[0] = 2.0 * iext->data[(int)nfcns - 1] * err + iext->data[(int)
                         nfcns - 2];
            j->data[1] = 2.0 * k1 * iext->data[(int)nfcns - 1];
            l->data[0] = iext->data[(int)nfcns - 3] - iext->data[(int)nfcns - 1];
            for (nut = 0; nut <= (int)nfcns - 3; nut++) {
                if (2 + nut == (int)nfcns - 1) {
                    k1 /= 2.0;
                    err /= 2.0;
                }

                j->data[nut + 2] = 0.0;
                i44 = x2->size[0] * x2->size[1];
                x2->size[0] = 1;
                x2->size[1] = (int)((2.0 + (double)nut) - 1.0) + 1;
                emxEnsureCapacity_real_T(x2, i44);
                loop_ub = (int)((2.0 + (double)nut) - 1.0);
                for (i44 = 0; i44 <= loop_ub; i44++) {
                    x2->data[x2->size[0] * i44] = 1.0 + (double)i44;
                }

                loop_ub = x2->size[0] * x2->size[1];
                for (i44 = 0; i44 < loop_ub; i44++) {
                    a->data[(int)x2->data[i44] - 1] = j->data[(int)x2->data[i44] - 1];
                }

                temp = 2.0 * err;
                loop_ub = x2->size[0] * x2->size[1];
                for (i44 = 0; i44 < loop_ub; i44++) {
                    j->data[(int)x2->data[i44] - 1] = temp * a->data[(int)x2->data[i44] -
                                                      1];
                }

                j->data[1] += 2.0 * a->data[0] * k1;
                i44 = x2->size[0] * x2->size[1];
                x2->size[0] = 1;
                x2->size[1] = (int)(((2.0 + (double)nut) - 1.0) - 1.0) + 1;
                emxEnsureCapacity_real_T(x2, i44);
                loop_ub = (int)(((2.0 + (double)nut) - 1.0) - 1.0);
                for (i44 = 0; i44 <= loop_ub; i44++) {
                    x2->data[x2->size[0] * i44] = 1.0 + (double)i44;
                }

                i44 = x->size[0] * x->size[1];
                x->size[0] = 1;
                x->size[1] = x2->size[1];
                emxEnsureCapacity_real_T(x, i44);
                loop_ub = x2->size[0] * x2->size[1];
                for (i44 = 0; i44 < loop_ub; i44++) {
                    x->data[i44] = a->data[(int)x2->data[i44]];
                }

                i44 = mtmp->size[0];
                mtmp->size[0] = x2->size[0] * x2->size[1];
                emxEnsureCapacity_real_T1(mtmp, i44);
                loop_ub = x2->size[0] * x2->size[1];
                for (i44 = 0; i44 < loop_ub; i44++) {
                    mtmp->data[i44] = (j->data[(int)x2->data[i44] - 1] + l->data[(int)
                                       x2->data[i44] - 1]) + k1 * x->data[i44];
                }

                loop_ub = mtmp->size[0];
                for (i44 = 0; i44 < loop_ub; i44++) {
                    j->data[(int)x2->data[i44] - 1] = mtmp->data[i44];
                }

                i44 = x2->size[0] * x2->size[1];
                x2->size[0] = 1;
                x2->size[1] = (int)(((2.0 + (double)nut) + 1.0) - 3.0) + 1;
                emxEnsureCapacity_real_T(x2, i44);
                loop_ub = (int)(((2.0 + (double)nut) + 1.0) - 3.0);
                for (i44 = 0; i44 <= loop_ub; i44++) {
                    x2->data[x2->size[0] * i44] = 3.0 + (double)i44;
                }

                i44 = x->size[0] * x->size[1];
                x->size[0] = 1;
                x->size[1] = x2->size[1];
                emxEnsureCapacity_real_T(x, i44);
                loop_ub = x2->size[0] * x2->size[1];
                for (i44 = 0; i44 < loop_ub; i44++) {
                    x->data[i44] = a->data[(int)x2->data[i44] - 2];
                }

                i44 = mtmp->size[0];
                mtmp->size[0] = x2->size[0] * x2->size[1];
                emxEnsureCapacity_real_T1(mtmp, i44);
                loop_ub = x2->size[0] * x2->size[1];
                for (i44 = 0; i44 < loop_ub; i44++) {
                    mtmp->data[i44] = j->data[(int)x2->data[i44] - 1] + k1 * x->data[i44];
                }

                loop_ub = mtmp->size[0];
                for (i44 = 0; i44 < loop_ub; i44++) {
                    j->data[(int)x2->data[i44] - 1] = mtmp->data[i44];
                }

                if (2 + nut != (int)nfcns - 1) {
                    i44 = x2->size[0] * x2->size[1];
                    x2->size[0] = 1;
                    x2->size[1] = (int)((2.0 + (double)nut) - 1.0) + 1;
                    emxEnsureCapacity_real_T(x2, i44);
                    loop_ub = (int)((2.0 + (double)nut) - 1.0);
                    for (i44 = 0; i44 <= loop_ub; i44++) {
                        x2->data[x2->size[0] * i44] = 1.0 + (double)i44;
                    }

                    loop_ub = x2->size[0] * x2->size[1];
                    for (i44 = 0; i44 < loop_ub; i44++) {
                        l->data[(int)x2->data[i44] - 1] = -a->data[(int)x2->data[i44] - 1];
                    }

                    l->data[0] += iext->data[((int)nfcns - nut) - 4];
                }
            }

            loop_ub = (int)nfcns;
            for (i44 = 0; i44 < loop_ub; i44++) {
                iext->data[i44] = j->data[i44];
            }
        }

        /*  alpha must be at lease >=3 */
        if (nfcns <= 3.0) {
            /* alpha(nfcns + 1) = 0; */
            /* alpha(nfcns + 2) = 0; */
            i44 = j->size[0] * j->size[1];
            j->size[0] = 1;
            j->size[1] = iext->size[0] + 2;
            emxEnsureCapacity_real_T(j, i44);
            loop_ub = iext->size[0];
            for (i44 = 0; i44 < loop_ub; i44++) {
                j->data[j->size[0] * i44] = iext->data[i44];
            }

            j->data[j->size[0] * iext->size[0]] = 0.0;
            j->data[j->size[0] * (iext->size[0] + 1)] = 0.0;
        } else {
            i44 = j->size[0] * j->size[1];
            j->size[0] = 1;
            j->size[1] = iext->size[0];
            emxEnsureCapacity_real_T(j, i44);
            loop_ub = iext->size[0];
            for (i44 = 0; i44 < loop_ub; i44++) {
                j->data[j->size[0] * i44] = iext->data[i44];
            }
        }

        /* alpha=alpha'; */
        /*  now that's done! */
        if (nodd != 0.0) {
            i44 = x->size[0] * x->size[1];
            x->size[0] = 1;
            x->size[1] = (int)floor(-((0.0 - (nfcns - 1.0)) - -1.0)) + 1;
            emxEnsureCapacity_real_T(x, i44);
            loop_ub = (int)floor(-((0.0 - (nfcns - 1.0)) - -1.0));
            for (i44 = 0; i44 <= loop_ub; i44++) {
                x->data[x->size[0] * i44] = j->data[(int)((nfcns + 1.0) + (-1.0 -
                                                    (double)i44)) - 1];
            }

            i44 = h->size[0] * h->size[1];
            h->size[0] = 1;
            h->size[1] = x->size[1] + 1;
            emxEnsureCapacity_real_T(h, i44);
            loop_ub = x->size[1];
            for (i44 = 0; i44 < loop_ub; i44++) {
                h->data[h->size[0] * i44] = 0.5 * x->data[x->size[0] * i44];
            }

            h->data[h->size[0] * x->size[1]] = j->data[0];
        } else {
            if ((nfcns - (nfcns - 1.0)) + 2.0 > nfcns) {
                i44 = 0;
                i45 = 1;
            } else {
                i44 = (int)nfcns - 1;
                i45 = -1;
            }

            i46 = x2->size[0] * x2->size[1];
            x2->size[0] = 1;
            x2->size[1] = (int)floor(-((0.0 - (nfcns - 1.0)) - -2.0)) + 1;
            emxEnsureCapacity_real_T(x2, i46);
            loop_ub = (int)floor(-((0.0 - (nfcns - 1.0)) - -2.0));
            for (i46 = 0; i46 <= loop_ub; i46++) {
                x2->data[x2->size[0] * i46] = j->data[(int)((nfcns + 1.0) + (double)(int)
                                                      (-2.0 - (double)i46)) - 1];
            }

            i46 = h->size[0] * h->size[1];
            h->size[0] = 1;
            h->size[1] = 2 + x2->size[1];
            emxEnsureCapacity_real_T(h, i46);
            h->data[0] = 0.25 * j->data[(int)nfcns - 1];
            loop_ub = x2->size[1];
            for (i46 = 0; i46 < loop_ub; i46++) {
                h->data[h->size[0] * (i46 + 1)] = 0.25 * (x2->data[x2->size[0] * i46] +
                                                  j->data[i44 + i45 * i46]);
            }

            h->data[h->size[0] * (1 + x2->size[1])] = 0.25 * (2.0 * j->data[0] +
                    j->data[1]);
        }
    }

    emxFree_real_T(&mtmp);
    emxFree_int32_T(&b_iext);
    emxFree_int8_T(&r22);
    emxFree_real_T(&b_wt);
    emxFree_real_T(&a);
    emxFree_real_T(&b);
    emxFree_int32_T(&r21);
    emxFree_real_T(&x2);
    emxFree_real_T(&l);
    emxFree_real_T(&ad);
    emxFree_real_T(&x);
    emxFree_real_T(&y);
    emxFree_real_T(&iext);
    emxFree_real_T(&j);
    *dev = b_dev;
    *valid = b_valid;
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
    int ii_size_idx_1;
    int ii;
    boolean_T exitg1;
    int ii_data[1];
    int b_ii_size_idx_1;
    boolean_T b_a;
    int b_ii_data[1];
    signed char varargin_1_data[1];
    int varargin_2_data[1];
    int c_ii_size_idx_1;
    int trz_len_data_idx_0;
    double d0;

    /*  end freqs */
    idx = 0;
    ii_size_idx_1 = 1;
    ii = b_size[1];
    exitg1 = false;
    while ((!exitg1) && (ii > 0)) {
        if (b_data[ii - 1] != 0.0) {
            idx = 1;
            ii_data[0] = ii;
            exitg1 = true;
        } else {
            ii--;
        }
    }

    if (idx == 0) {
        ii_size_idx_1 = 0;
    }

    idx = 0;
    b_ii_size_idx_1 = 1;
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
        b_ii_size_idx_1 = 0;
    }

    for (ii = 0; ii < ii_size_idx_1; ii++) {
        varargin_1_data[ii] = (signed char)((signed char)b_size[1] - (signed char)
                                            ii_data[ii]);
    }

    idx = a->size[1];
    for (ii = 0; ii < b_ii_size_idx_1; ii++) {
        varargin_2_data[ii] = idx - b_ii_data[ii];
    }

    if (ii_size_idx_1 <= b_ii_size_idx_1) {
        c_ii_size_idx_1 = ii_size_idx_1;
    } else {
        c_ii_size_idx_1 = 0;
    }

    if (1 <= c_ii_size_idx_1) {
        idx = varargin_1_data[0];
        trz_len_data_idx_0 = varargin_2_data[0];
        if (idx < trz_len_data_idx_0) {
            trz_len_data_idx_0 = idx;
        }
    }

    if (trz_len_data_idx_0 > 0) {
        d0 = (double)b_size[1] - (double)trz_len_data_idx_0;
        if (1.0 > d0) {
            idx = 0;
        } else {
            idx = (int)d0;
        }

        bR_size[0] = 1;
        bR_size[1] = idx;
        for (ii = 0; ii < idx; ii++) {
            bR_data[bR_size[0] * ii] = b_data[ii];
        }

        d0 = (double)a->size[1] - (double)trz_len_data_idx_0;
        if (1.0 > d0) {
            idx = 0;
        } else {
            idx = (int)d0;
        }

        ii = aR->size[0] * aR->size[1];
        aR->size[0] = 1;
        aR->size[1] = idx;
        emxEnsureCapacity_creal_T(aR, ii);
        for (ii = 0; ii < idx; ii++) {
            aR->data[aR->size[0] * ii] = a->data[ii];
        }
    } else {
        bR_size[0] = 1;
        bR_size[1] = b_size[1];
        idx = b_size[0] * b_size[1];
        for (ii = 0; ii < idx; ii++) {
            bR_data[ii] = b_data[ii];
        }

        ii = aR->size[0] * aR->size[1];
        aR->size[0] = 1;
        aR->size[1] = a->size[1];
        emxEnsureCapacity_creal_T(aR, ii);
        idx = a->size[0] * a->size[1];
        for (ii = 0; ii < idx; ii++) {
            aR->data[ii] = a->data[ii];
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
    double b_u1;
    double q;
    if (!((!rtIsNaN(u0)) && (!rtIsInf(u0)) && ((!rtIsNaN(u1)) && (!rtIsInf(u1))))) {
        y = rtNaN;
    } else {
        if (u1 < 0.0) {
            b_u1 = ceil(u1);
        } else {
            b_u1 = floor(u1);
        }

        if ((u1 != 0.0) && (u1 != b_u1)) {
            q = fabs(u0 / u1);
            if (fabs(q - floor(q + 0.5)) <= DBL_EPSILON * q) {
                y = 0.0 * u0;
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
    emxArray_real_T *r15;
    int i38;
    int loop_ub;
    emxArray_real_T *digw;
    emxArray_creal_T *s;
    emxArray_creal_T *y;
    boolean_T b9;
    int k;
    double s_re;
    double s_im;
    emxInit_real_T(&r15, 2);

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
    i38 = r15->size[0] * r15->size[1];
    r15->size[0] = 1;
    r15->size[1] = w->size[1];
    emxEnsureCapacity_real_T(r15, i38);
    loop_ub = w->size[0] * w->size[1];
    for (i38 = 0; i38 < loop_ub; i38++) {
        r15->data[i38] = 6.2831853071795862 * w->data[i38];
    }

    emxInit_real_T(&digw, 2);
    emxInit_creal_T(&s, 2);
    rdivide(r15, Fs, digw);

    /*  Convert from Hz to rad/sample for computational purposes */
    i38 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = digw->size[1];
    emxEnsureCapacity_creal_T(s, i38);
    loop_ub = digw->size[0] * digw->size[1];
    emxFree_real_T(&r15);
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
    b9 = (y->size[1] == 0);
    if (!b9) {
        i38 = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity_creal_T(y, i38);
        loop_ub = y->size[1];
        for (i38 = 0; i38 < loop_ub; i38++) {
            y->data[y->size[0] * i38].re = b[0];
            y->data[y->size[0] * i38].im = 0.0;
        }

        for (k = 0; k < 18; k++) {
            i38 = y->size[0] * y->size[1];
            y->size[0] = 1;
            y->size[1] = s->size[1];
            emxEnsureCapacity_creal_T(y, i38);
            loop_ub = s->size[0] * s->size[1];
            for (i38 = 0; i38 < loop_ub; i38++) {
                s_re = s->data[i38].re * y->data[i38].re - s->data[i38].im * y->data[i38]
                       .im;
                s_im = s->data[i38].re * y->data[i38].im + s->data[i38].im * y->data[i38]
                       .re;
                y->data[i38].re = s_re + b[k + 1];
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
        s_re = digw->data[i38] * 0.0;
        s_im = digw->data[i38];
        s->data[i38].re = 18.0 * s_re;
        s->data[i38].im = 18.0 * s_im;
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
    emxArray_real_T *r16;
    int i39;
    int loop_ub;
    emxArray_real_T *digw;
    emxArray_creal_T *s;
    emxArray_creal_T *y;
    boolean_T b10;
    int k;
    double s_re;
    double s_im;
    emxInit_real_T(&r16, 2);

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
    i39 = r16->size[0] * r16->size[1];
    r16->size[0] = 1;
    r16->size[1] = w->size[1];
    emxEnsureCapacity_real_T(r16, i39);
    loop_ub = w->size[0] * w->size[1];
    for (i39 = 0; i39 < loop_ub; i39++) {
        r16->data[i39] = 6.2831853071795862 * w->data[i39];
    }

    emxInit_real_T(&digw, 2);
    emxInit_creal_T(&s, 2);
    rdivide(r16, Fs, digw);

    /*  Convert from Hz to rad/sample for computational purposes */
    i39 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = digw->size[1];
    emxEnsureCapacity_creal_T(s, i39);
    loop_ub = digw->size[0] * digw->size[1];
    emxFree_real_T(&r16);
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
    b10 = (y->size[1] == 0);
    if (!b10) {
        i39 = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity_creal_T(y, i39);
        loop_ub = y->size[1];
        for (i39 = 0; i39 < loop_ub; i39++) {
            y->data[y->size[0] * i39].re = b[0];
            y->data[y->size[0] * i39].im = 0.0;
        }

        for (k = 0; k < 84; k++) {
            i39 = y->size[0] * y->size[1];
            y->size[0] = 1;
            y->size[1] = s->size[1];
            emxEnsureCapacity_creal_T(y, i39);
            loop_ub = s->size[0] * s->size[1];
            for (i39 = 0; i39 < loop_ub; i39++) {
                s_re = s->data[i39].re * y->data[i39].re - s->data[i39].im * y->data[i39]
                       .im;
                s_im = s->data[i39].re * y->data[i39].im + s->data[i39].im * y->data[i39]
                       .re;
                y->data[i39].re = s_re + b[k + 1];
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
        s_re = digw->data[i39] * 0.0;
        s_im = digw->data[i39];
        s->data[i39].re = 84.0 * s_re;
        s->data[i39].im = 84.0 * s_im;
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
    if ((a > b) || rtIsNaN(b)) {
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
    int i59;
    int i;
    int im1n;
    int in;
    creal_T alpha1;
    int c;
    double alpha1_re;
    double alpha1_im;
    double xnorm;
    double beta1;
    int jy;
    int knt;
    boolean_T b_tau;
    double ai;
    int lastv;
    int lastc;
    int k;
    boolean_T exitg1;
    creal_T b_c;
    int ix;
    int exitg2;
    n = a->size[0];
    if (a->size[0] < 1) {
        ntau = 0;
    } else {
        ntau = a->size[0] - 1;
    }

    emxInit_creal_T1(&tau, 1);
    emxInit_creal_T1(&work, 1);
    i59 = tau->size[0];
    tau->size[0] = ntau;
    emxEnsureCapacity_creal_T1(tau, i59);
    ntau = a->size[0];
    i59 = work->size[0];
    work->size[0] = ntau;
    emxEnsureCapacity_creal_T1(work, i59);
    for (i59 = 0; i59 < ntau; i59++) {
        work->data[i59].re = 0.0;
        work->data[i59].im = 0.0;
    }

    for (i = 0; i + 1 < n; i++) {
        im1n = i * n + 2;
        in = (i + 1) * n;
        alpha1 = a->data[(i + a->size[0] * i) + 1];
        ntau = i + 3;
        if (!(ntau < n)) {
            ntau = n;
        }

        ntau += i * n;
        c = (n - i) - 2;
        alpha1_re = 0.0;
        alpha1_im = 0.0;
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
                    i59 = (ntau + c) - 1;
                    do {
                        knt++;
                        for (k = ntau; k <= i59; k++) {
                            xnorm = a->data[k - 1].re;
                            ai = a->data[k - 1].im;
                            a->data[k - 1].re = 9.9792015476736E+291 * xnorm - 0.0 * ai;
                            a->data[k - 1].im = 9.9792015476736E+291 * ai + 0.0 * xnorm;
                        }

                        beta1 *= 9.9792015476736E+291;
                        alpha1.re *= 9.9792015476736E+291;
                        alpha1.im *= 9.9792015476736E+291;
                    } while (!(fabs(beta1) >= 1.0020841800044864E-292));

                    beta1 = xdlapy3(alpha1.re, alpha1.im, xnrm2(c, a, ntau));
                    if (alpha1.re >= 0.0) {
                        beta1 = -beta1;
                    }

                    xnorm = beta1 - alpha1.re;
                    if (0.0 - alpha1.im == 0.0) {
                        alpha1_re = xnorm / beta1;
                        alpha1_im = 0.0;
                    } else if (xnorm == 0.0) {
                        alpha1_re = 0.0;
                        alpha1_im = (0.0 - alpha1.im) / beta1;
                    } else {
                        alpha1_re = xnorm / beta1;
                        alpha1_im = (0.0 - alpha1.im) / beta1;
                    }

                    b_c.re = alpha1.re - beta1;
                    b_c.im = alpha1.im;
                    xscal(c, recip(b_c), a, ntau);
                    for (k = 1; k <= knt; k++) {
                        beta1 *= 1.0020841800044864E-292;
                    }

                    alpha1.re = beta1;
                    alpha1.im = 0.0;
                } else {
                    xnorm = beta1 - a->data[(i + a->size[0] * i) + 1].re;
                    ai = 0.0 - a->data[(i + a->size[0] * i) + 1].im;
                    if (ai == 0.0) {
                        alpha1_re = xnorm / beta1;
                        alpha1_im = 0.0;
                    } else if (xnorm == 0.0) {
                        alpha1_re = 0.0;
                        alpha1_im = ai / beta1;
                    } else {
                        alpha1_re = xnorm / beta1;
                        alpha1_im = ai / beta1;
                    }

                    b_c.re = a->data[(i + a->size[0] * i) + 1].re - beta1;
                    b_c.im = a->data[(i + a->size[0] * i) + 1].im;
                    xscal(c, recip(b_c), a, ntau);
                    alpha1.re = beta1;
                    alpha1.im = 0.0;
                }
            }
        }

        tau->data[i].re = alpha1_re;
        tau->data[i].im = alpha1_im;
        a->data[(i + a->size[0] * i) + 1].re = 1.0;
        a->data[(i + a->size[0] * i) + 1].im = 0.0;
        c = (n - i) - 3;
        jy = (i + im1n) - 1;
        b_tau = ((tau->data[i].re != 0.0) || (tau->data[i].im != 0.0));
        if (b_tau) {
            lastv = c + 2;
            ntau = jy + c;
            exitg1 = false;
            while ((!exitg1) && (lastv > 0)) {
                b_tau = ((a->data[ntau + 1].re == 0.0) && (a->data[ntau + 1].im == 0.0));
                if (b_tau) {
                    lastv--;
                    ntau--;
                } else {
                    exitg1 = true;
                }
            }

            lastc = n;
            exitg1 = false;
            while ((!exitg1) && (lastc > 0)) {
                ntau = in + lastc;
                c = ntau;
                do {
                    exitg2 = 0;
                    if ((n > 0) && (c <= ntau + (lastv - 1) * n)) {
                        b_tau = ((a->data[c - 1].re != 0.0) || (a->data[c - 1].im != 0.0));
                        if (b_tau) {
                            exitg2 = 1;
                        } else {
                            c += n;
                        }
                    } else {
                        lastc--;
                        exitg2 = 2;
                    }
                } while (exitg2 == 0);

                if (exitg2 == 1) {
                    exitg1 = true;
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
                i59 = (in + n * (lastv - 1)) + 1;
                knt = in + 1;
                while ((n > 0) && (knt <= i59)) {
                    b_c.re = a->data[ix].re - 0.0 * a->data[ix].im;
                    b_c.im = a->data[ix].im + 0.0 * a->data[ix].re;
                    ntau = 0;
                    k = (knt + lastc) - 1;
                    for (c = knt; c <= k; c++) {
                        xnorm = a->data[c - 1].re * b_c.re - a->data[c - 1].im * b_c.im;
                        ai = a->data[c - 1].re * b_c.im + a->data[c - 1].im * b_c.re;
                        work->data[ntau].re += xnorm;
                        work->data[ntau].im += ai;
                        ntau++;
                    }

                    ix++;
                    knt += n;
                }
            }

            alpha1_re = -tau->data[i].re;
            alpha1_im = -tau->data[i].im;
            if (!((alpha1_re == 0.0) && (alpha1_im == 0.0))) {
                ntau = in;
                for (knt = 1; knt <= lastv; knt++) {
                    b_tau = ((a->data[jy].re != 0.0) || (a->data[jy].im != 0.0));
                    if (b_tau) {
                        b_c.re = a->data[jy].re * alpha1_re + a->data[jy].im * alpha1_im;
                        b_c.im = a->data[jy].re * alpha1_im - a->data[jy].im * alpha1_re;
                        ix = 0;
                        i59 = lastc + ntau;
                        for (k = ntau; k + 1 <= i59; k++) {
                            xnorm = work->data[ix].re * b_c.re - work->data[ix].im * b_c.im;
                            ai = work->data[ix].re * b_c.im + work->data[ix].im * b_c.re;
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
        alpha1_re = tau->data[i].re;
        alpha1_im = -tau->data[i].im;
        if ((alpha1_re != 0.0) || (alpha1_im != 0.0)) {
            lastv = c + 2;
            ntau = im1n + c;
            exitg1 = false;
            while ((!exitg1) && (lastv > 0)) {
                b_tau = ((a->data[ntau + 1].re == 0.0) && (a->data[ntau + 1].im == 0.0));
                if (b_tau) {
                    lastv--;
                    ntau--;
                } else {
                    exitg1 = true;
                }
            }

            lastc = (n - i) - 1;
            exitg1 = false;
            while ((!exitg1) && (lastc > 0)) {
                ntau = jy + (lastc - 1) * n;
                c = ntau;
                do {
                    exitg2 = 0;
                    if (c <= (ntau + lastv) - 1) {
                        b_tau = ((a->data[c - 1].re != 0.0) || (a->data[c - 1].im != 0.0));
                        if (b_tau) {
                            exitg2 = 1;
                        } else {
                            c++;
                        }
                    } else {
                        lastc--;
                        exitg2 = 2;
                    }
                } while (exitg2 == 0);

                if (exitg2 == 1) {
                    exitg1 = true;
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
                i59 = jy + n * (lastc - 1);
                knt = jy;
                while ((n > 0) && (knt <= i59)) {
                    ix = im1n;
                    b_c.re = 0.0;
                    b_c.im = 0.0;
                    k = (knt + lastv) - 1;
                    for (c = knt - 1; c + 1 <= k; c++) {
                        b_c.re += a->data[c].re * a->data[ix].re + a->data[c].im * a->
                                  data[ix].im;
                        b_c.im += a->data[c].re * a->data[ix].im - a->data[c].im * a->
                                  data[ix].re;
                        ix++;
                    }

                    work->data[ntau].re += b_c.re - 0.0 * b_c.im;
                    work->data[ntau].im += b_c.im + 0.0 * b_c.re;
                    ntau++;
                    knt += n;
                }
            }

            alpha1_re = -alpha1_re;
            alpha1_im = -alpha1_im;
            if (!((alpha1_re == 0.0) && (alpha1_im == 0.0))) {
                ntau = jy - 1;
                jy = 0;
                for (knt = 1; knt <= lastc; knt++) {
                    b_tau = ((work->data[jy].re != 0.0) || (work->data[jy].im != 0.0));
                    if (b_tau) {
                        b_c.re = work->data[jy].re * alpha1_re + work->data[jy].im *
                                 alpha1_im;
                        b_c.im = work->data[jy].re * alpha1_im - work->data[jy].im *
                                 alpha1_re;
                        ix = im1n;
                        i59 = lastv + ntau;
                        for (k = ntau; k + 1 <= i59; k++) {
                            xnorm = a->data[ix].re * b_c.re - a->data[ix].im * b_c.im;
                            ai = a->data[ix].re * b_c.im + a->data[ix].im * b_c.re;
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
    int i60;
    int k;
    double x_re;
    double x_im;
    i60 = (ix0 + n) - 1;
    for (k = ix0; k <= i60; k++) {
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
    int jcol;
    int ii;
    int nzcount;
    double anrm;
    boolean_T exitg1;
    double absxk;
    boolean_T ilascl;
    double anrmto;
    int ilo;
    double ctoc;
    int ihi;
    boolean_T notdone;
    int exitg3;
    double cfrom1;
    int i;
    double cto1;
    int j;
    double mul;
    creal_T b_At;
    creal_T c_At;
    double c;
    creal_T atmp;
    boolean_T exitg4;
    int exitg2;
    boolean_T d_At;
    double stemp_re;
    emxInit_creal_T(&At, 2);
    jcol = At->size[0] * At->size[1];
    At->size[0] = A->size[0];
    At->size[1] = A->size[1];
    emxEnsureCapacity_creal_T(At, jcol);
    ii = A->size[0] * A->size[1];
    for (jcol = 0; jcol < ii; jcol++) {
        At->data[jcol] = A->data[jcol];
    }

    nzcount = 0;
    jcol = alpha1->size[0];
    alpha1->size[0] = At->size[0];
    emxEnsureCapacity_creal_T1(alpha1, jcol);
    ii = At->size[0];
    for (jcol = 0; jcol < ii; jcol++) {
        alpha1->data[jcol].re = 0.0;
        alpha1->data[jcol].im = 0.0;
    }

    jcol = beta1->size[0];
    beta1->size[0] = At->size[0];
    emxEnsureCapacity_creal_T1(beta1, jcol);
    ii = At->size[0];
    for (jcol = 0; jcol < ii; jcol++) {
        beta1->data[jcol].re = 0.0;
        beta1->data[jcol].im = 0.0;
    }

    if (!((At->size[0] == 0) || (At->size[1] == 0))) {
        anrm = 0.0;
        jcol = 0;
        exitg1 = false;
        while ((!exitg1) && (jcol <= At->size[0] * At->size[1] - 1)) {
            absxk = rt_hypotd_snf(At->data[jcol].re, At->data[jcol].im);
            if (rtIsNaN(absxk)) {
                anrm = rtNaN;
                exitg1 = true;
            } else {
                if (absxk > anrm) {
                    anrm = absxk;
                }

                jcol++;
            }
        }

        if (!((!rtIsInf(anrm)) && (!rtIsNaN(anrm)))) {
            jcol = alpha1->size[0];
            alpha1->size[0] = At->size[0];
            emxEnsureCapacity_creal_T1(alpha1, jcol);
            ii = At->size[0];
            for (jcol = 0; jcol < ii; jcol++) {
                alpha1->data[jcol].re = rtNaN;
                alpha1->data[jcol].im = 0.0;
            }

            jcol = beta1->size[0];
            beta1->size[0] = At->size[0];
            emxEnsureCapacity_creal_T1(beta1, jcol);
            ii = At->size[0];
            for (jcol = 0; jcol < ii; jcol++) {
                beta1->data[jcol].re = rtNaN;
                beta1->data[jcol].im = 0.0;
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

                    jcol = At->size[0] * At->size[1];
                    emxEnsureCapacity_creal_T(At, jcol);
                    jcol = At->size[0];
                    ii = At->size[1];
                    ii *= jcol;
                    for (jcol = 0; jcol < ii; jcol++) {
                        At->data[jcol].re *= mul;
                        At->data[jcol].im *= mul;
                    }
                }
            }

            ilo = 0;
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
                        jcol = 1;
                        exitg4 = false;
                        while ((!exitg4) && (jcol <= ihi)) {
                            d_At = ((At->data[(ii + At->size[0] * (jcol - 1)) - 1].re != 0.0) ||
                                    (At->data[(ii + At->size[0] * (jcol - 1)) - 1].im != 0.0));
                            if (d_At || (ii == jcol)) {
                                if (nzcount == 0) {
                                    j = jcol;
                                    nzcount = 1;
                                    jcol++;
                                } else {
                                    nzcount = 2;
                                    exitg4 = true;
                                }
                            } else {
                                jcol++;
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
                        nzcount = At->size[0];
                        if (i != ihi) {
                            for (jcol = 0; jcol + 1 <= nzcount; jcol++) {
                                atmp = At->data[(i + At->size[0] * jcol) - 1];
                                At->data[(i + At->size[0] * jcol) - 1] = At->data[(ihi +
                                        At->size[0] * jcol) - 1];
                                At->data[(ihi + At->size[0] * jcol) - 1] = atmp;
                            }
                        }

                        if (j != ihi) {
                            for (jcol = 0; jcol + 1 <= ihi; jcol++) {
                                atmp = At->data[jcol + At->size[0] * (j - 1)];
                                At->data[jcol + At->size[0] * (j - 1)] = At->data[jcol +
                                        At->size[0] * (ihi - 1)];
                                At->data[jcol + At->size[0] * (ihi - 1)] = atmp;
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
                        jcol = ilo + 1;
                        exitg1 = false;
                        while ((!exitg1) && (jcol <= ihi)) {
                            nzcount = 0;
                            i = ihi;
                            j = jcol;
                            ii = ilo + 1;
                            exitg4 = false;
                            while ((!exitg4) && (ii <= ihi)) {
                                d_At = ((At->data[(ii + At->size[0] * (jcol - 1)) - 1].re != 0.0)
                                        || (At->data[(ii + At->size[0] * (jcol - 1)) - 1].im !=
                                            0.0));
                                if (d_At || (ii == jcol)) {
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
                                jcol++;
                            }
                        }

                        if (!notdone) {
                            exitg2 = 1;
                        } else {
                            nzcount = At->size[0];
                            if (i != ilo + 1) {
                                for (jcol = ilo; jcol + 1 <= nzcount; jcol++) {
                                    atmp = At->data[(i + At->size[0] * jcol) - 1];
                                    At->data[(i + At->size[0] * jcol) - 1] = At->data[ilo +
                                            At->size[0] * jcol];
                                    At->data[ilo + At->size[0] * jcol] = atmp;
                                }
                            }

                            if (j != ilo + 1) {
                                for (jcol = 0; jcol + 1 <= ihi; jcol++) {
                                    atmp = At->data[jcol + At->size[0] * (j - 1)];
                                    At->data[jcol + At->size[0] * (j - 1)] = At->data[jcol +
                                            At->size[0] * ilo];
                                    At->data[jcol + At->size[0] * ilo] = atmp;
                                }
                            }

                            ilo++;
                            if (ilo + 1 == ihi) {
                                exitg2 = 1;
                            }
                        }
                    } while (exitg2 == 0);
                }
            }

            nzcount = At->size[0];
            if ((!(At->size[0] <= 1)) && (!(ihi < ilo + 3))) {
                for (jcol = ilo; jcol + 1 < ihi - 1; jcol++) {
                    for (ii = ihi - 1; ii + 1 > jcol + 2; ii--) {
                        b_At = At->data[(ii + At->size[0] * jcol) - 1];
                        c_At = At->data[ii + At->size[0] * jcol];
                        xzlartg(b_At, c_At, &c, &atmp, &At->data[(ii + At->size[0] * jcol) -
                                                          1]);
                        At->data[ii + At->size[0] * jcol].re = 0.0;
                        At->data[ii + At->size[0] * jcol].im = 0.0;
                        for (j = jcol + 1; j + 1 <= nzcount; j++) {
                            absxk = atmp.re * At->data[ii + At->size[0] * j].re - atmp.im *
                                    At->data[ii + At->size[0] * j].im;
                            ctoc = atmp.re * At->data[ii + At->size[0] * j].im + atmp.im *
                                   At->data[ii + At->size[0] * j].re;
                            stemp_re = c * At->data[(ii + At->size[0] * j) - 1].re + absxk;
                            absxk = c * At->data[(ii + At->size[0] * j) - 1].im + ctoc;
                            ctoc = At->data[(ii + At->size[0] * j) - 1].re;
                            cfrom1 = At->data[(ii + At->size[0] * j) - 1].im;
                            cto1 = At->data[(ii + At->size[0] * j) - 1].im;
                            mul = At->data[(ii + At->size[0] * j) - 1].re;
                            At->data[ii + At->size[0] * j].re = c * At->data[ii + At->size[0] *
                                                                j].re - (atmp.re * ctoc + atmp.im * cfrom1);
                            At->data[ii + At->size[0] * j].im = c * At->data[ii + At->size[0] *
                                                                j].im - (atmp.re * cto1 - atmp.im * mul);
                            At->data[(ii + At->size[0] * j) - 1].re = stemp_re;
                            At->data[(ii + At->size[0] * j) - 1].im = absxk;
                        }

                        atmp.re = -atmp.re;
                        atmp.im = -atmp.im;
                        for (i = 0; i + 1 <= ihi; i++) {
                            absxk = atmp.re * At->data[i + At->size[0] * (ii - 1)].re -
                                    atmp.im * At->data[i + At->size[0] * (ii - 1)].im;
                            ctoc = atmp.re * At->data[i + At->size[0] * (ii - 1)].im + atmp.im
                                   * At->data[i + At->size[0] * (ii - 1)].re;
                            stemp_re = c * At->data[i + At->size[0] * ii].re + absxk;
                            absxk = c * At->data[i + At->size[0] * ii].im + ctoc;
                            ctoc = At->data[i + At->size[0] * ii].re;
                            cfrom1 = At->data[i + At->size[0] * ii].im;
                            cto1 = At->data[i + At->size[0] * ii].im;
                            mul = At->data[i + At->size[0] * ii].re;
                            At->data[i + At->size[0] * (ii - 1)].re = c * At->data[i +
                                    At->size[0] * (ii - 1)].re - (atmp.re * ctoc + atmp.im * cfrom1);
                            At->data[i + At->size[0] * (ii - 1)].im = c * At->data[i +
                                    At->size[0] * (ii - 1)].im - (atmp.re * cto1 - atmp.im * mul);
                            At->data[i + At->size[0] * ii].re = stemp_re;
                            At->data[i + At->size[0] * ii].im = absxk;
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

                    jcol = alpha1->size[0];
                    emxEnsureCapacity_creal_T1(alpha1, jcol);
                    ii = alpha1->size[0];
                    for (jcol = 0; jcol < ii; jcol++) {
                        alpha1->data[jcol].re *= mul;
                        alpha1->data[jcol].im *= mul;
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
    int b_info;
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
    boolean_T exitg2;
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

    b_info = -1;
    if ((A->size[0] == 1) && (A->size[1] == 1)) {
        ihi = 1;
    }

    jm1 = alpha1->size[0];
    alpha1->size[0] = A->size[0];
    emxEnsureCapacity_creal_T1(alpha1, jm1);
    jp1 = A->size[0];
    for (jm1 = 0; jm1 < jp1; jm1++) {
        alpha1->data[jm1].re = 0.0;
        alpha1->data[jm1].im = 0.0;
    }

    jm1 = beta1->size[0];
    beta1->size[0] = A->size[0];
    emxEnsureCapacity_creal_T1(beta1, jm1);
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
                    exitg2 = false;
                    while ((!exitg2) && (j + 1 >= ilo)) {
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
                            exitg2 = true;
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
                    emxEnsureCapacity_creal_T1(alpha1, jm1);
                    for (jm1 = 0; jm1 < jp1; jm1++) {
                        alpha1->data[jm1].re = rtNaN;
                        alpha1->data[jm1].im = 0.0;
                    }

                    jp1 = beta1->size[0];
                    jm1 = beta1->size[0];
                    beta1->size[0] = jp1;
                    emxEnsureCapacity_creal_T1(beta1, jm1);
                    for (jm1 = 0; jm1 < jp1; jm1++) {
                        beta1->data[jm1].re = rtNaN;
                        beta1->data[jm1].im = 0.0;
                    }

                    b_info = 0;
                    exitg1 = 1;
                } else if (goto60) {
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
                            ctemp.re = ascale * b_A->data[j + b_A->size[0] * j].re - shift.re *
                                       bscale;
                            ctemp.im = ascale * b_A->data[j + b_A->size[0] * j].im - shift.im *
                                       bscale;
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
                                 (b_A->data[j + b_A->size[0] * (j - 1)].im)) * temp2 <= anorm *
                                b_atol) {
                                goto90 = true;
                                exitg2 = true;
                            } else {
                                jp1 = j;
                                j--;
                            }
                        }

                        if (!goto90) {
                            istart = ifirst;
                            ctemp.re = ascale * b_A->data[(ifirst + b_A->size[0] * (ifirst - 1))
                                                          - 1].re - shift.re * bscale;
                            ctemp.im = ascale * b_A->data[(ifirst + b_A->size[0] * (ifirst - 1))
                                                          - 1].im - shift.im * bscale;
                            goto90 = true;
                        }
                    }

                    if (goto90) {
                        goto90 = false;
                        b_ascale.re = ascale * b_A->data[istart + b_A->size[0] * (istart - 1)]
                                      .re;
                        b_ascale.im = ascale * b_A->data[istart + b_A->size[0] * (istart - 1)]
                                      .im;
                        b_xzlartg(ctemp, b_ascale, &imAij, &shift);
                        j = istart;
                        jm1 = istart - 2;
                        while (j < ilast + 1) {
                            if (j > istart) {
                                b_ascale = b_A->data[(j + b_A->size[0] * jm1) - 1];
                                c_A = b_A->data[j + b_A->size[0] * jm1];
                                xzlartg(b_ascale, c_A, &imAij, &shift, &b_A->data[(j + b_A->
                                        size[0] * jm1) - 1]);
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
                                                                       b_A->size[0] * jp1].re - (shift.re * anorm + shift.im * reAij);
                                b_A->data[j + b_A->size[0] * jp1].im = imAij * b_A->data[j +
                                                                       b_A->size[0] * jp1].im - (shift.re * scale - shift.im * sumsq);
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
                                b_A->data[i + b_A->size[0] * (j - 1)].re = imAij * b_A->data[i +
                                        b_A->size[0] * (j - 1)].re - (shift.re * anorm + shift.im *
                                                                      reAij);
                                b_A->data[i + b_A->size[0] * (j - 1)].im = imAij * b_A->data[i +
                                        b_A->size[0] * (j - 1)].im - (shift.re * scale - shift.im *
                                                                      sumsq);
                                b_A->data[i + b_A->size[0] * j].re = ad22_re;
                                b_A->data[i + b_A->size[0] * j].im = ad22_im;
                            }

                            jm1 = j - 1;
                            j++;
                        }
                    }

                    jiter++;
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
            b_info = ilast;
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
    *info = b_info + 1;
}

/*
 * Arguments    : creal_T *alpha1
 *                creal_T *x
 * Return Type  : creal_T
 */
static creal_T xzlarfg(creal_T *alpha1, creal_T *x)
{
    creal_T tau;
    double xnorm;
    double beta1;
    int knt;
    double ai;
    creal_T b_alpha1;
    double x_re;
    double x_im;
    int k;
    tau.re = 0.0;
    tau.im = 0.0;
    xnorm = rt_hypotd_snf(x->re, x->im);
    if ((xnorm != 0.0) || (alpha1->im != 0.0)) {
        beta1 = xdlapy3(alpha1->re, alpha1->im, xnorm);
        if (alpha1->re >= 0.0) {
            beta1 = -beta1;
        }

        if (fabs(beta1) < 1.0020841800044864E-292) {
            knt = 0;
            do {
                knt++;
                x->re *= 9.9792015476736E+291;
                x->im *= 9.9792015476736E+291;
                beta1 *= 9.9792015476736E+291;
                alpha1->re *= 9.9792015476736E+291;
                alpha1->im *= 9.9792015476736E+291;
            } while (!(fabs(beta1) >= 1.0020841800044864E-292));

            beta1 = xdlapy3(alpha1->re, alpha1->im, rt_hypotd_snf(x->re, x->im));
            if (alpha1->re >= 0.0) {
                beta1 = -beta1;
            }

            xnorm = beta1 - alpha1->re;
            ai = 0.0 - alpha1->im;
            if (ai == 0.0) {
                tau.re = xnorm / beta1;
                tau.im = 0.0;
            } else if (xnorm == 0.0) {
                tau.re = 0.0;
                tau.im = ai / beta1;
            } else {
                tau.re = xnorm / beta1;
                tau.im = ai / beta1;
            }

            b_alpha1.re = alpha1->re - beta1;
            b_alpha1.im = alpha1->im;
            *alpha1 = recip(b_alpha1);
            xnorm = alpha1->re;
            ai = alpha1->im;
            x_re = x->re;
            x_im = x->im;
            x->re = xnorm * x_re - ai * x_im;
            x->im = xnorm * x_im + ai * x_re;
            for (k = 1; k <= knt; k++) {
                beta1 *= 1.0020841800044864E-292;
            }

            alpha1->re = beta1;
            alpha1->im = 0.0;
        } else {
            xnorm = beta1 - alpha1->re;
            ai = 0.0 - alpha1->im;
            if (ai == 0.0) {
                tau.re = xnorm / beta1;
                tau.im = 0.0;
            } else if (xnorm == 0.0) {
                tau.re = 0.0;
                tau.im = ai / beta1;
            } else {
                tau.re = xnorm / beta1;
                tau.im = ai / beta1;
            }

            b_alpha1.re = alpha1->re - beta1;
            b_alpha1.im = alpha1->im;
            *alpha1 = recip(b_alpha1);
            xnorm = alpha1->re;
            ai = alpha1->im;
            x_re = x->re;
            x_im = x->im;
            x->re = xnorm * x_re - ai * x_im;
            x->im = xnorm * x_im + ai * x_re;
            alpha1->re = beta1;
            alpha1->im = 0.0;
        }
    }

    return tau;
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
 *                emxArray_real_T *c
 *                double *d
 * Return Type  : void
 */
static void zp2ss_cg(emxArray_creal_T *a, emxArray_real_T *b,
                     emxArray_real_T *c,
                     double *d)
{
    int i6;
    creal_T b_c[3];
    int r1;
    double wn;
    double a22;
    int k;
    double den[3];
    double b1[2];
    double t[4];
    double B[4];
    int r2;
    double Y[4];
    emxArray_int8_T *reshapes_f2;
    emxArray_cint8_T *b_a;
    emxArray_cint8_T *result;
    int i7;
    signed char a_re;
    signed char a_im;
    emxArray_real_T *b_b1;
    emxArray_real_T *varargin_2;
    int result_im;
    emxArray_real_T *r3;

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
    emxEnsureCapacity_real_T1(b, i6);
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
    b_c[0].re = 1.0;
    b_c[0].im = 0.0;
    for (r1 = 0; r1 < 2; r1++) {
        wn = -(0.86602540378443871 + -1.7320508075688774 * (double)r1);
        a22 = b_c[r1].re;
        b_c[r1 + 1].re = 0.49999999999999978 * b_c[r1].re - wn * b_c[r1].im;
        b_c[r1 + 1].im = 0.49999999999999978 * b_c[r1].im + wn * a22;
        k = r1;
        while (k + 1 > 1) {
            wn = 0.86602540378443871 + -1.7320508075688774 * (double)r1;
            b_c[1].re -= -0.49999999999999978 * b_c[0].re - wn * b_c[0].im;
            b_c[1].im -= -0.49999999999999978 * b_c[0].im + wn * b_c[0].re;
            k = 0;
        }
    }

    for (i6 = 0; i6 < 3; i6++) {
        den[i6] = b_c[i6].re;
    }

    for (k = 0; k < 2; k++) {
        b1[k] = rt_hypotd_snf(-0.49999999999999978, 0.86602540378443871 +
                              -1.7320508075688774 * (double)k);
    }

    wn = sqrt(b1[0] * b1[1]);
    b1[0] = 1.0;
    b1[1] = 1.0 / wn;
    for (i6 = 0; i6 < 4; i6++) {
        t[i6] = 0.0;
    }

    /*  Balancing transformation */
    B[0] = -den[1];
    B[2] = -den[2];
    for (r1 = 0; r1 < 2; r1++) {
        t[r1 + (r1 << 1)] = b1[r1];
        B[1 + (r1 << 1)] = 1.0 - (double)r1;
    }

    if (t[1] > t[0]) {
        r1 = 1;
        r2 = 0;
    } else {
        r1 = 0;
        r2 = 1;
    }

    wn = t[r2] / t[r1];
    a22 = t[2 + r2] - wn * t[2 + r1];
    for (k = 0; k < 2; k++) {
        Y[1 + (k << 1)] = (B[r2 + (k << 1)] - B[r1 + (k << 1)] * wn) / a22;
        Y[k << 1] = (B[r1 + (k << 1)] - Y[1 + (k << 1)] * t[2 + r1]) / t[r1];
    }

    if (t[1] > t[0]) {
        r1 = 1;
        r2 = 0;
    } else {
        r1 = 0;
        r2 = 1;
    }

    emxInit_int8_T(&reshapes_f2, 2);
    wn = t[r2] / t[r1];
    b1[1] = ((1.0 - (double)r2) - (1.0 - (double)r1) * wn) / (t[2 + r2] - wn * t[2
            + r1]);
    b1[0] = ((1.0 - (double)r1) - b1[1] * t[2 + r1]) / t[r1];

    /*  [a,b,c,d] = series(a,b,c,d,a1,b1,c1,d1); */
    /*  Next lines perform series connection */
    i6 = reshapes_f2->size[0] * reshapes_f2->size[1];
    reshapes_f2->size[0] = 1;
    reshapes_f2->size[1] = 2;
    emxEnsureCapacity_int8_T(reshapes_f2, i6);
    for (i6 = 0; i6 < 2; i6++) {
        reshapes_f2->data[i6] = 0;
    }

    emxInit_cint8_T(&b_a, 2);
    i6 = b_a->size[0] * b_a->size[1];
    b_a->size[0] = a->size[0];
    b_a->size[1] = a->size[1];
    emxEnsureCapacity_cint8_T(b_a, i6);
    r1 = a->size[1];
    for (i6 = 0; i6 < r1; i6++) {
        r2 = a->size[0];
        for (i7 = 0; i7 < r2; i7++) {
            a_re = (signed char)a->data[i7 + a->size[0] * i6].re;
            a_im = (signed char)a->data[i7 + a->size[0] * i6].im;
            b_a->data[i7 + b_a->size[0] * i6].re = a_re;
            b_a->data[i7 + b_a->size[0] * i6].im = a_im;
        }
    }

    emxInit_cint8_T(&result, 2);
    i6 = result->size[0] * result->size[1];
    result->size[0] = 1;
    result->size[1] = 1 + reshapes_f2->size[1];
    emxEnsureCapacity_cint8_T(result, i6);
    for (i6 = 0; i6 < 1; i6++) {
        for (i7 = 0; i7 < 1; i7++) {
            result->data[0] = b_a->data[0];
        }
    }

    emxFree_cint8_T(&b_a);
    r1 = reshapes_f2->size[1];
    for (i6 = 0; i6 < r1; i6++) {
        r2 = reshapes_f2->size[0];
        for (i7 = 0; i7 < r2; i7++) {
            result->data[i7 + result->size[0] * (i6 + 1)].re = reshapes_f2->data[i7 +
                    reshapes_f2->size[0] * i6];
            result->data[i7 + result->size[0] * (i6 + 1)].im = 0;
        }
    }

    emxFree_int8_T(&reshapes_f2);
    emxInit_real_T(&b_b1, 2);
    i6 = b_b1->size[0] * b_b1->size[1];
    b_b1->size[0] = 2;
    b_b1->size[1] = c->size[1];
    emxEnsureCapacity_real_T(b_b1, i6);
    for (i6 = 0; i6 < 2; i6++) {
        r1 = c->size[1];
        for (i7 = 0; i7 < r1; i7++) {
            b_b1->data[i6 + b_b1->size[0] * i7] = b1[i6] * c->data[c->size[0] * i7];
        }
    }

    for (i6 = 0; i6 < 2; i6++) {
        for (i7 = 0; i7 < 2; i7++) {
            B[i6 + (i7 << 1)] = 0.0;
            for (r1 = 0; r1 < 2; r1++) {
                B[i6 + (i7 << 1)] += Y[i6 + (r1 << 1)] * t[r1 + (i7 << 1)];
            }
        }
    }

    emxInit_real_T(&varargin_2, 2);
    i6 = varargin_2->size[0] * varargin_2->size[1];
    varargin_2->size[0] = 2;
    varargin_2->size[1] = b_b1->size[1] + 2;
    emxEnsureCapacity_real_T(varargin_2, i6);
    r1 = b_b1->size[1];
    for (i6 = 0; i6 < r1; i6++) {
        for (i7 = 0; i7 < 2; i7++) {
            varargin_2->data[i7 + varargin_2->size[0] * i6] = b_b1->data[i7 +
                    b_b1->size[0] * i6];
        }
    }

    for (i6 = 0; i6 < 2; i6++) {
        for (i7 = 0; i7 < 2; i7++) {
            varargin_2->data[i7 + varargin_2->size[0] * (i6 + b_b1->size[1])] = B[i7 +
                    (i6 << 1)];
        }
    }

    emxFree_real_T(&b_b1);
    if (!(result->size[1] == 0)) {
        r1 = result->size[1];
    } else {
        r1 = 3;
    }

    r2 = !(result->size[1] == 0);
    i6 = a->size[0] * a->size[1];
    a->size[0] = r2 + 2;
    a->size[1] = r1;
    emxEnsureCapacity_creal_T(a, i6);
    for (i6 = 0; i6 < r1; i6++) {
        for (i7 = 0; i7 < r2; i7++) {
            k = result->data[i7 + r2 * i6].re;
            result_im = result->data[i7 + r2 * i6].im;
            a->data[i7 + a->size[0] * i6].re = k;
            a->data[i7 + a->size[0] * i6].im = result_im;
        }
    }

    emxFree_cint8_T(&result);
    for (i6 = 0; i6 < r1; i6++) {
        for (i7 = 0; i7 < 2; i7++) {
            a->data[(i7 + r2) + a->size[0] * i6].re = varargin_2->data[i7 + (i6 << 1)];
            a->data[(i7 + r2) + a->size[0] * i6].im = 0.0;
        }
    }

    emxFree_real_T(&varargin_2);
    r1 = b->size[0];
    i6 = b->size[0];
    b->size[0] = r1 + 2;
    emxEnsureCapacity_real_T1(b, i6);
    for (i6 = 0; i6 < 2; i6++) {
        b->data[r1 + i6] = b1[i6] * 0.0;
    }

    for (i6 = 0; i6 < 2; i6++) {
        b1[i6] = 0.0;
        for (i7 = 0; i7 < 2; i7++) {
            b1[i6] += (double)i7 * t[i7 + (i6 << 1)];
        }
    }

    emxInit_real_T(&r3, 2);
    i6 = r3->size[0] * r3->size[1];
    r3->size[0] = 1;
    r3->size[1] = c->size[1] + 2;
    emxEnsureCapacity_real_T(r3, i6);
    r1 = c->size[1];
    for (i6 = 0; i6 < r1; i6++) {
        r3->data[r3->size[0] * i6] = 0.0 * c->data[c->size[0] * i6];
    }

    for (i6 = 0; i6 < 2; i6++) {
        r3->data[r3->size[0] * (i6 + c->size[1])] = b1[i6];
    }

    i6 = c->size[0] * c->size[1];
    c->size[0] = 1;
    c->size[1] = r3->size[1];
    emxEnsureCapacity_real_T(c, i6);
    r1 = r3->size[1];
    for (i6 = 0; i6 < r1; i6++) {
        c->data[c->size[0] * i6] = r3->data[r3->size[0] * i6];
    }

    emxFree_real_T(&r3);

    /*  Apply gain k: */
    i6 = c->size[0] * c->size[1];
    c->size[0] = 1;
    emxEnsureCapacity_real_T(c, i6);
    r1 = c->size[0];
    r2 = c->size[1];
    r1 *= r2;
    for (i6 = 0; i6 < r1; i6++) {
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
 * Return Type  : void
 */
void internal_design_filter_cg(double Rdata, double Fpass, double Fstop, double
                               caldiv, double FIR, double HB1, double PLL_mult, double Apass, double Astop,
                               double phEQ, double HB2, double HB3, const char Type[7], const char RxTx[2],
                               double RFbw, double DAC_div, double converter_rate, double PLL_rate, double
                               Fcenter, double wnom, double FIRdBmin, double int_FIR, double maxTaps, short
                               outputTaps[128], double *numOutputTaps, double *filterGain)
{
    emxArray_creal_T *a1;
    emxArray_creal_T *a2;
    double b1[4];
    double b2[2];
    creal_T a2_data[2];
    int a2_size[2];
    int b1_size[2];
    int i0;
    double b1_data[4];
    int b2_size[2];
    double b2_data[4];
    double hb1_coeff[15];
    static const double dv0[15] = { -0.00323486328125, 0.0, 0.01910400390625, 0.0,
                                    -0.07049560546875, 0.0, 0.30450439453125, 0.5, 0.30450439453125, 0.0,
                                    -0.07049560546875, 0.0, 0.01910400390625, 0.0, -0.00323486328125
                                  };

    static const double dv1[15] = { -0.00390625, 0.0, 0.0205078125, 0.0,
                                    -0.07177734375, 0.0, 0.30224609375, 0.49462890625, 0.30224609375, 0.0,
                                    -0.07177734375, 0.0, 0.0205078125, 0.0, -0.00390625
                                  };

    int hb3_coeff_size[2];
    int dec_int3_coeff_size[2];
    double hb3_coeff_data[5];
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
                                    -0.02398681640625, 0.0, 0.00506591796875, 0.00335693359375
                                  };

    int hb3;
    int dec_int3;
    double sigma;
    signed char i1;
    char enables[4];
    double w[2048];
    double phi[2048];
    static const double hb2_coeff[7] = { -0.03515625, 0.0, 0.28515625, 0.5,
                                         0.28515625, 0.0, -0.03515625
                                       };

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
    int loop_ub;
    emxArray_real_T *omega;
    emxArray_creal_T *rg1;
    emxArray_creal_T *c_combinedResponse;
    emxArray_creal_T *r0;
    double apnd;
    double absa;
    double rg1_im;
    double re;
    double im;
    emxArray_real_T *sw;
    emxArray_real_T *fg2;
    int n;
    emxArray_real_T *omega2;
    emxArray_creal_T *rgN;
    emxArray_real_T *b_omega2;
    emxArray_creal_T *d_combinedResponse;
    emxArray_real_T *a;
    emxArray_real_T *wg;
    emxArray_real_T *weight;
    boolean_T exitg1;
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
    emxArray_real_T *F4;
    emxArray_real_T *b_W1;
    int exitg2;
    int i2;
    int b_loop_ub;
    boolean_T valid;
    int i3;
    double firTapsPreScale[128];
    double b_firTapsPreScale[128];
    short i4;
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
        for (i0 = 0; i0 < 2; i0++) {
            b1_data[i0] = b2[i0];
        }

        i0 = a1->size[0] * a1->size[1];
        a1->size[0] = 1;
        a1->size[1] = 2;
        emxEnsureCapacity_creal_T(a1, i0);
        for (i0 = 0; i0 < 2; i0++) {
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
        a2->size[1] = 2;
        emxEnsureCapacity_creal_T(a2, i0);
        for (i0 = 0; i0 < 2; i0++) {
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
        hb3 = 50;
        dec_int3 = 49;
    } else if (HB3 == 3.0) {
        hb3 = 49;
        dec_int3 = 51;
    } else {
        hb3 = 49;
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
    enables[2] = (signed char)hb3;
    enables[3] = (signed char)dec_int3;

    /*  Find out the best fit delay on passband */
    memset(&w[0], 0, sizeof(double) << 11);
    memset(&phi[0], 0, sizeof(double) << 11);
    w[0] = -Fpass;
    for (hb3 = 0; hb3 < 2047; hb3++) {
        w[hb3 + 1] = w[0] - 2.0 * w[0] * (2.0 + (double)hb3) / 2048.0;
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

    b_power(b_combinedResponse, invariance);
    for (i0 = 0; i0 < 2048; i0++) {
        b_combinedResponse[i0] = combinedResponse[i0].im;
    }

    b_power(b_combinedResponse, b_phi);
    for (i0 = 0; i0 < 2048; i0++) {
        invariance[i0] += b_phi[i0];
    }

    phi[0] = rt_atan2d_snf(combinedResponse[0].im, combinedResponse[0].re);
    for (hb3 = 0; hb3 < 2047; hb3++) {
        sigma = rt_atan2d_snf(combinedResponse[hb3 + 1].im, combinedResponse[hb3 + 1]
                              .re) - phi[hb3];
        phi[hb3 + 1] = phi[hb3] + (sigma - 6.2831853071795862 * floor(sigma /
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
    if (!((Gpass < Gstop - 1.0) || rtIsNaN(Gstop - 1.0))) {
        Gpass = Gstop - 1.0;
    }

    emxInit_real_T(&fg, 2);
    i0 = fg->size[0] * fg->size[1];
    fg->size[0] = 1;
    fg->size[1] = (int)(Gpass + 1.0);
    emxEnsureCapacity_real_T(fg, i0);
    loop_ub = (int)(Gpass + 1.0);
    for (i0 = 0; i0 < loop_ub; i0++) {
        fg->data[i0] = 0.0;
    }

    emxInit_real_T(&omega, 2);
    i0 = omega->size[0] * omega->size[1];
    omega->size[0] = 1;
    omega->size[1] = (int)(Gpass + 1.0);
    emxEnsureCapacity_real_T(omega, i0);
    loop_ub = (int)(Gpass + 1.0);
    for (i0 = 0; i0 < loop_ub; i0++) {
        omega->data[i0] = 0.0;
    }

    /*  passband */
    for (hb3 = 0; hb3 < (int)(Gpass + 1.0); hb3++) {
        fg->data[hb3] = ((1.0 + (double)hb3) - 1.0) / 16384.0;
        omega->data[hb3] = fg->data[hb3] * clkFIR;
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
    i0 = rg1->size[0] * rg1->size[1];
    rg1->size[0] = 1;
    rg1->size[1] = c_combinedResponse->size[1];
    emxEnsureCapacity_creal_T(rg1, i0);
    loop_ub = c_combinedResponse->size[0] * c_combinedResponse->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
        sigmax = c_combinedResponse->data[i0].re;
        apnd = c_combinedResponse->data[i0].im;
        absa = rg1->data[i0].re;
        rg1_im = rg1->data[i0].im;
        rg1->data[i0].re = sigmax * absa - apnd * rg1_im;
        rg1->data[i0].im = sigmax * rg1_im + apnd * absa;
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
    b_rdivide(r0, rg1, c_combinedResponse);
    sigma = Gpass + 1.0;

    /*  Expand memory correctly */
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
        sw->size[0] = 1;
        sw->size[1] = (int)(8192.0 - Gstop) + 1;
        emxEnsureCapacity_real_T(sw, i0);
        loop_ub = (int)(8192.0 - Gstop);
        for (i0 = 0; i0 <= loop_ub; i0++) {
            sw->data[sw->size[0] * i0] = Gstop + (double)i0;
        }
    } else {
        sigmax = floor((8192.0 - Gstop) + 0.5);
        apnd = Gstop + sigmax;
        absa = fabs(Gstop);
        if (!(absa > 8192.0)) {
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
        emxEnsureCapacity_real_T(sw, i0);
        if (n > 0) {
            sw->data[0] = Gstop;
            if (n > 1) {
                sw->data[n - 1] = apnd;
                hb3 = (n - 1) / 2;
                for (dec_int3 = 1; dec_int3 < hb3; dec_int3++) {
                    sw->data[dec_int3] = Gstop + (double)dec_int3;
                    sw->data[(n - dec_int3) - 1] = apnd - (double)dec_int3;
                }

                if (hb3 << 1 == n - 1) {
                    sw->data[hb3] = (Gstop + apnd) / 2.0;
                } else {
                    sw->data[hb3] = Gstop + (double)hb3;
                    sw->data[hb3 + 1] = apnd - (double)hb3;
                }
            }
        }
    }

    emxInit_real_T(&fg2, 2);
    i0 = fg2->size[0] * fg2->size[1];
    fg2->size[0] = 1;
    fg2->size[1] = (int)((unsigned int)sw->size[1] + fg->size[1]);
    emxEnsureCapacity_real_T(fg2, i0);
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
        sw->size[0] = 1;
        sw->size[1] = (int)(8192.0 - Gstop) + 1;
        emxEnsureCapacity_real_T(sw, i0);
        loop_ub = (int)(8192.0 - Gstop);
        for (i0 = 0; i0 <= loop_ub; i0++) {
            sw->data[sw->size[0] * i0] = Gstop + (double)i0;
        }
    } else {
        sigmax = floor((8192.0 - Gstop) + 0.5);
        apnd = Gstop + sigmax;
        absa = fabs(Gstop);
        if (!(absa > 8192.0)) {
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
        emxEnsureCapacity_real_T(sw, i0);
        if (n > 0) {
            sw->data[0] = Gstop;
            if (n > 1) {
                sw->data[n - 1] = apnd;
                hb3 = (n - 1) / 2;
                for (dec_int3 = 1; dec_int3 < hb3; dec_int3++) {
                    sw->data[dec_int3] = Gstop + (double)dec_int3;
                    sw->data[(n - dec_int3) - 1] = apnd - (double)dec_int3;
                }

                if (hb3 << 1 == n - 1) {
                    sw->data[hb3] = (Gstop + apnd) / 2.0;
                } else {
                    sw->data[hb3] = Gstop + (double)hb3;
                    sw->data[hb3 + 1] = apnd - (double)hb3;
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
        omega2->data[i0] = omega->data[omega->size[0] * i0];
    }

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
        sw->size[0] = 1;
        sw->size[1] = (int)(8192.0 - Gstop) + 1;
        emxEnsureCapacity_real_T(sw, i0);
        loop_ub = (int)(8192.0 - Gstop);
        for (i0 = 0; i0 <= loop_ub; i0++) {
            sw->data[sw->size[0] * i0] = Gstop + (double)i0;
        }
    } else {
        sigmax = floor((8192.0 - Gstop) + 0.5);
        apnd = Gstop + sigmax;
        absa = fabs(Gstop);
        if (!(absa > 8192.0)) {
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
        emxEnsureCapacity_real_T(sw, i0);
        if (n > 0) {
            sw->data[0] = Gstop;
            if (n > 1) {
                sw->data[n - 1] = apnd;
                hb3 = (n - 1) / 2;
                for (dec_int3 = 1; dec_int3 < hb3; dec_int3++) {
                    sw->data[dec_int3] = Gstop + (double)dec_int3;
                    sw->data[(n - dec_int3) - 1] = apnd - (double)dec_int3;
                }

                if (hb3 << 1 == n - 1) {
                    sw->data[hb3] = (Gstop + apnd) / 2.0;
                } else {
                    sw->data[hb3] = Gstop + (double)hb3;
                    sw->data[hb3 + 1] = apnd - (double)hb3;
                }
            }
        }
    }

    emxInit_creal_T(&rgN, 2);
    i0 = rgN->size[0] * rgN->size[1];
    rgN->size[0] = 1;
    rgN->size[1] = (int)((unsigned int)sw->size[1] + c_combinedResponse->size[1]);
    emxEnsureCapacity_creal_T(rgN, i0);
    loop_ub = (int)((unsigned int)sw->size[1] + c_combinedResponse->size[1]);
    for (i0 = 0; i0 < loop_ub; i0++) {
        rgN->data[i0].re = 0.0;
        rgN->data[i0].im = 0.0;
    }

    loop_ub = c_combinedResponse->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
        rgN->data[i0] = c_combinedResponse->data[c_combinedResponse->size[0] * i0];
    }

    /*  stop band */
    for (hb3 = 0; hb3 < (int)(8192.0 + (1.0 - Gstop)); hb3++) {
        sigma++;
        fg2->data[(int)sigma - 1] = (Gstop + (double)hb3) / 16384.0;
        omega2->data[(int)sigma - 1] = fg2->data[(int)sigma - 1] * clkFIR;
        rgN->data[(int)sigma - 1].re = 0.0;
        rgN->data[(int)sigma - 1].im = 0.0;
    }

    /*  Generate responses then convolve */
    if (Gpass + 2.0 > omega2->size[1]) {
        i0 = 0;
        n = 0;
    } else {
        i0 = (int)(Gpass + 2.0) - 1;
        n = omega2->size[1];
    }

    emxInit_real_T(&b_omega2, 2);
    hb3 = b_omega2->size[0] * b_omega2->size[1];
    b_omega2->size[0] = 1;
    b_omega2->size[1] = n - i0;
    emxEnsureCapacity_real_T(b_omega2, hb3);
    loop_ub = n - i0;
    for (n = 0; n < loop_ub; n++) {
        b_omega2->data[b_omega2->size[0] * n] = omega2->data[i0 + n];
    }

    b_generateCascadedResponseRx(enables, b_omega2, converter_rate, hb1_coeff,
                                 hb2_coeff, hb3_coeff_data, hb3_coeff_size, dec_int3_coeff_data,
                                 dec_int3_coeff_size, c_combinedResponse);
    if (Gpass + 2.0 > omega2->size[1]) {
        i0 = 0;
        n = 0;
    } else {
        i0 = (int)(Gpass + 2.0) - 1;
        n = omega2->size[1];
    }

    hb3 = b_omega2->size[0] * b_omega2->size[1];
    b_omega2->size[0] = 1;
    b_omega2->size[1] = n - i0;
    emxEnsureCapacity_real_T(b_omega2, hb3);
    loop_ub = n - i0;
    for (n = 0; n < loop_ub; n++) {
        b_omega2->data[b_omega2->size[0] * n] = omega2->data[i0 + n];
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
    emxInit_real_T(&a, 2);
    if (b_strcmp(RxTx)) {
        rdivide(fg, dBinv(-Astop), omega);
    } else {
        sigma = FIR;
        b_sqrt(&sigma);
        i0 = a->size[0] * a->size[1];
        a->size[0] = 1;
        a->size[1] = fg->size[1];
        emxEnsureCapacity_real_T(a, i0);
        loop_ub = fg->size[0] * fg->size[1];
        for (i0 = 0; i0 < loop_ub; i0++) {
            a->data[i0] = sigma * fg->data[i0];
        }

        rdivide(a, dBinv(-Astop), omega);
    }

    emxInit_real_T(&wg, 2);
    sigma = dBinv(FIRdBmin);
    hb3 = omega->size[1];
    i0 = wg->size[0] * wg->size[1];
    wg->size[0] = 1;
    wg->size[1] = omega->size[1];
    emxEnsureCapacity_real_T(wg, i0);
    for (dec_int3 = 0; dec_int3 + 1 <= hb3; dec_int3++) {
        if ((omega->data[dec_int3] > sigma) || rtIsNaN(sigma)) {
            sigmax = omega->data[dec_int3];
        } else {
            sigmax = sigma;
        }

        wg->data[dec_int3] = sigmax;
    }

    if (phEQ == -1.0) {
        b_abs(rgN, sw);
        i0 = rgN->size[0] * rgN->size[1];
        rgN->size[0] = 1;
        rgN->size[1] = sw->size[1];
        emxEnsureCapacity_creal_T(rgN, i0);
        loop_ub = sw->size[0] * sw->size[1];
        for (i0 = 0; i0 < loop_ub; i0++) {
            rgN->data[i0].re = sw->data[i0];
            rgN->data[i0].im = 0.0;
        }
    }

    emxInit_real_T(&weight, 2);
    b_abs(rg1, sw);
    rdivide(sw, dBinv(Apass / 2.0) - 1.0, fg);
    i0 = weight->size[0] * weight->size[1];
    weight->size[0] = 1;
    weight->size[1] = fg->size[1] + wg->size[1];
    emxEnsureCapacity_real_T(weight, i0);
    loop_ub = fg->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
        weight->data[weight->size[0] * i0] = fg->data[fg->size[0] * i0];
    }

    loop_ub = wg->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
        weight->data[weight->size[0] * (i0 + fg->size[1])] = wg->data[wg->size[0] *
                i0];
    }

    hb3 = 1;
    n = weight->size[1];
    sigma = weight->data[0];
    if (weight->size[1] > 1) {
        if (rtIsNaN(weight->data[0])) {
            dec_int3 = 2;
            exitg1 = false;
            while ((!exitg1) && (dec_int3 <= n)) {
                hb3 = dec_int3;
                if (!rtIsNaN(weight->data[dec_int3 - 1])) {
                    sigma = weight->data[dec_int3 - 1];
                    exitg1 = true;
                } else {
                    dec_int3++;
                }
            }
        }

        if (hb3 < weight->size[1]) {
            while (hb3 + 1 <= n) {
                if (weight->data[hb3] > sigma) {
                    sigma = weight->data[hb3];
                }

                hb3++;
            }
        }
    }

    i0 = a->size[0] * a->size[1];
    a->size[0] = 1;
    a->size[1] = weight->size[1];
    emxEnsureCapacity_real_T(a, i0);
    loop_ub = weight->size[0] * weight->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
        a->data[i0] = weight->data[i0];
    }

    rdivide(a, sigma, weight);

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
        F1->data[F1->size[0] * i0] = fg2->data[i0] * 2.0;
    }

    if (Gpass + 2.0 > fg2->size[1]) {
        i0 = 0;
        n = 0;
    } else {
        i0 = (int)(Gpass + 2.0) - 1;
        n = fg2->size[1];
    }

    emxInit_real_T(&F2, 2);
    hb3 = F2->size[0] * F2->size[1];
    F2->size[0] = 1;
    F2->size[1] = n - i0;
    emxEnsureCapacity_real_T(F2, hb3);
    loop_ub = n - i0;
    for (n = 0; n < loop_ub; n++) {
        F2->data[F2->size[0] * n] = fg2->data[i0 + n] * 2.0;
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
        A1->data[A1->size[0] * i0] = fg->data[i0];
    }

    if (Gpass + 2.0 > fg->size[1]) {
        i0 = 0;
        n = 0;
    } else {
        i0 = (int)(Gpass + 2.0) - 1;
        n = fg->size[1];
    }

    emxInit_real_T(&A2, 2);
    hb3 = A2->size[0] * A2->size[1];
    A2->size[0] = 1;
    A2->size[1] = n - i0;
    emxEnsureCapacity_real_T(A2, hb3);
    loop_ub = n - i0;
    for (n = 0; n < loop_ub; n++) {
        A2->data[A2->size[0] * n] = fg->data[i0 + n];
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
        W1->data[W1->size[0] * i0] = weight->data[i0];
    }

    if (Gpass + 2.0 > weight->size[1]) {
        i0 = 0;
        n = 0;
    } else {
        i0 = (int)(Gpass + 2.0) - 1;
        n = weight->size[1];
    }

    emxInit_real_T(&W2, 2);
    hb3 = W2->size[0] * W2->size[1];
    W2->size[0] = 1;
    W2->size[1] = n - i0;
    emxEnsureCapacity_real_T(W2, hb3);
    loop_ub = n - i0;
    for (n = 0; n < loop_ub; n++) {
        W2->data[W2->size[0] * n] = weight->data[i0 + n];
    }

    emxInit_real_T(&tap_store, 2);

    /*  % Determine the number of taps for FIR */
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
    sigma = maxTaps / 16.0;
    i0 = tap_store->size[0] * tap_store->size[1];
    tap_store->size[0] = (int)sigma;
    tap_store->size[1] = (int)maxTaps;
    emxEnsureCapacity_real_T(tap_store, i0);
    loop_ub = (int)sigma * (int)maxTaps;
    for (i0 = 0; i0 < loop_ub; i0++) {
        tap_store->data[i0] = 0.0;
    }

    emxInit_real_T1(&Apass_actual_vector, 1);
    sigma = maxTaps / 16.0;
    i0 = Apass_actual_vector->size[0];
    Apass_actual_vector->size[0] = (int)sigma;
    emxEnsureCapacity_real_T1(Apass_actual_vector, i0);
    loop_ub = (int)sigma;
    for (i0 = 0; i0 < loop_ub; i0++) {
        Apass_actual_vector->data[i0] = 0.0;
    }

    emxInit_real_T1(&Astop_actual_vector, 1);
    sigma = maxTaps / 16.0;
    i0 = Astop_actual_vector->size[0];
    Astop_actual_vector->size[0] = (int)sigma;
    emxEnsureCapacity_real_T1(Astop_actual_vector, i0);
    loop_ub = (int)sigma;
    for (i0 = 0; i0 < loop_ub; i0++) {
        Astop_actual_vector->data[i0] = 0.0;
    }

    i = 1U;

    /*  Design filter */
    emxInit_real_T(&ccoef, 2);
    emxInit_real_T(&F4, 2);
    emxInit_real_T(&b_W1, 2);
    do {
        exitg2 = 0;
        if (int_FIR != 0.0) {
            b1[0] = F1->data[0];
            b1[1] = F1->data[F1->size[1] - 1];
            b1[2] = F2->data[0];
            b1[3] = F2->data[F2->size[1] - 1];
            i0 = a->size[0] * a->size[1];
            a->size[0] = 1;
            a->size[1] = A1->size[1] + A2->size[1];
            emxEnsureCapacity_real_T(a, i0);
            loop_ub = A1->size[1];
            for (i0 = 0; i0 < loop_ub; i0++) {
                a->data[a->size[0] * i0] = A1->data[A1->size[0] * i0];
            }

            loop_ub = A2->size[1];
            for (i0 = 0; i0 < loop_ub; i0++) {
                a->data[a->size[0] * (i0 + A1->size[1])] = A2->data[A2->size[0] * i0];
            }

            i0 = b_omega2->size[0] * b_omega2->size[1];
            b_omega2->size[0] = 1;
            b_omega2->size[1] = F1->size[1] + F2->size[1];
            emxEnsureCapacity_real_T(b_omega2, i0);
            loop_ub = F1->size[1];
            for (i0 = 0; i0 < loop_ub; i0++) {
                b_omega2->data[b_omega2->size[0] * i0] = F1->data[F1->size[0] * i0];
            }

            loop_ub = F2->size[1];
            for (i0 = 0; i0 < loop_ub; i0++) {
                b_omega2->data[b_omega2->size[0] * (i0 + F1->size[1])] = F2->data
                        [F2->size[0] * i0];
            }

            i0 = b_W1->size[0] * b_W1->size[1];
            b_W1->size[0] = 1;
            b_W1->size[1] = W1->size[1] + W2->size[1];
            emxEnsureCapacity_real_T(b_W1, i0);
            loop_ub = W1->size[1];
            for (i0 = 0; i0 < loop_ub; i0++) {
                b_W1->data[b_W1->size[0] * i0] = W1->data[W1->size[0] * i0];
            }

            loop_ub = W2->size[1];
            for (i0 = 0; i0 < loop_ub; i0++) {
                b_W1->data[b_W1->size[0] * (i0 + W1->size[1])] = W2->data[W2->size[0] *
                        i0];
            }

            firpm_cg(clkFIR - 1.0, b1, a, b_omega2, b_W1, ccoef);
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
            dec_int3 = 0;
            exitg1 = false;
            while ((!exitg1) && (dec_int3 < 126)) {
                b1[0] = F1->data[0];
                b1[1] = F1->data[F1->size[1] - 1];
                b1[2] = F2->data[0];
                b1[3] = F2->data[F2->size[1] - 1];
                i0 = a->size[0] * a->size[1];
                a->size[0] = 1;
                a->size[1] = A1->size[1] + A2->size[1];
                emxEnsureCapacity_real_T(a, i0);
                loop_ub = A1->size[1];
                for (i0 = 0; i0 < loop_ub; i0++) {
                    a->data[a->size[0] * i0] = A1->data[A1->size[0] * i0];
                }

                loop_ub = A2->size[1];
                for (i0 = 0; i0 < loop_ub; i0++) {
                    a->data[a->size[0] * (i0 + A1->size[1])] = A2->data[A2->size[0] * i0];
                }

                i0 = b_omega2->size[0] * b_omega2->size[1];
                b_omega2->size[0] = 1;
                b_omega2->size[1] = F1->size[1] + F2->size[1];
                emxEnsureCapacity_real_T(b_omega2, i0);
                loop_ub = F1->size[1];
                for (i0 = 0; i0 < loop_ub; i0++) {
                    b_omega2->data[b_omega2->size[0] * i0] = F1->data[F1->size[0] * i0];
                }

                loop_ub = F2->size[1];
                for (i0 = 0; i0 < loop_ub; i0++) {
                    b_omega2->data[b_omega2->size[0] * (i0 + F1->size[1])] = F2->data
                            [F2->size[0] * i0];
                }

                i0 = b_W1->size[0] * b_W1->size[1];
                b_W1->size[0] = 1;
                b_W1->size[1] = W1->size[1] + W2->size[1];
                emxEnsureCapacity_real_T(b_W1, i0);
                loop_ub = W1->size[1];
                for (i0 = 0; i0 < loop_ub; i0++) {
                    b_W1->data[b_W1->size[0] * i0] = W1->data[W1->size[0] * i0];
                }

                loop_ub = W2->size[1];
                for (i0 = 0; i0 < loop_ub; i0++) {
                    b_W1->data[b_W1->size[0] * (i0 + W1->size[1])] = W2->data[W2->size[0] *
                            i0];
                }

                b_firpm_cg(3.0 + (double)dec_int3, b1, a, b_omega2, b_W1, ccoef, &valid,
                           &sigmax);

                /*  Check if design meets specs */
                if ((sigmax < sigma) && valid) {
                    exitg1 = true;
                } else {
                    dec_int3++;
                }
            }
        }

        /*  Enable phase equalization and apply update to taps */
        if (phEQ != -1.0) {
            if (1 > fg2->size[1]) {
                i0 = 1;
                n = 1;
                hb3 = 0;
            } else {
                i0 = fg2->size[1];
                n = -1;
                hb3 = 1;
            }

            i2 = fg->size[0] * fg->size[1];
            fg->size[0] = 1;
            fg->size[1] = div_s32_floor(hb3 - i0, n) + 1;
            emxEnsureCapacity_real_T(fg, i2);
            loop_ub = div_s32_floor(hb3 - i0, n);
            for (hb3 = 0; hb3 <= loop_ub; hb3++) {
                fg->data[fg->size[0] * hb3] = 0.5 - fg2->data[(i0 + n * hb3) - 1];
            }

            if (1 > rgN->size[1]) {
                i0 = 1;
                n = 1;
                hb3 = 0;
            } else {
                i0 = rgN->size[1];
                n = -1;
                hb3 = 1;
            }

            i2 = omega->size[0] * omega->size[1];
            omega->size[0] = 1;
            omega->size[1] = div_s32_floor(hb3 - i0, n) + 1;
            emxEnsureCapacity_real_T(omega, i2);
            loop_ub = div_s32_floor(hb3 - i0, n);
            for (hb3 = 0; hb3 <= loop_ub; hb3++) {
                omega->data[omega->size[0] * hb3] = rgN->data[(i0 + n * hb3) - 1].im;
            }

            if (1 > weight->size[1]) {
                i0 = 1;
                n = 1;
                hb3 = 0;
            } else {
                i0 = weight->size[1];
                n = -1;
                hb3 = 1;
            }

            i2 = sw->size[0] * sw->size[1];
            sw->size[0] = 1;
            sw->size[1] = div_s32_floor(hb3 - i0, n) + 1;
            emxEnsureCapacity_real_T(sw, i2);
            loop_ub = div_s32_floor(hb3 - i0, n);
            for (i2 = 0; i2 <= loop_ub; i2++) {
                sw->data[sw->size[0] * i2] = weight->data[(i0 + n * i2) - 1];
            }

            if (1.0 > (8192.0 - Gstop) + 1.0) {
                loop_ub = 0;
            } else {
                loop_ub = (int)((8192.0 - Gstop) + 1.0);
            }

            i2 = wg->size[0] * wg->size[1];
            wg->size[0] = 1;
            wg->size[1] = loop_ub;
            emxEnsureCapacity_real_T(wg, i2);
            for (i2 = 0; i2 < loop_ub; i2++) {
                wg->data[wg->size[0] * i2] = fg->data[i2] * 2.0;
            }

            if ((8192.0 - Gstop) + 2.0 > fg->size[1]) {
                i2 = 0;
                dec_int3 = 0;
            } else {
                i2 = (int)((8192.0 - Gstop) + 2.0) - 1;
                dec_int3 = fg->size[1];
            }

            i3 = F4->size[0] * F4->size[1];
            F4->size[0] = 1;
            F4->size[1] = dec_int3 - i2;
            emxEnsureCapacity_real_T(F4, i3);
            loop_ub = dec_int3 - i2;
            for (dec_int3 = 0; dec_int3 < loop_ub; dec_int3++) {
                F4->data[F4->size[0] * dec_int3] = fg->data[i2 + dec_int3] * 2.0;
            }

            if (1.0 > (8192.0 - Gstop) + 1.0) {
                loop_ub = 0;
            } else {
                loop_ub = (int)((8192.0 - Gstop) + 1.0);
            }

            if ((8192.0 - Gstop) + 2.0 > omega->size[1]) {
                i2 = 0;
                dec_int3 = 0;
            } else {
                i2 = (int)((8192.0 - Gstop) + 2.0) - 1;
                dec_int3 = omega->size[1];
            }

            if (1.0 > (8192.0 - Gstop) + 1.0) {
                b_loop_ub = 0;
            } else {
                b_loop_ub = (int)((8192.0 - Gstop) + 1.0);
            }

            if ((8192.0 - Gstop) + 2.0 > div_s32_floor(hb3 - i0, n) + 1) {
                i3 = 0;
                i0 = -1;
            } else {
                i3 = (int)((8192.0 - Gstop) + 2.0) - 1;
                i0 = div_s32_floor(hb3 - i0, n);
            }

            if (int_FIR != 0.0) {
                sigma = clkFIR - 1.0;
            } else {
                sigma = (double)ccoef->size[1] - 1.0;
            }

            b1[0] = wg->data[0];
            b1[1] = wg->data[wg->size[1] - 1];
            b1[2] = F4->data[0];
            b1[3] = F4->data[F4->size[1] - 1];
            n = b_W1->size[0] * b_W1->size[1];
            b_W1->size[0] = 1;
            b_W1->size[1] = (loop_ub + dec_int3) - i2;
            emxEnsureCapacity_real_T(b_W1, n);
            for (n = 0; n < loop_ub; n++) {
                b_W1->data[b_W1->size[0] * n] = omega->data[n];
            }

            hb3 = dec_int3 - i2;
            for (n = 0; n < hb3; n++) {
                b_W1->data[b_W1->size[0] * (n + loop_ub)] = omega->data[i2 + n];
            }

            n = a->size[0] * a->size[1];
            a->size[0] = 1;
            a->size[1] = wg->size[1] + F4->size[1];
            emxEnsureCapacity_real_T(a, n);
            loop_ub = wg->size[1];
            for (n = 0; n < loop_ub; n++) {
                a->data[a->size[0] * n] = wg->data[wg->size[0] * n];
            }

            loop_ub = F4->size[1];
            for (n = 0; n < loop_ub; n++) {
                a->data[a->size[0] * (n + wg->size[1])] = F4->data[F4->size[0] * n];
            }

            n = b_omega2->size[0] * b_omega2->size[1];
            b_omega2->size[0] = 1;
            b_omega2->size[1] = ((b_loop_ub + i0) - i3) + 1;
            emxEnsureCapacity_real_T(b_omega2, n);
            for (n = 0; n < b_loop_ub; n++) {
                b_omega2->data[b_omega2->size[0] * n] = sw->data[n];
            }

            loop_ub = i0 - i3;
            for (i0 = 0; i0 <= loop_ub; i0++) {
                b_omega2->data[b_omega2->size[0] * (i0 + b_loop_ub)] = sw->data[i3 + i0];
            }

            firpm_cg(sigma, b1, b_W1, a, b_omega2, fg);
            i0 = fg->size[1];
            for (dec_int3 = 0; dec_int3 < i0; dec_int3++) {
                fg->data[dec_int3] = -fg->data[dec_int3] * mpower(-1.0, (1.0 + (double)
                                     dec_int3) - 1.0);
            }
        } else {
            for (i0 = 0; i0 < 2; i0++) {
                b2[i0] = ccoef->size[i0];
            }

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
            tap_store->data[((int)i + tap_store->size[0] * i0) - 1] = ccoef->
                    data[ccoef->size[0] * i0] + fg->data[fg->size[0] * i0];
        }

        /*  scoef ==0 when no EQ */
        determineBestFractionLength(tap_store, i, ccoef->size[1], a);
        loop_ub = a->size[1];
        for (i0 = 0; i0 < loop_ub; i0++) {
            tap_store->data[((int)i + tap_store->size[0] * i0) - 1] = a->data[a->size
                    [0] * i0];
        }

        if (b_strcmp(RxTx)) {
            if (1.0 > Gpass + 1.0) {
                loop_ub = 0;
            } else {
                loop_ub = (int)(Gpass + 1.0);
            }

            if (1 > ccoef->size[1]) {
                b_loop_ub = 0;
            } else {
                b_loop_ub = ccoef->size[1];
            }

            i0 = b_omega2->size[0] * b_omega2->size[1];
            b_omega2->size[0] = 1;
            b_omega2->size[1] = loop_ub;
            emxEnsureCapacity_real_T(b_omega2, i0);
            for (i0 = 0; i0 < loop_ub; i0++) {
                b_omega2->data[b_omega2->size[0] * i0] = omega2->data[i0];
            }

            i0 = a->size[0] * a->size[1];
            a->size[0] = 1;
            a->size[1] = b_loop_ub;
            emxEnsureCapacity_real_T(a, i0);
            for (i0 = 0; i0 < b_loop_ub; i0++) {
                a->data[a->size[0] * i0] = tap_store->data[((int)i + tap_store->size[0] *
                                           i0) - 1];
            }

            c_generateCascadedResponseRx(enables, b_omega2, converter_rate, hb1_coeff,
                                         hb2_coeff, hb3_coeff_data, hb3_coeff_size, dec_int3_coeff_data,
                                         dec_int3_coeff_size, a, c_combinedResponse);
            if (Gpass + 2.0 > omega2->size[1]) {
                i0 = 0;
                n = 0;
            } else {
                i0 = (int)(Gpass + 2.0) - 1;
                n = omega2->size[1];
            }

            if (1 > ccoef->size[1]) {
                loop_ub = 0;
            } else {
                loop_ub = ccoef->size[1];
            }

            hb3 = b_omega2->size[0] * b_omega2->size[1];
            b_omega2->size[0] = 1;
            b_omega2->size[1] = n - i0;
            emxEnsureCapacity_real_T(b_omega2, hb3);
            b_loop_ub = n - i0;
            for (n = 0; n < b_loop_ub; n++) {
                b_omega2->data[b_omega2->size[0] * n] = omega2->data[i0 + n];
            }

            i0 = a->size[0] * a->size[1];
            a->size[0] = 1;
            a->size[1] = loop_ub;
            emxEnsureCapacity_real_T(a, i0);
            for (i0 = 0; i0 < loop_ub; i0++) {
                a->data[a->size[0] * i0] = tap_store->data[((int)i + tap_store->size[0] *
                                           i0) - 1];
            }

            c_generateCascadedResponseRx(enables, b_omega2, converter_rate, hb1_coeff,
                                         hb2_coeff, hb3_coeff_data, hb3_coeff_size, dec_int3_coeff_data,
                                         dec_int3_coeff_size, a, rg1);
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
                b_omega2->data[b_omega2->size[0] * i0] = omega2->data[i0];
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
                n = 0;
            } else {
                i0 = (int)(Gpass + 2.0) - 1;
                n = omega2->size[1];
            }

            hb3 = b_omega2->size[0] * b_omega2->size[1];
            b_omega2->size[0] = 1;
            b_omega2->size[1] = n - i0;
            emxEnsureCapacity_real_T(b_omega2, hb3);
            loop_ub = n - i0;
            for (n = 0; n < loop_ub; n++) {
                b_omega2->data[b_omega2->size[0] * n] = omega2->data[i0 + n];
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
                absa = rg1->data[i0].re;
                rg1_im = rg1->data[i0].im;
                d_combinedResponse->data[i0].re = re * absa - im * rg1_im;
                d_combinedResponse->data[i0].im = re * rg1_im + im * absa;
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
                b_loop_ub = 0;
            } else {
                b_loop_ub = ccoef->size[1];
            }

            i0 = b_omega2->size[0] * b_omega2->size[1];
            b_omega2->size[0] = 1;
            b_omega2->size[1] = loop_ub;
            emxEnsureCapacity_real_T(b_omega2, i0);
            for (i0 = 0; i0 < loop_ub; i0++) {
                b_omega2->data[b_omega2->size[0] * i0] = omega2->data[i0];
            }

            i0 = a->size[0] * a->size[1];
            a->size[0] = 1;
            a->size[1] = b_loop_ub;
            emxEnsureCapacity_real_T(a, i0);
            for (i0 = 0; i0 < b_loop_ub; i0++) {
                a->data[a->size[0] * i0] = tap_store->data[((int)i + tap_store->size[0] *
                                           i0) - 1];
            }

            c_generateCascadedResponseRx(enables, b_omega2, converter_rate, hb1_coeff,
                                         hb2_coeff, hb3_coeff_data, hb3_coeff_size, dec_int3_coeff_data,
                                         dec_int3_coeff_size, a, c_combinedResponse);
            if (Gpass + 2.0 > omega2->size[1]) {
                i0 = 0;
                n = 0;
            } else {
                i0 = (int)(Gpass + 2.0) - 1;
                n = omega2->size[1];
            }

            if (1 > ccoef->size[1]) {
                loop_ub = 0;
            } else {
                loop_ub = ccoef->size[1];
            }

            hb3 = b_omega2->size[0] * b_omega2->size[1];
            b_omega2->size[0] = 1;
            b_omega2->size[1] = n - i0;
            emxEnsureCapacity_real_T(b_omega2, hb3);
            b_loop_ub = n - i0;
            for (n = 0; n < b_loop_ub; n++) {
                b_omega2->data[b_omega2->size[0] * n] = omega2->data[i0 + n];
            }

            i0 = a->size[0] * a->size[1];
            a->size[0] = 1;
            a->size[1] = loop_ub;
            emxEnsureCapacity_real_T(a, i0);
            for (i0 = 0; i0 < loop_ub; i0++) {
                a->data[a->size[0] * i0] = tap_store->data[((int)i + tap_store->size[0] *
                                           i0) - 1];
            }

            c_generateCascadedResponseRx(enables, b_omega2, converter_rate, hb1_coeff,
                                         hb2_coeff, hb3_coeff_data, hb3_coeff_size, dec_int3_coeff_data,
                                         dec_int3_coeff_size, a, rg1);
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
                b_omega2->data[b_omega2->size[0] * i0] = omega2->data[i0];
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
                n = 0;
            } else {
                i0 = (int)(Gpass + 2.0) - 1;
                n = omega2->size[1];
            }

            hb3 = b_omega2->size[0] * b_omega2->size[1];
            b_omega2->size[0] = 1;
            b_omega2->size[1] = n - i0;
            emxEnsureCapacity_real_T(b_omega2, hb3);
            loop_ub = n - i0;
            for (n = 0; n < loop_ub; n++) {
                b_omega2->data[b_omega2->size[0] * n] = omega2->data[i0 + n];
            }

            d_analogresp(b_omega2, converter_rate, b1_data, b1_size, a1, b2_data,
                         b2_size, a2, r0);
            i0 = d_combinedResponse->size[0] * d_combinedResponse->size[1];
            d_combinedResponse->size[0] = 1;
            d_combinedResponse->size[1] = rg1->size[1];
            emxEnsureCapacity_creal_T(d_combinedResponse, i0);
            loop_ub = rg1->size[0] * rg1->size[1];
            for (i0 = 0; i0 < loop_ub; i0++) {
                absa = rg1->data[i0].re;
                rg1_im = rg1->data[i0].im;
                re = r0->data[i0].re;
                im = r0->data[i0].im;
                d_combinedResponse->data[i0].re = absa * re - rg1_im * im;
                d_combinedResponse->data[i0].im = absa * im + rg1_im * re;
            }

            b_abs(d_combinedResponse, omega);
        }

        /*  quantitative values about actual passband and stopband */
        hb3 = 1;
        n = fg->size[1];
        sigma = fg->data[0];
        if (fg->size[1] > 1) {
            if (rtIsNaN(fg->data[0])) {
                dec_int3 = 2;
                exitg1 = false;
                while ((!exitg1) && (dec_int3 <= n)) {
                    hb3 = dec_int3;
                    if (!rtIsNaN(fg->data[dec_int3 - 1])) {
                        sigma = fg->data[dec_int3 - 1];
                        exitg1 = true;
                    } else {
                        dec_int3++;
                    }
                }
            }

            if (hb3 < fg->size[1]) {
                while (hb3 + 1 <= n) {
                    if (fg->data[hb3] > sigma) {
                        sigma = fg->data[hb3];
                    }

                    hb3++;
                }
            }
        }

        hb3 = 1;
        n = fg->size[1];
        sigmax = fg->data[0];
        if (fg->size[1] > 1) {
            if (rtIsNaN(fg->data[0])) {
                dec_int3 = 2;
                exitg1 = false;
                while ((!exitg1) && (dec_int3 <= n)) {
                    hb3 = dec_int3;
                    if (!rtIsNaN(fg->data[dec_int3 - 1])) {
                        sigmax = fg->data[dec_int3 - 1];
                        exitg1 = true;
                    } else {
                        dec_int3++;
                    }
                }
            }

            if (hb3 < fg->size[1]) {
                while (hb3 + 1 <= n) {
                    if (fg->data[hb3] < sigmax) {
                        sigmax = fg->data[hb3];
                    }

                    hb3++;
                }
            }
        }

        Apass_actual_vector->data[(int)i - 1] = mag2db(sigma) - mag2db(sigmax);
        hb3 = 1;
        n = omega->size[1];
        sigma = omega->data[0];
        if (omega->size[1] > 1) {
            if (rtIsNaN(omega->data[0])) {
                dec_int3 = 2;
                exitg1 = false;
                while ((!exitg1) && (dec_int3 <= n)) {
                    hb3 = dec_int3;
                    if (!rtIsNaN(omega->data[dec_int3 - 1])) {
                        sigma = omega->data[dec_int3 - 1];
                        exitg1 = true;
                    } else {
                        dec_int3++;
                    }
                }
            }

            if (hb3 < omega->size[1]) {
                while (hb3 + 1 <= n) {
                    if (omega->data[hb3] > sigma) {
                        sigma = omega->data[hb3];
                    }

                    hb3++;
                }
            }
        }

        Astop_actual_vector->data[(int)i - 1] = -mag2db(sigma);
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
                fg->data[fg->size[0] * i0] = tap_store->data[tap_store->size[0] * i0];
            }

            /* Apass_actual = Apass_actual_vector(1); */
            /* Astop_actual = Astop_actual_vector(1); */
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
                fg->data[fg->size[0] * i0] = tap_store->data[tap_store->size[0] * i0];
            }

            /* Apass_actual = Apass_actual_vector(1); */
            /* Astop_actual = Astop_actual_vector(1); */
            exitg2 = 1;
        } else if ((Apass_actual_vector->data[(int)i - 1] > Apass) ||
                   (Astop_actual_vector->data[(int)i - 1] < Astop)) {
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
                fg->data[fg->size[0] * i0] = tap_store->data[((int)i + tap_store->size[0]
                                             * i0) - 2];
            }

            /* Apass_actual = Apass_actual_vector(i-1); */
            /* Astop_actual = Astop_actual_vector(i-1); */
            exitg2 = 1;
        } else {
            clkFIR -= 16.0;
            i++;
        }
    } while (exitg2 == 0);

    emxFree_real_T(&b_omega2);
    emxFree_creal_T(&d_combinedResponse);
    emxFree_real_T(&b_W1);
    emxFree_creal_T(&r0);
    emxFree_creal_T(&c_combinedResponse);
    emxFree_real_T(&F4);
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
    emxFree_real_T(&wg);
    emxFree_creal_T(&rgN);
    emxFree_real_T(&omega2);
    emxFree_real_T(&fg2);
    emxFree_creal_T(&rg1);
    emxFree_real_T(&omega);
    emxFree_creal_T(&a2);
    emxFree_creal_T(&a1);
    if (c_strcmp(RxTx)) {
        if ((int_FIR == 1.0) && (FIR == 2.0)) {
            if (rt_remd_snf(fg->size[1], 32.0) != 0.0) {
                i0 = a->size[0] * a->size[1];
                a->size[0] = 1;
                a->size[1] = 16 + fg->size[1];
                emxEnsureCapacity_real_T(a, i0);
                for (i0 = 0; i0 < 8; i0++) {
                    a->data[a->size[0] * i0] = 0.0;
                }

                loop_ub = fg->size[1];
                for (i0 = 0; i0 < loop_ub; i0++) {
                    a->data[a->size[0] * (i0 + 8)] = fg->data[fg->size[0] * i0];
                }

                for (i0 = 0; i0 < 8; i0++) {
                    a->data[a->size[0] * ((i0 + fg->size[1]) + 8)] = 0.0;
                }

                i0 = fg->size[0] * fg->size[1];
                fg->size[0] = 1;
                fg->size[1] = a->size[1];
                emxEnsureCapacity_real_T(fg, i0);
                loop_ub = a->size[1];
                for (i0 = 0; i0 < loop_ub; i0++) {
                    fg->data[fg->size[0] * i0] = a->data[a->size[0] * i0];
                }
            }
        } else {
            if ((int_FIR == 1.0) && (FIR == 4.0) && (rt_remd_snf(fg->size[1], 64.0) !=
                    0.0)) {
                sigma = (ceil((double)fg->size[1] / 64.0) * 64.0 - (double)fg->size[1]) /
                        2.0;
                i0 = a->size[0] * a->size[1];
                a->size[0] = 1;
                a->size[1] = ((int)sigma + fg->size[1]) + (int)sigma;
                emxEnsureCapacity_real_T(a, i0);
                loop_ub = (int)sigma;
                for (i0 = 0; i0 < loop_ub; i0++) {
                    a->data[a->size[0] * i0] = 0.0;
                }

                loop_ub = fg->size[1];
                for (i0 = 0; i0 < loop_ub; i0++) {
                    a->data[a->size[0] * (i0 + (int)sigma)] = fg->data[fg->size[0] * i0];
                }

                loop_ub = (int)sigma;
                for (i0 = 0; i0 < loop_ub; i0++) {
                    a->data[a->size[0] * ((i0 + (int)sigma) + fg->size[1])] = 0.0;
                }

                i0 = fg->size[0] * fg->size[1];
                fg->size[0] = 1;
                fg->size[1] = a->size[1];
                emxEnsureCapacity_real_T(fg, i0);
                loop_ub = a->size[1];
                for (i0 = 0; i0 < loop_ub; i0++) {
                    fg->data[fg->size[0] * i0] = a->data[a->size[0] * i0];
                }
            }
        }
    }

    emxFree_real_T(&a);

    /*  There will always be 128 taps output */
    memset(&firTapsPreScale[0], 0, sizeof(double) << 7);
    loop_ub = fg->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
        firTapsPreScale[i0] = fg->data[fg->size[0] * i0];
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
    hb3 = 1;
    sigma = firTapsPreScale[0];
    if (rtIsNaN(firTapsPreScale[0])) {
        dec_int3 = 2;
        exitg1 = false;
        while ((!exitg1) && (dec_int3 < 129)) {
            hb3 = dec_int3;
            if (!rtIsNaN(firTapsPreScale[dec_int3 - 1])) {
                sigma = firTapsPreScale[dec_int3 - 1];
                exitg1 = true;
            } else {
                dec_int3++;
            }
        }
    }

    if (hb3 < 128) {
        while (hb3 + 1 < 129) {
            if (firTapsPreScale[hb3] > sigma) {
                sigma = firTapsPreScale[hb3];
            }

            hb3++;
        }
    }

    sigma = b_log2(sigma);
    sigma = ceil(sigma);
    switch ((int)(1.0 + sigma)) {
    case 2:
        hb3 = 6;
        break;

    case 1:
        hb3 = 0;
        break;

    case 0:
        hb3 = -6;
        break;

    default:
        hb3 = -12;
        break;
    }

    if (b_strcmp(RxTx)) {
        if (1.0 + sigma > 2.0) {
            hb3 = 6;
        }
    } else {
        if (FIR == 2.0) {
            hb3 += 6;
        } else {
            if (FIR == 4.0) {
                hb3 += 12;
            }
        }

        if (hb3 > 0) {
            hb3 = 0;
        } else {
            if (hb3 < -6) {
                hb3 = -6;
            }
        }
    }

    /*  Scale taps */
    memcpy(&b_firTapsPreScale[0], &firTapsPreScale[0], sizeof(double) << 7);
    b_determineBestFractionLength(b_firTapsPreScale, firTapsPreScale);
    sigma = mpower(2.0, 16.0 - (1.0 + sigma));
    for (i0 = 0; i0 < 128; i0++) {
        sigmax = rt_roundd_snf(firTapsPreScale[i0] * sigma);
        if (sigmax < 32768.0) {
            if (sigmax >= -32768.0) {
                i4 = (short)sigmax;
            } else {
                i4 = MIN_int16_T;
            }
        } else if (sigmax >= 32768.0) {
            i4 = MAX_int16_T;
        } else {
            i4 = 0;
        }

        outputTaps[i0] = i4;
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
    *filterGain = hb3;
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
