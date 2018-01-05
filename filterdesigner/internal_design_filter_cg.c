/*
 * Sponsored Third Party Support License -- for use only to support
 * products interfaced to MathWorks software under terms specified in your
 * company's restricted use license agreement.
 * File: internal_design_filter_cg.c
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 05-Jan-2018 14:44:11
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "internal_design_filter_cg.h"
#include "internal_design_filter_cg_emxutil.h"

/* Type Definitions */
#include <stdio.h>
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
static void b_power(const emxArray_real_T *a, emxArray_real_T *y);
static void b_rdivide(const emxArray_creal_T *x, const emxArray_creal_T *y,
                      emxArray_creal_T *z);
static void b_sinc(emxArray_real_T *x);
static void b_sqrt(creal_T *x);
static boolean_T b_strcmp(const char a[2]);
static double b_sum(const emxArray_real_T *x);
static void b_us(const double o[7], double u[7]);
static void b_xscal(int n, const creal_T a, emxArray_creal_T *x, int ix0, int
                    incx);
static void b_xzlartg(const creal_T f, const creal_T g, double *cs,
                      creal_T *sn);
static void butter_cg(double Wn, double num[2], creal_T den_data[], int
                      den_size[2]);
static void c_abs(const double x_data[], const int x_size[2], double y_data[],
                  int y_size[2]);
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
static double c_sum(const double x[128]);
static void c_us(const double o_data[], const int o_size[2], double u_data[],
                 int u_size[2]);
static int cfprintf(const char * varargin_1);
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
static void determineBestFractionLength(const double tap_store_data[], const int
                                        tap_store_size[2], double i, double M, emxArray_real_T *taps);
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
static void j_polyval(const double p_data[], const int p_size[2], const
                      emxArray_creal_T *x, emxArray_creal_T *y);
static void k_freqz_cg(const emxArray_real_T *w, double Fs,
                       emxArray_creal_T *hh);
static void l_freqz_cg(const double b[15], const emxArray_real_T *w, double Fs,
                       emxArray_creal_T *hh);
static void lp2lp_cg(const emxArray_creal_T *a, const emxArray_real_T *b, double
                     wo, emxArray_creal_T *at, emxArray_real_T *bt, double *dt);
static void m_freqz_cg(const double b[7], const emxArray_real_T *w, double Fs,
                       emxArray_creal_T *hh);
static double mag2db(double y);
static void n_freqz_cg(const double b_data[], const int b_size[2], const
                       emxArray_real_T *w, double Fs, emxArray_creal_T *hh);
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
static void upsample(const double x_data[], const int x_size[2], double N,
                     double y_data[], int y_size[2]);
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

    double a[2048];
    double dv14[2048];
    double y[2048];
    static creal_T dcv1[2048];
    double abc_im;
    double a_im;
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
            a[kstr] = f[kstr] / Fconverter;
        }

        sinc(a);
        for (kstr = 0; kstr < 2048; kstr++) {
            dv14[kstr] = 6.2831853071795862 * f[kstr];
        }

        freqs_cg(b1_data, b1_size, a1, dv14, abc);
        for (kstr = 0; kstr < 2048; kstr++) {
            dv14[kstr] = 6.2831853071795862 * f[kstr];
        }

        freqs_cg(b2_data, b2_size, a2, dv14, dcv1);
        for (kstr = 0; kstr < 2048; kstr++) {
            abc_im = a[kstr] * abc[kstr].re;
            a_im = a[kstr] * abc[kstr].im;
            abc[kstr].re = abc_im * dcv1[kstr].re - a_im * dcv1[kstr].im;
            abc[kstr].im = abc_im * dcv1[kstr].im + a_im * dcv1[kstr].re;
        }
        break;

    case 1:
        for (kstr = 0; kstr < 2048; kstr++) {
            a[kstr] = f[kstr] / Fconverter;
        }

        sinc(a);
        for (kstr = 0; kstr < 2048; kstr++) {
            y[kstr] = rt_powd_snf(a[kstr], 3.0);
            dv14[kstr] = 6.2831853071795862 * f[kstr];
        }

        freqs_cg(b1_data, b1_size, a1, dv14, abc);
        for (kstr = 0; kstr < 2048; kstr++) {
            dv14[kstr] = 6.2831853071795862 * f[kstr];
        }

        freqs_cg(b2_data, b2_size, a2, dv14, dcv1);
        for (kstr = 0; kstr < 2048; kstr++) {
            abc_im = abc[kstr].re * dcv1[kstr].im + abc[kstr].im * dcv1[kstr].re;
            abc[kstr].re = y[kstr] * (abc[kstr].re * dcv1[kstr].re - abc[kstr].im *
                                      dcv1[kstr].im);
            abc[kstr].im = y[kstr] * abc_im;
        }
        break;

    default:
        /*  Default to Rx */
        for (kstr = 0; kstr < 2048; kstr++) {
            a[kstr] = f[kstr] / Fconverter;
        }

        sinc(a);
        for (kstr = 0; kstr < 2048; kstr++) {
            y[kstr] = rt_powd_snf(a[kstr], 3.0);
            dv14[kstr] = 6.2831853071795862 * f[kstr];
        }

        freqs_cg(b1_data, b1_size, a1, dv14, abc);
        for (kstr = 0; kstr < 2048; kstr++) {
            dv14[kstr] = 6.2831853071795862 * f[kstr];
        }

        freqs_cg(b2_data, b2_size, a2, dv14, dcv1);
        for (kstr = 0; kstr < 2048; kstr++) {
            abc_im = abc[kstr].re * dcv1[kstr].im + abc[kstr].im * dcv1[kstr].re;
            abc[kstr].re = y[kstr] * (abc[kstr].re * dcv1[kstr].re - abc[kstr].im *
                                      dcv1[kstr].im);
            abc[kstr].im = y[kstr] * abc_im;
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
    p = false;
    for (k = 0; k + 1 <= nx; k++) {
        if (p || (rtIsInf(x->data[k].re) || rtIsInf(x->data[k].im)) || (rtIsNaN
                (x->data[k].re) || rtIsNaN(x->data[k].im))) {
            p = true;
        } else {
            p = false;
        }
    }

    return p;
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
    emxEnsureCapacity((emxArray__common *)y, k, sizeof(double));
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
    int exitg1;
    static const char cv46[2] = { 'T', 'x' };

    emxArray_creal_T *r19;
    emxArray_real_T *r20;
    emxArray_real_T *r21;
    static const char cv47[2] = { 'R', 'x' };

    emxArray_real_T *r22;
    int loop_ub;
    emxArray_real_T *r23;
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

    emxInit_creal_T(&r19, 2);
    emxInit_real_T(&r20, 2);
    emxInit_real_T(&r21, 2);
    switch (kstr) {
    case 0:
        emxInit_real_T(&r22, 2);
        rdivide(f, Fconverter, r21);
        b_sinc(r21);
        kstr = r22->size[0] * r22->size[1];
        r22->size[0] = 1;
        r22->size[1] = f->size[1];
        emxEnsureCapacity((emxArray__common *)r22, kstr, sizeof(double));
        loop_ub = f->size[0] * f->size[1];
        for (kstr = 0; kstr < loop_ub; kstr++) {
            r22->data[kstr] = 6.2831853071795862 * f->data[kstr];
        }

        emxInit_real_T(&r23, 2);
        b_freqs_cg(b1_data, b1_size, a1, r22, abc);
        kstr = r23->size[0] * r23->size[1];
        r23->size[0] = 1;
        r23->size[1] = f->size[1];
        emxEnsureCapacity((emxArray__common *)r23, kstr, sizeof(double));
        loop_ub = f->size[0] * f->size[1];
        emxFree_real_T(&r22);
        for (kstr = 0; kstr < loop_ub; kstr++) {
            r23->data[kstr] = 6.2831853071795862 * f->data[kstr];
        }

        b_freqs_cg(b2_data, b2_size, a2, r23, r19);
        kstr = abc->size[0] * abc->size[1];
        abc->size[0] = 1;
        abc->size[1] = r21->size[1];
        emxEnsureCapacity((emxArray__common *)abc, kstr, sizeof(creal_T));
        loop_ub = r21->size[0] * r21->size[1];
        emxFree_real_T(&r23);
        for (kstr = 0; kstr < loop_ub; kstr++) {
            abc_re = r21->data[kstr] * abc->data[kstr].re;
            abc_im = r21->data[kstr] * abc->data[kstr].im;
            re = r19->data[kstr].re;
            im = r19->data[kstr].im;
            abc->data[kstr].re = abc_re * re - abc_im * im;
            abc->data[kstr].im = abc_re * im + abc_im * re;
        }
        break;

    case 1:
        emxInit_real_T(&r22, 2);
        kstr = r22->size[0] * r22->size[1];
        r22->size[0] = 1;
        r22->size[1] = f->size[1];
        emxEnsureCapacity((emxArray__common *)r22, kstr, sizeof(double));
        loop_ub = f->size[0] * f->size[1];
        for (kstr = 0; kstr < loop_ub; kstr++) {
            r22->data[kstr] = 6.2831853071795862 * f->data[kstr];
        }

        emxInit_real_T(&r23, 2);
        b_freqs_cg(b1_data, b1_size, a1, r22, abc);
        kstr = r23->size[0] * r23->size[1];
        r23->size[0] = 1;
        r23->size[1] = f->size[1];
        emxEnsureCapacity((emxArray__common *)r23, kstr, sizeof(double));
        loop_ub = f->size[0] * f->size[1];
        emxFree_real_T(&r22);
        for (kstr = 0; kstr < loop_ub; kstr++) {
            r23->data[kstr] = 6.2831853071795862 * f->data[kstr];
        }

        b_freqs_cg(b2_data, b2_size, a2, r23, r19);
        rdivide(f, Fconverter, r21);
        b_sinc(r21);
        b_power(r21, r20);
        kstr = abc->size[0] * abc->size[1];
        abc->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)abc, kstr, sizeof(creal_T));
        kstr = abc->size[0];
        loop_ub = abc->size[1];
        loop_ub *= kstr;
        emxFree_real_T(&r23);
        for (kstr = 0; kstr < loop_ub; kstr++) {
            abc_re = abc->data[kstr].re * r19->data[kstr].re - abc->data[kstr].im *
                     r19->data[kstr].im;
            abc_im = abc->data[kstr].re * r19->data[kstr].im + abc->data[kstr].im *
                     r19->data[kstr].re;
            abc->data[kstr].re = r20->data[kstr] * abc_re;
            abc->data[kstr].im = r20->data[kstr] * abc_im;
        }
        break;

    default:
        emxInit_real_T(&r22, 2);

        /*  Default to Rx */
        kstr = r22->size[0] * r22->size[1];
        r22->size[0] = 1;
        r22->size[1] = f->size[1];
        emxEnsureCapacity((emxArray__common *)r22, kstr, sizeof(double));
        loop_ub = f->size[0] * f->size[1];
        for (kstr = 0; kstr < loop_ub; kstr++) {
            r22->data[kstr] = 6.2831853071795862 * f->data[kstr];
        }

        emxInit_real_T(&r23, 2);
        b_freqs_cg(b1_data, b1_size, a1, r22, abc);
        kstr = r23->size[0] * r23->size[1];
        r23->size[0] = 1;
        r23->size[1] = f->size[1];
        emxEnsureCapacity((emxArray__common *)r23, kstr, sizeof(double));
        loop_ub = f->size[0] * f->size[1];
        emxFree_real_T(&r22);
        for (kstr = 0; kstr < loop_ub; kstr++) {
            r23->data[kstr] = 6.2831853071795862 * f->data[kstr];
        }

        b_freqs_cg(b2_data, b2_size, a2, r23, r19);
        rdivide(f, Fconverter, r21);
        b_sinc(r21);
        b_power(r21, r20);
        kstr = abc->size[0] * abc->size[1];
        abc->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)abc, kstr, sizeof(creal_T));
        kstr = abc->size[0];
        loop_ub = abc->size[1];
        loop_ub *= kstr;
        emxFree_real_T(&r23);
        for (kstr = 0; kstr < loop_ub; kstr++) {
            abc_re = abc->data[kstr].re * r19->data[kstr].re - abc->data[kstr].im *
                     r19->data[kstr].im;
            abc_im = abc->data[kstr].re * r19->data[kstr].im + abc->data[kstr].im *
                     r19->data[kstr].re;
            abc->data[kstr].re = r20->data[kstr] * abc_re;
            abc->data[kstr].im = r20->data[kstr] * abc_im;
        }
        break;
    }

    emxFree_real_T(&r21);
    emxFree_real_T(&r20);
    emxFree_creal_T(&r19);
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
    double dv15[128];
    double mtmp;
    double e[16];
    double v;
    double dv16[128];
    short i58;
    double dv17[128];
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
    int itmp;
    int ix;
    boolean_T exitg1;
    memset(&r[0], 0, sizeof(double) << 11);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        mtmp = tap_store[ixstart] * 2.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i58 = (short)mtmp;
            } else {
                i58 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i58 = MAX_int16_T;
        } else {
            i58 = 0;
        }

        r[ixstart << 4] = (double)i58 * 0.5;
        b_r[ixstart] = r[ixstart << 4] - tap_store[ixstart];
        mtmp = tap_store[ixstart] * 4.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i58 = (short)mtmp;
            } else {
                i58 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i58 = MAX_int16_T;
        } else {
            i58 = 0;
        }

        r[1 + (ixstart << 4)] = (double)i58 * 0.25;
    }

    d_abs(b_r, dv15);
    e[0] = c_sum(dv15);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        b_r[ixstart] = r[1 + (ixstart << 4)] - tap_store[ixstart];
        mtmp = tap_store[ixstart] * 8.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i58 = (short)mtmp;
            } else {
                i58 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i58 = MAX_int16_T;
        } else {
            i58 = 0;
        }

        r[2 + (ixstart << 4)] = (double)i58 * 0.125;
    }

    d_abs(b_r, dv16);
    e[1] = c_sum(dv16);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        b_r[ixstart] = r[2 + (ixstart << 4)] - tap_store[ixstart];
        mtmp = tap_store[ixstart] * 16.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i58 = (short)mtmp;
            } else {
                i58 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i58 = MAX_int16_T;
        } else {
            i58 = 0;
        }

        r[3 + (ixstart << 4)] = (double)i58 * 0.0625;
    }

    d_abs(b_r, dv17);
    e[2] = c_sum(dv17);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        b_r[ixstart] = r[3 + (ixstart << 4)] - tap_store[ixstart];
        mtmp = tap_store[ixstart] * 32.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i58 = (short)mtmp;
            } else {
                i58 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i58 = MAX_int16_T;
        } else {
            i58 = 0;
        }

        r[4 + (ixstart << 4)] = (double)i58 * 0.03125;
    }

    d_abs(b_r, dv18);
    e[3] = c_sum(dv18);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        b_r[ixstart] = r[4 + (ixstart << 4)] - tap_store[ixstart];
        mtmp = tap_store[ixstart] * 64.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i58 = (short)mtmp;
            } else {
                i58 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i58 = MAX_int16_T;
        } else {
            i58 = 0;
        }

        r[5 + (ixstart << 4)] = (double)i58 * 0.015625;
    }

    d_abs(b_r, dv19);
    e[4] = c_sum(dv19);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        b_r[ixstart] = r[5 + (ixstart << 4)] - tap_store[ixstart];
        mtmp = tap_store[ixstart] * 128.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i58 = (short)mtmp;
            } else {
                i58 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i58 = MAX_int16_T;
        } else {
            i58 = 0;
        }

        r[6 + (ixstart << 4)] = (double)i58 * 0.0078125;
    }

    d_abs(b_r, dv20);
    e[5] = c_sum(dv20);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        b_r[ixstart] = r[6 + (ixstart << 4)] - tap_store[ixstart];
        mtmp = tap_store[ixstart] * 256.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i58 = (short)mtmp;
            } else {
                i58 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i58 = MAX_int16_T;
        } else {
            i58 = 0;
        }

        r[7 + (ixstart << 4)] = (double)i58 * 0.00390625;
    }

    d_abs(b_r, dv21);
    e[6] = c_sum(dv21);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        b_r[ixstart] = r[7 + (ixstart << 4)] - tap_store[ixstart];
        mtmp = tap_store[ixstart] * 512.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i58 = (short)mtmp;
            } else {
                i58 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i58 = MAX_int16_T;
        } else {
            i58 = 0;
        }

        r[8 + (ixstart << 4)] = (double)i58 * 0.001953125;
    }

    d_abs(b_r, dv22);
    e[7] = c_sum(dv22);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        b_r[ixstart] = r[8 + (ixstart << 4)] - tap_store[ixstart];
        mtmp = tap_store[ixstart] * 1024.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i58 = (short)mtmp;
            } else {
                i58 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i58 = MAX_int16_T;
        } else {
            i58 = 0;
        }

        r[9 + (ixstart << 4)] = (double)i58 * 0.0009765625;
    }

    d_abs(b_r, dv23);
    e[8] = c_sum(dv23);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        b_r[ixstart] = r[9 + (ixstart << 4)] - tap_store[ixstart];
        mtmp = tap_store[ixstart] * 2048.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i58 = (short)mtmp;
            } else {
                i58 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i58 = MAX_int16_T;
        } else {
            i58 = 0;
        }

        r[10 + (ixstart << 4)] = (double)i58 * 0.00048828125;
    }

    d_abs(b_r, dv24);
    e[9] = c_sum(dv24);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        b_r[ixstart] = r[10 + (ixstart << 4)] - tap_store[ixstart];
        mtmp = tap_store[ixstart] * 4096.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i58 = (short)mtmp;
            } else {
                i58 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i58 = MAX_int16_T;
        } else {
            i58 = 0;
        }

        r[11 + (ixstart << 4)] = (double)i58 * 0.000244140625;
    }

    d_abs(b_r, dv25);
    e[10] = c_sum(dv25);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        b_r[ixstart] = r[11 + (ixstart << 4)] - tap_store[ixstart];
        mtmp = tap_store[ixstart] * 8192.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i58 = (short)mtmp;
            } else {
                i58 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i58 = MAX_int16_T;
        } else {
            i58 = 0;
        }

        r[12 + (ixstart << 4)] = (double)i58 * 0.0001220703125;
    }

    d_abs(b_r, dv26);
    e[11] = c_sum(dv26);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        b_r[ixstart] = r[12 + (ixstart << 4)] - tap_store[ixstart];
        mtmp = tap_store[ixstart] * 16384.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i58 = (short)mtmp;
            } else {
                i58 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i58 = MAX_int16_T;
        } else {
            i58 = 0;
        }

        r[13 + (ixstart << 4)] = (double)i58 * 6.103515625E-5;
    }

    d_abs(b_r, dv27);
    e[12] = c_sum(dv27);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        b_r[ixstart] = r[13 + (ixstart << 4)] - tap_store[ixstart];
        mtmp = tap_store[ixstart] * 32768.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i58 = (short)mtmp;
            } else {
                i58 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i58 = MAX_int16_T;
        } else {
            i58 = 0;
        }

        r[14 + (ixstart << 4)] = (double)i58 * 3.0517578125E-5;
    }

    d_abs(b_r, dv28);
    e[13] = c_sum(dv28);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        b_r[ixstart] = r[14 + (ixstart << 4)] - tap_store[ixstart];
        mtmp = tap_store[ixstart] * 65536.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i58 = (short)mtmp;
            } else {
                i58 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i58 = MAX_int16_T;
        } else {
            i58 = 0;
        }

        r[15 + (ixstart << 4)] = (double)i58 * 1.52587890625E-5;
    }

    d_abs(b_r, dv29);
    e[14] = c_sum(dv29);
    for (ixstart = 0; ixstart < 128; ixstart++) {
        b_r[ixstart] = r[15 + (ixstart << 4)] - tap_store[ixstart];
    }

    d_abs(b_r, dv30);
    e[15] = c_sum(dv30);
    ixstart = 1;
    mtmp = e[0];
    itmp = 0;
    if (rtIsNaN(e[0])) {
        ix = 2;
        exitg1 = false;
        while ((!exitg1) && (ix < 17)) {
            ixstart = ix;
            if (!rtIsNaN(e[ix - 1])) {
                mtmp = e[ix - 1];
                itmp = ix - 1;
                exitg1 = true;
            } else {
                ix++;
            }
        }
    }

    if (ixstart < 16) {
        while (ixstart + 1 < 17) {
            if (e[ixstart] < mtmp) {
                mtmp = e[ixstart];
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
    double bim;
    double digw[2048];
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
        bim = 6.2831853071795862 * options->w[i64] / options->Fs;
        dcv2[i64].re = bim * 0.0;
        dcv2[i64].im = bim;
        digw[i64] = bim;
    }

    b_exp(dcv2);
    polyval(b, dcv2, h);
    for (i64 = 0; i64 < 2048; i64++) {
        dcv2[i64].re = 14.0 * (digw[i64] * 0.0);
        dcv2[i64].im = 14.0 * digw[i64];
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
            bim = fabs(dcv2[i64].im);
            if (brm > bim) {
                bim = dcv2[i64].im / dcv2[i64].re;
                d = dcv2[i64].re + bim * dcv2[i64].im;
                h[i64].re = (h[i64].re + bim * h[i64].im) / d;
                h[i64].im = (h[i64].im - bim * h_re) / d;
            } else if (bim == brm) {
                if (dcv2[i64].re > 0.0) {
                    bim = 0.5;
                } else {
                    bim = -0.5;
                }

                if (dcv2[i64].im > 0.0) {
                    d = 0.5;
                } else {
                    d = -0.5;
                }

                h[i64].re = (h[i64].re * bim + h[i64].im * d) / brm;
                h[i64].im = (h[i64].im * bim - h_re * d) / brm;
            } else {
                bim = dcv2[i64].re / dcv2[i64].im;
                d = dcv2[i64].im + bim * dcv2[i64].re;
                h[i64].re = (bim * h[i64].re + h[i64].im) / d;
                h[i64].im = (bim * h[i64].im - h_re) / d;
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
    emxArray_real_T *r30;
    int i53;
    emxArray_real_T *b_h;
    double b_ff[4];
    double x;
    boolean_T b_valid;
    int h_idx_0;
    double d2;
    int i54;
    emxArray_real_T *c_h;
    int i55;
    int loop_ub;
    emxArray_real_T *d_h;
    emxInit_real_T(&grid, 2);
    emxInit_real_T(&des, 2);
    emxInit_real_T(&wt, 2);
    emxInit_real_T(&r30, 2);

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
    rdivide(grid, 2.0, r30);
    emxFree_real_T(&grid);
    for (i53 = 0; i53 < 4; i53++) {
        b_ff[i53] = ff[i53] / 2.0;
    }

    emxInit_real_T(&b_h, 2);
    remezm(order + 1.0, b_ff, r30, des, wt, b_h, &x, &b_valid);
    h_idx_0 = b_h->size[0] * b_h->size[1];
    i53 = h->size[0] * h->size[1];
    h->size[0] = 1;
    h->size[1] = h_idx_0;
    emxEnsureCapacity((emxArray__common *)h, i53, sizeof(double));
    emxFree_real_T(&r30);
    emxFree_real_T(&wt);
    emxFree_real_T(&des);
    for (i53 = 0; i53 < h_idx_0; i53++) {
        h->data[h->size[0] * i53] = b_h->data[i53];
    }

    emxFree_real_T(&b_h);

    /*  make it a row */
    d2 = (double)h->size[1] - rt_remd_snf(order + 1.0, 2.0);
    if (1.0 > d2) {
        i53 = 1;
        h_idx_0 = 1;
        i54 = 0;
    } else {
        i53 = (int)d2;
        h_idx_0 = -1;
        i54 = 1;
    }

    emxInit_real_T(&c_h, 2);
    i55 = c_h->size[0] * c_h->size[1];
    c_h->size[0] = 1;
    c_h->size[1] = (h->size[1] + div_s32_floor(i54 - i53, h_idx_0)) + 1;
    emxEnsureCapacity((emxArray__common *)c_h, i55, sizeof(double));
    loop_ub = h->size[1];
    for (i55 = 0; i55 < loop_ub; i55++) {
        c_h->data[c_h->size[0] * i55] = h->data[h->size[0] * i55];
    }

    loop_ub = div_s32_floor(i54 - i53, h_idx_0);
    for (i54 = 0; i54 <= loop_ub; i54++) {
        c_h->data[c_h->size[0] * (i54 + h->size[1])] = h->data[(i53 + h_idx_0 * i54)
                - 1];
    }

    i53 = h->size[0] * h->size[1];
    h->size[0] = 1;
    h->size[1] = c_h->size[1];
    emxEnsureCapacity((emxArray__common *)h, i53, sizeof(double));
    loop_ub = c_h->size[1];
    for (i53 = 0; i53 < loop_ub; i53++) {
        h->data[h->size[0] * i53] = c_h->data[c_h->size[0] * i53];
    }

    emxFree_real_T(&c_h);
    if (1 > h->size[1]) {
        i53 = 1;
        h_idx_0 = 1;
        i54 = 0;
    } else {
        i53 = h->size[1];
        h_idx_0 = -1;
        i54 = 1;
    }

    emxInit_real_T(&d_h, 2);
    i55 = d_h->size[0] * d_h->size[1];
    d_h->size[0] = 1;
    d_h->size[1] = div_s32_floor(i54 - i53, h_idx_0) + 1;
    emxEnsureCapacity((emxArray__common *)d_h, i55, sizeof(double));
    loop_ub = div_s32_floor(i54 - i53, h_idx_0);
    for (i54 = 0; i54 <= loop_ub; i54++) {
        d_h->data[d_h->size[0] * i54] = h->data[(i53 + h_idx_0 * i54) - 1];
    }

    i53 = h->size[0] * h->size[1];
    h->size[0] = 1;
    h->size[1] = d_h->size[1];
    emxEnsureCapacity((emxArray__common *)h, i53, sizeof(double));
    loop_ub = d_h->size[1];
    for (i53 = 0; i53 < loop_ub; i53++) {
        h->data[h->size[0] * i53] = d_h->data[d_h->size[0] * i53];
    }

    emxFree_real_T(&d_h);
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
    int i41;
    int loop_ub;
    emxArray_creal_T *y;
    boolean_T b11;
    int k;
    double a_re;
    double a_im;
    double s_re;
    double s_im;
    emxInit_creal_T(&s, 2);
    emxInit_creal_T(&b_a, 2);
    removeTrailingZero(b_data, b_size, a, b_b_data, b_b_size, b_a);
    i41 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = w->size[1];
    emxEnsureCapacity((emxArray__common *)s, i41, sizeof(creal_T));
    loop_ub = w->size[0] * w->size[1];
    for (i41 = 0; i41 < loop_ub; i41++) {
        s->data[i41].re = w->data[i41] * 0.0;
        s->data[i41].im = w->data[i41];
    }

    emxInit_creal_T(&y, 2);
    i41 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = s->size[1];
    emxEnsureCapacity((emxArray__common *)y, i41, sizeof(creal_T));
    if ((y->size[1] == 0) || (b_a->size[1] == 0)) {
        b11 = true;
    } else {
        b11 = false;
    }

    if (!b11) {
        i41 = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)y, i41, sizeof(creal_T));
        loop_ub = y->size[1];
        for (i41 = 0; i41 < loop_ub; i41++) {
            y->data[y->size[0] * i41] = b_a->data[0];
        }

        for (k = 0; k <= b_a->size[1] - 2; k++) {
            i41 = y->size[0] * y->size[1];
            y->size[0] = 1;
            y->size[1] = s->size[1];
            emxEnsureCapacity((emxArray__common *)y, i41, sizeof(creal_T));
            a_re = b_a->data[k + 1].re;
            a_im = b_a->data[k + 1].im;
            loop_ub = s->size[0] * s->size[1];
            for (i41 = 0; i41 < loop_ub; i41++) {
                s_re = s->data[i41].re * y->data[i41].re - s->data[i41].im * y->data[i41]
                       .im;
                s_im = s->data[i41].re * y->data[i41].im + s->data[i41].im * y->data[i41]
                       .re;
                y->data[i41].re = s_re + a_re;
                y->data[i41].im = s_im + a_im;
            }
        }
    }

    j_polyval(b_b_data, b_b_size, s, b_a);
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
    int iy;
    double h_u[7];
    int k;
    static const char cv36[4] = { '1', '2', '1', '1' };

    double combinedResponse_re;
    double combinedResponse_im;
    double d2_re;
    static const char cv37[4] = { '1', '1', '2', '1' };

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
        n_freqz_cg(tmp_data, tmp_size, w, Fs, combinedResponse);
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
        emxEnsureCapacity((emxArray__common *)combinedResponse, ix, sizeof(creal_T));
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
        n_freqz_cg(tmp_data, tmp_size, w, Fs, d2);
        ix = combinedResponse->size[0] * combinedResponse->size[1];
        combinedResponse->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)combinedResponse, ix, sizeof(creal_T));
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
        n_freqz_cg(tmp_data, tmp_size, w, Fs, d2);
        ix = combinedResponse->size[0] * combinedResponse->size[1];
        combinedResponse->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)combinedResponse, ix, sizeof(creal_T));
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
        n_freqz_cg(tmp_data, tmp_size, w, Fs, d3);
        ix = combinedResponse->size[0] * combinedResponse->size[1];
        combinedResponse->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)combinedResponse, ix, sizeof(creal_T));
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
        n_freqz_cg(tmp_data, tmp_size, w, Fs, combinedResponse);

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
        n_freqz_cg(tmp_data, tmp_size, w, Fs, d2);
        ix = combinedResponse->size[0] * combinedResponse->size[1];
        combinedResponse->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)combinedResponse, ix, sizeof(creal_T));
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
        n_freqz_cg(tmp_data, tmp_size, w, Fs, d2);
        ix = combinedResponse->size[0] * combinedResponse->size[1];
        combinedResponse->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)combinedResponse, ix, sizeof(creal_T));
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
        n_freqz_cg(tmp_data, tmp_size, w, Fs, d3);
        ix = combinedResponse->size[0] * combinedResponse->size[1];
        combinedResponse->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)combinedResponse, ix, sizeof(creal_T));
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
    emxEnsureCapacity((emxArray__common *)y, k, sizeof(double));
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
    emxEnsureCapacity((emxArray__common *)z, i30, sizeof(creal_T));
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
 * Arguments    : emxArray_real_T *x
 * Return Type  : void
 */
static void b_sinc(emxArray_real_T *x)
{
    int i73;
    int k;
    i73 = x->size[1];
    for (k = 0; k < i73; k++) {
        if (fabs(x->data[k]) < 1.0020841800044864E-292) {
            x->data[k] = 1.0;
        } else {
            x->data[k] *= 3.1415926535897931;
            x->data[k] = sin(x->data[k]) / x->data[k];
        }
    }
}

/*
 * Arguments    : creal_T *x
 * Return Type  : void
 */
static void b_sqrt(creal_T *x)
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
    emxArray_creal_T *r5;
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
    emxEnsureCapacity((emxArray__common *)c, i5, sizeof(double));
    c->data[0] = 1.0;
    tmp_size[0] = 1;
    tmp_size[1] = 1;
    tmp_data[0].re = -1.0;
    tmp_data[0].im = 0.0;
    b_tmp_size[0] = 1;
    b_tmp_data[0] = 1.0;
    emxInit_creal_T(&a, 2);
    emxInit_real_T1(&b, 1);
    emxInit_creal_T(&r5, 2);
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
    poly(a, r5);
    den_size[0] = 1;
    den_size[1] = r5->size[1];
    loop_ub = r5->size[0] * r5->size[1];
    emxFree_real_T(&c);
    emxFree_real_T(&b);
    emxFree_creal_T(&a);
    for (i5 = 0; i5 < loop_ub; i5++) {
        den_data[i5] = r5->data[i5];
    }

    emxFree_creal_T(&r5);

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
 * Arguments    : const double x_data[]
 *                const int x_size[2]
 *                double y_data[]
 *                int y_size[2]
 * Return Type  : void
 */
static void c_abs(const double x_data[], const int x_size[2], double y_data[],
                  int y_size[2])
{
    int k;
    y_size[0] = 1;
    y_size[1] = (unsigned char)x_size[1];
    for (k = 0; k + 1 <= x_size[1]; k++) {
        y_data[k] = fabs(x_data[k]);
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
    emxArray_real_T *r31;
    int b_abc;
    int loop_ub;
    emxArray_real_T *r32;
    emxArray_creal_T *r33;
    emxArray_real_T *r34;
    emxArray_real_T *r35;
    double abc_re;
    double abc_im;
    emxInit_real_T(&r31, 2);
    b_abc = r31->size[0] * r31->size[1];
    r31->size[0] = 1;
    r31->size[1] = f->size[1];
    emxEnsureCapacity((emxArray__common *)r31, b_abc, sizeof(double));
    loop_ub = f->size[0] * f->size[1];
    for (b_abc = 0; b_abc < loop_ub; b_abc++) {
        r31->data[b_abc] = 6.2831853071795862 * f->data[b_abc];
    }

    emxInit_real_T(&r32, 2);
    b_freqs_cg(b1_data, b1_size, a1, r31, abc);
    b_abc = r32->size[0] * r32->size[1];
    r32->size[0] = 1;
    r32->size[1] = f->size[1];
    emxEnsureCapacity((emxArray__common *)r32, b_abc, sizeof(double));
    loop_ub = f->size[0] * f->size[1];
    emxFree_real_T(&r31);
    for (b_abc = 0; b_abc < loop_ub; b_abc++) {
        r32->data[b_abc] = 6.2831853071795862 * f->data[b_abc];
    }

    emxInit_creal_T(&r33, 2);
    emxInit_real_T(&r34, 2);
    emxInit_real_T(&r35, 2);
    b_freqs_cg(b2_data, b2_size, a2, r32, r33);
    rdivide(f, Fconverter, r35);
    b_sinc(r35);
    b_power(r35, r34);
    b_abc = abc->size[0] * abc->size[1];
    abc->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)abc, b_abc, sizeof(creal_T));
    b_abc = abc->size[0];
    loop_ub = abc->size[1];
    loop_ub *= b_abc;
    emxFree_real_T(&r32);
    emxFree_real_T(&r35);
    for (b_abc = 0; b_abc < loop_ub; b_abc++) {
        abc_re = abc->data[b_abc].re * r33->data[b_abc].re - abc->data[b_abc].im *
                 r33->data[b_abc].im;
        abc_im = abc->data[b_abc].re * r33->data[b_abc].im + abc->data[b_abc].im *
                 r33->data[b_abc].re;
        abc->data[b_abc].re = r34->data[b_abc] * abc_re;
        abc->data[b_abc].im = r34->data[b_abc] * abc_im;
    }

    emxFree_real_T(&r34);
    emxFree_creal_T(&r33);
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
    double bim;
    double digw[2048];
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
        bim = 6.2831853071795862 * options->w[i65] / options->Fs;
        dcv3[i65].re = bim * 0.0;
        dcv3[i65].im = bim;
        digw[i65] = bim;
    }

    b_exp(dcv3);
    b_polyval(b, dcv3, h);
    for (i65 = 0; i65 < 2048; i65++) {
        dcv3[i65].re = 6.0 * (digw[i65] * 0.0);
        dcv3[i65].im = 6.0 * digw[i65];
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
            bim = fabs(dcv3[i65].im);
            if (brm > bim) {
                bim = dcv3[i65].im / dcv3[i65].re;
                d = dcv3[i65].re + bim * dcv3[i65].im;
                h[i65].re = (h[i65].re + bim * h[i65].im) / d;
                h[i65].im = (h[i65].im - bim * h_re) / d;
            } else if (bim == brm) {
                if (dcv3[i65].re > 0.0) {
                    bim = 0.5;
                } else {
                    bim = -0.5;
                }

                if (dcv3[i65].im > 0.0) {
                    d = 0.5;
                } else {
                    d = -0.5;
                }

                h[i65].re = (h[i65].re * bim + h[i65].im * d) / brm;
                h[i65].im = (h[i65].im * bim - h_re * d) / brm;
            } else {
                bim = dcv3[i65].re / dcv3[i65].im;
                d = dcv3[i65].im + bim * dcv3[i65].re;
                h[i65].re = (bim * h[i65].re + h[i65].im) / d;
                h[i65].im = (bim * h[i65].im - h_re) / d;
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
    int exitg1;
    emxArray_creal_T *d2;
    emxArray_creal_T *d3;
    static const char cv48[4] = { '2', '1', '1', '1' };

    double u[15];
    double tmp_data[29];
    int tmp_size[2];
    double b_u[30];
    double c_u[14];
    double d_u[60];
    double e_u[45];
    double f_u[21];
    double g_u[90];
    int stages;
    double h_u[7];
    int k;
    double u_data[1024];
    static const char cv49[4] = { '1', '2', '1', '1' };

    int u_size[2];
    double b_u_data[1023];
    double combinedResponse_re;
    double combinedResponse_im;
    double d2_re;
    static const char cv50[4] = { '1', '1', '2', '1' };

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
                if (ix + 1 < 5) {
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
                    if (ix + 1 < 5) {
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
                        if (ix + 1 < 5) {
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
                            ix = 5;
                        } else {
                            b_bool = false;
                            ix = 0;
                            do {
                                exitg1 = 0;
                                if (ix + 1 < 5) {
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
                                    if (ix + 1 < 5) {
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
                                        if (ix + 1 < 5) {
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
                                            if (ix + 1 < 5) {
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
                                                if (ix + 1 < 5) {
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
                                                    if (ix + 1 < 5) {
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

        l_freqz_cg(u, w, Fs, combinedResponse);
        stages = 1;
        break;

    case 2:
        /*  Hb2 */
        for (ix = 0; ix < 7; ix++) {
            h_u[ix] = 0.0;
        }

        ix = 0;
        stages = 0;
        for (k = 0; k < 7; k++) {
            h_u[stages] = hb2_coeff[ix];
            ix++;
            stages++;
        }

        m_freqz_cg(h_u, w, Fs, combinedResponse);
        stages = 1;
        break;

    case 3:
        /*  Hb3 */
        c_us(hb3_coeff_data, hb3_coeff_size, tmp_data, tmp_size);
        n_freqz_cg(tmp_data, tmp_size, w, Fs, combinedResponse);
        stages = 1;
        break;

    case 4:
        /*  Hb2,Hb1 */
        memset(&b_u[0], 0, 30U * sizeof(double));
        ix = 0;
        stages = 0;
        for (k = 0; k < 15; k++) {
            b_u[stages] = hb1_coeff[ix];
            ix++;
            stages += 2;
        }

        o_freqz_cg(*(double (*)[29])&b_u[0], w, Fs, combinedResponse);
        for (ix = 0; ix < 7; ix++) {
            h_u[ix] = 0.0;
        }

        ix = 0;
        stages = 0;
        for (k = 0; k < 7; k++) {
            h_u[stages] = hb2_coeff[ix];
            ix++;
            stages++;
        }

        m_freqz_cg(h_u, w, Fs, d2);
        ix = combinedResponse->size[0] * combinedResponse->size[1];
        combinedResponse->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)combinedResponse, ix, sizeof(creal_T));
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
        memset(&b_u[0], 0, 30U * sizeof(double));
        ix = 0;
        stages = 0;
        for (k = 0; k < 15; k++) {
            b_u[stages] = hb1_coeff[ix];
            ix++;
            stages += 2;
        }

        o_freqz_cg(*(double (*)[29])&b_u[0], w, Fs, combinedResponse);
        c_us(hb3_coeff_data, hb3_coeff_size, tmp_data, tmp_size);
        n_freqz_cg(tmp_data, tmp_size, w, Fs, d2);
        ix = combinedResponse->size[0] * combinedResponse->size[1];
        combinedResponse->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)combinedResponse, ix, sizeof(creal_T));
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
        memset(&c_u[0], 0, 14U * sizeof(double));
        ix = 0;
        stages = 0;
        for (k = 0; k < 7; k++) {
            c_u[stages] = hb2_coeff[ix];
            ix++;
            stages += 2;
        }

        p_freqz_cg(*(double (*)[13])&c_u[0], w, Fs, combinedResponse);
        c_us(hb3_coeff_data, hb3_coeff_size, tmp_data, tmp_size);
        n_freqz_cg(tmp_data, tmp_size, w, Fs, d2);
        ix = combinedResponse->size[0] * combinedResponse->size[1];
        combinedResponse->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)combinedResponse, ix, sizeof(creal_T));
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

        q_freqz_cg(*(double (*)[57])&d_u[0], w, Fs, combinedResponse);
        memset(&c_u[0], 0, 14U * sizeof(double));
        ix = 0;
        stages = 0;
        for (k = 0; k < 7; k++) {
            c_u[stages] = hb2_coeff[ix];
            ix++;
            stages += 2;
        }

        p_freqz_cg(*(double (*)[13])&c_u[0], w, Fs, d2);
        c_us(hb3_coeff_data, hb3_coeff_size, tmp_data, tmp_size);
        n_freqz_cg(tmp_data, tmp_size, w, Fs, d3);
        ix = combinedResponse->size[0] * combinedResponse->size[1];
        combinedResponse->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)combinedResponse, ix, sizeof(creal_T));
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
        n_freqz_cg(tmp_data, tmp_size, w, Fs, combinedResponse);
        stages = 1;

        /*  RECHECK ALL DEC BY 3     */
        break;

    case 9:
        /*  Dec/Int3,Hb1 */
        memset(&e_u[0], 0, 45U * sizeof(double));
        ix = 0;
        stages = 0;
        for (k = 0; k < 15; k++) {
            e_u[stages] = hb1_coeff[ix];
            ix++;
            stages += 3;
        }

        r_freqz_cg(*(double (*)[43])&e_u[0], w, Fs, combinedResponse);
        c_us(dec_int3_coeff_data, dec_int3_coeff_size, tmp_data, tmp_size);
        n_freqz_cg(tmp_data, tmp_size, w, Fs, d2);
        ix = combinedResponse->size[0] * combinedResponse->size[1];
        combinedResponse->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)combinedResponse, ix, sizeof(creal_T));
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
        memset(&f_u[0], 0, 21U * sizeof(double));
        ix = 0;
        stages = 0;
        for (k = 0; k < 7; k++) {
            f_u[stages] = hb2_coeff[ix];
            ix++;
            stages += 3;
        }

        s_freqz_cg(*(double (*)[19])&f_u[0], w, Fs, combinedResponse);
        c_us(dec_int3_coeff_data, dec_int3_coeff_size, tmp_data, tmp_size);
        n_freqz_cg(tmp_data, tmp_size, w, Fs, d2);
        ix = combinedResponse->size[0] * combinedResponse->size[1];
        combinedResponse->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)combinedResponse, ix, sizeof(creal_T));
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
        /*  Dec/Int3,Hb2,Hb1 {Hm4,Hm2c34,Hm1} */
        memset(&g_u[0], 0, 90U * sizeof(double));
        ix = 0;
        stages = 0;
        for (k = 0; k < 15; k++) {
            g_u[stages] = hb1_coeff[ix];
            ix++;
            stages += 6;
        }

        t_freqz_cg(*(double (*)[85])&g_u[0], w, Fs, combinedResponse);
        memset(&f_u[0], 0, 21U * sizeof(double));
        ix = 0;
        stages = 0;
        for (k = 0; k < 7; k++) {
            f_u[stages] = hb2_coeff[ix];
            ix++;
            stages += 3;
        }

        s_freqz_cg(*(double (*)[19])&f_u[0], w, Fs, d2);
        c_us(dec_int3_coeff_data, dec_int3_coeff_size, tmp_data, tmp_size);
        n_freqz_cg(tmp_data, tmp_size, w, Fs, d3);
        ix = combinedResponse->size[0] * combinedResponse->size[1];
        combinedResponse->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)combinedResponse, ix, sizeof(creal_T));
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
        ix = (int)rt_powd_snf(2.0, stages);
        upsample(extraTaps_data, extraTaps_size, ix, u_data, tmp_size);
        ix = (tmp_size[1] - ix) + 1;
        if (1 > ix) {
            stages = 0;
        } else {
            stages = ix;
        }

        u_size[0] = 1;
        u_size[1] = stages;
        for (ix = 0; ix < stages; ix++) {
            b_u_data[ix] = u_data[ix];
        }

        n_freqz_cg(b_u_data, u_size, w, Fs, d2);
        ix = combinedResponse->size[0] * combinedResponse->size[1];
        combinedResponse->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)combinedResponse, ix, sizeof(creal_T));
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
 * Arguments    : const emxArray_real_T *x
 *                const emxArray_real_T *y
 *                emxArray_real_T *z
 * Return Type  : void
 */
static void c_rdivide(const emxArray_real_T *x, const emxArray_real_T *y,
                      emxArray_real_T *z)
{
    int i51;
    int loop_ub;
    i51 = z->size[0] * z->size[1];
    z->size[0] = 1;
    z->size[1] = x->size[1];
    emxEnsureCapacity((emxArray__common *)z, i51, sizeof(double));
    loop_ub = x->size[0] * x->size[1];
    for (i51 = 0; i51 < loop_ub; i51++) {
        z->data[i51] = x->data[i51] / y->data[i51];
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
    double b_u_data[29];
    int k;
    ix = (signed char)o_size[1];
    for (iy = 0; iy < ix; iy++) {
        b_u_data[iy] = 0.0;
    }

    ix = 0;
    iy = 0;
    for (k = 1; k <= o_size[1]; k++) {
        b_u_data[iy] = o_data[ix];
        ix++;
        iy++;
    }

    ix = (signed char)o_size[1];
    u_size[0] = 1;
    u_size[1] = (signed char)o_size[1];
    for (iy = 0; iy < ix; iy++) {
        u_data[iy] = b_u_data[iy];
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
    emxArray_real_T *r36;
    emxArray_real_T *r37;
    int i57;
    int loop_ub;
    emxArray_real_T *r38;
    emxArray_creal_T *r39;
    double re;
    double im;
    double b_re;
    double b_im;
    emxInit_real_T(&r36, 2);
    emxInit_real_T(&r37, 2);
    rdivide(f, Fconverter, r36);
    b_sinc(r36);
    i57 = r37->size[0] * r37->size[1];
    r37->size[0] = 1;
    r37->size[1] = f->size[1];
    emxEnsureCapacity((emxArray__common *)r37, i57, sizeof(double));
    loop_ub = f->size[0] * f->size[1];
    for (i57 = 0; i57 < loop_ub; i57++) {
        r37->data[i57] = 6.2831853071795862 * f->data[i57];
    }

    emxInit_real_T(&r38, 2);
    b_freqs_cg(b1_data, b1_size, a1, r37, abc);
    i57 = r38->size[0] * r38->size[1];
    r38->size[0] = 1;
    r38->size[1] = f->size[1];
    emxEnsureCapacity((emxArray__common *)r38, i57, sizeof(double));
    loop_ub = f->size[0] * f->size[1];
    emxFree_real_T(&r37);
    for (i57 = 0; i57 < loop_ub; i57++) {
        r38->data[i57] = 6.2831853071795862 * f->data[i57];
    }

    emxInit_creal_T(&r39, 2);
    b_freqs_cg(b2_data, b2_size, a2, r38, r39);
    i57 = abc->size[0] * abc->size[1];
    abc->size[0] = 1;
    abc->size[1] = r36->size[1];
    emxEnsureCapacity((emxArray__common *)abc, i57, sizeof(creal_T));
    loop_ub = r36->size[0] * r36->size[1];
    emxFree_real_T(&r38);
    for (i57 = 0; i57 < loop_ub; i57++) {
        re = r36->data[i57] * abc->data[i57].re;
        im = r36->data[i57] * abc->data[i57].im;
        b_re = r39->data[i57].re;
        b_im = r39->data[i57].im;
        abc->data[i57].re = re * b_re - im * b_im;
        abc->data[i57].im = re * b_im + im * b_re;
    }

    emxFree_real_T(&r36);
    emxFree_creal_T(&r39);
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
    int i66;
    double b_b_data[29];
    creal_T dcv4[2048];
    double bim;
    double digw[2048];
    double h_re;
    double brm;
    double d;

    /* -------------------------------------------------------------------------- */
    b_idx_0 = b_size[1];
    loop_ub = b_size[1];
    for (i66 = 0; i66 < loop_ub; i66++) {
        b_b_data[i66] = b_data[i66];
    }

    b_size[0] = 1;
    for (i66 = 0; i66 < b_idx_0; i66++) {
        b_data[i66] = b_b_data[i66];
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
    for (i66 = 0; i66 < 2048; i66++) {
        w[i66] = options->w[i66];
        bim = 6.2831853071795862 * options->w[i66] / options->Fs;
        dcv4[i66].re = bim * 0.0;
        dcv4[i66].im = bim;
        digw[i66] = bim;
    }

    b_exp(dcv4);
    c_polyval(b_data, b_size, dcv4, h);
    for (i66 = 0; i66 < 2048; i66++) {
        dcv4[i66].re = ((double)b_size[1] - 1.0) * (digw[i66] * 0.0);
        dcv4[i66].im = ((double)b_size[1] - 1.0) * digw[i66];
    }

    b_exp(dcv4);
    for (i66 = 0; i66 < 2048; i66++) {
        h_re = h[i66].re;
        if (dcv4[i66].im == 0.0) {
            if (h[i66].im == 0.0) {
                h[i66].re /= dcv4[i66].re;
                h[i66].im = 0.0;
            } else if (h[i66].re == 0.0) {
                h[i66].re = 0.0;
                h[i66].im /= dcv4[i66].re;
            } else {
                h[i66].re /= dcv4[i66].re;
                h[i66].im /= dcv4[i66].re;
            }
        } else if (dcv4[i66].re == 0.0) {
            if (h[i66].re == 0.0) {
                h[i66].re = h[i66].im / dcv4[i66].im;
                h[i66].im = 0.0;
            } else if (h[i66].im == 0.0) {
                h[i66].re = 0.0;
                h[i66].im = -(h_re / dcv4[i66].im);
            } else {
                h[i66].re = h[i66].im / dcv4[i66].im;
                h[i66].im = -(h_re / dcv4[i66].im);
            }
        } else {
            brm = fabs(dcv4[i66].re);
            bim = fabs(dcv4[i66].im);
            if (brm > bim) {
                bim = dcv4[i66].im / dcv4[i66].re;
                d = dcv4[i66].re + bim * dcv4[i66].im;
                h[i66].re = (h[i66].re + bim * h[i66].im) / d;
                h[i66].im = (h[i66].im - bim * h_re) / d;
            } else if (bim == brm) {
                if (dcv4[i66].re > 0.0) {
                    bim = 0.5;
                } else {
                    bim = -0.5;
                }

                if (dcv4[i66].im > 0.0) {
                    d = 0.5;
                } else {
                    d = -0.5;
                }

                h[i66].re = (h[i66].re * bim + h[i66].im * d) / brm;
                h[i66].im = (h[i66].im * bim - h_re * d) / brm;
            } else {
                bim = dcv4[i66].re / dcv4[i66].im;
                d = dcv4[i66].im + bim * dcv4[i66].re;
                h[i66].re = (bim * h[i66].re + h[i66].im) / d;
                h[i66].im = (bim * h[i66].im - h_re) / d;
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
    int i14;
    static const char cv18[8] = { 'o', 'n', 'e', 's', 'i', 'd', 'e', 'd' };

    int b_b_size[2];
    static const char cv19[7] = { 'o', 'm', 'i', 't', 't', 'e', 'd' };

    int loop_ub;
    double b_b_data[29];
    double b_w[2048];

    /*  Cast to enforce precision rules */
    options.Fs = Fs;
    memcpy(&options.nfft[0], &w[0], sizeof(double) << 11);
    memcpy(&options.w[0], &w[0], sizeof(double) << 11);

    /*  Remaining are default or for advanced use */
    options.fvflag = 1.0;
    for (i14 = 0; i14 < 8; i14++) {
        options.range[i14] = cv18[i14];
    }

    options.centerdc = 0.0;
    for (i14 = 0; i14 < 7; i14++) {
        options.configlevel[i14] = cv19[i14];
    }

    b_b_size[0] = 1;
    b_b_size[1] = b_size[1];
    loop_ub = b_size[0] * b_size[1];
    for (i14 = 0; i14 < loop_ub; i14++) {
        b_b_data[i14] = b_data[i14];
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
    int i15;
    int k;
    double x_im;
    for (i15 = 0; i15 < 2048; i15++) {
        y[i15].re = p[0];
        y[i15].im = 0.0;
    }

    for (k = 0; k < 28; k++) {
        for (i15 = 0; i15 < 2048; i15++) {
            x_im = x[i15].re * y[i15].im + x[i15].im * y[i15].re;
            y[i15].re = (x[i15].re * y[i15].re - x[i15].im * y[i15].im) + p[k + 1];
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
 * Codegen workaround for fixed fi call requirements
 * Arguments    : const double tap_store_data[]
 *                const int tap_store_size[2]
 *                double i
 *                double M
 *                emxArray_real_T *taps
 * Return Type  : void
 */
static void determineBestFractionLength(const double tap_store_data[], const int
                                        tap_store_size[2], double i, double M, emxArray_real_T *taps)
{
    int ixstart;
    int ix;
    emxArray_real_T *r;
    double org_data[128];
    double mtmp;
    double v;
    int r_size[2];
    short i56;
    double r_data[128];
    double tmp_data[128];
    int tmp_size[2];
    emxArray_real_T b_tmp_data;
    double e[16];
    int b_r_size[2];
    emxArray_real_T c_tmp_data;
    int c_r_size[2];
    emxArray_real_T d_tmp_data;
    int d_r_size[2];
    emxArray_real_T e_tmp_data;
    int e_r_size[2];
    emxArray_real_T f_tmp_data;
    int f_r_size[2];
    emxArray_real_T g_tmp_data;
    int g_r_size[2];
    emxArray_real_T h_tmp_data;
    int h_r_size[2];
    emxArray_real_T i_tmp_data;
    int i_r_size[2];
    emxArray_real_T j_tmp_data;
    int j_r_size[2];
    emxArray_real_T k_tmp_data;
    int k_r_size[2];
    emxArray_real_T l_tmp_data;
    int l_r_size[2];
    emxArray_real_T m_tmp_data;
    int m_r_size[2];
    emxArray_real_T n_tmp_data;
    int n_r_size[2];
    emxArray_real_T o_tmp_data;
    int o_r_size[2];
    emxArray_real_T p_tmp_data;
    int p_r_size[2];
    emxArray_real_T q_tmp_data;
    int itmp;
    boolean_T exitg1;
    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    for (ix = 0; ix < ixstart; ix++) {
        org_data[ix] = tap_store_data[((int)i + tap_store_size[0] * ix) - 1];
    }

    emxInit_real_T(&r, 2);
    ix = r->size[0] * r->size[1];
    r->size[0] = 16;
    r->size[1] = (int)M;
    emxEnsureCapacity((emxArray__common *)r, ix, sizeof(double));
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
        mtmp = tap_store_data[((int)i + tap_store_size[0] * ix) - 1] * 2.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i56 = (short)mtmp;
            } else {
                i56 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i56 = MAX_int16_T;
        } else {
            i56 = 0;
        }

        r->data[r->size[0] * ix] = (double)i56 * 0.5;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    r_size[0] = 1;
    r_size[1] = ixstart;
    for (ix = 0; ix < ixstart; ix++) {
        r_data[ix] = r->data[r->size[0] * ix] - org_data[ix];
    }

    c_abs(r_data, r_size, tmp_data, tmp_size);
    b_tmp_data.data = (double *)&tmp_data;
    b_tmp_data.size = (int *)&tmp_size;
    b_tmp_data.allocatedSize = 128;
    b_tmp_data.numDimensions = 2;
    b_tmp_data.canFreeData = false;
    e[0] = b_sum(&b_tmp_data);
    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        mtmp = tap_store_data[((int)i + tap_store_size[0] * ix) - 1] * 4.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i56 = (short)mtmp;
            } else {
                i56 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i56 = MAX_int16_T;
        } else {
            i56 = 0;
        }

        r->data[1 + r->size[0] * ix] = (double)i56 * 0.25;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    b_r_size[0] = 1;
    b_r_size[1] = ixstart;
    for (ix = 0; ix < ixstart; ix++) {
        r_data[ix] = r->data[1 + r->size[0] * ix] - org_data[ix];
    }

    c_abs(r_data, b_r_size, tmp_data, tmp_size);
    c_tmp_data.data = (double *)&tmp_data;
    c_tmp_data.size = (int *)&tmp_size;
    c_tmp_data.allocatedSize = 128;
    c_tmp_data.numDimensions = 2;
    c_tmp_data.canFreeData = false;
    e[1] = b_sum(&c_tmp_data);
    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        mtmp = tap_store_data[((int)i + tap_store_size[0] * ix) - 1] * 8.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i56 = (short)mtmp;
            } else {
                i56 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i56 = MAX_int16_T;
        } else {
            i56 = 0;
        }

        r->data[2 + r->size[0] * ix] = (double)i56 * 0.125;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    c_r_size[0] = 1;
    c_r_size[1] = ixstart;
    for (ix = 0; ix < ixstart; ix++) {
        r_data[ix] = r->data[2 + r->size[0] * ix] - org_data[ix];
    }

    c_abs(r_data, c_r_size, tmp_data, tmp_size);
    d_tmp_data.data = (double *)&tmp_data;
    d_tmp_data.size = (int *)&tmp_size;
    d_tmp_data.allocatedSize = 128;
    d_tmp_data.numDimensions = 2;
    d_tmp_data.canFreeData = false;
    e[2] = b_sum(&d_tmp_data);
    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        mtmp = tap_store_data[((int)i + tap_store_size[0] * ix) - 1] * 16.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i56 = (short)mtmp;
            } else {
                i56 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i56 = MAX_int16_T;
        } else {
            i56 = 0;
        }

        r->data[3 + r->size[0] * ix] = (double)i56 * 0.0625;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    d_r_size[0] = 1;
    d_r_size[1] = ixstart;
    for (ix = 0; ix < ixstart; ix++) {
        r_data[ix] = r->data[3 + r->size[0] * ix] - org_data[ix];
    }

    c_abs(r_data, d_r_size, tmp_data, tmp_size);
    e_tmp_data.data = (double *)&tmp_data;
    e_tmp_data.size = (int *)&tmp_size;
    e_tmp_data.allocatedSize = 128;
    e_tmp_data.numDimensions = 2;
    e_tmp_data.canFreeData = false;
    e[3] = b_sum(&e_tmp_data);
    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        mtmp = tap_store_data[((int)i + tap_store_size[0] * ix) - 1] * 32.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i56 = (short)mtmp;
            } else {
                i56 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i56 = MAX_int16_T;
        } else {
            i56 = 0;
        }

        r->data[4 + r->size[0] * ix] = (double)i56 * 0.03125;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    e_r_size[0] = 1;
    e_r_size[1] = ixstart;
    for (ix = 0; ix < ixstart; ix++) {
        r_data[ix] = r->data[4 + r->size[0] * ix] - org_data[ix];
    }

    c_abs(r_data, e_r_size, tmp_data, tmp_size);
    f_tmp_data.data = (double *)&tmp_data;
    f_tmp_data.size = (int *)&tmp_size;
    f_tmp_data.allocatedSize = 128;
    f_tmp_data.numDimensions = 2;
    f_tmp_data.canFreeData = false;
    e[4] = b_sum(&f_tmp_data);
    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        mtmp = tap_store_data[((int)i + tap_store_size[0] * ix) - 1] * 64.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i56 = (short)mtmp;
            } else {
                i56 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i56 = MAX_int16_T;
        } else {
            i56 = 0;
        }

        r->data[5 + r->size[0] * ix] = (double)i56 * 0.015625;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    f_r_size[0] = 1;
    f_r_size[1] = ixstart;
    for (ix = 0; ix < ixstart; ix++) {
        r_data[ix] = r->data[5 + r->size[0] * ix] - org_data[ix];
    }

    c_abs(r_data, f_r_size, tmp_data, tmp_size);
    g_tmp_data.data = (double *)&tmp_data;
    g_tmp_data.size = (int *)&tmp_size;
    g_tmp_data.allocatedSize = 128;
    g_tmp_data.numDimensions = 2;
    g_tmp_data.canFreeData = false;
    e[5] = b_sum(&g_tmp_data);
    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        mtmp = tap_store_data[((int)i + tap_store_size[0] * ix) - 1] * 128.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i56 = (short)mtmp;
            } else {
                i56 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i56 = MAX_int16_T;
        } else {
            i56 = 0;
        }

        r->data[6 + r->size[0] * ix] = (double)i56 * 0.0078125;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    g_r_size[0] = 1;
    g_r_size[1] = ixstart;
    for (ix = 0; ix < ixstart; ix++) {
        r_data[ix] = r->data[6 + r->size[0] * ix] - org_data[ix];
    }

    c_abs(r_data, g_r_size, tmp_data, tmp_size);
    h_tmp_data.data = (double *)&tmp_data;
    h_tmp_data.size = (int *)&tmp_size;
    h_tmp_data.allocatedSize = 128;
    h_tmp_data.numDimensions = 2;
    h_tmp_data.canFreeData = false;
    e[6] = b_sum(&h_tmp_data);
    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        mtmp = tap_store_data[((int)i + tap_store_size[0] * ix) - 1] * 256.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i56 = (short)mtmp;
            } else {
                i56 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i56 = MAX_int16_T;
        } else {
            i56 = 0;
        }

        r->data[7 + r->size[0] * ix] = (double)i56 * 0.00390625;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    h_r_size[0] = 1;
    h_r_size[1] = ixstart;
    for (ix = 0; ix < ixstart; ix++) {
        r_data[ix] = r->data[7 + r->size[0] * ix] - org_data[ix];
    }

    c_abs(r_data, h_r_size, tmp_data, tmp_size);
    i_tmp_data.data = (double *)&tmp_data;
    i_tmp_data.size = (int *)&tmp_size;
    i_tmp_data.allocatedSize = 128;
    i_tmp_data.numDimensions = 2;
    i_tmp_data.canFreeData = false;
    e[7] = b_sum(&i_tmp_data);
    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        mtmp = tap_store_data[((int)i + tap_store_size[0] * ix) - 1] * 512.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i56 = (short)mtmp;
            } else {
                i56 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i56 = MAX_int16_T;
        } else {
            i56 = 0;
        }

        r->data[8 + r->size[0] * ix] = (double)i56 * 0.001953125;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    i_r_size[0] = 1;
    i_r_size[1] = ixstart;
    for (ix = 0; ix < ixstart; ix++) {
        r_data[ix] = r->data[8 + r->size[0] * ix] - org_data[ix];
    }

    c_abs(r_data, i_r_size, tmp_data, tmp_size);
    j_tmp_data.data = (double *)&tmp_data;
    j_tmp_data.size = (int *)&tmp_size;
    j_tmp_data.allocatedSize = 128;
    j_tmp_data.numDimensions = 2;
    j_tmp_data.canFreeData = false;
    e[8] = b_sum(&j_tmp_data);
    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        mtmp = tap_store_data[((int)i + tap_store_size[0] * ix) - 1] * 1024.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i56 = (short)mtmp;
            } else {
                i56 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i56 = MAX_int16_T;
        } else {
            i56 = 0;
        }

        r->data[9 + r->size[0] * ix] = (double)i56 * 0.0009765625;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    j_r_size[0] = 1;
    j_r_size[1] = ixstart;
    for (ix = 0; ix < ixstart; ix++) {
        r_data[ix] = r->data[9 + r->size[0] * ix] - org_data[ix];
    }

    c_abs(r_data, j_r_size, tmp_data, tmp_size);
    k_tmp_data.data = (double *)&tmp_data;
    k_tmp_data.size = (int *)&tmp_size;
    k_tmp_data.allocatedSize = 128;
    k_tmp_data.numDimensions = 2;
    k_tmp_data.canFreeData = false;
    e[9] = b_sum(&k_tmp_data);
    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        mtmp = tap_store_data[((int)i + tap_store_size[0] * ix) - 1] * 2048.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i56 = (short)mtmp;
            } else {
                i56 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i56 = MAX_int16_T;
        } else {
            i56 = 0;
        }

        r->data[10 + r->size[0] * ix] = (double)i56 * 0.00048828125;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    k_r_size[0] = 1;
    k_r_size[1] = ixstart;
    for (ix = 0; ix < ixstart; ix++) {
        r_data[ix] = r->data[10 + r->size[0] * ix] - org_data[ix];
    }

    c_abs(r_data, k_r_size, tmp_data, tmp_size);
    l_tmp_data.data = (double *)&tmp_data;
    l_tmp_data.size = (int *)&tmp_size;
    l_tmp_data.allocatedSize = 128;
    l_tmp_data.numDimensions = 2;
    l_tmp_data.canFreeData = false;
    e[10] = b_sum(&l_tmp_data);
    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        mtmp = tap_store_data[((int)i + tap_store_size[0] * ix) - 1] * 4096.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i56 = (short)mtmp;
            } else {
                i56 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i56 = MAX_int16_T;
        } else {
            i56 = 0;
        }

        r->data[11 + r->size[0] * ix] = (double)i56 * 0.000244140625;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    l_r_size[0] = 1;
    l_r_size[1] = ixstart;
    for (ix = 0; ix < ixstart; ix++) {
        r_data[ix] = r->data[11 + r->size[0] * ix] - org_data[ix];
    }

    c_abs(r_data, l_r_size, tmp_data, tmp_size);
    m_tmp_data.data = (double *)&tmp_data;
    m_tmp_data.size = (int *)&tmp_size;
    m_tmp_data.allocatedSize = 128;
    m_tmp_data.numDimensions = 2;
    m_tmp_data.canFreeData = false;
    e[11] = b_sum(&m_tmp_data);
    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        mtmp = tap_store_data[((int)i + tap_store_size[0] * ix) - 1] * 8192.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i56 = (short)mtmp;
            } else {
                i56 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i56 = MAX_int16_T;
        } else {
            i56 = 0;
        }

        r->data[12 + r->size[0] * ix] = (double)i56 * 0.0001220703125;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    m_r_size[0] = 1;
    m_r_size[1] = ixstart;
    for (ix = 0; ix < ixstart; ix++) {
        r_data[ix] = r->data[12 + r->size[0] * ix] - org_data[ix];
    }

    c_abs(r_data, m_r_size, tmp_data, tmp_size);
    n_tmp_data.data = (double *)&tmp_data;
    n_tmp_data.size = (int *)&tmp_size;
    n_tmp_data.allocatedSize = 128;
    n_tmp_data.numDimensions = 2;
    n_tmp_data.canFreeData = false;
    e[12] = b_sum(&n_tmp_data);
    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        mtmp = tap_store_data[((int)i + tap_store_size[0] * ix) - 1] * 16384.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i56 = (short)mtmp;
            } else {
                i56 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i56 = MAX_int16_T;
        } else {
            i56 = 0;
        }

        r->data[13 + r->size[0] * ix] = (double)i56 * 6.103515625E-5;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    n_r_size[0] = 1;
    n_r_size[1] = ixstart;
    for (ix = 0; ix < ixstart; ix++) {
        r_data[ix] = r->data[13 + r->size[0] * ix] - org_data[ix];
    }

    c_abs(r_data, n_r_size, tmp_data, tmp_size);
    o_tmp_data.data = (double *)&tmp_data;
    o_tmp_data.size = (int *)&tmp_size;
    o_tmp_data.allocatedSize = 128;
    o_tmp_data.numDimensions = 2;
    o_tmp_data.canFreeData = false;
    e[13] = b_sum(&o_tmp_data);
    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        mtmp = tap_store_data[((int)i + tap_store_size[0] * ix) - 1] * 32768.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i56 = (short)mtmp;
            } else {
                i56 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i56 = MAX_int16_T;
        } else {
            i56 = 0;
        }

        r->data[14 + r->size[0] * ix] = (double)i56 * 3.0517578125E-5;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    o_r_size[0] = 1;
    o_r_size[1] = ixstart;
    for (ix = 0; ix < ixstart; ix++) {
        r_data[ix] = r->data[14 + r->size[0] * ix] - org_data[ix];
    }

    c_abs(r_data, o_r_size, tmp_data, tmp_size);
    p_tmp_data.data = (double *)&tmp_data;
    p_tmp_data.size = (int *)&tmp_size;
    p_tmp_data.allocatedSize = 128;
    p_tmp_data.numDimensions = 2;
    p_tmp_data.canFreeData = false;
    e[14] = b_sum(&p_tmp_data);
    if (1.0 > M) {
        ixstart = -1;
    } else {
        ixstart = (int)M - 1;
    }

    for (ix = 0; ix <= ixstart; ix++) {
        mtmp = tap_store_data[((int)i + tap_store_size[0] * ix) - 1] * 65536.0;
        v = fabs(mtmp);
        if (v < 4.503599627370496E+15) {
            if (v >= 0.5) {
                mtmp = floor(mtmp + 0.5);
            } else {
                mtmp *= 0.0;
            }
        }

        if (mtmp < 32768.0) {
            if (mtmp >= -32768.0) {
                i56 = (short)mtmp;
            } else {
                i56 = MIN_int16_T;
            }
        } else if (mtmp >= 32768.0) {
            i56 = MAX_int16_T;
        } else {
            i56 = 0;
        }

        r->data[15 + r->size[0] * ix] = (double)i56 * 1.52587890625E-5;
    }

    if (1.0 > M) {
        ixstart = 0;
    } else {
        ixstart = (int)M;
    }

    p_r_size[0] = 1;
    p_r_size[1] = ixstart;
    for (ix = 0; ix < ixstart; ix++) {
        r_data[ix] = r->data[15 + r->size[0] * ix] - org_data[ix];
    }

    c_abs(r_data, p_r_size, tmp_data, tmp_size);
    q_tmp_data.data = (double *)&tmp_data;
    q_tmp_data.size = (int *)&tmp_size;
    q_tmp_data.allocatedSize = 128;
    q_tmp_data.numDimensions = 2;
    q_tmp_data.canFreeData = false;
    e[15] = b_sum(&q_tmp_data);
    ixstart = 1;
    mtmp = e[0];
    itmp = 0;
    if (rtIsNaN(e[0])) {
        ix = 2;
        exitg1 = false;
        while ((!exitg1) && (ix < 17)) {
            ixstart = ix;
            if (!rtIsNaN(e[ix - 1])) {
                mtmp = e[ix - 1];
                itmp = ix - 1;
                exitg1 = true;
            } else {
                ix++;
            }
        }
    }

    if (ixstart < 16) {
        while (ixstart + 1 < 17) {
            if (e[ixstart] < mtmp) {
                mtmp = e[ixstart];
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
    emxEnsureCapacity((emxArray__common *)taps, ix, sizeof(double));
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
    int i67;
    creal_T dcv5[2048];
    double bim;
    double digw[2048];
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
        bim = 6.2831853071795862 * options->w[i67] / options->Fs;
        dcv5[i67].re = bim * 0.0;
        dcv5[i67].im = bim;
        digw[i67] = bim;
    }

    b_exp(dcv5);
    d_polyval(b, dcv5, h);
    for (i67 = 0; i67 < 2048; i67++) {
        dcv5[i67].re = 28.0 * (digw[i67] * 0.0);
        dcv5[i67].im = 28.0 * digw[i67];
    }

    b_exp(dcv5);
    for (i67 = 0; i67 < 2048; i67++) {
        h_re = h[i67].re;
        if (dcv5[i67].im == 0.0) {
            if (h[i67].im == 0.0) {
                h[i67].re /= dcv5[i67].re;
                h[i67].im = 0.0;
            } else if (h[i67].re == 0.0) {
                h[i67].re = 0.0;
                h[i67].im /= dcv5[i67].re;
            } else {
                h[i67].re /= dcv5[i67].re;
                h[i67].im /= dcv5[i67].re;
            }
        } else if (dcv5[i67].re == 0.0) {
            if (h[i67].re == 0.0) {
                h[i67].re = h[i67].im / dcv5[i67].im;
                h[i67].im = 0.0;
            } else if (h[i67].im == 0.0) {
                h[i67].re = 0.0;
                h[i67].im = -(h_re / dcv5[i67].im);
            } else {
                h[i67].re = h[i67].im / dcv5[i67].im;
                h[i67].im = -(h_re / dcv5[i67].im);
            }
        } else {
            brm = fabs(dcv5[i67].re);
            bim = fabs(dcv5[i67].im);
            if (brm > bim) {
                bim = dcv5[i67].im / dcv5[i67].re;
                d = dcv5[i67].re + bim * dcv5[i67].im;
                h[i67].re = (h[i67].re + bim * h[i67].im) / d;
                h[i67].im = (h[i67].im - bim * h_re) / d;
            } else if (bim == brm) {
                if (dcv5[i67].re > 0.0) {
                    bim = 0.5;
                } else {
                    bim = -0.5;
                }

                if (dcv5[i67].im > 0.0) {
                    d = 0.5;
                } else {
                    d = -0.5;
                }

                h[i67].re = (h[i67].re * bim + h[i67].im * d) / brm;
                h[i67].im = (h[i67].im * bim - h_re * d) / brm;
            } else {
                bim = dcv5[i67].re / dcv5[i67].im;
                d = dcv5[i67].im + bim * dcv5[i67].re;
                h[i67].re = (bim * h[i67].re + h[i67].im) / d;
                h[i67].im = (bim * h[i67].im - h_re) / d;
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
static void e_polyval(const double p[13], const creal_T x[2048],
                      creal_T y[2048])
{
    int i17;
    int k;
    double x_im;
    for (i17 = 0; i17 < 2048; i17++) {
        y[i17].re = p[0];
        y[i17].im = 0.0;
    }

    for (k = 0; k < 12; k++) {
        for (i17 = 0; i17 < 2048; i17++) {
            x_im = x[i17].re * y[i17].im + x[i17].im * y[i17].re;
            y[i17].re = (x[i17].re * y[i17].re - x[i17].im * y[i17].im) + p[k + 1];
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
        emxEnsureCapacity((emxArray__common *)V, info, sizeof(creal_T));
        i = A->size[0];
        for (info = 0; info < i; info++) {
            V->data[info].re = 0.0;
            V->data[info].im = 0.0;
        }
    } else if (anyNonFinite(A)) {
        if ((A->size[0] == 1) && (A->size[1] == 1)) {
            info = V->size[0];
            V->size[0] = 1;
            emxEnsureCapacity((emxArray__common *)V, info, sizeof(creal_T));
            V->data[0].re = rtNaN;
            V->data[0].im = 0.0;
        } else {
            info = V->size[0];
            V->size[0] = A->size[0];
            emxEnsureCapacity((emxArray__common *)V, info, sizeof(creal_T));
            i = A->size[0];
            for (info = 0; info < i; info++) {
                V->data[info].re = rtNaN;
                V->data[info].im = 0.0;
            }
        }
    } else if ((A->size[0] == 1) && (A->size[1] == 1)) {
        info = V->size[0];
        V->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)V, info, sizeof(creal_T));
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
                emxEnsureCapacity((emxArray__common *)T, info, sizeof(creal_T));
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
                emxEnsureCapacity((emxArray__common *)T, info, sizeof(creal_T));
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
            emxEnsureCapacity((emxArray__common *)V, info, sizeof(creal_T));
            for (info = 0; info + 1 <= T->size[0]; info++) {
                V->data[info] = T->data[info + T->size[0] * info];
            }

            emxFree_creal_T(&T);
        } else {
            emxInit_creal_T1(&T, 1);
            xzgeev(A, &info, V, T);
            info = V->size[0];
            emxEnsureCapacity((emxArray__common *)V, info, sizeof(creal_T));
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
    int ldh;
    int knt;
    int i;
    double SMLNUM;
    double aa;
    boolean_T exitg1;
    double tst;
    double ba;
    int L;
    creal_T u2;
    boolean_T goto140;
    int its;
    boolean_T exitg2;
    int k;
    boolean_T exitg3;
    creal_T b_u2;
    double htmp1;
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
    int c_k;
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
                aa = h->data[i + h->size[0] * (i - 1)].re;
                tst = h->data[i + h->size[0] * (i - 1)].im;
                ba = fabs(h->data[i + h->size[0] * (i - 1)].re) + fabs(h->data[i +
                        h->size[0] * (i - 1)].im);
                if (tst == 0.0) {
                    u2.re = aa / ba;
                    u2.im = 0.0;
                } else if (aa == 0.0) {
                    u2.re = 0.0;
                    u2.im = tst / ba;
                } else {
                    u2.re = aa / ba;
                    u2.im = tst / ba;
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
                aa = h->data[i + h->size[0] * (i - 1)].im;
                h->data[i + h->size[0] * (i - 1)].re = rt_hypotd_snf(tst, aa);
                h->data[i + h->size[0] * (i - 1)].im = 0.0;
                b_xscal(n - i, u2, h, (i + i * ldh) + 1, ldh);
                b_u2.re = u2.re;
                b_u2.im = -u2.im;
                knt = i + 2;
                if (n < knt) {
                    knt = n;
                }

                xscal(knt, b_u2, h, 1 + i * ldh);
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
                        b_u2 = h->data[(i + h->size[0] * i) - 1];
                        b_sqrt(&b_u2);
                        u2 = h->data[i + h->size[0] * (i - 1)];
                        b_sqrt(&u2);
                        u_re = b_u2.re * u2.re - b_u2.im * u2.im;
                        u_im = b_u2.re * u2.im + b_u2.im * u2.re;
                        s = fabs(u_re) + fabs(u_im);
                        if (s != 0.0) {
                            tst = h->data[(i + h->size[0] * (i - 1)) - 1].re - h->data[i +
                                    h->size[0] * i].re;
                            aa = h->data[(i + h->size[0] * (i - 1)) - 1].im - h->data[i +
                                    h->size[0] * i].im;
                            x_re = 0.5 * tst;
                            x_im = 0.5 * aa;
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

                            htmp1 = ba;
                            ba = ba * ba - ab * ab;
                            ab = htmp1 * ab + ab * htmp1;
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

                            b_u2.re = ba + (u2.re * u2.re - u2.im * u2.im);
                            b_u2.im = ab + (u2.re * u2.im + u2.im * u2.re);
                            b_sqrt(&b_u2);
                            u2.re = s * b_u2.re;
                            u2.im = s * b_u2.im;
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

                                if (ba * u2.re + ab * u2.im < 0.0) {
                                    u2.re = -u2.re;
                                    u2.im = -u2.im;
                                }
                            }

                            ba = x_re + u2.re;
                            htmp1 = x_im + u2.im;
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

                        u2 = v[0];
                        x_re = v[1].re;
                        x_im = v[1].im;
                        ba = 0.0;
                        ab = 0.0;
                        tst = rt_hypotd_snf(v[1].re, v[1].im);
                        if ((tst != 0.0) || (v[0].im != 0.0)) {
                            htmp1 = xdlapy3(v[0].re, v[0].im, tst);
                            if (v[0].re >= 0.0) {
                                htmp1 = -htmp1;
                            }

                            if (fabs(htmp1) < 1.0020841800044864E-292) {
                                knt = 0;
                                do {
                                    knt++;
                                    x_re *= 9.9792015476736E+291;
                                    x_im *= 9.9792015476736E+291;
                                    htmp1 *= 9.9792015476736E+291;
                                    u2.re *= 9.9792015476736E+291;
                                    u2.im *= 9.9792015476736E+291;
                                } while (!(fabs(htmp1) >= 1.0020841800044864E-292));

                                htmp1 = xdlapy3(u2.re, u2.im, rt_hypotd_snf(x_re, x_im));
                                if (u2.re >= 0.0) {
                                    htmp1 = -htmp1;
                                }

                                aa = htmp1 - u2.re;
                                if (0.0 - u2.im == 0.0) {
                                    ba = aa / htmp1;
                                    ab = 0.0;
                                } else if (aa == 0.0) {
                                    ba = 0.0;
                                    ab = (0.0 - u2.im) / htmp1;
                                } else {
                                    ba = aa / htmp1;
                                    ab = (0.0 - u2.im) / htmp1;
                                }

                                b_u2.re = u2.re - htmp1;
                                b_u2.im = u2.im;
                                u2 = recip(b_u2);
                                tst = x_re;
                                x_re = u2.re * x_re - u2.im * x_im;
                                x_im = u2.re * x_im + u2.im * tst;
                                for (c_k = 1; c_k <= knt; c_k++) {
                                    htmp1 *= 1.0020841800044864E-292;
                                }

                                u2.re = htmp1;
                                u2.im = 0.0;
                            } else {
                                aa = htmp1 - v[0].re;
                                if (0.0 - v[0].im == 0.0) {
                                    ba = aa / htmp1;
                                    ab = 0.0;
                                } else if (aa == 0.0) {
                                    ba = 0.0;
                                    ab = (0.0 - v[0].im) / htmp1;
                                } else {
                                    ba = aa / htmp1;
                                    ab = (0.0 - v[0].im) / htmp1;
                                }

                                b_u2.re = v[0].re - htmp1;
                                b_u2.im = v[0].im;
                                b_u2 = recip(b_u2);
                                x_re = b_u2.re * v[1].re - b_u2.im * v[1].im;
                                x_im = b_u2.re * v[1].im + b_u2.im * v[1].re;
                                u2.re = htmp1;
                                u2.im = 0.0;
                            }
                        }

                        v[0] = u2;
                        v[1].re = x_re;
                        v[1].im = x_im;
                        if (b_k > m) {
                            h->data[(b_k + h->size[0] * (b_k - 2)) - 1] = u2;
                            h->data[b_k + h->size[0] * (b_k - 2)].re = 0.0;
                            h->data[b_k + h->size[0] * (b_k - 2)].im = 0.0;
                        }

                        htmp1 = ba * x_re - ab * x_im;
                        for (knt = b_k - 1; knt + 1 <= n; knt++) {
                            tst = ba * h->data[(b_k + h->size[0] * knt) - 1].re - -ab *
                                  h->data[(b_k + h->size[0] * knt) - 1].im;
                            aa = ba * h->data[(b_k + h->size[0] * knt) - 1].im + -ab * h->
                                 data[(b_k + h->size[0] * knt) - 1].re;
                            u2.re = tst + htmp1 * h->data[b_k + h->size[0] * knt].re;
                            u2.im = aa + htmp1 * h->data[b_k + h->size[0] * knt].im;
                            h->data[(b_k + h->size[0] * knt) - 1].re -= u2.re;
                            h->data[(b_k + h->size[0] * knt) - 1].im -= u2.im;
                            h->data[b_k + h->size[0] * knt].re -= u2.re * x_re - u2.im * x_im;
                            h->data[b_k + h->size[0] * knt].im -= u2.re * x_im + u2.im * x_re;
                        }

                        if (b_k + 2 < i + 1) {
                            c_k = b_k;
                        } else {
                            c_k = i - 1;
                        }

                        for (knt = 0; knt + 1 <= c_k + 2; knt++) {
                            tst = ba * h->data[knt + h->size[0] * (b_k - 1)].re - ab * h->
                                  data[knt + h->size[0] * (b_k - 1)].im;
                            aa = ba * h->data[knt + h->size[0] * (b_k - 1)].im + ab * h->
                                 data[knt + h->size[0] * (b_k - 1)].re;
                            u2.re = tst + htmp1 * h->data[knt + h->size[0] * b_k].re;
                            u2.im = aa + htmp1 * h->data[knt + h->size[0] * b_k].im;
                            h->data[knt + h->size[0] * (b_k - 1)].re -= u2.re;
                            h->data[knt + h->size[0] * (b_k - 1)].im -= u2.im;
                            h->data[knt + h->size[0] * b_k].re -= u2.re * x_re - u2.im * -x_im;
                            h->data[knt + h->size[0] * b_k].im -= u2.re * -x_im + u2.im * x_re;
                        }

                        if ((b_k == m) && (m > k + 1)) {
                            u2.re = 1.0 - ba;
                            u2.im = 0.0 - ab;
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

                                    b_u2.re = u2.re;
                                    b_u2.im = -u2.im;
                                    xscal(knt - 1, b_u2, h, 1 + (knt - 1) * ldh);
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
                            b_u2.re = u2.re;
                            b_u2.im = -u2.im;
                            b_xscal((n - i) - 1, b_u2, h, (i + (i + 1) * ldh) + 1, ldh);
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
    int i68;
    creal_T dcv6[2048];
    double bim;
    double digw[2048];
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
        bim = 6.2831853071795862 * options->w[i68] / options->Fs;
        dcv6[i68].re = bim * 0.0;
        dcv6[i68].im = bim;
        digw[i68] = bim;
    }

    b_exp(dcv6);
    e_polyval(b, dcv6, h);
    for (i68 = 0; i68 < 2048; i68++) {
        dcv6[i68].re = 12.0 * (digw[i68] * 0.0);
        dcv6[i68].im = 12.0 * digw[i68];
    }

    b_exp(dcv6);
    for (i68 = 0; i68 < 2048; i68++) {
        h_re = h[i68].re;
        if (dcv6[i68].im == 0.0) {
            if (h[i68].im == 0.0) {
                h[i68].re /= dcv6[i68].re;
                h[i68].im = 0.0;
            } else if (h[i68].re == 0.0) {
                h[i68].re = 0.0;
                h[i68].im /= dcv6[i68].re;
            } else {
                h[i68].re /= dcv6[i68].re;
                h[i68].im /= dcv6[i68].re;
            }
        } else if (dcv6[i68].re == 0.0) {
            if (h[i68].re == 0.0) {
                h[i68].re = h[i68].im / dcv6[i68].im;
                h[i68].im = 0.0;
            } else if (h[i68].im == 0.0) {
                h[i68].re = 0.0;
                h[i68].im = -(h_re / dcv6[i68].im);
            } else {
                h[i68].re = h[i68].im / dcv6[i68].im;
                h[i68].im = -(h_re / dcv6[i68].im);
            }
        } else {
            brm = fabs(dcv6[i68].re);
            bim = fabs(dcv6[i68].im);
            if (brm > bim) {
                bim = dcv6[i68].im / dcv6[i68].re;
                d = dcv6[i68].re + bim * dcv6[i68].im;
                h[i68].re = (h[i68].re + bim * h[i68].im) / d;
                h[i68].im = (h[i68].im - bim * h_re) / d;
            } else if (bim == brm) {
                if (dcv6[i68].re > 0.0) {
                    bim = 0.5;
                } else {
                    bim = -0.5;
                }

                if (dcv6[i68].im > 0.0) {
                    d = 0.5;
                } else {
                    d = -0.5;
                }

                h[i68].re = (h[i68].re * bim + h[i68].im * d) / brm;
                h[i68].im = (h[i68].im * bim - h_re * d) / brm;
            } else {
                bim = dcv6[i68].re / dcv6[i68].im;
                d = dcv6[i68].im + bim * dcv6[i68].re;
                h[i68].re = (bim * h[i68].re + h[i68].im) / d;
                h[i68].im = (bim * h[i68].im - h_re) / d;
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
static void f_polyval(const double p[57], const creal_T x[2048],
                      creal_T y[2048])
{
    int i19;
    int k;
    double x_im;
    for (i19 = 0; i19 < 2048; i19++) {
        y[i19].re = p[0];
        y[i19].im = 0.0;
    }

    for (k = 0; k < 56; k++) {
        for (i19 = 0; i19 < 2048; i19++) {
            x_im = x[i19].re * y[i19].im + x[i19].im * y[i19].re;
            y[i19].re = (x[i19].re * y[i19].re - x[i19].im * y[i19].im) + p[k + 1];
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
    emxArray_real_T *r24;
    int i42;
    emxArray_real_T *b_h;
    double b_ff[4];
    double err;
    boolean_T valid;
    int h_idx_0;
    int i43;
    emxArray_real_T *c_h;
    int i44;
    int loop_ub;
    emxArray_real_T *d_h;
    emxInit_real_T(&grid, 2);
    emxInit_real_T(&des, 2);
    emxInit_real_T(&wt, 2);
    emxInit_real_T(&r24, 2);

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
    rdivide(grid, 2.0, r24);
    emxFree_real_T(&grid);
    for (i42 = 0; i42 < 4; i42++) {
        b_ff[i42] = ff[i42] / 2.0;
    }

    emxInit_real_T(&b_h, 2);
    remezm(order + 1.0, b_ff, r24, des, wt, b_h, &err, &valid);
    h_idx_0 = b_h->size[0] * b_h->size[1];
    i42 = h->size[0] * h->size[1];
    h->size[0] = 1;
    h->size[1] = h_idx_0;
    emxEnsureCapacity((emxArray__common *)h, i42, sizeof(double));
    emxFree_real_T(&r24);
    emxFree_real_T(&wt);
    emxFree_real_T(&des);
    for (i42 = 0; i42 < h_idx_0; i42++) {
        h->data[h->size[0] * i42] = b_h->data[i42];
    }

    emxFree_real_T(&b_h);

    /*  make it a row */
    err = (double)h->size[1] - rt_remd_snf(order + 1.0, 2.0);
    if (1.0 > err) {
        i42 = 1;
        h_idx_0 = 1;
        i43 = 0;
    } else {
        i42 = (int)err;
        h_idx_0 = -1;
        i43 = 1;
    }

    emxInit_real_T(&c_h, 2);
    i44 = c_h->size[0] * c_h->size[1];
    c_h->size[0] = 1;
    c_h->size[1] = (h->size[1] + div_s32_floor(i43 - i42, h_idx_0)) + 1;
    emxEnsureCapacity((emxArray__common *)c_h, i44, sizeof(double));
    loop_ub = h->size[1];
    for (i44 = 0; i44 < loop_ub; i44++) {
        c_h->data[c_h->size[0] * i44] = h->data[h->size[0] * i44];
    }

    loop_ub = div_s32_floor(i43 - i42, h_idx_0);
    for (i43 = 0; i43 <= loop_ub; i43++) {
        c_h->data[c_h->size[0] * (i43 + h->size[1])] = h->data[(i42 + h_idx_0 * i43)
                - 1];
    }

    i42 = h->size[0] * h->size[1];
    h->size[0] = 1;
    h->size[1] = c_h->size[1];
    emxEnsureCapacity((emxArray__common *)h, i42, sizeof(double));
    loop_ub = c_h->size[1];
    for (i42 = 0; i42 < loop_ub; i42++) {
        h->data[h->size[0] * i42] = c_h->data[c_h->size[0] * i42];
    }

    emxFree_real_T(&c_h);
    if (1 > h->size[1]) {
        i42 = 1;
        h_idx_0 = 1;
        i43 = 0;
    } else {
        i42 = h->size[1];
        h_idx_0 = -1;
        i43 = 1;
    }

    emxInit_real_T(&d_h, 2);
    i44 = d_h->size[0] * d_h->size[1];
    d_h->size[0] = 1;
    d_h->size[1] = div_s32_floor(i43 - i42, h_idx_0) + 1;
    emxEnsureCapacity((emxArray__common *)d_h, i44, sizeof(double));
    loop_ub = div_s32_floor(i43 - i42, h_idx_0);
    for (i43 = 0; i43 <= loop_ub; i43++) {
        d_h->data[d_h->size[0] * i43] = h->data[(i42 + h_idx_0 * i43) - 1];
    }

    i42 = h->size[0] * h->size[1];
    h->size[0] = 1;
    h->size[1] = d_h->size[1];
    emxEnsureCapacity((emxArray__common *)h, i42, sizeof(double));
    loop_ub = d_h->size[1];
    for (i42 = 0; i42 < loop_ub; i42++) {
        h->data[h->size[0] * i42] = d_h->data[d_h->size[0] * i42];
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
    emxArray_int32_T *r25;
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
    emxInit_int32_T(&r25, 2);
    while (l + 1.0 <= 4.0) {
        a = grid[(int)j - 1] + delf;
        ngrid = ff[(int)(l + 1.0) - 1] + delf;
        if (rtIsNaN(a) || rtIsNaN(delf) || rtIsNaN(ngrid)) {
            k = newgrid->size[0] * newgrid->size[1];
            newgrid->size[0] = 1;
            newgrid->size[1] = 1;
            emxEnsureCapacity((emxArray__common *)newgrid, k, sizeof(double));
            newgrid->data[0] = rtNaN;
        } else if ((delf == 0.0) || ((a < ngrid) && (delf < 0.0)) || ((ngrid < a) &&
                   (delf > 0.0))) {
            k = newgrid->size[0] * newgrid->size[1];
            newgrid->size[0] = 1;
            newgrid->size[1] = 0;
            emxEnsureCapacity((emxArray__common *)newgrid, k, sizeof(double));
        } else if ((rtIsInf(a) || rtIsInf(ngrid)) && (rtIsInf(delf) || (a == ngrid))) {
            k = newgrid->size[0] * newgrid->size[1];
            newgrid->size[0] = 1;
            newgrid->size[1] = 1;
            emxEnsureCapacity((emxArray__common *)newgrid, k, sizeof(double));
            newgrid->data[0] = rtNaN;
        } else if (rtIsInf(delf)) {
            k = newgrid->size[0] * newgrid->size[1];
            newgrid->size[0] = 1;
            newgrid->size[1] = 1;
            emxEnsureCapacity((emxArray__common *)newgrid, k, sizeof(double));
            newgrid->data[0] = a;
        } else if ((floor(a) == a) && (floor(delf) == delf)) {
            k = newgrid->size[0] * newgrid->size[1];
            newgrid->size[0] = 1;
            newgrid->size[1] = (int)floor((ngrid - a) / delf) + 1;
            emxEnsureCapacity((emxArray__common *)newgrid, k, sizeof(double));
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
            emxEnsureCapacity((emxArray__common *)newgrid, k, sizeof(double));
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
                emxEnsureCapacity((emxArray__common *)newgrid, k, sizeof(double));
                newgrid->data[0] = rtNaN;
            } else if ((delf1 == 0.0) || ((a < ngrid) && (delf1 < 0.0)) || ((ngrid < a)
                       && (delf1 > 0.0))) {
                k = newgrid->size[0] * newgrid->size[1];
                newgrid->size[0] = 1;
                newgrid->size[1] = 0;
                emxEnsureCapacity((emxArray__common *)newgrid, k, sizeof(double));
            } else if ((rtIsInf(a) || rtIsInf(ngrid)) && (rtIsInf(delf1) || (a ==
                       ngrid))) {
                k = newgrid->size[0] * newgrid->size[1];
                newgrid->size[0] = 1;
                newgrid->size[1] = 1;
                emxEnsureCapacity((emxArray__common *)newgrid, k, sizeof(double));
                newgrid->data[0] = rtNaN;
            } else if (rtIsInf(delf1)) {
                k = newgrid->size[0] * newgrid->size[1];
                newgrid->size[0] = 1;
                newgrid->size[1] = 1;
                emxEnsureCapacity((emxArray__common *)newgrid, k, sizeof(double));
                newgrid->data[0] = a;
            } else if ((floor(a) == a) && (floor(delf1) == delf1)) {
                k = newgrid->size[0] * newgrid->size[1];
                newgrid->size[0] = 1;
                newgrid->size[1] = (int)floor((ngrid - a) / delf1) + 1;
                emxEnsureCapacity((emxArray__common *)newgrid, k, sizeof(double));
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
                emxEnsureCapacity((emxArray__common *)newgrid, k, sizeof(double));
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
        nm1d2 = r25->size[0] * r25->size[1];
        r25->size[0] = 1;
        r25->size[1] = (int)((double)k - 1.0) + 1;
        emxEnsureCapacity((emxArray__common *)r25, nm1d2, sizeof(int));
        nm1d2 = (int)((double)k - 1.0);
        for (k = 0; k <= nm1d2; k++) {
            r25->data[r25->size[0] * k] = (int)(gridSize + (1.0 + (double)k));
        }

        nm1d2 = newgrid->size[0] * newgrid->size[1];
        for (k = 0; k < nm1d2; k++) {
            grid[r25->data[k] - 1] = newgrid->data[k];
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

    emxFree_int32_T(&r25);
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
    emxEnsureCapacity((emxArray__common *)gridactual, k, sizeof(double));
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
    int i27;
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
    for (i27 = 0; i27 < 2048; i27++) {
        s[i27].re = w[i27] * 0.0;
        s[i27].im = w[i27];
    }

    b0 = (b_a->size[1] == 0);
    if (!b0) {
        for (i27 = 0; i27 < 2048; i27++) {
            y[i27] = b_a->data[0];
        }

        for (k = 0; k <= b_a->size[1] - 2; k++) {
            bim = b_a->data[k + 1].re;
            d = b_a->data[k + 1].im;
            for (i27 = 0; i27 < 2048; i27++) {
                brm = s[i27].re * y[i27].im + s[i27].im * y[i27].re;
                y[i27].re = (s[i27].re * y[i27].re - s[i27].im * y[i27].im) + bim;
                y[i27].im = brm + d;
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
            brm = fabs(y[i27].re);
            bim = fabs(y[i27].im);
            if (brm > bim) {
                bim = y[i27].im / y[i27].re;
                d = y[i27].re + bim * y[i27].im;
                h[i27].re = (h[i27].re + bim * h[i27].im) / d;
                h[i27].im = (h[i27].im - bim * h_re) / d;
            } else if (bim == brm) {
                if (y[i27].re > 0.0) {
                    bim = 0.5;
                } else {
                    bim = -0.5;
                }

                if (y[i27].im > 0.0) {
                    d = 0.5;
                } else {
                    d = -0.5;
                }

                h[i27].re = (h[i27].re * bim + h[i27].im * d) / brm;
                h[i27].im = (h[i27].im * bim - h_re * d) / brm;
            } else {
                bim = y[i27].re / y[i27].im;
                d = y[i27].im + bim * y[i27].re;
                h[i27].re = (bim * h[i27].re + h[i27].im) / d;
                h[i27].im = (bim * h[i27].im - h_re) / d;
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
    int i69;
    creal_T dcv7[2048];
    double bim;
    double digw[2048];
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
        bim = 6.2831853071795862 * options->w[i69] / options->Fs;
        dcv7[i69].re = bim * 0.0;
        dcv7[i69].im = bim;
        digw[i69] = bim;
    }

    b_exp(dcv7);
    f_polyval(b, dcv7, h);
    for (i69 = 0; i69 < 2048; i69++) {
        dcv7[i69].re = 56.0 * (digw[i69] * 0.0);
        dcv7[i69].im = 56.0 * digw[i69];
    }

    b_exp(dcv7);
    for (i69 = 0; i69 < 2048; i69++) {
        h_re = h[i69].re;
        if (dcv7[i69].im == 0.0) {
            if (h[i69].im == 0.0) {
                h[i69].re /= dcv7[i69].re;
                h[i69].im = 0.0;
            } else if (h[i69].re == 0.0) {
                h[i69].re = 0.0;
                h[i69].im /= dcv7[i69].re;
            } else {
                h[i69].re /= dcv7[i69].re;
                h[i69].im /= dcv7[i69].re;
            }
        } else if (dcv7[i69].re == 0.0) {
            if (h[i69].re == 0.0) {
                h[i69].re = h[i69].im / dcv7[i69].im;
                h[i69].im = 0.0;
            } else if (h[i69].im == 0.0) {
                h[i69].re = 0.0;
                h[i69].im = -(h_re / dcv7[i69].im);
            } else {
                h[i69].re = h[i69].im / dcv7[i69].im;
                h[i69].im = -(h_re / dcv7[i69].im);
            }
        } else {
            brm = fabs(dcv7[i69].re);
            bim = fabs(dcv7[i69].im);
            if (brm > bim) {
                bim = dcv7[i69].im / dcv7[i69].re;
                d = dcv7[i69].re + bim * dcv7[i69].im;
                h[i69].re = (h[i69].re + bim * h[i69].im) / d;
                h[i69].im = (h[i69].im - bim * h_re) / d;
            } else if (bim == brm) {
                if (dcv7[i69].re > 0.0) {
                    bim = 0.5;
                } else {
                    bim = -0.5;
                }

                if (dcv7[i69].im > 0.0) {
                    d = 0.5;
                } else {
                    d = -0.5;
                }

                h[i69].re = (h[i69].re * bim + h[i69].im * d) / brm;
                h[i69].im = (h[i69].im * bim - h_re * d) / brm;
            } else {
                bim = dcv7[i69].re / dcv7[i69].im;
                d = dcv7[i69].im + bim * dcv7[i69].re;
                h[i69].re = (bim * h[i69].re + h[i69].im) / d;
                h[i69].im = (bim * h[i69].im - h_re) / d;
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
static void g_polyval(const double p[43], const creal_T x[2048],
                      creal_T y[2048])
{
    int i21;
    int k;
    double x_im;
    for (i21 = 0; i21 < 2048; i21++) {
        y[i21].re = p[0];
        y[i21].im = 0.0;
    }

    for (k = 0; k < 42; k++) {
        for (i21 = 0; i21 < 2048; i21++) {
            x_im = x[i21].re * y[i21].im + x[i21].im * y[i21].re;
            y[i21].re = (x[i21].re * y[i21].re - x[i21].im * y[i21].im) + p[k + 1];
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
    int i70;
    creal_T dcv8[2048];
    double bim;
    double digw[2048];
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
        bim = 6.2831853071795862 * options->w[i70] / options->Fs;
        dcv8[i70].re = bim * 0.0;
        dcv8[i70].im = bim;
        digw[i70] = bim;
    }

    b_exp(dcv8);
    g_polyval(b, dcv8, h);
    for (i70 = 0; i70 < 2048; i70++) {
        dcv8[i70].re = 42.0 * (digw[i70] * 0.0);
        dcv8[i70].im = 42.0 * digw[i70];
    }

    b_exp(dcv8);
    for (i70 = 0; i70 < 2048; i70++) {
        h_re = h[i70].re;
        if (dcv8[i70].im == 0.0) {
            if (h[i70].im == 0.0) {
                h[i70].re /= dcv8[i70].re;
                h[i70].im = 0.0;
            } else if (h[i70].re == 0.0) {
                h[i70].re = 0.0;
                h[i70].im /= dcv8[i70].re;
            } else {
                h[i70].re /= dcv8[i70].re;
                h[i70].im /= dcv8[i70].re;
            }
        } else if (dcv8[i70].re == 0.0) {
            if (h[i70].re == 0.0) {
                h[i70].re = h[i70].im / dcv8[i70].im;
                h[i70].im = 0.0;
            } else if (h[i70].im == 0.0) {
                h[i70].re = 0.0;
                h[i70].im = -(h_re / dcv8[i70].im);
            } else {
                h[i70].re = h[i70].im / dcv8[i70].im;
                h[i70].im = -(h_re / dcv8[i70].im);
            }
        } else {
            brm = fabs(dcv8[i70].re);
            bim = fabs(dcv8[i70].im);
            if (brm > bim) {
                bim = dcv8[i70].im / dcv8[i70].re;
                d = dcv8[i70].re + bim * dcv8[i70].im;
                h[i70].re = (h[i70].re + bim * h[i70].im) / d;
                h[i70].im = (h[i70].im - bim * h_re) / d;
            } else if (bim == brm) {
                if (dcv8[i70].re > 0.0) {
                    bim = 0.5;
                } else {
                    bim = -0.5;
                }

                if (dcv8[i70].im > 0.0) {
                    d = 0.5;
                } else {
                    d = -0.5;
                }

                h[i70].re = (h[i70].re * bim + h[i70].im * d) / brm;
                h[i70].im = (h[i70].im * bim - h_re * d) / brm;
            } else {
                bim = dcv8[i70].re / dcv8[i70].im;
                d = dcv8[i70].im + bim * dcv8[i70].re;
                h[i70].re = (bim * h[i70].re + h[i70].im) / d;
                h[i70].im = (bim * h[i70].im - h_re) / d;
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
static void h_polyval(const double p[19], const creal_T x[2048],
                      creal_T y[2048])
{
    int i23;
    int k;
    double x_im;
    for (i23 = 0; i23 < 2048; i23++) {
        y[i23].re = p[0];
        y[i23].im = 0.0;
    }

    for (k = 0; k < 18; k++) {
        for (i23 = 0; i23 < 2048; i23++) {
            x_im = x[i23].re * y[i23].im + x[i23].im * y[i23].re;
            y[i23].re = (x[i23].re * y[i23].re - x[i23].im * y[i23].im) + p[k + 1];
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
    int i71;
    creal_T dcv9[2048];
    double bim;
    double digw[2048];
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
        bim = 6.2831853071795862 * options->w[i71] / options->Fs;
        dcv9[i71].re = bim * 0.0;
        dcv9[i71].im = bim;
        digw[i71] = bim;
    }

    b_exp(dcv9);
    h_polyval(b, dcv9, h);
    for (i71 = 0; i71 < 2048; i71++) {
        dcv9[i71].re = 18.0 * (digw[i71] * 0.0);
        dcv9[i71].im = 18.0 * digw[i71];
    }

    b_exp(dcv9);
    for (i71 = 0; i71 < 2048; i71++) {
        h_re = h[i71].re;
        if (dcv9[i71].im == 0.0) {
            if (h[i71].im == 0.0) {
                h[i71].re /= dcv9[i71].re;
                h[i71].im = 0.0;
            } else if (h[i71].re == 0.0) {
                h[i71].re = 0.0;
                h[i71].im /= dcv9[i71].re;
            } else {
                h[i71].re /= dcv9[i71].re;
                h[i71].im /= dcv9[i71].re;
            }
        } else if (dcv9[i71].re == 0.0) {
            if (h[i71].re == 0.0) {
                h[i71].re = h[i71].im / dcv9[i71].im;
                h[i71].im = 0.0;
            } else if (h[i71].im == 0.0) {
                h[i71].re = 0.0;
                h[i71].im = -(h_re / dcv9[i71].im);
            } else {
                h[i71].re = h[i71].im / dcv9[i71].im;
                h[i71].im = -(h_re / dcv9[i71].im);
            }
        } else {
            brm = fabs(dcv9[i71].re);
            bim = fabs(dcv9[i71].im);
            if (brm > bim) {
                bim = dcv9[i71].im / dcv9[i71].re;
                d = dcv9[i71].re + bim * dcv9[i71].im;
                h[i71].re = (h[i71].re + bim * h[i71].im) / d;
                h[i71].im = (h[i71].im - bim * h_re) / d;
            } else if (bim == brm) {
                if (dcv9[i71].re > 0.0) {
                    bim = 0.5;
                } else {
                    bim = -0.5;
                }

                if (dcv9[i71].im > 0.0) {
                    d = 0.5;
                } else {
                    d = -0.5;
                }

                h[i71].re = (h[i71].re * bim + h[i71].im * d) / brm;
                h[i71].im = (h[i71].im * bim - h_re * d) / brm;
            } else {
                bim = dcv9[i71].re / dcv9[i71].im;
                d = dcv9[i71].im + bim * dcv9[i71].re;
                h[i71].re = (bim * h[i71].re + h[i71].im) / d;
                h[i71].im = (bim * h[i71].im - h_re) / d;
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
static void i_polyval(const double p[85], const creal_T x[2048],
                      creal_T y[2048])
{
    int i25;
    int k;
    double x_im;
    for (i25 = 0; i25 < 2048; i25++) {
        y[i25].re = p[0];
        y[i25].im = 0.0;
    }

    for (k = 0; k < 84; k++) {
        for (i25 = 0; i25 < 2048; i25++) {
            x_im = x[i25].re * y[i25].im + x[i25].im * y[i25].re;
            y[i25].re = (x[i25].re * y[i25].re - x[i25].im * y[i25].im) + p[k + 1];
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
    emxEnsureCapacity((emxArray__common *)y, low_ip1, sizeof(double));
    nd2 = varargin_2->size[0] * varargin_2->size[1];
    for (low_ip1 = 0; low_ip1 < nd2; low_ip1++) {
        y->data[low_ip1] = varargin_2->data[low_ip1];
    }

    emxInit_real_T(&x, 2);
    low_ip1 = x->size[0] * x->size[1];
    x->size[0] = 1;
    x->size[1] = varargin_1->size[1];
    emxEnsureCapacity((emxArray__common *)x, low_ip1, sizeof(double));
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
    emxEnsureCapacity((emxArray__common *)Vq, low_ip1, sizeof(double));
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
    int i72;
    creal_T dcv10[2048];
    double bim;
    double digw[2048];
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
    for (i72 = 0; i72 < 2048; i72++) {
        w[i72] = options->w[i72];
        bim = 6.2831853071795862 * options->w[i72] / options->Fs;
        dcv10[i72].re = bim * 0.0;
        dcv10[i72].im = bim;
        digw[i72] = bim;
    }

    b_exp(dcv10);
    i_polyval(b, dcv10, h);
    for (i72 = 0; i72 < 2048; i72++) {
        dcv10[i72].re = 84.0 * (digw[i72] * 0.0);
        dcv10[i72].im = 84.0 * digw[i72];
    }

    b_exp(dcv10);
    for (i72 = 0; i72 < 2048; i72++) {
        h_re = h[i72].re;
        if (dcv10[i72].im == 0.0) {
            if (h[i72].im == 0.0) {
                h[i72].re /= dcv10[i72].re;
                h[i72].im = 0.0;
            } else if (h[i72].re == 0.0) {
                h[i72].re = 0.0;
                h[i72].im /= dcv10[i72].re;
            } else {
                h[i72].re /= dcv10[i72].re;
                h[i72].im /= dcv10[i72].re;
            }
        } else if (dcv10[i72].re == 0.0) {
            if (h[i72].re == 0.0) {
                h[i72].re = h[i72].im / dcv10[i72].im;
                h[i72].im = 0.0;
            } else if (h[i72].im == 0.0) {
                h[i72].re = 0.0;
                h[i72].im = -(h_re / dcv10[i72].im);
            } else {
                h[i72].re = h[i72].im / dcv10[i72].im;
                h[i72].im = -(h_re / dcv10[i72].im);
            }
        } else {
            brm = fabs(dcv10[i72].re);
            bim = fabs(dcv10[i72].im);
            if (brm > bim) {
                bim = dcv10[i72].im / dcv10[i72].re;
                d = dcv10[i72].re + bim * dcv10[i72].im;
                h[i72].re = (h[i72].re + bim * h[i72].im) / d;
                h[i72].im = (h[i72].im - bim * h_re) / d;
            } else if (bim == brm) {
                if (dcv10[i72].re > 0.0) {
                    bim = 0.5;
                } else {
                    bim = -0.5;
                }

                if (dcv10[i72].im > 0.0) {
                    d = 0.5;
                } else {
                    d = -0.5;
                }

                h[i72].re = (h[i72].re * bim + h[i72].im * d) / brm;
                h[i72].im = (h[i72].im * bim - h_re * d) / brm;
            } else {
                bim = dcv10[i72].re / dcv10[i72].im;
                d = dcv10[i72].im + bim * dcv10[i72].re;
                h[i72].re = (bim * h[i72].re + h[i72].im) / d;
                h[i72].im = (bim * h[i72].im - h_re) / d;
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
 * Arguments    : const double p_data[]
 *                const int p_size[2]
 *                const emxArray_creal_T *x
 *                emxArray_creal_T *y
 * Return Type  : void
 */
static void j_polyval(const double p_data[], const int p_size[2], const
                      emxArray_creal_T *x, emxArray_creal_T *y)
{
    int i34;
    boolean_T b4;
    int loop_ub;
    int k;
    double x_re;
    double x_im;
    i34 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = x->size[1];
    emxEnsureCapacity((emxArray__common *)y, i34, sizeof(creal_T));
    if ((y->size[1] == 0) || (p_size[1] == 0)) {
        b4 = true;
    } else {
        b4 = false;
    }

    if (!b4) {
        i34 = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)y, i34, sizeof(creal_T));
        loop_ub = y->size[1];
        for (i34 = 0; i34 < loop_ub; i34++) {
            y->data[y->size[0] * i34].re = p_data[0];
            y->data[y->size[0] * i34].im = 0.0;
        }

        for (k = 0; k <= p_size[1] - 2; k++) {
            i34 = y->size[0] * y->size[1];
            y->size[0] = 1;
            y->size[1] = x->size[1];
            emxEnsureCapacity((emxArray__common *)y, i34, sizeof(creal_T));
            loop_ub = x->size[0] * x->size[1];
            for (i34 = 0; i34 < loop_ub; i34++) {
                x_re = x->data[i34].re * y->data[i34].re - x->data[i34].im * y->data[i34]
                       .im;
                x_im = x->data[i34].re * y->data[i34].im + x->data[i34].im * y->data[i34]
                       .re;
                y->data[i34].re = x_re + p_data[k + 1];
                y->data[i34].im = x_im;
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
    emxArray_real_T *r8;
    int i29;
    int loop_ub;
    emxArray_real_T *digw;
    emxArray_creal_T *s;
    emxArray_creal_T *y;
    boolean_T b1;
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
    i29 = r8->size[0] * r8->size[1];
    r8->size[0] = 1;
    r8->size[1] = w->size[1];
    emxEnsureCapacity((emxArray__common *)r8, i29, sizeof(double));
    loop_ub = w->size[0] * w->size[1];
    for (i29 = 0; i29 < loop_ub; i29++) {
        r8->data[i29] = 6.2831853071795862 * w->data[i29];
    }

    emxInit_real_T(&digw, 2);
    emxInit_creal_T(&s, 2);
    rdivide(r8, Fs, digw);

    /*  Convert from Hz to rad/sample for computational purposes */
    i29 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = digw->size[1];
    emxEnsureCapacity((emxArray__common *)s, i29, sizeof(creal_T));
    loop_ub = digw->size[0] * digw->size[1];
    emxFree_real_T(&r8);
    for (i29 = 0; i29 < loop_ub; i29++) {
        s->data[i29].re = digw->data[i29] * 0.0;
        s->data[i29].im = digw->data[i29];
    }

    emxInit_creal_T(&y, 2);
    c_exp(s);

    /*  Digital frequency must be used for this calculation */
    i29 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = s->size[1];
    emxEnsureCapacity((emxArray__common *)y, i29, sizeof(creal_T));
    b1 = (y->size[1] == 0);
    if (!b1) {
        i29 = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)y, i29, sizeof(creal_T));
        loop_ub = y->size[1];
        for (i29 = 0; i29 < loop_ub; i29++) {
            y->data[y->size[0] * i29].re = 1.0;
            y->data[y->size[0] * i29].im = 0.0;
        }
    }

    i29 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = digw->size[1];
    emxEnsureCapacity((emxArray__common *)s, i29, sizeof(creal_T));
    loop_ub = digw->size[0] * digw->size[1];
    for (i29 = 0; i29 < loop_ub; i29++) {
        re = digw->data[i29] * 0.0;
        im = digw->data[i29];
        s->data[i29].re = 0.0 * re;
        s->data[i29].im = 0.0 * im;
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
static void l_freqz_cg(const double b[15], const emxArray_real_T *w, double Fs,
                       emxArray_creal_T *hh)
{
    emxArray_real_T *r9;
    int i31;
    int loop_ub;
    emxArray_real_T *digw;
    emxArray_creal_T *s;
    emxArray_creal_T *y;
    boolean_T b2;
    int k;
    double re;
    double im;
    emxInit_real_T(&r9, 2);

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
    i31 = r9->size[0] * r9->size[1];
    r9->size[0] = 1;
    r9->size[1] = w->size[1];
    emxEnsureCapacity((emxArray__common *)r9, i31, sizeof(double));
    loop_ub = w->size[0] * w->size[1];
    for (i31 = 0; i31 < loop_ub; i31++) {
        r9->data[i31] = 6.2831853071795862 * w->data[i31];
    }

    emxInit_real_T(&digw, 2);
    emxInit_creal_T(&s, 2);
    rdivide(r9, Fs, digw);

    /*  Convert from Hz to rad/sample for computational purposes */
    i31 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = digw->size[1];
    emxEnsureCapacity((emxArray__common *)s, i31, sizeof(creal_T));
    loop_ub = digw->size[0] * digw->size[1];
    emxFree_real_T(&r9);
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
    emxEnsureCapacity((emxArray__common *)y, i31, sizeof(creal_T));
    b2 = (y->size[1] == 0);
    if (!b2) {
        i31 = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)y, i31, sizeof(creal_T));
        loop_ub = y->size[1];
        for (i31 = 0; i31 < loop_ub; i31++) {
            y->data[y->size[0] * i31].re = b[0];
            y->data[y->size[0] * i31].im = 0.0;
        }

        for (k = 0; k < 14; k++) {
            i31 = y->size[0] * y->size[1];
            y->size[0] = 1;
            y->size[1] = s->size[1];
            emxEnsureCapacity((emxArray__common *)y, i31, sizeof(creal_T));
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
    emxEnsureCapacity((emxArray__common *)s, i31, sizeof(creal_T));
    loop_ub = digw->size[0] * digw->size[1];
    for (i31 = 0; i31 < loop_ub; i31++) {
        re = digw->data[i31] * 0.0;
        im = digw->data[i31];
        s->data[i31].re = 14.0 * re;
        s->data[i31].im = 14.0 * im;
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
    int i59;
    int loop_ub;

    /*  Transform lowpass to lowpass */
    i59 = at->size[0] * at->size[1];
    at->size[0] = a->size[0];
    at->size[1] = a->size[1];
    emxEnsureCapacity((emxArray__common *)at, i59, sizeof(creal_T));
    loop_ub = a->size[0] * a->size[1];
    for (i59 = 0; i59 < loop_ub; i59++) {
        at->data[i59].re = wo * a->data[i59].re;
        at->data[i59].im = wo * a->data[i59].im;
    }

    i59 = bt->size[0];
    bt->size[0] = b->size[0];
    emxEnsureCapacity((emxArray__common *)bt, i59, sizeof(double));
    loop_ub = b->size[0];
    for (i59 = 0; i59 < loop_ub; i59++) {
        bt->data[i59] = wo * b->data[i59];
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
    emxArray_real_T *r10;
    int i32;
    int loop_ub;
    emxArray_real_T *digw;
    emxArray_creal_T *s;
    emxArray_creal_T *y;
    boolean_T b3;
    int k;
    double re;
    double im;
    emxInit_real_T(&r10, 2);

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
    i32 = r10->size[0] * r10->size[1];
    r10->size[0] = 1;
    r10->size[1] = w->size[1];
    emxEnsureCapacity((emxArray__common *)r10, i32, sizeof(double));
    loop_ub = w->size[0] * w->size[1];
    for (i32 = 0; i32 < loop_ub; i32++) {
        r10->data[i32] = 6.2831853071795862 * w->data[i32];
    }

    emxInit_real_T(&digw, 2);
    emxInit_creal_T(&s, 2);
    rdivide(r10, Fs, digw);

    /*  Convert from Hz to rad/sample for computational purposes */
    i32 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = digw->size[1];
    emxEnsureCapacity((emxArray__common *)s, i32, sizeof(creal_T));
    loop_ub = digw->size[0] * digw->size[1];
    emxFree_real_T(&r10);
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
    emxEnsureCapacity((emxArray__common *)y, i32, sizeof(creal_T));
    b3 = (y->size[1] == 0);
    if (!b3) {
        i32 = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)y, i32, sizeof(creal_T));
        loop_ub = y->size[1];
        for (i32 = 0; i32 < loop_ub; i32++) {
            y->data[y->size[0] * i32].re = b[0];
            y->data[y->size[0] * i32].im = 0.0;
        }

        for (k = 0; k < 6; k++) {
            i32 = y->size[0] * y->size[1];
            y->size[0] = 1;
            y->size[1] = s->size[1];
            emxEnsureCapacity((emxArray__common *)y, i32, sizeof(creal_T));
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
    emxEnsureCapacity((emxArray__common *)s, i32, sizeof(creal_T));
    loop_ub = digw->size[0] * digw->size[1];
    for (i32 = 0; i32 < loop_ub; i32++) {
        re = digw->data[i32] * 0.0;
        im = digw->data[i32];
        s->data[i32].re = 6.0 * re;
        s->data[i32].im = 6.0 * im;
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
static void n_freqz_cg(const double b_data[], const int b_size[2], const
                       emxArray_real_T *w, double Fs, emxArray_creal_T *hh)
{
    int y_size[2];
    int loop_ub;
    int i33;
    emxArray_real_T *r11;
    double y_data[1023];
    emxArray_real_T *digw;
    emxArray_creal_T *s;
    emxArray_creal_T *r12;
    double re;
    double im;

    /*  Cast to enforce precision rules */
    /*  Remaining are default or for advanced use */
    /*  Make b a row */
    /* -------------------------------------------------------------------------- */
    y_size[0] = 1;
    y_size[1] = b_size[1];
    loop_ub = b_size[1];
    for (i33 = 0; i33 < loop_ub; i33++) {
        y_data[i33] = b_data[i33];
    }

    emxInit_real_T(&r11, 2);

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
    i33 = r11->size[0] * r11->size[1];
    r11->size[0] = 1;
    r11->size[1] = w->size[1];
    emxEnsureCapacity((emxArray__common *)r11, i33, sizeof(double));
    loop_ub = w->size[0] * w->size[1];
    for (i33 = 0; i33 < loop_ub; i33++) {
        r11->data[i33] = 6.2831853071795862 * w->data[i33];
    }

    emxInit_real_T(&digw, 2);
    emxInit_creal_T(&s, 2);
    rdivide(r11, Fs, digw);

    /*  Convert from Hz to rad/sample for computational purposes */
    i33 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = digw->size[1];
    emxEnsureCapacity((emxArray__common *)s, i33, sizeof(creal_T));
    loop_ub = digw->size[0] * digw->size[1];
    emxFree_real_T(&r11);
    for (i33 = 0; i33 < loop_ub; i33++) {
        s->data[i33].re = digw->data[i33] * 0.0;
        s->data[i33].im = digw->data[i33];
    }

    emxInit_creal_T(&r12, 2);
    c_exp(s);

    /*  Digital frequency must be used for this calculation */
    j_polyval(y_data, y_size, s, r12);
    i33 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = digw->size[1];
    emxEnsureCapacity((emxArray__common *)s, i33, sizeof(creal_T));
    loop_ub = digw->size[0] * digw->size[1];
    for (i33 = 0; i33 < loop_ub; i33++) {
        re = digw->data[i33] * 0.0;
        im = digw->data[i33];
        s->data[i33].re = ((double)y_size[1] - 1.0) * re;
        s->data[i33].im = ((double)y_size[1] - 1.0) * im;
    }

    emxFree_real_T(&digw);
    c_exp(s);
    b_rdivide(r12, s, hh);

    /*  Generate the default structure to pass to freqzplot */
    /*  If rad/sample, Fs is empty */
    emxFree_creal_T(&r12);
    emxFree_creal_T(&s);
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
    emxArray_real_T *r13;
    int i35;
    int loop_ub;
    emxArray_real_T *digw;
    emxArray_creal_T *s;
    emxArray_creal_T *y;
    boolean_T b5;
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
    i35 = r13->size[0] * r13->size[1];
    r13->size[0] = 1;
    r13->size[1] = w->size[1];
    emxEnsureCapacity((emxArray__common *)r13, i35, sizeof(double));
    loop_ub = w->size[0] * w->size[1];
    for (i35 = 0; i35 < loop_ub; i35++) {
        r13->data[i35] = 6.2831853071795862 * w->data[i35];
    }

    emxInit_real_T(&digw, 2);
    emxInit_creal_T(&s, 2);
    rdivide(r13, Fs, digw);

    /*  Convert from Hz to rad/sample for computational purposes */
    i35 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = digw->size[1];
    emxEnsureCapacity((emxArray__common *)s, i35, sizeof(creal_T));
    loop_ub = digw->size[0] * digw->size[1];
    emxFree_real_T(&r13);
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
    emxEnsureCapacity((emxArray__common *)y, i35, sizeof(creal_T));
    b5 = (y->size[1] == 0);
    if (!b5) {
        i35 = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)y, i35, sizeof(creal_T));
        loop_ub = y->size[1];
        for (i35 = 0; i35 < loop_ub; i35++) {
            y->data[y->size[0] * i35].re = b[0];
            y->data[y->size[0] * i35].im = 0.0;
        }

        for (k = 0; k < 28; k++) {
            i35 = y->size[0] * y->size[1];
            y->size[0] = 1;
            y->size[1] = s->size[1];
            emxEnsureCapacity((emxArray__common *)y, i35, sizeof(creal_T));
            loop_ub = s->size[0] * s->size[1];
            for (i35 = 0; i35 < loop_ub; i35++) {
                re = s->data[i35].re * y->data[i35].re - s->data[i35].im * y->data[i35].
                     im;
                im = s->data[i35].re * y->data[i35].im + s->data[i35].im * y->data[i35].
                     re;
                y->data[i35].re = re + b[k + 1];
                y->data[i35].im = im;
            }
        }
    }

    i35 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = digw->size[1];
    emxEnsureCapacity((emxArray__common *)s, i35, sizeof(creal_T));
    loop_ub = digw->size[0] * digw->size[1];
    for (i35 = 0; i35 < loop_ub; i35++) {
        re = digw->data[i35] * 0.0;
        im = digw->data[i35];
        s->data[i35].re = 28.0 * re;
        s->data[i35].im = 28.0 * im;
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
    emxArray_real_T *r14;
    int i36;
    int loop_ub;
    emxArray_real_T *digw;
    emxArray_creal_T *s;
    emxArray_creal_T *y;
    boolean_T b6;
    int k;
    double re;
    double im;
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
    i36 = r14->size[0] * r14->size[1];
    r14->size[0] = 1;
    r14->size[1] = w->size[1];
    emxEnsureCapacity((emxArray__common *)r14, i36, sizeof(double));
    loop_ub = w->size[0] * w->size[1];
    for (i36 = 0; i36 < loop_ub; i36++) {
        r14->data[i36] = 6.2831853071795862 * w->data[i36];
    }

    emxInit_real_T(&digw, 2);
    emxInit_creal_T(&s, 2);
    rdivide(r14, Fs, digw);

    /*  Convert from Hz to rad/sample for computational purposes */
    i36 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = digw->size[1];
    emxEnsureCapacity((emxArray__common *)s, i36, sizeof(creal_T));
    loop_ub = digw->size[0] * digw->size[1];
    emxFree_real_T(&r14);
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
    emxEnsureCapacity((emxArray__common *)y, i36, sizeof(creal_T));
    b6 = (y->size[1] == 0);
    if (!b6) {
        i36 = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)y, i36, sizeof(creal_T));
        loop_ub = y->size[1];
        for (i36 = 0; i36 < loop_ub; i36++) {
            y->data[y->size[0] * i36].re = b[0];
            y->data[y->size[0] * i36].im = 0.0;
        }

        for (k = 0; k < 12; k++) {
            i36 = y->size[0] * y->size[1];
            y->size[0] = 1;
            y->size[1] = s->size[1];
            emxEnsureCapacity((emxArray__common *)y, i36, sizeof(creal_T));
            loop_ub = s->size[0] * s->size[1];
            for (i36 = 0; i36 < loop_ub; i36++) {
                re = s->data[i36].re * y->data[i36].re - s->data[i36].im * y->data[i36].
                     im;
                im = s->data[i36].re * y->data[i36].im + s->data[i36].im * y->data[i36].
                     re;
                y->data[i36].re = re + b[k + 1];
                y->data[i36].im = im;
            }
        }
    }

    i36 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = digw->size[1];
    emxEnsureCapacity((emxArray__common *)s, i36, sizeof(creal_T));
    loop_ub = digw->size[0] * digw->size[1];
    for (i36 = 0; i36 < loop_ub; i36++) {
        re = digw->data[i36] * 0.0;
        im = digw->data[i36];
        s->data[i36].re = 12.0 * re;
        s->data[i36].im = 12.0 * im;
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
    emxArray_creal_T *r6;
    emxInit_creal_T1(&r6, 1);
    eig(x, r6);
    vector_poly(r6, c);
    emxFree_creal_T(&r6);
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
        y[k] = a[k] * a[k];
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
    emxArray_real_T *r15;
    int i37;
    int loop_ub;
    emxArray_real_T *digw;
    emxArray_creal_T *s;
    emxArray_creal_T *y;
    boolean_T b7;
    int k;
    double re;
    double im;
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
    i37 = r15->size[0] * r15->size[1];
    r15->size[0] = 1;
    r15->size[1] = w->size[1];
    emxEnsureCapacity((emxArray__common *)r15, i37, sizeof(double));
    loop_ub = w->size[0] * w->size[1];
    for (i37 = 0; i37 < loop_ub; i37++) {
        r15->data[i37] = 6.2831853071795862 * w->data[i37];
    }

    emxInit_real_T(&digw, 2);
    emxInit_creal_T(&s, 2);
    rdivide(r15, Fs, digw);

    /*  Convert from Hz to rad/sample for computational purposes */
    i37 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = digw->size[1];
    emxEnsureCapacity((emxArray__common *)s, i37, sizeof(creal_T));
    loop_ub = digw->size[0] * digw->size[1];
    emxFree_real_T(&r15);
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
    emxEnsureCapacity((emxArray__common *)y, i37, sizeof(creal_T));
    b7 = (y->size[1] == 0);
    if (!b7) {
        i37 = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)y, i37, sizeof(creal_T));
        loop_ub = y->size[1];
        for (i37 = 0; i37 < loop_ub; i37++) {
            y->data[y->size[0] * i37].re = b[0];
            y->data[y->size[0] * i37].im = 0.0;
        }

        for (k = 0; k < 56; k++) {
            i37 = y->size[0] * y->size[1];
            y->size[0] = 1;
            y->size[1] = s->size[1];
            emxEnsureCapacity((emxArray__common *)y, i37, sizeof(creal_T));
            loop_ub = s->size[0] * s->size[1];
            for (i37 = 0; i37 < loop_ub; i37++) {
                re = s->data[i37].re * y->data[i37].re - s->data[i37].im * y->data[i37].
                     im;
                im = s->data[i37].re * y->data[i37].im + s->data[i37].im * y->data[i37].
                     re;
                y->data[i37].re = re + b[k + 1];
                y->data[i37].im = im;
            }
        }
    }

    i37 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = digw->size[1];
    emxEnsureCapacity((emxArray__common *)s, i37, sizeof(creal_T));
    loop_ub = digw->size[0] * digw->size[1];
    for (i37 = 0; i37 < loop_ub; i37++) {
        re = digw->data[i37] * 0.0;
        im = digw->data[i37];
        s->data[i37].re = 56.0 * re;
        s->data[i37].im = 56.0 * im;
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
    emxArray_real_T *r16;
    int i38;
    int loop_ub;
    emxArray_real_T *digw;
    emxArray_creal_T *s;
    emxArray_creal_T *y;
    boolean_T b8;
    int k;
    double re;
    double im;
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
    i38 = r16->size[0] * r16->size[1];
    r16->size[0] = 1;
    r16->size[1] = w->size[1];
    emxEnsureCapacity((emxArray__common *)r16, i38, sizeof(double));
    loop_ub = w->size[0] * w->size[1];
    for (i38 = 0; i38 < loop_ub; i38++) {
        r16->data[i38] = 6.2831853071795862 * w->data[i38];
    }

    emxInit_real_T(&digw, 2);
    emxInit_creal_T(&s, 2);
    rdivide(r16, Fs, digw);

    /*  Convert from Hz to rad/sample for computational purposes */
    i38 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = digw->size[1];
    emxEnsureCapacity((emxArray__common *)s, i38, sizeof(creal_T));
    loop_ub = digw->size[0] * digw->size[1];
    emxFree_real_T(&r16);
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
    emxEnsureCapacity((emxArray__common *)y, i38, sizeof(creal_T));
    b8 = (y->size[1] == 0);
    if (!b8) {
        i38 = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)y, i38, sizeof(creal_T));
        loop_ub = y->size[1];
        for (i38 = 0; i38 < loop_ub; i38++) {
            y->data[y->size[0] * i38].re = b[0];
            y->data[y->size[0] * i38].im = 0.0;
        }

        for (k = 0; k < 42; k++) {
            i38 = y->size[0] * y->size[1];
            y->size[0] = 1;
            y->size[1] = s->size[1];
            emxEnsureCapacity((emxArray__common *)y, i38, sizeof(creal_T));
            loop_ub = s->size[0] * s->size[1];
            for (i38 = 0; i38 < loop_ub; i38++) {
                re = s->data[i38].re * y->data[i38].re - s->data[i38].im * y->data[i38].
                     im;
                im = s->data[i38].re * y->data[i38].im + s->data[i38].im * y->data[i38].
                     re;
                y->data[i38].re = re + b[k + 1];
                y->data[i38].im = im;
            }
        }
    }

    i38 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = digw->size[1];
    emxEnsureCapacity((emxArray__common *)s, i38, sizeof(creal_T));
    loop_ub = digw->size[0] * digw->size[1];
    for (i38 = 0; i38 < loop_ub; i38++) {
        re = digw->data[i38] * 0.0;
        im = digw->data[i38];
        s->data[i38].re = 42.0 * re;
        s->data[i38].im = 42.0 * im;
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
    int i28;
    int loop_ub;
    i28 = z->size[0] * z->size[1];
    z->size[0] = 1;
    z->size[1] = x->size[1];
    emxEnsureCapacity((emxArray__common *)z, i28, sizeof(double));
    loop_ub = x->size[0] * x->size[1];
    for (i28 = 0; i28 < loop_ub; i28++) {
        z->data[i28] = x->data[i28] / y;
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
    emxArray_int32_T *r29;
    int i52;
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
    emxInit_int32_T(&r29, 2);
    while (l <= (int)m - 1) {
        if ((m == 0.0) || (((m > 0.0) && (1.0 + (double)l > n)) || ((0.0 > m) && (n >
                           1.0 + (double)l)))) {
            i52 = 1;
            i = 1;
            end = 0;
        } else {
            i52 = l + 1;
            i = (int)m;
            end = (int)n;
        }

        b_x = x->data[(int)k - 1];
        loop_ub = xx->size[0] * xx->size[1];
        xx->size[0] = 1;
        xx->size[1] = div_s32_floor(end - i52, i) + 1;
        emxEnsureCapacity((emxArray__common *)xx, loop_ub, sizeof(double));
        loop_ub = div_s32_floor(end - i52, i);
        for (end = 0; end <= loop_ub; end++) {
            xx->data[xx->size[0] * end] = 2.0 * (b_x - x->data[(i52 + i * end) - 1]);
        }

        end = xx->size[1] - 1;
        loop_ub = 0;
        for (i = 0; i <= end; i++) {
            if (xx->data[i] != 0.0) {
                loop_ub++;
            }
        }

        i52 = r29->size[0] * r29->size[1];
        r29->size[0] = 1;
        r29->size[1] = loop_ub;
        emxEnsureCapacity((emxArray__common *)r29, i52, sizeof(int));
        loop_ub = 0;
        for (i = 0; i <= end; i++) {
            if (xx->data[i] != 0.0) {
                r29->data[loop_ub] = i + 1;
                loop_ub++;
            }
        }

        if (r29->size[1] == 0) {
            b_x = 1.0;
        } else {
            b_x = xx->data[r29->data[0] - 1];
            for (loop_ub = 2; loop_ub <= r29->size[1]; loop_ub++) {
                b_x *= xx->data[r29->data[r29->size[0] * (loop_ub - 1)] - 1];
            }
        }

        y *= b_x;
        l++;
    }

    emxFree_int32_T(&r29);
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
    int i45;
    double temp;
    int loop_ub;
    emxArray_real_T *iext;
    emxArray_real_T *r26;
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
    double b_dev;
    double devl;
    emxArray_real_T *ad;
    int niter;
    double jchnge;
    double d1;
    emxArray_real_T *l;
    emxArray_real_T *a;
    emxArray_int32_T *r27;
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
    emxArray_int8_T *r28;
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
    boolean_T guard1 = false;
    int exitg1;
    double d_j;
    double c_k1;
    int nut;
    double knz;
    double b_l;
    int nu;
    double varargin_1[2];
    double dnum;
    int i46;
    int i47;
    int i48;
    int i49;
    int i50;
    boolean_T guard2 = false;
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
        i45 = x2->size[0] * x2->size[1];
        x2->size[0] = 1;
        x2->size[1] = grid->size[1];
        emxEnsureCapacity((emxArray__common *)x2, i45, sizeof(double));
        loop_ub = grid->size[0] * grid->size[1];
        for (i45 = 0; i45 < loop_ub; i45++) {
            x2->data[i45] = 3.1415926535897931 * grid->data[i45];
        }

        b_cos(x2);
        c_rdivide(des, x2, j);
        i45 = des->size[0] * des->size[1];
        des->size[0] = 1;
        des->size[1] = j->size[1];
        emxEnsureCapacity((emxArray__common *)des, i45, sizeof(double));
        loop_ub = j->size[0] * j->size[1];
        for (i45 = 0; i45 < loop_ub; i45++) {
            des->data[i45] = j->data[i45];
        }

        i45 = x2->size[0] * x2->size[1];
        x2->size[0] = 1;
        x2->size[1] = grid->size[1];
        emxEnsureCapacity((emxArray__common *)x2, i45, sizeof(double));
        loop_ub = grid->size[0] * grid->size[1];
        for (i45 = 0; i45 < loop_ub; i45++) {
            x2->data[i45] = 3.1415926535897931 * grid->data[i45];
        }

        emxInit_real_T(&r26, 2);
        b_cos(x2);
        i45 = r26->size[0] * r26->size[1];
        r26->size[0] = 1;
        r26->size[1] = x2->size[1];
        emxEnsureCapacity((emxArray__common *)r26, i45, sizeof(double));
        loop_ub = x2->size[0] * x2->size[1];
        for (i45 = 0; i45 < loop_ub; i45++) {
            r26->data[i45] = x2->data[i45];
        }

        i45 = wt->size[0] * wt->size[1];
        wt->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)wt, i45, sizeof(double));
        ixstart = wt->size[0];
        flag = wt->size[1];
        loop_ub = ixstart * flag;
        for (i45 = 0; i45 < loop_ub; i45++) {
            wt->data[i45] *= r26->data[i45];
        }

        emxFree_real_T(&r26);
    }

    temp = ((double)grid->size[1] - 1.0) / nfcns;
    if (rtIsNaN(nfcns)) {
        i45 = j->size[0] * j->size[1];
        j->size[0] = 1;
        j->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)j, i45, sizeof(double));
        j->data[0] = rtNaN;
    } else if (nfcns < 1.0) {
        i45 = j->size[0] * j->size[1];
        j->size[0] = 1;
        j->size[1] = 0;
        emxEnsureCapacity((emxArray__common *)j, i45, sizeof(double));
    } else if (rtIsInf(nfcns) && (1.0 == nfcns)) {
        i45 = j->size[0] * j->size[1];
        j->size[0] = 1;
        j->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)j, i45, sizeof(double));
        j->data[0] = rtNaN;
    } else {
        i45 = j->size[0] * j->size[1];
        j->size[0] = 1;
        j->size[1] = (int)floor(nfcns - 1.0) + 1;
        emxEnsureCapacity((emxArray__common *)j, i45, sizeof(double));
        loop_ub = (int)floor(nfcns - 1.0);
        for (i45 = 0; i45 <= loop_ub; i45++) {
            j->data[j->size[0] * i45] = 1.0 + (double)i45;
        }
    }

    i45 = x2->size[0] * x2->size[1];
    x2->size[0] = 1;
    x2->size[1] = j->size[1] + 1;
    emxEnsureCapacity((emxArray__common *)x2, i45, sizeof(double));
    loop_ub = j->size[1];
    for (i45 = 0; i45 < loop_ub; i45++) {
        x2->data[x2->size[0] * i45] = temp * (j->data[j->size[0] * i45] - 1.0) + 1.0;
    }

    emxInit_real_T1(&iext, 1);
    x2->data[x2->size[0] * j->size[1]] = grid->size[1];
    c_fix(x2);
    i45 = iext->size[0];
    iext->size[0] = x2->size[1] + 1;
    emxEnsureCapacity((emxArray__common *)iext, i45, sizeof(double));
    loop_ub = x2->size[1];
    for (i45 = 0; i45 < loop_ub; i45++) {
        iext->data[i45] = x2->data[x2->size[0] * i45];
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
    i45 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)y, i45, sizeof(double));
    y->data[0] = -1.0;
    b_dev = -1.0;
    devl = -1.0;
    i45 = x->size[0] * x->size[1];
    x->size[0] = 1;
    x->size[1] = (int)(nfcns + 1.0);
    emxEnsureCapacity((emxArray__common *)x, i45, sizeof(double));
    loop_ub = (int)(nfcns + 1.0);
    for (i45 = 0; i45 < loop_ub; i45++) {
        x->data[i45] = 0.0;
    }

    emxInit_real_T(&ad, 2);
    niter = 0;
    jchnge = 1.0;
    d1 = (nfcns - 1.0) / 15.0;
    b_fix(&d1);
    i45 = ad->size[0] * ad->size[1];
    ad->size[0] = 1;
    ad->size[1] = (int)(nfcns + 1.0);
    emxEnsureCapacity((emxArray__common *)ad, i45, sizeof(double));
    loop_ub = (int)(nfcns + 1.0);
    for (i45 = 0; i45 < loop_ub; i45++) {
        ad->data[i45] = 0.0;
    }

    /*  index manager(s) */
    emxInit_real_T(&l, 2);
    emxInit_real_T(&a, 2);
    emxInit_int32_T(&r27, 2);
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
    emxInit_int8_T(&r28, 2);
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
    emxInit_real_T1(&c_j, 1);
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

                i45 = l->size[0] * l->size[1];
                l->size[0] = 1;
                l->size[1] = loop_ub;
                emxEnsureCapacity((emxArray__common *)l, i45, sizeof(double));
                for (i45 = 0; i45 < loop_ub; i45++) {
                    l->data[l->size[0] * i45] = iext->data[i45];
                }

                i45 = x->size[0] * x->size[1];
                x->size[0] = 1;
                x->size[1] = l->size[1];
                emxEnsureCapacity((emxArray__common *)x, i45, sizeof(double));
                nut = l->size[0] * l->size[1];
                for (i45 = 0; i45 < nut; i45++) {
                    x->data[i45] = 6.2831853071795862 * grid->data[(int)l->data[i45] - 1];
                }

                b_cos(x);
                for (ixstart = 0; ixstart < (int)(nfcns + 1.0); ixstart++) {
                    ad->data[ixstart] = remezdd(1.0 + (double)ixstart, nfcns + 1.0, d1 +
                                                1.0, x);
                }

                for (i45 = 0; i45 < 2; i45++) {
                    varargin_1[i45] = ad->size[i45];
                }

                i45 = j->size[0] * j->size[1];
                j->size[0] = 1;
                j->size[1] = (int)varargin_1[1];
                emxEnsureCapacity((emxArray__common *)j, i45, sizeof(double));
                nut = (int)varargin_1[1];
                for (i45 = 0; i45 < nut; i45++) {
                    j->data[i45] = 1.0;
                }

                if (2.0 > nfcns + 1.0) {
                    i45 = 0;
                    i46 = 1;
                    i47 = 0;
                    i48 = 0;
                    i49 = 1;
                } else {
                    i45 = 1;
                    i46 = 2;
                    i47 = (int)(nfcns + 1.0);
                    i48 = 1;
                    i49 = 2;
                }

                i50 = r28->size[0] * r28->size[1];
                r28->size[0] = 1;
                r28->size[1] = (int)varargin_1[1];
                emxEnsureCapacity((emxArray__common *)r28, i50, sizeof(signed char));
                nut = (int)varargin_1[1];
                for (i50 = 0; i50 < nut; i50++) {
                    r28->data[r28->size[0] * i50] = 1;
                }

                nut = div_s32_floor((i47 - i45) - 1, i46);
                for (i47 = 0; i47 <= nut; i47++) {
                    j->data[i48 + i49 * i47] = -(double)r28->data[i45 + i46 * i47];
                }

                i45 = b->size[0];
                b->size[0] = l->size[1];
                emxEnsureCapacity((emxArray__common *)b, i45, sizeof(double));
                nut = l->size[1];
                for (i45 = 0; i45 < nut; i45++) {
                    b->data[i45] = des->data[(int)l->data[l->size[0] * i45] - 1];
                }

                guard2 = false;
                if (ad->size[1] == 1) {
                    guard2 = true;
                } else {
                    i45 = b_iext->size[0];
                    b_iext->size[0] = loop_ub;
                    emxEnsureCapacity((emxArray__common *)b_iext, i45, sizeof(int));
                    for (i45 = 0; i45 < loop_ub; i45++) {
                        b_iext->data[i45] = (int)iext->data[i45];
                    }

                    if (b_iext->size[0] == 1) {
                        guard2 = true;
                    } else {
                        dnum = 0.0;
                        for (i45 = 0; i45 < ad->size[1]; i45++) {
                            dnum += ad->data[ad->size[0] * i45] * b->data[i45];
                        }
                    }
                }

                if (guard2) {
                    dnum = 0.0;
                    for (i45 = 0; i45 < ad->size[1]; i45++) {
                        dnum += ad->data[ad->size[0] * i45] * b->data[i45];
                    }
                }

                i45 = c_wt->size[0] * c_wt->size[1];
                c_wt->size[0] = 1;
                c_wt->size[1] = l->size[1];
                emxEnsureCapacity((emxArray__common *)c_wt, i45, sizeof(double));
                loop_ub = l->size[0] * l->size[1];
                for (i45 = 0; i45 < loop_ub; i45++) {
                    c_wt->data[i45] = wt->data[(int)l->data[i45] - 1];
                }

                c_rdivide(ad, c_wt, x2);
                i45 = b->size[0];
                b->size[0] = x2->size[1];
                emxEnsureCapacity((emxArray__common *)b, i45, sizeof(double));
                loop_ub = x2->size[1];
                for (i45 = 0; i45 < loop_ub; i45++) {
                    b->data[i45] = x2->data[x2->size[0] * i45];
                }

                if ((j->size[1] == 1) || (b->size[0] == 1)) {
                    temp = 0.0;
                    for (i45 = 0; i45 < j->size[1]; i45++) {
                        temp += j->data[j->size[0] * i45] * b->data[i45];
                    }
                } else {
                    temp = 0.0;
                    for (i45 = 0; i45 < j->size[1]; i45++) {
                        temp += j->data[j->size[0] * i45] * b->data[i45];
                    }
                }

                b_dev = dnum / temp;
                nu = 1;
                if (b_dev > 0.0) {
                    nu = -1;
                }

                b_dev *= -(double)nu;
                temp = (double)nu * b_dev;
                i45 = b_a->size[0] * b_a->size[1];
                b_a->size[0] = 1;
                b_a->size[1] = j->size[1];
                emxEnsureCapacity((emxArray__common *)b_a, i45, sizeof(double));
                loop_ub = j->size[0] * j->size[1];
                for (i45 = 0; i45 < loop_ub; i45++) {
                    b_a->data[i45] = temp * j->data[i45];
                }

                i45 = b_wt->size[0] * b_wt->size[1];
                b_wt->size[0] = 1;
                b_wt->size[1] = l->size[1];
                emxEnsureCapacity((emxArray__common *)b_wt, i45, sizeof(double));
                loop_ub = l->size[0] * l->size[1];
                for (i45 = 0; i45 < loop_ub; i45++) {
                    b_wt->data[i45] = wt->data[(int)l->data[i45] - 1];
                }

                c_rdivide(b_a, b_wt, x2);
                i45 = y->size[0] * y->size[1];
                y->size[0] = 1;
                y->size[1] = l->size[1];
                emxEnsureCapacity((emxArray__common *)y, i45, sizeof(double));
                loop_ub = l->size[0] * l->size[1];
                for (i45 = 0; i45 < loop_ub; i45++) {
                    y->data[i45] = des->data[(int)l->data[i45] - 1] + x2->data[i45];
                }

                if (b_dev <= devl) {
                    /* warning(message('signal:firpm:DidNotConverge',niter)) */
                    cfprintf("DidNotConverge");
                    i45 = h->size[0] * h->size[1];
                    h->size[0] = (int)nfilt;
                    h->size[1] = 1;
                    emxEnsureCapacity((emxArray__common *)h, i45, sizeof(double));
                    loop_ub = (int)nfilt;
                    for (i45 = 0; i45 < loop_ub; i45++) {
                        h->data[i45] = 0.0;
                    }

                    b_dev = -1.0;

                    /* iext */
                    b_valid = false;
                    exitg1 = 1;
                } else {
                    devl = b_dev;
                    jchnge = 0.0;
                    c_k1 = iext->data[0];
                    knz = iext->data[(int)(nfcns + 1.0) - 1];
                    temp = 0.0;
                    nut = -nu;
                    d_j = 1.0;
                    flag34 = 1;
                    while (d_j < (nfcns + 1.0) + 1.0) {
                        dnum = iext->data[(int)(unsigned int)d_j];
                        b_l = iext->data[(int)d_j - 1] + 1.0;
                        nut = -nut;
                        if (d_j == 2.0) {
                            b_y1 = comp;
                        }

                        comp = b_dev;
                        flag = 1;
                        if (iext->data[(int)d_j - 1] + 1.0 < iext->data[(int)(d_j + 1.0) - 1]) {
                            /*  gee */
                            err = cos(6.2831853071795862 * grid->data[(int)(iext->data[(int)
                                      d_j - 1] + 1.0) - 1]);
                            i45 = n_x->size[0] * n_x->size[1];
                            n_x->size[0] = 1;
                            n_x->size[1] = x->size[1];
                            emxEnsureCapacity((emxArray__common *)n_x, i45, sizeof(double));
                            loop_ub = x->size[0] * x->size[1];
                            for (i45 = 0; i45 < loop_ub; i45++) {
                                n_x->data[i45] = err - x->data[i45];
                            }

                            c_rdivide(ad, n_x, j);
                            i45 = b->size[0];
                            b->size[0] = y->size[1];
                            emxEnsureCapacity((emxArray__common *)b, i45, sizeof(double));
                            loop_ub = y->size[1];
                            for (i45 = 0; i45 < loop_ub; i45++) {
                                b->data[i45] = y->data[y->size[0] * i45];
                            }

                            if ((j->size[1] == 1) || (b->size[0] == 1)) {
                                dtemp = 0.0;
                                for (i45 = 0; i45 < j->size[1]; i45++) {
                                    dtemp += j->data[j->size[0] * i45] * b->data[i45];
                                }
                            } else {
                                dtemp = 0.0;
                                for (i45 = 0; i45 < j->size[1]; i45++) {
                                    dtemp += j->data[j->size[0] * i45] * b->data[i45];
                                }
                            }

                            err = b_sum(j);
                            err = (dtemp / err - des->data[(int)(iext->data[(int)d_j - 1] +
                                                                 1.0) - 1]) * wt->data[(int)(iext->data[(int)d_j - 1] + 1.0)
                                                                         - 1];
                            dtemp = (double)nut * err - b_dev;
                            if (dtemp > 0.0) {
                                comp = (double)nut * err;
                                b_l = (iext->data[(int)d_j - 1] + 1.0) + 1.0;
                                exitg2 = false;
                                while ((!exitg2) && (b_l < dnum)) {
                                    /*  gee */
                                    err = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                                    i45 = m_x->size[0] * m_x->size[1];
                                    m_x->size[0] = 1;
                                    m_x->size[1] = x->size[1];
                                    emxEnsureCapacity((emxArray__common *)m_x, i45, sizeof(double));
                                    loop_ub = x->size[0] * x->size[1];
                                    for (i45 = 0; i45 < loop_ub; i45++) {
                                        m_x->data[i45] = err - x->data[i45];
                                    }

                                    c_rdivide(ad, m_x, j);
                                    i45 = b->size[0];
                                    b->size[0] = y->size[1];
                                    emxEnsureCapacity((emxArray__common *)b, i45, sizeof(double));
                                    loop_ub = y->size[1];
                                    for (i45 = 0; i45 < loop_ub; i45++) {
                                        b->data[i45] = y->data[y->size[0] * i45];
                                    }

                                    if ((j->size[1] == 1) || (b->size[0] == 1)) {
                                        dtemp = 0.0;
                                        for (i45 = 0; i45 < j->size[1]; i45++) {
                                            dtemp += j->data[j->size[0] * i45] * b->data[i45];
                                        }
                                    } else {
                                        dtemp = 0.0;
                                        for (i45 = 0; i45 < j->size[1]; i45++) {
                                            dtemp += j->data[j->size[0] * i45] * b->data[i45];
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
                                        exitg2 = true;
                                    }
                                }

                                iext->data[(int)d_j - 1] = b_l - 1.0;
                                d_j++;
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
                                err = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                                i45 = l_x->size[0] * l_x->size[1];
                                l_x->size[0] = 1;
                                l_x->size[1] = x->size[1];
                                emxEnsureCapacity((emxArray__common *)l_x, i45, sizeof(double));
                                loop_ub = x->size[0] * x->size[1];
                                for (i45 = 0; i45 < loop_ub; i45++) {
                                    l_x->data[i45] = err - x->data[i45];
                                }

                                c_rdivide(ad, l_x, j);
                                i45 = b->size[0];
                                b->size[0] = y->size[1];
                                emxEnsureCapacity((emxArray__common *)b, i45, sizeof(double));
                                loop_ub = y->size[1];
                                for (i45 = 0; i45 < loop_ub; i45++) {
                                    b->data[i45] = y->data[y->size[0] * i45];
                                }

                                if ((j->size[1] == 1) || (b->size[0] == 1)) {
                                    dtemp = 0.0;
                                    for (i45 = 0; i45 < j->size[1]; i45++) {
                                        dtemp += j->data[j->size[0] * i45] * b->data[i45];
                                    }
                                } else {
                                    dtemp = 0.0;
                                    for (i45 = 0; i45 < j->size[1]; i45++) {
                                        dtemp += j->data[j->size[0] * i45] * b->data[i45];
                                    }
                                }

                                err = b_sum(j);
                                err = (dtemp / err - des->data[(int)b_l - 1]) * wt->data[(int)
                                        b_l - 1];
                                dtemp = (double)nut * err - comp;
                                if ((dtemp > 0.0) || (jchnge > 0.0)) {
                                    exitg2 = true;
                                } else {
                                    b_l--;
                                }
                            }

                            if (b_l <= temp) {
                                b_l = iext->data[(int)d_j - 1] + 1.0;
                                if (jchnge > 0.0) {
                                    iext->data[(int)d_j - 1] = (iext->data[(int)d_j - 1] + 1.0) -
                                                               1.0;
                                    d_j++;
                                    temp = b_l - 1.0;
                                    jchnge++;
                                } else {
                                    b_l = (iext->data[(int)d_j - 1] + 1.0) + 1.0;
                                    exitg2 = false;
                                    while ((!exitg2) && (b_l < dnum)) {
                                        /*  gee */
                                        err = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                                        i45 = k_x->size[0] * k_x->size[1];
                                        k_x->size[0] = 1;
                                        k_x->size[1] = x->size[1];
                                        emxEnsureCapacity((emxArray__common *)k_x, i45, sizeof
                                                          (double));
                                        loop_ub = x->size[0] * x->size[1];
                                        for (i45 = 0; i45 < loop_ub; i45++) {
                                            k_x->data[i45] = err - x->data[i45];
                                        }

                                        c_rdivide(ad, k_x, j);
                                        i45 = b->size[0];
                                        b->size[0] = y->size[1];
                                        emxEnsureCapacity((emxArray__common *)b, i45, sizeof(double));
                                        loop_ub = y->size[1];
                                        for (i45 = 0; i45 < loop_ub; i45++) {
                                            b->data[i45] = y->data[y->size[0] * i45];
                                        }

                                        if ((j->size[1] == 1) || (b->size[0] == 1)) {
                                            dtemp = 0.0;
                                            for (i45 = 0; i45 < j->size[1]; i45++) {
                                                dtemp += j->data[j->size[0] * i45] * b->data[i45];
                                            }
                                        } else {
                                            dtemp = 0.0;
                                            for (i45 = 0; i45 < j->size[1]; i45++) {
                                                dtemp += j->data[j->size[0] * i45] * b->data[i45];
                                            }
                                        }

                                        err = b_sum(j);
                                        err = (dtemp / err - des->data[(int)b_l - 1]) * wt->data
                                              [(int)b_l - 1];
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
                                            err = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                                            i45 = j_x->size[0] * j_x->size[1];
                                            j_x->size[0] = 1;
                                            j_x->size[1] = x->size[1];
                                            emxEnsureCapacity((emxArray__common *)j_x, i45, sizeof
                                                              (double));
                                            loop_ub = x->size[0] * x->size[1];
                                            for (i45 = 0; i45 < loop_ub; i45++) {
                                                j_x->data[i45] = err - x->data[i45];
                                            }

                                            c_rdivide(ad, j_x, j);
                                            i45 = b->size[0];
                                            b->size[0] = y->size[1];
                                            emxEnsureCapacity((emxArray__common *)b, i45, sizeof
                                                              (double));
                                            loop_ub = y->size[1];
                                            for (i45 = 0; i45 < loop_ub; i45++) {
                                                b->data[i45] = y->data[y->size[0] * i45];
                                            }

                                            if ((j->size[1] == 1) || (b->size[0] == 1)) {
                                                dtemp = 0.0;
                                                for (i45 = 0; i45 < j->size[1]; i45++) {
                                                    dtemp += j->data[j->size[0] * i45] * b->data[i45];
                                                }
                                            } else {
                                                dtemp = 0.0;
                                                for (i45 = 0; i45 < j->size[1]; i45++) {
                                                    dtemp += j->data[j->size[0] * i45] * b->data[i45];
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
                                                exitg2 = true;
                                            }
                                        }

                                        iext->data[(int)d_j - 1] = b_l - 1.0;
                                        d_j++;
                                        temp = b_l - 1.0;
                                        jchnge = 1.0;
                                    } else {
                                        temp = iext->data[(int)d_j - 1];
                                        d_j++;
                                    }
                                }
                            } else if (dtemp > 0.0) {
                                comp = (double)nut * err;
                                b_l--;
                                exitg2 = false;
                                while ((!exitg2) && (b_l > temp)) {
                                    /*  gee */
                                    err = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                                    i45 = i_x->size[0] * i_x->size[1];
                                    i_x->size[0] = 1;
                                    i_x->size[1] = x->size[1];
                                    emxEnsureCapacity((emxArray__common *)i_x, i45, sizeof(double));
                                    loop_ub = x->size[0] * x->size[1];
                                    for (i45 = 0; i45 < loop_ub; i45++) {
                                        i_x->data[i45] = err - x->data[i45];
                                    }

                                    c_rdivide(ad, i_x, j);
                                    i45 = b->size[0];
                                    b->size[0] = y->size[1];
                                    emxEnsureCapacity((emxArray__common *)b, i45, sizeof(double));
                                    loop_ub = y->size[1];
                                    for (i45 = 0; i45 < loop_ub; i45++) {
                                        b->data[i45] = y->data[y->size[0] * i45];
                                    }

                                    if ((j->size[1] == 1) || (b->size[0] == 1)) {
                                        dtemp = 0.0;
                                        for (i45 = 0; i45 < j->size[1]; i45++) {
                                            dtemp += j->data[j->size[0] * i45] * b->data[i45];
                                        }
                                    } else {
                                        dtemp = 0.0;
                                        for (i45 = 0; i45 < j->size[1]; i45++) {
                                            dtemp += j->data[j->size[0] * i45] * b->data[i45];
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
                                        exitg2 = true;
                                    }
                                }

                                temp = iext->data[(int)d_j - 1];
                                iext->data[(int)d_j - 1] = b_l + 1.0;
                                d_j++;
                                jchnge++;
                            } else {
                                temp = iext->data[(int)d_j - 1];
                                d_j++;
                            }
                        }
                    }

                    do {
                        exitg3 = 0;
                        if (d_j == (nfcns + 1.0) + 1.0) {
                            varargin_1[1] = iext->data[0];
                            ixstart = 1;
                            if (rtIsNaN(c_k1)) {
                                flag = 2;
                                exitg2 = false;
                                while ((!exitg2) && (flag < 3)) {
                                    ixstart = 2;
                                    if (!rtIsNaN(varargin_1[1])) {
                                        c_k1 = varargin_1[1];
                                        exitg2 = true;
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
                                exitg2 = false;
                                while ((!exitg2) && (flag < 3)) {
                                    ixstart = 2;
                                    if (!rtIsNaN(varargin_1[1])) {
                                        knz = varargin_1[1];
                                        exitg2 = true;
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
                            exitg2 = false;
                            while ((!exitg2) && (b_l < c_k1)) {
                                /*  gee */
                                err = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                                i45 = h_x->size[0] * h_x->size[1];
                                h_x->size[0] = 1;
                                h_x->size[1] = x->size[1];
                                emxEnsureCapacity((emxArray__common *)h_x, i45, sizeof(double));
                                loop_ub = x->size[0] * x->size[1];
                                for (i45 = 0; i45 < loop_ub; i45++) {
                                    h_x->data[i45] = err - x->data[i45];
                                }

                                c_rdivide(ad, h_x, j);
                                i45 = b->size[0];
                                b->size[0] = y->size[1];
                                emxEnsureCapacity((emxArray__common *)b, i45, sizeof(double));
                                loop_ub = y->size[1];
                                for (i45 = 0; i45 < loop_ub; i45++) {
                                    b->data[i45] = y->data[y->size[0] * i45];
                                }

                                if ((j->size[1] == 1) || (b->size[0] == 1)) {
                                    dtemp = 0.0;
                                    for (i45 = 0; i45 < j->size[1]; i45++) {
                                        dtemp += j->data[j->size[0] * i45] * b->data[i45];
                                    }
                                } else {
                                    dtemp = 0.0;
                                    for (i45 = 0; i45 < j->size[1]; i45++) {
                                        dtemp += j->data[j->size[0] * i45] * b->data[i45];
                                    }
                                }

                                err = b_sum(j);
                                err = (dtemp / err - des->data[(int)b_l - 1]) * wt->data[(int)
                                        b_l - 1];
                                dtemp = err * -(double)nu - comp;
                                if (dtemp > 0.0) {
                                    comp = -(double)nu * err;
                                    b_l++;
                                    exitg4 = false;
                                    while ((!exitg4) && (b_l < c_k1)) {
                                        /*  gee */
                                        err = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                                        i45 = g_x->size[0] * g_x->size[1];
                                        g_x->size[0] = 1;
                                        g_x->size[1] = x->size[1];
                                        emxEnsureCapacity((emxArray__common *)g_x, i45, sizeof
                                                          (double));
                                        loop_ub = x->size[0] * x->size[1];
                                        for (i45 = 0; i45 < loop_ub; i45++) {
                                            g_x->data[i45] = err - x->data[i45];
                                        }

                                        c_rdivide(ad, g_x, j);
                                        i45 = b->size[0];
                                        b->size[0] = y->size[1];
                                        emxEnsureCapacity((emxArray__common *)b, i45, sizeof(double));
                                        loop_ub = y->size[1];
                                        for (i45 = 0; i45 < loop_ub; i45++) {
                                            b->data[i45] = y->data[y->size[0] * i45];
                                        }

                                        if ((j->size[1] == 1) || (b->size[0] == 1)) {
                                            dtemp = 0.0;
                                            for (i45 = 0; i45 < j->size[1]; i45++) {
                                                dtemp += j->data[j->size[0] * i45] * b->data[i45];
                                            }
                                        } else {
                                            dtemp = 0.0;
                                            for (i45 = 0; i45 < j->size[1]; i45++) {
                                                dtemp += j->data[j->size[0] * i45] * b->data[i45];
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
                                            exitg4 = true;
                                        }
                                    }

                                    iext->data[(int)((nfcns + 1.0) + 1.0) - 1] = b_l - 1.0;
                                    d_j = ((nfcns + 1.0) + 1.0) + 1.0;
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
                                    err = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                                    i45 = f_x->size[0] * f_x->size[1];
                                    f_x->size[0] = 1;
                                    f_x->size[1] = x->size[1];
                                    emxEnsureCapacity((emxArray__common *)f_x, i45, sizeof(double));
                                    loop_ub = x->size[0] * x->size[1];
                                    for (i45 = 0; i45 < loop_ub; i45++) {
                                        f_x->data[i45] = err - x->data[i45];
                                    }

                                    c_rdivide(ad, f_x, j);
                                    i45 = b->size[0];
                                    b->size[0] = y->size[1];
                                    emxEnsureCapacity((emxArray__common *)b, i45, sizeof(double));
                                    loop_ub = y->size[1];
                                    for (i45 = 0; i45 < loop_ub; i45++) {
                                        b->data[i45] = y->data[y->size[0] * i45];
                                    }

                                    if ((j->size[1] == 1) || (b->size[0] == 1)) {
                                        dtemp = 0.0;
                                        for (i45 = 0; i45 < j->size[1]; i45++) {
                                            dtemp += j->data[j->size[0] * i45] * b->data[i45];
                                        }
                                    } else {
                                        dtemp = 0.0;
                                        for (i45 = 0; i45 < j->size[1]; i45++) {
                                            dtemp += j->data[j->size[0] * i45] * b->data[i45];
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
                                        exitg4 = false;
                                        while ((!exitg4) && (b_l > knz)) {
                                            /*  gee */
                                            err = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                                            i45 = e_x->size[0] * e_x->size[1];
                                            e_x->size[0] = 1;
                                            e_x->size[1] = x->size[1];
                                            emxEnsureCapacity((emxArray__common *)e_x, i45, sizeof
                                                              (double));
                                            loop_ub = x->size[0] * x->size[1];
                                            for (i45 = 0; i45 < loop_ub; i45++) {
                                                e_x->data[i45] = err - x->data[i45];
                                            }

                                            c_rdivide(ad, e_x, j);
                                            i45 = b->size[0];
                                            b->size[0] = y->size[1];
                                            emxEnsureCapacity((emxArray__common *)b, i45, sizeof
                                                              (double));
                                            loop_ub = y->size[1];
                                            for (i45 = 0; i45 < loop_ub; i45++) {
                                                b->data[i45] = y->data[y->size[0] * i45];
                                            }

                                            if ((j->size[1] == 1) || (b->size[0] == 1)) {
                                                dtemp = 0.0;
                                                for (i45 = 0; i45 < j->size[1]; i45++) {
                                                    dtemp += j->data[j->size[0] * i45] * b->data[i45];
                                                }
                                            } else {
                                                dtemp = 0.0;
                                                for (i45 = 0; i45 < j->size[1]; i45++) {
                                                    dtemp += j->data[j->size[0] * i45] * b->data[i45];
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
                                                exitg4 = true;
                                            }
                                        }

                                        iext->data[(int)((nfcns + 1.0) + 1.0) - 1] = b_l + 1.0;
                                        d_j = ((nfcns + 1.0) + 1.0) + 1.0;
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
                                            i45 = -2;
                                            i46 = 0;
                                        } else {
                                            i45 = -1;
                                            i46 = (int)temp;
                                        }

                                        temp = (nfcns + 1.0) - nfcns;
                                        if (temp > (nfcns + 1.0) - 1.0) {
                                            i47 = 1;
                                            i48 = 0;
                                        } else {
                                            i47 = (int)temp;
                                            i48 = (int)((nfcns + 1.0) - 1.0);
                                        }

                                        /*  Update index */
                                        temp = (nfcns + 1.0) - nfcns;
                                        if (2.0 > temp) {
                                            i49 = -2;
                                            i50 = 0;
                                        } else {
                                            i49 = -1;
                                            i50 = (int)temp;
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
                                        mtmp->size[1] = ((i46 - i45) + i48) - i47;
                                        emxEnsureCapacity((emxArray__common *)mtmp, nu, sizeof
                                                          (double));
                                        mtmp->data[0] = c_k1;
                                        loop_ub = i46 - i45;
                                        for (nu = 0; nu <= loop_ub - 3; nu++) {
                                            mtmp->data[mtmp->size[0] * (nu + 1)] = iext->data[(i45 +
                                                                                   nu) + 2];
                                        }

                                        loop_ub = i48 - i47;
                                        for (i48 = 0; i48 <= loop_ub; i48++) {
                                            mtmp->data[mtmp->size[0] * (((i48 + i46) - i45) - 1)] =
                                                iext->data[(i47 + i48) - 1];
                                        }

                                        loop_ub = mtmp->size[1];
                                        i45 = r27->size[0] * r27->size[1];
                                        r27->size[0] = 1;
                                        r27->size[1] = loop_ub;
                                        emxEnsureCapacity((emxArray__common *)r27, i45, sizeof(int));
                                        for (i45 = 0; i45 < loop_ub; i45++) {
                                            r27->data[r27->size[0] * i45] = i45;
                                        }

                                        i45 = b_mtmp->size[0];
                                        b_mtmp->size[0] = ((i50 - i49) + flag) - ixstart;
                                        emxEnsureCapacity((emxArray__common *)b_mtmp, i45, sizeof
                                                          (double));
                                        b_mtmp->data[0] = c_k1;
                                        loop_ub = i50 - i49;
                                        for (i45 = 0; i45 <= loop_ub - 3; i45++) {
                                            b_mtmp->data[i45 + 1] = iext->data[(i49 + i45) + 2];
                                        }

                                        loop_ub = flag - ixstart;
                                        for (i45 = 0; i45 <= loop_ub; i45++) {
                                            b_mtmp->data[((i45 + i50) - i49) - 1] = iext->data
                                                                                    [(ixstart + i45) - 1];
                                        }

                                        loop_ub = r27->size[1];
                                        for (i45 = 0; i45 < loop_ub; i45++) {
                                            iext->data[r27->data[r27->size[0] * i45]] = b_mtmp->data[(*
                                                    (int (*)[2])r27->size)[0] * i45];
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

                    if ((flag34 != 0) && (d_j > (nfcns + 1.0) + 1.0)) {
                        if (luck > 9) {
                            if (2.0 > nfcns + 1.0) {
                                i45 = 0;
                                i46 = 0;
                            } else {
                                i45 = 1;
                                i46 = (int)(nfcns + 1.0);
                            }

                            if (nfcns + 1.0 > (nfcns + 1.0) - 1.0) {
                                i47 = 0;
                                i48 = 0;
                            } else {
                                i47 = (int)(nfcns + 1.0) - 1;
                                i48 = (int)((nfcns + 1.0) - 1.0);
                            }

                            /*  Update index */
                            if (2.0 > nfcns + 1.0) {
                                i49 = 0;
                                i50 = 0;
                            } else {
                                i49 = 1;
                                i50 = (int)(nfcns + 1.0);
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
                            e_iext->size[1] = (((i46 - i45) + i48) - i47) + 2;
                            emxEnsureCapacity((emxArray__common *)e_iext, nu, sizeof(double));
                            loop_ub = i46 - i45;
                            for (nu = 0; nu < loop_ub; nu++) {
                                e_iext->data[e_iext->size[0] * nu] = iext->data[i45 + nu];
                            }

                            loop_ub = i48 - i47;
                            for (nu = 0; nu < loop_ub; nu++) {
                                e_iext->data[e_iext->size[0] * ((nu + i46) - i45)] = iext->
                                        data[i47 + nu];
                            }

                            e_iext->data[e_iext->size[0] * (((i46 - i45) + i48) - i47)] =
                                iext->data[(int)((nfcns + 1.0) + 1.0) - 1];
                            e_iext->data[e_iext->size[0] * ((((i46 - i45) + i48) - i47) + 1)] =
                                iext->data[(int)((nfcns + 1.0) + 1.0) - 1];
                            loop_ub = e_iext->size[1];
                            i45 = r27->size[0] * r27->size[1];
                            r27->size[0] = 1;
                            r27->size[1] = loop_ub;
                            emxEnsureCapacity((emxArray__common *)r27, i45, sizeof(int));
                            for (i45 = 0; i45 < loop_ub; i45++) {
                                r27->data[r27->size[0] * i45] = i45;
                            }

                            temp = iext->data[(int)((nfcns + 1.0) + 1.0) - 1];
                            dnum = iext->data[(int)((nfcns + 1.0) + 1.0) - 1];
                            i45 = f_iext->size[0];
                            f_iext->size[0] = (((i50 - i49) + flag) - ixstart) + 2;
                            emxEnsureCapacity((emxArray__common *)f_iext, i45, sizeof(double));
                            loop_ub = i50 - i49;
                            for (i45 = 0; i45 < loop_ub; i45++) {
                                f_iext->data[i45] = iext->data[i49 + i45];
                            }

                            loop_ub = flag - ixstart;
                            for (i45 = 0; i45 < loop_ub; i45++) {
                                f_iext->data[(i45 + i50) - i49] = iext->data[ixstart + i45];
                            }

                            f_iext->data[((i50 - i49) + flag) - ixstart] = temp;
                            f_iext->data[(((i50 - i49) + flag) - ixstart) + 1] = dnum;
                            loop_ub = r27->size[1];
                            for (i45 = 0; i45 < loop_ub; i45++) {
                                iext->data[r27->data[r27->size[0] * i45]] = f_iext->data[(*(int
                                        (*)[2])r27->size)[0] * i45];
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

                            c_k1 = iext->data[(int)((nfcns + 1.0) + 1.0) - 1];
                            comp = b_y1 * 1.00001;
                            b_l = ((double)ngrid + 1.0) - 1.0;
                            exitg2 = false;
                            while ((!exitg2) && (b_l > knz)) {
                                /*  gee */
                                err = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                                i45 = d_x->size[0] * d_x->size[1];
                                d_x->size[0] = 1;
                                d_x->size[1] = x->size[1];
                                emxEnsureCapacity((emxArray__common *)d_x, i45, sizeof(double));
                                loop_ub = x->size[0] * x->size[1];
                                for (i45 = 0; i45 < loop_ub; i45++) {
                                    d_x->data[i45] = err - x->data[i45];
                                }

                                c_rdivide(ad, d_x, j);
                                i45 = b->size[0];
                                b->size[0] = y->size[1];
                                emxEnsureCapacity((emxArray__common *)b, i45, sizeof(double));
                                loop_ub = y->size[1];
                                for (i45 = 0; i45 < loop_ub; i45++) {
                                    b->data[i45] = y->data[y->size[0] * i45];
                                }

                                if ((j->size[1] == 1) || (b->size[0] == 1)) {
                                    dtemp = 0.0;
                                    for (i45 = 0; i45 < j->size[1]; i45++) {
                                        dtemp += j->data[j->size[0] * i45] * b->data[i45];
                                    }
                                } else {
                                    dtemp = 0.0;
                                    for (i45 = 0; i45 < j->size[1]; i45++) {
                                        dtemp += j->data[j->size[0] * i45] * b->data[i45];
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
                                    exitg4 = false;
                                    while ((!exitg4) && (b_l > knz)) {
                                        /*  gee */
                                        err = cos(6.2831853071795862 * grid->data[(int)b_l - 1]);
                                        i45 = c_x->size[0] * c_x->size[1];
                                        c_x->size[0] = 1;
                                        c_x->size[1] = x->size[1];
                                        emxEnsureCapacity((emxArray__common *)c_x, i45, sizeof
                                                          (double));
                                        loop_ub = x->size[0] * x->size[1];
                                        for (i45 = 0; i45 < loop_ub; i45++) {
                                            c_x->data[i45] = err - x->data[i45];
                                        }

                                        c_rdivide(ad, c_x, j);
                                        i45 = b->size[0];
                                        b->size[0] = y->size[1];
                                        emxEnsureCapacity((emxArray__common *)b, i45, sizeof(double));
                                        loop_ub = y->size[1];
                                        for (i45 = 0; i45 < loop_ub; i45++) {
                                            b->data[i45] = y->data[y->size[0] * i45];
                                        }

                                        if ((j->size[1] == 1) || (b->size[0] == 1)) {
                                            dtemp = 0.0;
                                            for (i45 = 0; i45 < j->size[1]; i45++) {
                                                dtemp += j->data[j->size[0] * i45] * b->data[i45];
                                            }
                                        } else {
                                            dtemp = 0.0;
                                            for (i45 = 0; i45 < j->size[1]; i45++) {
                                                dtemp += j->data[j->size[0] * i45] * b->data[i45];
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
                                            exitg4 = true;
                                        }
                                    }

                                    iext->data[(int)((nfcns + 1.0) + 1.0) - 1] = b_l + 1.0;
                                    jchnge++;
                                    if (2.0 > nfcns + 1.0) {
                                        i45 = 0;
                                        i46 = 0;
                                    } else {
                                        i45 = 1;
                                        i46 = (int)(nfcns + 1.0);
                                    }

                                    if (nfcns + 1.0 > (nfcns + 1.0) - 1.0) {
                                        i47 = 0;
                                        i48 = 0;
                                    } else {
                                        i47 = (int)(nfcns + 1.0) - 1;
                                        i48 = (int)((nfcns + 1.0) - 1.0);
                                    }

                                    /*  Update index */
                                    if (2.0 > nfcns + 1.0) {
                                        i49 = 0;
                                        i50 = 0;
                                    } else {
                                        i49 = 1;
                                        i50 = (int)(nfcns + 1.0);
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
                                    c_iext->size[1] = (((i46 - i45) + i48) - i47) + 2;
                                    emxEnsureCapacity((emxArray__common *)c_iext, nu, sizeof
                                                      (double));
                                    loop_ub = i46 - i45;
                                    for (nu = 0; nu < loop_ub; nu++) {
                                        c_iext->data[c_iext->size[0] * nu] = iext->data[i45 + nu];
                                    }

                                    loop_ub = i48 - i47;
                                    for (nu = 0; nu < loop_ub; nu++) {
                                        c_iext->data[c_iext->size[0] * ((nu + i46) - i45)] =
                                            iext->data[i47 + nu];
                                    }

                                    c_iext->data[c_iext->size[0] * (((i46 - i45) + i48) - i47)] =
                                        iext->data[(int)((nfcns + 1.0) + 1.0) - 1];
                                    c_iext->data[c_iext->size[0] * ((((i46 - i45) + i48) - i47) +
                                                                    1)] = iext->data[(int)((nfcns + 1.0) + 1.0) - 1];
                                    loop_ub = c_iext->size[1];
                                    i45 = r27->size[0] * r27->size[1];
                                    r27->size[0] = 1;
                                    r27->size[1] = loop_ub;
                                    emxEnsureCapacity((emxArray__common *)r27, i45, sizeof(int));
                                    for (i45 = 0; i45 < loop_ub; i45++) {
                                        r27->data[r27->size[0] * i45] = i45;
                                    }

                                    temp = iext->data[(int)((nfcns + 1.0) + 1.0) - 1];
                                    dnum = iext->data[(int)((nfcns + 1.0) + 1.0) - 1];
                                    i45 = d_iext->size[0];
                                    d_iext->size[0] = (((i50 - i49) + flag) - ixstart) + 2;
                                    emxEnsureCapacity((emxArray__common *)d_iext, i45, sizeof
                                                      (double));
                                    loop_ub = i50 - i49;
                                    for (i45 = 0; i45 < loop_ub; i45++) {
                                        d_iext->data[i45] = iext->data[i49 + i45];
                                    }

                                    loop_ub = flag - ixstart;
                                    for (i45 = 0; i45 < loop_ub; i45++) {
                                        d_iext->data[(i45 + i50) - i49] = iext->data[ixstart + i45];
                                    }

                                    d_iext->data[((i50 - i49) + flag) - ixstart] = temp;
                                    d_iext->data[(((i50 - i49) + flag) - ixstart) + 1] = dnum;
                                    loop_ub = r27->size[1];
                                    for (i45 = 0; i45 < loop_ub; i45++) {
                                        iext->data[r27->data[r27->size[0] * i45]] = d_iext->data
                                                [(*(int (*)[2])r27->size)[0] * i45];
                                    }

                                    exitg2 = true;
                                } else {
                                    b_l--;
                                }
                            }

                            if (luck != 6) {
                                temp = (nfcns + 1.0) - nfcns;
                                if (2.0 > temp) {
                                    i45 = -2;
                                    i46 = 0;
                                } else {
                                    i45 = -1;
                                    i46 = (int)temp;
                                }

                                temp = (nfcns + 1.0) - nfcns;
                                if (temp > (nfcns + 1.0) - 1.0) {
                                    i47 = 1;
                                    i48 = 0;
                                } else {
                                    i47 = (int)temp;
                                    i48 = (int)((nfcns + 1.0) - 1.0);
                                }

                                /*  Update index */
                                temp = (nfcns + 1.0) - nfcns;
                                if (2.0 > temp) {
                                    i49 = -2;
                                    i50 = 0;
                                } else {
                                    i49 = -1;
                                    i50 = (int)temp;
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
                                k1->size[1] = ((i46 - i45) + i48) - i47;
                                emxEnsureCapacity((emxArray__common *)k1, nu, sizeof(double));
                                k1->data[0] = c_k1;
                                loop_ub = i46 - i45;
                                for (nu = 0; nu <= loop_ub - 3; nu++) {
                                    k1->data[k1->size[0] * (nu + 1)] = iext->data[(i45 + nu) + 2];
                                }

                                loop_ub = i48 - i47;
                                for (i48 = 0; i48 <= loop_ub; i48++) {
                                    k1->data[k1->size[0] * (((i48 + i46) - i45) - 1)] = iext->
                                            data[(i47 + i48) - 1];
                                }

                                loop_ub = k1->size[1];
                                i45 = r27->size[0] * r27->size[1];
                                r27->size[0] = 1;
                                r27->size[1] = loop_ub;
                                emxEnsureCapacity((emxArray__common *)r27, i45, sizeof(int));
                                for (i45 = 0; i45 < loop_ub; i45++) {
                                    r27->data[r27->size[0] * i45] = i45;
                                }

                                i45 = b_k1->size[0];
                                b_k1->size[0] = ((i50 - i49) + flag) - ixstart;
                                emxEnsureCapacity((emxArray__common *)b_k1, i45, sizeof(double));
                                b_k1->data[0] = c_k1;
                                loop_ub = i50 - i49;
                                for (i45 = 0; i45 <= loop_ub - 3; i45++) {
                                    b_k1->data[i45 + 1] = iext->data[(i49 + i45) + 2];
                                }

                                loop_ub = flag - ixstart;
                                for (i45 = 0; i45 <= loop_ub; i45++) {
                                    b_k1->data[((i45 + i50) - i49) - 1] = iext->data[(ixstart +
                                                                          i45) - 1];
                                }

                                loop_ub = r27->size[1];
                                for (i45 = 0; i45 < loop_ub; i45++) {
                                    iext->data[r27->data[r27->size[0] * i45]] = b_k1->data[(*(int
                                            (*)[2])r27->size)[0] * i45];
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
        d_j = -1.0;
        c_k1 = -1.0;

        /*  initialize memory */
        /* x(nzz) = -2; */
        i45 = x2->size[0] * x2->size[1];
        x2->size[0] = 1;
        x2->size[1] = x->size[1] + 1;
        emxEnsureCapacity((emxArray__common *)x2, i45, sizeof(double));
        loop_ub = x->size[1];
        for (i45 = 0; i45 < loop_ub; i45++) {
            x2->data[x2->size[0] * i45] = x->data[x->size[0] * i45];
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
            d_j = -(dtemp + dnum) / (dtemp - dnum);
        }

        i45 = a->size[0] * a->size[1];
        a->size[0] = 1;
        a->size[1] = (int)nfcns;
        emxEnsureCapacity((emxArray__common *)a, i45, sizeof(double));
        for (nut = 0; nut < (int)nfcns; nut++) {
            temp = ((1.0 + (double)nut) - 1.0) * knz;
            dnum = cos(6.2831853071795862 * temp);
            if (nu != 1) {
                dnum = (dnum - d_j) / c_k1;
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
                i45 = b_x->size[0] * b_x->size[1];
                b_x->size[0] = 1;
                b_x->size[1] = (int)(nfcns + 1.0);
                emxEnsureCapacity((emxArray__common *)b_x, i45, sizeof(double));
                for (i45 = 0; i45 < loop_ub; i45++) {
                    b_x->data[b_x->size[0] * i45] = err - x2->data[i45];
                }

                c_rdivide(ad, b_x, j);
                i45 = b->size[0];
                b->size[0] = y->size[1];
                emxEnsureCapacity((emxArray__common *)b, i45, sizeof(double));
                loop_ub = y->size[1];
                for (i45 = 0; i45 < loop_ub; i45++) {
                    b->data[i45] = y->data[y->size[0] * i45];
                }

                if ((j->size[1] == 1) || (b->size[0] == 1)) {
                    dtemp = 0.0;
                    for (i45 = 0; i45 < j->size[1]; i45++) {
                        dtemp += j->data[j->size[0] * i45] * b->data[i45];
                    }
                } else {
                    dtemp = 0.0;
                    for (i45 = 0; i45 < j->size[1]; i45++) {
                        dtemp += j->data[j->size[0] * i45] * b->data[i45];
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
        i45 = j->size[0] * j->size[1];
        j->size[0] = 1;
        j->size[1] = (int)nfcns;
        emxEnsureCapacity((emxArray__common *)j, i45, sizeof(double));
        for (nut = 0; nut < (int)nfcns; nut++) {
            dnum = ((1.0 + (double)nut) - 1.0) * temp;
            if (nfcns - 1.0 < 1.0) {
                j->data[nut] = a->data[0];
            } else {
                if (2.0 > nfcns) {
                    i45 = 1;
                    i46 = 1;
                } else {
                    i45 = 2;
                    i46 = (int)nfcns + 1;
                }

                if (nfcns - 1.0 < 1.0) {
                    i47 = y->size[0] * y->size[1];
                    y->size[0] = 1;
                    y->size[1] = 0;
                    emxEnsureCapacity((emxArray__common *)y, i47, sizeof(double));
                } else {
                    i47 = y->size[0] * y->size[1];
                    y->size[0] = 1;
                    y->size[1] = (int)floor((nfcns - 1.0) - 1.0) + 1;
                    emxEnsureCapacity((emxArray__common *)y, i47, sizeof(double));
                    loop_ub = (int)floor((nfcns - 1.0) - 1.0);
                    for (i47 = 0; i47 <= loop_ub; i47++) {
                        y->data[y->size[0] * i47] = 1.0 + (double)i47;
                    }
                }

                i47 = y->size[0] * y->size[1];
                y->size[0] = 1;
                emxEnsureCapacity((emxArray__common *)y, i47, sizeof(double));
                ixstart = y->size[0];
                flag = y->size[1];
                loop_ub = ixstart * flag;
                for (i47 = 0; i47 < loop_ub; i47++) {
                    y->data[i47] *= dnum;
                }

                b_cos(y);
                i47 = y->size[0] * y->size[1];
                y->size[0] = 1;
                emxEnsureCapacity((emxArray__common *)y, i47, sizeof(double));
                ixstart = y->size[0];
                flag = y->size[1];
                loop_ub = ixstart * flag;
                for (i47 = 0; i47 < loop_ub; i47++) {
                    y->data[i47] *= 2.0;
                }

                i47 = b->size[0];
                b->size[0] = i46 - i45;
                emxEnsureCapacity((emxArray__common *)b, i47, sizeof(double));
                loop_ub = i46 - i45;
                for (i47 = 0; i47 < loop_ub; i47++) {
                    b->data[i47] = a->data[(i45 + i47) - 1];
                }

                if ((y->size[1] == 1) || (i46 - i45 == 1)) {
                    dtemp = 0.0;
                    for (i45 = 0; i45 < y->size[1]; i45++) {
                        dtemp += y->data[y->size[0] * i45] * b->data[i45];
                    }
                } else {
                    dtemp = 0.0;
                    for (i45 = 0; i45 < y->size[1]; i45++) {
                        dtemp += y->data[y->size[0] * i45] * b->data[i45];
                    }
                }

                j->data[nut] = a->data[0] + dtemp;
            }
        }

        if (2.0 > nfcns) {
            i45 = -1;
            i46 = 0;
        } else {
            i45 = 0;
            i46 = (int)nfcns;
        }

        i47 = iext->size[0];
        iext->size[0] = i46 - i45;
        emxEnsureCapacity((emxArray__common *)iext, i47, sizeof(double));
        iext->data[0] = j->data[0] / jchnge;
        loop_ub = i46 - i45;
        for (i46 = 0; i46 <= loop_ub - 2; i46++) {
            iext->data[i46 + 1] = 2.0 * j->data[(i45 + i46) + 1] / jchnge;
        }

        i45 = j->size[0] * j->size[1];
        j->size[0] = 1;
        j->size[1] = (int)nfcns;
        emxEnsureCapacity((emxArray__common *)j, i45, sizeof(double));
        loop_ub = (int)nfcns;
        for (i45 = 0; i45 < loop_ub; i45++) {
            j->data[i45] = 0.0;
        }

        i45 = l->size[0] * l->size[1];
        l->size[0] = 1;
        l->size[1] = (int)nfcns - 1;
        emxEnsureCapacity((emxArray__common *)l, i45, sizeof(double));
        loop_ub = (int)nfcns - 1;
        for (i45 = 0; i45 < loop_ub; i45++) {
            l->data[i45] = 0.0;
        }

        if (nu != 1) {
            j->data[0] = 2.0 * iext->data[(int)nfcns - 1] * d_j + iext->data[(int)
                         nfcns - 2];
            j->data[1] = 2.0 * c_k1 * iext->data[(int)nfcns - 1];
            l->data[0] = iext->data[(int)nfcns - 3] - iext->data[(int)nfcns - 1];
            for (nut = 0; nut <= (int)nfcns - 3; nut++) {
                if (2 + nut == (int)nfcns - 1) {
                    c_k1 /= 2.0;
                    d_j /= 2.0;
                }

                j->data[nut + 2] = 0.0;
                i45 = x2->size[0] * x2->size[1];
                x2->size[0] = 1;
                x2->size[1] = (int)((2.0 + (double)nut) - 1.0) + 1;
                emxEnsureCapacity((emxArray__common *)x2, i45, sizeof(double));
                loop_ub = (int)((2.0 + (double)nut) - 1.0);
                for (i45 = 0; i45 <= loop_ub; i45++) {
                    x2->data[x2->size[0] * i45] = 1.0 + (double)i45;
                }

                loop_ub = x2->size[0] * x2->size[1];
                for (i45 = 0; i45 < loop_ub; i45++) {
                    a->data[(int)x2->data[i45] - 1] = j->data[(int)x2->data[i45] - 1];
                }

                temp = 2.0 * d_j;
                loop_ub = x2->size[0] * x2->size[1];
                for (i45 = 0; i45 < loop_ub; i45++) {
                    j->data[(int)x2->data[i45] - 1] = temp * a->data[(int)x2->data[i45] -
                                                      1];
                }

                j->data[1] += 2.0 * a->data[0] * c_k1;
                i45 = x2->size[0] * x2->size[1];
                x2->size[0] = 1;
                x2->size[1] = (int)(((2.0 + (double)nut) - 1.0) - 1.0) + 1;
                emxEnsureCapacity((emxArray__common *)x2, i45, sizeof(double));
                loop_ub = (int)(((2.0 + (double)nut) - 1.0) - 1.0);
                for (i45 = 0; i45 <= loop_ub; i45++) {
                    x2->data[x2->size[0] * i45] = 1.0 + (double)i45;
                }

                i45 = x->size[0] * x->size[1];
                x->size[0] = 1;
                x->size[1] = x2->size[1];
                emxEnsureCapacity((emxArray__common *)x, i45, sizeof(double));
                loop_ub = x2->size[0] * x2->size[1];
                for (i45 = 0; i45 < loop_ub; i45++) {
                    x->data[i45] = a->data[(int)x2->data[i45]];
                }

                i45 = b_j->size[0];
                b_j->size[0] = x2->size[0] * x2->size[1];
                emxEnsureCapacity((emxArray__common *)b_j, i45, sizeof(double));
                loop_ub = x2->size[0] * x2->size[1];
                for (i45 = 0; i45 < loop_ub; i45++) {
                    b_j->data[i45] = (j->data[(int)x2->data[i45] - 1] + l->data[(int)
                                      x2->data[i45] - 1]) + c_k1 * x->data[i45];
                }

                loop_ub = b_j->size[0];
                for (i45 = 0; i45 < loop_ub; i45++) {
                    j->data[(int)x2->data[i45] - 1] = b_j->data[i45];
                }

                i45 = x2->size[0] * x2->size[1];
                x2->size[0] = 1;
                x2->size[1] = (int)(((2.0 + (double)nut) + 1.0) - 3.0) + 1;
                emxEnsureCapacity((emxArray__common *)x2, i45, sizeof(double));
                loop_ub = (int)(((2.0 + (double)nut) + 1.0) - 3.0);
                for (i45 = 0; i45 <= loop_ub; i45++) {
                    x2->data[x2->size[0] * i45] = 3.0 + (double)i45;
                }

                i45 = x->size[0] * x->size[1];
                x->size[0] = 1;
                x->size[1] = x2->size[1];
                emxEnsureCapacity((emxArray__common *)x, i45, sizeof(double));
                loop_ub = x2->size[0] * x2->size[1];
                for (i45 = 0; i45 < loop_ub; i45++) {
                    x->data[i45] = a->data[(int)x2->data[i45] - 2];
                }

                i45 = c_j->size[0];
                c_j->size[0] = x2->size[0] * x2->size[1];
                emxEnsureCapacity((emxArray__common *)c_j, i45, sizeof(double));
                loop_ub = x2->size[0] * x2->size[1];
                for (i45 = 0; i45 < loop_ub; i45++) {
                    c_j->data[i45] = j->data[(int)x2->data[i45] - 1] + c_k1 * x->data[i45];
                }

                loop_ub = c_j->size[0];
                for (i45 = 0; i45 < loop_ub; i45++) {
                    j->data[(int)x2->data[i45] - 1] = c_j->data[i45];
                }

                if (2 + nut != (int)nfcns - 1) {
                    i45 = x2->size[0] * x2->size[1];
                    x2->size[0] = 1;
                    x2->size[1] = (int)((2.0 + (double)nut) - 1.0) + 1;
                    emxEnsureCapacity((emxArray__common *)x2, i45, sizeof(double));
                    loop_ub = (int)((2.0 + (double)nut) - 1.0);
                    for (i45 = 0; i45 <= loop_ub; i45++) {
                        x2->data[x2->size[0] * i45] = 1.0 + (double)i45;
                    }

                    loop_ub = x2->size[0] * x2->size[1];
                    for (i45 = 0; i45 < loop_ub; i45++) {
                        l->data[(int)x2->data[i45] - 1] = -a->data[(int)x2->data[i45] - 1];
                    }

                    l->data[0] += iext->data[((int)nfcns - nut) - 4];
                }
            }

            loop_ub = (int)nfcns;
            for (i45 = 0; i45 < loop_ub; i45++) {
                iext->data[i45] = j->data[i45];
            }
        }

        /*  alpha must be at lease >=3 */
        if (nfcns <= 3.0) {
            /* alpha(nfcns + 1) = 0; */
            /* alpha(nfcns + 2) = 0; */
            i45 = j->size[0] * j->size[1];
            j->size[0] = 1;
            j->size[1] = iext->size[0] + 2;
            emxEnsureCapacity((emxArray__common *)j, i45, sizeof(double));
            loop_ub = iext->size[0];
            for (i45 = 0; i45 < loop_ub; i45++) {
                j->data[j->size[0] * i45] = iext->data[i45];
            }

            j->data[j->size[0] * iext->size[0]] = 0.0;
            j->data[j->size[0] * (iext->size[0] + 1)] = 0.0;
        } else {
            i45 = j->size[0] * j->size[1];
            j->size[0] = 1;
            j->size[1] = iext->size[0];
            emxEnsureCapacity((emxArray__common *)j, i45, sizeof(double));
            loop_ub = iext->size[0];
            for (i45 = 0; i45 < loop_ub; i45++) {
                j->data[j->size[0] * i45] = iext->data[i45];
            }
        }

        /* alpha=alpha'; */
        /*  now that's done! */
        if (nodd != 0.0) {
            i45 = x->size[0] * x->size[1];
            x->size[0] = 1;
            x->size[1] = (int)floor(-((0.0 - (nfcns - 1.0)) - -1.0)) + 1;
            emxEnsureCapacity((emxArray__common *)x, i45, sizeof(double));
            loop_ub = (int)floor(-((0.0 - (nfcns - 1.0)) - -1.0));
            for (i45 = 0; i45 <= loop_ub; i45++) {
                x->data[x->size[0] * i45] = j->data[(int)((nfcns + 1.0) + (-1.0 -
                                                    (double)i45)) - 1];
            }

            i45 = h->size[0] * h->size[1];
            h->size[0] = 1;
            h->size[1] = x->size[1] + 1;
            emxEnsureCapacity((emxArray__common *)h, i45, sizeof(double));
            loop_ub = x->size[1];
            for (i45 = 0; i45 < loop_ub; i45++) {
                h->data[h->size[0] * i45] = 0.5 * x->data[x->size[0] * i45];
            }

            h->data[h->size[0] * x->size[1]] = j->data[0];
        } else {
            if ((nfcns - (nfcns - 1.0)) + 2.0 > nfcns) {
                i45 = 0;
                i46 = 1;
            } else {
                i45 = (int)nfcns - 1;
                i46 = -1;
            }

            i47 = x2->size[0] * x2->size[1];
            x2->size[0] = 1;
            x2->size[1] = (int)floor(-((0.0 - (nfcns - 1.0)) - -2.0)) + 1;
            emxEnsureCapacity((emxArray__common *)x2, i47, sizeof(double));
            loop_ub = (int)floor(-((0.0 - (nfcns - 1.0)) - -2.0));
            for (i47 = 0; i47 <= loop_ub; i47++) {
                x2->data[x2->size[0] * i47] = j->data[(int)((nfcns + 1.0) + (double)(int)
                                                      (-2.0 - (double)i47)) - 1];
            }

            i47 = h->size[0] * h->size[1];
            h->size[0] = 1;
            h->size[1] = 2 + x2->size[1];
            emxEnsureCapacity((emxArray__common *)h, i47, sizeof(double));
            h->data[0] = 0.25 * j->data[(int)nfcns - 1];
            loop_ub = x2->size[1];
            for (i47 = 0; i47 < loop_ub; i47++) {
                h->data[h->size[0] * (i47 + 1)] = 0.25 * (x2->data[x2->size[0] * i47] +
                                                  j->data[i45 + i46 * i47]);
            }

            h->data[h->size[0] * (1 + x2->size[1])] = 0.25 * (2.0 * j->data[0] +
                    j->data[1]);
        }
    }

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
    emxFree_int8_T(&r28);
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
    emxFree_int32_T(&r27);
    emxFree_real_T(&a);
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
        emxEnsureCapacity((emxArray__common *)aR, ii, sizeof(creal_T));
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
        emxEnsureCapacity((emxArray__common *)aR, ii, sizeof(creal_T));
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
    emxArray_real_T *r17;
    int i39;
    int loop_ub;
    emxArray_real_T *digw;
    emxArray_creal_T *s;
    emxArray_creal_T *y;
    boolean_T b9;
    int k;
    double re;
    double im;
    emxInit_real_T(&r17, 2);

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
    i39 = r17->size[0] * r17->size[1];
    r17->size[0] = 1;
    r17->size[1] = w->size[1];
    emxEnsureCapacity((emxArray__common *)r17, i39, sizeof(double));
    loop_ub = w->size[0] * w->size[1];
    for (i39 = 0; i39 < loop_ub; i39++) {
        r17->data[i39] = 6.2831853071795862 * w->data[i39];
    }

    emxInit_real_T(&digw, 2);
    emxInit_creal_T(&s, 2);
    rdivide(r17, Fs, digw);

    /*  Convert from Hz to rad/sample for computational purposes */
    i39 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = digw->size[1];
    emxEnsureCapacity((emxArray__common *)s, i39, sizeof(creal_T));
    loop_ub = digw->size[0] * digw->size[1];
    emxFree_real_T(&r17);
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
    emxEnsureCapacity((emxArray__common *)y, i39, sizeof(creal_T));
    b9 = (y->size[1] == 0);
    if (!b9) {
        i39 = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)y, i39, sizeof(creal_T));
        loop_ub = y->size[1];
        for (i39 = 0; i39 < loop_ub; i39++) {
            y->data[y->size[0] * i39].re = b[0];
            y->data[y->size[0] * i39].im = 0.0;
        }

        for (k = 0; k < 18; k++) {
            i39 = y->size[0] * y->size[1];
            y->size[0] = 1;
            y->size[1] = s->size[1];
            emxEnsureCapacity((emxArray__common *)y, i39, sizeof(creal_T));
            loop_ub = s->size[0] * s->size[1];
            for (i39 = 0; i39 < loop_ub; i39++) {
                re = s->data[i39].re * y->data[i39].re - s->data[i39].im * y->data[i39].
                     im;
                im = s->data[i39].re * y->data[i39].im + s->data[i39].im * y->data[i39].
                     re;
                y->data[i39].re = re + b[k + 1];
                y->data[i39].im = im;
            }
        }
    }

    i39 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = digw->size[1];
    emxEnsureCapacity((emxArray__common *)s, i39, sizeof(creal_T));
    loop_ub = digw->size[0] * digw->size[1];
    for (i39 = 0; i39 < loop_ub; i39++) {
        re = digw->data[i39] * 0.0;
        im = digw->data[i39];
        s->data[i39].re = 18.0 * re;
        s->data[i39].im = 18.0 * im;
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
    emxArray_real_T *r18;
    int i40;
    int loop_ub;
    emxArray_real_T *digw;
    emxArray_creal_T *s;
    emxArray_creal_T *y;
    boolean_T b10;
    int k;
    double re;
    double im;
    emxInit_real_T(&r18, 2);

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
    i40 = r18->size[0] * r18->size[1];
    r18->size[0] = 1;
    r18->size[1] = w->size[1];
    emxEnsureCapacity((emxArray__common *)r18, i40, sizeof(double));
    loop_ub = w->size[0] * w->size[1];
    for (i40 = 0; i40 < loop_ub; i40++) {
        r18->data[i40] = 6.2831853071795862 * w->data[i40];
    }

    emxInit_real_T(&digw, 2);
    emxInit_creal_T(&s, 2);
    rdivide(r18, Fs, digw);

    /*  Convert from Hz to rad/sample for computational purposes */
    i40 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = digw->size[1];
    emxEnsureCapacity((emxArray__common *)s, i40, sizeof(creal_T));
    loop_ub = digw->size[0] * digw->size[1];
    emxFree_real_T(&r18);
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
    emxEnsureCapacity((emxArray__common *)y, i40, sizeof(creal_T));
    b10 = (y->size[1] == 0);
    if (!b10) {
        i40 = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)y, i40, sizeof(creal_T));
        loop_ub = y->size[1];
        for (i40 = 0; i40 < loop_ub; i40++) {
            y->data[y->size[0] * i40].re = b[0];
            y->data[y->size[0] * i40].im = 0.0;
        }

        for (k = 0; k < 84; k++) {
            i40 = y->size[0] * y->size[1];
            y->size[0] = 1;
            y->size[1] = s->size[1];
            emxEnsureCapacity((emxArray__common *)y, i40, sizeof(creal_T));
            loop_ub = s->size[0] * s->size[1];
            for (i40 = 0; i40 < loop_ub; i40++) {
                re = s->data[i40].re * y->data[i40].re - s->data[i40].im * y->data[i40].
                     im;
                im = s->data[i40].re * y->data[i40].im + s->data[i40].im * y->data[i40].
                     re;
                y->data[i40].re = re + b[k + 1];
                y->data[i40].im = im;
            }
        }
    }

    i40 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = digw->size[1];
    emxEnsureCapacity((emxArray__common *)s, i40, sizeof(creal_T));
    loop_ub = digw->size[0] * digw->size[1];
    for (i40 = 0; i40 < loop_ub; i40++) {
        re = digw->data[i40] * 0.0;
        im = digw->data[i40];
        s->data[i40].re = 84.0 * re;
        s->data[i40].im = 84.0 * im;
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
 * Arguments    : const double x_data[]
 *                const int x_size[2]
 *                double N
 *                double y_data[]
 *                int y_size[2]
 * Return Type  : void
 */
static void upsample(const double x_data[], const int x_size[2], double N,
                     double y_data[], int y_size[2])
{
    int vleny;
    int iy;
    int k;
    vleny = (int)N * x_size[1];
    y_size[0] = 1;
    y_size[1] = vleny;
    for (iy = 0; iy < vleny; iy++) {
        y_data[iy] = 0.0;
    }

    vleny = 0;
    iy = 0;
    for (k = 1; k <= x_size[1]; k++) {
        y_data[iy] = x_data[vleny];
        vleny++;
        iy += (int)N;
    }
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
    emxEnsureCapacity((emxArray__common *)c, k, sizeof(creal_T));
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
    int i60;
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
    i60 = tau->size[0];
    tau->size[0] = ntau;
    emxEnsureCapacity((emxArray__common *)tau, i60, sizeof(creal_T));
    ntau = a->size[0];
    i60 = work->size[0];
    work->size[0] = ntau;
    emxEnsureCapacity((emxArray__common *)work, i60, sizeof(creal_T));
    for (i60 = 0; i60 < ntau; i60++) {
        work->data[i60].re = 0.0;
        work->data[i60].im = 0.0;
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
                    i60 = (ntau + c) - 1;
                    do {
                        knt++;
                        for (k = ntau; k <= i60; k++) {
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
                i60 = (in + n * (lastv - 1)) + 1;
                knt = in + 1;
                while ((n > 0) && (knt <= i60)) {
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
                        i60 = lastc + ntau;
                        for (k = ntau; k + 1 <= i60; k++) {
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
                i60 = jy + n * (lastc - 1);
                knt = jy;
                while ((n > 0) && (knt <= i60)) {
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
                        i60 = lastv + ntau;
                        for (k = ntau; k + 1 <= i60; k++) {
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
    int i61;
    int k;
    double x_re;
    double x_im;
    i61 = (ix0 + n) - 1;
    for (k = ix0; k <= i61; k++) {
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
    emxEnsureCapacity((emxArray__common *)At, jcol, sizeof(creal_T));
    ii = A->size[0] * A->size[1];
    for (jcol = 0; jcol < ii; jcol++) {
        At->data[jcol] = A->data[jcol];
    }

    nzcount = 0;
    jcol = alpha1->size[0];
    alpha1->size[0] = At->size[0];
    emxEnsureCapacity((emxArray__common *)alpha1, jcol, sizeof(creal_T));
    ii = At->size[0];
    for (jcol = 0; jcol < ii; jcol++) {
        alpha1->data[jcol].re = 0.0;
        alpha1->data[jcol].im = 0.0;
    }

    jcol = beta1->size[0];
    beta1->size[0] = At->size[0];
    emxEnsureCapacity((emxArray__common *)beta1, jcol, sizeof(creal_T));
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
            emxEnsureCapacity((emxArray__common *)alpha1, jcol, sizeof(creal_T));
            ii = At->size[0];
            for (jcol = 0; jcol < ii; jcol++) {
                alpha1->data[jcol].re = rtNaN;
                alpha1->data[jcol].im = 0.0;
            }

            jcol = beta1->size[0];
            beta1->size[0] = At->size[0];
            emxEnsureCapacity((emxArray__common *)beta1, jcol, sizeof(creal_T));
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
                    emxEnsureCapacity((emxArray__common *)At, jcol, sizeof(creal_T));
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
                    emxEnsureCapacity((emxArray__common *)alpha1, jcol, sizeof(creal_T));
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
    emxEnsureCapacity((emxArray__common *)b_A, jm1, sizeof(creal_T));
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
    emxEnsureCapacity((emxArray__common *)alpha1, jm1, sizeof(creal_T));
    jp1 = A->size[0];
    for (jm1 = 0; jm1 < jp1; jm1++) {
        alpha1->data[jm1].re = 0.0;
        alpha1->data[jm1].im = 0.0;
    }

    jm1 = beta1->size[0];
    beta1->size[0] = A->size[0];
    emxEnsureCapacity((emxArray__common *)beta1, jm1, sizeof(creal_T));
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
                    emxEnsureCapacity((emxArray__common *)alpha1, jm1, sizeof(creal_T));
                    for (jm1 = 0; jm1 < jp1; jm1++) {
                        alpha1->data[jm1].re = rtNaN;
                        alpha1->data[jm1].im = 0.0;
                    }

                    jp1 = beta1->size[0];
                    jm1 = beta1->size[0];
                    beta1->size[0] = jp1;
                    emxEnsureCapacity((emxArray__common *)beta1, jm1, sizeof(creal_T));
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
    double v[2];
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
    emxArray_real_T *r7;

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
    emxEnsureCapacity((emxArray__common *)a, i6, sizeof(creal_T));
    a->data[0].re = -1.0;
    a->data[0].im = 0.0;
    i6 = b->size[0];
    b->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)b, i6, sizeof(double));
    b->data[0] = 1.0;
    i6 = c->size[0] * c->size[1];
    c->size[0] = 1;
    c->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)c, i6, sizeof(double));
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
    v[0] = 1.0;
    v[1] = 1.0 / wn;
    for (i6 = 0; i6 < 4; i6++) {
        t[i6] = 0.0;
    }

    /*  Balancing transformation */
    B[0] = -den[1];
    B[2] = -den[2];
    for (r1 = 0; r1 < 2; r1++) {
        t[r1 + (r1 << 1)] = v[r1];
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
    emxEnsureCapacity((emxArray__common *)reshapes_f2, i6, sizeof(signed char));
    for (i6 = 0; i6 < 2; i6++) {
        reshapes_f2->data[i6] = 0;
    }

    emxInit_cint8_T(&b_a, 2);
    i6 = b_a->size[0] * b_a->size[1];
    b_a->size[0] = a->size[0];
    b_a->size[1] = a->size[1];
    emxEnsureCapacity((emxArray__common *)b_a, i6, sizeof(cint8_T));
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
    emxEnsureCapacity((emxArray__common *)result, i6, sizeof(cint8_T));
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
    emxEnsureCapacity((emxArray__common *)b_b1, i6, sizeof(double));
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
    emxEnsureCapacity((emxArray__common *)varargin_2, i6, sizeof(double));
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
    emxEnsureCapacity((emxArray__common *)a, i6, sizeof(creal_T));
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
    emxEnsureCapacity((emxArray__common *)b, i6, sizeof(double));
    for (i6 = 0; i6 < 2; i6++) {
        b->data[r1 + i6] = b1[i6] * 0.0;
    }

    for (i6 = 0; i6 < 2; i6++) {
        b1[i6] = 0.0;
        for (i7 = 0; i7 < 2; i7++) {
            b1[i6] += (double)i7 * t[i7 + (i6 << 1)];
        }
    }

    emxInit_real_T(&r7, 2);
    i6 = r7->size[0] * r7->size[1];
    r7->size[0] = 1;
    r7->size[1] = c->size[1] + 2;
    emxEnsureCapacity((emxArray__common *)r7, i6, sizeof(double));
    r1 = c->size[1];
    for (i6 = 0; i6 < r1; i6++) {
        r7->data[r7->size[0] * i6] = 0.0 * c->data[c->size[0] * i6];
    }

    for (i6 = 0; i6 < 2; i6++) {
        r7->data[r7->size[0] * (i6 + c->size[1])] = b1[i6];
    }

    i6 = c->size[0] * c->size[1];
    c->size[0] = 1;
    c->size[1] = r7->size[1];
    emxEnsureCapacity((emxArray__common *)c, i6, sizeof(double));
    r1 = r7->size[1];
    for (i6 = 0; i6 < r1; i6++) {
        c->data[c->size[0] * i6] = r7->data[r7->size[0] * i6];
    }

    emxFree_real_T(&r7);

    /*  Apply gain k: */
    i6 = c->size[0] * c->size[1];
    c->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)c, i6, sizeof(double));
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
 *                short outputTaps[128]
 *                double *numOutputTaps
 *                double *filterGain
 * Return Type  : void
 */
void internal_design_filter_cg(double Rdata, double Fpass, double Fstop, double
                               caldiv, double FIR, double HB1, double PLL_mult, double Apass, double Astop,
                               double phEQ, double HB2, double HB3, const char Type[7], const char RxTx[2],
                               double RFbw, double DAC_div, double converter_rate, double PLL_rate, double
                               Fcenter, double wnom, double FIRdBmin, double int_FIR, short outputTaps[128],
                               double *numOutputTaps, double *filterGain)
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
    int i;
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
    emxArray_creal_T *c_combinedResponse;
    emxArray_creal_T *r0;
    int d_combinedResponse;
    emxArray_creal_T *rg2;
    double apnd;
    double absa;
    double im;
    emxArray_creal_T *rg;
    emxArray_real_T *sw;
    emxArray_real_T *F3;
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
    boolean_T exitg1;
    emxArray_real_T *F1;
    emxArray_real_T *F2;
    emxArray_real_T *A1;
    emxArray_real_T *A2;
    emxArray_real_T *W1;
    emxArray_real_T *W2;
    int Nmax;
    int tap_store_size[2];
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
    int b_tap_store_size[2];
    boolean_T valid;
    int c_tap_store_size[2];
    double firTapsPreScale[128];
    int d_tap_store_size[2];
    int e_tap_store_size[2];
    emxArray_real_T *r3;
    emxArray_real_T *r4;
    double b_firTapsPreScale[128];
    short i4;
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
        a1->size[1] = 2;
        emxEnsureCapacity((emxArray__common *)a1, i0, sizeof(creal_T));
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
        emxEnsureCapacity((emxArray__common *)a2, i0, sizeof(creal_T));
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
    if (!((Gpass < Gstop - 1.0) || rtIsNaN(Gstop - 1.0))) {
        Gpass = Gstop - 1.0;
    }

    emxInit_real_T(&fg, 2);
    i0 = fg->size[0] * fg->size[1];
    fg->size[0] = 1;
    fg->size[1] = (int)(Gpass + 1.0);
    emxEnsureCapacity((emxArray__common *)fg, i0, sizeof(double));
    loop_ub = (int)(Gpass + 1.0);
    for (i0 = 0; i0 < loop_ub; i0++) {
        fg->data[i0] = 0.0;
    }

    emxInit_real_T(&omega, 2);
    i0 = omega->size[0] * omega->size[1];
    omega->size[0] = 1;
    omega->size[1] = (int)(Gpass + 1.0);
    emxEnsureCapacity((emxArray__common *)omega, i0, sizeof(double));
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
    emxEnsureCapacity((emxArray__common *)c_combinedResponse, i0, sizeof(creal_T));
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
    emxEnsureCapacity((emxArray__common *)rg2, i0, sizeof(creal_T));
    loop_ub = omega->size[0] * omega->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
        absa = omega->data[i0] * -0.0;
        im = omega->data[i0] * -6.2831853071795862;
        rg2->data[i0].re = sigma * absa;
        rg2->data[i0].im = sigma * im;
    }

    emxInit_creal_T(&rg, 2);
    emxInit_real_T(&sw, 2);
    c_exp(rg2);
    b_rdivide(rg2, c_combinedResponse, rg);
    b_abs(c_combinedResponse, sw);
    sigma = Gpass + 1.0;

    /*  Expand memory correctly */
    emxInit_real_T(&F3, 2);
    if (rtIsNaN(Gstop)) {
        i0 = F3->size[0] * F3->size[1];
        F3->size[0] = 1;
        F3->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)F3, i0, sizeof(double));
        F3->data[0] = rtNaN;
    } else if (8192.0 < Gstop) {
        i0 = F3->size[0] * F3->size[1];
        F3->size[0] = 1;
        F3->size[1] = 0;
        emxEnsureCapacity((emxArray__common *)F3, i0, sizeof(double));
    } else if (rtIsInf(Gstop) && (Gstop == 8192.0)) {
        i0 = F3->size[0] * F3->size[1];
        F3->size[0] = 1;
        F3->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)F3, i0, sizeof(double));
        F3->data[0] = rtNaN;
    } else if (Gstop == Gstop) {
        i0 = F3->size[0] * F3->size[1];
        F3->size[0] = 1;
        F3->size[1] = (int)(8192.0 - Gstop) + 1;
        emxEnsureCapacity((emxArray__common *)F3, i0, sizeof(double));
        loop_ub = (int)(8192.0 - Gstop);
        for (i0 = 0; i0 <= loop_ub; i0++) {
            F3->data[F3->size[0] * i0] = Gstop + (double)i0;
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

        i0 = F3->size[0] * F3->size[1];
        F3->size[0] = 1;
        F3->size[1] = n;
        emxEnsureCapacity((emxArray__common *)F3, i0, sizeof(double));
        if (n > 0) {
            F3->data[0] = Gstop;
            if (n > 1) {
                F3->data[n - 1] = apnd;
                dec_int3 = (n - 1) / 2;
                for (d_combinedResponse = 1; d_combinedResponse < dec_int3;
                     d_combinedResponse++) {
                    F3->data[d_combinedResponse] = Gstop + (double)d_combinedResponse;
                    F3->data[(n - d_combinedResponse) - 1] = apnd - (double)
                            d_combinedResponse;
                }

                if (dec_int3 << 1 == n - 1) {
                    F3->data[dec_int3] = (Gstop + apnd) / 2.0;
                } else {
                    F3->data[dec_int3] = Gstop + (double)dec_int3;
                    F3->data[dec_int3 + 1] = apnd - (double)dec_int3;
                }
            }
        }
    }

    emxInit_real_T(&fg2, 2);
    i0 = fg2->size[0] * fg2->size[1];
    fg2->size[0] = 1;
    fg2->size[1] = (int)((unsigned int)F3->size[1] + fg->size[1]);
    emxEnsureCapacity((emxArray__common *)fg2, i0, sizeof(double));
    loop_ub = (int)((unsigned int)F3->size[1] + fg->size[1]);
    for (i0 = 0; i0 < loop_ub; i0++) {
        fg2->data[i0] = 0.0;
    }

    loop_ub = fg->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
        fg2->data[i0] = fg->data[fg->size[0] * i0];
    }

    if (rtIsNaN(Gstop)) {
        i0 = F3->size[0] * F3->size[1];
        F3->size[0] = 1;
        F3->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)F3, i0, sizeof(double));
        F3->data[0] = rtNaN;
    } else if (8192.0 < Gstop) {
        i0 = F3->size[0] * F3->size[1];
        F3->size[0] = 1;
        F3->size[1] = 0;
        emxEnsureCapacity((emxArray__common *)F3, i0, sizeof(double));
    } else if (rtIsInf(Gstop) && (Gstop == 8192.0)) {
        i0 = F3->size[0] * F3->size[1];
        F3->size[0] = 1;
        F3->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)F3, i0, sizeof(double));
        F3->data[0] = rtNaN;
    } else if (Gstop == Gstop) {
        i0 = F3->size[0] * F3->size[1];
        F3->size[0] = 1;
        F3->size[1] = (int)(8192.0 - Gstop) + 1;
        emxEnsureCapacity((emxArray__common *)F3, i0, sizeof(double));
        loop_ub = (int)(8192.0 - Gstop);
        for (i0 = 0; i0 <= loop_ub; i0++) {
            F3->data[F3->size[0] * i0] = Gstop + (double)i0;
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

        i0 = F3->size[0] * F3->size[1];
        F3->size[0] = 1;
        F3->size[1] = n;
        emxEnsureCapacity((emxArray__common *)F3, i0, sizeof(double));
        if (n > 0) {
            F3->data[0] = Gstop;
            if (n > 1) {
                F3->data[n - 1] = apnd;
                dec_int3 = (n - 1) / 2;
                for (d_combinedResponse = 1; d_combinedResponse < dec_int3;
                     d_combinedResponse++) {
                    F3->data[d_combinedResponse] = Gstop + (double)d_combinedResponse;
                    F3->data[(n - d_combinedResponse) - 1] = apnd - (double)
                            d_combinedResponse;
                }

                if (dec_int3 << 1 == n - 1) {
                    F3->data[dec_int3] = (Gstop + apnd) / 2.0;
                } else {
                    F3->data[dec_int3] = Gstop + (double)dec_int3;
                    F3->data[dec_int3 + 1] = apnd - (double)dec_int3;
                }
            }
        }
    }

    emxInit_real_T(&omega2, 2);
    i0 = omega2->size[0] * omega2->size[1];
    omega2->size[0] = 1;
    omega2->size[1] = (int)((unsigned int)F3->size[1] + omega->size[1]);
    emxEnsureCapacity((emxArray__common *)omega2, i0, sizeof(double));
    loop_ub = (int)((unsigned int)F3->size[1] + omega->size[1]);
    for (i0 = 0; i0 < loop_ub; i0++) {
        omega2->data[i0] = 0.0;
    }

    loop_ub = omega->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
        omega2->data[i0] = omega->data[omega->size[0] * i0];
    }

    if (rtIsNaN(Gstop)) {
        i0 = F3->size[0] * F3->size[1];
        F3->size[0] = 1;
        F3->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)F3, i0, sizeof(double));
        F3->data[0] = rtNaN;
    } else if (8192.0 < Gstop) {
        i0 = F3->size[0] * F3->size[1];
        F3->size[0] = 1;
        F3->size[1] = 0;
        emxEnsureCapacity((emxArray__common *)F3, i0, sizeof(double));
    } else if (rtIsInf(Gstop) && (Gstop == 8192.0)) {
        i0 = F3->size[0] * F3->size[1];
        F3->size[0] = 1;
        F3->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)F3, i0, sizeof(double));
        F3->data[0] = rtNaN;
    } else if (Gstop == Gstop) {
        i0 = F3->size[0] * F3->size[1];
        F3->size[0] = 1;
        F3->size[1] = (int)(8192.0 - Gstop) + 1;
        emxEnsureCapacity((emxArray__common *)F3, i0, sizeof(double));
        loop_ub = (int)(8192.0 - Gstop);
        for (i0 = 0; i0 <= loop_ub; i0++) {
            F3->data[F3->size[0] * i0] = Gstop + (double)i0;
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

        i0 = F3->size[0] * F3->size[1];
        F3->size[0] = 1;
        F3->size[1] = n;
        emxEnsureCapacity((emxArray__common *)F3, i0, sizeof(double));
        if (n > 0) {
            F3->data[0] = Gstop;
            if (n > 1) {
                F3->data[n - 1] = apnd;
                dec_int3 = (n - 1) / 2;
                for (d_combinedResponse = 1; d_combinedResponse < dec_int3;
                     d_combinedResponse++) {
                    F3->data[d_combinedResponse] = Gstop + (double)d_combinedResponse;
                    F3->data[(n - d_combinedResponse) - 1] = apnd - (double)
                            d_combinedResponse;
                }

                if (dec_int3 << 1 == n - 1) {
                    F3->data[dec_int3] = (Gstop + apnd) / 2.0;
                } else {
                    F3->data[dec_int3] = Gstop + (double)dec_int3;
                    F3->data[dec_int3 + 1] = apnd - (double)dec_int3;
                }
            }
        }
    }

    emxInit_creal_T(&rgN, 2);
    i0 = rgN->size[0] * rgN->size[1];
    rgN->size[0] = 1;
    rgN->size[1] = (int)((unsigned int)F3->size[1] + rg->size[1]);
    emxEnsureCapacity((emxArray__common *)rgN, i0, sizeof(creal_T));
    loop_ub = (int)((unsigned int)F3->size[1] + rg->size[1]);
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
    emxEnsureCapacity((emxArray__common *)b_omega2, i3, sizeof(double));
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
    emxEnsureCapacity((emxArray__common *)c_omega2, i3, sizeof(double));
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
    emxEnsureCapacity((emxArray__common *)e_combinedResponse, i0, sizeof(creal_T));
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
        emxEnsureCapacity((emxArray__common *)x, i0, sizeof(double));
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
    emxEnsureCapacity((emxArray__common *)fg, i0, sizeof(double));
    for (d_combinedResponse = 0; d_combinedResponse + 1 <= dec_int3;
         d_combinedResponse++) {
        if ((omega->data[d_combinedResponse] > sigma) || rtIsNaN(sigma)) {
            sigmax = omega->data[d_combinedResponse];
        } else {
            sigmax = sigma;
        }

        fg->data[d_combinedResponse] = sigmax;
    }

    if (phEQ == -1.0) {
        b_abs(rgN, F3);
        i0 = rgN->size[0] * rgN->size[1];
        rgN->size[0] = 1;
        rgN->size[1] = F3->size[1];
        emxEnsureCapacity((emxArray__common *)rgN, i0, sizeof(creal_T));
        loop_ub = F3->size[0] * F3->size[1];
        for (i0 = 0; i0 < loop_ub; i0++) {
            rgN->data[i0].re = F3->data[i0];
            rgN->data[i0].im = 0.0;
        }
    }

    emxInit_real_T(&weight, 2);
    rdivide(sw, dBinv(Apass / 2.0) - 1.0, F3);
    i0 = weight->size[0] * weight->size[1];
    weight->size[0] = 1;
    weight->size[1] = F3->size[1] + fg->size[1];
    emxEnsureCapacity((emxArray__common *)weight, i0, sizeof(double));
    loop_ub = F3->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
        weight->data[weight->size[0] * i0] = F3->data[F3->size[0] * i0];
    }

    loop_ub = fg->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
        weight->data[weight->size[0] * (i0 + F3->size[1])] = fg->data[fg->size[0] *
                i0];
    }

    dec_int3 = 1;
    n = weight->size[1];
    sigma = weight->data[0];
    if (weight->size[1] > 1) {
        if (rtIsNaN(weight->data[0])) {
            d_combinedResponse = 2;
            exitg1 = false;
            while ((!exitg1) && (d_combinedResponse <= n)) {
                dec_int3 = d_combinedResponse;
                if (!rtIsNaN(weight->data[d_combinedResponse - 1])) {
                    sigma = weight->data[d_combinedResponse - 1];
                    exitg1 = true;
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
    emxEnsureCapacity((emxArray__common *)b_weight, i0, sizeof(double));
    loop_ub = weight->size[0] * weight->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
        b_weight->data[i0] = weight->data[i0];
    }

    rdivide(b_weight, sigma, weight);

    /*  Set up design for FIR filter */
    i0 = fg->size[0] * fg->size[1];
    fg->size[0] = 1;
    fg->size[1] = rgN->size[1];
    emxEnsureCapacity((emxArray__common *)fg, i0, sizeof(double));
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
    emxEnsureCapacity((emxArray__common *)F1, i0, sizeof(double));
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
    emxEnsureCapacity((emxArray__common *)F2, i3, sizeof(double));
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
    emxEnsureCapacity((emxArray__common *)A1, i0, sizeof(double));
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
    emxEnsureCapacity((emxArray__common *)A2, i3, sizeof(double));
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
    emxEnsureCapacity((emxArray__common *)W1, i0, sizeof(double));
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
    emxEnsureCapacity((emxArray__common *)W2, i3, sizeof(double));
    loop_ub = i2 - i0;
    for (i2 = 0; i2 < loop_ub; i2++) {
        W2->data[W2->size[0] * i2] = weight->data[i0 + i2];
    }

    /*  Determine the number of taps for FIR */
    if (b_strcmp(RxTx)) {
        if (hb3 == 1) {
            clkFIR = 16.0 * floor(converter_rate / Rdata);
            if (!(clkFIR < 128.0)) {
                clkFIR = 128.0;
            }
        } else {
            clkFIR = 16.0 * floor(converter_rate / (2.0 * Rdata));
            if (!(clkFIR < 128.0)) {
                clkFIR = 128.0;
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
        clkFIR = Nmax;
        if (sigma < clkFIR) {
            clkFIR = sigma;
        }
    }

    apnd = clkFIR / 16.0;
    tap_store_size[0] = (int)apnd;
    tap_store_size[1] = (int)clkFIR;
    loop_ub = (int)apnd * (int)clkFIR;
    for (i0 = 0; i0 < loop_ub; i0++) {
        tap_store_data[i0] = 0.0;
    }

    loop_ub = (int)(clkFIR / 16.0);
    for (i0 = 0; i0 < loop_ub; i0++) {
        Apass_actual_vector_data[i0] = 0.0;
    }

    loop_ub = (int)(clkFIR / 16.0);
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
            b1[0] = F1->data[0];
            b1[1] = F1->data[F1->size[1] - 1];
            b1[2] = F2->data[0];
            b1[3] = F2->data[F2->size[1] - 1];
            i0 = b_A1->size[0] * b_A1->size[1];
            b_A1->size[0] = 1;
            b_A1->size[1] = A1->size[1] + A2->size[1];
            emxEnsureCapacity((emxArray__common *)b_A1, i0, sizeof(double));
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
            emxEnsureCapacity((emxArray__common *)b_F1, i0, sizeof(double));
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
            emxEnsureCapacity((emxArray__common *)b_W1, i0, sizeof(double));
            loop_ub = W1->size[1];
            for (i0 = 0; i0 < loop_ub; i0++) {
                b_W1->data[b_W1->size[0] * i0] = W1->data[W1->size[0] * i0];
            }

            loop_ub = W2->size[1];
            for (i0 = 0; i0 < loop_ub; i0++) {
                b_W1->data[b_W1->size[0] * (i0 + W1->size[1])] = W2->data[W2->size[0] *
                        i0];
            }

            firpm_cg(clkFIR - 1.0, b1, b_A1, b_F1, b_W1, ccoef);
        } else {
            /*  Check different designs until we reach required ripple condition */
            sigma = rt_powd_snf(10.0, -Astop / 20.0);

            /*  Peak Ripple */
            i0 = ccoef->size[0] * ccoef->size[1];
            ccoef->size[0] = 1;
            ccoef->size[1] = 1;
            emxEnsureCapacity((emxArray__common *)ccoef, i0, sizeof(double));
            ccoef->data[0] = 0.0;

            /*  Predef type */
            d_combinedResponse = 0;
            exitg1 = false;
            while ((!exitg1) && (d_combinedResponse < 126)) {
                b1[0] = F1->data[0];
                b1[1] = F1->data[F1->size[1] - 1];
                b1[2] = F2->data[0];
                b1[3] = F2->data[F2->size[1] - 1];
                i0 = c_A1->size[0] * c_A1->size[1];
                c_A1->size[0] = 1;
                c_A1->size[1] = A1->size[1] + A2->size[1];
                emxEnsureCapacity((emxArray__common *)c_A1, i0, sizeof(double));
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
                emxEnsureCapacity((emxArray__common *)c_F1, i0, sizeof(double));
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
                emxEnsureCapacity((emxArray__common *)c_W1, i0, sizeof(double));
                loop_ub = W1->size[1];
                for (i0 = 0; i0 < loop_ub; i0++) {
                    c_W1->data[c_W1->size[0] * i0] = W1->data[W1->size[0] * i0];
                }

                loop_ub = W2->size[1];
                for (i0 = 0; i0 < loop_ub; i0++) {
                    c_W1->data[c_W1->size[0] * (i0 + W1->size[1])] = W2->data[W2->size[0] *
                            i0];
                }

                b_firpm_cg(3.0 + (double)d_combinedResponse, b1, c_A1, c_F1, c_W1, ccoef,
                           &valid, &sigmax);

                /*  Check if design meets specs */
                if ((sigmax < sigma) && valid) {
                    exitg1 = true;
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
            emxEnsureCapacity((emxArray__common *)fg, d_combinedResponse, sizeof
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
            emxEnsureCapacity((emxArray__common *)omega, d_combinedResponse, sizeof
                              (double));
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
            emxEnsureCapacity((emxArray__common *)sw, d_combinedResponse, sizeof
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
            emxEnsureCapacity((emxArray__common *)F3, d_combinedResponse, sizeof
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
            emxEnsureCapacity((emxArray__common *)F4, n, sizeof(double));
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
                sigma = clkFIR - 1.0;
            } else {
                sigma = (double)ccoef->size[1] - 1.0;
            }

            b1[0] = F3->data[0];
            b1[1] = F3->data[F3->size[1] - 1];
            b1[2] = F4->data[0];
            b1[3] = F4->data[F4->size[1] - 1];
            i2 = b_omega->size[0] * b_omega->size[1];
            b_omega->size[0] = 1;
            b_omega->size[1] = (loop_ub + dec_int3) - d_combinedResponse;
            emxEnsureCapacity((emxArray__common *)b_omega, i2, sizeof(double));
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
            emxEnsureCapacity((emxArray__common *)b_F3, i2, sizeof(double));
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
            emxEnsureCapacity((emxArray__common *)b_sw, i2, sizeof(double));
            for (i2 = 0; i2 < hb3; i2++) {
                b_sw->data[b_sw->size[0] * i2] = sw->data[i2];
            }

            loop_ub = i0 - n;
            for (i0 = 0; i0 <= loop_ub; i0++) {
                b_sw->data[b_sw->size[0] * (i0 + hb3)] = sw->data[n + i0];
            }

            firpm_cg(sigma, b1, b_omega, b_F3, b_sw, fg);
            i0 = fg->size[1];
            for (d_combinedResponse = 0; d_combinedResponse < i0; d_combinedResponse++) {
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
            emxEnsureCapacity((emxArray__common *)fg, i0, sizeof(double));
            loop_ub = (int)b2[1];
            for (i0 = 0; i0 < loop_ub; i0++) {
                fg->data[i0] = 0.0;
            }
        }

        loop_ub = ccoef->size[1];
        for (i0 = 0; i0 < loop_ub; i0++) {
            tap_store_data[i + (int)apnd * i0] = ccoef->data[ccoef->size[0] * i0] +
                                                 fg->data[fg->size[0] * i0];
        }

        /*  scoef ==0 when no EQ */
        determineBestFractionLength(tap_store_data, tap_store_size, i + 1,
                                    ccoef->size[1], F3);
        loop_ub = F3->size[1];
        for (i0 = 0; i0 < loop_ub; i0++) {
            tap_store_data[i + (int)apnd * i0] = F3->data[F3->size[0] * i0];
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
            emxEnsureCapacity((emxArray__common *)k_omega2, i0, sizeof(double));
            for (i0 = 0; i0 < loop_ub; i0++) {
                k_omega2->data[k_omega2->size[0] * i0] = omega2->data[i0];
            }

            b_tap_store_size[0] = 1;
            b_tap_store_size[1] = hb3;
            for (i0 = 0; i0 < hb3; i0++) {
                firTapsPreScale[i0] = tap_store_data[i + (int)apnd * i0];
            }

            c_generateCascadedResponseRx(enables, k_omega2, converter_rate, hb1_coeff,
                                         hb2_coeff, hb3_coeff_data, hb3_coeff_size, dec_int3_coeff_data,
                                         dec_int3_coeff_size, firTapsPreScale, b_tap_store_size, rg2);
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
            emxEnsureCapacity((emxArray__common *)j_omega2, i3, sizeof(double));
            hb3 = i2 - i0;
            for (i2 = 0; i2 < hb3; i2++) {
                j_omega2->data[j_omega2->size[0] * i2] = omega2->data[i0 + i2];
            }

            d_tap_store_size[0] = 1;
            d_tap_store_size[1] = loop_ub;
            for (i0 = 0; i0 < loop_ub; i0++) {
                firTapsPreScale[i0] = tap_store_data[i + (int)apnd * i0];
            }

            c_generateCascadedResponseRx(enables, j_omega2, converter_rate, hb1_coeff,
                                         hb2_coeff, hb3_coeff_data, hb3_coeff_size, dec_int3_coeff_data,
                                         dec_int3_coeff_size, firTapsPreScale, d_tap_store_size, rg);
            if (1.0 > Gpass + 1.0) {
                loop_ub = 0;
            } else {
                loop_ub = (int)(Gpass + 1.0);
            }

            i0 = i_omega2->size[0] * i_omega2->size[1];
            i_omega2->size[0] = 1;
            i_omega2->size[1] = loop_ub;
            emxEnsureCapacity((emxArray__common *)i_omega2, i0, sizeof(double));
            for (i0 = 0; i0 < loop_ub; i0++) {
                i_omega2->data[i_omega2->size[0] * i0] = omega2->data[i0];
            }

            c_analogresp(i_omega2, converter_rate, b1_data, b1_size, a1, b2_data,
                         b2_size, a2, r0);
            i0 = r2->size[0] * r2->size[1];
            r2->size[0] = 1;
            r2->size[1] = r0->size[1];
            emxEnsureCapacity((emxArray__common *)r2, i0, sizeof(creal_T));
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
            emxEnsureCapacity((emxArray__common *)h_omega2, i3, sizeof(double));
            loop_ub = i2 - i0;
            for (i2 = 0; i2 < loop_ub; i2++) {
                h_omega2->data[h_omega2->size[0] * i2] = omega2->data[i0 + i2];
            }

            c_analogresp(h_omega2, converter_rate, b1_data, b1_size, a1, b2_data,
                         b2_size, a2, r0);
            i0 = r1->size[0] * r1->size[1];
            r1->size[0] = 1;
            r1->size[1] = r0->size[1];
            emxEnsureCapacity((emxArray__common *)r1, i0, sizeof(creal_T));
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
            emxEnsureCapacity((emxArray__common *)g_omega2, i0, sizeof(double));
            for (i0 = 0; i0 < loop_ub; i0++) {
                g_omega2->data[g_omega2->size[0] * i0] = omega2->data[i0];
            }

            c_tap_store_size[0] = 1;
            c_tap_store_size[1] = hb3;
            for (i0 = 0; i0 < hb3; i0++) {
                firTapsPreScale[i0] = tap_store_data[i + (int)apnd * i0];
            }

            c_generateCascadedResponseRx(enables, g_omega2, converter_rate, hb1_coeff,
                                         hb2_coeff, hb3_coeff_data, hb3_coeff_size, dec_int3_coeff_data,
                                         dec_int3_coeff_size, firTapsPreScale, c_tap_store_size, rg2);
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
            emxEnsureCapacity((emxArray__common *)f_omega2, i3, sizeof(double));
            hb3 = i2 - i0;
            for (i2 = 0; i2 < hb3; i2++) {
                f_omega2->data[f_omega2->size[0] * i2] = omega2->data[i0 + i2];
            }

            e_tap_store_size[0] = 1;
            e_tap_store_size[1] = loop_ub;
            for (i0 = 0; i0 < loop_ub; i0++) {
                firTapsPreScale[i0] = tap_store_data[i + (int)apnd * i0];
            }

            c_generateCascadedResponseRx(enables, f_omega2, converter_rate, hb1_coeff,
                                         hb2_coeff, hb3_coeff_data, hb3_coeff_size, dec_int3_coeff_data,
                                         dec_int3_coeff_size, firTapsPreScale, e_tap_store_size, rg);
            if (1.0 > Gpass + 1.0) {
                loop_ub = 0;
            } else {
                loop_ub = (int)(Gpass + 1.0);
            }

            i0 = e_omega2->size[0] * e_omega2->size[1];
            e_omega2->size[0] = 1;
            e_omega2->size[1] = loop_ub;
            emxEnsureCapacity((emxArray__common *)e_omega2, i0, sizeof(double));
            for (i0 = 0; i0 < loop_ub; i0++) {
                e_omega2->data[e_omega2->size[0] * i0] = omega2->data[i0];
            }

            d_analogresp(e_omega2, converter_rate, b1_data, b1_size, a1, b2_data,
                         b2_size, a2, r0);
            i0 = b_rg2->size[0] * b_rg2->size[1];
            b_rg2->size[0] = 1;
            b_rg2->size[1] = rg2->size[1];
            emxEnsureCapacity((emxArray__common *)b_rg2, i0, sizeof(creal_T));
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
            emxEnsureCapacity((emxArray__common *)d_omega2, i3, sizeof(double));
            loop_ub = i2 - i0;
            for (i2 = 0; i2 < loop_ub; i2++) {
                d_omega2->data[d_omega2->size[0] * i2] = omega2->data[i0 + i2];
            }

            d_analogresp(d_omega2, converter_rate, b1_data, b1_size, a1, b2_data,
                         b2_size, a2, r0);
            i0 = b_rg->size[0] * b_rg->size[1];
            b_rg->size[0] = 1;
            b_rg->size[1] = rg->size[1];
            emxEnsureCapacity((emxArray__common *)b_rg, i0, sizeof(creal_T));
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
                exitg1 = false;
                while ((!exitg1) && (d_combinedResponse <= n)) {
                    dec_int3 = d_combinedResponse;
                    if (!rtIsNaN(fg->data[d_combinedResponse - 1])) {
                        sigma = fg->data[d_combinedResponse - 1];
                        exitg1 = true;
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
                exitg1 = false;
                while ((!exitg1) && (d_combinedResponse <= n)) {
                    dec_int3 = d_combinedResponse;
                    if (!rtIsNaN(fg->data[d_combinedResponse - 1])) {
                        sigmax = fg->data[d_combinedResponse - 1];
                        exitg1 = true;
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
                exitg1 = false;
                while ((!exitg1) && (d_combinedResponse <= n)) {
                    dec_int3 = d_combinedResponse;
                    if (!rtIsNaN(omega->data[d_combinedResponse - 1])) {
                        sigma = omega->data[d_combinedResponse - 1];
                        exitg1 = true;
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
            emxEnsureCapacity((emxArray__common *)h, i0, sizeof(double));
            for (i0 = 0; i0 < loop_ub; i0++) {
                h->data[h->size[0] * i0] = tap_store_data[(int)apnd * i0];
            }

            /* Apass_actual = Apass_actual_vector(1); */
            /* Astop_actual = Astop_actual_vector(1); */
            exitg2 = 1;
        } else if ((Apass_actual_vector_data[0] > Apass) ||
                   (Astop_actual_vector_data[0] < Astop)) {
            if (1.0 > clkFIR) {
                loop_ub = 0;
            } else {
                loop_ub = (int)clkFIR;
            }

            i0 = h->size[0] * h->size[1];
            h->size[0] = 1;
            h->size[1] = loop_ub;
            emxEnsureCapacity((emxArray__common *)h, i0, sizeof(double));
            for (i0 = 0; i0 < loop_ub; i0++) {
                h->data[h->size[0] * i0] = tap_store_data[(int)apnd * i0];
            }

            /* Apass_actual = Apass_actual_vector(1); */
            /* Astop_actual = Astop_actual_vector(1); */
            exitg2 = 1;
        } else if ((Apass_actual_vector_data[i] > Apass) ||
                   (Astop_actual_vector_data[i] < Astop)) {
            if (1.0 > clkFIR + 16.0) {
                loop_ub = 0;
            } else {
                loop_ub = (int)(clkFIR + 16.0);
            }

            i0 = h->size[0] * h->size[1];
            h->size[0] = 1;
            h->size[1] = loop_ub;
            emxEnsureCapacity((emxArray__common *)h, i0, sizeof(double));
            for (i0 = 0; i0 < loop_ub; i0++) {
                h->data[h->size[0] * i0] = tap_store_data[(i + (int)apnd * i0) - 1];
            }

            /* Apass_actual = Apass_actual_vector(i-1); */
            /* Astop_actual = Astop_actual_vector(i-1); */
            exitg2 = 1;
        } else {
            clkFIR -= 16.0;
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
                emxEnsureCapacity((emxArray__common *)r4, i0, sizeof(double));
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
                emxEnsureCapacity((emxArray__common *)h, i0, sizeof(double));
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
                emxEnsureCapacity((emxArray__common *)r3, i0, sizeof(double));
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
                emxEnsureCapacity((emxArray__common *)h, i0, sizeof(double));
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

    sigma = ceil(b_log2(sigma));
    switch ((int)(1.0 + sigma)) {
    case 2:
        dec_int3 = 6;
        break;

    case 1:
        dec_int3 = 0;
        break;

    case 0:
        dec_int3 = -6;
        break;

    default:
        dec_int3 = -12;
        break;
    }

    if (b_strcmp(RxTx)) {
        if (1.0 + sigma > 2.0) {
            dec_int3 = 6;
        }
    } else {
        if (FIR == 2.0) {
            dec_int3 += 6;
        } else {
            if (FIR == 4.0) {
                dec_int3 += 12;
            }
        }

        if (dec_int3 > 0) {
            dec_int3 = 0;
        } else {
            if (dec_int3 < -6) {
                dec_int3 = -6;
            }
        }
    }

    /*  Scale taps */
    memcpy(&b_firTapsPreScale[0], &firTapsPreScale[0], sizeof(double) << 7);
    b_determineBestFractionLength(b_firTapsPreScale, firTapsPreScale);
    sigmax = rt_powd_snf(2.0, 16.0 - (1.0 + sigma));
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
    *numOutputTaps = h->size[1];
    *filterGain = dec_int3;
    emxFree_real_T(&h);
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
