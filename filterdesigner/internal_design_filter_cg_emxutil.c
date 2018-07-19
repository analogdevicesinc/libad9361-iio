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

/* Function Definitions */

/*
 * Arguments    : emxArray_cint8_T *emxArray
 *                int oldNumel
 * Return Type  : void
 */
void emxEnsureCapacity_cint8_T(emxArray_cint8_T *emxArray, int oldNumel)
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

        newData = calloc((unsigned int)i, sizeof(cint8_T));
        if (emxArray->data != NULL) {
            memcpy(newData, (void *)emxArray->data, sizeof(cint8_T) * oldNumel);
            if (emxArray->canFreeData) {
                free((void *)emxArray->data);
            }
        }

        emxArray->data = (cint8_T *)newData;
        emxArray->allocatedSize = i;
        emxArray->canFreeData = true;
    }
}

/*
 * Arguments    : emxArray_creal_T *emxArray
 *                int oldNumel
 * Return Type  : void
 */
void emxEnsureCapacity_creal_T(emxArray_creal_T *emxArray, int oldNumel)
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

        newData = calloc((unsigned int)i, sizeof(creal_T));
        if (emxArray->data != NULL) {
            memcpy(newData, (void *)emxArray->data, sizeof(creal_T) * oldNumel);
            if (emxArray->canFreeData) {
                free((void *)emxArray->data);
            }
        }

        emxArray->data = (creal_T *)newData;
        emxArray->allocatedSize = i;
        emxArray->canFreeData = true;
    }
}

/*
 * Arguments    : emxArray_creal_T *emxArray
 *                int oldNumel
 * Return Type  : void
 */
void emxEnsureCapacity_creal_T1(emxArray_creal_T *emxArray, int oldNumel)
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

        newData = calloc((unsigned int)i, sizeof(creal_T));
        if (emxArray->data != NULL) {
            memcpy(newData, (void *)emxArray->data, sizeof(creal_T) * oldNumel);
            if (emxArray->canFreeData) {
                free((void *)emxArray->data);
            }
        }

        emxArray->data = (creal_T *)newData;
        emxArray->allocatedSize = i;
        emxArray->canFreeData = true;
    }
}

/*
 * Arguments    : emxArray_int32_T *emxArray
 *                int oldNumel
 * Return Type  : void
 */
void emxEnsureCapacity_int32_T(emxArray_int32_T *emxArray, int oldNumel)
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

        newData = calloc((unsigned int)i, sizeof(int));
        if (emxArray->data != NULL) {
            memcpy(newData, (void *)emxArray->data, sizeof(int) * oldNumel);
            if (emxArray->canFreeData) {
                free((void *)emxArray->data);
            }
        }

        emxArray->data = (int *)newData;
        emxArray->allocatedSize = i;
        emxArray->canFreeData = true;
    }
}

/*
 * Arguments    : emxArray_int32_T *emxArray
 *                int oldNumel
 * Return Type  : void
 */
void emxEnsureCapacity_int32_T1(emxArray_int32_T *emxArray, int oldNumel)
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

        newData = calloc((unsigned int)i, sizeof(int));
        if (emxArray->data != NULL) {
            memcpy(newData, (void *)emxArray->data, sizeof(int) * oldNumel);
            if (emxArray->canFreeData) {
                free((void *)emxArray->data);
            }
        }

        emxArray->data = (int *)newData;
        emxArray->allocatedSize = i;
        emxArray->canFreeData = true;
    }
}

/*
 * Arguments    : emxArray_int8_T *emxArray
 *                int oldNumel
 * Return Type  : void
 */
void emxEnsureCapacity_int8_T(emxArray_int8_T *emxArray, int oldNumel)
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

        newData = calloc((unsigned int)i, sizeof(signed char));
        if (emxArray->data != NULL) {
            memcpy(newData, (void *)emxArray->data, sizeof(signed char) * oldNumel);
            if (emxArray->canFreeData) {
                free((void *)emxArray->data);
            }
        }

        emxArray->data = (signed char *)newData;
        emxArray->allocatedSize = i;
        emxArray->canFreeData = true;
    }
}

/*
 * Arguments    : emxArray_real_T *emxArray
 *                int oldNumel
 * Return Type  : void
 */
void emxEnsureCapacity_real_T(emxArray_real_T *emxArray, int oldNumel)
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

        newData = calloc((unsigned int)i, sizeof(double));
        if (emxArray->data != NULL) {
            memcpy(newData, (void *)emxArray->data, sizeof(double) * oldNumel);
            if (emxArray->canFreeData) {
                free((void *)emxArray->data);
            }
        }

        emxArray->data = (double *)newData;
        emxArray->allocatedSize = i;
        emxArray->canFreeData = true;
    }
}

/*
 * Arguments    : emxArray_real_T *emxArray
 *                int oldNumel
 * Return Type  : void
 */
void emxEnsureCapacity_real_T1(emxArray_real_T *emxArray, int oldNumel)
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

        newData = calloc((unsigned int)i, sizeof(double));
        if (emxArray->data != NULL) {
            memcpy(newData, (void *)emxArray->data, sizeof(double) * oldNumel);
            if (emxArray->canFreeData) {
                free((void *)emxArray->data);
            }
        }

        emxArray->data = (double *)newData;
        emxArray->allocatedSize = i;
        emxArray->canFreeData = true;
    }
}

/*
 * Arguments    : emxArray_cint8_T **pEmxArray
 * Return Type  : void
 */
void emxFree_cint8_T(emxArray_cint8_T **pEmxArray)
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
void emxFree_creal_T(emxArray_creal_T **pEmxArray)
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
void emxFree_int32_T(emxArray_int32_T **pEmxArray)
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
void emxFree_int8_T(emxArray_int8_T **pEmxArray)
{
    if (*pEmxArray != (emxArray_int8_T *)NULL) {
        if (((*pEmxArray)->data != (signed char *)NULL) && (*pEmxArray)->canFreeData) {
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
void emxFree_real_T(emxArray_real_T **pEmxArray)
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
 * Arguments    : emxArray_cint8_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
void emxInit_cint8_T(emxArray_cint8_T **pEmxArray, int numDimensions)
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
void emxInit_creal_T(emxArray_creal_T **pEmxArray, int numDimensions)
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
void emxInit_creal_T1(emxArray_creal_T **pEmxArray, int numDimensions)
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
void emxInit_int32_T(emxArray_int32_T **pEmxArray, int numDimensions)
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
void emxInit_int32_T1(emxArray_int32_T **pEmxArray, int numDimensions)
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
void emxInit_int8_T(emxArray_int8_T **pEmxArray, int numDimensions)
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
void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions)
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
void emxInit_real_T1(emxArray_real_T **pEmxArray, int numDimensions)
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
 * File trailer for internal_design_filter_cg_emxutil.c
 *
 * [EOF]
 */
