/*
 * Sponsored Third Party Support License -- for use only to support
 * products interfaced to MathWorks software under terms specified in your
 * company's restricted use license agreement.
 * File: internal_design_filter_cg_emxutil.c
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 05-Jan-2018 14:44:11
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "internal_design_filter_cg.h"
#include "internal_design_filter_cg_emxutil.h"

/* Function Definitions */

/*
 * Arguments    : emxArray__common *emxArray
 *                int oldNumel
 *                unsigned int elementSize
 * Return Type  : void
 */
void emxEnsureCapacity(emxArray__common *emxArray, int oldNumel, unsigned int
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

        newData = calloc((unsigned int)i, elementSize);
        if (emxArray->data != NULL) {
            memcpy(newData, emxArray->data, elementSize * oldNumel);
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
