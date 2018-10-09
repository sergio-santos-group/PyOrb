//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdint.h>
#include "cube.h"

/* Docstrings */
static char   module_docstring[] = "module_docstring.";
static char get_cube_docstring[] = "get_cube_docstring.";

/* Available functions */
static PyObject *cube_get_cube(PyObject *self, PyObject *args);

/* Module specification */
static PyMethodDef module_methods[] = {
    {"get_cube", cube_get_cube, METH_VARARGS, get_cube_docstring},
    {NULL, NULL, 0, NULL}
};

/* Initialize the module */
PyMODINIT_FUNC init_cube(void) {

  PyObject *m = Py_InitModule3("_cube", module_methods, module_docstring);
  if (m == NULL)
    return;
  import_array();   // Load numpy functionality
}

static PyObject *cube_get_cube(PyObject *self, PyObject *args)
{
  PyObject *obj1, *obj2, *obj3, *obj4, *obj5, *obj6, *obj7, *obj8, *obj9, *obj10, *obj11;
  
  if (!PyArg_ParseTuple(args, "OOOOOOOOOOO", &obj1, &obj2, &obj3, &obj4,
                                             &obj5, &obj6, &obj7, &obj8,
                                             &obj9, &obj10, &obj11)) {
    return NULL;
  }
  
  /* Interpret the input objects as numpy arrays. */
  PyObject *darr1 = PyArray_FROM_OTF(obj1, NPY_DOUBLE, NPY_IN_ARRAY);
  PyObject *darr2 = PyArray_FROM_OTF(obj2, NPY_DOUBLE, NPY_IN_ARRAY);
  PyObject *darr3 = PyArray_FROM_OTF(obj3, NPY_DOUBLE, NPY_IN_ARRAY);
  PyObject *darr4 = PyArray_FROM_OTF(obj4, NPY_DOUBLE, NPY_IN_ARRAY);
  PyObject *darr5 = PyArray_FROM_OTF(obj5, NPY_DOUBLE, NPY_IN_ARRAY);
  PyObject *darr6 = PyArray_FROM_OTF(obj6, NPY_DOUBLE, NPY_IN_ARRAY);
  
  PyObject *iarr7 = PyArray_FROM_OTF(obj7, NPY_INT64, NPY_IN_ARRAY);
  PyObject *iarr8 = PyArray_FROM_OTF(obj8, NPY_INT64, NPY_IN_ARRAY);
  PyObject *iarr9 = PyArray_FROM_OTF(obj9, NPY_INT64, NPY_IN_ARRAY);
  PyObject *iarr10= PyArray_FROM_OTF(obj10,NPY_INT64, NPY_IN_ARRAY);
  PyObject *iarr11= PyArray_FROM_OTF(obj11,NPY_INT64, NPY_IN_ARRAY);
  
  /* If that didn't work, throw an exception. */
  if (darr1==NULL || darr2==NULL || darr3==NULL || darr4==NULL || darr5==NULL || darr6==NULL ||
      iarr7==NULL || iarr8==NULL || iarr9==NULL || iarr10==NULL || iarr11==NULL) {
    Py_XDECREF(darr1);
    Py_XDECREF(darr2);
    Py_XDECREF(darr3);
    Py_XDECREF(darr4);
    Py_XDECREF(darr5);
    Py_XDECREF(darr6);
    Py_XDECREF(iarr7);
    Py_XDECREF(iarr8);
    Py_XDECREF(iarr9);
    Py_XDECREF(iarr10);
    Py_XDECREF(iarr11);
    return NULL;
  }
  
  int64_t nCubePts = (int64_t)PyArray_DIM(darr2, 0);  // number of grid points
  double *grid = (double*)PyArray_DATA(darr1);
  double *cube = (double*)PyArray_DATA(darr2);
  
  int64_t nCoords = (int64_t)PyArray_DIM(darr3, 0)/3; // number of TOMS
  double *coords  = (double*)PyArray_DATA(darr3);
  
  double *coefs = (double*)PyArray_DATA(darr4);
  double *gtoA  = (double*)PyArray_DATA(darr5);
  double *gtoCN = (double*)PyArray_DATA(darr6);
  
  int64_t nSymmetry = (int64_t)PyArray_DIM(iarr7, 0);
  
  int64_t *symmetry  = (int64_t*)PyArray_DATA(iarr7);
  int64_t *atomIdxs  = (int64_t*)PyArray_DATA(iarr8);
  int64_t *moIdxs  = (int64_t*)PyArray_DATA(iarr9);
  int64_t *gtoIdxs = (int64_t*)PyArray_DATA(iarr10);
  int64_t *cIdxs   = (int64_t*)PyArray_DATA(iarr11);
  
  /* Call the external C function. */
  get_cube(grid, cube, nCubePts,
           coords, nCoords,
           symmetry, nSymmetry,
           atomIdxs, coefs, gtoA, gtoCN,
           moIdxs, gtoIdxs, cIdxs);
  
  /* Clean up. */
  Py_DECREF(darr1);
  Py_DECREF(darr2);
  Py_DECREF(darr3);
  Py_DECREF(darr4);
  Py_DECREF(darr5);
  Py_DECREF(darr6);
  Py_DECREF(iarr7);
  Py_DECREF(iarr8);
  Py_DECREF(iarr9);
  Py_DECREF(iarr10);
  Py_DECREF(iarr11);
  
  /* Build the output tuple */
  PyObject *ret = Py_BuildValue("i", 1);
  return ret;
  
  //   if (value < 0.0) {
  //       PyErr_SetString(PyExc_RuntimeError,
  //                   "Chi-squared returned an impossible value.");
  //       return NULL;
  //   }

}
