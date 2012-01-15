/* **************************************************************************

   Some helper functions to convert between Python lists and C arrays.
 Thanks to Leo 

   ************************************************************************** */

/* PYLIST_TO_VECTOR */
double *PyList_to_vector(PyObject *input, int *size)
{
    /* Check if the input is a list */
    if (!PyList_Check(input))
    {
        PyErr_SetString(PyExc_TypeError,"Input must be a list.");
        return NULL;
    }
    else
    {
        int i;
        PyObject *item;
        double *result;

        /* Get the size of the vector from the number of lines in the PyList */
        *size = PyList_Size(input);

        /* Malloc some memory for the c vector */
        result = (double *)malloc((*size)*sizeof(double));

        /* Check if each item in input is a float and put them in the array */
        for(i=0; i<*size; i++)
        {
            item = PyList_GetItem(input, i);

            /* Check is the item is a float or int */
            if(!PyNumber_Check(item))
            {
                PyErr_SetString(PyExc_TypeError,"Input must be a list of floats or ints.");
                free(result);
                return NULL;
            }

            /* Put it in the right position in the C array */
            result[i] = PyFloat_AsDouble(item);
        }

        /* Return the C array */
        return result;
    }
}

/* VECTOR_TO_PYLIST */
PyObject *vector_to_PyList(double *vec, int size)
{
    int i;
    PyObject *result;

    /* Create a list to put the vector in */
    result = PyList_New(size);

    /* Now set the elements in it with the elements of the vector */
    for(i=0; i<size; i++)
    {
        PyList_SetItem(result, i, PyFloat_FromDouble(vec[i]));
    }

    /* Return the Python List */
    return result;
}

/*****************************************************************************/

/* 
Create a pyList with 2 collums
from two double arrays
*/
PyObject *Arrays2_to_PyList(double *a, double *b, int length)
{
    int i;
    PyObject *tmp, *result;

    /* Create the new list object */
    result = PyList_New(length);

    /* Now set the elements in it with the 2 elements */
    for(i=0; i<length; i++)
    {
        /* Create a new buffer row */
        tmp = PyList_New(2);
        
        PyList_SetItem(tmp, 0, PyFloat_FromDouble(a[i]));
        PyList_SetItem(tmp, 1, PyFloat_FromDouble(b[i]));

        /* Now put the row in the $result */
        PyList_SetItem(result, i, tmp);
    }

    /* Return the Python List */
    return result;
}
