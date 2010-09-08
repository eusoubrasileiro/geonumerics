%module dsprocessing

%{
#define SWIG_FILE_WITH_INIT
#include "../c/dsprocessing.h"
#include "typeconversions.c"
%}

/* temp variable to store the size of the input pyList, needed when creating a list back */
%{
unsigned int size_pyList;
%}
/* 
as its a mult type map for each block of variables it will only accept one input
In typemap for converting pyList inputs into array and also get 
the size of the pyList and call correctly the function  
$1 is the first argument of the function that will be called 
$2 is the second and so on..
$input is the input object expected here a List
*/
%typemap(in) (double *y, unsigned int N)
{
	 /* Check if is a list */
	  if (PyList_Check($input)) {
		/* convert the list to an array and also get's its size */
		$1 = PyList_to_vector($input, &size_pyList);
		$2 = size_pyList;   
	  }
	  else {
	    PyErr_SetString(PyExc_TypeError,"not a list");
	    return NULL;
	  }		
}

// %typemap(in, numinputs=0) double*
// {

// }

/* 
Out typemap for double array  back to pyList
$1 is the first result from the function called
size_pyList contains the initial list size
*/
%typemap(out) double *
{
	Py_XDECREF($result);   /* Blow away any previous result */
    $result = vector_to_PyList(result, 2);

    /* Check if there was any error in the conversion */
    if(!$result)
    {
        return NULL;
    }

}

/*
get the $1 input argument for the C code and free it
the wrapper guys used
*/ 
%typemap(freearg) (double *y, unsigned int )
{
	if($1)
    {
        free($1);
    }
}

/* 
dspLinear get a linear least squares over the input List
return the coeficients A+B*x as an list
list[0] = A 
list[1] = B 
*/
%feature("autodoc", "1");
extern double *dspLinear(double *y, unsigned int N);


/***********************************************************************/

/* 
as its a mult type map for each block of variables it will only accept one input
In typemap for converting pyList inputs into array and also get 
the size of the pyList and call correctly the function  
$1 is the first argument of the function that will be called 
$2 is the second and so on..
$input is the input object expected here a List
*/
%typemap(in) (double *x, unsigned int N)
{
	 /* Check if is a list */
	  if (PyList_Check($input)) {
		/* convert the list to an array and also get's its size */
		$1 = PyList_to_vector($input, &size_pyList);
		$2 = size_pyList;   
	  }
	  else {
	    PyErr_SetString(PyExc_TypeError,"not a list");
	    return NULL;
	  }		
}

/* 
Out typemap for double array  
result is the  result from the function called, or "result" contains that
or $1 contains that
size_pyList contains the initial list size, but doesnt matter because its always 2
so this Out TYPEMAP  here has its signature a litle diferent from the one before
the arguments names are taking in consideration if they are placed
*/
%typemap(argout) (double *x, unsigned int N)
{
	Py_XDECREF($input);   /* Blow away any previous data in the input vector*/
    $result = vector_to_PyList($1, size_pyList);

    /* Check if there was any error in the conversion */
    if(!$result)
    {
        return NULL;
    }

}


/* 
Finally dspDetrendLinear its output is not gonna be void 
Also return the coeficients A+B*x as an list
list[0] = A 
list[1] = B 
Those are not yet exposed to python.
*/
%feature("autodoc", "1");
extern double* dspDetrendLinear(double *x, unsigned int N);



/***********************************************************************/

/* 
as its a mult type map for each block of variables it will only accept one input
In typemap for converting pyList inputs into array and also get 
the size of the pyList and call correctly the function  
$1 is the first argument of the function that will be called 
$2 is the second and so on..
$input is the input object expected here a List
*/
%typemap(in) (double *x, unsigned int ns, double *a, double *b)
{
	 /* Check if is a list */
	  if (PyList_Check($input)) {
		/* convert the list to an array and also get's its size */
		$1 = PyList_to_vector($input, &size_pyList);
		$2 = size_pyList;  
		$3 = (double*) malloc(sizeof(double)*size_pyList);
		$4 = (double*) malloc(sizeof(double)*size_pyList);
	  }
	  else {
	    PyErr_SetString(PyExc_TypeError,"not a list");
	    return NULL;
	  }		
}

/* 
Out typemap for 2 double array  
result is the  result from the function called, or "result" contains that
size_pyList contains the initial list size,
*/
%typemap(argout) (double *x, unsigned int ns, double *a, double *b)
{
	/* check return status */
	if(result == -1)
		return NULL;
	else 
	{
		$result = Arrays2_to_PyList($3, $4, size_pyList);
		/* Check if there was any error in the conversion */
		if(!$result)
	    {
	        return NULL;
	    }
	}
}

%feature("autodoc", "1");
extern int dspDft(double *x, unsigned int ns, double *a, double *b);

// %pythonprepend dspDft(double *x, unsigned int ns, double **a, double **b)
// %{
	// printf("Xxx");
// %}

// %pythonappend dspDft(double*, unsigned int, double**, double**)
// %{
	// printf("Xxx");
// %}

%feature("autodoc", "1");
extern int dspFft(double *x, unsigned int ns, double *a, double *b);



/* 
Out typemap for 2 double array  
result is the  result from the function called, or "result" contains that
size_pyList contains the initial list size,
*/
%typemap(argout) (double *signal, unsigned int ns, double dt, double fc, double ramp)
{
	/* check return status */
	if(result == NULL)
		return NULL;
	else 
	{
		$result = vector_to_PyList(result, size_pyList);
		/* Check if there was any error in the conversion */
		if(!$result)
	    {
	        return NULL;
	    }
	}
}

/* 
as its a mult type map for each block of variables it will only accept one input
In typemap for converting pyList inputs into array and also get 
the size of the pyList and call correctly the function  
$1 is the first argument of the function that will be called 
$2 is the second and so on..
$input is the input object expected here a List
*/
%typemap(in) (double *signal, unsigned int ns)
{
	 /* Check if is a list */
	  if (PyList_Check($input)) {
		/* convert the list to an array and also get's its size */
		$1 = PyList_to_vector($input, &size_pyList);
		$2 = size_pyList;  		
	  }
	  else {
	    PyErr_SetString(PyExc_TypeError,"not a list");
	    return NULL;
	  }		
}

%typemap(in) double
{
	$1 = PyFloat_AsDouble($input);
}

%feature("autodoc", "1");
extern 
double* dspIIR_SincTrapezLowPassApply(double *signal, unsigned int ns, double dt, double fc, double ramp);



