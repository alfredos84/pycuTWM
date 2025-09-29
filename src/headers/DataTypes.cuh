#ifndef _DTYPESCONSTSCUH
#define _DTYPESCONSTSCUH

 
//  * Complex data type: a set of datatypes are
//  * defined to make the code more readable.
//  *
//  * Definitions for numbers
//  * real_t    : datatype for real numbers
//  * complex_t : datatype for complex numbers
//  * 
//  * Definitions for vectors:
//  * 
//  * rVech_t   : real vector host
//  * rVecd_t   : real vector device
//  * cVech_t   : complex vector host
//  * cVecd_t   : complex vector device

using real_t = float;
using complex_t = cufftComplex;
using rVech_t = thrust::host_vector<real_t>;
using rVecd_t = thrust::device_vector<real_t>;
using cVech_t = thrust::host_vector<complex_t>;
using cVecd_t = thrust::device_vector<complex_t>;

using json	  =  nlohmann::json;

	// Define memory size for vectors and matrices
const size_t nBytes1Dr = sizeof(real_t) * NT;		    // real 1D
const size_t nBytes1Dc = sizeof(complex_t) * NT;		// complex 1D
const size_t nBytes2Dr = sizeof(real_t) * NX * NY;		// real 2D 
const size_t nBytes2Dc = sizeof(complex_t) * NX * NY;	// complex 2D
const size_t nBytes3Dr = sizeof(real_t) * SIZE;	        // real 3D
const size_t nBytes3Dc = sizeof(complex_t) * SIZE;	    // complex 3D

// Define global constnts
const real_t PI = 3.14159265358979323846;		// pi
const real_t C  = 299792458*1E6/1E12;			// speed of ligth in vacuum [um/ps]
const real_t EPS0   = 8.8541878128E-12*1E12/1E6;	// vacuum pertivity [W.ps/V²μm] 



#endif // -> #ifdef _DTYPESCONSTS