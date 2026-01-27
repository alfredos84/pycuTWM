/*---------------------------------------------------------------------------*/
// * This file contains functions to solve the Split-Step Fourier method (SSMF)
// * needed to calculate the electric fields evolution along the nonlinear placed
// * in a ring a cavity using the function `run_opo()`.
// * 
// * For single pass simulations use the function `run_single_pass()`. In particular, 
// * this file should be used when only three equation describes the 
// * problem, i.e., sum or difference frequency generation (SFG or DFG).
// *
// * For any specific process, please check the form of the Couple wave equations
// * in the first function called dAdz().
/*---------------------------------------------------------------------------*/



#ifndef _SOLVERCUH
#define _SOLVERCUH


// Scales a vector after Fourier transforms in the space-reciprocal space
// (CUFFT_INVERSE mode)	
__global__ void cuFFT2DScaleAllFieldsDisp( complex_t *Ap, complex_t *As, complex_t *Ai )
{	
	// real_t size = static_cast<real_t>(SIZE);
	real_t size = NT;
	uint32_t idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint32_t idy = threadIdx.y + blockDim.y*blockIdx.y;

	for (uint32_t idt = 0; idt < NT; idt++){
		if( idx < NX and idy < NY ){
			Ap[IDX(idt, idx + NX * idy)] /= size;
			As[IDX(idt, idx + NX * idy)] /= size;
			Ai[IDX(idt, idx + NX * idy)] /= size;
		}
	}
	return ;

}


// Scales a vector after Fourier transforms in the space-reciprocal space
// (CUFFT_INVERSE mode)	
__global__ void cuFFT2DScaleAllFieldsDiff( complex_t *Ap, complex_t *As, complex_t *Ai )
{	
	real_t size = static_cast<real_t>(NX*NY);

	uint32_t idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint32_t idy = threadIdx.y + blockDim.y*blockIdx.y;

	for (uint32_t idt = 0; idt < NT; idt++){
		if( idx < NX and idy < NY ){
			Ap[IDX(idx,idy,idt)] /= size;
			As[IDX(idx,idy,idt)] /= size;
			Ai[IDX(idx,idy,idt)] /= size;
		}
	}
	return ;

}


// This kernel set the linear propagator operators useful in dispersion calculations
__global__ void dispersionPropagators ( complex_t *eiLz_p, complex_t *eiLz_s, complex_t *eiLz_i, 
					real_t vp, real_t vs, real_t vi, 
					real_t b2p, real_t b2s, real_t b2i,
					real_t dz, real_t *w, uint32_t pass)
{	
			
	uint32_t idw = threadIdx.x + blockDim.x*blockIdx.x;
	
	complex_t Im; Im.x = 0; Im.y = 1;
	int dir = powf(-1, pass);

	if( idw < NT ){
		eiLz_p[idw] = CpxExp(dir*dz*w[idw]*((1/vp-1/vi)+0.5*w[idw]*b2p ));
		eiLz_s[idw] = CpxExp(dir*dz*w[idw]*((1/vs-1/vs)+0.5*w[idw]*b2s ));
		eiLz_i[idw] = CpxExp(dir*dz*w[idw]*((1/vi-1/vs)+0.5*w[idw]*b2i )); 
	}

	return ;	
}


// Reorganize 3D tensor to 2D array
__global__ void kernelTensor3DTo2D_3Fileds(	complex_t *Ap2D, complex_t *Ap3D,
					complex_t *As2D, complex_t *As3D,
					complex_t *Ai2D, complex_t *Ai3D)
{
	uint32_t idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint32_t idy = threadIdx.y + blockDim.y*blockIdx.y;

	for (uint32_t idt = 0; idt < NT; idt++){
		if( idx < NX and idy < NY ){
			Ap2D[IDX(idt, idx + NX * idy)] = Ap3D[IDX(idx, idy, idt)];
			As2D[IDX(idt, idx + NX * idy)] = As3D[IDX(idx, idy, idt)];
			Ai2D[IDX(idt, idx + NX * idy)] = Ai3D[IDX(idx, idy, idt)];
		}
	}

	return;
}


// Reorganize 2D array back to 3D tensor
__global__ void kernelTensor2DTo3D_3Fileds(	complex_t *Ap3D, complex_t *Ap2D,
					complex_t *As3D, complex_t *As2D,
					complex_t *Ai3D, complex_t *Ai2D)
{
	uint32_t idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint32_t idy = threadIdx.y + blockDim.y*blockIdx.y;

	for (uint32_t idt = 0; idt < NT; idt++){
		if( idx < NX and idy < NY ){
			Ap3D[IDX(idx, idy, idt)] = Ap2D[IDX(idt, idx + NX * idy)];
			As3D[IDX(idx, idy, idt)] = As2D[IDX(idt, idx + NX * idy)];
			Ai3D[IDX(idx, idy, idt)] = Ai2D[IDX(idt, idx + NX * idy)];
		}
	}

	return;
}


// Product of complex numbers in GPU for Dispersion propagator
__global__ void kernelDispersionPropagatorProduct(	complex_t *Awp_prop, complex_t *Awp, complex_t *eiLz_p,
							complex_t *Aws_prop, complex_t *Aws, complex_t *eiLz_s,
							complex_t *Awi_prop, complex_t *Awi, complex_t *eiLz_i )
{
	uint32_t idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint32_t idy = threadIdx.y + blockDim.y*blockIdx.y;

	for (uint32_t idt = 0; idt < NT; idt++){
		if( idx < NX and idy < NY ){
			Awp_prop[IDX(idt, idx + NX * idy)] = Awp[IDX(idt, idx + NX * idy)] * eiLz_p[IDX(idt)];
			Aws_prop[IDX(idt, idx + NX * idy)] = Aws[IDX(idt, idx + NX * idy)] * eiLz_s[IDX(idt)];
			Awi_prop[IDX(idt, idx + NX * idy)] = Awi[IDX(idt, idx + NX * idy)] * eiLz_i[IDX(idt)];
		}
    }

	return ;
}


// This kernel set the beam propagator operators useful in diffraction calculations
__global__ void kernelSetDiffractionPropagator ( complex_t *eiQz_p, complex_t *eiQz_s, complex_t *eiQz_i, 
						real_t uX, real_t uY,
						real_t dx, real_t dy, real_t dz,
						real_t kp, real_t ks, real_t ki,
						real_t ap, real_t as, real_t ai, uint32_t pass )
{	
	real_t dfX  = 1/dx/NX;	real_t dfY  = 1/dy/NY;
		
	uint32_t idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint32_t idy = threadIdx.y + blockDim.y*blockIdx.y;
	
	int dir = powf(-1, pass);
	
	if( idx < NX and idy < NY){
		eiQz_p[IDX(idx,idy,0)] = CpxExp(dir*dz*(2*powf(PI,2)/kp * ( dfX*dfX*powf(idx - uX,2) + dfY*dfY*powf(idy - uY,2))))*expf(-0.5f*ap*dz); 
		eiQz_s[IDX(idx,idy,0)] = CpxExp(dir*dz*(2*powf(PI,2)/ks * ( dfX*dfX*powf(idx - uX,2) + dfY*dfY*powf(idy - uY,2))))*expf(-0.5f*as*dz); 
		eiQz_i[IDX(idx,idy,0)] = CpxExp(dir*dz*(2*powf(PI,2)/ki * ( dfX*dfX*powf(idx - uX,2) + dfY*dfY*powf(idy - uY,2))))*expf(-0.5f*ai*dz); 
	}
	
	
	return ;	
}


// Product of complex numbers in GPU for Diffraction propagator
__global__ void kernelDiffractionPropagatorProduct (complex_t *AQp_propagated, complex_t *eiQz_p, complex_t *AQp,
						complex_t *AQs_propagated, complex_t *eiQz_s, complex_t *AQs,
						complex_t *AQi_propagated, complex_t *eiQz_i, complex_t *AQi )
{	
	uint32_t idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint32_t idy = threadIdx.y + blockDim.y*blockIdx.y;

	for (uint32_t idt = 0; idt < NT; idt++){
		if( idx < NX and idy < NY ){
			AQp_propagated[IDX(idx,idy,idt)] = eiQz_p[IDX(idx,idy,0)] * AQp[IDX(idx,idy,idt)] ;
			AQs_propagated[IDX(idx,idy,idt)] = eiQz_s[IDX(idx,idy,0)] * AQs[IDX(idx,idy,idt)] ;
			AQi_propagated[IDX(idx,idy,idt)] = eiQz_i[IDX(idx,idy,0)] * AQi[IDX(idx,idy,idt)] ;
		}
	}

	return ;
}


/** Computes the nonlinear part: dA/dz=i.κ.Ax.Ay.exp(i.Δk.L) and saves the result in dAx (x represents the different fields) */
__global__ void dAdz( complex_t *dAp, complex_t *dAs,  complex_t *dAi, complex_t *Ap, complex_t *As, complex_t *Ai, 
	real_t kappa_p, real_t kappa_s, real_t kappa_i, real_t dk, real_t z, uint32_t pass )
{
	complex_t Im; Im.x = 0; Im.y = 1;
	real_t dir = powf(-1, pass);

	uint32_t idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint32_t idy = threadIdx.y + blockDim.y*blockIdx.y;

	for (uint32_t idt = 0; idt < NT; idt++){
		if( idx < NX and idy < NY ){
			
			#ifdef OPG
			dAp[IDX(idx,idy,idt)]  = Im * kappa_p * As[IDX(idx,idy,idt)] * Ai[IDX(idx,idy,idt)] * CpxExp(-dk*z) ;
			dAs[IDX(idx,idy,idt)]  = Im * kappa_s * Ap[IDX(idx,idy,idt)] * CpxConj(Ai[IDX(idx,idy,idt)]) * CpxExp(+dk*z);
			dAi[IDX(idx,idy,idt)]  = Im * kappa_i * Ap[IDX(idx,idy,idt)] * CpxConj(As[IDX(idx,idy,idt)]) * CpxExp(+dk*z);
			#endif

			#ifdef SFG
			dAp[IDX(idx,idy,idt)]  = Im * kappa_p * As[IDX(idx,idy,idt)] * Ai[IDX(idx,idy,idt)] * CpxExp(-dk*z) ;
			dAs[IDX(idx,idy,idt)]  = Im * kappa_s * Ap[IDX(idx,idy,idt)] * CpxConj(Ai[IDX(idx,idy,idt)]) * CpxExp(+dk*z);
			dAi[IDX(idx,idy,idt)]  = Im * kappa_i * Ap[IDX(idx,idy,idt)] * CpxConj(As[IDX(idx,idy,idt)]) * CpxExp(+dk*z);
			#endif

			#ifdef SHG
			dAp[IDX(idx,idy,idt)]  = 1.0f*Im * kappa_p * As[IDX(idx,idy,idt)] * CpxConj(Ap[IDX(idx,idy,idt)]) * CpxExp(+dir*dk*z) ;
			dAs[IDX(idx,idy,idt)]  = 0.5f*Im * kappa_s * Ap[IDX(idx,idy,idt)] * Ap[IDX(idx,idy,idt)] * CpxExp(-dir*dk*z);
			dAi[IDX(idx,idy,idt)]  = 0.5f*Im * kappa_s * Ap[IDX(idx,idy,idt)] * Ap[IDX(idx,idy,idt)] * CpxExp(-dir*dk*z);
			#endif			

		}
	}
	
	return ;
}


/** Computes a linear combination Ax + s.kx and saves the result in aux_x */
__global__ void kernelLinealCombination(complex_t *auxp, complex_t *auxs, complex_t *auxi,
					complex_t *Ap, complex_t *As, complex_t *Ai, 
					complex_t *kp, complex_t *ks, complex_t *ki, real_t s )
{
	uint32_t idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint32_t idy = threadIdx.y + blockDim.y*blockIdx.y;

	for (uint32_t idt = 0; idt < NT; idt++){
		if( idx < NX and idy < NY ){	
			auxp[IDX(idx,idy,idt)] = Ap[IDX(idx,idy,idt)] + kp[IDX(idx,idy,idt)] * s;
			auxs[IDX(idx,idy,idt)] = As[IDX(idx,idy,idt)] + ks[IDX(idx,idy,idt)] * s;
			auxi[IDX(idx,idy,idt)] = Ai[IDX(idx,idy,idt)] + ki[IDX(idx,idy,idt)] * s;
		}
	}

	return ;
}


/** This kernel computes the final sum after appling the Rounge-Kutta algorithm */
__global__ void rk4(complex_t *Ap, complex_t *As,  complex_t *Ai, 
		complex_t *k1p, complex_t *k1s, complex_t *k1i, 
		complex_t *k2p, complex_t *k2s, complex_t *k2i, 
		complex_t *k3p, complex_t *k3s, complex_t *k3i,
		complex_t *k4p, complex_t *k4s, complex_t *k4i, real_t dz )
{
	
	uint32_t idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint32_t idy = threadIdx.y + blockDim.y*blockIdx.y;

	for (uint32_t idt = 0; idt < NT; idt++){
		if( idx < NX and idy < NY ){
			Ap[IDX(idx,idy,idt)] += (k1p[IDX(idx,idy,idt)] + 2.0f*k2p[IDX(idx,idy,idt)] + 2.0f*k3p[IDX(idx,idy,idt)] + k4p[IDX(idx,idy,idt)]) * dz / 6.0f;
			As[IDX(idx,idy,idt)] += (k1s[IDX(idx,idy,idt)] + 2.0f*k2s[IDX(idx,idy,idt)] + 2.0f*k3s[IDX(idx,idy,idt)] + k4s[IDX(idx,idy,idt)]) * dz / 6.0f;
			Ai[IDX(idx,idy,idt)] += (k1i[IDX(idx,idy,idt)] + 2.0f*k2i[IDX(idx,idy,idt)] + 2.0f*k3i[IDX(idx,idy,idt)] + k4i[IDX(idx,idy,idt)]) * dz / 6.0f;
		}
	}
	
	return ;
}



// Difine the class Solver

class Solver
{	
public:	

	cVecd_t k1p, k2p, k3p, k4p;
	cVecd_t k1s, k2s, k3s, k4s;
	cVecd_t k1i, k2i, k3i, k4i;
	cVecd_t auxp, auxs, auxi;
	uint32_t npasses;
	rVech_t mirror_distances;
	bool add_phase_air;

	Crystal *Cr;	EFields *A; json config;

	Solver(Crystal *_Cr, EFields *_A, json _config) : Cr(_Cr), A(_A), config(_config)
	{	// Constructor

		k1p.resize(SIZE); k2p.resize(SIZE); k3p.resize(SIZE); k4p.resize(SIZE);
		k1s.resize(SIZE); k2s.resize(SIZE); k3s.resize(SIZE); k4s.resize(SIZE);
		k1i.resize(SIZE); k2i.resize(SIZE); k3i.resize(SIZE); k4i.resize(SIZE);
		auxp.resize(SIZE); auxs.resize(SIZE); auxi.resize(SIZE);

        // Pointers for dispersion()
        

		print_line_on_screen();
		printf("\nInitialize Solver...\n");
	}

	
	~Solver(){ printf("Solver finished.\n"); }

	// Methods definition	
	bool check_courant_stability(real_t dz , real_t kpdxdy2);
	void set_disp_propagators(uint32_t pass = 1);
	void dispersion();
	void set_diff_propagators(uint32_t pass = 1);
	void diffraction();
	void solver_rk4( real_t z, uint32_t pass = 1);
	void ssmf (uint32_t pass);
	void set_npasses(uint32_t _npasses);
	void set_mirror_distances(rVech_t _mirror_distances);
	void set_air_phase();
	void run_single_pass();
	void run_multipass();

};


// Methods declaration

bool Solver::check_courant_stability(real_t dz , real_t kpdxdy)
{
	bool condition;
	
	if(dz < kpdxdy){condition = true;}
	else{condition = false;}

	return condition;
}


void Solver::set_npasses(uint32_t _npasses)
{
	this->npasses = _npasses;
	return ;
}


void Solver::set_mirror_distances(rVech_t _mirror_distances)
{
	this->mirror_distances = _mirror_distances;
	for(int i = 0; i < mirror_distances.size(); i++){
		std::cout << "        ---> Set mirror distances:" << std::endl;	
		std::cout << "        ---> Mirror distances = " << mirror_distances[i] << " cm" << std::endl;
	}

	return;
}


void Solver::set_air_phase()
{
	this->add_phase_air = true;
	return;
}


// Set vectors for dispersion propagators
void Solver::set_disp_propagators( uint32_t pass )
{	
	// Parameters for kernels 1D
	dim3 block1D(BLKT);	dim3 grid1D((NT+BLKT-1)/BLKT);
	
	real_t *w_ptr = thrust::raw_pointer_cast(A->w.data());
	complex_t *eiLz_p_ptr = thrust::raw_pointer_cast(A->eiLz_p.data());
	complex_t *eiLz_s_ptr = thrust::raw_pointer_cast(A->eiLz_s.data());
	complex_t *eiLz_i_ptr = thrust::raw_pointer_cast(A->eiLz_i.data());

	dispersionPropagators<<<grid1D, block1D>>> (eiLz_p_ptr, eiLz_s_ptr, eiLz_i_ptr, 
												A->vp, A->vs, A->vi, A->b2p, A->b2s, A->b2i,
												A->dz, w_ptr, pass );
	CHECK(cudaDeviceSynchronize());
	
	return ;
}


// Applies the dispersion term to the electric fields
void Solver::dispersion ( )
{
	// Parameters for kernels: 3D-tensors converted into 2D
	dim3 block2D(BLKX, BLKY);
	dim3 grid2D((NX+BLKX-1)/BLKX, (NY+BLKY-1)/BLKY);

	complex_t *Ap_ptr = thrust::raw_pointer_cast(A->Ap.data());
	complex_t *As_ptr = thrust::raw_pointer_cast(A->As.data());
	complex_t *Ai_ptr = thrust::raw_pointer_cast(A->Ai.data());
	complex_t *eiLz_p_ptr = thrust::raw_pointer_cast(A->eiLz_p.data());
	complex_t *eiLz_s_ptr = thrust::raw_pointer_cast(A->eiLz_s.data());
	complex_t *eiLz_i_ptr = thrust::raw_pointer_cast(A->eiLz_i.data());
	
    complex_t *Aux_Ap_ptr = thrust::raw_pointer_cast(A->Aux_Ap.data());
	complex_t *Aux_As_ptr = thrust::raw_pointer_cast(A->Aux_As.data());
	complex_t *Aux_Ai_ptr = thrust::raw_pointer_cast(A->Aux_Ai.data());
    complex_t *Aux_Awp_ptr = thrust::raw_pointer_cast(A->Aux_Awp.data());
	complex_t *Aux_Aws_ptr = thrust::raw_pointer_cast(A->Aux_Aws.data());
	complex_t *Aux_Awi_ptr = thrust::raw_pointer_cast(A->Aux_Awi.data());    
    complex_t *Aux_Awp_prop_ptr = thrust::raw_pointer_cast(A->Aux_Awp_prop.data());
	complex_t *Aux_Aws_prop_ptr = thrust::raw_pointer_cast(A->Aux_Aws_prop.data());
	complex_t *Aux_Awi_prop_ptr = thrust::raw_pointer_cast(A->Aux_Awi_prop.data());        

	// Convert into 2D tensor to paralelize FFTs
	kernelTensor3DTo2D_3Fileds<<<grid2D,block2D>>>( Aux_Ap_ptr, Ap_ptr,
        Aux_As_ptr, As_ptr,
        Aux_Ai_ptr, Ai_ptr );

	// CHECK(cudaDeviceSynchronize());

	cufftExecC2C(A->planDispersion, (complex_t *)Aux_Ap_ptr, (complex_t *)Aux_Awp_ptr, CUFFT_INVERSE);
	// CHECK(cudaDeviceSynchronize());
	cufftExecC2C(A->planDispersion, (complex_t *)Aux_As_ptr, (complex_t *)Aux_Aws_ptr, CUFFT_INVERSE);
	// CHECK(cudaDeviceSynchronize());
	cufftExecC2C(A->planDispersion, (complex_t *)Aux_Ai_ptr, (complex_t *)Aux_Awi_ptr, CUFFT_INVERSE);
	// CHECK(cudaDeviceSynchronize());

	cuFFT2DScaleAllFieldsDisp<<<grid2D,block2D>>>(	Aux_Awp_ptr, Aux_Aws_ptr, Aux_Awi_ptr);
	// CHECK(cudaDeviceSynchronize());

	kernelDispersionPropagatorProduct<<<grid2D,block2D>>>(	Aux_Awp_prop_ptr, Aux_Awp_ptr, eiLz_p_ptr,
															Aux_Aws_prop_ptr, Aux_Aws_ptr, eiLz_s_ptr,
															Aux_Awi_prop_ptr, Aux_Awi_ptr, eiLz_i_ptr );
	// CHECK(cudaDeviceSynchronize());

	cufftExecC2C(A->planDispersion, (complex_t *)Aux_Awp_prop_ptr, (complex_t *)Aux_Ap_ptr, CUFFT_FORWARD);
	// CHECK(cudaDeviceSynchronize());
	cufftExecC2C(A->planDispersion, (complex_t *)Aux_Aws_prop_ptr, (complex_t *)Aux_As_ptr, CUFFT_FORWARD);
	// CHECK(cudaDeviceSynchronize());
	cufftExecC2C(A->planDispersion, (complex_t *)Aux_Awi_prop_ptr, (complex_t *)Aux_Ai_ptr, CUFFT_FORWARD);
	// CHECK(cudaDeviceSynchronize());

	// Convert into 3D tensor
	kernelTensor2DTo3D_3Fileds<<<grid2D,block2D>>>( Ap_ptr, Aux_Ap_ptr,
        As_ptr, Aux_As_ptr,
        Ai_ptr, Aux_Ai_ptr );
	// CHECK(cudaDeviceSynchronize());
	
	return ;
}


// Set vectors for diffraction propagators
void Solver::set_diff_propagators(uint32_t pass )
{	
	
	// Parameters for kernels 2D
	dim3 block2D(BLKX, BLKY);	dim3 grid2D((NX+BLKX-1)/BLKX, (NY+BLKY-1)/BLKY);
	// 		* Diffraction operation ∇²_{xy}: Beam propagator for Pump and Signal
	
	complex_t *eiQz_p_ptr = thrust::raw_pointer_cast(A->eiQz_p.data());
	complex_t *eiQz_s_ptr = thrust::raw_pointer_cast(A->eiQz_s.data());
	complex_t *eiQz_i_ptr = thrust::raw_pointer_cast(A->eiQz_i.data());

	kernelSetDiffractionPropagator<<<grid2D, block2D>>> (eiQz_p_ptr, eiQz_s_ptr, eiQz_i_ptr, 
												0.5f*real_t(NX),  0.5f*real_t(NY), 
												A->dx, A->dy, A->dz,
												A->kp, A->ks, A->ki,
												Cr->alpha_crp, Cr->alpha_crs, Cr->alpha_cri, pass);
	CHECK(cudaDeviceSynchronize());
	
	A->fftShift2D ( A->eiQz_p );	A->fftShift2D ( A->eiQz_s ); A->fftShift2D ( A->eiQz_i );

	return ;
}


// Applies the diffraction term to the electric fields
void Solver::diffraction ()
{	
	// Parameters for kernels 2D
	dim3 block2D(BLKX, BLKY);	dim3 grid2D((NX+BLKX-1)/BLKX, (NY+BLKY-1)/BLKY);

	complex_t *AuxQp_ptr = thrust::raw_pointer_cast(A->AuxQp.data());
    complex_t *Ap_ptr = thrust::raw_pointer_cast(A->Ap.data());
	complex_t *AQp_ptr = thrust::raw_pointer_cast(A->AQp.data());
	complex_t *eiQz_p_ptr = thrust::raw_pointer_cast(A->eiQz_p.data());

	complex_t *AuxQs_ptr = thrust::raw_pointer_cast(A->AuxQs.data());
	complex_t *As_ptr = thrust::raw_pointer_cast(A->As.data());
	complex_t *AQs_ptr = thrust::raw_pointer_cast(A->AQs.data());
	complex_t *eiQz_s_ptr = thrust::raw_pointer_cast(A->eiQz_s.data());

	complex_t *AuxQi_ptr = thrust::raw_pointer_cast(A->AuxQi.data());
	complex_t *Ai_ptr = thrust::raw_pointer_cast(A->Ai.data());
	complex_t *AQi_ptr = thrust::raw_pointer_cast(A->AQi.data());
	complex_t *eiQz_i_ptr = thrust::raw_pointer_cast(A->eiQz_i.data());


	cufftExecC2C(A->planDiffraction, (complex_t *)Ap_ptr, (complex_t *)AQp_ptr, CUFFT_FORWARD);
	// CHECK(cudaDeviceSynchronize());
	cufftExecC2C(A->planDiffraction, (complex_t *)As_ptr, (complex_t *)AQs_ptr, CUFFT_FORWARD);
	// CHECK(cudaDeviceSynchronize());
	cufftExecC2C(A->planDiffraction, (complex_t *)Ai_ptr, (complex_t *)AQi_ptr, CUFFT_FORWARD);
	// CHECK(cudaDeviceSynchronize());
	
	kernelDiffractionPropagatorProduct<<<grid2D, block2D>>> (	AuxQp_ptr, eiQz_p_ptr, AQp_ptr,
																AuxQs_ptr, eiQz_s_ptr, AQs_ptr,	
																AuxQi_ptr, eiQz_i_ptr, AQi_ptr	);
	// CHECK(cudaDeviceSynchronize());
	
	CHECK(cudaMemcpy(AQp_ptr, AuxQp_ptr, nBytes3Dc, cudaMemcpyDeviceToDevice));	
	CHECK(cudaMemcpy(AQs_ptr, AuxQs_ptr, nBytes3Dc, cudaMemcpyDeviceToDevice));
	CHECK(cudaMemcpy(AQi_ptr, AuxQi_ptr, nBytes3Dc, cudaMemcpyDeviceToDevice));
	
	cufftExecC2C(A->planDiffraction, (complex_t *)AQp_ptr, (complex_t *)Ap_ptr, CUFFT_INVERSE);
	// CHECK(cudaDeviceSynchronize());
	cufftExecC2C(A->planDiffraction, (complex_t *)AQs_ptr, (complex_t *)As_ptr, CUFFT_INVERSE);
	// CHECK(cudaDeviceSynchronize());
	cufftExecC2C(A->planDiffraction, (complex_t *)AQi_ptr, (complex_t *)Ai_ptr, CUFFT_INVERSE);
	// CHECK(cudaDeviceSynchronize());

	cuFFT2DScaleAllFieldsDiff<<<grid2D, block2D>>>(Ap_ptr, As_ptr, Ai_ptr);
	// CHECK(cudaDeviceSynchronize());		
		
	return ;
}


void Solver::solver_rk4( real_t z, uint32_t pass)
{	
	// A->Applies the Fourh-order Runge-Kutta Method with fixed step size dz
	// This function apply the fourth-order Runge-Kutta method	
	
	// Parameters for kernels
	dim3 block2D(BLKX, BLKY);	dim3 grid2D((NX+BLKX-1)/BLKX, (NY+BLKY-1)/BLKY);


	// Define pointers to use them in kernels
	complex_t * k1p_ptr = thrust::raw_pointer_cast(this->k1p.data());
	complex_t * k2p_ptr = thrust::raw_pointer_cast(this->k2p.data());
	complex_t * k3p_ptr = thrust::raw_pointer_cast(this->k3p.data());
	complex_t * k4p_ptr = thrust::raw_pointer_cast(this->k4p.data());
	
	complex_t * k1s_ptr = thrust::raw_pointer_cast(this->k1s.data());
	complex_t * k2s_ptr = thrust::raw_pointer_cast(this->k2s.data());
	complex_t * k3s_ptr = thrust::raw_pointer_cast(this->k3s.data());
	complex_t * k4s_ptr = thrust::raw_pointer_cast(this->k4s.data());

	complex_t * k1i_ptr = thrust::raw_pointer_cast(this->k1i.data());
	complex_t * k2i_ptr = thrust::raw_pointer_cast(this->k2i.data());
	complex_t * k3i_ptr = thrust::raw_pointer_cast(this->k3i.data());
	complex_t * k4i_ptr = thrust::raw_pointer_cast(this->k4i.data());
	
	complex_t * Ap_ptr  = thrust::raw_pointer_cast(A->Ap.data());
	complex_t * As_ptr  = thrust::raw_pointer_cast(A->As.data());
	complex_t * Ai_ptr  = thrust::raw_pointer_cast(A->Ai.data());

	complex_t * auxp_ptr = thrust::raw_pointer_cast(this->auxp.data());
	complex_t * auxs_ptr = thrust::raw_pointer_cast(this->auxs.data());
	complex_t * auxi_ptr = thrust::raw_pointer_cast(this->auxi.data());
	
	real_t dz = Cr->dz;
	real_t dk = A->dk; 
	real_t kp = A->kappa_p, ks = A->kappa_s, ki = A->kappa_i;

	//k1 = dAdz(kappas,dk,z,A)
	dAdz<<<grid2D,block2D>>>( k1p_ptr, k1s_ptr, k1i_ptr, Ap_ptr, As_ptr, Ai_ptr, kp, ks, ki, dk, z, pass );
	// CHECK(cudaDeviceSynchronize()); 

	//k2 = dAdz(kappas,dk,z+dz/2,A+k1/2) -> aux = A+k1/2
	kernelLinealCombination<<<grid2D,block2D>>>( auxp_ptr, auxs_ptr, auxi_ptr, Ap_ptr, As_ptr, Ai_ptr, k1p_ptr, k1s_ptr, k1i_ptr, 0.5f );
	// CHECK(cudaDeviceSynchronize());   
	dAdz<<<grid2D,block2D>>>( k2p_ptr, k2s_ptr, k2i_ptr, auxp_ptr, auxs_ptr, auxi_ptr, kp, ks, ki, dk, z+dz/4.0f, pass );
	// CHECK(cudaDeviceSynchronize());

	// k3 = dAdz(kappas,dk,z+dz/2,A+k2/2)
	kernelLinealCombination<<<grid2D,block2D>>>( auxp_ptr, auxs_ptr, auxi_ptr, Ap_ptr, As_ptr, Ai_ptr, k2p_ptr, k2s_ptr, k2i_ptr, 0.5f );
	// CHECK(cudaDeviceSynchronize());   
	dAdz<<<grid2D,block2D>>>( k3p_ptr, k3s_ptr, k3i_ptr, auxp_ptr, auxs_ptr, auxi_ptr, kp, ks, ki, dk, z+dz/4.0f, pass );
	// CHECK(cudaDeviceSynchronize());

	// k4 = dAdz(kappas,dk,z+dz,A+k3)
	kernelLinealCombination<<<grid2D,block2D>>>( auxp_ptr, auxs_ptr, auxi_ptr, Ap_ptr, As_ptr, Ai_ptr, k3p_ptr, k3s_ptr, k3i_ptr, 1.0f );
	// CHECK(cudaDeviceSynchronize());   
	dAdz<<<grid2D,block2D>>>( k4p_ptr, k4s_ptr, k4i_ptr, auxp_ptr, auxs_ptr, auxi_ptr, kp, ks, ki, dk, z+dz/2.0f, pass );
	// CHECK(cudaDeviceSynchronize());

	// A = A + (k1 + 2*k2 + 2*k3 + k4)/6
	rk4<<<grid2D,block2D>>>(Ap_ptr, As_ptr, Ai_ptr, 
							k1p_ptr, k1s_ptr, k1i_ptr,
							k2p_ptr, k2s_ptr, k2i_ptr, 
							k3p_ptr, k3s_ptr, k3i_ptr,
							k4p_ptr, k4s_ptr, k4i_ptr, 
							dz/2.0f );
	// CHECK(cudaDeviceSynchronize());

	return ;
	
}


void Solver::ssmf (uint32_t pass )
{
	bool save_inside = this->config["save_mode"]["save_full_field_inside_crystal"].get<bool>();
	bool save_p      = this->config["save_mode"]["save_pump"].get<bool>();
	bool save_s      = this->config["save_mode"]["save_signal"].get<bool>();
	bool save_i      = this->config["save_mode"]["save_idler"].get<bool>();

	real_t z = 0;
	for (uint32_t sz = 0; sz < NZ; sz++)
	{	
		solver_rk4(z, pass); // RK4 in dz/2
		#ifdef DISPERSION
		dispersion(); // Dispersion in dz
		// Show status
    	static int last_percent = -1; // Solo inicializa una vez
    	int percent = static_cast<int>(sz * 100.0 / (NZ - 1) + 0.5);
    	if (percent % 10 == 0 and percent != last_percent) {
        	std::cout << "\r        ---> Completed: " << percent << " %" << std::flush;
        	last_percent = percent;
        	if (percent == 100) std::cout << std::endl;
    	}
		#endif
		diffraction(); // Diffraction in dz
		solver_rk4(z, pass); // RK4 in dz/2
		z+=(Cr->dz);
		#ifndef DISPERSION
		if (save_inside && save_p)
    		write_zslice_h5(A->Ap, "pump_inside_crystal.h5", sz + NZ*(pass - 1));
		if (save_inside && save_s)
    		write_zslice_h5(A->As, "signal_inside_crystal.h5", sz + NZ*(pass - 1));
		if (save_inside && save_i)
    		write_zslice_h5(A->Ai, "idler_inside_crystal.h5", sz + NZ*(pass - 1));			
		#endif
	}

	return ;
}


void Solver::run_single_pass()
{
	if (check_courant_stability((Cr->dz) , 0.2f*(A->kp)*(Cr->dx)*(Cr->dy)))
	{
		bool save_inside = this->config["save_mode"]["save_full_field_inside_crystal"].get<bool>();
    	bool save_p      = this->config["save_mode"]["save_pump"].get<bool>();
	    bool save_s      = this->config["save_mode"]["save_signal"].get<bool>();
	    bool save_i      = this->config["save_mode"]["save_idler"].get<bool>();
		
		if (save_inside){
			if(save_p)
				init_full_field_h5("pump_inside_crystal.h5", NZ) ;
			if(save_s)
				init_full_field_h5("signal_inside_crystal.h5", NZ) ;
			if(save_i)
				init_full_field_h5("idler_inside_crystal.h5", NZ) ;
		}

		double iStart = seconds();	// Timing code
        A->set_plan_diffraction();
        A->set_plan_dispersion();
		A->Ap = A->Api;
		#ifdef DISPERSION
		set_disp_propagators();
		#endif
		set_diff_propagators();
		ssmf(1);
        A->destroy_cufft_plans();
		double iElaps = seconds() - iStart;	// finish timing
		timing_code( iElaps); // print time
		log_simulation_timing(config["save_mode"]["time_performance_filename"].get<std::string>(), iElaps);
       
	}
	
	else {
		std::cout << "        ---> Courant convergence criteria: failed." << std::endl;
		std::cout << "        ---> Please change step sizes. Courant factor = " << (0.2f*A->kp*Cr->dx*Cr->dy)/Cr->dz << " < 1" << std::endl;
	}

	return ;
}


void Solver::run_multipass()
{
	if (check_courant_stability((Cr->dz) , 0.2f*(A->kp)*(Cr->dx)*(Cr->dy)))
	{
		bool save_inside = this->config["save_mode"]["save_full_field_inside_crystal"].get<bool>();
    	bool save_p      = this->config["save_mode"]["save_pump"].get<bool>();
	    bool save_s      = this->config["save_mode"]["save_signal"].get<bool>();
	    bool save_i      = this->config["save_mode"]["save_idler"].get<bool>();
		
		if (save_inside){
			if(save_p)
				init_full_field_h5("pump_inside_crystal.h5", NZ*this->npasses) ;
			if(save_s)
				init_full_field_h5("signal_inside_crystal.h5", NZ*this->npasses) ;
			if(save_i)
				init_full_field_h5("idler_inside_crystal.h5", NZ*this->npasses) ;
		}
	
		double iStart = seconds();	// Timing code
        A->set_plan_diffraction();
        A->set_plan_dispersion();
		A->Ap = A->Api;
		for (int pass = 1; pass <= this->npasses; pass++){
			std::cout << "        ---> Computing pass #" << pass << std::endl;	
			#ifdef DISPERSION
			set_disp_propagators(pass);
			#endif
			set_diff_propagators(pass);
			ssmf(pass);
			if (pass < npasses){
				A->add_phase_in_air(A->lp, A->ls, this->mirror_distances[pass-1]);
				std::cout << "        ---> Phase added to fields." << std::endl;
			}
		}
		A->destroy_cufft_plans();
		double iElaps = seconds() - iStart;	// finish timing
		timing_code( iElaps); // print time
		log_simulation_timing(config["save_mode"]["time_performance_filename"].get<std::string>(), iElaps);
	}
	
	else {
		std::cout << "        ---> Courant convergence criteria: failed." << std::endl;
		std::cout << "        ---> Please change step sizes. Courant factor = " << (0.2f*A->kp*Cr->dx*Cr->dy)/Cr->dz << " < 1" << std::endl;
	}
	return ;
}


#endif // -> #ifdef _SOLVERCUH