/*---------------------------------------------------------------------------*/
// * This file contains the class Twm which models the three-wave mixing procces
// * in a nonlinear Cr
/*---------------------------------------------------------------------------*/


#ifndef _EFIELDSCUH
#define _EFIELDSCUH

#pragma once

///////////////////////////////////////////////////////////////////////////////////////////////////
// Scales a vector after Fourier transforms in the time-frequency domaing
// (CUFFT_INVERSE mode)
__global__ void cuFFT1D_Scale(complex_t *A)
{		
	
	real_t size = static_cast<real_t>(NT);
	
	uint32_t idt = threadIdx.x + blockDim.x*blockIdx.x;
			
	if( idt < SIZE ) {A[idt] = A[idt] / size;}
	
	return ;
	
}


// CUDA kernel to add accumulated phase to both fields
__global__ void addPhase_kernel( complex_t *Ap_ptr, complex_t *As_ptr,
				complex_t *auxp_ptr, complex_t *auxs_ptr,
	 			real_t phase_p, real_t phase_s)
{
    uint32_t idx = threadIdx.x + blockDim.x * blockIdx.x;
    uint32_t idy = threadIdx.y + blockDim.y * blockIdx.y;
    uint32_t idz = 0;

    if (idx < NX && idy < NY) {
        // Add phase for pump and signal using complex exponential
        Ap_ptr[IDX(idx, idy, idz)] = CpxExp(phase_p) * auxp_ptr[IDX(idx, idy, idz)];
        As_ptr[IDX(idx, idy, idz)] = CpxExp(phase_s) * auxs_ptr[IDX(idx, idy, idz)];
    }
    return;
}


// Copy a slice from a 3D tensor to a 2D matrix of complex numbers in GPU
__global__ void kernelGetSlice ( complex_t *Aux_Ax, complex_t *Ax_ptr, int slice )
{	
	uint32_t idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint32_t idy = threadIdx.y + blockDim.y*blockIdx.y;
		
	if( idx < NX and idy < NY ){
		Aux_Ax[IDX(idx,idy,0)] = Ax_ptr[IDX(idx,idy,slice)] ;
	}
	
	return ;
}


// Copy from a 2D matrix to a specific slice of a 3D tensor of complex numbers in GPU
__global__ void kernelSetSlice ( complex_t *Ax_ptr, complex_t *Aux_Ax, int slice )
{		
	uint32_t idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint32_t idy = threadIdx.y + blockDim.y*blockIdx.y;
		
	if( idx < NX and idy < NY ){
		Ax_ptr[IDX(idx,idy,slice)] = Aux_Ax[IDX(idx,idy,0)] ;
	}
	
	return ;
}


// Flips a vector for Fourier transforms
template<typename T>
void fftshift( T& V_flip, T V )
{
	int i, c = V.size()/2;
	for ( i = 0; i < V.size()/2; i++ ){
		V_flip[i+c] = V[i];
		V_flip[i]   = V[i+c];
	}
	
	return ;
}


// Set initial pump as a Gaussian beam in time and plane in XY cordinates
__global__ void setWavePlaneCW( complex_t *Ap_ptr, real_t Ap0, real_t waist, real_t uX, real_t uY, real_t dx, real_t dy )
{		
	uint32_t idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint32_t idy = threadIdx.y + blockDim.y*blockIdx.y;

	for (uint32_t idt = 0; idt < NT; idt++){
		if( idx < NX and idy < NY ){			
			Ap_ptr[IDX(idx,idy,idt)].x = (Ap0) * expf(((-powf((idx-uX)*dx,2)-powf((idy-uY)*dy,2))/(waist*waist)));
			Ap_ptr[IDX(idx,idy,idt)].y = 0.0f;
		}
	}
	
	return;
}


// Set initial pump as a Gaussian beam in time and plane in XY cordinates
__global__ void setWavePlanePulsed( complex_t *Ap_ptr, real_t *t_ptr, real_t Ap0, real_t waist, real_t tau, real_t uX, real_t uY, real_t dx, real_t dy )
{		
	uint32_t idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint32_t idy = threadIdx.y + blockDim.y*blockIdx.y;

	for (uint32_t idt = 0; idt < NT; idt++){
		if( idx < NX and idy < NY ){			
			Ap_ptr[IDX(idx,idy,idt)].x = (Ap0) * expf( -powf(t_ptr[idt]/tau, 2) - (powf((idx-uX)*dx,2)+powf((idy-uY)*dy,2))/(waist*waist) );
			Ap_ptr[IDX(idx,idy,idt)].y = 0.0f;
		}
		
	}
	
	return;
}


// Set initial pump as a CW beam in time and XY cordinates
__global__ void setFocusedCW( complex_t *Ap_ptr, real_t Ap0, complex_t MX, real_t waist, real_t uX, real_t uY, real_t dx, real_t dy )
{	
	
	uint32_t idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint32_t idy = threadIdx.y + blockDim.y*blockIdx.y;

	for (uint32_t idt = 0; idt < NT; idt++){
		if( idx < NX and idy < NY ){			
			Ap_ptr[IDX(idx,idy,idt)] = (Ap0/MX) * CpxExp(((-powf((idx-uX)*dx,2)-powf((idy-uY)*dy,2))/(waist*waist*MX)));
		}
	}
	
	return;
}


// Set initial pump as a Gaussian beam in time and XY cordinates
__global__ void setFocusedPulsed( complex_t *Ap_ptr, real_t *t_ptr, real_t Ap0, complex_t MX, real_t waist, real_t tau, real_t uX, real_t uY, real_t dx, real_t dy )
{	
	
	uint32_t idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint32_t idy = threadIdx.y + blockDim.y*blockIdx.y;
	
	for (uint32_t idt = 0; idt < NT; idt++){
		if( idx < NX and idy < NY ){			
			Ap_ptr[IDX(idx,idy,idt)] = (Ap0/MX) * expf(-powf(t_ptr[idt]/tau, 2)) * CpxExp(((-powf((idx-uX)*dx,2)-powf((idy-uY)*dy,2))/(waist*waist*MX)));
		}
	}
	
	return;
}


// Swap horizontally the values un a matrix
__global__ void fftShift2DH( complex_t *Field, complex_t *aux)
{	 
	uint32_t idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint32_t idy = threadIdx.y + blockDim.y*blockIdx.y;

	uint32_t c = (int) floor((real_t)NX/2);
	uint32_t idt = 0;
	if (idx < c and idy < NY){
		Field[IDX(idx+c,idy,idt)]  =  aux[IDX(idx,idy,idt)];
		Field[IDX(idx,idy,idt)]    =  aux[IDX(idx+c,idy,idt)];
	}

	return ;
}


// Swap vertically the values un a matrix
__global__ void fftShift2DV( complex_t *Field, complex_t *aux)
{	// Swap vertically the values un a matrix
	uint32_t idx = threadIdx.x + blockDim.x*blockIdx.x;
	uint32_t idy = threadIdx.y + blockDim.y*blockIdx.y;
	
	uint32_t r = (int) floor((real_t)NY/2);
	uint32_t idt = 0;
	if (idy < r and idx < NX){
		Field[IDX(idx,idy+r,idt)]  =  aux[IDX(idx,idy,idt)];
		Field[IDX(idx,idy,idt)]  =  aux[IDX(idx,idy+r,idt)];
	}  
	
	return ;
}


// Class Efields
// 
class EFields
{
public:
	cVecd_t Api;
    cVecd_t Ap, As, Ai;
	cVecd_t Awp, Aws, Awi;
	cVecd_t AQp, AQs, AQi;
	cVecd_t eiQz_p, eiQz_s, eiQz_i;
	cVecd_t eiLz_p, eiLz_s, eiLz_i;
    rVecd_t t, F, w;

    // Diffraction propagators (auxiliary)
    cVecd_t AuxQp, AuxQs, AuxQi;
    // Dispersion propagators (auxiliary)
    cVecd_t Aux_Ap, Aux_As, Aux_Ai;
    cVecd_t Aux_Awp, Aux_Aws, Aux_Awi;
    cVecd_t Aux_Awp_prop, Aux_Aws_prop, Aux_Awi_prop;


	real_t lp, ls, li, waist, Power;
	real_t np, ns, ni;
	real_t vp, vs, vi;
	real_t b2p, b2s, b2i;
	real_t alpha_crp, alpha_crs, alpha_cri;
	real_t kp, ks, ki; 
	real_t kappa_p, kappa_s, kappa_i;
	real_t dx, dy, dz;
	real_t dk, dkp;	// mismatch and group-velocity mismatch

    // cuFFT plans for diffraction and dispersion
    cufftHandle planDiffraction, planDispersion; 
    
    // Constructor
    EFields(real_t _lp, real_t _ls, real_t _li, real_t _Power, real_t _waist, Crystal *Cr) :
			lp(_lp), ls(_ls), li(_li), Power(_Power), waist(_waist)
    {
        // Initialization or other constructor logic if needed
		this->Api.resize(SIZE); // initial pump efield
		this->Ap.resize(SIZE); this->As.resize(SIZE); this->Ai.resize(SIZE);
		this->Awp.resize(SIZE); this->Aws.resize(SIZE); this->Awi.resize(SIZE);
		this->AQp.resize(SIZE); this->AQs.resize(SIZE); this->AQi.resize(SIZE);
        this->AuxQp.resize(SIZE); this->AuxQs.resize(SIZE); this->AuxQi.resize(SIZE);
        this->Aux_Ap.resize(SIZE); this->Aux_As.resize(SIZE); this->Aux_Ai.resize(SIZE);
        this->Aux_Awp.resize(SIZE); this->Aux_Aws.resize(SIZE); this->Aux_Awi.resize(SIZE);
        this->Aux_Awp_prop.resize(SIZE); this->Aux_Aws_prop.resize(SIZE); this->Aux_Awi_prop.resize(SIZE);
		this->eiQz_p.resize(NX*NY); this->eiQz_s.resize(NX*NY); this->eiQz_i.resize(NX*NY);
		this->eiLz_p.resize(NT); this->eiLz_s.resize(NT); this->eiLz_i.resize(NT);
		
		this->t.resize(NT); this->F.resize(NT); this->w.resize(NT);

		print_line_on_screen();
		printf("\nInitialize Electric fields.\n");
		np  = Cr->np; ns = Cr->ns; ni = Cr->ni;
		vp  = Cr->vp; vs = Cr->vs; vi = Cr->vi;
		b2p  = Cr->b2p; b2s = Cr->b2s; b2i = Cr->b2i;
		alpha_crp = Cr->alpha_crp; alpha_crs = Cr->alpha_crs; alpha_cri = Cr->alpha_cri;
		dx  = Cr->dx; dy = Cr->dy; dz = Cr->dz;
		kappa_p  = 2*PI*Cr->dQ/(np*lp);   // pump   kappa [1/V]
		kappa_s  = 2*PI*Cr->dQ/(ns*ls);   // signal kappa [1/V]
		kappa_i  = 2*PI*Cr->dQ/(ni*li);   // idler  kappa [1/V]
		kp  = 2.0f*PI*ns/lp;			  // pump   wavevector
		ks  = 2.0f*PI*ns/ls; 			  // signal wavevector
		ki  = 2.0f*PI*ni/ls;			  // idler  wavevector
		dk = Cr->dk;
		dkp = 1/vp-1/vs;
    }

    // Destructor
	~EFields(){	printf("Electric fields switched off.\n"); }

	// Methods definition
	void set_pump_field( real_t Power, real_t FWHM, real_t waist, real_t focalpoint, std::string mode );
	void noise_generator ( cVecd_t& Vec );
	void set_time_freq_vectors(real_t t_window);
	void fftShift2D ( cVecd_t& propagator );
	real_t ref_index_air(real_t lambda_um);
	void add_phase_in_air(real_t lp, real_t ls, real_t distance_cm);
    void set_plan_diffraction();
    void set_plan_dispersion();
    void destroy_cufft_plans();
};


// Methods declaration

void EFields::set_pump_field( real_t Power, real_t FWHM, real_t waist, real_t focalpoint, std::string mode )
{
	complex_t Im; Im.x = 0; Im.y = 1;
	real_t tau    = FWHM*sqrtf(2)/(2*sqrtf(2*logf(2)));
	real_t Inten  = Power/(PI*waist*waist);
	real_t w02    = waist*waist;
	real_t zR     = PI*(this->np)*w02/lp;
	real_t eta    = focalpoint/zR;
	real_t Ap0    = sqrtf(4*Power/(EPS0*C*PI*(this->np)*w02));
	complex_t MX  = (1-Im*eta);

	complex_t *Ap_ptr = thrust::raw_pointer_cast(this->Api.data());
	real_t *t_ptr = thrust::raw_pointer_cast(this->t.data());
	dim3 block2D(BLKX, BLKY);	dim3 grid2D((NX+BLKX-1)/BLKX, (NY+BLKY-1)/BLKY);

	std::string mode1 = "waveplane-cw"; std::string mode2 = "waveplane-pulsed";
	std::string mode3 = "focused-cw"; std::string mode4 = "focused-pulsed";

	if (mode.compare(mode1) == 0){
		setWavePlaneCW<<<grid2D,block2D>>>( Ap_ptr, Ap0, waist, 0.5*NX, 0.5*NY, this->dx, this->dy );
		// CHECK(cudaDeviceSynchronize()); 
		std::cout << "        ---> Pump field mode: " + mode << std::endl;
		std::cout << "        ---> Pump Power: " << Power << " W" << std::endl;
		std::cout << "        ---> Beam waist: " << waist << " \u03BCm\n" << std::endl;
	}
	else if (mode.compare(mode2) == 0){
		setWavePlanePulsed<<<grid2D,block2D>>>( Ap_ptr, t_ptr, Ap0, waist, tau, 0.5*NX, 0.5*NY, this->dx, this->dy );
		// CHECK(cudaDeviceSynchronize()); 
		std::cout << "        ---> Pump field mode: " + mode << std::endl;
		std::cout << "        ---> Pump Power: " << Power << " W" << std::endl;
		std::cout << "        ---> Beam waist: " << waist << " \u03BCm" << std::endl;
		std::cout << "        ---> FWHM: " << FWHM << " ps\n" << std::endl;
	}
	else if (mode.compare(mode3) == 0){
		setFocusedCW<<<grid2D,block2D>>>( Ap_ptr, Ap0, MX, waist, 0.5*NX, 0.5*NY, this->dx, this->dy );
		// CHECK(cudaDeviceSynchronize());
		std::cout << "        ---> Pump field mode: " + mode << std::endl;
		std::cout << "        ---> Setting Pump e-field with \u03BE = " << eta << std::endl;
		std::cout << "        ---> Pump Power: " << Power << " W" << std::endl;
		std::cout << "        ---> Beam waist: " << waist << " \u03BCm" << std::endl;
		std::cout << "        ---> Focal position: " << focalpoint*1e-3 << " mm\n" << std::endl;
	}
	else if (mode.compare(mode4) == 0){
		setFocusedPulsed<<<grid2D,block2D>>>( Ap_ptr, t_ptr, Ap0, MX, waist, tau, 0.5*NX, 0.5*NY, this->dx, this->dy );
		// CHECK(cudaDeviceSynchronize()); 
		std::cout << "        ---> Pump field mode: " + mode << std::endl;
		std::cout << "        ---> Setting Pump e-field with \u03BE = " << eta << std::endl;
		std::cout << "        ---> Pump Power: " << Power << " W" << std::endl;
		std::cout << "        ---> Beam waist: " << waist << " \u03BCm" << std::endl;
		std::cout << "        ---> Focal position: " << focalpoint*1e-3 << " mm" << std::endl;
		std::cout << "        ---> FWHM: " << FWHM << " ps\n" << std::endl;
		
	}
	else{std::cout << "Invalid option for pump electric field!!\n";}

	return ;
}



void EFields::noise_generator ( cVecd_t& Vec )
{	// Noise generator for initial signal/idler vectors 
	uint32_t seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::normal_distribution<real_t> distribution(0.0,1.0e-15);
	
	cVech_t A(Vec.size());
	real_t nsx, nsy;    
	for (int i=0; i<Vec.size(); ++i) {
		nsx = distribution(generator); A[i].x = static_cast<real_t>(nsx);
		nsy = distribution(generator); A[i].y = static_cast<real_t>(nsy);
	}
	thrust::copy(A.begin(), A.end(), Vec.begin());

	return ;	
}



void EFields::set_time_freq_vectors(real_t t_window)
{
	std::cout << "        ---> Set time and frequency vectors" << std::endl;
	this->t = linspace<decltype(this->t)>( -t_window*0.5, t_window*0.5, this->t.size() );
	this->F = linspace<decltype(this->F)>( -0.5*this->F.size()/t_window, +0.5*this->F.size()/t_window, this->F.size() );
	fftshift<decltype(this->F)>(this->w, this->F) ;
	this->w *= (2.0*PI);
	
	return ;
}



void EFields::fftShift2D ( cVecd_t& propagator )
{	// Standard fftshift in 2D

	complex_t *propagator_ptr = thrust::raw_pointer_cast(propagator.data());
	complex_t *aux;	CHECK(cudaMalloc((void **)&aux, nBytes2Dc));
	
	dim3 block2D(BLKX, BLKY);	dim3 grid2D((NX+BLKX-1)/BLKX, (NY+BLKY-1)/BLKY);
	
	CHECK(cudaMemcpy(aux, propagator_ptr, nBytes2Dc, cudaMemcpyDeviceToDevice));
	fftShift2DV<<<grid2D, block2D>>>(propagator_ptr, aux);
	cudaDeviceSynchronize();
	
	CHECK(cudaMemcpy(aux, propagator_ptr, nBytes2Dc, cudaMemcpyDeviceToDevice));
	fftShift2DH<<<grid2D, block2D>>>(propagator_ptr, aux);
	cudaDeviceSynchronize();
	
	CHECK(cudaFree(aux));
	
	return ;	
}


// Function to estimate the refractive index of air for a given wavelength (in um)
real_t EFields::ref_index_air(real_t lambda_um) 
{
    // Simple formula from Ciddor 1996, valid in visible/NIR
	// https://refractiveindex.info/?shelf=other&book=air&page=Ciddor
	
    real_t s = 1.0 / (lambda_um * lambda_um);
    real_t n_minus_1 = 0.05792105 / (238.0185 - s) + 0.00167917 / (57.362 - s);
    return 1.0 + n_minus_1;
}


// Add this function to your class (or as a free function if needed)
void EFields::add_phase_in_air(real_t lp, real_t ls, real_t distance_cm)
{
    // --- Compute refractive indices for each wavelength
    real_t n_p = ref_index_air(lp);
    real_t n_s = ref_index_air(ls);

    // --- Convert distance to meters
    real_t L_m = distance_cm * 1e-2;

    // --- Convert wavelengths to meters
    real_t lambda_p_m = lp * 1e-6;
    real_t lambda_s_m = ls * 1e-6;

    // --- Calculate accumulated phase in radians
    real_t phase_p = std::fmod(-4.0 * PI * n_p * L_m / lambda_p_m, 2.0 * PI);
    real_t phase_s = std::fmod(+4.0 * PI * n_s * L_m / lambda_s_m, 2.0 * PI);
	
	// --- Setup pointers for device arrays
    cVecd_t auxp = this->Ap;
    cVecd_t auxs = this->As;

    complex_t *Ap_ptr = thrust::raw_pointer_cast(this->Ap.data());
    complex_t *As_ptr = thrust::raw_pointer_cast(this->As.data());
    complex_t *auxp_ptr = thrust::raw_pointer_cast(auxp.data());
    complex_t *auxs_ptr = thrust::raw_pointer_cast(auxs.data());

    dim3 block2D(BLKX, BLKY);
    dim3 grid2D((NX + BLKX - 1) / BLKX, (NY + BLKY - 1) / BLKY);

	addPhase_kernel<<<grid2D, block2D>>>(Ap_ptr, As_ptr, auxp_ptr, auxs_ptr, phase_p, phase_s );
    cudaDeviceSynchronize();

    return;
}


void EFields::set_plan_diffraction()
{
    // This function can be used to set up cuFFT plans needed for diffraction
	int batch = NT;	int rank = 2;
	int nRows = NY;	int nCols = NX; int n[2] = {NY, NX};
	int idist = NX * NY; int odist = NX * NY;
	int inembed[] = {NY , NX}; 	int onembed[] = {NY, NX};
	int istride = 1; int ostride = 1;

	cufftPlanMany(&this->planDiffraction, rank, n, inembed, istride,
        idist, onembed, ostride, odist, CUFFT_C2C, batch);
}


void EFields::set_plan_dispersion()
{
    // This function can be used to set up cuFFT plans needed for dispersion
	int NT_cufft = NT;
	cufftPlanMany(&this->planDispersion, 1, &NT_cufft, nullptr, 1, NT_cufft,
        nullptr, 1, NT_cufft, CUFFT_C2C, NX * NY );
}


void EFields::destroy_cufft_plans()
{
    // This function destroys cuFFT plans
    cufftDestroy(this->planDiffraction);
    cufftDestroy(this->planDispersion);
}
#endif // -> #ifdef _EFIELDSCUH