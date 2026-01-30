// * Author name: Alfredo Daniel Sanchez
// * email:       alfredo.daniel.sanchez@gmail.com


#include "headers/Libraries.cuh"		// Required libraries
#include "headers/PackageLibraries.cuh"	// Required package libraries
#include <memory>


int main(int argc, char *argv[]){

	///////////////////////////////////////////////////////////////////////////////////
	// 1. Read config.json data set. It contains all the required parameters 
	
    if (argc < 2) {
        std::cerr << "Use: " << argv[0] << " config.json" << std::endl;
        return 1;
    }

    // Read JSON file
    std::ifstream config_file(argv[1]);
    json config;
    try {
        config_file >> config;
    } catch (const nlohmann::json::parse_error& e) {
        std::cerr << "Error by passing JSON file: " << e.what() << std::endl;
        return 1;
    }

	// GPU availability check
	is_gpu_available();

	print_grid();

	///////////////////////////////////////////////////////////////////////////////////
	// 2. Build the system

	// set the nonlinear crystal
	real_t LX = config["crystal"]["dimensions"]["LX"].get<real_t>();
	real_t LY = config["crystal"]["dimensions"]["LY"].get<real_t>();
	real_t Lcr = config["crystal"]["dimensions"]["Lcr"].get<real_t>();

	real_t lp = config["crystal"]["wavelengths"]["pump"].get<real_t>();
	real_t ls = config["crystal"]["wavelengths"]["signal"].get<real_t>();
	real_t li = config["crystal"]["wavelengths"]["idler"].get<real_t>();
    std::string crystal_type = config["crystal"]["type"].get<std::string>();
    
    // Instantiate the correct crystal based on type
	// Base class for crystals: MgOsPPLT, MgOPPLN, PPLN and ZGP
    std::unique_ptr<Crystal> cr1;
	if (crystal_type == "MgO:sPPLT") {
        real_t T = config["crystal"]["properties_pp"]["temperature"].get<real_t>();
        real_t Lambda = config["crystal"]["properties_pp"]["grating_period"].get<real_t>();
        cr1 = std::make_unique<MgOsPPLT>(LX, LY, Lcr, T, Lambda, lp, ls, li);
	} else if (crystal_type == "PPLN") {
		real_t T = config["crystal"]["properties_pp"]["temperature"].get<real_t>();
		real_t Lambda = config["crystal"]["properties_pp"]["grating_period"].get<real_t>();
		cr1 = std::make_unique<PPLN>(LX, LY, Lcr, T, Lambda, lp, ls, li);
    } else if (crystal_type == "MgO:PPLN") {
		real_t T = config["crystal"]["properties_pp"]["temperature"].get<real_t>();
		real_t Lambda = config["crystal"]["properties_pp"]["grating_period"].get<real_t>();
		cr1 = std::make_unique<MgOPPLN>(LX, LY, Lcr, T, Lambda, lp, ls, li);
	} else if (crystal_type == "ZGP") {
        // Assuming polarizations are passed as an array in JSON
        // Example: "polarizations": ["e", "e", "o"]
        auto pol_p = config["crystal"]["properties_birref"]["polarization_p"].get<char>();
		auto pol_s = config["crystal"]["properties_birref"]["polarization_s"].get<char>();
		auto pol_i = config["crystal"]["properties_birref"]["polarization_i"].get<char>();
        std::tuple<char, char, char> pol = std::make_tuple(pol_p, pol_s, pol_i);
        cr1 = std::make_unique<ZGP>(LX, LY, Lcr, pol, lp, ls, li); 
    } else {
        std::cerr << "Error: Unsupported crystal type specified in JSON: " << crystal_type << std::endl;
        return ;
    }
	real_t dk = config["crystal"]["wavelengths"]["dk"].get<real_t>();
	cr1->set_dk(dk);	cr1->getCrystalProp();
	

	// b. set electric fields	
	// real_t alphas = 0.5*( (1-cav1->Rs)+cr1->alpha_crs*cr1->Lcr ), alphai = alphas; 
	real_t Power = config["fields"]["pump"]["pump_power_W"].get<real_t>();
	real_t waist = config["fields"]["pump"]["waist_um"].get<real_t>();
	real_t FWHM = config["fields"]["pump"]["fwhm"].get<real_t>();
	real_t focalpoint = (config["fields"]["pump"]["focal_point_factor"].get<real_t>()) * Lcr;
	
	
	std::unique_ptr<EFields> A = std::make_unique<EFields>(lp, ls, li, Power, waist, cr1.get());
	real_t t_window =  config["fields"]["time_freq_vect"]["time_window"].get<real_t>();
	A->set_time_freq_vectors(t_window);
	
	// set pump mode: 
	// (1) "waveplane-cw"; (2) "waveplane-pulsed"; (3) "focused-cw"; (4) "focused-pulsed"
	std::string mode = config["fields"]["pump"]["mode"].get<std::string>();
	A->set_pump_field( Power, FWHM, waist, focalpoint, mode ); A->Ap = A->Api;
	
	// saveMatrixComplex_TimeSlice ( A->Ap, 0, "Input_pump_power" );
	bool signal_input = config["fields"]["signal"]["noise_generation"].get<bool>();
	bool idler_input = config["fields"]["idler"]["noise_generation"].get<bool>();
	bool degenerate = config["fields"]["degenerate"]["enabled"].get<bool>();
	
	if (signal_input and degenerate){
		A->noise_generator( A->As ); 
		A->Ai = A->As;
	}
	else if (signal_input and idler_input) {
		A->noise_generator( A->As ); 
		A->noise_generator( A->Ai );
	} 
	else {std::cout << "        ---> Error: Check signal and idler electric field definition." << std::endl;}



	///////////////////////////////////////////////////////////////////////////////////
	// 3. run package
	
	std::unique_ptr<Solver> solv1 = std::make_unique<Solver>(cr1.get(), A.get(), config);
	
	bool multi_pass = config["mul_pass_scheme"]["multipass"].get<bool>();
	uint32_t npasses = config["mul_pass_scheme"]["npasses"].get<int>();
	bool add_phase_air = config["mul_pass_scheme"]["add_phase_air"].get<bool>();	
	std::vector<real_t> tmp = config["mul_pass_scheme"]["distances_in_air"].get<std::vector<real_t>>();
	rVech_t mirror_distances(tmp.begin(), tmp.end());
	
	if (multi_pass and (npasses > 1) ){
		solv1->set_npasses(npasses);
		solv1->set_mirror_distances(mirror_distances);
		// std::cout << "        ---> Set mirror distances:" << std::endl;
		// for (int m = 1; m < npasses; m++)
		// 	std::cout << "        ---> Mirror #" << m << " at distance " << mirror_distances[m-1] << std::endl;

		if (add_phase_air){solv1->set_air_phase();}

		if(mode == "waveplane-cw" or mode == "focused-cw"){
			std::cout << "        ---> Running Solver: multipass mode" << std::endl;
			solv1->run_multipass();
		}
		else{std::cout << "        ---> Change the pump mode to waveplane-cw or focused-cw" << std::endl;}
	}
	else if (npasses == 1 and !multi_pass){
		std::cout << "        ---> Running Solver: single-pass mode" << std::endl;
		solv1->run_single_pass();
	}
	else
		{std::cout << "        ---> Error. Number of passes must be equal or grater than 1" << std::endl;}

	///////////////////////////////////////////////////////////////////////////////////
	// 4. Save output data
	
	save_input_pump_slices_XY (A.get(), config);
	save_output_slices_XY (A.get(), config);
	save_time_and_frequency_vectors_h5(A.get(), config);

	// /////////////////////////////////////////////////////////////////////////////////
	// // 5. Delete object instances
	
	///////////////////////////////////////////////////////////////////////////////////
	
	return 0;
	
}
