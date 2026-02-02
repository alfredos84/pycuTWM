#ifndef _CFGOBJECTS
#define _CFGOBJECTS


inline std::unique_ptr<Crystal> config_crystal(const json& config)
{
    /* 
        This function configures and returns a unique pointer to a 
        Crystal object based on the provided JSON configuration.
        It supports different crystal types including ZGP, MgOsPPLT,
        PPLN, and MgOPPLN.
    */

    // set the nonlinear crystal
    real_t LX = config["crystal"]["dimensions"]["LX"].get<real_t>();
    real_t LY = config["crystal"]["dimensions"]["LY"].get<real_t>();
    real_t Lcr = config["crystal"]["dimensions"]["Lcr"].get<real_t>();

    real_t lp = config["crystal"]["wavelengths"]["pump"].get<real_t>();
    real_t ls = config["crystal"]["wavelengths"]["signal"].get<real_t>();
    real_t li = config["crystal"]["wavelengths"]["idler"].get<real_t>();
    std::string crystal_type = config["crystal"]["type"].get<std::string>();

    std::unique_ptr<Crystal> cr1;
    using MakerPP = std::unique_ptr<Crystal>(*)(real_t, real_t, real_t, real_t, real_t, real_t, real_t, real_t);
    static const std::unordered_map<std::string, MakerPP> makers_pp = {
        {"MgO:sPPLT", [](real_t LX, real_t LY, real_t Lcr, real_t T, real_t Lambda, real_t lp, real_t ls, real_t li) -> std::unique_ptr<Crystal> {
            return std::make_unique<MgOsPPLT>(LX, LY, Lcr, T, Lambda, lp, ls, li);
        }},
        {"PPLN", [](real_t LX, real_t LY, real_t Lcr, real_t T, real_t Lambda, real_t lp, real_t ls, real_t li) -> std::unique_ptr<Crystal> {
            return std::make_unique<PPLN>(LX, LY, Lcr, T, Lambda, lp, ls, li);
        }},
        {"MgO:PPLN", [](real_t LX, real_t LY, real_t Lcr, real_t T, real_t Lambda, real_t lp, real_t ls, real_t li) -> std::unique_ptr<Crystal> {
            return std::make_unique<MgOPPLN>(LX, LY, Lcr, T, Lambda, lp, ls, li);
        }},
    };

    if (crystal_type == "ZGP") {
        auto pol_p = config["crystal"]["properties_birref"]["polarization_p"].get<char>();
        auto pol_s = config["crystal"]["properties_birref"]["polarization_s"].get<char>();
        auto pol_i = config["crystal"]["properties_birref"]["polarization_i"].get<char>();
        std::tuple<char, char, char> pol = std::make_tuple(pol_p, pol_s, pol_i);
        cr1 = std::make_unique<ZGP>(LX, LY, Lcr, pol, lp, ls, li);
    } else {
        const auto& props_pp = config["crystal"]["properties_pp"];
        real_t T = props_pp["temperature"].get<real_t>();
        real_t Lambda = props_pp["grating_period"].get<real_t>();
        auto it = makers_pp.find(crystal_type);
        if (it == makers_pp.end()) {
            std::cerr << "[Error]: Unsupported crystal type specified in JSON: " << crystal_type << std::endl;
            return nullptr;
        }
        cr1 = it->second(LX, LY, Lcr, T, Lambda, lp, ls, li);
    }

    real_t dk = config["crystal"]["wavelengths"]["dk"].get<real_t>();
    cr1->set_dk(dk);
    cr1->getCrystalProp();

    return cr1;
}


inline std::unique_ptr<EFields> config_efields(Crystal* Cr, const json& config)
{
    /* 
        This function configures and returns a unique pointer to an 
        EFields object based on the provided JSON configuration and
        the given Crystal object.
    */
    real_t Lcr = config["crystal"]["dimensions"]["Lcr"].get<real_t>();
    real_t lp = config["crystal"]["wavelengths"]["pump"].get<real_t>();
    real_t ls = config["crystal"]["wavelengths"]["signal"].get<real_t>();
    real_t li = config["crystal"]["wavelengths"]["idler"].get<real_t>();    
    real_t Power = config["fields"]["pump"]["pump_power_W"].get<real_t>();
    real_t waist = config["fields"]["pump"]["waist_um"].get<real_t>();
    real_t FWHM = config["fields"]["pump"]["fwhm"].get<real_t>();
    real_t focalpoint = config["fields"]["pump"]["focal_point_factor"].get<real_t>() * Lcr;

    auto A = std::make_unique<EFields>(lp, ls, li, Power, waist, Cr);

    real_t t_window = config["fields"]["time_freq_vect"]["time_window"].get<real_t>();
    A->set_time_freq_vectors(t_window);

    std::string mode = config["fields"]["pump"]["mode"].get<std::string>();
    A->set_pump_field(Power, FWHM, waist, focalpoint, mode);
    A->Ap = A->Api;

    bool signal_input = config["fields"]["signal"]["noise_generation"].get<bool>();
    bool idler_input = config["fields"]["idler"]["noise_generation"].get<bool>();
    bool degenerate = config["fields"]["degenerate"]["enabled"].get<bool>();

    if (signal_input && degenerate) {
        A->noise_generator(A->As);
        A->Ai = A->As;
    } else if (signal_input && idler_input) {
        A->noise_generator(A->As);
        A->noise_generator(A->Ai);
    } else {
        throw std::runtime_error("[ERROR]: Invalid signal/idler/degenerate configuration");
    }

    return A;
}



inline void run_TWM(Solver *solver, const json& config)
{
    /*
        This function runs the TWM solver based on the provided
        Solver object and JSON configuration. It supports both
        single-pass and multi-pass simulation modes.
    */

    // set pump mode: 
	// (1) "waveplane-cw"; (2) "waveplane-pulsed"; (3) "focused-cw"; (4) "focused-pulsed"
	std::string mode = config["fields"]["pump"]["mode"].get<std::string>();
	bool multi_pass = config["mul_pass_scheme"]["multipass"].get<bool>();
	uint32_t npasses = config["mul_pass_scheme"]["npasses"].get<int>();
	bool add_phase_air = config["mul_pass_scheme"]["add_phase_air"].get<bool>();	
	std::vector<real_t> tmp = config["mul_pass_scheme"]["distances_in_air"].get<std::vector<real_t>>();
	rVech_t mirror_distances(tmp.begin(), tmp.end());
	
	if (multi_pass and (npasses > 1) ){
		solver->set_npasses(npasses);
		solver->set_mirror_distances(mirror_distances);

		if (add_phase_air){solver->set_air_phase();}

		if(mode == "waveplane-cw" or mode == "focused-cw"){
			std::cout << "        ---> Running Solver: multipass mode" << std::endl;
			solver->run_multipass();
		}
		else{std::cout << "        ---> Change the pump mode to waveplane-cw or focused-cw" << std::endl;}
	}
	else if (npasses == 1 and !multi_pass){
		std::cout << "        ---> Running Solver: single-pass mode" << std::endl;
		solver->run_single_pass();
	}
	else
		{std::cout << "        ---> Error. Number of passes must be equal or grater than 1" << std::endl;}
    
    return ;

}

#endif
