// * Author name: Alfredo Daniel Sanchez
// * email:       alfredo.daniel.sanchez@gmail.com


#include "headers/Libraries.cuh"		// Required libraries
#include "headers/PackageLibraries.cuh"	// Required package libraries



int main(int argc, char *argv[]){

	///////////////////////////////////////////////////////////////////////////////////
	// 1. Read config.json data set that contains all the required parameters 
	
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

	// Print grid information: tensors size, blocks size, etc.
	print_grid();

	///////////////////////////////////////////////////////////////////////////////////
	// 2. Build the system

	// a. set the nonlinear crystal
	auto cr1 = config_crystal(config);
	if (!cr1) {return 1;}  // Exit if crystal configuration failed
	
	// b. set electric fields	
	auto A = config_efields(cr1.get(), config);

	// c. set solver
	auto solver = std::make_unique<Solver>(cr1.get(), A.get(), config);

	///////////////////////////////////////////////////////////////////////////////////
	// 3. run package
	run_TWM(solver.get(), config);
	
	
	///////////////////////////////////////////////////////////////////////////////////
	// 4. Save output data
	
	save_all_outputs (A.get(), config);
	
	return 0;
	
}
