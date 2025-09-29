/*---------------------------------------------------------------------------*/
// * This file set the policy for output files
// * the optical cavity as well as some intracavity elements.
/*---------------------------------------------------------------------------*/


#ifndef _SAVEOUTPUTS
#define _SAVEOUTPUTS


void init_full_field_h5(const std::string& filename, int NZtot) {
    using namespace H5;
    hsize_t dims[3] = { NX, NY, NZtot };

    H5File file(filename, H5F_ACC_TRUNC); // Crea o sobreescribe
    DataSpace dataspace(3, dims);

    file.createDataSet("real", PredType::NATIVE_FLOAT, dataspace);
    file.createDataSet("imag", PredType::NATIVE_FLOAT, dataspace);

    file.close();
}


void write_zslice_h5(const cVech_t& field, const std::string& filename, int z_slice) {
    using namespace H5;

    // Open file in read-write mode
    H5File file(filename, H5F_ACC_RDWR);
    DataSet dset_real = file.openDataSet("real");
    DataSet dset_imag = file.openDataSet("imag");

    hsize_t dims[3] = { NX, NY, 1 }; 

    // Define "hyperslab": [NX, NY, 1] in position z_slice
    hsize_t offset[3] = { 0, 0, static_cast<hsize_t>(z_slice) };
    DataSpace filespace = dset_real.getSpace();
    filespace.selectHyperslab(H5S_SELECT_SET, dims, offset);

    // Prepare slice real/imag plane [NX, NY]
    std::vector<real_t> arr_real(NX * NY), arr_imag(NX * NY);
    for (int y = 0; y < NY; ++y)
        for (int x = 0; x < NX; ++x) {
            int idx = x + NX * y;
            arr_real[idx] = field[idx].x;
            arr_imag[idx] = field[idx].y;
        }

    // Memory dataspace for 1 slice
    DataSpace memspace(3, dims);

    dset_real.write(arr_real.data(), PredType::NATIVE_FLOAT, memspace, filespace);
    dset_imag.write(arr_imag.data(), PredType::NATIVE_FLOAT, memspace, filespace);

    file.close();
}


void save_matrix_complex_h5_time(const cVech_t& field, const std::string& filename)
{
    using namespace H5;
    hsize_t dims[3] = { static_cast<hsize_t>(NX), static_cast<hsize_t>(NY), static_cast<hsize_t>(NT) }; // <-- ¡Orden deseado!

    H5File file(filename, H5F_ACC_TRUNC);
    DataSpace dataspace(3, dims);

    DataSet dset_real = file.createDataSet("real", PredType::NATIVE_FLOAT, dataspace);
    DataSet dset_imag = file.createDataSet("imag", PredType::NATIVE_FLOAT, dataspace);

    // Reordena el buffer para que coincida con el orden [NX, NY, NT]
    std::vector<real_t> arr_real(NX * NY * NT);
    std::vector<real_t> arr_imag(NX * NY * NT);

    for (int t = 0; t < NT; ++t) {
        for (int y = 0; y < NY; ++y) {
            for (int x = 0; x < NX; ++x) {
                int idx_flat = x + NX * (y + NY * t);    // Acceso en tu memoria original
                int idx_h5   = x + NX * (y + NY * t);    // Igual, porque vamos a guardar en ese orden
                arr_real[idx_h5] = field[idx_flat].x;
                arr_imag[idx_h5] = field[idx_flat].y;
            }
        }
    }

    dset_real.write(arr_real.data(), PredType::NATIVE_FLOAT);
    dset_imag.write(arr_imag.data(), PredType::NATIVE_FLOAT);

    file.close();
    return ;
}


void save_output_slices_XY(EFields *A, json config) {
    bool save_XY = config["save_mode"]["save_fields_XY_times_slides"].get<bool>();
    bool save_p  = config["save_mode"]["save_pump"].get<bool>();
    bool save_s  = config["save_mode"]["save_signal"].get<bool>();
    bool save_i  = config["save_mode"]["save_idler"].get<bool>();

    if (save_XY && save_p)
        save_matrix_complex_h5_time(A->Ap, "pump_output_XY.h5");
    if (save_XY && save_s)
        save_matrix_complex_h5_time(A->As, "signal_output_XY.h5");
    if (save_XY && save_i)
        save_matrix_complex_h5_time(A->Ai, "idler_output_XY.h5");
}


void save_input_pump_slices_XY (EFields *A, json config)
{
    // Define grid size
    bool save_p = config["save_mode"]["save_in_pump"].get<bool>();
    
    if (save_p)
        save_matrix_complex_h5_time(A->Api, "pump_input_XY.h5");
    
    return ;
}


void save_timedomain_fields (EFields *A, json config)
{
    // Save output signal electric field in the temporal domain.
	// This subroutine take each `x,y` pair and save the corresponding the time vector.

    // Define grid size
    bool save_time   = config["save_mode"]["save_timedomain_fields"].get<bool>();
    bool save_p      = config["save_mode"]["save_pump"].get<bool>();
    bool save_s      = config["save_mode"]["save_signal"].get<bool>();
    bool save_i      = config["save_mode"]["save_idler"].get<bool>();
    
    if (save_time and save_p)
    {
        for(uint x = 0; x < NX; x++){
            for(uint y = 0; y < NY; y++){
                save_electric_field3D_time (A->Ap, x, y, "pump_output_time_vector_x_" + std::to_string(x) + "_y_" + std::to_string(y));
            }
        }
    }

    if (save_time and save_s)
    {
        for(uint x = 0; x < NX; x++){
            for(uint y = 0; y < NY; y++){
                save_electric_field3D_time (A->As, x, y, "signal_output_time_vector_x_" + std::to_string(x) + "_y_" + std::to_string(y));
            }
        }
    }

    if (save_time and save_i)
    {
        for(uint x = 0; x < NX; x++){
            for(uint y = 0; y < NY; y++){
                save_electric_field3D_time (A->Ai, x, y, "idler_output_time_vector_x_" + std::to_string(x) + "_y_" + std::to_string(y));
            }
        }
    }   

    return ;
}


void save_timedomain_fields_center (EFields *A, json config)
{
    // Save output signal electric field in the temporal domain.
	// This subroutine take only `x=NX/2,y=NY/2` pair and save the corresponding the time vector.

    // Define grid size
    bool save_time   = config["save_mode"]["save_timedomain_fields_center"].get<bool>();
    bool save_p      = config["save_mode"]["save_pump"].get<bool>();
    bool save_s      = config["save_mode"]["save_signal"].get<bool>();
    bool save_i      = config["save_mode"]["save_idler"].get<bool>();

    if (save_time and save_p)
    {
        int x=NX/2,y=NY/2;
        save_electric_field3D_time (A->Ap, x, y, "pump_output_time_vector_x_" + std::to_string(x) + "_y_" + std::to_string(y));
    }

    if (save_time and save_s)
    {
        int x=NX/2,y=NY/2;
        save_electric_field3D_time (A->As, x, y, "signal_output_time_vector_x_" + std::to_string(x) + "_y_" + std::to_string(y));
    }

    if (save_time and save_i)
    {
        int x=NX/2,y=NY/2;
        save_electric_field3D_time (A->Ai, x, y, "idler_output_time_vector_x_" + std::to_string(x) + "_y_" + std::to_string(y));
    }   

    return ;
}


void save_full_field_inside_crystal(EFields *A, json config, uint32_t z_slice)
{

    bool save_inside = config["save_mode"]["save_full_field_inside_crystal"].get<bool>();
    bool save_p      = config["save_mode"]["save_pump"].get<bool>();
    bool save_s      = config["save_mode"]["save_signal"].get<bool>();
    bool save_i      = config["save_mode"]["save_idler"].get<bool>();

    if (save_inside and save_p)
    {
        save_matrix_complex_gpu ( A->Ap, "pump_output_XY_z_slice_" + std::to_string(z_slice) );
    }

    if (save_inside and save_s)
    {
        save_matrix_complex_gpu ( A->As, "signal_output_XY_z_slice_" + std::to_string(z_slice) );
    }

    if (save_inside and save_i)
    {
        save_matrix_complex_gpu ( A->Ai, "idler_output_XY_z_slice_" + std::to_string(z_slice) );
    }
    return ;
}


void save_time_and_frequency_vectors(EFields *A, json config)
{
    // Define grid size
    bool save = config["save_mode"]["save_time_and_freq_vectors"].get<bool>();
    
    if (save)
    {
        save_vector_real_gpu (A->t, "time");
        save_vector_real_gpu (A->F, "frequency");
    }

    return ;
}


#endif // -> #ifdef _SAVEOUTPUTS







// void save_output_slices_XY (EFields *A, json config)
// {
//     // Define grid size
//     bool save_XY     = config["save_mode"]["save_fields_XY_times_slides"].get<bool>();
//     bool save_p      = config["save_mode"]["save_pump"].get<bool>();
//     bool save_s      = config["save_mode"]["save_signal"].get<bool>();
//     bool save_i      = config["save_mode"]["save_idler"].get<bool>();
    
//     if (save_XY and save_p)
//     {
//         for (uint slice = 0; slice < NT; slice ++)
//         {
//             if (NT > 1)
//                 save_matrix_complex_time_slice_gpu ( A->Ap, slice, "pump_output_XY_timeslice_" + std::to_string(slice) );
//             else
//                 save_matrix_complex_time_slice_gpu ( A->Ap, slice, "pump_output_XY");
//         }
//     }

//     if (save_XY and save_s)
//     {
//         for (uint slice = 0; slice < NT; slice ++)
//         {
//             if (NT > 1)
//                 save_matrix_complex_time_slice_gpu ( A->As, slice, "signal_output_XY_timeslice_" + std::to_string(slice) );
//             else
//                 save_matrix_complex_time_slice_gpu ( A->As, slice, "signal_output_XY");    
//         }
//     }

//     if (save_XY and save_i)
//     {
//         for (uint slice = 0; slice < NT; slice ++)
//         {
//             if (NT > 1)
//                 save_matrix_complex_time_slice_gpu ( A->Ai, slice, "idler_output_XY_timeslice_" + std::to_string(slice) );
//             else
//                 save_matrix_complex_time_slice_gpu ( A->Ai, slice, "idler_output_XY");
//         }
//     }   

//     return ;
// }


// void save_input_pump_slices_XY (EFields *A, json config)
// {
//     // Define grid size
//     bool save_p = config["save_mode"]["save_in_pump"].get<bool>();
    
//     if (save_p)
//     {
//         for (uint slice = 0; slice < NT; slice ++)
//         {
//             if (NT > 1)
//                 save_matrix_complex_time_slice_gpu ( A->Api, slice, "pump_input_XY_timeslice_" + std::to_string(slice) );
//             else
//                 save_matrix_complex_time_slice_gpu ( A->Api, slice, "pump_input_XY");
//         }
//     }
//     return ;
// }