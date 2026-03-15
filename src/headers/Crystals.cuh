
#ifndef _CRYSTALS
#define _CRYSTALS


// Abstract base class for all nonlinear crystals
class Crystal {
    public:
        // Crystal geometry and material properties
        real_t Lcr;
        real_t LX;
        real_t LY;
        real_t T;
        real_t Lambda;
        real_t lp, ls, li;
        // Nonlinear and absorption constants
        real_t d33;
        real_t dQ;
        real_t alpha_crp;
        real_t alpha_crs;
        real_t alpha_cri;
        real_t beta_crs;
        real_t chi3p;
        real_t chi3s;
        real_t chi3i;
        real_t np, ns, ni;
        real_t vp, vs, vi;
        real_t b2p, b2s, b2i;
        real_t b3p, b3s, b3i;
        real_t dx, dy, dz;
        real_t kp, ks, ki, dk;
        std::string name;

        // Crystal() {}
        virtual ~Crystal() {}
    
        // Print the properties of the crystal (must be implemented by children)
        virtual void getCrystalProp() = 0;
    
        // Return the crystal type as a string (must be implemented by children)
        virtual std::string get_crystal_type() const = 0;

        virtual void set_dk(real_t dk) = 0;
    };

#include "CrystalsDataBase/PPLN.cuh"
#include "CrystalsDataBase/MgOPPLN.cuh"
#include "CrystalsDataBase/MgOsPPLT.cuh"
#include "CrystalsDataBase/ZGP.cuh"

#endif // _CRYSTALS