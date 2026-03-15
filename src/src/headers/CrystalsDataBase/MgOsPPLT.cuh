#ifndef _MGOSPPLT
#define _MGOSPPLT



// This class contains the Sellmeier equations and properties for the MgO:sPPLT nonlinear crystal and other χ^(2) properties.
// Sellmeier equations from: O. Gayer et al., Appl. Phys. B 91, 343–348 (2008).

class MgOsPPLT : public Crystal
{
public:
    bool QPM;

    MgOsPPLT(real_t _LX, real_t _LY, real_t _Lcr, real_t _T, real_t _Lambda, real_t _lp, real_t _ls, real_t _li)
    {
        LX = _LX; LY = _LY; Lcr = _Lcr;
        T = _T; Lambda = _Lambda;
        lp = _lp; ls = _ls; li = _li;

        print_line_on_screen();
        printf("\nInitialize crystal MgOsPPLT.\n");

        name = "MgO:sPPLT";
        QPM = true;
        d33 = 13.7e-6;      dQ = 2.0 * d33 / PI;
        alpha_crp = 0.021e-4; alpha_crs = 0.002e-4; alpha_cri = 0.002e-4;

        dx = static_cast<real_t>(LX / (NX - 1.0f));
        dy = static_cast<real_t>(LY / (NY - 1.0f));
        dz = static_cast<real_t>(Lcr / (NZ - 1.0f));

        np = this->n(lp, T); ns = this->n(ls, T); ni = this->n(li, T);
        vp = this->GV(lp, T); vs = this->GV(ls, T); vi = this->GV(li, T);
        b2p = this->GVD(lp, T); b2s = this->GVD(ls, T); b2i = this->GVD(li, T);
        b3p = this->TOD(lp, T); b3s = this->TOD(ls, T); b3i = this->TOD(li, T);

        kp = 2.0f * PI * ns / lp;
        ks = 2.0f * PI * ns / ls;
        ki = 2.0f * PI * ni / ls;
    }

    ~MgOsPPLT() override { printf("MgOsPPLT deleted\n"); }
    
    void set_dk(real_t dk) override;
    real_t n(real_t L, real_t T);
    real_t resonances(real_t L, real_t T, int p);
    real_t dndl(real_t L, real_t T);
    real_t d2ndl2(real_t L, real_t T);
    real_t d3ndl3(real_t L, real_t T);
    real_t GV(real_t L, real_t T);
    real_t GVD(real_t L, real_t T);
    real_t TOD(real_t L, real_t T);
    void getCrystalProp() override;
    std::string get_crystal_type() const override { return "MgOsPPLT"; }
};

// ---------- Methods  ----------

void MgOsPPLT::set_dk(real_t dk)
{
    this->dk = dk;
    return ;
}


/** Returns the MgO:sPPLT extraordinary refractive index */
real_t MgOsPPLT::n(real_t L, real_t T)
{
    // The following values are taken from O. Gayer et al.
    real_t f = (T - 24.5) * (T + 570.82);
    real_t a1 = 5.113;
    real_t a2 = 0.0996;
    real_t a3 = 0.2102;
    real_t a4 = 189.69;
    real_t a5 = 12.48;
    real_t a6 = 1.32e-2;
    real_t b1 = 2.767e-6;
    real_t b2 = 3.728e-8;
    real_t b3 = 5.290e-8;
    real_t b4 = 1.275e-4;
    real_t G1 = a1 + b1 * f;
    real_t G2 = a2 + b2 * f;
    real_t G3 = a3 + b3 * f;
    real_t G4 = a4 + b4 * f;
    return sqrtf(G1 + G2 / (powf(L, 2) - powf(G3, 2)) + G4 / (powf(L, 2) - powf(a5, 2)) - a6 * L * L);
}

/** Auxiliary function related with the resonances */
real_t MgOsPPLT::resonances(real_t L, real_t T, int p)
{
    real_t f = (T - 24.5) * (T + 570.82);
    real_t a1 = 5.113;
    real_t a2 = 0.0996;
    real_t a3 = 0.2102;
    real_t a4 = 189.69;
    real_t a5 = 12.48;
    real_t a6 = 1.32e-2;
    real_t b1 = 2.767e-6;
    real_t b2 = 3.728e-8;
    real_t b3 = 5.290e-8;
    real_t b4 = 1.275e-4;
    real_t G1 = a1 + b1 * f;
    real_t G2 = a2 + b2 * f;
    real_t G3 = a3 + b3 * f;
    real_t G4 = a4 + b4 * f;
    return G2 / powf((powf(L, 2) - powf(G3, 2)), p) + G4 / powf((powf(L, 2) - powf(a5, 2)), p);
}

/** Returns dn/dλ */
real_t MgOsPPLT::dndl(real_t L, real_t T)
{
    real_t f = (T - 24.5) * (T + 570.82);
    real_t a1 = 5.113;
    real_t a2 = 0.0996;
    real_t a3 = 0.2102;
    real_t a4 = 189.69;
    real_t a5 = 12.48;
    real_t a6 = 1.32e-2;
    real_t b1 = 2.767e-6;
    real_t b2 = 3.728e-8;
    real_t b3 = 5.290e-8;
    real_t b4 = 1.275e-4;
    real_t G1 = a1 + b1 * f;
    real_t G2 = a2 + b2 * f;
    real_t G3 = a3 + b3 * f;
    real_t G4 = a4 + b4 * f;

    return -L * (resonances(L, T, 2) + a6) / n(L, T);
}

/** Returns d²n/dλ² */
real_t MgOsPPLT::d2ndl2(real_t L, real_t T)
{
    real_t f = (T - 24.5) * (T + 570.82);
    real_t a1 = 5.113;
    real_t a2 = 0.0996;
    real_t a3 = 0.2102;
    real_t a4 = 189.69;
    real_t a5 = 12.48;
    real_t a6 = 1.32e-2;
    real_t b1 = 2.767e-6;
    real_t b2 = 3.728e-8;
    real_t b3 = 5.290e-8;
    real_t b4 = 1.275e-4;
    real_t G1 = a1 + b1 * f;
    real_t G2 = a2 + b2 * f;
    real_t G3 = a3 + b3 * f;
    real_t G4 = a4 + b4 * f;

    real_t A = (L * dndl(L, T) / powf(n(L, T), 2) - 1 / n(L, T)) * (resonances(L, T, 2) + a6);
    real_t B = 4 * L * L / n(L, T) * resonances(L, T, 3);

    return A + B;
}

/** Returns d³n/dλ³ */
real_t MgOsPPLT::d3ndl3(real_t L, real_t T)
{
    real_t f = (T - 24.5) * (T + 570.82);
    real_t a1 = 5.113;
    real_t a2 = 0.0996;
    real_t a3 = 0.2102;
    real_t a4 = 189.69;
    real_t a5 = 12.48;
    real_t a6 = 1.32e-2;
    real_t b1 = 2.767e-6;
    real_t b2 = 3.728e-8;
    real_t b3 = 5.290e-8;
    real_t b4 = 1.275e-4;
    real_t G1 = a1 + b1 * f;
    real_t G2 = a2 + b2 * f;
    real_t G3 = a3 + b3 * f;
    real_t G4 = a4 + b4 * f;

    real_t A1 = (2 * dndl(L, T) + L * d2ndl2(L, T)) / powf(n(L, T), 2);
    real_t A2 = -2 * L * powf(dndl(L, T), 2) / powf(n(L, T), 3);
    real_t AA = (A1 + A2) * (resonances(L, T, 2) + a6);
    real_t B1 = 12 * L / n(L, T);
    real_t B2 = -4 * L * L / n(L, T) * d2ndl2(L, T) * (1 - 1 / n(L, T));
    real_t BB = (B1 + B2) * resonances(L, T, 3);
    real_t CC = -24 * L * L * L / n(L, T) * resonances(L, T, 4);

    return AA + BB + CC;
}

/** Returns the group velocity vg(λ) = c / (n(λ) - λ dn/dλ) */
real_t MgOsPPLT::GV(real_t L, real_t T)
{
    return C / (n(L, T) - L * dndl(L, T));
}

/** Returns the group velocity dispersion β2(λ) = λ^3/(2πc²) * d²n/dλ² */
real_t MgOsPPLT::GVD(real_t L, real_t T)
{
    return powf(L, 3) * d2ndl2(L, T) / (2 * PI * C * C);
}

/** Returns the third-order dispersion β3(λ) = -λ^4/(4π²c³) * [3 d²n/dλ² + λ d³n/dλ³] */
real_t MgOsPPLT::TOD(real_t L, real_t T)
{
    return -powf(L, 4) / (4 * PI * PI * C * C * C) * (3 * d2ndl2(L, T) + L * d3ndl3(L, T));
}

void MgOsPPLT::getCrystalProp()
{
    std::cout << "Crystal name = " << this->name << std::endl;
    std::cout << "        ---> Temp            = " << T << " ºC" << std::endl;
    std::cout << "        ---> lp              = " << this->lp << " um" << std::endl;
    std::cout << "        ---> ls              = " << this->ls << " um" << std::endl;
    std::cout << "        ---> li              = " << this->li << " um" << std::endl;
    std::cout << "        ---> np              = " << this->np << std::endl;
    std::cout << "        ---> ns              = " << this->ns << std::endl;
    std::cout << "        ---> ni              = " << this->ni << std::endl;
    std::cout << "        ---> vgp             = " << this->vp << " um/ps" << std::endl;
    std::cout << "        ---> vgs             = " << this->vs << " um/ps" << std::endl;
    std::cout << "        ---> vgi             = " << this->vi << " um/ps" << std::endl;
    std::cout << "        ---> b2p             = " << this->b2p << " ps²/um" << std::endl;
    std::cout << "        ---> b2s             = " << this->b2s << " ps²/um" << std::endl;
    std::cout << "        ---> b2i             = " << this->b2i << " ps²/um" << std::endl;
    std::cout << "        ---> dQ              = " << this->dQ * 1e6 << " pm/V" << std::endl;
    std::cout << "        ---> Λ               = " << this->Lambda << " μm" << std::endl;
    std::cout << "        ---> αcp             = " << this->alpha_crp << " μm⁻¹" << std::endl;
    std::cout << "        ---> αcs             = " << this->alpha_crs << " μm⁻¹" << std::endl;
    std::cout << "        ---> αci             = " << this->alpha_cri << " μm⁻¹" << std::endl;
    std::cout << "        ---> dx              = " << this->dx << " μm" << std::endl;
    std::cout << "        ---> dy              = " << this->dy << " μm" << std::endl;
    std::cout << "        ---> dz              = " << this->dz << " μm" << std::endl;
    std::cout << "        ---> \u0394k              = " << this->dk << " \u03BCm⁻¹"  << std::endl;        
    std::cout << "        ---> Crystal length  = " << this->Lcr * 1e-3 << " mm\n" << std::endl;
    return;
}

#endif // _MGOSPPLT
