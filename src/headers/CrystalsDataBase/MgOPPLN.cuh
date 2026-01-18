#ifndef _MGOPPLN 
#define _MGOPPLN 


// MgOPPLN nonlinear crystal: Sellmeier equations and χ² material properties.
// Reference: O. Gayer et al.

class MgOPPLN : public Crystal
{
public:
    bool QPM;

    MgOPPLN(real_t _LX, real_t _LY, real_t _Lcr, real_t _T, real_t _Lambda, real_t _lp, real_t _ls, real_t _li)
    {
        LX = _LX; LY = _LY; Lcr = _Lcr;
        T = _T; Lambda = _Lambda;
        lp = _lp; ls = _ls; li = _li;

        print_line_on_screen();
        printf("\nInitialize crystal MgOPPLN.\n");

        name = "MgO:PPLN";
        QPM = true;
        d33 = 25.20e-6;    dQ = 2.0 * d33 / PI;
        alpha_crp = 0.002e-4; alpha_crs = 0.025e-4; alpha_cri = 0.025e-4;
        
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

    ~MgOPPLN() override { printf("MgOPPLN deleted\n"); }
    
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
    std::string get_crystal_type() const override { return "MgOPPLN"; }
};

// ---------- Methods  ----------

void MgOPPLN::set_dk(real_t dk)
{
    this->dk = dk;
    return ;
}


/** This function returns the MgO:PPLN extraordinary refractive index */
real_t MgOPPLN::n(real_t L, real_t T)
{
    real_t f = (T - 24.5) * (T + 570.82);
    real_t a1 = 5.756;
    real_t a2 = 0.0983;
    real_t a3 = 0.2020;
    real_t a4 = 189.32;
    real_t a5 = 12.52;
    real_t a6 = 1.32e-2;
    real_t b1 = 2.860e-6;
    real_t b2 = 4.700e-8;
    real_t b3 = 6.113e-8;
    real_t b4 = 1.516e-4;
    real_t G1 = a1 + b1 * f;
    real_t G2 = a2 + b2 * f;
    real_t G3 = a3 + b3 * f;
    real_t G4 = a4 + b4 * f;
    return sqrtf(G1 + G2 / (powf(L, 2) - powf(G3, 2)) + G4 / (powf(L, 2) - powf(a5, 2)) - a6 * L * L);
}

/** This function is an auxiliary function related with the resonances */
real_t MgOPPLN::resonances(real_t L, real_t T, int p)
{
    real_t f = (T - 24.5) * (T + 570.82);
    real_t a1 = 5.756;
    real_t a2 = 0.0983;
    real_t a3 = 0.2020;
    real_t a4 = 189.32;
    real_t a5 = 12.52;
    real_t a6 = 1.32e-2;
    real_t b1 = 2.860e-6;
    real_t b2 = 4.700e-8;
    real_t b3 = 6.113e-8;
    real_t b4 = 1.516e-4;
    real_t G1 = a1 + b1 * f;
    real_t G2 = a2 + b2 * f;
    real_t G3 = a3 + b3 * f;
    real_t G4 = a4 + b4 * f;
    return G2 / powf((powf(L, 2) - powf(G3, 2)), p) + G4 / powf((powf(L, 2) - powf(a5, 2)), p);
}

/** Returns the first-order derivative of the refractive index respect to the wavelength dn/dλ. */
real_t MgOPPLN::dndl(real_t L, real_t T)
{
    real_t f = (T - 24.5) * (T + 570.82);
    real_t a1 = 5.756;
    real_t a2 = 0.0983;
    real_t a3 = 0.2020;
    real_t a4 = 189.32;
    real_t a5 = 12.52;
    real_t a6 = 1.32e-2;
    real_t b1 = 2.860e-6;
    real_t b2 = 4.700e-8;
    real_t b3 = 6.113e-8;
    real_t b4 = 1.516e-4;
    real_t G1 = a1 + b1 * f;
    real_t G2 = a2 + b2 * f;
    real_t G3 = a3 + b3 * f;
    real_t G4 = a4 + b4 * f;
    return -L * (resonances(L, T, 2) + a6) / n(L, T);
}

/** Returns the second-order derivative of the refractive index respect to the wavelength d²n/dλ². */
real_t MgOPPLN::d2ndl2(real_t L, real_t T)
{
    real_t f = (T - 24.5) * (T + 570.82);
    real_t a1 = 5.756;
    real_t a2 = 0.0983;
    real_t a3 = 0.2020;
    real_t a4 = 189.32;
    real_t a5 = 12.52;
    real_t a6 = 1.32e-2;
    real_t b1 = 2.860e-6;
    real_t b2 = 4.700e-8;
    real_t b3 = 6.113e-8;
    real_t b4 = 1.516e-4;
    real_t G1 = a1 + b1 * f;
    real_t G2 = a2 + b2 * f;
    real_t G3 = a3 + b3 * f;
    real_t G4 = a4 + b4 * f;
    real_t A = (L * dndl(L, T) / powf(n(L, T), 2) - 1 / n(L, T)) * (resonances(L, T, 2) + a6);
    real_t B = 4 * L * L / n(L, T) * resonances(L, T, 3);
    return A + B;
}

/** Returns the third-order derivative of the refractive index respect to the wavelength d³n/dλ³. */
real_t MgOPPLN::d3ndl3(real_t L, real_t T)
{
    real_t f = (T - 24.5) * (T + 570.82);
    real_t a1 = 5.756;
    real_t a2 = 0.0983;
    real_t a3 = 0.2020;
    real_t a4 = 189.32;
    real_t a5 = 12.52;
    real_t a6 = 1.32e-2;
    real_t b1 = 2.860e-6;
    real_t b2 = 4.700e-8;
    real_t b3 = 6.113e-8;
    real_t b4 = 1.516e-4;
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

/** Returns the group-velocity vg(λ) = c/(n(λ)-λdn/dλ). */
real_t MgOPPLN::GV(real_t L, real_t T)
{
    return C / (n(L, T) - L * dndl(L, T));
}

/** Returns the group-velocity β2(λ)=λ^3/(2πc²)(d²n/dλ²). */
real_t MgOPPLN::GVD(real_t L, real_t T)
{
    return powf(L, 3) * d2ndl2(L, T) / (2 * PI * C * C);
}

/** Returns the TOD β3(λ)=-λ^4/(4π²c³)[3.d²n/dλ² + λ.d³n/dλ³]. */
real_t MgOPPLN::TOD(real_t L, real_t T)
{
    return -powf(L, 4) / (4 * PI * PI * C * C * C) * (3 * d2ndl2(L, T) + L * d3ndl3(L, T));
}

void MgOPPLN::getCrystalProp()
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

#endif // _MGOPPLN
