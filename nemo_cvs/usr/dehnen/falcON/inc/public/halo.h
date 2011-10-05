// -*- C++ -*-                                                                  
////////////////////////////////////////////////////////////////////////////////
///                                                                             
/// \file   inc/public/halo.h                                                   
///                                                                             
/// \brief  classes for generating spherical halo models                        
///                                                                             
/// \author Paul McMillan                                                       
/// \author Walter Dehnen                                                       
/// \date   2000-2008                                                           
///                                                                             
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// Copyright (C) 2000-2008  Walter Dehnen, Paul McMillan                        
//                                                                              
// This program is free software; you can redistribute it and/or modify         
// it under the terms of the GNU General Public License as published by         
// the Free Software Foundation; either version 2 of the License, or (at        
// your option) any later version.                                              
//                                                                              
// This program is distributed in the hope that it will be useful, but          
// WITHOUT ANY WARRANTY; without even the implied warranty of                   
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU            
// General Public License for more details.                                     
//                                                                              
// You should have received a copy of the GNU General Public License            
// along with this program; if not, write to the Free Software                  
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                    
//                                                                              
////////////////////////////////////////////////////////////////////////////////
#ifndef falcON_included_halo_h
#define falcON_included_halo_h

#ifndef falcON_included_externacc_h
#  include <externacc.h>
#endif
#ifndef falcON_included_sample_h
#  include <public/sample.h>
#endif
////////////////////////////////////////////////////////////////////////////////
namespace WDutils { template<class C, class B> class Pspline; }
namespace falcON {
  using namespace WDutils;
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::HaloDensity                                                  
  //                                                                            
  /// abstract base class for a spherical halo model                            
  /// used in the construction of halo equilibrium and initial conditions       
  ///                                                                           
  // ///////////////////////////////////////////////////////////////////////////
  class HaloDensity {
  public:
    /// pure virtual function: density at given radius
    /// \return density
    /// \param r (input) radius
    virtual double operator()(double r) const = 0;
    /// pure virtual function: negative logarithmic density slope at r->0
    /// \return gamma(r->0)
    virtual double inner_gamma() const = 0;
    /// pure virtual function: scale radius
    /// \return scale radius
    virtual double scale_radius() const = 0;
    /// pure virtual function: truncation radius (if any)
    /// \return truncation radius, zero if density decays like power law
    virtual double trunc_radius() const = 0;
    /// pure virtual function: negative logarithmic density slope at r->oo
    /// \return gamma(r->oo)
    virtual double outer_gamma() const = 0;
    /// pure virtual function: density and its first derivative
    /// \return density
    /// \param r (input) radius
    /// \param rh1 (output) drho/dr
    virtual double operator()(double r, double&rh1) const = 0;
    /// pure virtual function: density and its first two derivatives
    /// \return density
    /// \param r (input) radius
    /// \param rh1 (output) drho/dr
    /// \param rh2 (output) d^2rho/dr^2
    virtual double operator()(double r, double&rh1, double&rh2) const = 0;
    /// noop dtor
    virtual~HaloDensity() {}
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::HaloPotential                                                
  //                                                                            
  /// given a HaloDensity, we find the potential generated by it and an optional
  /// external monopole potential.                                              
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class HaloPotential {
  protected:
    const HaloDensity &DEN;                        ///< density of halo         
    const acceleration*MON;                        ///< external monopole       
    double             Ah,At;                      ///< gamma_0 for rho, psi    
    double             Ch;                         ///< gamma_00 for rho        
    double             ps0;                        ///< Psi(0) (if finite)      
    int                n,n0,n1,nm;                 ///< size of tables         
    bool               tr;                         ///< truncated?
    double             Mht,Mtt,go,g3,g4,fmt,fmh;   ///< for r>r[n1]
    Array<double,1>    r,lr;                       ///< r_i, ln(r_i)            
    Array<double,1>    mh;                         ///< M_halo(<r_i)            
    Array<double,1>    mt;                         ///< M_tot(<r_i)             
    Array<double,1>    ps;                         ///< Psi_tot(r_i)            
    Array<double,1>    rh;                         ///< rho_tot(r_i)            
    Array<double,1>    ec;                         ///< E_c(r_i)                
    Array<double,1>    dp;                         ///< -v_c^2(r_i) = dPsi/dlnr 
    Pspline<double,double>   *PS;                  ///< penta spline: Psi(r)    
  public:
    /// constructor
    /// \param halo  density model for halo
    /// \param mono  external monopole potential
    /// \param r_max maximum radius for tables
    /// \note When using an external potential the density of which extends
    ///       well beyond that of \a halo, it is advisable to give \a r_max
    ///       to be a radius beyond which the \b total potential is well
    ///       approximated by GM/r.
    HaloPotential(HaloDensity const&halo, const acceleration*mono,
		  double r_max=0.) falcON_THROWING;
    /// destructor
    ~HaloPotential();
    /// potential Psi(r) using polynomial interpolation
    double Ps(double R) const;
    /// potential Phi(r), and (-dPhi/dr)/r using penta spline
    /// \return Phi(r)
    /// \param  rq  (input) radius squared
    /// \param  acx (output) (-dPhi/dr)/r, so that acc = pos * acx
    double PotAcc(double rq, double&acx) const;
    /// total cumulative mass
    double Mt(double R) const;
    /// halo cumulative mass
    double Mh(double R) const;
    /// total halo mass
    double Mh_tot() const { return Mht; }
    /// total mass density
    double rhot(double R) const;
    /// v_circ^2(r)
    double vcq(double R) const;
    /// Omega^2(r)
    double omq(double R) const;
    /// kappa^2(r)
    double kpq(double R) const;
    /// gamma := 2*Omega/kappa 
    double gam(double R) const;
    // ln R_psi(E)
    double lnRPsi(double P) const;
    // R_psi(E)
    double RPsi(double P) const;
    /// Eps_c(r)
    double Epc(double R) const;
    /// R_circ(Eps)
    double RcE(double E) const;
    /// estimate for R(Eps, L^2, cos[eta])
    double Rap(double E, double Lq, double ce) const;
    /// estimate for R_peri(Eps, L^2)
    double Rp(double E, double Lq) const;
    /// estimate for R_apo(Eps, L^2)
    double Ra(double E, double Lq) const;
    /// radius given cumulative halo mass
    double RMh(double M) const falcON_THROWING;
    //--------------------------------------------------------------------------
    /// \name constant access to some data
    //{@
    /// size of tables
    int Ntab() const { return n; }
    /// r_i
    double const&rad(int i) const { return r[i]; }
    /// ln(r_i)
    double const&lnr(int i) const { return lr[i]; }
    /// M_halo(r<r_i)
    double const&mhr(int i) const { return mh[i]; }
    /// Psi_tot(r_i)
    double const&psi(int i) const { return ps[i]; }
    /// M_tot(r<r_i)
    double const&mtr(int i) const { return mt[i]; }
    /// Eps_circ(r_i)
    double const&epc(int i) const { return ec[i]; }
    //@}
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::HaloModel                                                    
  //                                                                            
  /// given a HaloDensity, we construct a full spherical halo model with        
  /// Cuddeford (1991) distribution function (default: isotropic) in the        
  /// presence of an external spherical potential.                              
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class HaloModel : public HaloPotential {
  private:
    const double       B;                          ///< beta parameter of DF    
    const double       RA;                         ///< anisotropy radius of DF 
    Array<double,1>    lg;                         ///< table: log(g(Q))        
    bool               g_nonmon;
  public:
    /// constructor
    /// \param halo density model for halo
    /// \param beta anisotropy parameter beta 
    /// \param r_a anisotropy radius (0 maps to infinity)
    /// \param mono external monopole potential
    /// \param r_max maximum radius for tables
    /// \note When using an external potential the density of which extends
    ///       well beyond that of \a halo, it is advisable to give \a r_max
    ///       to be a radius beyond which the \b total potential is well
    ///       approximated by GM/r.
    HaloModel(HaloDensity const&halo, double beta, double r_a,
	      const acceleration*mono, double r_max=0.) falcON_THROWING;
    /// ln g(Q)
    double lnG(double Q) const falcON_THROWING;
    /// distribution function f(Eps,L^2)
    double fEL(double E, double Lq) const;
    /// ln g(Q=Psi_tot(r_i))
    double const&lng(int i) const { return lg[i]; }
    /// is g(Q) non-monotinic?
    bool const& non_monotonic() const { return g_nonmon; }
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  /// class HaloSampler                                                         
  //                                                                            
  /// For generating N-body realisations of a halo model                        
  // ///////////////////////////////////////////////////////////////////////////
  class HaloSampler :
    private HaloModel,
    public  SphericalSampler {
  public:
#ifdef falcON_PROPER
    /// constructor
    /// \param halo model for halo density
    /// \param beta beta of distritubtion function
    /// \param r_a  Ossipkov-Merritt anisotropy radius
    /// \param mono external monopole potential
    /// \param MArs mass adaption: scale radiu
    /// \param MAmm mass adaption: m_max/m_min
    /// \param MAet mass adaption: shape parameter
    /// \param MAnm mass adaption: n_max per (E,L) (default: mm)
    /// \param MApr mass adaption: using R_peri (or R_circ) for mass adaption ?
    /// \param care care for non-monotonic DF?
    /// \param r_max maximum radius for tables
    /// \note When using an external potential the density of which extends
    ///       well beyond that of \a halo, it is advisable to give \a r_max
    ///       to be a radius beyond which the \b total potential is well
    ///       approximated by GM/r.
    HaloSampler(HaloDensity const&halo, double bet, double r_a,
		const acceleration*mono,
		double MArs, double MAmm, double MAet, double MAnm, bool MApr,
		double r_max=0) :
      HaloModel(halo,bet,r_a,mono,r_max),
      SphericalSampler(Mh_tot(),r_a,bet,MArs,MAmm,MAet,MAnm,MApr,
		       non_monotonic())
    {}
#endif
    /// constructor
    /// \param halo model for halo density
    /// \param beta beta of distritubtion function
    /// \param r_a  Ossipkov-Merritt anisotropy radius
    /// \param mono external monopole potential
    /// \param care care for non-monotonic DF?
    /// \param r_max maximum radius for tables
    /// \note When using an external potential the density of which extends
    ///       well beyond that of \a halo, it is advisable to give \a r_max
    ///       to be a radius beyond which the \b total potential is well
    ///       approximated by GM/r.
    HaloSampler(HaloDensity const&halo, double bet, double r_a,
		const acceleration*mono, double r_max=0) :
      HaloModel(halo,bet,r_a,mono,r_max),
      SphericalSampler(Mh_tot(),r_a,bet,non_monotonic())
    {}
    //--------------------------------------------------------------------------
    HaloModel const&Model() const { return *this; }
    //--------------------------------------------------------------------------
    double Mt() const { return HaloModel::Mh_tot(); }
    double DF(double q)           const { return exp(HaloModel::lnG(q)); }
    double Ps(double _r)          const { return HaloModel::Ps(_r); }
    double rM(double m)           const { return HaloModel::RMh(m); }
    double Re(double e)           const { return HaloModel::RcE(e); } 
    double Rp(double e, double l) const { return HaloModel::Rp(e,l*l); }
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::HaloModifier                                                 
  //                                                                            
  /// given rho(r), get Rho(r) := rho(sqrt[r^2+c^2]) trunc(r/|rt|)              
  /// where trunc(x) = sech(x) for rt>0 and 2/(sech(x)+1/sech(x)) for rt<0      
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class HaloModifier {
    const double rc,rcq;              ///< core radius, core radius squared
    const double rt,irt;              ///< truncation radius & its inverse
    const bool   sechtr;              ///< sech(x) truncation?
  public:
    /// sqrt(r^2+rc^2)
    double core(double r) const;
    /// sqrt(r^2+rc^2) and its first derivative
    double core(double r, double&d1) const;
    /// sqrt(r^2+rc^2) and its first two derivatives
    double core(double r, double&d1, double&d2) const;
    /// cored but untrancated density at given radius
    double cored(HaloDensity const&m, double r) const;
    /// cored but untrancated density and its first derivative
    double cored(HaloDensity const&m, double r, double&rh1) const;
    /// cored but untrancated density and its first two derivatives
    double cored(HaloDensity const&m, double r, double&rh1, double&rh2) const;
    /// truncation factor
    double trunc(double r) const;
    /// truncation factor and its first derivative
    double trunc(double r, double&t1) const;
    /// truncation factor and its first two derivatives
    double trunc(double r, double&t1, double&t2) const;
    /// core radius
    double const&r_c() const { return rc; }
    /// truncation radius
    double const&r_t() const { return rt; }
    /// constructor
    /// \param c core radius
    /// \param t truncation radius
    HaloModifier(double c, double t) falcON_THROWING;
    /// modified density at given radius
    double operator()(HaloDensity const&m, double r) const;
    /// modified density and its first derivative
    double operator()(HaloDensity const&m, double r, double&rh1) const;
    /// modified density and its first two derivatives
    double operator()(HaloDensity const&m, double r, double&rh1, double&rh2)
      const;
    /// is the core radius actually non-zero?
    bool cored() const { return rcq; }
    /// is the truncation radius finite?
    bool truncated() const { return irt; }
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::DoublePowerLawHalo                                           
  //                                                                            
  /// implements a halo density model with density\n                            
  /// rho(r) = r^(-gamma_i) * (1 + r^eta)^((gamma_o-gamma_i)/eta)               
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class DoublePowerLawHalo : public HaloDensity {
    const double go,gi,et;
    const double gg,al;
  public:
    /// specifying the model
    enum Model {
      Default   = 0,  ///< Zhao (1996)               : gi=???  go=???   et=???
      Plummer   = 1,  ///< Plummer (1911)            : gi=0    go=5     et=2
      Jaffe     = 2,  ///< Jaffe (1983)              : gi=2    go=4     et=1
      Hernquist = 3,  ///< Hernquist (1990)          : gi=1    go=4     et=1
      Dehnen    = 4,  ///< Dehnen (1993)             : gi=???  go=4     et=1
      NFW       = 5,  ///< NFW (1996)                : gi=1    go=3     et=1
      Moore     = 6,  ///< Moore (1999)              : gi=3/2  go=3     et=3/2
      DM        = 7   ///< Dehnen & McLaughlin (2005): gi=7/9  go=31/9  et=4/9
    };
    /// find model from string
    /// \param[in] model "Plummer", "Jaffe", "Hernquist", "Dehnen", "Zhao", 
    ///                  "NFW", "Moore", or "DM"
    /// \param[in] warn  if true, a non-matching @a model triggers a warning
    static Model model(const char*model, bool warn=true);
    /// get name string from Model
    /// \param[in] model model to get name for
    static const char*name(Model model);
    static double null_value() { return -1.; }
    static double inner_value(Model, double val=-1.) falcON_THROWING;
    static double outer_value(Model, double val=-1.) falcON_THROWING;
    static double trans_value(Model, double val=-1.) falcON_THROWING;
    static double inner_default(Model);
    static double outer_default(Model);
    static double trans_default(Model);
    /// constructor from parameters
    /// \param[in] inner inner power law slope gamma_i of density
    /// \param[in] outer outer power law slope gamma_o of density
    /// \param[in] trans transition steepness eta
    DoublePowerLawHalo(double inner, double outer, double trans);
    /// constructor from model and parameters as required
    /// \param[in] model model
    /// \param[in] inner optional: inner power law slope gamma_i of density
    /// \param[in] outer optional: iouter power law slope gamma_o of density
    /// \param[in] trans optional: itransition steepness eta
    DoublePowerLawHalo(Model  model,
		       double inner=-1., double outer=-1., double trans=-1.)
    falcON_THROWING;
    /// negative logarithmic density slope at r->0
    double inner_gamma() const { return gi; }
    /// scale radius: one
    double scale_radius() const { return 1.; }
    /// transition steepness
    double transition() const { return et; }
    /// no truncation radius
    double trunc_radius() const { return 0.; }
    /// negative logarithmic density slope at r->oo
    double outer_gamma() const { return go; }
    /// density at given (scaled) radius
    double operator()(double x) const;
    /// density and its first derivative
    double operator()(double x, double&rh1) const;
    /// density and its first two derivatives
    double operator()(double x, double&rh1, double&rh2) const;
    /// total mass if multiplied with a truncation factor and used with core
    /// \param M modifier, specifying truncation and core
    double Mtot(const HaloModifier&M) const falcON_THROWING;
  };
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  // class falcON::ModifiedDoublePowerLawHalo                                   
  //                                                                            
  /// implements a halo with truncated double power-law density                 
  //                                                                            
  // ///////////////////////////////////////////////////////////////////////////
  class ModifiedDoublePowerLawHalo : public HaloDensity {
    const DoublePowerLawHalo Model;
    const HaloModifier       Modif;
    const double             rsc,irs,mt,rh0,fc1,fc2;
  public:
    /// constructor from parameters
    /// It is required that the total mass is finite; otherwise, we error out.
    /// The total mass is finite if outer > 3 or trunc > 0.
    /// \param scale scale radius
    /// \param core  core radius
    /// \param trunc truncation radius (zero -> infinity, ie. no truncation)
    /// \param mtot  total mass
    /// \param inner inner power law slope of density
    /// \param outer outer power law slope of density
    /// \param trans transition steepness
    ModifiedDoublePowerLawHalo(double scale, double core, double trunc,
			       double mtot,
			       double inner, double outer, double trans)
      falcON_THROWING :
      Model(inner,outer,trans),
      Modif(core/scale,trunc/scale),
      rsc  (scale),
      irs  (1/scale),
      mt   (mtot),
      rh0  (cube(irs)*mtot/Model.Mtot(Modif)),
      fc1  (rh0*irs),
      fc2  (fc1*irs) 
    {
      if(isinf(rsc)) falcON_THROW("ModifiedDoublePowerLawHalo: r_s=inf\n");
      if(rsc   ==0.) falcON_THROW("ModifiedDoublePowerLawHalo: r_s==0\n");
      if(rsc   < 0.) falcON_THROW("ModifiedDoublePowerLawHalo: r_s=%g<0\n",rsc);
      if(isinf(mt) ) falcON_THROW("ModifiedDoublePowerLawHalo: M=inf\n");
      if(mt    ==0.) falcON_THROW("ModifiedDoublePowerLawHalo: M==0\n");
      if(mt     <0.) falcON_THROW("ModifiedDoublePowerLawHalo: M=%g<0\n",mt);
      DebugInfo(2,"ModifiedDoublePowerLawHalo: rh0=%f\n",rh0);
    }
    /// constructor from model and parameters
    /// It is required that the total mass is finite; otherwise, we error out.
    /// The total mass is finite if outer > 3 or trunc > 0.
    /// \param scale scale radius
    /// \param core  core radius
    /// \param trunc truncation radius (zero -> infinity, ie. no truncation)
    /// \param mtot  total mass
    /// \param inner inner power law slope of density
    /// \param outer outer power law slope of density
    /// \param trans transition steepness
    ModifiedDoublePowerLawHalo(double scale, double core, double trunc,
			       double mtot,
			       DoublePowerLawHalo::Model model,
			       double inner=-1., double outer=-1.,
			       double trans=-1.)
      falcON_THROWING :
      Model(model,inner,outer,trans),
      Modif(core/scale,trunc/scale),
      rsc  (scale),
      irs  (1/scale),
      mt   (mtot),
      rh0  (cube(irs)*mtot/Model.Mtot(Modif)),
      fc1  (rh0*irs),
      fc2  (fc1*irs) 
    {
      if(isinf(rsc)) falcON_THROW("ModifiedDoublePowerLawHalo: r_s=inf\n");
      if(rsc   ==0.) falcON_THROW("ModifiedDoublePowerLawHalo: r_s==0\n");
      if(rsc   < 0.) falcON_THROW("ModifiedDoublePowerLawHalo: r_s=%g<0\n",rsc);
      if(isinf(mt) ) falcON_THROW("ModifiedDoublePowerLawHalo: M=inf\n");
      if(mt    ==0.) falcON_THROW("ModifiedDoublePowerLawHalo: M==0\n");
      if(mt     <0.) falcON_THROW("ModifiedDoublePowerLawHalo: M=%g<0\n",mt);
      DebugInfo(2,"ModifiedDoublePowerLawHalo: rh0=%f\n",rh0);
    }
    /// re-set the model's total mass
    void reset_mass(double Mt) {
      const double fac = Mt/mt;
      if(fac != 1.0) {
	const_cast<double&>(mt)   = Mt;
	const_cast<double&>(rh0) *= fac;
	const_cast<double&>(fc1) *= fac;
	const_cast<double&>(fc2) *= fac;
      }
    }
    /// total mass
    double const&total_mass() const {
      return mt;
    }
    /// density normalisation
    double const&rho0() const {
      return rh0;
    }
    /// negative logarithmic density slope at r->0
    double inner_gamma() const {
      return Modif.cored()? 0. : Model.inner_gamma();
    }
    /// scale radius
    double scale_radius() const {
      return rsc;
    }
    /// core radius
    double core_radius() const {
      return rsc*Modif.r_c();
    }
    /// truncation radius
    double trunc_radius() const {
      return rsc*Modif.r_t();
    }
    /// negative logarithmic density slope at r->oo
    double outer_gamma() const {
      return Modif.truncated()? 100. : Model.outer_gamma();
    }
    /// transition steepness
    double transition() const {
      return Model.transition();
    }
    /// density at given (scaled) radius
    double operator()(double r) const {
      return rh0*Modif(Model,r*irs);
    }
    /// density and its first derivative
    double operator()(double r, double&rh1) const {
      double rho=rh0*Modif(Model,r*irs,rh1);
      rh1 *= fc1;
      return rho;
    }
    /// density and its first two derivatives
    double operator()(double r, double&rh1, double&rh2) const {
      double rho=rh0*Modif(Model,r*irs,rh1,rh2);
      rh1 *= fc1;
      rh2 *= fc2;
      return rho;
    }
  };
////////////////////////////////////////////////////////////////////////////////
} // namespace falcON
falcON_TRAITS(falcON::HaloDensity,"falcON::HaloDensity");
falcON_TRAITS(falcON::HaloModel,"falcON::HaloModel");
falcON_TRAITS(falcON::HaloModifier,"falcON::HaloModifier");
falcON_TRAITS(falcON::DoublePowerLawHalo,"falcON::DoublePowerLawHalo");
falcON_TRAITS(falcON::ModifiedDoublePowerLawHalo,"falcON::ModifiedDoublePowerLawHalo");
////////////////////////////////////////////////////////////////////////////////
#endif
