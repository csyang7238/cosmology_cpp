#include "BackgroundCosmology.h"

//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double OmegaB, 
    double OmegaCDM, 
    double OmegaK,
    double Neff, 
    double TCMB) :
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
  OmegaK(OmegaK),
  Neff(Neff), 
  TCMB(TCMB)
{

  //=============================================================================
  // TODO: Compute OmegaR, OmegaNu, OmegaLambda, H0, ...
  //=============================================================================

  H0          = h * Constants.H0_over_h;
  OmegaR      = 2 * (pow(Constants.pi, 2.0) / 30) * (pow(Constants.k_b * TCMB, 4.0) / (pow(Constants.hbar, 3.0) * pow(Constants.c, 5.0))) * (8 * Constants.pi * Constants.G / (3*pow(H0,2.0)));
  OmegaNu     = Neff * OmegaR * 0.875 * (pow(0.3636,1.3333));
  OmegaLambda = 1 - OmegaB - OmegaCDM - OmegaK - OmegaR - OmegaNu;
}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(){
  Utils::StartTiming("Eta");
    
  //=============================================================================
  // TODO: Set the range of x and the number of points for the splines
  // For this Utils::linspace(x_start, x_end, npts) is useful
  //=============================================================================
  const int n_pts = 100;
  Vector x_array = Utils::linspace(Constants.x_start, Constants.x_end, n_pts);

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){

    //=============================================================================
    // TODO: Set the rhs of the detadx and detadx ODE 
    //=============================================================================

    detadx[0] = Constants.c/Hp_of_x(x);

    return GSL_SUCCESS;
  };

  // The ODE for dt/dx
  ODEFunction dtdx = [&](double x, const double *t, double *dtdx) {

    //=============================================================================
    // TODO: Set the rhs of the detadx and dtdx ODE 
    //=============================================================================

    dtdx[0] = 1 / H_of_x(x);

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set the initial condition, set up the ODE system, solve and make
  // the spline eta_of_x_spline 
  //=============================================================================
 
  // Set the IC 
  // for eta
  double eta_ini = Constants.c / Hp_of_x(Constants.x_start);
  Vector eta_ic{ eta_ini };

  // for time
  double t_ini = 1 / (2*H_of_x(Constants.x_start));
  Vector t_ic{ t_ini };

  // Solve the ODE
  ODESolver ode;
  ODESolver ode2;
  ode.solve(detadx, x_array, eta_ic);
  ode2.solve(dtdx, x_array, t_ic);

  // Get the solution (we only have one component so index 0 holds the solution)
  auto eta_array = ode.get_data_by_component(0);
  auto t_array = ode2.get_data_by_component(0);

  eta_of_x_spline.set_out_of_bounds_warning(true);
  eta_of_x_spline.create(x_array, eta_array, "eta_of_x");

  t_of_x_spline.set_out_of_bounds_warning(true);
  t_of_x_spline.create(x_array, t_array, "t_of_x");

  Utils::EndTiming("Eta");
}

//====================================================
// Get methods
//====================================================

double BackgroundCosmology::H_of_x(double x) const{

  //=============================================================================================
  // TODO: Calc. Hubble parameter based on the present time density parameters and scale factors
  //=============================================================================================

  return H0 * sqrt((OmegaB+OmegaCDM)*exp(-3*x) + (OmegaR+OmegaNu)*exp(-4*x) + OmegaK*exp(-2*x) + OmegaLambda);
}

double BackgroundCosmology::Hp_of_x(double x) const{

  //=============================================================================
  // TODO: Calc. scaled Hubble parameter 
  //=============================================================================

  return exp(x) * H_of_x(x);
}

double BackgroundCosmology::dHpdx_of_x(double x) const{

  //=============================================================================
  // TODO: Calc. first derivative of scaled Hubble parameter.
  //=============================================================================

  return Hp_of_x(x) + (exp(3*x)/Hp_of_x(x)) * (pow(H0,2.0)/2) * (-3 * (OmegaB + OmegaCDM) * exp(-4 * x) - 4 * (OmegaR + OmegaNu) * exp(-5 * x) - 2 * OmegaK * exp(-3 * x));
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{

  //=============================================================================
  // TODO: Calc. second derivative of scaled Hubble parameter.
  //=============================================================================

  return 5 * dHpdx_of_x(x) - (1 / Hp_of_x(x)) * (pow(dHpdx_of_x(x), 2.0)) - 3 * Hp_of_x(x) + (exp(4 * x) / Hp_of_x(x)) * (pow(H0, 2.0) / 2) * (12 * (OmegaB + OmegaCDM) * exp(-5 * x) + 20 * (OmegaR + OmegaNu) * exp(-6 * x) + 6 * OmegaK * exp(-4 * x));
}

double BackgroundCosmology::get_OmegaB(double x) const{ 
  if(x == 0.0) return OmegaB;

  //=============================================================================
  // TODO: Calc. density parameter (baryon) at a certain x.
  //=============================================================================

  return OmegaB / (exp(3 * x) * pow(H_of_x(x)/H0, 2.0));
}

double BackgroundCosmology::get_OmegaR(double x) const{ 
  if(x == 0.0) return OmegaR;

  //=============================================================================
  // TODO: Calc. density parameter (radiation) at a certain x.
  //=============================================================================
  
  return OmegaR / (exp(4 * x) * pow(H_of_x(x) / H0, 2.0));
}

double BackgroundCosmology::get_OmegaNu(double x) const{ 
  if(x == 0.0) return OmegaNu;

  //=============================================================================
  // TODO: Calc. density parameter (neutrino) at a certain x.
  //=============================================================================

  return OmegaNu / (exp(4 * x) * pow(H_of_x(x) / H0, 2.0));
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  if(x == 0.0) return OmegaCDM;

  //=============================================================================
  // TODO: Calc. density parameter (CDM) at a certain x.
  //=============================================================================

  return OmegaCDM / (exp(3 * x) * pow(H_of_x(x) / H0, 2.0));
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  if(x == 0.0) return OmegaLambda;

  //=============================================================================
  // TODO: Calc. density parameter (cosmological constant) at a certain x.
  //=============================================================================

  return OmegaLambda / pow(H_of_x(x) / H0, 2.0);
}

double BackgroundCosmology::get_OmegaK(double x) const{ 
  if(x == 0.0) return OmegaK;

  //=============================================================================
  // TODO: Calc. density parameter (curvature) at a certain x.
  //=============================================================================

  return OmegaK / (exp(2 * x) * pow(H_of_x(x) / H0, 2.0));
}
    
double BackgroundCosmology::get_luminosity_distance_of_x(double x) const{
  //=============================================================================
  // TODO: Implement...
  //=============================================================================

  return get_comoving_distance_of_x(x)/exp(x);
}
double BackgroundCosmology::get_comoving_distance_of_x(double x) const{
  //=============================================================================
  // TODO: Implement...
  //=============================================================================

  return eta_of_x(0.0) - eta_of_x(x);
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::t_of_x(double x) const {
  return t_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const{ 
  return H0; 
}

double BackgroundCosmology::get_h() const{ 
  return h; 
}

double BackgroundCosmology::get_Neff() const{ 
  return Neff; 
}

double BackgroundCosmology::get_TCMB(double x) const{ 
  if(x == 0.0) return TCMB;
  return TCMB * exp(-x); 
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:                     " << OmegaB                         << "\n";
  std::cout << "OmegaCDM:                   " << OmegaCDM                       << "\n";
  std::cout << "OmegaLambda:                " << OmegaLambda                    << "\n";
  std::cout << "OmegaK:                     " << OmegaK                         << "\n";
  std::cout << "OmegaNu:                    " << OmegaNu                        << "\n";
  std::cout << "OmegaR:                     " << OmegaR                         << "\n";
  std::cout << "Neff:                       " << Neff                           << "\n";
  std::cout << "h:                          " << h                              << "\n";
  std::cout << "TCMB:                       " << TCMB                           << "\n";
  std::cout << "Age of the Universe (Gy):   " << t_of_x(0)/(1e9*365*24*3600)    << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = -15.0;
  const double x_max =  5.0;
  const int    n_pts =  100;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";
    fp << eta_of_x(x)        << " ";
    fp << t_of_x(x)          << " ";
    fp << Hp_of_x(x)         << " ";
    fp << dHpdx_of_x(x)      << " ";
    fp << ddHpddx_of_x(x)    << " ";
    fp << get_OmegaB(x)      << " ";
    fp << get_OmegaCDM(x)    << " ";
    fp << get_OmegaLambda(x) << " ";
    fp << get_OmegaR(x)      << " ";
    fp << get_OmegaNu(x)     << " ";
    fp << get_OmegaK(x)      << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

