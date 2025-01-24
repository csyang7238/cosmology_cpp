#include"RecombinationHistory.h"

//====================================================
// Constructors
//====================================================
   
RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo, 
    double Yp) :
  cosmo(cosmo),
  Yp(Yp)
{}

//====================================================
// Do all the solving we need to do
//====================================================

void RecombinationHistory::solve(){
    
  // Compute and spline Xe, ne
  solve_number_density_electrons();
   
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  //solve_for_optical_depth_tau();
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");
  
  //=============================================================================
  // TODO: Set up x-array and make arrays to store X_e(x) and n_e(x) on
  //=============================================================================
  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
  Vector Xe_arr = Vector(npts_rec_arrays);
  Vector ne_arr = Vector(npts_rec_arrays);

  // Calculate recombination history
  bool saha_regime = true;
  for(int i = 0; i < npts_rec_arrays; i++){

    //==============================================================
    // TODO: Get X_e from solving the Saha equation so
    // implement the function electron_fraction_from_saha_equation
    //==============================================================
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit)
      saha_regime = false;

    if(saha_regime){
      
      //=============================================================================
      // TODO: Store the result we got from the Saha equation
      //=============================================================================
      Xe_arr[i] = Xe_current;
      ne_arr[i] = ne_current;

    } else {

      //==============================================================
      // TODO: Compute X_e from current time til today by solving 
      // the Peebles equation (NB: if you solve all in one go remember to
      // exit the for-loop!)
      // Implement rhs_peebles_ode
      //==============================================================
      Vector x_peebles(x_array.begin() + (i-1), x_array.end());
      Vector ne_peebles_arr = Vector(npts_rec_arrays-i+1);

      // The Peebles ODE equation
      ODESolver peebles_Xe_ode;
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        return rhs_peebles_ode(x, Xe, dXedx);
      };
      
      //=============================================================================
      // TODO: Set up IC, solve the ODE and fetch the result 
      //=============================================================================
      Vector Xe_ic{Xe_arr[i-1]};
      peebles_Xe_ode.solve(dXedx, x_peebles, Xe_ic);
      auto Xe_peebles_arr = peebles_Xe_ode.get_data_by_component(0);
      // To calculate the ne from Xe derived from Peebles equation
      for (int j = 0; j < x_peebles.size(); ++j) {
          ne_peebles_arr[j] = Xe_peebles_arr[j] * compute_nb(x_peebles[j]);
      }
      // Fill the vector from the i-th index onwards from ode solution array
      std::copy(Xe_peebles_arr.begin(), Xe_peebles_arr.end(), Xe_arr.begin() + (i-1));
      std::copy(ne_peebles_arr.begin(), ne_peebles_arr.end(), ne_arr.begin() + (i-1));

      break;
    
    }
  }

  //=============================================================================
  // TODO: Spline the result. Implement and make sure the Xe_of_x, ne_of_x 
  // functions are working
  //=============================================================================
  // 
  // Apply logarithm to Xe and ne
  auto log_Xe_arr = log(Xe_arr);
  auto log_ne_arr = log(ne_arr);
  
  //log_Xe_of_x_spline.set_out_of_bounds_warning(true);
  //log_Xe_of_x_spline.create(x_array, log_Xe_arr, "Xe");  //later need to be exponentiated
  //log_ne_of_x_spline.set_out_of_bounds_warning(true);
  //log_ne_of_x_spline.create(x_array, log_ne_arr, "ne");

  Utils::EndTiming("Xe");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  const double a           = exp(x);
 
  // Physical constants
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double epsilon_0   = Constants.epsilon_0;
  const double H0_over_h   = Constants.H0_over_h;

  // Fetch cosmological parameters
  const double OmegaB      = cosmo->get_OmegaB(x);
  const double TCMB        = cosmo->get_TCMB(x);
  const double H_of_x      = cosmo->H_of_x(x);

  // Electron fraction and number density
  double Xe = 0.0;
  double ne = 0.0;
  
  //=============================================================================
  // TODO: Compute Xe and ne from the Saha equation
  //=============================================================================
  const double T_b = TCMB / a;
  const double nb = compute_nb(x);
  const double saha_rhs  = (1.0 / nb) * pow(m_e * k_b * T_b / (2.0 * Constants.pi * pow(hbar, 2.0)), 1.5) * exp(-epsilon_0 / (k_b * T_b));
  
  if ((4.0 / saha_rhs) < 1e-8) {
      Xe = 1.0;
  }
  else {
      Xe = -(saha_rhs / 2.0) * (1.0 - sqrt(1.0 + (4.0 / saha_rhs)));
  }
  
  ne = Xe * nb;

  return std::pair<double,double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of a and X_e
  const double X_e         = Xe[0];
  const double a           = exp(x);

  // Physical constants in SI units
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double c           = Constants.c;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double sigma_T     = Constants.sigma_T;
  const double lambda_2s1s = Constants.lambda_2s1s;
  const double epsilon_0   = Constants.epsilon_0;

  // Cosmological parameters
  const double OmegaB = cosmo->get_OmegaB(x);
  const double TCMB   = cosmo->get_TCMB(x);
  const double H_of_x = cosmo->H_of_x(x);
  const double H0     = cosmo->get_H0();

  //=============================================================================
  // TODO: Write the expression for dXedx
  //=============================================================================
  const double T_b = TCMB / a;
  const double n_H = (1 - Yp) * 3 * pow(H0, 2.0) * OmegaB / (8.0 * Constants.pi * G * m_H * pow(a, 3.0));
  const double n_1s = (1 - X_e) * n_H;
  const double lambda_alpha = H_of_x * (pow(3.0*epsilon_0,3.0)/(pow(8.0*Constants.pi,2.0)*pow(c*hbar,3.0)*n_1s));
  const double phi_2 = 0.448 * log(epsilon_0 / (k_b * T_b));
  const double alpha_2 = (8.0 / sqrt(3.0 * Constants.pi)) * c * sigma_T * (sqrt(epsilon_0 / k_b * T_b)) * phi_2;
  const double beta = alpha_2 * pow((m_e * k_b * T_b) / (2.0 * Constants.pi * pow(hbar, 2.0)), 1.5) * exp(-epsilon_0 / (k_b * T_b));
  const double beta_combined = alpha_2 * pow((m_e * k_b * T_b) / (2.0 * Constants.pi * pow(hbar, 2.0)), 1.5) * exp(-epsilon_0 / (4.0 * k_b * T_b));
  const double Cr = (lambda_2s1s + lambda_alpha)/ (lambda_2s1s + lambda_alpha + beta_combined);
  
  dXedx[0] = (Cr/H_of_x)*(beta*(1.0-X_e)-n_H*alpha_2*pow(X_e,2.0));

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation
  const int npts = 1000;
  Vector x_array_rev = Utils::linspace(x_end, x_start, npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){

    //=============================================================================
    // TODO: Write the expression for dtaudx
    //=============================================================================
    const double sigma_T = Constants.sigma_T;
    const double c = Constants.c;

    const double H_of_x = cosmo->H_of_x(x);

    // Set the derivative for photon optical depth
    dtaudx[0] = -c*ne_of_x(x)*sigma_T/H_of_x;

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set up and solve the ODE and make tau splines
  //=============================================================================
  double tau_ini = 0.0;
  Vector tau_ic{ tau_ini };
  ODESolver tau_ode;
  tau_ode.solve(dtaudx, x_array_rev, tau_ic);
  auto tau_arr = tau_ode.get_data_by_component(0);
  //tau_of_x_spline.create(x_array_rev, tau_arr, "tau");

  //=============================================================================
  // TODO: Compute visibility functions and spline everything
  //=============================================================================
  Vector g_tilde_arr(npts);
  for (int i = 0; i < g_tilde_arr.size(); ++i) {
      g_tilde_arr[i] = -dtaudx_of_x(x_array_rev[i]) * exp(-tau_of_x(x_array_rev[i]));
  }
  //g_tilde_of_x_spline.create(x_array_rev, g_tilde_arr, "g");

  Utils::EndTiming("opticaldepth");
}

//====================================================
// Get methods
//====================================================
double RecombinationHistory::compute_nb(double x) const {
    const double a      = exp(x);
    const double G      = Constants.G;
    const double m_H    = Constants.m_H;
    const double OmegaB = cosmo->get_OmegaB(x);
    const double H0     = cosmo->get_H0();
    
    return (1 - Yp) * 3 * pow(H0, 2.0) * OmegaB / (8 * Constants.pi * G * m_H * pow(a, 3.0));
}

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement. Either from the tau-spline tau_of_x_spline.deriv_(x) or 
  // from a separate spline if you choose to do this
  //=============================================================================

  return tau_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================

  return tau_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================

  return g_tilde_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================

  return g_tilde_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::Xe_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================

  return exp(log_Xe_of_x_spline(x));
}

double RecombinationHistory::ne_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================

  return exp(log_ne_of_x_spline(x));
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = 5000;
  const double x_min   = x_start;
  const double x_max   = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                    << " ";
    fp << compute_nb(x)        << " ";
    fp << Xe_of_x(x)           << " ";
    fp << ne_of_x(x)           << " ";
    fp << tau_of_x(x)          << " ";
    fp << dtaudx_of_x(x)       << " ";
    fp << ddtauddx_of_x(x)     << " ";
    fp << g_tilde_of_x(x)      << " ";
    fp << dgdx_tilde_of_x(x)   << " ";
    fp << ddgddx_tilde_of_x(x) << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

