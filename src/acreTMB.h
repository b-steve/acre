#ifndef acreTMB_hpp
#define acreTMB_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

#include <map>

//DIY functions
template<class Type>
Type bessi0(Type x){
  Type ax;
  Type ans;
  Type y;
  if ((ax = fabs(x)) < 3.75){
    y = x/3.75;
    y *= y;
    ans = 1.0 + y*(3.5156229 + y*(3.0899424 + y*(1.2067492 + y*(0.2659732 + y*(0.360768e-1 + y*0.45813e-2)))));
  } else {
    y = 3.75/ax;
    ans = (exp(ax)/sqrt(ax))*(0.39894228 + y*(0.1328592e-1 + y*(0.225319e-2 + 
      y*(-0.157565e-2 + y*(0.916281e-2 + y*(-0.2057706e-1 + y*(0.2635537e-1 + 
      y*(-0.1647633e-1 + y*0.392377e-2))))))));
  }
  return ans;
}


template<class Type>
Type Gamma(Type z){
  if(z < 0.5){
    z = 1 - z;
    Type y = sqrt(2 * 3.14159265358979323846) * pow(z + 6.5, (z - 0.5)) * exp((-1) * z - 6.5) * 
      (0.99999999999980993 + 676.5203681218851/(z) + (-1259.1392167224028)/(z + 1) + 771.32342877765313/(z + 2) +
      (-176.61502916214059)/(z + 3) + 12.507343278686905/(z + 4) + (-0.13857109526572012)/(z + 5));
    //+ 9.9843695780195716e-6/(z + 6) + 1.5056327351493116e-7/(z + 7));
    
    return 3.14159265358979323846 / (sin(3.14159265358979323846 * (z - 1)) * y);
  } else {
    Type y = sqrt(2 * 3.14159265358979323846) * pow(z + 6.5, (z - 0.5)) * exp((-1) * z - 6.5) * 
      (0.99999999999980993 + 676.5203681218851/(z) + (-1259.1392167224028)/(z + 1) + 771.32342877765313/(z + 2) +
      (-176.61502916214059)/(z + 3) + 12.507343278686905/(z + 4) + (-0.13857109526572012)/(z + 5)); 
    //+ 9.9843695780195716e-6/(z + 6) + 1.5056327351493116e-7/(z + 7));
    
    return y;
  }
}

int incheck_scalar(const int &a, const vector<int> &b){
  int lengthb = b.size();
  int ans = 0;
  for(int i = 0; i < lengthb; i++){
    if(a == b(i)){
      ans = 1;
      break;
    }
  }
  return ans;
}


template<class Type>
void trans(Type *x, const int &link){
  
  if(link == 2){
    *x = exp(*x);
  }
  if(link == 3){
    *x = exp(*x)/(1 + exp(*x));
  }
  
}

//The linear predictors for one parameter, on the link scale: at each session-trap row, from
//DX, and at each mask point, from DX_mask. The coefficients in beta are those for the columns
//of DX followed by those for DX_mask, with numbers of columns in row k of par_n_col. A
//coefficient outside its bounds makes the likelihood zero.
template<class Type>
void linear_predictor(const vector<Type> &beta, const matrix<Type> &DX, const matrix<Type> &DX_mask,
                      const matrix<Type> &bound, const matrix<int> &par_n_col, int k, int nrow_mask,
                      Type &nll, vector<Type> &vec_full, vector<Type> &vec_mask){
  int n_full = par_n_col(k, 0);
  int n_mask = par_n_col(k, 1);
  for(int i = 0; i < n_full + n_mask; i++){
    if(beta(i) < bound(0, i) || beta(i) > bound(1, i)) nll += std::numeric_limits<double>::infinity();
  }
  vector<Type> beta_full = beta.head(n_full);
  vec_full = DX * beta_full;
  vec_mask.resize(nrow_mask);
  if(n_mask == 0){
    vec_mask.setZero();
  } else {
    vector<Type> beta_mask = beta.tail(n_mask);
    vec_mask = DX_mask * beta_mask;
  }
}

//the value of x back-transformed by the link function
template<class Type>
Type linked(Type x, const int &link){
  trans(&x, link);
  return x;
}







//detect functions---------------------------------------------------------------------

//hn
template<class Type>
Type det_hn(const Type &dx, const vector<Type> &param){
  Type ans = 0.0;
  Type g0 = param(0);
  Type sigma = param(1);
  ans = g0 * exp(-1 * pow(dx, 2) / (2 * pow(sigma, 2)));
  return ans;
}

//hhn
template<class Type>
Type det_hhn(const Type &dx, const vector<Type> &param){
  Type ans = 0.0;
  Type lambda0 = param(0);
  Type sigma = param(1);
  Type lambda_d = lambda0 * exp(pow(dx, 2) * (-0.5) / pow(sigma, 2));
  ans = 1 - exp((-1) * lambda_d);
  return ans;
}

//hr
template<class Type>
Type det_hr(const Type &dx, const vector<Type> &param){
  Type ans = 0.0;
  Type g0 = param(0);
  Type sigma = param(1);
  Type z = param(2);
  ans = g0 * (1 - exp((-1) * pow(dx / sigma, (-z))));
  return ans;
}

//lth
template<class Type>
Type det_lth(const Type &dx, const vector<Type> &param){
  Type ans = 0.0;
  Type shape_1 = param(0);
  Type shape_2 = param(1);
  Type scale = param(2);
  ans = 0.5 - 0.5 * erf(shape_1 - exp(shape_2 - scale * dx));
  return ans;
}

//th
template<class Type>
Type det_th(const Type &dx, const vector<Type> &param){
  Type ans = 0.0;
  Type shape = param(0);
  Type scale = param(1);
  ans = 0.5 - 0.5 * erf(dx / scale - shape);
  return ans;
}

//ss
//in ss, we do something different, here we use
//the same input, but output is back-transformed essx
//instead of the probability be detected
template<class Type>
Type essx_ss_identical(const Type &dx, const vector<Type> &param){
  Type b0_ss = param(0);
  Type b1_ss = param(1);
  Type essx = b0_ss - b1_ss * dx;
  return essx;
}

template<class Type>
Type essx_ss_log(const Type &dx, const vector<Type> &param){
  Type b0_ss = param(0);
  Type b1_ss = param(1);
  Type essx = exp(b0_ss - b1_ss * dx);
  return essx;
}

template<class Type>
Type essx_ss_spherical(const Type &dx, const vector<Type> &param){
  Type b0_ss = param(0);
  Type b1_ss = param(1);
  Type essx = 0.0;
  if(dx > 1){
    essx += b0_ss - 20 * log10(dx) - b1_ss * (dx - 1);
  } else {
    essx += b0_ss;
  }
  return essx;
}

//ss_dir
//currently not in use
template<class Type>
Type essx_ss_dir_identical(const Type &dx, const vector<Type> &param){
  Type b0_ss = param(0);
  Type b1_ss = param(1);
  Type b2_ss = param(2);
  Type theta = param(3);
  Type essx = b0_ss - (b1_ss - b2_ss * (cos(theta) - 1)) * dx;
  return essx;
}

template<class Type>
Type essx_ss_dir_log(const Type &dx, const vector<Type> &param){
  Type b0_ss = param(0);
  Type b1_ss = param(1);
  Type b2_ss = param(2);
  Type theta = param(3);
  Type essx = exp(b0_ss - (b1_ss - b2_ss * (cos(theta) - 1)) * dx);
  return essx;
}

template<class Type>
Type essx_ss_dir_spherical(const Type &dx, const vector<Type> &param){
  Type b0_ss = param(0);
  Type b1_ss = param(1);
  Type b2_ss = param(2);
  Type theta = param(3);
  Type essx = 0.0;
  if(dx > 1){
    essx += b0_ss;
  } else {
    essx += b0_ss - 20 * log10(dx) - (b1_ss - b2_ss * (cos(theta) - 1)) * (dx - 1);
  }
  return essx;
}

//ss_het
//currently not in use
template<class Type>
Type essx_ss_het(const Type &dx, const vector<Type> &param){
  Type b0_ss = param(0);
  Type b1_ss = param(1);
  Type u = param(2);
  Type essx = b0_ss - b1_ss * dx + u;
  return essx;
}

//end of detection functions------------------------------------------------------------------------

template<class Type>
Type acreTMB(objective_function<Type>* obj)
{
  Type nll = Type(0.0);
  DATA_INTEGER(n_sessions);
  DATA_IVECTOR(n_animals);
  DATA_IVECTOR(n_IDs);
  DATA_IVECTOR(n_traps);
  DATA_IVECTOR(n_masks);
  DATA_IVECTOR(n_calls_each_animal);
  
  DATA_INTEGER(nrow_data_mask);
  DATA_INTEGER(nrow_dx);
  
  DATA_VECTOR(A);
  DATA_VECTOR(survey_length);
  DATA_SCALAR(cue_rates);
  
  // Conditional likelihood 
  DATA_INTEGER(is_conditional);
  DATA_INTEGER(is_two_stage);
  
  DATA_SCALAR(cutoff);
  
  
  DATA_INTEGER(detfn_index);
  DATA_IVECTOR(param_og);
  DATA_IMATRIX(par_n_col);
  DATA_IVECTOR(par_link);
  DATA_INTEGER(ss_link);
  
  
  DATA_INTEGER(is_animalID);
  DATA_INTEGER(is_ss);
  DATA_INTEGER(is_bearing);
  DATA_INTEGER(is_toa);
  DATA_INTEGER(is_dist);
  DATA_INTEGER(is_local);
  
  
  DATA_VECTOR(capt_bin);
  DATA_VECTOR(capt_bearing);
  DATA_VECTOR(capt_dist);
  DATA_VECTOR(capt_ss);
  DATA_VECTOR(dx);
  DATA_VECTOR(theta);
  DATA_VECTOR(index_local);
  DATA_VECTOR(toa_ssq);
  
  
  DATA_MATRIX(g0_DX);
  DATA_MATRIX(sigma_DX);
  DATA_MATRIX(lambda0_DX);
  DATA_MATRIX(z_DX);
  DATA_MATRIX(shape_1_DX);
  DATA_MATRIX(shape_2_DX);
  DATA_MATRIX(shape_DX);
  DATA_MATRIX(scale_DX);
  DATA_MATRIX(b0_ss_DX);
  DATA_MATRIX(b1_ss_DX);
  DATA_MATRIX(b2_ss_DX);
  DATA_MATRIX(sigma_ss_DX);
  DATA_MATRIX(kappa_DX);
  DATA_MATRIX(alpha_DX);
  DATA_MATRIX(sigma_toa_DX);
  //sigma_b0_ss is not extentable
  DATA_MATRIX(D_DX);
  DATA_MATRIX(mu_DX);
  
  
  DATA_MATRIX(g0_DX_mask);
  DATA_MATRIX(sigma_DX_mask);
  DATA_MATRIX(lambda0_DX_mask);
  DATA_MATRIX(z_DX_mask);
  DATA_MATRIX(shape_1_DX_mask);
  DATA_MATRIX(shape_2_DX_mask);
  DATA_MATRIX(shape_DX_mask);
  DATA_MATRIX(scale_DX_mask);
  DATA_MATRIX(b0_ss_DX_mask);
  DATA_MATRIX(b1_ss_DX_mask);
  DATA_MATRIX(b2_ss_DX_mask);
  DATA_MATRIX(sigma_ss_DX_mask);
  DATA_MATRIX(kappa_DX_mask);
  DATA_MATRIX(alpha_DX_mask);
  DATA_MATRIX(sigma_toa_DX_mask);
  //sigma_b0_ss is not extentable
  DATA_MATRIX(D_DX_mask);
  DATA_MATRIX(mu_DX_mask);
  
  DATA_MATRIX(g0_bound);
  DATA_MATRIX(sigma_bound);
  DATA_MATRIX(lambda0_bound);
  DATA_MATRIX(z_bound);
  DATA_MATRIX(shape_1_bound);
  DATA_MATRIX(shape_2_bound);
  DATA_MATRIX(shape_bound);
  DATA_MATRIX(scale_bound);
  DATA_MATRIX(b0_ss_bound);
  DATA_MATRIX(b1_ss_bound);
  DATA_MATRIX(b2_ss_bound);
  DATA_MATRIX(sigma_ss_bound);
  DATA_MATRIX(kappa_bound);
  DATA_MATRIX(alpha_bound);
  DATA_MATRIX(sigma_toa_bound);
  DATA_MATRIX(sigma_b0_ss_bound);
  DATA_MATRIX(D_bound);
  DATA_MATRIX(mu_bound);
  

  //define a pointer variable which point to a detect function
  //based on our input "detfn_index"
  //and also point out the number of parameters for that detfn
  
  Type (*detfn)(const Type &dx, const vector<Type> &param);
  int n_detfn_param;
  if(detfn_index == 1){
    detfn = det_hn;
    n_detfn_param = 2;
  } else if(detfn_index == 2){
    detfn = det_hhn;
    n_detfn_param = 2;
  } else if(detfn_index == 3){
    detfn = det_hr;
    n_detfn_param = 3;
  } else if(detfn_index == 5){
    detfn = det_lth;
    n_detfn_param = 3;
  } else if(detfn_index == 4){
    detfn = det_th;
    n_detfn_param = 2;
  } else if(detfn_index == 6){
    if(ss_link == 1){
      detfn = essx_ss_identical;
    } else if(ss_link == 2){
      detfn = essx_ss_log;
    } else if(ss_link == 4){
      detfn = essx_ss_spherical;
    }
    n_detfn_param = 2;
  }
  //define the parameter vector which will be used later
  vector<Type> detfn_param(n_detfn_param);
  
  //Each parameter's coefficients are those for the columns of its design matrix for session-trap
  //rows (*_DX) followed by those for its mask-level design matrix (*_DX_mask), on the scale of
  //its link function. Parameters not used by the model are fixed in MakeADFun() via 'map'.
  PARAMETER_VECTOR(g0);
  vector<Type> g0_vec_full, g0_vec_mask;
  linear_predictor(g0, g0_DX, g0_DX_mask, g0_bound, par_n_col, 0, nrow_data_mask, nll, g0_vec_full, g0_vec_mask);
  if(incheck_scalar(1, param_og) == 1) ADREPORT(g0);
  PARAMETER_VECTOR(sigma);
  vector<Type> sigma_vec_full, sigma_vec_mask;
  linear_predictor(sigma, sigma_DX, sigma_DX_mask, sigma_bound, par_n_col, 1, nrow_data_mask, nll, sigma_vec_full, sigma_vec_mask);
  if(incheck_scalar(2, param_og) == 1) ADREPORT(sigma);
  PARAMETER_VECTOR(lambda0);
  vector<Type> lambda0_vec_full, lambda0_vec_mask;
  linear_predictor(lambda0, lambda0_DX, lambda0_DX_mask, lambda0_bound, par_n_col, 2, nrow_data_mask, nll, lambda0_vec_full, lambda0_vec_mask);
  if(incheck_scalar(3, param_og) == 1) ADREPORT(lambda0);
  PARAMETER_VECTOR(z);
  vector<Type> z_vec_full, z_vec_mask;
  linear_predictor(z, z_DX, z_DX_mask, z_bound, par_n_col, 3, nrow_data_mask, nll, z_vec_full, z_vec_mask);
  if(incheck_scalar(4, param_og) == 1) ADREPORT(z);
  PARAMETER_VECTOR(shape_1);
  vector<Type> shape_1_vec_full, shape_1_vec_mask;
  linear_predictor(shape_1, shape_1_DX, shape_1_DX_mask, shape_1_bound, par_n_col, 4, nrow_data_mask, nll, shape_1_vec_full, shape_1_vec_mask);
  if(incheck_scalar(5, param_og) == 1) ADREPORT(shape_1);
  PARAMETER_VECTOR(shape_2);
  vector<Type> shape_2_vec_full, shape_2_vec_mask;
  linear_predictor(shape_2, shape_2_DX, shape_2_DX_mask, shape_2_bound, par_n_col, 5, nrow_data_mask, nll, shape_2_vec_full, shape_2_vec_mask);
  if(incheck_scalar(6, param_og) == 1) ADREPORT(shape_2);
  PARAMETER_VECTOR(shape);
  vector<Type> shape_vec_full, shape_vec_mask;
  linear_predictor(shape, shape_DX, shape_DX_mask, shape_bound, par_n_col, 6, nrow_data_mask, nll, shape_vec_full, shape_vec_mask);
  if(incheck_scalar(7, param_og) == 1) ADREPORT(shape);
  PARAMETER_VECTOR(scale);
  vector<Type> scale_vec_full, scale_vec_mask;
  linear_predictor(scale, scale_DX, scale_DX_mask, scale_bound, par_n_col, 7, nrow_data_mask, nll, scale_vec_full, scale_vec_mask);
  if(incheck_scalar(8, param_og) == 1) ADREPORT(scale);
  PARAMETER_VECTOR(b0_ss);
  vector<Type> b0_ss_vec_full, b0_ss_vec_mask;
  linear_predictor(b0_ss, b0_ss_DX, b0_ss_DX_mask, b0_ss_bound, par_n_col, 8, nrow_data_mask, nll, b0_ss_vec_full, b0_ss_vec_mask);
  if(incheck_scalar(9, param_og) == 1) ADREPORT(b0_ss);
  PARAMETER_VECTOR(b1_ss);
  vector<Type> b1_ss_vec_full, b1_ss_vec_mask;
  linear_predictor(b1_ss, b1_ss_DX, b1_ss_DX_mask, b1_ss_bound, par_n_col, 9, nrow_data_mask, nll, b1_ss_vec_full, b1_ss_vec_mask);
  if(incheck_scalar(10, param_og) == 1) ADREPORT(b1_ss);
  PARAMETER_VECTOR(b2_ss);
  vector<Type> b2_ss_vec_full, b2_ss_vec_mask;
  linear_predictor(b2_ss, b2_ss_DX, b2_ss_DX_mask, b2_ss_bound, par_n_col, 10, nrow_data_mask, nll, b2_ss_vec_full, b2_ss_vec_mask);
  if(incheck_scalar(11, param_og) == 1) ADREPORT(b2_ss);
  PARAMETER_VECTOR(sigma_ss);
  vector<Type> sigma_ss_vec_full, sigma_ss_vec_mask;
  linear_predictor(sigma_ss, sigma_ss_DX, sigma_ss_DX_mask, sigma_ss_bound, par_n_col, 11, nrow_data_mask, nll, sigma_ss_vec_full, sigma_ss_vec_mask);
  if(incheck_scalar(12, param_og) == 1) ADREPORT(sigma_ss);
  PARAMETER_VECTOR(kappa);
  vector<Type> kappa_vec_full, kappa_vec_mask;
  linear_predictor(kappa, kappa_DX, kappa_DX_mask, kappa_bound, par_n_col, 12, nrow_data_mask, nll, kappa_vec_full, kappa_vec_mask);
  if(incheck_scalar(13, param_og) == 1) ADREPORT(kappa);
  PARAMETER_VECTOR(alpha);
  vector<Type> alpha_vec_full, alpha_vec_mask;
  linear_predictor(alpha, alpha_DX, alpha_DX_mask, alpha_bound, par_n_col, 13, nrow_data_mask, nll, alpha_vec_full, alpha_vec_mask);
  if(incheck_scalar(14, param_og) == 1) ADREPORT(alpha);
  PARAMETER_VECTOR(sigma_toa);
  vector<Type> sigma_toa_vec_full, sigma_toa_vec_mask;
  linear_predictor(sigma_toa, sigma_toa_DX, sigma_toa_DX_mask, sigma_toa_bound, par_n_col, 14, nrow_data_mask, nll, sigma_toa_vec_full, sigma_toa_vec_mask);
  if(incheck_scalar(15, param_og) == 1) ADREPORT(sigma_toa);
  //sigma_b0_ss is not extendable, so it is a scalar
  PARAMETER(sigma_b0_ss);
  if(sigma_b0_ss < sigma_b0_ss_bound(0, 0) || sigma_b0_ss > sigma_b0_ss_bound(1, 0)) nll += std::numeric_limits<double>::infinity();
  if(incheck_scalar(16, param_og) == 1) ADREPORT(sigma_b0_ss);
  //D is used by every model, so it is always reported
  PARAMETER_VECTOR(D);
  vector<Type> D_vec_full, D_vec_mask;
  ADREPORT(D);
  linear_predictor(D, D_DX, D_DX_mask, D_bound, par_n_col, 16, nrow_data_mask, nll, D_vec_full, D_vec_mask);
  //mu is the call rate in the animal model
  PARAMETER_VECTOR(mu);
  vector<Type> mu_vec_full, mu_vec_mask;
  linear_predictor(mu, mu_DX, mu_DX_mask, mu_bound, par_n_col, 17, nrow_data_mask, nll, mu_vec_full, mu_vec_mask);
  if(incheck_scalar(18, param_og) == 1) ADREPORT(mu);
  //a latent variable, currently unused
  PARAMETER_VECTOR(u);
  
  //---------------------------------------------------------------------------------------
  //The likelihood.
  //
  //Both models are written in terms of groups of detected calls. Under the animal model
  //(Stevenson et al., 2021), a group is a detected animal and its detected calls; under the
  //call model, each detected call is a group of its own. Each group contributes an integral
  //over mask points of a group-level factor (the density and, for the animal model, the
  //Poisson probability of its number of detected calls) times, for each of its calls, the
  //probability of the call's binary capture history and the density of its auxiliary data.
  //The models differ only in the expected number of detected groups (lambda_theta), the
  //group-level factor, and the final normalisation.
  //
  //Data layouts, by row:
  //  capture data (capt_*): session - call - trap (calls ordered by group)
  //  *_DX, so *_vec_full:   session - trap (there are no call-level covariates)
  //  data_mask:             session - mask
  //  dx, theta and essx:    session - mask - trap
  //  data_IDmask:           session - call - mask
  //---------------------------------------------------------------------------------------

  vector<Type> essx(nrow_dx);
  essx.setZero();
  vector<Type> esa(n_sessions);

  //first rows of the current session in each data layout
  int off_full = 0;
  int off_trap = 0;
  int off_mask = 0;
  int off_dt = 0;
  int off_idmask = 0;
  //first element of n_calls_each_animal for the current session
  int off_group = 0;

  //with no auxiliary data, calls in the call model that share a binary capture history make
  //identical contributions, so each unique history is evaluated once
  int no_aux = (is_toa == 0 && is_bearing == 0 && is_dist == 0 && is_ss == 0);
  int share_histories = (is_animalID == 0 && no_aux == 1 && is_local == 0);

  for(int s = 0; s < n_sessions; s++){
    int n_m = n_masks(s);
    int n_t = n_traps(s);
    int n_c = n_IDs(s);
    int n_g = (is_animalID == 1) ? n_animals(s) : n_c;
    Type area_unit = A(s);
    Type survey_len = survey_length(s);

    //density, call rate and detection probabilities at each mask point; D, mu and
    //sigma_toa have no trap-level covariates, so they come from the session's first trap
    vector<Type> D_m(n_m);
    vector<Type> mu_m(n_m);
    vector<Type> p_dot(n_m);
    matrix<Type> p_k(n_m, n_t);
    Type lambda_theta = Type(0.0);
    esa(s) = Type(0.0);
    for(int m = 0; m < n_m; m++){
      int mrow = off_mask + m;
      D_m(m) = linked(D_vec_full(off_trap) + D_vec_mask(mrow), par_link(16));
      p_dot(m) = Type(1.0);
      for(int t = 0; t < n_t; t++){
        int row = off_trap + t;
        int dt = off_dt + m * n_t + t;
        if(detfn_index == 1){
          detfn_param(0) = linked(g0_vec_full(row) + g0_vec_mask(mrow), par_link(0));
          detfn_param(1) = linked(sigma_vec_full(row) + sigma_vec_mask(mrow), par_link(1));
        } else if(detfn_index == 2){
          detfn_param(0) = linked(lambda0_vec_full(row) + lambda0_vec_mask(mrow), par_link(2));
          detfn_param(1) = linked(sigma_vec_full(row) + sigma_vec_mask(mrow), par_link(1));
        } else if(detfn_index == 3){
          detfn_param(0) = linked(g0_vec_full(row) + g0_vec_mask(mrow), par_link(0));
          detfn_param(1) = linked(sigma_vec_full(row) + sigma_vec_mask(mrow), par_link(1));
          detfn_param(2) = linked(z_vec_full(row) + z_vec_mask(mrow), par_link(3));
        } else if(detfn_index == 4){
          detfn_param(0) = linked(shape_vec_full(row) + shape_vec_mask(mrow), par_link(6));
          detfn_param(1) = linked(scale_vec_full(row) + scale_vec_mask(mrow), par_link(7));
        } else if(detfn_index == 5){
          detfn_param(0) = linked(shape_1_vec_full(row) + shape_1_vec_mask(mrow), par_link(4));
          detfn_param(1) = linked(shape_2_vec_full(row) + shape_2_vec_mask(mrow), par_link(5));
          detfn_param(2) = linked(scale_vec_full(row) + scale_vec_mask(mrow), par_link(7));
        } else if(detfn_index == 6){
          detfn_param(0) = linked(b0_ss_vec_full(row) + b0_ss_vec_mask(mrow), par_link(8));
          detfn_param(1) = linked(b1_ss_vec_full(row) + b1_ss_vec_mask(mrow), par_link(9));
        }

        if(is_ss == 0){
          p_k(m, t) = (*detfn)(dx(dt), detfn_param);
        } else {
          //for signal strengths, detfn gives the expected received signal strength
          Type sigma_ss_tem = linked(sigma_ss_vec_full(row) + sigma_ss_vec_mask(mrow), par_link(11));
          essx(dt) = (*detfn)(dx(dt), detfn_param);
          p_k(m, t) = 1 - pnorm((cutoff - essx(dt)) / sigma_ss_tem);
        }
        p_dot(m) *= 1 - p_k(m, t);
      }
      p_dot(m) = 1 - p_dot(m);

      //probability that a group at m is detected: for a call, p_dot; for an animal, that
      //at least one of its calls is detected
      if(is_animalID == 0){
        esa(s) += p_dot(m);
        lambda_theta += D_m(m) * p_dot(m);
      } else {
        mu_m(m) = linked(mu_vec_full(off_trap) + mu_vec_mask(mrow), par_link(17));
        Type p_animal = 1 - exp((-1) * mu_m(m) * survey_len * p_dot(m));
        esa(s) += p_animal;
        lambda_theta += D_m(m) * p_animal;
      }
    }

    esa(s) *= area_unit;
    if(is_animalID == 0){
      lambda_theta *= area_unit * survey_len * cue_rates;
    } else {
      lambda_theta *= area_unit;
    }

    //the Poisson term for the number of detected groups; the conditional likelihood
    //conditions on this number instead
    if(is_conditional == 0){
      nll += lambda_theta;
    }

    //unique binary capture histories in this session, and for each the probability of the
    //history at every mask point (it depends only on the history, as there are no
    //call-level covariates)
    std::map<std::vector<int>, int> uid_of;
    std::vector<int> call_uid(n_c);
    std::vector<int> uid_n_calls;
    std::vector<int> uid_n_det;
    std::vector<std::vector<int> > uid_hist;
    for(int c = 0; c < n_c; c++){
      std::vector<int> hist(n_t);
      int n_det = 0;
      for(int t = 0; t < n_t; t++){
        hist[t] = (asDouble(capt_bin(off_full + c * n_t + t)) == 1) ? 1 : 0;
        n_det += hist[t];
      }
      std::map<std::vector<int>, int>::iterator it = uid_of.find(hist);
      if(it == uid_of.end()){
        int u = uid_hist.size();
        uid_of[hist] = u;
        uid_hist.push_back(hist);
        uid_n_calls.push_back(1);
        uid_n_det.push_back(n_det);
        call_uid[c] = u;
      } else {
        call_uid[c] = it->second;
        uid_n_calls[it->second] += 1;
      }
    }
    int n_uid = uid_hist.size();
    matrix<Type> fw(n_uid, n_m);
    for(int u = 0; u < n_uid; u++){
      for(int m = 0; m < n_m; m++){
        fw(u, m) = Type(1.0);
        for(int t = 0; t < n_t; t++){
          Type w = Type(uid_hist[u][t]);
          if(is_ss == 0){
            fw(u, m) *= pow(p_k(m, t), w) * pow((1 - p_k(m, t)), (1 - w));
          } else {
            //the detection terms are replaced by the signal strength densities
            fw(u, m) *= pow((1 - p_k(m, t)), (1 - w));
          }
        }
      }
    }
    std::vector<int> uid_done(n_uid, 0);

    //the contribution of each group
    int c0 = 0;
    for(int g = 0; g < n_g; g++){
      int n_cg = (is_animalID == 1) ? n_calls_each_animal(off_group + g) : 1;
      Type multiplicity = Type(1.0);
      if(share_histories == 1){
        int u = call_uid[c0];
        if(uid_done[u] == 1){
          c0 += n_cg;
          continue;
        }
        uid_done[u] = 1;
        multiplicity = Type(uid_n_calls[u]);
      }

      Type l_i = Type(0.0);
      for(int m = 0; m < n_m; m++){
        int mrow = off_mask + m;
        //with local integration, only mask points near the group's detections are used
        if(is_local == 1 && index_local(off_idmask + c0 * n_m + m) != 1) continue;

        Type l_w = Type(1.0);
        Type fy_toa_log = Type(0.0);
        Type fy_bear_log = Type(0.0);
        Type fy_dist_log = Type(0.0);
        Type fy_ss_log = Type(0.0);
        for(int ic = 0; ic < n_cg; ic++){
          int c = c0 + ic;
          l_w *= fw(call_uid[c], m);

          //times of arrival
          if(is_toa == 1){
            if(!(std::isnan(asDouble(sigma_toa_vec_full(off_trap))) || std::isnan(asDouble(sigma_toa_vec_mask(mrow))))){
              Type sigma_toa_tem = linked(sigma_toa_vec_full(off_trap) + sigma_toa_vec_mask(mrow), par_link(14));
              int n_det = uid_n_det[call_uid[c]];
              fy_toa_log += (1 - n_det) * log(sigma_toa_tem) + (-0.5) * toa_ssq(off_idmask + c * n_m + m) / pow(sigma_toa_tem, 2);
            }
          }

          //bearings, distances and signal strengths exist only for detections
          for(int t = 0; t < n_t; t++){
            int row = off_full + c * n_t + t;
            int trow = off_trap + t;
            if(uid_hist[call_uid[c]][t] != 1) continue;
            int dt = off_dt + m * n_t + t;

            if(is_bearing == 1 && !std::isnan(asDouble(capt_bearing(row)))){
              Type kappa_tem = linked(kappa_vec_full(trow) + kappa_vec_mask(mrow), par_link(12));
              fy_bear_log += kappa_tem * cos(capt_bearing(row) - theta(dt)) - log(bessi0(kappa_tem));
            }

            if(is_dist == 1 && !std::isnan(asDouble(capt_dist(row)))){
              Type alpha_tem = linked(alpha_vec_full(trow) + alpha_vec_mask(mrow), par_link(13));
              fy_dist_log += ((-1) * (alpha_tem * (log(dx(dt)) - log(alpha_tem)) + log(Gamma(alpha_tem))) + (alpha_tem - 1) *
                log(capt_dist(row)) - alpha_tem * capt_dist(row) / dx(dt));
            }

            if(is_ss == 1){
              Type sigma_ss_tem = linked(sigma_ss_vec_full(trow) + sigma_ss_vec_mask(mrow), par_link(11));
              fy_ss_log += dnorm(capt_ss(row), essx(dt), sigma_ss_tem, true);
            }
          }
        }

        //the group-level factor; under the conditional likelihood, D does not enter
        if(is_animalID == 0){
          if(is_conditional == 0){
            l_i += D_m(m) * l_w * exp(fy_toa_log + fy_bear_log + fy_dist_log + fy_ss_log);
          } else {
            l_i += l_w * exp(fy_toa_log + fy_bear_log + fy_dist_log + fy_ss_log);
          }
        } else {
          l_i += (is_conditional == 0 ? D_m(m) : Type(1.0)) * dpois(Type(n_cg), Type(mu_m(m) * survey_len)) * l_w *
            exp(mu_m(m) * survey_len * (1 - p_dot(m))) *
            exp(fy_toa_log + fy_bear_log + fy_dist_log + fy_ss_log);
        }
      }

      //the call model's group-level factor includes the expected number of calls per
      //animal, cue_rates * survey_len
      if(is_animalID == 0){
        l_i *= survey_len * cue_rates;
      }
      //the conditional likelihood conditions on the number of detected groups in the
      //session, so each contribution is normalised by the session's ESA; DBL_MIN keeps
      //the log finite if the integral underflows
      if(is_conditional == 0){
        nll -= multiplicity * log((l_i * area_unit) + std::numeric_limits<double>::min());
      } else {
        nll -= multiplicity * log((l_i * area_unit / esa(s)) + std::numeric_limits<double>::min());
      }
      c0 += n_cg;
    }

    off_full += n_c * n_t;
    off_trap += n_t;
    off_mask += n_m;
    off_dt += n_m * n_t;
    off_idmask += n_c * n_m;
    if(is_animalID == 1) off_group += n_animals(s);
  }

  ADREPORT(esa);
  
  if (is_two_stage) {
    // Report the log_esa for each session
    vector<Type> log_esa(n_sessions);
    for(int s = 0; s < n_sessions; s++) log_esa(s) = log(esa(s));
    ADREPORT(log_esa);
  }
  
  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif

