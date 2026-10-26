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

vector<int> incheck(const vector<int> &a, const vector<int> &b){
  int lengtha = a.size();
  int lengthb = b.size();
  vector<int> ans(lengtha);
  for(int i = 0; i < lengtha; i++){
    ans(i) = 0;
    for(int j = 0; j < lengthb; j++){
      if(a(i) == b(j)){
        ans(i) = 1;
        break;
      }
    }
  }
  return ans;
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
  Type *pointer_nll = &nll;
  Type Inf = std::numeric_limits<double>::infinity();
  Type *pointer_Inf = &Inf;
  Type pi = 3.14159265358979323846;
  DATA_INTEGER(n_sessions);
  DATA_IVECTOR(n_animals);
  DATA_IVECTOR(n_IDs);
  DATA_IVECTOR(n_IDs_for_datafull); 
  DATA_IVECTOR(n_traps);
  DATA_IVECTOR(n_masks);
  DATA_IVECTOR(n_calls_each_animal);
  
  DATA_INTEGER(nrow_data_full);
  DATA_INTEGER(nrow_data_mask);
  DATA_INTEGER(nrow_dx);
  DATA_INTEGER(nrow_id_mask);
  
  DATA_VECTOR(A);
  DATA_VECTOR(survey_length);
  DATA_SCALAR(sound_speed);
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
  DATA_INTEGER(is_ss_origin);
  DATA_INTEGER(is_bearing);
  DATA_INTEGER(is_toa);
  DATA_INTEGER(is_dist);
  DATA_INTEGER(is_local);
  DATA_INTEGER(is_freqs);
  
  
  DATA_VECTOR(capt_bin);
  DATA_VECTOR(capt_bearing);
  DATA_VECTOR(capt_dist);
  DATA_VECTOR(capt_ss);
  DATA_VECTOR(capt_toa);
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
  //DATA_MATRIX(sigma_b0_ss_DX);
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
  //DATA_MATRIX(sigma_b0_ss_DX_mask);
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
  
  //declare all parameters and generate its corresponding
  //data vectors based on data_full and data_mask
  //and only report this parameter if it is used in the model
  //to avoid NaN in std, outside of this .cpp, the unused
  //parameters should be fixed in MakeADFun() by "map=list(...)"
  
  
  //in addition, very important point, although all parameters
  //are used here in their original names, such as g0, however,
  //all of them are transformed by their corresponding link function
  //in the likelihood, they need to be backtransformed
  
  int i;
  //g0
  PARAMETER_VECTOR(g0);
  
  for(i = 0; i < (par_n_col(0, 0) + par_n_col(0, 1)); i++){
    if(g0(i) < g0_bound(0, i) || g0(i) > g0_bound(1, i)) *pointer_nll += *pointer_Inf; //infinity is Inf or INFINITY?
  }
  
  vector<Type> g0_full = g0.head(par_n_col(0, 0));
  vector<Type> g0_mask = g0.tail(par_n_col(0, 1));
  vector<Type> g0_vec_full = g0_DX * g0_full;
  vector<Type> g0_vec_mask(nrow_data_mask);
  
  if(par_n_col(0, 1) == 0){
    g0_vec_mask.setZero();
  } else {
    g0_vec_mask = g0_DX_mask * g0_mask;
  }
  if(incheck_scalar(1, param_og) == 1) ADREPORT(g0);
  
  //sigma
  PARAMETER_VECTOR(sigma);
  
  for(i = 0; i < (par_n_col(1, 0) + par_n_col(1, 1)); i++){
    if(sigma(i) < sigma_bound(0, i) || sigma(i) > sigma_bound(1, i)) *pointer_nll += *pointer_Inf;
  }
  
  vector<Type> sigma_full = sigma.head(par_n_col(1, 0));
  vector<Type> sigma_mask = sigma.tail(par_n_col(1, 1));
  vector<Type> sigma_vec_full = sigma_DX * sigma_full;
  vector<Type> sigma_vec_mask(nrow_data_mask);
  
  if(par_n_col(1, 1) == 0){
    sigma_vec_mask.setZero();
  } else {
    sigma_vec_mask = sigma_DX_mask * sigma_mask;
  }
  
  if(incheck_scalar(2, param_og) == 1) ADREPORT(sigma);
  
  //lambda0
  PARAMETER_VECTOR(lambda0);
  
  for(i = 0; i < (par_n_col(2, 0) + par_n_col(2, 1)); i++){
    if(lambda0(i) < lambda0_bound(0, i) || lambda0(i) > lambda0_bound(1, i)) *pointer_nll += *pointer_Inf;
  }
  
  vector<Type> lambda0_full = lambda0.head(par_n_col(2, 0));
  vector<Type> lambda0_mask = lambda0.tail(par_n_col(2, 1));
  vector<Type> lambda0_vec_full = lambda0_DX * lambda0_full;
  vector<Type> lambda0_vec_mask(nrow_data_mask);
  
  if(par_n_col(2, 1) == 0){
    lambda0_vec_mask.setZero();
  } else {
    lambda0_vec_mask = lambda0_DX_mask * lambda0_mask;
  }
  if(incheck_scalar(3, param_og) == 1) ADREPORT(lambda0);
  
  //z
  PARAMETER_VECTOR(z);
  for(i = 0; i < (par_n_col(3, 0) + par_n_col(3, 1)); i++){
    if(z(i) < z_bound(0, i) || z(i) > z_bound(1, i)) *pointer_nll += *pointer_Inf;
  }
  
  vector<Type> z_full = z.head(par_n_col(3, 0));
  vector<Type> z_mask = z.tail(par_n_col(3, 1));
  vector<Type> z_vec_full = z_DX * z_full;
  vector<Type> z_vec_mask(nrow_data_mask);
  
  if(par_n_col(3, 1) == 0){
    z_vec_mask.setZero();
  } else {
    z_vec_mask = z_DX_mask * z_mask;
  }

  if(incheck_scalar(4, param_og) == 1) ADREPORT(z);
  
  //shape_1
  PARAMETER_VECTOR(shape_1);
  
  for(i = 0; i < (par_n_col(4, 0) + par_n_col(4, 1)); i++){
    if(shape_1(i) < shape_1_bound(0, i) || shape_1(i) > shape_1_bound(1, i)) *pointer_nll += *pointer_Inf;
  }
  
  vector<Type> shape_1_full = shape_1.head(par_n_col(4, 0));
  vector<Type> shape_1_mask = shape_1.tail(par_n_col(4, 1));
  vector<Type> shape_1_vec_full = shape_1_DX * shape_1_full;
  vector<Type> shape_1_vec_mask(nrow_data_mask);
  
  if(par_n_col(4, 1) == 0){
    shape_1_vec_mask.setZero();
  } else {
    shape_1_vec_mask = shape_1_DX_mask * shape_1_mask;
  }
  if(incheck_scalar(5, param_og) == 1) ADREPORT(shape_1);
  
  //shape_2
  PARAMETER_VECTOR(shape_2);
  
  for(i = 0; i < (par_n_col(5, 0) + par_n_col(5, 1)); i++){
    if(shape_2(i) < shape_2_bound(0, i) || shape_2(i) > shape_2_bound(1, i)) *pointer_nll += *pointer_Inf;
  }
  
  vector<Type> shape_2_full = shape_2.head(par_n_col(5, 0));
  vector<Type> shape_2_mask = shape_2.tail(par_n_col(5, 1));
  vector<Type> shape_2_vec_full = shape_2_DX * shape_2_full;
  vector<Type> shape_2_vec_mask(nrow_data_mask);
  
  if(par_n_col(5, 1) == 0){
    shape_2_vec_mask.setZero();
  } else {
    shape_2_vec_mask = shape_2_DX_mask * shape_2_mask;
  }
  if(incheck_scalar(6, param_og) == 1) ADREPORT(shape_2);
  
  //shape
  PARAMETER_VECTOR(shape);
  
  for(i = 0; i < (par_n_col(6, 0) + par_n_col(6, 1)); i++){
    if(shape(i) < shape_bound(0, i) || shape(i) > shape_bound(1, i)) *pointer_nll += *pointer_Inf;
  }
  
  vector<Type> shape_full = shape.head(par_n_col(6, 0));
  vector<Type> shape_mask = shape.tail(par_n_col(6, 1));
  vector<Type> shape_vec_full = shape_DX * shape_full;
  vector<Type> shape_vec_mask(nrow_data_mask);
  
  if(par_n_col(6, 1) == 0){
    shape_vec_mask.setZero();
  } else {
    shape_vec_mask = shape_DX_mask * shape_mask;
  }
  if(incheck_scalar(7, param_og) == 1) ADREPORT(shape);
  
  //scale
  PARAMETER_VECTOR(scale);
  
  for(i = 0; i < (par_n_col(7, 0) + par_n_col(7, 1)); i++){
    if(scale(i) < scale_bound(0, i) || scale(i) > scale_bound(1, i)) *pointer_nll += *pointer_Inf;
  }
  
  vector<Type> scale_full = scale.head(par_n_col(7, 0));
  vector<Type> scale_mask = scale.tail(par_n_col(7, 1));
  vector<Type> scale_vec_full = scale_DX * scale_full;
  vector<Type> scale_vec_mask(nrow_data_mask);
  
  if(par_n_col(7, 1) == 0){
    scale_vec_mask.setZero();
  } else {
    scale_vec_mask = scale_DX_mask * scale_mask;
  }
  if(incheck_scalar(8, param_og) == 1) ADREPORT(scale);
  
  //b0_ss
  PARAMETER_VECTOR(b0_ss);
  
  for(i = 0; i < (par_n_col(8, 0) + par_n_col(8, 1)); i++){
    if(b0_ss(i) < b0_ss_bound(0, i) || b0_ss(i) > b0_ss_bound(1, i)) *pointer_nll += *pointer_Inf;
  }
  
  vector<Type> b0_ss_full = b0_ss.head(par_n_col(8, 0));
  vector<Type> b0_ss_mask = b0_ss.tail(par_n_col(8, 1));
  vector<Type> b0_ss_vec_full = b0_ss_DX * b0_ss_full;
  vector<Type> b0_ss_vec_mask(nrow_data_mask);
  
  if(par_n_col(8, 1) == 0){
    b0_ss_vec_mask.setZero();
  } else {
    b0_ss_vec_mask = b0_ss_DX_mask * b0_ss_mask;
  }
  if(incheck_scalar(9, param_og) == 1) ADREPORT(b0_ss);
  
  //b1_ss
  PARAMETER_VECTOR(b1_ss);
  
  for(i = 0; i < (par_n_col(9, 0) + par_n_col(9, 1)); i++){
    if(b1_ss(i) < b1_ss_bound(0, i) || b1_ss(i) > b1_ss_bound(1, i)) *pointer_nll += *pointer_Inf;
  }
  
  vector<Type> b1_ss_full = b1_ss.head(par_n_col(9, 0));
  vector<Type> b1_ss_mask = b1_ss.tail(par_n_col(9, 1));
  vector<Type> b1_ss_vec_full = b1_ss_DX * b1_ss_full;
  vector<Type> b1_ss_vec_mask(nrow_data_mask);
  
  if(par_n_col(9, 1) == 0){
    b1_ss_vec_mask.setZero();
  } else {
    b1_ss_vec_mask = b1_ss_DX_mask * b1_ss_mask;
  }
  if(incheck_scalar(10, param_og) == 1) ADREPORT(b1_ss);
  
  //b2_ss
  PARAMETER_VECTOR(b2_ss);
  
  for(i = 0; i < (par_n_col(10, 0) + par_n_col(10, 1)); i++){
    if(b2_ss(i) < b2_ss_bound(0, i) || b2_ss(i) > b2_ss_bound(1, i)) *pointer_nll += *pointer_Inf;
  }
  
  vector<Type> b2_ss_full = b2_ss.head(par_n_col(10, 0));
  vector<Type> b2_ss_mask = b2_ss.tail(par_n_col(10, 1));
  vector<Type> b2_ss_vec_full = b2_ss_DX * b2_ss_full;
  vector<Type> b2_ss_vec_mask(nrow_data_mask);
  
  if(par_n_col(10, 1) == 0){
    b2_ss_vec_mask.setZero();
  } else {
    b2_ss_vec_mask = b2_ss_DX_mask * b2_ss_mask;
  }
  if(incheck_scalar(11, param_og) == 1) ADREPORT(b2_ss);
  
  
  //sigma_ss
  PARAMETER_VECTOR(sigma_ss);
  
  for(i = 0; i < (par_n_col(11, 0) + par_n_col(11, 1)); i++){
    if(sigma_ss(i) < sigma_ss_bound(0, i) || sigma_ss(i) > sigma_ss_bound(1, i)) *pointer_nll += *pointer_Inf;
  }
  
  vector<Type> sigma_ss_full = sigma_ss.head(par_n_col(11, 0));
  vector<Type> sigma_ss_mask = sigma_ss.tail(par_n_col(11, 1));
  vector<Type> sigma_ss_vec_full = sigma_ss_DX * sigma_ss_full;
  vector<Type> sigma_ss_vec_mask(nrow_data_mask);
  
  if(par_n_col(11, 1) == 0){
    sigma_ss_vec_mask.setZero();
  } else {
    sigma_ss_vec_mask = sigma_ss_DX_mask * sigma_ss_mask;
  }
  if(incheck_scalar(12, param_og) == 1) ADREPORT(sigma_ss);
  
  //kappa
  PARAMETER_VECTOR(kappa);
  
  for(i = 0; i < (par_n_col(12, 0) + par_n_col(12, 1)); i++){
    if(kappa(i) < kappa_bound(0, i) || kappa(i) > kappa_bound(1, i)) *pointer_nll += *pointer_Inf;
  }
  
  
  vector<Type> kappa_full = kappa.head(par_n_col(12, 0));
  vector<Type> kappa_mask = kappa.tail(par_n_col(12, 1));
  vector<Type> kappa_vec_full = kappa_DX * kappa_full;
  vector<Type> kappa_vec_mask(nrow_data_mask);
  
  if(par_n_col(12, 1) == 0){
    kappa_vec_mask.setZero();
  } else {
    kappa_vec_mask = kappa_DX_mask * kappa_mask;
  }
  if(incheck_scalar(13, param_og) == 1) ADREPORT(kappa);
  
  //alpha
  PARAMETER_VECTOR(alpha);
  
  for(i = 0; i < (par_n_col(13, 0) + par_n_col(13, 1)); i++){
    if(alpha(i) < alpha_bound(0, i) || alpha(i) > alpha_bound(1, i)) *pointer_nll += *pointer_Inf;
  }
  
  
  vector<Type> alpha_full = alpha.head(par_n_col(13, 0));
  vector<Type> alpha_mask = alpha.tail(par_n_col(13, 1));
  vector<Type> alpha_vec_full = alpha_DX * alpha_full;
  vector<Type> alpha_vec_mask(nrow_data_mask);
  
  if(par_n_col(13, 1) == 0){
    alpha_vec_mask.setZero();
  } else {
    alpha_vec_mask = alpha_DX_mask * alpha_mask;
  }
  if(incheck_scalar(14, param_og) == 1) ADREPORT(alpha);
  
  //sigma_toa
  PARAMETER_VECTOR(sigma_toa);
  
  for(i = 0; i < (par_n_col(14, 0) + par_n_col(14, 1)); i++){
    if(sigma_toa(i) < sigma_toa_bound(0, i) || sigma_toa(i) > sigma_toa_bound(1, i)) *pointer_nll += *pointer_Inf;
  }
  
  vector<Type> sigma_toa_full = sigma_toa.head(par_n_col(14, 0));
  vector<Type> sigma_toa_mask = sigma_toa.tail(par_n_col(14, 1));
  vector<Type> sigma_toa_vec_full = sigma_toa_DX * sigma_toa_full;
  vector<Type> sigma_toa_vec_mask(nrow_data_mask);
  
  if(par_n_col(14, 1) == 0){
    sigma_toa_vec_mask.setZero();
  } else {
    sigma_toa_vec_mask = sigma_toa_DX_mask * sigma_toa_mask;
  }
  if(incheck_scalar(15, param_og) == 1) ADREPORT(sigma_toa);
  
  //sigma_b0_ss, this is not extentable, so just declare it as a scalar
  PARAMETER(sigma_b0_ss);
  
  if(sigma_b0_ss < sigma_b0_ss_bound(0, 0) || sigma_b0_ss > sigma_b0_ss_bound(1, 0)) *pointer_nll += *pointer_Inf;
  
  
  if(incheck_scalar(16, param_og) == 1) ADREPORT(sigma_b0_ss);
  
  //settle with "D" as it will be used regardless type of model
  PARAMETER_VECTOR(D);
  
  for(i = 0; i < (par_n_col(16, 0) + par_n_col(16, 1)); i++){
    if(D(i) < D_bound(0, i) || D(i) > D_bound(1, i)) *pointer_nll += *pointer_Inf;
  }
  
  ADREPORT(D);
  
  vector<Type> D_full = D.head(par_n_col(16, 0));
  vector<Type> D_mask = D.tail(par_n_col(16, 1));
  vector<Type> D_vec_full = D_DX * D_full;
  vector<Type> D_vec_mask(nrow_data_mask);
  
  if(par_n_col(16, 1) == 0){
    D_vec_mask.setZero();
  } else {
    D_vec_mask = D_DX_mask * D_mask;
  }
  

  //deal with 'mu', which is the cue.rate for animal_ID included model
  PARAMETER_VECTOR(mu);
  
  for(i = 0; i < (par_n_col(17, 0) + par_n_col(17, 1)); i++){
    if(mu(i) < mu_bound(0, i) || mu(i) > mu_bound(1, i)) *pointer_nll += *pointer_Inf;
  }
  
  
  vector<Type> mu_full = mu.head(par_n_col(17, 0));
  vector<Type> mu_mask = mu.tail(par_n_col(17, 1));
  vector<Type> mu_vec_full = mu_DX * mu_full;
  vector<Type> mu_vec_mask(nrow_data_mask);
  
  if(par_n_col(17, 1) == 0){
    mu_vec_mask.setZero();
  } else {
    mu_vec_mask = mu_DX_mask * mu_mask;
  }
  if(incheck_scalar(18, param_og) == 1) ADREPORT(mu);

  //finally declare the latent variable "u"
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
  //  data_full:          session - call - trap (calls ordered by group; a session with no
  //                      detections has one placeholder call, see n_IDs_for_datafull)
  //  data_mask:          session - mask
  //  dx, theta and essx: session - mask - trap
  //  data_IDmask:        session - call - mask
  //---------------------------------------------------------------------------------------

  vector<Type> essx(nrow_dx);
  essx.setZero();
  vector<Type> esa(n_sessions);

  //first rows of the current session in each data layout
  int off_full = 0;
  int off_mask = 0;
  int off_dt = 0;
  int off_idmask = 0;
  //first element of n_calls_each_animal for the current session
  int off_group = 0;

  //with no auxiliary data, calls in the call model that share a binary capture history make
  //identical contributions, so each unique history is evaluated once
  int no_aux = (is_toa == 0 && is_bearing == 0 && is_dist == 0 && is_ss_origin == 0);
  int share_histories = (is_animalID == 0 && no_aux == 1 && is_local == 0);

  for(int s = 0; s < n_sessions; s++){
    int n_m = n_masks(s);
    int n_t = n_traps(s);
    int n_c = n_IDs(s);
    int n_g = (is_animalID == 1) ? n_animals(s) : n_c;
    Type area_unit = A(s);
    Type survey_len = survey_length(s);

    //density, call rate and detection probabilities at each mask point; D, mu and
    //sigma_toa have no trap-level covariates, so they come from the session's first row
    vector<Type> D_m(n_m);
    vector<Type> mu_m(n_m);
    vector<Type> p_dot(n_m);
    matrix<Type> p_k(n_m, n_t);
    Type lambda_theta = Type(0.0);
    esa(s) = Type(0.0);
    for(int m = 0; m < n_m; m++){
      int mrow = off_mask + m;
      D_m(m) = linked(D_vec_full(off_full) + D_vec_mask(mrow), par_link(16));
      p_dot(m) = Type(1.0);
      for(int t = 0; t < n_t; t++){
        int row = off_full + t;
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
        } else if(is_ss_origin == 1){
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
        mu_m(m) = linked(mu_vec_full(off_full) + mu_vec_mask(mrow), par_link(17));
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
            if(!(std::isnan(asDouble(sigma_toa_vec_full(off_full))) || std::isnan(asDouble(sigma_toa_vec_mask(mrow))))){
              Type sigma_toa_tem = linked(sigma_toa_vec_full(off_full) + sigma_toa_vec_mask(mrow), par_link(14));
              int n_det = uid_n_det[call_uid[c]];
              fy_toa_log += (1 - n_det) * log(sigma_toa_tem) + (-0.5) * toa_ssq(off_idmask + c * n_m + m) / pow(sigma_toa_tem, 2);
            }
          }

          //bearings, distances and signal strengths exist only for detections
          for(int t = 0; t < n_t; t++){
            int row = off_full + c * n_t + t;
            if(uid_hist[call_uid[c]][t] != 1) continue;
            int dt = off_dt + m * n_t + t;

            if(is_bearing == 1 && !std::isnan(asDouble(capt_bearing(row)))){
              Type kappa_tem = linked(kappa_vec_full(row) + kappa_vec_mask(mrow), par_link(12));
              fy_bear_log += kappa_tem * cos(capt_bearing(row) - theta(dt)) - log(bessi0(kappa_tem));
            }

            //the animal model also skips recorded values of exactly zero
            if(is_dist == 1 && !std::isnan(asDouble(capt_dist(row))) &&
               (is_animalID == 0 || asDouble(capt_dist(row)) != 0)){
              Type alpha_tem = linked(alpha_vec_full(row) + alpha_vec_mask(mrow), par_link(13));
              fy_dist_log += ((-1) * (alpha_tem * (log(dx(dt)) - log(alpha_tem)) + log(Gamma(alpha_tem))) + (alpha_tem - 1) *
                log(capt_dist(row)) - alpha_tem * capt_dist(row) / dx(dt));
            }

            if(is_ss_origin == 1 && (is_animalID == 0 || asDouble(capt_ss(row)) != 0)){
              Type sigma_ss_tem = linked(sigma_ss_vec_full(row) + sigma_ss_vec_mask(mrow), par_link(11));
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

      //the conditional likelihood conditions on the number of detected groups in the
      //session, so each contribution is normalised by the session's ESA
      if(is_animalID == 0){
        if(is_conditional == 0){
          nll -= multiplicity * log(l_i * area_unit * survey_len * cue_rates);
        } else {
          nll -= multiplicity * log(l_i * area_unit * survey_len * cue_rates * (1 / esa(s)));
        }
      } else {
        if(is_conditional == 0){
          nll -= log((l_i * area_unit) + std::numeric_limits<double>::min());
        } else {
          nll -= log((l_i * area_unit / esa(s)) + std::numeric_limits<double>::min());
        }
      }
      c0 += n_cg;
    }

    off_full += n_IDs_for_datafull(s) * n_t;
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

