functions {
  real[] det_SEIR_meta(real t,
                       real[] y,
                       real[] theta,
                       real[] x_r,
                       int[] x_i) {

    int nm = 17;
    int nc = 152;
    real beta1 = theta[1];
    real beta2 = 0.33;
    real beta3 = theta[2];
    matrix[nm, nm] M_tt = to_matrix(x_r[1:(nm*nm)], nm, nm);
    matrix[nm, nm] M_hh = to_matrix(x_r[((nm*nm)+1):(2*nm*nm)], nm, nm);
    matrix[nc, nc] C = to_matrix(x_r[((2*nm*nm)+1):((2*nm*nm)+(nc*nc))], nc, nc);
    real Kbar = x_r[((2*nm*nm)+(nc*nc)+1)];
    real TL = x_r[((2*nm*nm)+(nc*nc)+2)];
    vector[nc*nm] N = to_vector(x_r[((2*nm*nm)+(nc*nc)+3):((2*nm*nm)+(nc*nc)+2+(nc*nm))]);
    vector[nc*nm] N_sum = to_vector(x_r[((2*nm*nm)+(nc*nc)+2+(nc*nm)+1):((2*nm*nm)+(nc*nc)+2+(nc*nm)+(nc*nm))]);
    vector[49] W = to_vector(x_r[((2*nm*nm)+(nc*nc)+2+(nc*nm)+(nc*nm)+1):((2*nm*nm)+(nc*nc)+2+(nc*nm)+(nc*nm)+49)]);
    real nu = 0.25;
    real gamma = theta[3];
    int sidxs[nc*nm] = x_i[1:nc*nm]; // indices of susceptible compartments
    int iidxs[nc*nm] = x_i[(nc*nm)+1:(2*nc*nm)]; // indices of infectious compartments
    vector[nc] i_sum;
    vector[nc*nm] infec_rate;
    vector[nc*nm] tmp;
    real dydt[nc*nm*4];

    // Real to integer conversion workaround
    int t_idx = 1;

    while (t_idx < t) {
      t_idx += 1;
    }

    for (i in 1:nc) {
      int k = (i-1)*nm;
      i_sum[i] = sum(y[iidxs[k+1:k+nm]])/N_sum[k+1];
    }

    if (t < TL) {
      for (i in 1:nc) {
        int k = (i-1)*nm;
        tmp[k+1:k+nm] = M_tt*to_vector(y[iidxs[k+1:k+nm]]) + beta2*W[t_idx]*Kbar*dot_product(row(C,i), i_sum);
      }
      infec_rate = beta1 * tmp .* to_vector(y[sidxs]) ./ N;
    } else {
      for (i in 1:nc) {
        int k = (i-1)*nm;
        tmp[k+1:k+nm] = M_hh*to_vector(y[iidxs[k+1:k+nm]]) + beta2*W[t_idx]*Kbar*dot_product(row(C,i), i_sum);
      }
      infec_rate = beta1 * beta3 * tmp .* to_vector(y[sidxs]) ./ N;
    }

    for (i in 1:nc) {
      for (j in 1:nm) {
        int k = (i-1)*nm + j;
        int l = k * 4 - 3;

        dydt[l] = -infec_rate[k];
        dydt[l+1] = infec_rate[k] - nu * y[l+1];
        dydt[l+2] = nu * y[l+1] - gamma * y[l+2];
        dydt[l+3] = gamma * y[l+2];
      }
    }

    return dydt;
  }

  real[ , ] integrate_ode_explicit_trapezoidal(real[] initial_state, real initial_time, real[] times, real[] params, real[] real_data, int[] integer_data) {
    real h;
    vector[size(initial_state)] dstate_dt_initial_time;
    vector[size(initial_state)] dstate_dt_tidx;
    vector[size(initial_state)] k;
    real state_estimate[size(times),size(initial_state)];

    h = times[1] - initial_time;
    dstate_dt_initial_time = to_vector(det_SEIR_meta(initial_time, initial_state, params, real_data, integer_data));
    k = h*dstate_dt_initial_time;
    state_estimate[1,] = to_array_1d(to_vector(initial_state) + h*(dstate_dt_initial_time + to_vector(det_SEIR_meta(times[1], to_array_1d(to_vector(initial_state)+k), params, real_data, integer_data)))/2);

    for (tidx in 1:size(times)-1) {
      h = (times[tidx+1] - times[tidx]);
      dstate_dt_tidx = to_vector(det_SEIR_meta(times[tidx], state_estimate[tidx], params, real_data, integer_data));
      k = h*dstate_dt_tidx;
      state_estimate[tidx+1,] = to_array_1d(to_vector(state_estimate[tidx,]) + h*(dstate_dt_tidx + to_vector(det_SEIR_meta(times[tidx+1], to_array_1d(to_vector(state_estimate[tidx,])+k), params, real_data, integer_data)))/2);
    }

    return state_estimate;
  }
}
data {
  int<lower=1> T;
  real y0[10336];
  real t0;
  real ts[T];
  real x_r[28901];
  int sidxs[2584];
  int eidxs[2584];
  int iidxs[2584];
  int ridxs[2584];
  int<lower=0> reported_cases[T];
}
transformed data {
  int x_i[5168] = append_array(sidxs, iidxs);
}
parameters {
  real<lower=0> beta1;
  real<lower=0> beta3;
  real<lower=0> gamma;
  //real<lower=0> I0;
  real<lower=0> reciprocal_sqrt_r;
}
model {
  beta1 ~ gamma(1, 1);
  beta3 ~ gamma(20, 20);
  gamma ~ gamma(100, 400);
  //I0 ~ gamma(1.5, 0.05);
  reciprocal_sqrt_r ~ normal(0, 0.5);
  real r = 1. / square(reciprocal_sqrt_r);

  real theta[3];
  theta[1] = beta1;
  theta[2] = beta3;
  theta[3] = gamma;

  real y_hat[T,10336] = integrate_ode_explicit_trapezoidal(y0, t0, ts, theta, x_r, x_i);

  vector[T+1] removed_total; 
  removed_total[1] = sum(y0[ridxs]);

  for (t in 2:T+1) {
    removed_total[t] = sum(y_hat[t-1, ridxs]);
  }

  vector[T] removed_total_incr = removed_total[2:] - removed_total[:size(removed_total)-1];

  real phi = 0.1;
  vector[T] mu = removed_total_incr * phi;

  reported_cases ~ neg_binomial_2(mu, r);
}
