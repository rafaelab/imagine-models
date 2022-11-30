#include <vector>
#include <cassert>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <fftw3.h>
#include <random>

// #include "exceptions.h"



template<typename G, typename T>
class Field {
protected:
    // Fields

    // Constructors
    Field() {};
    //Field(const Field<G, T>& f) {};
    // Methods

public:
  // Constructors
  virtual ~Field() = default;
  // Fields

  // Methods
  virtual T at_position(const double &x, const double &y, const double &z) const = 0;
  virtual std::vector<double> on_grid(const G &grid_x, const G &grid_y, const G &grid_z) const = 0;
  virtual std::vector<double> on_grid(const int &nx, const int &ny, const int &nz, const double &dx, const double &dy, const double &dz) const = 0;

  // This is the interface function to CRPRopa
  std::vector<double> getField(const std::vector<double> &pos_vec) const {
    T evm = at_position(pos_vec[0], pos_vec[1], pos_vec[2]);
    return evm;
  }
    // Evaluate scalar valued functions on irregular grids
  std::vector<double> evaluate_function_on_grid(const int size_x, const int size_y, const int size_z,
                                                const G &ggx, const G &ggy, const G &ggz,
                                                std::function<double(double, double, double)> eval) const {
     std::vector<double> b(size_x*size_y*size_z);
     for (int i=0; i < size_x; i++) {
         int m = i*size_x;
         for (int j=0; j < size_y; j++) {
             int n = j*size_y;
             for (int k=0; k < size_z; k++) {
                 b[m + n + k] = eval(ggx.at(i), ggy.at(j), ggz.at(k));
         }   }   }
     return b;
     }


    // Evaluate vector valued functions on irregular grids
  std::vector<double> evaluate_function_on_grid(const int size_x, const int size_y, const int size_z,
                                                const G &ggx, const G &ggy, const G &ggz,
                                                std::function<std::vector<double>(double, double, double)> eval) const {
     std::vector<double> b;
     std::vector<double>::iterator bend = b.begin();
     for (int i=0; i < size_x; i++) {
         for (int j=0; j < size_y; j++) {
             for (int k=0; k < size_z; k++) {
                 std::vector<double> v = eval(ggx.at(i), ggy.at(j), ggz.at(k));
                 b.insert(b.end(), v.begin(), v.end());
         }   }   }
     return b;
   }

     // Evaluate vector valued functions on regular grids
   std::vector<double> evaluate_function_on_grid(const int size_x, const int size_y, const int size_z,
                                                 const double &dx, const double &dy, const double &dz,
                                                 const double &zpx, const double &zpy, const double &zpz,
                                                 std::function<std::vector<double>(double, double, double)> eval) const {
      std::vector<double> b;
      for (int i=0; i < size_x; i++) {
          for (int j=0; j < size_y; j++) {
              for (int k=0; k < size_z; k++) {
                  std::vector<double> v = eval(zpx + i*dx, zpy + j*dy, zpz + k*dz);
                  b.insert(b.end(), v.begin(), v.end());
          }   }   }
      return b;
    }

    // Evaluate scalar valued functions on irregular grids
  std::vector<double> evaluate_function_on_grid(const G &ggx, const G &ggy, const G &ggz,
                                                const double &dx, const double &dy, const double &dz,
                                                const double &zpx, const double &zpy, const double &zpz,
                                                std::function<double(double, double, double)> eval) const {
     std::vector<double> b(size_x*size_y*size_z);
     for (int i=0; i < size_x; i++) {
         int m = i*size_x;
         for (int j=0; j < size_y; j++) {
             int n = j*size_y;
             for (int k=0; k < size_z; k++) {
                 b[m + n + k] = eval(zpx + i*dx, zpy + j*dy, zpz + k*dz);
         }   }   }
     return b;
     }

};


template<typename G, typename T>
class RegularField : public Field<G, T>  {
protected:
    // Fields

    // Constructors
    using Field<G, T> :: Field;
    RegularField() : Field<G, T>() {};
    // Methods

public:
  // Constructors
  virtual ~RegularField() = default;
  // Fields

  // Methods
  virtual T at_position(const double &x, const double &y, const double &z) const = 0;

  std::vector<double> on_grid(const G &grid_x, const G &grid_y, const G &grid_z) const {
      int siz_x = grid_x.size();
      int siz_y = grid_y.size();
      int siz_z = grid_z.size();

      std::vector<double> b = RegularField::evaluate_function_on_grid(grid_x, grid_y, grid_z, siz_x, siz_y, siz_z,
        [this](double &xx, double &yy, double &zz) {return at_position(xx, yy, zz);});
    return b;
  }

  std::vector<double> on_grid(const std::vector<int> &n, const std::vector<double> &zeropoint, const std::vector<double> &increment) const {
    std::vector<double> b = RegularField::evaluate_function_on_grid(n, zeropoint, increment,
      [this](double &xx, double &yy, double &zz) {return at_position(xx, yy, zz);});
  return b;
  };


};


template<typename G, typename T>
class SumRegularField : public RegularField<G, T>{
public:
  SumRegularField(const RegularField<G, T>& summand1, const RegularField<G, T>& summand2) :
  RegularField<G, T>(summand1),
    summand1(summand1), summand2(summand2)  {
      }

  const RegularField<G, T> &summand1;
  const RegularField<G, T> &summand2;

  T at_position(const double& x, const double& y, const double& z) const {

  T sum;
  if constexpr (std::is_same<T, double>::value) {
    sum = sum_scalar_models(x, y, z);
    }
  else
    {
    sum =  sum_vector_models(x, y, z);
  }
  return sum;
}

  std::vector<double> sum_vector_models(const double& x, const double& y,
    const double& z) const {
      std::vector<double> sum1 = summand1.at_position(x, y, z);
      std::vector<double> sum2 = summand2.at_position(x, y, z);
      std::vector<double> sum(3);

      std::transform(sum1.begin(), sum1.end(), sum2.begin(), sum.begin(), std::plus<double>());
      return sum;
      }

  double sum_scalar_models(const double& x, const double& y, const double& z) const {
      return summand1.at_position(x, y, z) + summand2.at_position(x, y, z);
    }

};

template<typename G, typename T>
class RandomField  : public Field<G, T>  {
protected:
    // Fields
    // constructors
    using Field<G, T> :: Field;

    RandomField() : Field<G, T>() {};
    ~RandomField() {
      if (clean_switch) {
      fftw_destroy_plan(plan_c0_bw);
      fftw_destroy_plan(plan_c1_bw);
      fftw_destroy_plan(plan_c0_fw);
      fftw_destroy_plan(plan_c1_fw);
      fftw_free(c0);
      fftw_free(c1);
      #ifdef _OPENMP
            fftw_cleanup_threads();
      #else
            fftw_cleanup();
      #endif};
    }

    // methods

public:

  // Fields
  double rms, k0, k1, a0, a1 // Pspec params
  std::unique_ptr<ham_float[]> bx, by, bz;
  // Fourier domain magnetic field
  // for/backward FFT plans
  fftw_plan plan_c0_bw, plan_c1_bw, plan_c0_fw, plan_c1_fw;
  // for destructor
  bool clean_switch = false;
  // methods
  void at_position(const double &x, const double &y, const double &z) const = 0;

  std::vector<double> on_grid(const G &grid_x, const G &grid_y, const G &grid_z) = 0;

  std::vector<double> on_grid(const int &nx, const int &ny, const int &nz, const double &dx, const double &dy, const double &dz) const {
      std::vector<fftw_complex> random_numbers = random_numbers();
      fftw_execute_dft(plan_c0_bw, c0, c0);
      fftw_execute_dft(plan_c1_bw, c1, c1);
  };


  virtual spatial_profile(const double &x, const double &y, const double &z) const = 0;


  std::vector<fftw_complex> draw_random_numbers(const double &nx, const double &ny, const double &nz, const double &dx, const double &dy, const double &dz, const int seed) const {

  fftw_complex *c0, *c1;
  #ifdef _OPENMP
    random_device r;
    std::vector<std::default_random_engine> generators;
    gsl_rng **threadvec = new gsl_rng *[omp_get_max_threads()];
    for (int i = 0, N = omp_get_max_threads(); i < N; ++i) {
    generators.emplace_back(default_random_engine(r()));
      }
  #else
    gsl_rng *r{gsl_rng_alloc(gsl_rng_taus)};
    gsl_rng_set(r, toolkit::random_seed(seed));
  #endif
    // start Fourier space filling, physical k in 1/kpc dimension
    // physical dk^3
    const double  dk3 = 1. / (lx * ly * lz);
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
    for int i = 0; i < nx; ++i) {
  #ifdef _OPENMP
      auto seed_id = threadvec[omp_get_thread_num()];
  #else
      auto seed_id = r;
  #endif
      double kx = i / dx;
      if (i >= (nx + 1) / 2)
        kx -= 1. / dx;
      // it's faster to calculate indeces manually
      const int idx_lv1 = i * ny * nz;
      for (int j = 0; j < ny; ++j) {
        double ky = j / ly;
        if (j >= (ny + 1) / 2)
          ky -= 1. / dy;
        const int idx_lv2 = idx_lv1 + j * nz;
        for int l = 0; l < nz; ++l) {
          // 0th term is fixed to zero in allocation
          if (i == 0 and j == 0 and l == 0)
            continue;
          double kz = 1. / dz};
          if (l >= (nz + 1) / 2)
            kz -= 1. / dz;
          const double ks = std::sqrt(kx * kx + ky * ky + kz * kz);
          const int idx = idx_lv2 + l;
          // turbulent power is shared in following pattern
          // P ~ (bx^2 + by^2 + bz^2)
          // c0^2 ~ c1^2 ~ (bx^2 + by^2) ~ P*2/3
          // as renormalization comes in PHASE II,
          // 1/3, P0 in spectrum, dk3 are numerically redundant
          // while useful for precision check
          const double sigma = std::sqrt(0.33333333 * spectrum(ks, rms, k0, k1, a0, a1) * dk3);
          c0[idx][0] = sigma * gsl_ran_ugaussian(seed_id);
          c0[idx][1] = sigma * gsl_ran_ugaussian(seed_id);
          c1[idx][0] = sigma * gsl_ran_ugaussian(seed_id);
          c1[idx][1] = sigma * gsl_ran_ugaussian(seed_id); //Hermitian symmetry??
        } // l
      }   // j
    }     // i

    }

  // this function is adapted from https://github.com/hammurabi-dev/hammurabiX/blob/master/source/field/b/brnd_jf12.cc
  // original author: https://github.com/gioacchinowang
  double spectrum(const double &abs_k, const double rms, const double k0, const double k1, const double a0, const double a1) const {
  const double p0 = rms*rms;
  const double unit = 1. / (4 * cgs::pi * abs_k * abs_k);   // units fixing, wave vector in 1/kpc units
  // power laws
  const double band1 = double(abs_k < k1);
  const double band2 = double(abs_k > k1) * double(k < k0);
  const double band3 = double(abs_k > k0);
  const double P = band1 * std::pow(k0 / k1, a1) * std::pow(abs_k / k1, 6.0) +
                   band2 / std::pow(abs_k / k0, a1) +
                   band3 / std::pow(abs_k / k0, a0);
  return P * p0 * unit;
}

// this function is adapted from https://github.com/hammurabi-dev/hammurabiX/blob/master/source/field/b/brnd_jf12.cc
// original author: https://github.com/gioacchinowang
std::vector<double> divergence_cleaner(const std::vector<double> &k,
                                       const std::vector<double> &b) const {

  double k_length = 0;
  double b_length = 0;
  double b_dot_k = 0;
  for (int i = 0; i < 3; ++i) {
    k_length += static_cast<double>(k[i] * k[i]);
    b_length += static_cast<double>(b[i] * b[i]);
    b_dot_k += static_cast<double>(b[i] * k[i]);
    }
  if (k_length == 0 or b_length == 0) {
    return b;
    }
  k_length = std::sqrt(k_length);

  const double inv_k_mod = 1. / k_length;
  // multiply \sqrt(3) for preserving spectral power statistically
  return std::vector<double>{
      1.73205081 * (b[0] - k[0] * b_dot_k * inv_k_mod),
      1.73205081 * (b[1] - k[1] * b_dot_k * inv_k_mod),
      1.73205081 * (b[2] - k[2] * b_dot_k * inv_k_mod)};
}


};
