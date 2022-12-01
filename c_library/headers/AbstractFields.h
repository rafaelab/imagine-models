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
    int ndim = 0;
    std::vector<fftw_plan> c2r;
    std::vector<fftw_plan> r2c;
    // constructors
    //RandomField() : Field<G, T>() {};
    RandomField(int i) : Field<G, T>(), ndim(i) {};
    RandomField(int i, int shape_x, int shape_y, int shape_z) : Field<G, T>(), ndim(i) {};
    ~RandomField() {
      if (clean_switch) {
        for (int i=0; i< ndim; i++) {
          fftw_destroy_plan(c2r[i]);
          fftw_destroy_plan(r2c[i]);
        }
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
  int seed;
  // for destructor
  bool clean_switch = false;
  // methods
  void at_position(const double &x, const double &y, const double &z) const = 0;

  std::vector<double> on_grid(const G &grid_x, const G &grid_y, const G &grid_z) = 0;

  std::vector<double> on_grid(const std::vector<int> &n, const std::vector<double> &zeropoint, const std::vector<double> &increment) const {
      std::vector<fftw_complex> b = draw_random_numbers(n, increment, seed);
      for (int i=0; i< ndim; i++) {
          fftw_execute_dft(c2r[i], b[i], b[i]);
      }
      auto b = RandomField::evaluate_function_on_grid(n, zeropoint, increment,
        [this](double &xx, double &yy, double &zz) {
          int _nx = (int)((xx - zeropoint[0])/increment[0])
          int _ny = (int)((yy - zeropoint[1])/increment[1])
          int _nz = (int)((zz - zeropoint[2])/increment[2])
          return b[_nz + n[2]*(_ny + n[1]*_nx) ]*spatial_profile(xx, yy, zz);
        });

      for (int i=0; i< ndim; i++) {
          fftw_execute_dft(r2c[i], b[i], b[i]);
      }

      b = divergence_cleaner(b)

      for (int i=0; i< ndim; i++) {
          fftw_execute_dft(c2r[i], b[i], b[i]);
      }

      retrun b;
      fftw_free(c0);
      fftw_free(c1);
  };


  virtual spatial_profile(const double &x, const double &y, const double &z) const = 0;


  std::vector<fftw_complex> draw_random_numbers(const std::vector<int> &n, const std::vector<double> &increment, const int seed) const {

    lx = n[0]*increment[0]
    ly = n[1]*increment[1]
    lz = n[2]*increment[2]
    std::vector<*fftw_complex> out(ndim, fftw_malloc(sizeof(fftw_complex) * nx*ny*nz);

    #ifdef _OPENMP
      std::seed_seq seq{seed}
      std::vector<std::mt19937> generators;
      std::vector<std::uint32_t> seeds(omp_get_max_threads());
      seq.generate(seeds.begin(), seeds.end());
      for (int i = 0, N = omp_get_max_threads(); i < N; ++i) {
        generators.emplace_back(std::mt19937(seeds[i]);
        }
    #else
      auto gen = std::mt19937(seed)
    #endif
      // start Fourier space filling, physical k in 1/kpc dimension
      // physical dk^3
      const double  dk3 = 1. / (lx * ly * lz);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
      for (int i = 0; i < n[0]; ++i) {
    #ifdef _OPENMP
        auto gen = generators[omp_get_thread_num()];
    #else
        auto seed_id = r;
    #endif
        double kx = i / increment[0];
        if (i >= (n[0] + 1) / 2)
          kx -= 1. / increment[0];
        // it's faster to calculate indeces manually
        const int idx_lv1 = i * n[1] * n[2];
        for (int j = 0; j < n[1]; ++j) {
          double ky = j / increment[1];
          if (j >= (n[1] + 1) / 2)
            ky -= 1. / increment[1];
          const int idx_lv2 = idx_lv1 + j * n[2];
          for (int l = 0; l < n[2]/2 + 1; ++l) {
            // 0th term is fixed to zero in allocation, last loop runs only until nz/2 due to complex array and real outputs
            if (i == 0 and j == 0 and l == 0)
              continue;
            double kz = 1. / increment[2]};
            const double ks = std::sqrt(kx * kx + ky * ky + kz * kz);
            const int idx = idx_lv2 + l;
            const double sigma = std::sqrt(0.33333333 * spectrum(ks, rms, k0, k1, a0, a1) * dk3);
            std::normal_distribution<> nd{0,sigma};
            for (int m=0; m<ndim; m++) {
              out[m][idx][0] = nd(gen);
              out[m][idx][1] = nd(gen);
              }
          } // l
        }   // j
      }     // i
      //free random memory!!!!
      return out;
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
std::vector<double> divergence_cleaner(std::vector<fftw_complex> &b_complex) const {
  if (ndim !=3) {
    throw DivergenceException();
    }

  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
    #endif
    for (int i = 0; i < nx; ++i) {
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
        for (int l = 0; l < nz/2 + 1; ++l) {
          // 0th term is fixed to zero in allocation
          if (i == 0 and j == 0 and l == 0)
            continue;
          double kz = 1. / dz};
          const int idx = idx_lv2 + l;
          double k_length = 0;
          double b_length = 0;
          double b_dot_k = 0;
          for (int m=0; m<ndim; m++) {
            k_length += static_cast<double>(k[m] * k[m]);
            b_length += static_cast<double>(b[m] * b[m]);
            b_dot_k += static_cast<double>(b[m] * k[m]);
            }
          if (k_length == 0 or b_length == 0) {
            continue;
            }
          k_length = std::sqrt(k_length);

          const double inv_k_mod = 1. / k_length;
          for (int m=0; m<ndim; m++) {
          // multiply \sqrt(3) for preserving spectral power statistically
              b_complex[m][idx] = 1.73205081 * (b_complex[m][idx] - k[m] * b_dot_k * inv_k_mod)};
           }
       } // l
     }   // j
    }     // i
    return b_complex;


  }


};
