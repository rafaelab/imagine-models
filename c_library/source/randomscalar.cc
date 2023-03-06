#include "../headers/RandomField.h"

void RandomScalarField::draw_random_numbers(fftw_complex* vec,  const std::array<int, 3> &grid_shape, const std::array<double, 3> &grid_increment, const int seed)  {

  double lx = grid_shape[0]*grid_increment[0];
  double ly = grid_shape[1]*grid_increment[1];
  double lz = grid_shape[2]*grid_increment[2];

  #ifdef _OPENMP
    std::seed_seq seq{seed}
    std::vector<std::mt19937> generators;
    std::vector<std::uint32_t> seeds(omp_get_max_threads());
    seq.generate(seeds.begin(), seeds.end());
    for (int i = 0, grid_shape = omp_get_max_threads(); i < grid_shape; ++i) {
      generators.emplace_back(std::mt19937(seeds[i]);
      }
  #else
    auto gen = std::mt19937(seed);
  #endif
  // start Fourier space filling, physical k in 1/kpc dimension
  // physical dk^3
  const double  dk3 = 1. / (lx * ly * lz);
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
    for (int i = 0; i < grid_shape[0]; ++i) {
  #ifdef _OPENMP
      auto gen = generators[omp_get_thread_num()];
  #endif
      double kx = i / lx;
      if (i >= (grid_shape[0] + 1) / 2)
        kx -= 1. / grid_increment[0];
      // it's faster to calculate indeces manually
      const int idx_lv1 = i * grid_shape[1] * grid_shape[2];
      for (int j = 0; j < grid_shape[1]; ++j) {
        double ky = j / ly;
        if (j >= (grid_shape[1] + 1) / 2)
          ky -= 1. /grid_increment[1];
        const int idx_lv2 = idx_lv1 + j * grid_shape[2];
        for (int l = 0; l < grid_shape[2]/2 + 1; ++l) {
          // 0th term is fixed to zero in allocation, last loop runs only until nz/2 due to complex array and real outputs
          if (i == 0 and j == 0 and l == 0)
            continue;
          double kz = 1. / lz;
          const double ks = std::sqrt(kx * kx + ky * ky + kz * kz);
          const int idx = idx_lv2 + l;
          const double sigma = std::sqrt(0.33333333 * spectrum(ks, rms, k0, k1, a0, a1) * dk3);
          std::normal_distribution<> nd{0, sigma};
          vec[idx][0] = nd(gen);
          vec[idx][1] = nd(gen);
          //  c_field[m][idx][1] = nd(gen);
        
        } // l
      }   // j
    }     // i
  };

