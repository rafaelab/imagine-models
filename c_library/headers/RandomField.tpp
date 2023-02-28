template<typename POSTYPE, typename GRIDTYPE>
void RandomField<POSTYPE, GRIDTYPE>::draw_random_numbers(std::array<fftw_complex*, 3> vec,  const std::array<int, 3> &grid_shape, const std::array<double, 3> &grid_increment, const int seed)  {

  int lx = grid_shape[0]*grid_increment[0];
  int ly = grid_shape[1]*grid_increment[1];
  int lz = grid_shape[2]*grid_increment[2];

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
          for (auto element: vec) {
            element[idx][0] = nd(gen);
            element[idx][1] = nd(gen);
          //  c_field[m][idx][1] = nd(gen);
            }
        } // l
      }   // j
    }     // i
  };

  // this function is adapted from https://github.com/hammurabi-dev/hammurabiX/blob/master/source/field/b/brnd_jf12.cc
  // original author: https://github.com/gioacchinowang
template<typename POSTYPE, typename GRIDTYPE>
double RandomField<POSTYPE, GRIDTYPE>::spectrum(const double &abs_k, const double &rms, const double &k0, const double &k1, const double &a0, const double &a1) const {
const double p0 = rms*rms;
double pi = 3.141592653589793;
const double unit = 1. / (4 * pi * abs_k * abs_k);   // units fixing, wave vector in 1/kpc units
// power laws
const double band1 = double(abs_k < k1);
const double band2 = double(abs_k > k1) * double(abs_k < k0);
const double band3 = double(abs_k > k0);
const double P = band1 * std::pow(k0 / k1, a1) * std::pow(abs_k / k1, 6.0) +
                 band2 / std::pow(abs_k / k0, a1) +
                 band3 / std::pow(abs_k / k0, a0);
return P * p0 * unit;
}
