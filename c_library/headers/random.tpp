template<typename POSTYPE, typename GRIDTYPE>
void RandomField<POSTYPE, GRIDTYPE>::draw_random_numbers_complex(fftw_complex* vec,  const std::array<int, 3> &shp, const std::array<double, 3> &inc, const int seed)  {


  bool debug_random = false;
  auto gen = std::mt19937(seed);
  int no_of_nyqists = 0;
  int no_of_lc = 0;
  int no_of_pc = 0;
  int no_of_free = 0;
  int no_of_rand = 0;

  double lx = shp[0]*inc[0];
  double ly = shp[1]*inc[1];
  double lz = shp[2]*inc[2];

  int sz_half = shp[2]/2 + 1;
  
  float nyqind_x = shp[0]/2.;
  float nyqind_y = shp[1]/2.;
  float nyqind_z = shp[2]/2.;

  for (int i = 0; i < shp[0]; ++i) {
    double kx = (double)i / lx;
    if (i >= (shp[0] + 1) / 2)
      kx -= 1./ inc[0];
    const int idx_lv1 = i * shp[1] * sz_half;
    for (int j = 0; j < shp[1]; ++j) {
      double ky = (double)j / ly;
      if (j >= (shp[1] + 1) / 2)
        ky -= 1./ inc[1];
      const int idx_lv2 = idx_lv1 + j * sz_half;
      for (int l = 0; l < sz_half; ++l) {
        const int idx = idx_lv2 + l;
        if (debug_random) {
          std::cout << "At Index (i, j, k): (" << i << j << l << ")" << std::endl;
          std::cout << "flat array index " << idx <<  std::endl;
        }
        // at first we deal with the monopoles and nyqist terms in 3d, in order to ensure a real field
        if (l == 0 and j == 0 and i == 0) {
          // Full Monopole is set to zero, dealt with seperately in the outer scope
          vec[0][0] = 1.;
          vec[0][1] = 0.;
          if (debug_random) {
          std::cout << "Type: Monopole" << std::endl;
          std::cout << "Array val (real, imag): " << vec[idx][0] << ", " << vec[idx][1] << std::endl;
          std::cout <<  "\n";
          }
          continue;
        }

        double kz = (double)l / lz;
        const double ks = std::sqrt(kx * kx + ky * ky + kz * kz);
        // TODO: make parameters below free
        double sigma = calculate_fourier_sigma(ks);
        std::normal_distribution<double> nd{0, sigma};

        bool is_nyquist = (l == 0 and j == 0 and i == nyqind_x) or
                          (l == 0 and j == nyqind_y and i == 0) or
                          (l == nyqind_z and j == 0 and i == 0) or
                          (l == 0 and j == nyqind_y and i == nyqind_x) or
                          (l == nyqind_z and j == 0 and i == nyqind_x) or
                          (l == nyqind_z and j == nyqind_y and i == 0) or
                          (l == nyqind_z and j == nyqind_y and i == nyqind_x);

        bool is_line_conjugate = !is_nyquist and (
                                 (l == 0 and j == 0 and i > nyqind_x) or
                                 (l == 0 and j == nyqind_y and i > nyqind_x) or
                                 (l == 0 and i == 0 and j > nyqind_y) or
                                 (l == 0 and i == nyqind_x and j > nyqind_y) or
                                 (l == nyqind_z and j == 0 and i > nyqind_x) or
                                 (l == nyqind_z and j == nyqind_y and i > nyqind_x) or
                                 (l == nyqind_z and i == 0 and j > nyqind_y) or
                                 (l == nyqind_z and i == nyqind_x and j > nyqind_y)
                                 );

        bool is_plane_conjugate = !(is_nyquist or is_line_conjugate) and (
                                  (l == 0 and i > nyqind_x) or
                                  (l == nyqind_z and i > nyqind_x)
                                  );


        if (is_nyquist) {
          // enforcing real nyqist term
          vec[idx][0] = nd(gen);
          vec[idx][1] = 0.;

          if (debug_random) {
            no_of_nyqists += 1;
            no_of_rand += 1;
            std::cout << "Type: Nyquist" << std::endl;
            std::cout << "Array val (real/imag): " << vec[idx][0] << ", " << vec[idx][1] << std::endl;
            std::cout <<  "\n";
          }
        }
        else if (is_line_conjugate) {
          // enforcing hermitian symmetry in the edge lines
          int line_idx = shp[0]  * shp[1] * sz_half + 1; // undefined default to catch errors
          if (j == 0 or j == nyqind_y) {
            line_idx = (shp[0] - i) * shp[1] * sz_half +  j * sz_half + l;
            if (debug_random) {
              std::cout << "Type: Line Conjugate" << std::endl;
              std::cout << "Conjugate to (shape_i - i,j,k): (" << (shp[0] - i) << j  << l << ")" << std::endl;
              }
            }
          if (i == 0 or i == nyqind_x) {
            line_idx = i * shp[1] * sz_half + (shp[1] - j) * sz_half  + l;
            if (debug_random) {
              std::cout << "Type: Line Conjugate" << std::endl;
              std::cout << "Conjugate to (i, shape_j - j,k): (" << i << (shp[1] - j)   << l << ")" << std::endl;
              }
            }
          vec[idx][0] = vec[line_idx][0];
          vec[idx][1] = -vec[line_idx][1];

          if (debug_random) {
            no_of_lc += 1;
            std::cout << "Conjugate flat array index: " << line_idx  << std::endl;
            std::cout << "array val (real/imag)" << vec[idx][0] << ", " << vec[idx][1] << std::endl;
            std::cout << "conj array val (real/imag)" << vec[line_idx][0] << ", " << vec[line_idx][1] << std::endl;
            std::cout << "\n";
          }
        }
        else if (is_plane_conjugate) {
          // enforcing hermitian symmetry in the edge x-y planes
          int plane_idx =  (shp[0] - i) * shp[1] * sz_half + (shp[1] - j) * sz_half + l;
          vec[idx][1] = -vec[plane_idx][1];
          vec[idx][0] = vec[plane_idx][0];
          
          if (debug_random) {
            no_of_pc += 1;
            std::cout << "Type: Plane Conjugate" << std::endl;
            std::cout <<  "Conjugate to (i, j, k): (" << (shp[0] - i) << (shp[1] - j)  << l << ")" << std::endl;
            std::cout <<  "Conjugate flat array index: " << plane_idx << std::endl;
            std::cout <<  "array val (real/imag)" << vec[idx][0] << ", " << vec[idx][1] << std::endl;
            std::cout <<  "conj array val (real/imag)" << vec[plane_idx][0] << ", " << vec[plane_idx][1] << std::endl;
            std::cout <<  "\n";
          }
        }
        else {

          vec[idx][1] = nd(gen);
          vec[idx][0] = nd(gen);
          
          if (debug_random) {
            std::cout << "Type: Standard" << std::endl;
            std::cout <<  "array val (real/imag)" << vec[idx][0] << ", " << vec[idx][1] << std::endl;
            std::cout <<  "\n";
            no_of_rand += 2;
            no_of_free += 1;
            }
          }


      }
    }
  }
  if (debug_random) {
    std::cout <<  "\n" << std::endl;

    std::cout <<  "\nnumber of nyqists: " << no_of_nyqists << "\nnumber of lc: " << no_of_lc <<  "\nnumber of pc: " << no_of_pc << "\nnumber of free: " << no_of_free<<  std::endl;
    std::cout <<  "\ndegrees of freedom: " << (shp[1] * nyqind_z*2 * shp[0]) << "\nnumber of random numbers drawn: " << no_of_rand << std::endl;
    std::cout <<  "\n" << std::endl;
  }

}

/*
  double lx = shp[0]*inc[0];
  double ly = shp[1]*inc[1];
  double lz = shp[2]*inc[2];

  #ifdef _OPENMP
    std::seed_seq seq{seed}
    std::vector<std::mt19937> generators;
    std::vector<std::uint32_t> seeds(omp_get_max_threads());
    seq.generate(seeds.begin(), seeds.end());
    for (int i = 0, shp = omp_get_max_threads(); i < shp; ++i) {
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
    for (int i = 0; i < shp[0]; ++i) {
  #ifdef _OPENMP
      auto gen = generators[omp_get_thread_num()];
  #endif
      double kx = i / lx;
      if (i >= (shp[0] + 1) / 2)
        kx -= 1. / inc[0];
      // it's faster to calculate indeces manually
      const int idx_lv1 = i * shp[1] * shp[2];
      for (int j = 0; j < shp[1]; ++j) {
        double ky = j / ly;
        if (j >= (shp[1] + 1) / 2)
          ky -= 1. /inc[1];
        const int idx_lv2 = idx_lv1 + j * shp[2];
        for (int l = 0; l < shp[2]/2 + 1; ++l) {
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

*/