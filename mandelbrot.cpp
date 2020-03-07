#include <complex>
#include <cstdio>
#include <math.h>
#include <string>
#include "./bmp_writer.hpp"

typedef std::complex<double> complex;

using std::cout, std::endl;

struct vec3 {
    float x=0, y=0, z=0;
    vec3() {}
    vec3(float vx, float vy, float vz) : x{vx}, y{vy}, z{vz} {}
    void set(vec3 v) {
        x=v.x;
        y=v.y;
        z=v.z;
    }
};

vec3 operator*(const vec3& v, const float& f) {
  return vec3(v.x*f,v.y*f,v.z*f);
}

vec3 getCol(const float& V) {
    double K{log(2)};
    double x{log(V)/K};
    double a{1/log(2)};
    double b{1/(3*sqrt(2)*log(2))};
    double c{1/(7*pow(3,0.125)*log(2))};
    vec3 rel(1+cos(a*x), 1+cos(b*x), 1+cos(c*x));
    return rel * 255 * 0.5;
}

vec3 MandelbrotCalculate(const complex& c, const size_t& maxiter) {
    // iterates z = z + c until |z| >= 2 or maxiter is reached,
    // returns the number of iterations.
    complex z{c};
    size_t n{0};
    vec3 col;
    double power{1};
    double V;
    for(; n < maxiter; ++n)
    {
        if( std::abs(z) >= 1000.0) {
            //cout << "n: " << n << endl;
            V = log(pow(std::abs(z),2))/power;
            col.set(getCol(V));
            return col;
        }
        z = z*z + c;
        power *= 2;
    }
    return vec3(255,255,255);
}

int main() {
    float scale{1};

    const int width = 1920*scale, height = 1080*scale, num_pixels = width*height;

    // Create a 24 bits/pixel BMP image in memory, modify it, save it on disk
    BMP bmp3(width, height, false);

    size_t max_iter{5};
    double scal{1E-11};

    for(size_t iter{1}; iter <= max_iter; iter++) {
      scal = scal / 5;
      const complex center(
            -1.740062382579339905220844167065825638296641720436171866879862,
            0.028175339779211048992411521144319509687539076742990608570401
          ), span(scal*3.7, -scal*(4/3.0)*2.7*height/width);
      //const complex center(-.7, 0), span(3.7, -(4/3.0)*2.7*height/width);
      const complex begin = center-span/2.0;
      const int maxiter = 10000;

      #pragma omp parallel for ordered schedule(dynamic)
      for(int pix = 0; pix<num_pixels; ++pix)
      {
          const int x{pix%width}, y{pix/width};

          complex c{begin + complex(x * span.real() / (width +1.0),
                                    y * span.imag() / (height+1.0))};

          vec3 color{MandelbrotCalculate(c, maxiter)};
          //if(n == maxiter) n = 255;

          #pragma omp ordered
          {
              //if(x == 0)
              //    std::cout << "Progress: " << static_cast<float>(100*y)/height << "%" << std::endl;
              bmp3.write_pixel(x, y, color.x, color.y, color.z, 255);
          }
      }
      std::string filename{"frames/test_"};
      for(size_t i{0}; i < floor(log10(max_iter))-floor(log10(iter)); i++) {
        filename += "0";
      }
      filename += std::to_string(iter);
      filename += ".bmp";
      bmp3.write(filename);
      if(iter % 1 == 0)
          std::cout << "Progress: " << static_cast<float>(100*iter)/max_iter << "%" << std::endl;
    }
}
