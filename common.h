#ifndef _UROBORO_COMMON_H
#define _UROBORO_COMMON_H

#ifdef MML_REAL_DOUBLE
using real = double;
#else
using real = float;
#endif

namespace uroboro {


	inline real sqrt(real x) {

		#ifdef MSVC_ASM
		__asm {
			fld x
			fsqrt
			fst x
		}
		#else
		asm("fsqrt" : "+t"(x));
		#endif

		return x;
	}


	inline real exp(real x) {

		real sqrx = x * x;
		real sqrx2 = sqrx * sqrx;

		return 1.f + x
		+ (sqrx / 2.f)
		+ (sqrx * x / 6.f)
		+ (sqrx2 / 24.f)
		+ (sqrx2 * x / 120.f)
		+ (sqrx2 * sqrx / 720.f)
		+ (sqrx2 * sqrx * x / 5040.f)
		+ (sqrx2 * sqrx2 / 40320.f);
	}


	inline real ln(real x) {

		return 0;
	}


	inline real pow(real x, int n) {

		real res;
		if(n > 0) {

			res = x;
			for(int i = 1; i < n; ++i)
				res *= x;

		} else if(n < 0) {

			res = 1 / x;
			for(int i = 1; i < -n; ++i)
				res /= x;
		} else
			return x;

		return res;
	}


	inline real pow(real x, real n) {

		return ln(n * exp(x));
	}


	inline real sin(real x) {

		#ifdef MSVC_ASM
		__asm {
			fld x
			fsin
			fst x
		}
		#else
		asm("fsin" : "+t"(x));
		#endif

		return x;
	}


	inline real cos(real x) {

		#ifdef MSVC_ASM
		__asm {
			fld x
			fcos
			fst x
		}
		#else
		asm("fcos" : "+t"(x));
		#endif

		return x;
	}


	inline real tan(real x) {

		real cos = x;
		#ifdef MSVC_ASM
		__asm {
			fld x
			fsin
			fst x
			fld cos
			fcos
			fst cos
		}
		#else
		asm("fsin" : "+t"(x));
		asm("fcos" : "+t"(cos));
		#endif

		return x / cos;
	}


	inline real cot(real x) {

		real cos = x;
		#ifdef MSVC_ASM
		__asm {
			fld x
			fsin
			fst x
			fld cos
			fcos
			fst cos
		}
		#else
		asm("fsin" : "+t"(x));
		asm("fcos" : "+t"(cos));
		#endif

		return cos / x;
	}


	inline real asin(real x) {

		real sqrx = x * x;
		real sqrx2 = sqrx * sqrx;

		return x
		+ (sqrx * x / 6.f)
		+ (sqrx2 * x * 3.f / 40.f)
		+ (sqrx2 * sqrx * x * 5.f / 112.f)
		+ (sqrx2 * sqrx2 * x * 35.f / 1152.f);
	}


	inline real acos(real x) {

		return 0;
	}


	inline real atan(real x) {

		real sqrx = x * x;
		real sqrx2 = sqrx * sqrx;

		return x
		- (sqrx * x / 3.f)
		+ (sqrx2 * x / 5.f)
		- (sqrx2 * sqrx * x / 7.f)
		+ (sqrx2 * sqrx2 * sqrx * x / 11.f);
	}


	inline real abs(real x) {

		#ifdef MSVC_ASM
		__asm {
			fld x
			fabs
			fst x
		}
		#else
		asm("fabs" : "+t"(x));
		#endif


		return x;
	}


}

#endif
