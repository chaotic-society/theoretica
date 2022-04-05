#ifndef UROBORO_PSEUDO_RANDOM_H
#define UROBORO_PSEUDO_RANDOM_H

namespace uroboro {

  unsigned int rand_congruential(unsigned int x, unsigned int a, unsigned int c, unsigned int m) {

  	return (a * x + c) % m;
  }

  unsigned int rand_congruential(unsigned int x, const std::vector<unsigned int>& state) {

  	if(state.size() != 3) {
  		UMATH_ERROR("rand_congruential",state.size(),INVALID_ARGUMENT);
  		return 0;
  	}

  	const unsigned int a = state[0];
  	const unsigned int c = state[1];
  	const unsigned int m = state[2];

  	if(a > m || c > m) {
  		UMATH_ERROR("rand_congruential", max(a, c), INVALID_ARGUMENT);
  		return 0;
  	}

  	return rand_congruential(x, a, c, m);
  }

}


#endif