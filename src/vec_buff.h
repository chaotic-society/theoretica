#ifndef UROBORO_VEC_BUFF_H
#define UROBORO_VEC_BUFF_H
#include <vector>

#include "./constants.h"

namespace uroboro {


	using vec_buff = std::vector<real>;


	// Sum the products of two sets of values
	inline real product_sum(const vec_buff& X, const vec_buff& Y) {

		if(X.size() != Y.size())
			// throw...
			return 0;

		real res = 0;
		for(int i = 0; i < X.size(); i++)
			res += X[i] * Y[i];

		return res;
	}


    // Sum the products of three sets of values
    inline real product_sum(const vec_buff& X, const vec_buff& Y, const vec_buff& Z) {

		if(X.size() != Y.size() || X.size() != Z.size())
			// throw...
			return 0;

		real res = 0;
		for(int i = 0; i < X.size(); i++)
			res += X[i] * Y[i] * Z[i];

		return res;
    }


    // Sum the quotients of two sets of values
	inline real quotient_sum(const vec_buff& X, const vec_buff& Y) {

		if(X.size() != Y.size())
			// throw...
			return 0;

		real res = 0;
		for(int i = 0; i < X.size(); i++)
			res += X[i] / Y[i];

		return res;
	}


    // Sum the squares of a set of values
    inline real sum_squares(const vec_buff& X) {

		if(!X.size())
			// throw...
			return 0;

		real res = 0;
		for(int i = 0; i < X.size(); i++)
			res += X[i] * X[i];

		return res;
    }


    // Sum together a set of values
	inline real sum(const vec_buff& X) {

		if(!X.size())
			// throw...
			return 0;

		real res = 0;
		for(int i = 0; i < X.size(); i++)
			res += X[i];

		return res;
	}

}

#endif
