
///
/// @file algebra_types.h Linear algebra type definitions
///

#ifndef THEORETICA_ALGEBRA_TYPES_H
#define THEORETICA_ALGEBRA_TYPES_H

#include "./vec.h"
#include "./mat.h"


namespace theoretica {


	/// A 2x2 matrix with real entries
	using mat2 = mat<real, 2, 2>;

	/// A 3x3 matrix with real entries
	using mat3 = mat<real, 3, 3>;

	/// A 4x4 matrix with real entries
	using mat4 = mat<real, 4, 4>;

	/// A 2x2 matrix with real entries
	using cmat2 = mat<complex<real>, 2, 2>;

	/// A 3x3 matrix with real entries
	using cmat3 = mat<complex<real>, 3, 3>;

	/// A 4x4 matrix with real entries
	using cmat4 = mat<complex<real>, 4, 4>;


	/// A 2-dimensional vector with real elements
	using vec2 = vec<2, real>;

	/// A 3-dimensional vector with real elements
	using vec3 = vec<3, real>;

	/// A 4-dimensional vector with real elements
	using vec4 = vec<4, real>;

	/// A 2-dimensional vector with complex elements
	using cvec2 = vec<2, complex<>>;

	/// A 3-dimensional vector with complex elements
	using cvec3 = vec<3, complex<>>;

	/// A 4-dimensional vector with complex elements
	using cvec4 = vec<4, complex<>>;

}

#endif
