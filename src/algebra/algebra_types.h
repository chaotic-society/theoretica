
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

	/// A variable size matrix with complex entries
	using cmat = mat<complex<>>;

	/// A 2x2 matrix with real entries
	using cmat2 = mat<complex<real>, 2, 2>;

	/// A 3x3 matrix with real entries
	using cmat3 = mat<complex<real>, 3, 3>;

	/// A 4x4 matrix with real entries
	using cmat4 = mat<complex<real>, 4, 4>;


	/// A 2-dimensional vector with real elements
	using vec2 = vec<real, 2>;

	/// A 3-dimensional vector with real elements
	using vec3 = vec<real, 3>;

	/// A 4-dimensional vector with real elements
	using vec4 = vec<real, 4>;

	/// A variable size vector with complex elements
	using cvec = vec<complex<>>;

	/// A 2-dimensional vector with complex elements
	using cvec2 = vec<complex<>, 2>;

	/// A 3-dimensional vector with complex elements
	using cvec3 = vec<complex<>, 3>;

	/// A 4-dimensional vector with complex elements
	using cvec4 = vec<complex<>, 4>;

}

#endif
