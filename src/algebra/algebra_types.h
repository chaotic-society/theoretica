
///
/// @file algebra_types.h Linear algebra type definitions
///

#ifndef THEORETICA_ALGEBRA_TYPES_H
#define THEORETICA_ALGEBRA_TYPES_H

#include "./vec.h"
#include "./mat.h"


namespace theoretica {


	/// A 2x2 matrix with real entries
	using mat2 = mat<2, 2, real>;

	/// A 3x3 matrix with real entries
	using mat3 = mat<3, 3, real>;

	/// A 4x4 matrix with real entries
	using mat4 = mat<4, 4, real>;

	/// A 2x2 matrix with real entries
	using cmat2 = mat<2, 2, complex<>>;

	/// A 3x3 matrix with real entries
	using cmat3 = mat<3, 3, complex<>>;

	/// A 4x4 matrix with real entries
	using cmat4 = mat<4, 4, complex<>>;


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
