
///
/// @file geometry.h Geometrical transformations
///

#ifndef THEORETICA_GEOMETRY_H
#define THEORETICA_GEOMETRY_H


namespace theoretica {


	/// Sphere inversion of a point with respect to
	/// a sphere of radius r centered at a point c
	/// @param p The vector to transform
	/// @param c The center of the sphere
	/// @param r The radius of the sphere
	template<typename Vector1, typename Vector2>
	inline Vector1 sphere_inversion(
		const Vector1& p, const Vector2& c = Vector2(0), real r = 1) {

		Vector1 q = p - c;
		return c + q * square(r / q.norm());
	}


}

#endif
