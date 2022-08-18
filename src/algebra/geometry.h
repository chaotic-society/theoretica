
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
	template<typename Vector>
	inline Vector sphere_inversion(const Vector& p, const Vector& c = Vector(0), real r = 1) {
		Vector q = p - c;
		return c + q * square(r / q.length());
	}


	/// Sphere inversion of a point with respect to
	/// a sphere of radius r centered at a point c
	/// @param p The vector to transform
	/// @param c The center of the sphere
	/// @param r The radius of the sphere
	template<unsigned int N, typename T>
	inline vec<N, T> sphere_inversion(const vec<N, T>& p, const vec<N, T>& c = vec<N, T>(0), T r = T(1)) {
		vec<N> q = p - c;
		return c + q * square(r) / q.square_length();
	}


}

#endif
