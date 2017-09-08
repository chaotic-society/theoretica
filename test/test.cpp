#include "lightcpptest.h"
#include "uroboro.h"
#include "utility.h"
#include <sstream>

#define TEST_TOL(a, b) if(!good_enough(a, b)) {	\
	std::ostringstream stream; \
	stream << "Expected value: " << b << "\n\t\tFunction returned: " << a << "\n"; \
	TEST_ERROR(stream.str()); }

#define TEST_TOL_D(desc, a, b) if(!good_enough(a, b)) {	\
	std::ostringstream stream; \
	stream << desc << "Expected value: " << b << "\n\t\tFunction returned: " << a << "\n"; \
	TEST_ERROR(stream.str()); }

using namespace uroboro;

constexpr real TOLERANCE = 0.01;

bool good_enough(real a, real b) {
	return abs(b - a) < TOLERANCE;
}

void vec3_misc() {}


int main(int argc, char const *argv[]) {

	TEST_STARTUP();

	TEST_BEGIN_MODULE(uroboro);

	set_output_prec(8);

	TEST_BEGIN_VOID(sqrt);

		TEST_TOL(sqrt(4), 2);
		TEST_TOL(sqrt(2), SQRT2);

	TEST_END();


	TEST_BEGIN_OVERLOAD_VOID(log, real, real);

		TEST_TOL(log(E), 1);
		TEST_TOL(log(E * E), 2);

	TEST_END();


	TEST_BEGIN_OVERLOAD_VOID(log2, real, real);

		TEST_TOL(log2(2), 1);
		TEST_TOL(log2(4), 2);
		TEST_TOL(log2(8), 3);

	TEST_END();


	TEST_BEGIN_OVERLOAD_VOID(log10, real, real);

		// log10 is safe only for x in range 0-20
		TEST_TOL(log10(10), 1.f);
		TEST_TOL(log10(100), 2.f);
		TEST_TOL(log10(1000), 3.f);

	TEST_END();


	TEST_BEGIN_VOID(exp);

		TEST_TOL(exp(2), E * E);
		TEST_TOL(exp(1), E);

	TEST_END();


	TEST_BEGIN_VOID(powf);

		// powf is safe only for x in range 0-1
		TEST_TOL(powf(1.f, 0.5f), 1.f);

	TEST_END();


	TEST_BEGIN_VOID(sin);

		// sin, cos and tan are precise
		TEST_TOL(sin(0.5f), 0.4794255386);
		TEST_TOL(sin(3), 0.14112000806);

	TEST_END();


	TEST_BEGIN_VOID(cos);

		TEST_TOL(cos(0.5f), 0.87758256189);
		TEST_TOL(cos(3), -0.9899924966);

	TEST_END();


	TEST_BEGIN_VOID(tan);

		TEST_TOL(tan(0.5f), 0.54630248984);
		TEST_TOL(tan(3), -0.14254654307);

	TEST_END();


	TEST_BEGIN_VOID(asin);

		// After 0.9 gets less precise
		TEST_TOL(asin(0.5f), 0.5235987756);
		TEST_TOL(asin(0.9f), 1.119769515);

	TEST_END();


	TEST_BEGIN_VOID(acos);

		// After 0.9 gets less precise
		TEST_TOL(acos(0.5f), 1.0471975512);
		TEST_TOL(acos(0.9f), 0.4510268118);

	TEST_END();


	TEST_BEGIN_VOID(atan);

		// atan is really imprecise
		TEST_TOL(atan(0.5f), 0.54630248984);
		TEST_TOL(atan(0.9f), 0.78037308007);

	TEST_END();


	TEST_BEGIN_VOID(degrees);

		TEST_TOL(degrees(50), 2864.789f);

	TEST_END();


	TEST_BEGIN_VOID(degrees);

		TEST_TOL(radians(50), 0.8726646f);

	TEST_END();


	// Trick test framework with fake function
	TEST_BEGIN_VOID(vec3_misc);

		vec3 v3 = vec3();
		vec4 v4 = vec4();
		mat4 m = mat4();

		TEST_TOL(v3.x, 0);

		v3 = m * v3;

		TEST_TOL(v3.x, 0);

	TEST_END();


	TEST_END_MODULE();

	TEST_EXIT();
}
