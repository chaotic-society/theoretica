#include "uroboro.h"
#include "utility.h"

using namespace umath;

int main(int argc, char const *argv[]) {

	auto vec1 = vec<4>({1, 3, 1, 2});
	auto vec2 = vec<3>({-1, 2, 1});

	auto mat1 = mat<4, 3>(2);
	auto mat2 = mat1 * 3;

	print_mat(mat1);
	print_mat(mat2);
	print_vec(vec1);
	print_vec(vec2);

	mat1.set_row(0, vec1);
	print_mat(mat1);

	mat1.set_column(1, vec2);
	print_mat(mat1);

	print_mat(mat1.transposed());

	return 0;
}
