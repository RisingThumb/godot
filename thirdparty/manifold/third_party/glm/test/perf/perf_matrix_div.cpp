#define GLM_FORCE_INLINE
#include <glm/ext/matrix_float4x4.hpp>
#include <glm/ext/matrix_double4x4.hpp>
#include <glm/ext/matrix_transform.hpp>
#include <glm/ext/matrix_relational.hpp>
#include <glm/ext/vector_float4.hpp>
#if GLM_CONFIG_SIMD == GLM_ENABLE
#include <glm/gtc/type_aligned.hpp>
#include <vector>
#include <chrono>
#include <cstdio>

template <typename matType>
static void test_mat_div_mat(matType const& M, std::vector<matType> const& I, std::vector<matType>& O)
{
	for (std::size_t i = 0, n = I.size(); i < n; ++i)
		O[i] = M / I[i];
}

template <typename matType>
static int launch_mat_div_mat(std::vector<matType>& O, matType const& Transform, matType const& Scale, std::size_t Samples)
{
	typedef typename matType::value_type T;

	std::vector<matType> I(Samples);
	O.resize(Samples);

	for(std::size_t i = 0; i < Samples; ++i)
		I[i] = Scale * static_cast<T>(i) + Scale;

	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	test_mat_div_mat<matType>(Transform, I, O);
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

	return static_cast<int>(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count());
}

template <typename packedMatType, typename alignedMatType>
static int comp_mat2_div_mat2(std::size_t Samples)
{
	typedef typename packedMatType::value_type T;
	
	int Error = 0;

	packedMatType const Transform(1, 2, 3, 4);
	packedMatType const Scale(0.01, 0.02, 0.03, 0.05);

	std::vector<packedMatType> SISD;
	std::printf("- SISD: %d us\n", launch_mat_div_mat<packedMatType>(SISD, Transform, Scale, Samples));

	std::vector<alignedMatType> SIMD;
	std::printf("- SIMD: %d us\n", launch_mat_div_mat<alignedMatType>(SIMD, Transform, Scale, Samples));

	for(std::size_t i = 0; i < Samples; ++i)
	{
		packedMatType const A = SISD[i];
		packedMatType const B = SIMD[i];
		Error += glm::all(glm::equal(A, B, static_cast<T>(0.001))) ? 0 : 1;
		assert(!Error);
	}
	
	return Error;
}

template <typename packedMatType, typename alignedMatType>
static int comp_mat3_div_mat3(std::size_t Samples)
{
	typedef typename packedMatType::value_type T;
	
	int Error = 0;

	packedMatType const Transform(1, 2, 3, 4, 5, 6, 7, 8, 9);
	packedMatType const Scale(0.01, 0.02, 0.03, 0.05, 0.01, 0.02, 0.03, 0.05, 0.01);

	std::vector<packedMatType> SISD;
	std::printf("- SISD: %d us\n", launch_mat_div_mat<packedMatType>(SISD, Transform, Scale, Samples));

	std::vector<alignedMatType> SIMD;
	std::printf("- SIMD: %d us\n", launch_mat_div_mat<alignedMatType>(SIMD, Transform, Scale, Samples));

	for(std::size_t i = 0; i < Samples; ++i)
	{
		packedMatType const A = SISD[i];
		packedMatType const B = SIMD[i];
		Error += glm::all(glm::equal(A, B, static_cast<T>(0.001))) ? 0 : 1;
		assert(!Error);
	}
	
	return Error;
}

template <typename packedMatType, typename alignedMatType>
static int comp_mat4_div_mat4(std::size_t Samples)
{
	typedef typename packedMatType::value_type T;
	
	int Error = 0;

	packedMatType const Transform(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);
	packedMatType const Scale(0.01, 0.02, 0.05, 0.04, 0.02, 0.08, 0.05, 0.01, 0.08, 0.03, 0.05, 0.06, 0.02, 0.03, 0.07, 0.05);

	std::vector<packedMatType> SISD;
	std::printf("- SISD: %d us\n", launch_mat_div_mat<packedMatType>(SISD, Transform, Scale, Samples));

	std::vector<alignedMatType> SIMD;
	std::printf("- SIMD: %d us\n", launch_mat_div_mat<alignedMatType>(SIMD, Transform, Scale, Samples));

	for(std::size_t i = 0; i < Samples; ++i)
	{
		packedMatType const A = SISD[i];
		packedMatType const B = SIMD[i];
		Error += glm::all(glm::equal(A, B, static_cast<T>(0.001))) ? 0 : 1;
		assert(!Error);
	}
	
	return Error;
}

int main()
{
	std::size_t const Samples = 100000;

	int Error = 0;

	std::printf("mat2 / mat2:\n");
	Error += comp_mat2_div_mat2<glm::mat2, glm::aligned_mat2>(Samples);
	
	std::printf("dmat2 / dmat2:\n");
	Error += comp_mat2_div_mat2<glm::dmat2, glm::aligned_dmat2>(Samples);

	std::printf("mat3 / mat3:\n");
	Error += comp_mat3_div_mat3<glm::mat3, glm::aligned_mat3>(Samples);
	
	std::printf("dmat3 / dmat3:\n");
	Error += comp_mat3_div_mat3<glm::dmat3, glm::aligned_dmat3>(Samples);

	std::printf("mat4 / mat4:\n");
	Error += comp_mat4_div_mat4<glm::mat4, glm::aligned_mat4>(Samples);
	
	std::printf("dmat4 / dmat4:\n");
	Error += comp_mat4_div_mat4<glm::dmat4, glm::aligned_dmat4>(Samples);

	return Error;
}

#else

int main()
{
	return 0;
}

#endif
