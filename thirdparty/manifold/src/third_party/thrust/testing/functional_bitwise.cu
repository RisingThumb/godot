#include <unittest/unittest.h>
#include <thrust/functional.h>
#include <thrust/transform.h>

#include <functional>
#include <algorithm>
    
const size_t NUM_SAMPLES = 10000;

// STL doesn't necessarily have these available
namespace ref
{

template<typename T>
  struct bit_and
{
  T operator()(const T &lhs, const T &rhs) const
  {
    return lhs & rhs;
  }
};

template<typename T>
  struct bit_or
{
  T operator()(const T &lhs, const T &rhs) const
  {
    return lhs | rhs;
  }
};

template<typename T>
  struct bit_xor
{
  T operator()(const T &lhs, const T &rhs) const
  {
    return lhs ^ rhs;
  }
};

}

template <class InputVector, class OutputVector, class Operator, class ReferenceOperator>
void TestBinaryFunctional(void)
{
    typedef typename InputVector::value_type  InputType;
    typedef typename OutputVector::value_type OutputType;
    
    thrust::host_vector<InputType>  std_input1 = unittest::random_samples<InputType>(NUM_SAMPLES);
    thrust::host_vector<InputType>  std_input2 = unittest::random_samples<InputType>(NUM_SAMPLES);
    thrust::host_vector<OutputType> std_output(NUM_SAMPLES);

    // Replace zeros to avoid divide by zero exceptions
    std::replace(std_input2.begin(), std_input2.end(), (InputType) 0, (InputType) 1);

    InputVector input1 = std_input1; 
    InputVector input2 = std_input2; 
    OutputVector output(NUM_SAMPLES);

    thrust::transform(    input1.begin(),     input1.end(),      input2.begin(),     output.begin(),          Operator());
    thrust::transform(std_input1.begin(), std_input1.end(),  std_input2.begin(), std_output.begin(), ReferenceOperator());

    // Note: FP division is not bit-equal, even when nvcc is invoked with --prec-div
    ASSERT_ALMOST_EQUAL(output, std_output);
}



// XXX add bool to list
// Instantiate a macro for all integer-like data types
#define INSTANTIATE_INTEGER_TYPES(Macro, vector_type, operator_name)   \
Macro(vector_type, operator_name, unittest::int8_t  )                  \
Macro(vector_type, operator_name, unittest::uint8_t )                  \
Macro(vector_type, operator_name, unittest::int16_t )                  \
Macro(vector_type, operator_name, unittest::uint16_t)                  \
Macro(vector_type, operator_name, unittest::int32_t )                  \
Macro(vector_type, operator_name, unittest::uint32_t)                  \
Macro(vector_type, operator_name, unittest::int64_t )                  \
Macro(vector_type, operator_name, unittest::uint64_t)

// Instantiate a macro for all integer and floating point data types
#define INSTANTIATE_ALL_TYPES(Macro, vector_type, operator_name)       \
INSTANTIATE_INTEGER_TYPES(Macro, vector_type, operator_name)           \
Macro(vector_type, operator_name, float)

// bitwise_op(T,T) -> T
#define INSTANTIATE_BINARY_BITWISE_FUNCTIONAL_TEST(vector_type, operator_name, data_type) \
    TestBinaryFunctional< thrust::vector_type<data_type>,                                \
                          thrust::vector_type<data_type>,                                \
                          thrust::operator_name<data_type>,                              \
                          ref::operator_name<data_type> >();

// op(T,T) -> T (for bitwise op and integer T only)
#define DECLARE_BINARY_BITWISE_FUNCTIONAL_UNITTEST(operator_name, OperatorName)                         \
void Test##OperatorName##FunctionalHost(void)                                                           \
{                                                                                                       \
    INSTANTIATE_INTEGER_TYPES( INSTANTIATE_BINARY_BITWISE_FUNCTIONAL_TEST, host_vector,   operator_name);   \
}                                                                                                       \
DECLARE_UNITTEST(Test##OperatorName##FunctionalHost);                                                   \
void Test##OperatorName##FunctionalDevice(void)                                                         \
{                                                                                                       \
    INSTANTIATE_INTEGER_TYPES( INSTANTIATE_BINARY_BITWISE_FUNCTIONAL_TEST, device_vector, operator_name);   \
}                                                                                                       \
DECLARE_UNITTEST(Test##OperatorName##FunctionalDevice);


// Create the unit tests
DECLARE_BINARY_BITWISE_FUNCTIONAL_UNITTEST(bit_and,       BitAnd      );
DECLARE_BINARY_BITWISE_FUNCTIONAL_UNITTEST(bit_or,        BitOr       );
DECLARE_BINARY_BITWISE_FUNCTIONAL_UNITTEST(bit_xor,       BitXor      );

