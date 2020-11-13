/*******************************************************************************
* Copyright 2018 Intel Corporation
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*******************************************************************************/

#ifndef JIT_AVX512_CORE_GEMM_S8U8S32_HPP
#define JIT_AVX512_CORE_GEMM_S8U8S32_HPP

#include <cstdint>
#include "mkldnn_types.h"

namespace mkldnn {
namespace impl {
namespace cpu {

mkldnn_status_t jit_avx512_core_gemm_s8u8s32(
        const char *transA, const char *transB, const char *offsetC,
        const int *m, const int *n, const int *k,
        const float *alpha, const int8_t *a, const int *lda, const int8_t *oa,
        const uint8_t *b, const int *ldb, const int8_t *ob,
        const float *beta, int32_t *c, const int *ldc, const int32_t *oc);

}
}
}

#endif // JIT_AVX512_CORE_GEMM_S8U8S32_HPP
