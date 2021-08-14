/////////////////////////////////////////////////////////////////////////////////
// File changes (yyyy-mm-dd)
// 2021-28-7: jeremymoyes3@gmail.com: first commit
/////////////////////////////////////////////////////////////////////////////////

#[compute]

#version 450

#VERSION_DEFINES

#extension GL_GOOGLE_include_directive : require

#define A_GPU
#define A_GLSL

#ifdef MODE_FSR_UPSCALE_NORMAL

#define A_HALF

#endif

#include "ffx_a.h"

layout(local_size_x = 64, local_size_y = 1, local_size_z = 1) in;

layout(rgba16f, set = 1, binding = 0) uniform restrict writeonly image2D fsr_image;
layout(set = 0, binding = 0) uniform sampler2D source_image;

#define FSR_UPSCALE_PASS_TYPE_EASU 0
#define FSR_UPSCALE_PASS_TYPE_RCAS 1

layout(push_constant, binding = 1, std430) uniform Params {
	float resolution_width, resolution_height;
	float upscaled_width, upscaled_height;
	float sharpness;
	int pass;
}
params;

AU4 Const0, Const1, Const2, Const3;

#ifdef MODE_FSR_UPSCALE_FALLBACK

#define FSR_EASU_F
AF4 FsrEasuRF(AF2 p) {
	AF4 res = textureGather(source_image, p, 0);
	return res;
}
AF4 FsrEasuGF(AF2 p) {
	AF4 res = textureGather(source_image, p, 1);
	return res;
}
AF4 FsrEasuBF(AF2 p) {
	AF4 res = textureGather(source_image, p, 2);
	return res;
}

#define FSR_RCAS_F
AF4 FsrRcasLoadF(ASU2 p) {
	return AF4(texelFetch(source_image, ASU2(p), 0));
}
void FsrRcasInputF(inout AF1 r, inout AF1 g, inout AF1 b) {}

#else

#define FSR_EASU_H
AH4 FsrEasuRH(AF2 p) {
	AH4 res = AH4(textureGather(source_image, p, 0));
	return res;
}
AH4 FsrEasuGH(AF2 p) {
	AH4 res = AH4(textureGather(source_image, p, 1));
	return res;
}
AH4 FsrEasuBH(AF2 p) {
	AH4 res = AH4(textureGather(source_image, p, 2));
	return res;
}

#define FSR_RCAS_H
AH4 FsrRcasLoadH(ASW2 p) {
	return AH4(texelFetch(source_image, ASU2(p), 0));
}
void FsrRcasInputH(inout AH1 r, inout AH1 g, inout AH1 b) {}

#endif

#include "ffx_fsr1.h"

void fsr_easu_pass(AU2 pos) {
#ifdef MODE_FSR_UPSCALE_NORMAL

	AH3 Gamma2Color = AH3(0, 0, 0);
	FsrEasuH(Gamma2Color, pos, Const0, Const1, Const2, Const3);
	imageStore(fsr_image, ASU2(pos), AH4(Gamma2Color, 1));

#else

	AF3 Gamma2Color = AF3(0, 0, 0);
	FsrEasuF(Gamma2Color, pos, Const0, Const1, Const2, Const3);
	imageStore(fsr_image, ASU2(pos), AF4(Gamma2Color, 1));

#endif
}

void fsr_rcas_pass(AU2 pos) {
#ifdef MODE_FSR_UPSCALE_NORMAL

	AH3 Gamma2Color = AH3(0, 0, 0);
	FsrRcasH(Gamma2Color.r, Gamma2Color.g, Gamma2Color.b, pos, Const0);
	imageStore(fsr_image, ASU2(pos), AH4(Gamma2Color, 1));

#else

	AF3 Gamma2Color = AF3(0, 0, 0);
	FsrRcasF(Gamma2Color.r, Gamma2Color.g, Gamma2Color.b, pos, Const0);
	imageStore(fsr_image, ASU2(pos), AF4(Gamma2Color, 1));

#endif
}

void fsr_pass(AU2 pos) {
	if (params.pass == FSR_UPSCALE_PASS_TYPE_EASU) {
		fsr_easu_pass(pos);
	} else if (params.pass == FSR_UPSCALE_PASS_TYPE_RCAS) {
		fsr_rcas_pass(pos);
	}
}

void main() {
	// Clang does not like unused functions. If ffx_a.h is included in the binary, clang will throw a fit and not compile
	if (params.pass == FSR_UPSCALE_PASS_TYPE_EASU) {
		FsrEasuCon(Const0, Const1, Const2, Const3, params.resolution_width, params.resolution_height, params.resolution_width, params.resolution_height, params.upscaled_width, params.upscaled_height);
	} else if (params.pass == FSR_UPSCALE_PASS_TYPE_RCAS) {
		FsrRcasCon(Const0, params.sharpness);
	}

	AU2 gxy = ARmp8x8(gl_LocalInvocationID.x) + AU2(gl_WorkGroupID.x << 4u, gl_WorkGroupID.y << 4u);

	fsr_pass(gxy);
	gxy.x += 8u;
	fsr_pass(gxy);
	gxy.y += 8u;
	fsr_pass(gxy);
	gxy.x -= 8u;
	fsr_pass(gxy);
}
