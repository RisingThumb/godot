/// @ref ext_matrix_clip_space
/// @file glm/ext/matrix_clip_space.hpp
///
/// @defgroup ext_matrix_clip_space GLM_EXT_matrix_clip_space
/// @ingroup ext
///
/// Defines functions that generate clip space transformation matrices.
///
/// The matrices generated by this extension use standard OpenGL fixed-function
/// conventions. For example, the lookAt function generates a transform from world
/// space into the specific eye space that the projective matrix functions
/// (perspective, ortho, etc) are designed to expect. The OpenGL compatibility
/// specifications defines the particular layout of this eye space.
///
/// Include <glm/ext/matrix_clip_space.hpp> to use the features of this extension.
///
/// @see ext_matrix_transform
/// @see ext_matrix_projection

#pragma once

// Dependencies
#include "../ext/scalar_constants.hpp"
#include "../geometric.hpp"
#include "../trigonometric.hpp"

#if GLM_MESSAGES == GLM_ENABLE && !defined(GLM_EXT_INCLUDED)
#	pragma message("GLM: GLM_EXT_matrix_clip_space extension included")
#endif

namespace glm
{
	/// @addtogroup ext_matrix_clip_space
	/// @{

	/// Creates a matrix for projecting two-dimensional coordinates onto the screen.
	///
	/// @tparam T A floating-point scalar type
	///
	/// @see - glm::ortho(T const& left, T const& right, T const& bottom, T const& top, T const& zNear, T const& zFar)
	/// @see <a href="https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/gluOrtho2D.xml">gluOrtho2D man page</a>
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> ortho(
		T left, T right, T bottom, T top);

	/// Creates a matrix for an orthographic parallel viewing volume, using left-handed coordinates.
	/// The near and far clip planes correspond to z normalized device coordinates of 0 and +1 respectively. (Direct3D clip volume definition)
	///
	/// @tparam T A floating-point scalar type
	///
	/// @see - glm::ortho(T const& left, T const& right, T const& bottom, T const& top)
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> orthoLH_ZO(
		T left, T right, T bottom, T top, T zNear, T zFar);

	/// Creates a matrix for an orthographic parallel viewing volume using right-handed coordinates.
	/// The near and far clip planes correspond to z normalized device coordinates of -1 and +1 respectively. (OpenGL clip volume definition)
	///
	/// @tparam T A floating-point scalar type
	///
	/// @see - glm::ortho(T const& left, T const& right, T const& bottom, T const& top)
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> orthoLH_NO(
		T left, T right, T bottom, T top, T zNear, T zFar);

	/// Creates a matrix for an orthographic parallel viewing volume, using left-handed coordinates.
	/// The near and far clip planes correspond to z normalized device coordinates of 0 and +1 respectively. (Direct3D clip volume definition)
	///
	/// @tparam T A floating-point scalar type
	///
	/// @see - glm::ortho(T const& left, T const& right, T const& bottom, T const& top)
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> orthoRH_ZO(
		T left, T right, T bottom, T top, T zNear, T zFar);

	/// Creates a matrix for an orthographic parallel viewing volume, using right-handed coordinates.
	/// The near and far clip planes correspond to z normalized device coordinates of -1 and +1 respectively. (OpenGL clip volume definition)
	///
	/// @tparam T A floating-point scalar type
	///
	/// @see - glm::ortho(T const& left, T const& right, T const& bottom, T const& top)
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> orthoRH_NO(
		T left, T right, T bottom, T top, T zNear, T zFar);

	/// Creates a matrix for an orthographic parallel viewing volume, using left-handed coordinates.
	/// The near and far clip planes correspond to z normalized device coordinates of 0 and +1 respectively. (Direct3D clip volume definition)
	///
	/// @tparam T A floating-point scalar type
	///
	/// @see - glm::ortho(T const& left, T const& right, T const& bottom, T const& top)
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> orthoZO(
		T left, T right, T bottom, T top, T zNear, T zFar);

	/// Creates a matrix for an orthographic parallel viewing volume, using left-handed coordinates if GLM_FORCE_LEFT_HANDED if defined or right-handed coordinates otherwise.
	/// The near and far clip planes correspond to z normalized device coordinates of -1 and +1 respectively. (OpenGL clip volume definition)
	///
	/// @tparam T A floating-point scalar type
	///
	/// @see - glm::ortho(T const& left, T const& right, T const& bottom, T const& top)
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> orthoNO(
		T left, T right, T bottom, T top, T zNear, T zFar);

	/// Creates a matrix for an orthographic parallel viewing volume, using left-handed coordinates.
	/// If GLM_FORCE_DEPTH_ZERO_TO_ONE is defined, the near and far clip planes correspond to z normalized device coordinates of 0 and +1 respectively. (Direct3D clip volume definition)
	/// Otherwise, the near and far clip planes correspond to z normalized device coordinates of -1 and +1 respectively. (OpenGL clip volume definition)
	///
	/// @tparam T A floating-point scalar type
	///
	/// @see - glm::ortho(T const& left, T const& right, T const& bottom, T const& top)
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> orthoLH(
		T left, T right, T bottom, T top, T zNear, T zFar);

	/// Creates a matrix for an orthographic parallel viewing volume, using right-handed coordinates.
	/// If GLM_FORCE_DEPTH_ZERO_TO_ONE is defined, the near and far clip planes correspond to z normalized device coordinates of 0 and +1 respectively. (Direct3D clip volume definition)
	/// Otherwise, the near and far clip planes correspond to z normalized device coordinates of -1 and +1 respectively. (OpenGL clip volume definition)
	///
	/// @tparam T A floating-point scalar type
	///
	/// @see - glm::ortho(T const& left, T const& right, T const& bottom, T const& top)
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> orthoRH(
		T left, T right, T bottom, T top, T zNear, T zFar);

	/// Creates a matrix for an orthographic parallel viewing volume, using the default handedness and default near and far clip planes definition.
	/// To change default handedness use GLM_FORCE_LEFT_HANDED. To change default near and far clip planes definition use GLM_FORCE_DEPTH_ZERO_TO_ONE.
	///
	/// @tparam T A floating-point scalar type
	///
	/// @see - glm::ortho(T const& left, T const& right, T const& bottom, T const& top)
	/// @see <a href="https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/glOrtho.xml">glOrtho man page</a>
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> ortho(
		T left, T right, T bottom, T top, T zNear, T zFar);

	/// Creates a left handed frustum matrix.
	/// The near and far clip planes correspond to z normalized device coordinates of 0 and +1 respectively. (Direct3D clip volume definition)
	///
	/// @tparam T A floating-point scalar type
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> frustumLH_ZO(
		T left, T right, T bottom, T top, T near, T far);

	/// Creates a left handed frustum matrix.
	/// The near and far clip planes correspond to z normalized device coordinates of -1 and +1 respectively. (OpenGL clip volume definition)
	///
	/// @tparam T A floating-point scalar type
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> frustumLH_NO(
		T left, T right, T bottom, T top, T near, T far);

	/// Creates a right handed frustum matrix.
	/// The near and far clip planes correspond to z normalized device coordinates of 0 and +1 respectively. (Direct3D clip volume definition)
	///
	/// @tparam T A floating-point scalar type
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> frustumRH_ZO(
		T left, T right, T bottom, T top, T near, T far);

	/// Creates a right handed frustum matrix.
	/// The near and far clip planes correspond to z normalized device coordinates of -1 and +1 respectively. (OpenGL clip volume definition)
	///
	/// @tparam T A floating-point scalar type
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> frustumRH_NO(
		T left, T right, T bottom, T top, T near, T far);

	/// Creates a frustum matrix using left-handed coordinates if GLM_FORCE_LEFT_HANDED if defined or right-handed coordinates otherwise.
	/// The near and far clip planes correspond to z normalized device coordinates of 0 and +1 respectively. (Direct3D clip volume definition)
	///
	/// @tparam T A floating-point scalar type
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> frustumZO(
		T left, T right, T bottom, T top, T near, T far);

	/// Creates a frustum matrix using left-handed coordinates if GLM_FORCE_LEFT_HANDED if defined or right-handed coordinates otherwise.
	/// The near and far clip planes correspond to z normalized device coordinates of -1 and +1 respectively. (OpenGL clip volume definition)
	///
	/// @tparam T A floating-point scalar type
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> frustumNO(
		T left, T right, T bottom, T top, T near, T far);

	/// Creates a left handed frustum matrix.
	/// If GLM_FORCE_DEPTH_ZERO_TO_ONE is defined, the near and far clip planes correspond to z normalized device coordinates of 0 and +1 respectively. (Direct3D clip volume definition)
	/// Otherwise, the near and far clip planes correspond to z normalized device coordinates of -1 and +1 respectively. (OpenGL clip volume definition)
	///
	/// @tparam T A floating-point scalar type
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> frustumLH(
		T left, T right, T bottom, T top, T near, T far);

	/// Creates a right handed frustum matrix.
	/// If GLM_FORCE_DEPTH_ZERO_TO_ONE is defined, the near and far clip planes correspond to z normalized device coordinates of 0 and +1 respectively. (Direct3D clip volume definition)
	/// Otherwise, the near and far clip planes correspond to z normalized device coordinates of -1 and +1 respectively. (OpenGL clip volume definition)
	///
	/// @tparam T A floating-point scalar type
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> frustumRH(
		T left, T right, T bottom, T top, T near, T far);

	/// Creates a frustum matrix with default handedness, using the default handedness and default near and far clip planes definition.
	/// To change default handedness use GLM_FORCE_LEFT_HANDED. To change default near and far clip planes definition use GLM_FORCE_DEPTH_ZERO_TO_ONE.
	///
	/// @tparam T A floating-point scalar type
	/// @see <a href="https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/glFrustum.xml">glFrustum man page</a>
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> frustum(
		T left, T right, T bottom, T top, T near, T far);


	/// Creates a matrix for a right handed, symmetric perspective-view frustum.
	/// The near and far clip planes correspond to z normalized device coordinates of 0 and +1 respectively. (Direct3D clip volume definition)
	///
	/// @param fovy Specifies the field of view angle, in degrees, in the y direction. Expressed in radians.
	/// @param aspect Specifies the aspect ratio that determines the field of view in the x direction. The aspect ratio is the ratio of x (width) to y (height).
	/// @param near Specifies the distance from the viewer to the near clipping plane (always positive).
	/// @param far Specifies the distance from the viewer to the far clipping plane (always positive).
	///
	/// @tparam T A floating-point scalar type
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> perspectiveRH_ZO(
		T fovy, T aspect, T near, T far);

	/// Creates a matrix for a right handed, symmetric perspective-view frustum.
	/// The near and far clip planes correspond to z normalized device coordinates of -1 and +1 respectively. (OpenGL clip volume definition)
	///
	/// @param fovy Specifies the field of view angle, in degrees, in the y direction. Expressed in radians.
	/// @param aspect Specifies the aspect ratio that determines the field of view in the x direction. The aspect ratio is the ratio of x (width) to y (height).
	/// @param near Specifies the distance from the viewer to the near clipping plane (always positive).
	/// @param far Specifies the distance from the viewer to the far clipping plane (always positive).
	///
	/// @tparam T A floating-point scalar type
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> perspectiveRH_NO(
		T fovy, T aspect, T near, T far);

	/// Creates a matrix for a left handed, symmetric perspective-view frustum.
	/// The near and far clip planes correspond to z normalized device coordinates of 0 and +1 respectively. (Direct3D clip volume definition)
	///
	/// @param fovy Specifies the field of view angle, in degrees, in the y direction. Expressed in radians.
	/// @param aspect Specifies the aspect ratio that determines the field of view in the x direction. The aspect ratio is the ratio of x (width) to y (height).
	/// @param near Specifies the distance from the viewer to the near clipping plane (always positive).
	/// @param far Specifies the distance from the viewer to the far clipping plane (always positive).
	///
	/// @tparam T A floating-point scalar type
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> perspectiveLH_ZO(
		T fovy, T aspect, T near, T far);

	/// Creates a matrix for a left handed, symmetric perspective-view frustum.
	/// The near and far clip planes correspond to z normalized device coordinates of -1 and +1 respectively. (OpenGL clip volume definition)
	///
	/// @param fovy Specifies the field of view angle, in degrees, in the y direction. Expressed in radians.
	/// @param aspect Specifies the aspect ratio that determines the field of view in the x direction. The aspect ratio is the ratio of x (width) to y (height).
	/// @param near Specifies the distance from the viewer to the near clipping plane (always positive).
	/// @param far Specifies the distance from the viewer to the far clipping plane (always positive).
	///
	/// @tparam T A floating-point scalar type
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> perspectiveLH_NO(
		T fovy, T aspect, T near, T far);

	/// Creates a matrix for a symmetric perspective-view frustum using left-handed coordinates if GLM_FORCE_LEFT_HANDED if defined or right-handed coordinates otherwise.
	/// The near and far clip planes correspond to z normalized device coordinates of 0 and +1 respectively. (Direct3D clip volume definition)
	///
	/// @param fovy Specifies the field of view angle, in degrees, in the y direction. Expressed in radians.
	/// @param aspect Specifies the aspect ratio that determines the field of view in the x direction. The aspect ratio is the ratio of x (width) to y (height).
	/// @param near Specifies the distance from the viewer to the near clipping plane (always positive).
	/// @param far Specifies the distance from the viewer to the far clipping plane (always positive).
	///
	/// @tparam T A floating-point scalar type
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> perspectiveZO(
		T fovy, T aspect, T near, T far);

	/// Creates a matrix for a symmetric perspective-view frustum using left-handed coordinates if GLM_FORCE_LEFT_HANDED if defined or right-handed coordinates otherwise.
	/// The near and far clip planes correspond to z normalized device coordinates of -1 and +1 respectively. (OpenGL clip volume definition)
	///
	/// @param fovy Specifies the field of view angle, in degrees, in the y direction. Expressed in radians.
	/// @param aspect Specifies the aspect ratio that determines the field of view in the x direction. The aspect ratio is the ratio of x (width) to y (height).
	/// @param near Specifies the distance from the viewer to the near clipping plane (always positive).
	/// @param far Specifies the distance from the viewer to the far clipping plane (always positive).
	///
	/// @tparam T A floating-point scalar type
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> perspectiveNO(
		T fovy, T aspect, T near, T far);

	/// Creates a matrix for a right handed, symmetric perspective-view frustum.
	/// If GLM_FORCE_DEPTH_ZERO_TO_ONE is defined, the near and far clip planes correspond to z normalized device coordinates of 0 and +1 respectively. (Direct3D clip volume definition)
	/// Otherwise, the near and far clip planes correspond to z normalized device coordinates of -1 and +1 respectively. (OpenGL clip volume definition)
	///
	/// @param fovy Specifies the field of view angle, in degrees, in the y direction. Expressed in radians.
	/// @param aspect Specifies the aspect ratio that determines the field of view in the x direction. The aspect ratio is the ratio of x (width) to y (height).
	/// @param near Specifies the distance from the viewer to the near clipping plane (always positive).
	/// @param far Specifies the distance from the viewer to the far clipping plane (always positive).
	///
	/// @tparam T A floating-point scalar type
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> perspectiveRH(
		T fovy, T aspect, T near, T far);

	/// Creates a matrix for a left handed, symmetric perspective-view frustum.
	/// If GLM_FORCE_DEPTH_ZERO_TO_ONE is defined, the near and far clip planes correspond to z normalized device coordinates of 0 and +1 respectively. (Direct3D clip volume definition)
	/// Otherwise, the near and far clip planes correspond to z normalized device coordinates of -1 and +1 respectively. (OpenGL clip volume definition)
	///
	/// @param fovy Specifies the field of view angle, in degrees, in the y direction. Expressed in radians.
	/// @param aspect Specifies the aspect ratio that determines the field of view in the x direction. The aspect ratio is the ratio of x (width) to y (height).
	/// @param near Specifies the distance from the viewer to the near clipping plane (always positive).
	/// @param far Specifies the distance from the viewer to the far clipping plane (always positive).
	///
	/// @tparam T A floating-point scalar type
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> perspectiveLH(
		T fovy, T aspect, T near, T far);

	/// Creates a matrix for a symmetric perspective-view frustum based on the default handedness and default near and far clip planes definition.
	/// To change default handedness use GLM_FORCE_LEFT_HANDED. To change default near and far clip planes definition use GLM_FORCE_DEPTH_ZERO_TO_ONE.
	///
	/// @param fovy Specifies the field of view angle in the y direction. Expressed in radians.
	/// @param aspect Specifies the aspect ratio that determines the field of view in the x direction. The aspect ratio is the ratio of x (width) to y (height).
	/// @param near Specifies the distance from the viewer to the near clipping plane (always positive).
	/// @param far Specifies the distance from the viewer to the far clipping plane (always positive).
	///
	/// @tparam T A floating-point scalar type
	/// @see <a href="https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/gluPerspective.xml">gluPerspective man page</a>
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> perspective(
		T fovy, T aspect, T near, T far);

	/// Builds a perspective projection matrix based on a field of view using right-handed coordinates.
	/// The near and far clip planes correspond to z normalized device coordinates of 0 and +1 respectively. (Direct3D clip volume definition)
	///
	/// @param fov Expressed in radians.
	/// @param width Width of the viewport
	/// @param height Height of the viewport
	/// @param near Specifies the distance from the viewer to the near clipping plane (always positive).
	/// @param far Specifies the distance from the viewer to the far clipping plane (always positive).
	///
	/// @tparam T A floating-point scalar type
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> perspectiveFovRH_ZO(
		T fov, T width, T height, T near, T far);

	/// Builds a perspective projection matrix based on a field of view using right-handed coordinates.
	/// The near and far clip planes correspond to z normalized device coordinates of -1 and +1 respectively. (OpenGL clip volume definition)
	///
	/// @param fov Expressed in radians.
	/// @param width Width of the viewport
	/// @param height Height of the viewport
	/// @param near Specifies the distance from the viewer to the near clipping plane (always positive).
	/// @param far Specifies the distance from the viewer to the far clipping plane (always positive).
	///
	/// @tparam T A floating-point scalar type
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> perspectiveFovRH_NO(
		T fov, T width, T height, T near, T far);

	/// Builds a perspective projection matrix based on a field of view using left-handed coordinates.
	/// The near and far clip planes correspond to z normalized device coordinates of 0 and +1 respectively. (Direct3D clip volume definition)
	///
	/// @param fov Expressed in radians.
	/// @param width Width of the viewport
	/// @param height Height of the viewport
	/// @param near Specifies the distance from the viewer to the near clipping plane (always positive).
	/// @param far Specifies the distance from the viewer to the far clipping plane (always positive).
	///
	/// @tparam T A floating-point scalar type
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> perspectiveFovLH_ZO(
		T fov, T width, T height, T near, T far);

	/// Builds a perspective projection matrix based on a field of view using left-handed coordinates.
	/// The near and far clip planes correspond to z normalized device coordinates of -1 and +1 respectively. (OpenGL clip volume definition)
	///
	/// @param fov Expressed in radians.
	/// @param width Width of the viewport
	/// @param height Height of the viewport
	/// @param near Specifies the distance from the viewer to the near clipping plane (always positive).
	/// @param far Specifies the distance from the viewer to the far clipping plane (always positive).
	///
	/// @tparam T A floating-point scalar type
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> perspectiveFovLH_NO(
		T fov, T width, T height, T near, T far);

	/// Builds a perspective projection matrix based on a field of view using left-handed coordinates if GLM_FORCE_LEFT_HANDED if defined or right-handed coordinates otherwise.
	/// The near and far clip planes correspond to z normalized device coordinates of 0 and +1 respectively. (Direct3D clip volume definition)
	///
	/// @param fov Expressed in radians.
	/// @param width Width of the viewport
	/// @param height Height of the viewport
	/// @param near Specifies the distance from the viewer to the near clipping plane (always positive).
	/// @param far Specifies the distance from the viewer to the far clipping plane (always positive).
	///
	/// @tparam T A floating-point scalar type
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> perspectiveFovZO(
		T fov, T width, T height, T near, T far);

	/// Builds a perspective projection matrix based on a field of view using left-handed coordinates if GLM_FORCE_LEFT_HANDED if defined or right-handed coordinates otherwise.
	/// The near and far clip planes correspond to z normalized device coordinates of -1 and +1 respectively. (OpenGL clip volume definition)
	///
	/// @param fov Expressed in radians.
	/// @param width Width of the viewport
	/// @param height Height of the viewport
	/// @param near Specifies the distance from the viewer to the near clipping plane (always positive).
	/// @param far Specifies the distance from the viewer to the far clipping plane (always positive).
	///
	/// @tparam T A floating-point scalar type
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> perspectiveFovNO(
		T fov, T width, T height, T near, T far);

	/// Builds a right handed perspective projection matrix based on a field of view.
	/// If GLM_FORCE_DEPTH_ZERO_TO_ONE is defined, the near and far clip planes correspond to z normalized device coordinates of 0 and +1 respectively. (Direct3D clip volume definition)
	/// Otherwise, the near and far clip planes correspond to z normalized device coordinates of -1 and +1 respectively. (OpenGL clip volume definition)
	///
	/// @param fov Expressed in radians.
	/// @param width Width of the viewport
	/// @param height Height of the viewport
	/// @param near Specifies the distance from the viewer to the near clipping plane (always positive).
	/// @param far Specifies the distance from the viewer to the far clipping plane (always positive).
	///
	/// @tparam T A floating-point scalar type
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> perspectiveFovRH(
		T fov, T width, T height, T near, T far);

	/// Builds a left handed perspective projection matrix based on a field of view.
	/// If GLM_FORCE_DEPTH_ZERO_TO_ONE is defined, the near and far clip planes correspond to z normalized device coordinates of 0 and +1 respectively. (Direct3D clip volume definition)
	/// Otherwise, the near and far clip planes correspond to z normalized device coordinates of -1 and +1 respectively. (OpenGL clip volume definition)
	///
	/// @param fov Expressed in radians.
	/// @param width Width of the viewport
	/// @param height Height of the viewport
	/// @param near Specifies the distance from the viewer to the near clipping plane (always positive).
	/// @param far Specifies the distance from the viewer to the far clipping plane (always positive).
	///
	/// @tparam T A floating-point scalar type
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> perspectiveFovLH(
		T fov, T width, T height, T near, T far);

	/// Builds a perspective projection matrix based on a field of view and the default handedness and default near and far clip planes definition.
	/// To change default handedness use GLM_FORCE_LEFT_HANDED. To change default near and far clip planes definition use GLM_FORCE_DEPTH_ZERO_TO_ONE.
	///
	/// @param fov Expressed in radians.
	/// @param width Width of the viewport
	/// @param height Height of the viewport
	/// @param near Specifies the distance from the viewer to the near clipping plane (always positive).
	/// @param far Specifies the distance from the viewer to the far clipping plane (always positive).
	///
	/// @tparam T A floating-point scalar type
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> perspectiveFov(
		T fov, T width, T height, T near, T far);

	/// Creates a matrix for a left handed, symmetric perspective-view frustum with far plane at infinite.
	///
	/// @param fovy Specifies the field of view angle, in degrees, in the y direction. Expressed in radians.
	/// @param aspect Specifies the aspect ratio that determines the field of view in the x direction. The aspect ratio is the ratio of x (width) to y (height).
	/// @param near Specifies the distance from the viewer to the near clipping plane (always positive).
	///
	/// @tparam T A floating-point scalar type
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> infinitePerspectiveLH(
		T fovy, T aspect, T near);

	/// Creates a matrix for a right handed, symmetric perspective-view frustum with far plane at infinite.
	///
	/// @param fovy Specifies the field of view angle, in degrees, in the y direction. Expressed in radians.
	/// @param aspect Specifies the aspect ratio that determines the field of view in the x direction. The aspect ratio is the ratio of x (width) to y (height).
	/// @param near Specifies the distance from the viewer to the near clipping plane (always positive).
	///
	/// @tparam T A floating-point scalar type
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> infinitePerspectiveRH(
		T fovy, T aspect, T near);

	/// Creates a matrix for a symmetric perspective-view frustum with far plane at infinite with default handedness.
	///
	/// @param fovy Specifies the field of view angle, in degrees, in the y direction. Expressed in radians.
	/// @param aspect Specifies the aspect ratio that determines the field of view in the x direction. The aspect ratio is the ratio of x (width) to y (height).
	/// @param near Specifies the distance from the viewer to the near clipping plane (always positive).
	///
	/// @tparam T A floating-point scalar type
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> infinitePerspective(
		T fovy, T aspect, T near);

	/// Creates a matrix for a symmetric perspective-view frustum with far plane at infinite for graphics hardware that doesn't support depth clamping.
	///
	/// @param fovy Specifies the field of view angle, in degrees, in the y direction. Expressed in radians.
	/// @param aspect Specifies the aspect ratio that determines the field of view in the x direction. The aspect ratio is the ratio of x (width) to y (height).
	/// @param near Specifies the distance from the viewer to the near clipping plane (always positive).
	///
	/// @tparam T A floating-point scalar type
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> tweakedInfinitePerspective(
		T fovy, T aspect, T near);

	/// Creates a matrix for a symmetric perspective-view frustum with far plane at infinite for graphics hardware that doesn't support depth clamping.
	///
	/// @param fovy Specifies the field of view angle, in degrees, in the y direction. Expressed in radians.
	/// @param aspect Specifies the aspect ratio that determines the field of view in the x direction. The aspect ratio is the ratio of x (width) to y (height).
	/// @param near Specifies the distance from the viewer to the near clipping plane (always positive).
	/// @param ep Epsilon
	///
	/// @tparam T A floating-point scalar type
	template<typename T>
	GLM_FUNC_DECL mat<4, 4, T, defaultp> tweakedInfinitePerspective(
		T fovy, T aspect, T near, T ep);

	/// @}
}//namespace glm

#include "matrix_clip_space.inl"
