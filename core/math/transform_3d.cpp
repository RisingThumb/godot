/*************************************************************************/
/*  transform_3d.cpp                                                     */
/*************************************************************************/
/*                       This file is part of:                           */
/*                           GODOT ENGINE                                */
/*                      https://godotengine.org                          */
/*************************************************************************/
/* Copyright (c) 2007-2022 Juan Linietsky, Ariel Manzur.                 */
/* Copyright (c) 2014-2022 Godot Engine contributors (cf. AUTHORS.md).   */
/*                                                                       */
/* Permission is hereby granted, free of charge, to any person obtaining */
/* a copy of this software and associated documentation files (the       */
/* "Software"), to deal in the Software without restriction, including   */
/* without limitation the rights to use, copy, modify, merge, publish,   */
/* distribute, sublicense, and/or sell copies of the Software, and to    */
/* permit persons to whom the Software is furnished to do so, subject to */
/* the following conditions:                                             */
/*                                                                       */
/* The above copyright notice and this permission notice shall be        */
/* included in all copies or substantial portions of the Software.       */
/*                                                                       */
/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,       */
/* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF    */
/* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.*/
/* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY  */
/* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,  */
/* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE     */
/* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                */
/*************************************************************************/

#include "transform_3d.h"

#include "core/math/math_funcs.h"
#include "core/string/print_string.h"

void Transform3D::affine_invert() {
	basis.invert();
	origin = basis.xform(-origin);
}

Transform3D Transform3D::affine_inverse() const {
	Transform3D ret = *this;
	ret.affine_invert();
	return ret;
}

void Transform3D::invert() {
	basis.transpose();
	origin = basis.xform(-origin);
}

Transform3D Transform3D::inverse() const {
	// FIXME: this function assumes the basis is a rotation matrix, with no scaling.
	// Transform3D::affine_inverse can handle matrices with scaling, so GDScript should eventually use that.
	Transform3D ret = *this;
	ret.invert();
	return ret;
}

void Transform3D::rotate(const Vector3 &p_axis, real_t p_phi) {
	*this = rotated(p_axis, p_phi);
}

Transform3D Transform3D::rotated(const Vector3 &p_axis, real_t p_phi) const {
	return Transform3D(Basis(p_axis, p_phi), Vector3()) * (*this);
}

void Transform3D::rotate_basis(const Vector3 &p_axis, real_t p_phi) {
	basis.rotate(p_axis, p_phi);
}

Transform3D Transform3D::looking_at(const Vector3 &p_target, const Vector3 &p_up) const {
	Transform3D t = *this;
	t.basis = Basis::looking_at(p_target - origin, p_up);
	return t;
}

void Transform3D::set_look_at(const Vector3 &p_eye, const Vector3 &p_target, const Vector3 &p_up) {
	basis = Basis::looking_at(p_target - p_eye, p_up);
	origin = p_eye;
}

Transform3D Transform3D::sphere_interpolate_with(const Transform3D &p_transform, real_t p_c) const {
	// If you plan on calling sphere_interpolate_with with the same two transform but for different values of p_c.
	// The following terms can be cached for use by the last block of code.

	Transform3D inverse = this->affine_inverse();
	Transform3D h = p_transform * inverse;
	Basis s = h.basis.log();
	real_t s_0 = s.elements[2][1];
	real_t s_1 = s.elements[0][2];
	real_t s_2 = s.elements[1][0];
	real_t theta = Math::sqrt(s_0 * s_0 + s_1 * s_1 + s_2 * s_2);
	Basis inv_v_1 = s._compute_inverse_v_1(theta);
	Vector3 u = inv_v_1.xform(h.origin);
	// End section that can be cached.

	Basis interp_r = s.exp(p_c, theta);
	Basis interp_t_t_times_v = s._compute_t_times_v(theta, p_c);
	Transform3D interp_h;
	interp_h.basis = interp_r * basis;
	interp_h.origin = interp_r.xform(origin) + interp_t_t_times_v.xform(u);
	return interp_h;
}

Transform3D Transform3D::interpolate_with(const Transform3D &p_transform, real_t p_c) const {
	Transform3D interp;

	interp.basis = basis.lerp(p_transform.basis, p_c);
	interp.origin = origin.lerp(p_transform.origin, p_c);

	return interp;
}

void Transform3D::scale(const Vector3 &p_scale) {
	basis.scale(p_scale);
	origin *= p_scale;
}

Transform3D Transform3D::scaled(const Vector3 &p_scale) const {
	Transform3D t = *this;
	t.scale(p_scale);
	return t;
}

void Transform3D::scale_basis(const Vector3 &p_scale) {
	basis.scale(p_scale);
}

void Transform3D::translate(real_t p_tx, real_t p_ty, real_t p_tz) {
	translate(Vector3(p_tx, p_ty, p_tz));
}

void Transform3D::translate(const Vector3 &p_translation) {
	for (int i = 0; i < 3; i++) {
		origin[i] += basis[i].dot(p_translation);
	}
}

Transform3D Transform3D::translated(const Vector3 &p_translation) const {
	Transform3D t = *this;
	t.translate(p_translation);
	return t;
}

void Transform3D::orthonormalize() {
	basis.orthonormalize();
}

Transform3D Transform3D::orthonormalized() const {
	Transform3D _copy = *this;
	_copy.orthonormalize();
	return _copy;
}

void Transform3D::orthogonalize() {
	basis.orthogonalize();
}

Transform3D Transform3D::orthogonalized() const {
	Transform3D _copy = *this;
	_copy.orthogonalize();
	return _copy;
}

bool Transform3D::is_equal_approx(const Transform3D &p_transform) const {
	return basis.is_equal_approx(p_transform.basis) && origin.is_equal_approx(p_transform.origin);
}

bool Transform3D::operator==(const Transform3D &p_transform) const {
	return (basis == p_transform.basis && origin == p_transform.origin);
}

bool Transform3D::operator!=(const Transform3D &p_transform) const {
	return (basis != p_transform.basis || origin != p_transform.origin);
}

void Transform3D::operator*=(const Transform3D &p_transform) {
	origin = xform(p_transform.origin);
	basis *= p_transform.basis;
}

Transform3D Transform3D::operator*(const Transform3D &p_transform) const {
	Transform3D t = *this;
	t *= p_transform;
	return t;
}

void Transform3D::operator*=(const real_t p_val) {
	origin *= p_val;
	basis *= p_val;
}

Transform3D Transform3D::operator*(const real_t p_val) const {
	Transform3D ret(*this);
	ret *= p_val;
	return ret;
}

Transform3D::operator String() const {
	return "[X: " + basis.get_axis(0).operator String() +
			", Y: " + basis.get_axis(1).operator String() +
			", Z: " + basis.get_axis(2).operator String() +
			", O: " + origin.operator String() + "]";
}

Transform3D::Transform3D(const Basis &p_basis, const Vector3 &p_origin) :
		basis(p_basis),
		origin(p_origin) {
}

Transform3D::Transform3D(const Vector3 &p_x, const Vector3 &p_y, const Vector3 &p_z, const Vector3 &p_origin) :
		origin(p_origin) {
	basis.set_axis(0, p_x);
	basis.set_axis(1, p_y);
	basis.set_axis(2, p_z);
}

Transform3D::Transform3D(real_t xx, real_t xy, real_t xz, real_t yx, real_t yy, real_t yz, real_t zx, real_t zy, real_t zz, real_t ox, real_t oy, real_t oz) {
	basis = Basis(xx, xy, xz, yx, yy, yz, zx, zy, zz);
	origin = Vector3(ox, oy, oz);
}

Transform3D Transform3D::cubic_interpolate(const Transform3D &p_b, const Transform3D &p_pre_a, const Transform3D &p_post_b, const real_t &p_weight) const {
	real_t t2 = (1.0 - p_weight) * p_weight * 2;
	Transform3D sp = p_pre_a.sphere_interpolate_with(p_post_b, p_weight);
	Transform3D sq = this->sphere_interpolate_with(p_b, p_weight);
	return sp.sphere_interpolate_with(sq, t2);
}