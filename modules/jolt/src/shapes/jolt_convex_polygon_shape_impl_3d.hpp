#pragma once

#include "shapes/jolt_shape_impl_3d.hpp"

class JoltConvexPolygonShapeImpl3D final : public JoltShapeImpl3D {
public:
	ShapeType get_type() const override { return ShapeType::SHAPE_CONVEX_POLYGON; }

	bool is_convex() const override { return true; }

	Variant get_data() const override;

	void set_data(const Variant& p_data) override;

	float get_margin() const override { return margin; }

	void set_margin(float p_margin) override;

	String to_string() const;

private:
	JPH::ShapeRefC _build() const override;

	PackedVector3Array vertices;

	float margin = 0.04f;
};
