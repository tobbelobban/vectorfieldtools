/*********************************************************************************
 *
 * Inviwo - Interactive Visualization Workshop
 *
 * Copyright (c) 2021 Inviwo Foundation
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *********************************************************************************/

#include <KTH/vectorfieldtools/algorithm/jacobiancompute.h>

namespace inviwo {

vec3 JacobianCompute::forward_difference(const size3_t p1, const size3_t p2, const float h) {
	const size3_t vec_dims = curr_vector_field_->getDimensions();
	const size_t xy_ = vec_dims.x * vec_dims.y;
	const size_t p1_index = xy_ * p1.z + vec_dims.x * p1.y + p1.x;
	const size_t p2_index = xy_ * p2.z + vec_dims.x * p2.y + p2.x;
	return curr_vector_field_->getRepresentation<VolumeRAM>()->dispatch<vec3, dispatching::filter::Float3s>(
		[h, p1_index, p2_index](auto vector_field_pr) {
			const auto vector_field_data = vector_field_pr->getDataTyped();
			return vec3(vector_field_data[p2_index] - vector_field_data[p1_index]) / h;
        }
    );
}

vec3 JacobianCompute::backward_difference(const size3_t p1, const size3_t p2, const float h) {
	const size3_t vec_dims = curr_vector_field_->getDimensions();
	const size_t xy_ = vec_dims.x * vec_dims.y;
	const size_t p1_index = xy_ * p1.z + vec_dims.x * p1.y + p1.x;
	const size_t p2_index = xy_ * p2.z + vec_dims.x * p2.y + p2.x;
	return curr_vector_field_->getRepresentation<VolumeRAM>()->dispatch<vec3, dispatching::filter::Float3s>(
		[h, p1_index, p2_index](auto vector_field_pr) {
			const auto vector_field_data = vector_field_pr->getDataTyped();
			return vec3(vector_field_data[p2_index] - vector_field_data[p1_index]) / h;
        }
    );
}

vec3 JacobianCompute::central_difference(const size3_t p1, const size3_t p2, const float h) {
	const size3_t vec_dims = curr_vector_field_->getDimensions();
	const size_t xy_ = vec_dims.x * vec_dims.y;
	const size_t p1_index = xy_ * p1.z + vec_dims.x * p1.y + p1.x;
	const size_t p2_index = xy_ * p2.z + vec_dims.x * p2.y + p2.x;
	return curr_vector_field_->getRepresentation<VolumeRAM>()->dispatch<vec3, dispatching::filter::Float3s>(
		[h, p1_index, p2_index](auto vector_field_pr) {
			const auto vector_field_data = vector_field_pr->getDataTyped();
			return vec3(vector_field_data[p2_index] - vector_field_data[p1_index]) / (2.0f * h);
        }
    );
}

mat3 JacobianCompute::get(const std::shared_ptr<const Volume> vector_field, const size3_t pos) {
	curr_vector_field_ = vector_field;
	curr_pos_ = pos;
	const size3_t field_dims = curr_vector_field_->getDimensions();
	const vec3 spacing = curr_vector_field_->getWorldSpaceGradientSpacing();
	// dx
	vec3 Fx(0,0,0);
	if(curr_pos_.x > 0 && curr_pos_.x < field_dims.x-1) {
		// central diff in x
		Fx = central_difference(curr_pos_ + size3_t(1,0,0), curr_pos_ - size3_t(1,0,0), spacing.x);
	} else if(curr_pos_.x == 0) {
		// forward diff in x
		Fx = forward_difference(curr_pos_, curr_pos_ + size3_t(1,0,0), spacing.x);
	} else {
		// backward diff in x
		Fx = backward_difference(curr_pos_ - size3_t(1,0,0), curr_pos_, spacing.x);
	}
	// dy
	vec3 Fy(0,0,0);
	if(curr_pos_.y > 0 && curr_pos_.y < field_dims.y-1) {
		// central diff in y
		Fy = central_difference(curr_pos_ + size3_t(0,1,0), curr_pos_ - size3_t(0,1,0), spacing.y);
	} else if(curr_pos_.y == 0) {
		// forward diff in y
		Fy = forward_difference(curr_pos_, curr_pos_ + size3_t(0,1,0), spacing.y);
	} else {
		// backward diff in y
		Fy = backward_difference(curr_pos_ - size3_t(0,1,0), curr_pos_, spacing.y);
	}
	// dz
	vec3 Fz(0,0,0);
	if(curr_pos_.z > 0 && curr_pos_.z < field_dims.z-1) {
		// central diff in z
		Fz = central_difference(curr_pos_ + size3_t(0,0,1), curr_pos_ - size3_t(0,0,1), spacing.z);
	} else if(curr_pos_.z == 0) {
		// forward diff in z
		Fz = forward_difference(curr_pos_, curr_pos_ + size3_t(0,0,1), spacing.z);
	} else {
		// backward diff in z
		Fz = backward_difference(curr_pos_ - size3_t(0,0,1), curr_pos_, spacing.z);
	}
	// store jacobian in column-major order...
	mat3 jacobian_;
	// column 0
	jacobian_[0] = vec3(Fx.x, Fy.x, Fz.x);
	// column 1
	jacobian_[1] = vec3(Fx.y, Fy.y, Fz.y);
	// column 2
	jacobian_[2] = vec3(Fx.z, Fy.z, Fz.z);

	return jacobian_;
}

}  // namespace inviwo
