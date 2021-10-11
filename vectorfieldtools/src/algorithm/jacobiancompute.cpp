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

/*
* Does (p2 - p1)/h
*/
vec3 JacobianCompute::forward_difference(const size3_t p2, const size3_t p1, const float h) {
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

/*
* Does (p2 - p1)/h
*/
vec3 JacobianCompute::backward_difference(const size3_t p2, const size3_t p1, const float h) {
	return forward_difference(p2,p1,h);
}

/*
* Does (p2 - p1)/(2h)
*/
vec3 JacobianCompute::central_difference(const size3_t p2, const size3_t p1, const float h) {
	return forward_difference(p2,p1, 2.0f*h);
}

/*
* gets the Jacobian at volume grid point pos = [x,y,z]
*/
mat3 JacobianCompute::get(const std::shared_ptr<const Volume> vector_field, const size3_t pos) {
	
	// store pointer to current vector field 
	curr_vector_field_ = vector_field;
	const size3_t field_dims = curr_vector_field_->getDimensions();
	const vec3 spacing = curr_vector_field_->getWorldSpaceGradientSpacing();
	// dx
	vec3 Fx(0,0,0);
	if(pos.x > 0 && pos.x < field_dims.x-1) {
		// central diff in x
		Fx = central_difference(pos + size3_t(1,0,0), pos - size3_t(1,0,0), spacing.x);
	} else if(pos.x == 0) {
		// forward diff in x
		Fx = forward_difference(pos + size3_t(1,0,0), pos, spacing.x);
	} else {
		// backward diff in x
		Fx = backward_difference(pos, pos - size3_t(1,0,0), spacing.x);
	}
	// dy
	vec3 Fy(0,0,0);
	if(pos.y > 0 && pos.y < field_dims.y-1) {
		// central diff in y
		Fy = central_difference(pos + size3_t(0,1,0), pos - size3_t(0,1,0), spacing.y);
	} else if(pos.y == 0) {
		// forward diff in y
		Fy = forward_difference(pos + size3_t(0,1,0), pos, spacing.y);
	} else {
		// backward diff in y
		Fy = backward_difference(pos, pos - size3_t(0,1,0), spacing.y);
	}
	// dz
	vec3 Fz(0,0,0);
	if(pos.z > 0 && pos.z < field_dims.z-1) {
		// central diff in z
		Fz = central_difference(pos + size3_t(0,0,1), pos - size3_t(0,0,1), spacing.z);
	} else if(pos.z == 0) {
		// forward diff in z
		Fz = forward_difference(pos + size3_t(0,0,1), pos, spacing.z);
	} else {
		// backward diff in z
		Fz = backward_difference(pos, pos - size3_t(0,0,1), spacing.z);
	}
	// store jacobian in column-major order... J = [Fx Fy Fz]
	mat3 jacobian_;
	// column 0
	jacobian_[0] = Fx;
	// column 1
	jacobian_[1] = Fy;
	// column 2
	jacobian_[2] = Fz;

	return jacobian_;
}

}  // namespace inviwo
