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

#include <KTH/vectorfieldtools/algorithm/omega2s2.h>

namespace inviwo {

mat3 Omega2S2::get(const std::shared_ptr<const Volume> vector_field, const size3_t pos) {	
	// compute jacobian at pos
	const mat3 jacobian = jacobian_computer_.get(vector_field, pos);
	
	// S = .5(J + transpose(J))
	mat3 S;
	
	// column 0
	S[0] = vec3(	jacobian[0][0],						// u_x
			.5f * (jacobian[0][1] + jacobian[1][0]),		// .5(v_x + u_y)
			.5f * (jacobian[0][2] + jacobian[2][0])		);	// .5(w_x + u_z)
	// column 1
	S[1] = vec3(	.5f * (jacobian[1][0] + jacobian[0][1]),		// .5(u_y + v_x)
			jacobian[1][1],						// v_y
			.5f * (jacobian[1][2] + jacobian[2][1])		);	// .5(w_y + v_z)
	// column 2
	S[2] = vec3(	.5f * (jacobian[2][0] + jacobian[0][2]),		// .5(u_z + w_x)
			.5f * (jacobian[2][1] + jacobian[1][2]),		// .5(v_z + w_y)
			jacobian[2][2]					);	// w_z
	
	// Omega = .5(J - transpose(J))
	mat3 Omega;		
	
	// column 0
	Omega[0] = vec3(0.0f,							// u_x - u_x = 0
			.5f * (jacobian[0][1] - jacobian[1][0]),		// .5(v_x - u_y)
			.5f * (jacobian[0][2] - jacobian[2][0])		);	// .5(w_x - u_z)
	// column 1
	Omega[1] = vec3(.5f * (jacobian[1][0] - jacobian[0][1]),		// .5(u_y - v_x)
			0.0f,							// v_y - v_y = 0
			.5f * (jacobian[1][2] - jacobian[2][1])		);	// .5(w_y - v_z)
	// column 2
	Omega[2] = vec3(.5f * (jacobian[2][0] - jacobian[0][2]),		// .5(u_z - w_x)
			.5f * (jacobian[2][1] - jacobian[1][2]),		// .5(v_z - w_y)
			0.0f						);	// w_z - w_z = 0

	// Omega^2 + S^2
	return Omega*Omega + S*S;
}

}  // namespace inviwo
