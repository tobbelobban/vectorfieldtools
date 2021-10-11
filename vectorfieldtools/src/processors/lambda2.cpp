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

#include <KTH/vectorfieldtools/processors/lambda2.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo Lambda2::processorInfo_{
    "org.inviwo.Lambda2",      // Class identifier
    "Lambda2",                // Display name
    "Undefined",              // Category
    CodeState::Experimental,  // Code state
    Tags::None,               // Tags
};
const ProcessorInfo Lambda2::getProcessorInfo() const { return processorInfo_; }

Lambda2::Lambda2()
    : Processor()
    , volume_in_("Vector_field_volume_inport")
	, volume_out_("Scalar_field_volume_outport") {

    addPort(volume_in_);
	addPort(volume_out_);
}

Eigen::Matrix3f Lambda2::glmToEigenMat3FLOAT(const glm::mat3 glm_mat3) {
	Eigen::Matrix3f Eigen_mat3;
	for(size_t col = 0; col < 3; ++col) {
		for(size_t row = 0; row < 3; ++row) {
			Eigen_mat3(row,col) = glm_mat3[col][row];
		}
	}
	return Eigen_mat3;
}

void Lambda2::process() {
    const std::shared_ptr<const Volume> vector_field = volume_in_.getData();
	const size3_t dims = vector_field->getDimensions();
	// make dest volume
	auto L2_vol_repr = std::make_shared<VolumeRAMPrecision<float>>(dims);
    float* L2_raw_ptr = L2_vol_repr->getDataTyped();
	// iterate over vector field and compute L2
	float max_val = std::numeric_limits<float>::min();
	float min_val = std::numeric_limits<float>::max();
	size_t dst_index = 0;
	for(size_t iz = 0; iz < dims.z; ++iz) {
		for(size_t iy = 0; iy < dims.y; ++iy) {
			for(size_t ix = 0; ix < dims.x; ++ix) {				
				const auto eigenvalues = eigen_3f_solver.compute(glmToEigenMat3FLOAT(omega2s2.get(vector_field, size3_t(ix,iy,iz))),false).eigenvalues();
				float reals[3] = {eigenvalues.col(0)[0].real(), eigenvalues.col(0)[1].real(), eigenvalues.col(0)[2].real()};
				std::sort(reals, reals+3);
				float res = reals[1];
				L2_raw_ptr[dst_index++] = res;
				if(res > max_val) max_val = res;
				if(res < min_val) min_val = res;
			}
		}	
	}
	std::shared_ptr<Volume> Lambda2 = std::make_shared<Volume>(L2_vol_repr);
    Lambda2->setBasis(vector_field->getBasis());
	Lambda2->setOffset(vector_field->getOffset());
	Lambda2->copyMetaDataFrom(*vector_field);	
	Lambda2->dataMap_.valueRange = vec2(min_val, max_val);
	Lambda2->dataMap_.dataRange = vec2(min_val, max_val);
	Lambda2->setModelMatrix(vector_field->getModelMatrix());
	Lambda2->setWorldMatrix(vector_field->getWorldMatrix());
    volume_out_.setData(Lambda2);
}

}  // namespace inviwo
