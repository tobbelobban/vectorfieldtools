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

#include <KTH/vectorfieldtools/processors/okuboweiss.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo OkuboWeiss::processorInfo_{
    "org.inviwo.OkuboWeiss",      // Class identifier
    "Okubo Weiss",                // Display name
    "Vector Field Visualization",              // Category
    CodeState::Experimental,  // Code state
    Tags::None,               // Tags
};
const ProcessorInfo OkuboWeiss::getProcessorInfo() const { return processorInfo_; }

OkuboWeiss::OkuboWeiss()
    : Processor()
    , vol_inport_("Vector_field_volume_inport")
	, vol_outport_("Scalar_field_volume_outport") {
    addPort(vol_inport_);
	addPort(vol_outport_);
}

void OkuboWeiss::process() {
    const std::shared_ptr<const Volume> vector_field = vol_inport_.getData();
	const size3_t dims = vector_field->getDimensions();
	// make dest volume
	auto okubo_weiss_vol_repr = std::make_shared<VolumeRAMPrecision<float>>(dims);
    float* okubo_weiss_raw_ptr = okubo_weiss_vol_repr->getDataTyped();
	// iterate over vector field and compute OW
	size_t dst_index = 0;
	std::vector<float> j;
	float max_val = std::numeric_limits<float>::min();
	float min_val = std::numeric_limits<float>::max();
	for(size_t iz = 0; iz < dims.z; ++iz) {
		for(size_t iy = 0; iy < dims.y; ++iy) {
			for(size_t ix = 0; ix < dims.x; ++ix) {
				// compute jacobian at (ix, iy, iz)
				j = jacobian_computer.get(vector_field, size3_t(ix,iy,iz));
				// compute & store OW
				// -2(u_y*v_x + u_z*w_x + w_y*v_z) - u_x^2 - v_y^2 - w_z^2
				float res = -2.0f * (j[1] * j[3] + j[2] * j[6] + j[7] * j[5]) - j[0] * j[0] - j[4] * j[4] - j[8] * j[8];
				okubo_weiss_raw_ptr[dst_index++] = res;
				if(res > max_val) max_val = res;
				if(res < min_val) min_val = res;
			}
		}	
	}
	std::shared_ptr<Volume> OW = std::make_shared<Volume>(okubo_weiss_vol_repr);
    OW->setBasis(vector_field->getBasis());
	OW->setOffset(vector_field->getOffset());
	OW->copyMetaDataFrom(*vector_field);	
	OW->dataMap_.valueRange = vec2(min_val, max_val);
	OW->dataMap_.dataRange = vec2(min_val, max_val);
	OW->setModelMatrix(vector_field->getModelMatrix());
	OW->setWorldMatrix(vector_field->getWorldMatrix());
    vol_outport_.setData(OW);
}

}  // namespace inviwo
