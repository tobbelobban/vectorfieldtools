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

#include <KTH/vectorfieldtools/processors/vectorfield3dcurl.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo VectorField3DCurl::processorInfo_{
    "org.inviwo.VectorField3DCurl",      // Class identifier
    "Vector Field 3D Curl",                // Display name
    "Vector Field Visualization",              // Category
    CodeState::Experimental,  // Code state
    Tags::None,               // Tags
};
const ProcessorInfo VectorField3DCurl::getProcessorInfo() const { return processorInfo_; }

VectorField3DCurl::VectorField3DCurl()
    : Processor()
	, vol_inport_("volume_inport")
    , vol_outport_("volume_outport") {

    addPort(vol_inport_);
	addPort(vol_outport_);
}

void VectorField3DCurl::process() {
    auto vector_field = vol_inport_.getData();
    const size3_t vol_dims = vector_field->getDimensions();
    auto dst_ram = std::make_shared<VolumeRAMPrecision<vec3>>(vol_dims);
    vec3* dst_data = dst_ram->getDataTyped();
    // iterate over each vector
	size_t dst_index = 0;
    for(int iz = 0; iz < vol_dims.z; ++iz) {
		for(int iy = 0; iy < vol_dims.y; ++iy) {
			for(int ix = 0; ix < vol_dims.x; ++ix) {
				// jacobian in column-major order
				mat3 j = jacobian_computer.get(vector_field, size3_t(ix,iy,iz));
				 //compute & store curl 
				dst_data[dst_index++] = vec3(	j[1][2] - j[2][1],	
												j[2][0] - j[0][2],		
												j[0][1] - j[1][0]	);	
			}
		}
	}
	std::shared_ptr<Volume> curl = std::make_shared<Volume>(dst_ram);
    curl->setBasis(vector_field->getBasis());
	curl->setOffset(vector_field->getOffset());
	curl->copyMetaDataFrom(*vector_field);	
	curl->setModelMatrix(vector_field->getModelMatrix());
	curl->setWorldMatrix(vector_field->getWorldMatrix());
    vol_outport_.setData(curl);
}

}  // namespace inviwo
