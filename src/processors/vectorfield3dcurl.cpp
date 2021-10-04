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
	const size_t xy_ = vol_dims.x * vol_dims.y;
	const vec3 spacing = vector_field->getWorldSpaceGradientSpacing();
	float max_curl_mag = 0.0f;
	JacobianCompute jacobian_computer;
	std::vector<float> j(9, 1.0f);
	std::shared_ptr<Volume> curl_;
    auto dst_ram = std::make_shared<VolumeRAMPrecision<vec3>>(vol_dims);
    vec3* dst_data = dst_ram->getDataTyped();
    // iterate over each vector
    for(int iz = 0; iz < vol_dims.z; ++iz) {
		size_t z_offset = iz * xy_;
		for(int iy = 0; iy < vol_dims.y; ++iy) {
			size_t yz_offset = z_offset + vol_dims.x * iy;
			for(int ix = 0; ix < vol_dims.x; ++ix) {
				j = jacobian_computer.get(vector_field, size3_t(ix,iy,iz));
				 //store curl magnitude in output volume
				dst_data[yz_offset + ix] = vec3(	j[5] - j[7],	// v.z - w.y
													j[6] - j[2],	// w.x - u.z 
													j[1] - j[3] );	// u.y - v.x	
				//if(max_curl_mag < dst_data[yz_offset + ix]) max_curl_mag = dst_data[yz_offset + ix];
			}
		}
	}
	//LogProcessorInfo(max_curl_mag);
	std::shared_ptr<Volume> curl = std::make_shared<Volume>(dst_ram);
    curl->setBasis(vector_field->getBasis());
	curl->setOffset(vector_field->getOffset());
	curl->copyMetaDataFrom(*vector_field);
	//curl->dataMap_.valueRange = vec2(0, max_curl_mag);
	//curl->dataMap_.dataRange = vec2(0, max_curl_mag);
	curl->setModelMatrix(vector_field->getModelMatrix());
	curl->setWorldMatrix(vector_field->getWorldMatrix());
    vol_outport_.setData(curl);
}

}  // namespace inviwo
