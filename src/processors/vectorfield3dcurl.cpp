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
    "Vector Field3DCurl",                // Display name
    "Undefined",              // Category
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
	const int xy_ = vol_dims.x * vol_dims.y;
	const vec3 spacing = vector_field->getWorldSpaceGradientSpacing();

	auto curl_ = vector_field->getRepresentation<VolumeRAM>()->dispatch<std::shared_ptr<Volume>, dispatching::filter::Float3s>(
		[&](auto vector_field_pr) {
			using VolVecType = util::PrecisionValueType<decltype(vector_field_pr)>;
			const VolVecType* vector_field_data = vector_field_pr->getDataTyped();
            auto dstRam = std::make_shared<VolumeRAMPrecision<float>>(vol_dims);
            float* dstData = dstRam->getDataTyped();
            // iterate over each vector
            for(int i = 0; i < vol_dims.z; ++i) {
				int z_offset = i * xy_;
				for(int j = 0; j < vol_dims.y; ++j) {
					int yz_offset = z_offset + vol_dims.y * j;
					for(int k = 0; k < vol_dims.z; ++k) {
                        // if on domain boundary - skip
                        if( i == 0 || i == vol_dims.x-1 ||
                            j == 0 || j == vol_dims.y-1 ||
                            k == 0 || k == vol_dims.z-1) {
							dstData[yz_offset + k] = 0.0f;
							continue;
						}
                        // otherwise, use central differences to compute Jacobian
						const vec3 u = vec3(vector_field_data[yz_offset + k + 1] - vector_field_data[yz_offset + k - 1]) / (2.0f * spacing.x);
						const vec3 v = vec3(vector_field_data[yz_offset + k + vol_dims.y] - vector_field_data[yz_offset + k - vol_dims.y]) / (2.0f * spacing.y);
						const vec3 w = vec3(vector_field_data[yz_offset + k + xy_] - vector_field_data[yz_offset + k - xy_]) / (2.0f * spacing.z);
                        // store curl magnitude in output volume
                        dstData[yz_offset + k] = length(vec3(v.z - w.y, w.x - u.z, u.y - v.x));
					}
				}
			}
            return std::make_shared<Volume>(dstRam);
        }
    );
    curl_.setBasis(vector_field->getBasis());
    vol_outport_.setData(curl_);
}

}  // namespace inviwo
