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

#include <KTH/vectorfieldtools/processors/qhunt.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo QHunt::processorInfo_{
    "org.inviwo.QHunt",      // Class identifier
    "QHunt",                // Display name
    "Undefined",              // Category
    CodeState::Experimental,  // Code state
    Tags::None,               // Tags
};
const ProcessorInfo QHunt::getProcessorInfo() const { return processorInfo_; }

QHunt::QHunt()
    : Processor()
    , volume_in_("Vector_field_volume_inport")
	, volume_out_("Scalar_field_volume_outport") {

    addPort(volume_in_);
	addPort(volume_out_);
}

void QHunt::process() {
	const std::shared_ptr<const Volume> vector_field = volume_in_.getData();
	const size3_t dims = vector_field->getDimensions();
	// make dest volume
	auto QHunt_vol_repr = std::make_shared<VolumeRAMPrecision<float>>(dims);
    float* QHunt_raw_ptr = QHunt_vol_repr->getDataTyped();
	// iterate over vector field and compute OW
	size_t dst_index = 0;
	float max_val = std::numeric_limits<float>::min();
	float min_val = std::numeric_limits<float>::max();
	for(size_t iz = 0; iz < dims.z; ++iz) {
		for(size_t iy = 0; iy < dims.y; ++iy) {
			for(size_t ix = 0; ix < dims.x; ++ix) {
				// compute Omega^2 + S^2
				const mat3 Omega2S2 = omega2s2_.get(vector_field, size3_t(ix,iy,iz));
				const float res = -Omega2S2[0][0] - Omega2S2[1][1] - Omega2S2[2][2];
				QHunt_raw_ptr[dst_index++] = res;
				if(res > max_val) max_val = res;
				if(res < min_val) min_val = res;
			}
		}	
	}
	std::shared_ptr<Volume> QHunt = std::make_shared<Volume>(QHunt_vol_repr);
    QHunt->setBasis(vector_field->getBasis());
	QHunt->setOffset(vector_field->getOffset());
	QHunt->copyMetaDataFrom(*vector_field);	
	QHunt->dataMap_.valueRange = vec2(min_val, max_val);
	QHunt->dataMap_.dataRange = vec2(min_val, max_val);
	QHunt->setModelMatrix(vector_field->getModelMatrix());
	QHunt->setWorldMatrix(vector_field->getWorldMatrix());
    volume_out_.setData(QHunt);
}

}  // namespace inviwo
