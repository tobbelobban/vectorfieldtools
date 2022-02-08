/*********************************************************************************
 *
 * Inviwo - Interactive Visualization Workshop
 *
 * Copyright (c) 2015-2021 Inviwo Foundation
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

#include <KTH/vectorfieldtools/processors/qhuntgl.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo QHuntGL::processorInfo_{
	"org.inviwo.QHuntGL",  // Class identifier
	"QHunt GL",           // Display name
	"Vector Field Visualization",     // Category
	CodeState::Stable,                // Code state
	Tags::GL,                         // Tags
};
const ProcessorInfo QHuntGL::getProcessorInfo() const { return processorInfo_; }

QHuntGL::QHuntGL() : VolumeGLProcessor("qhunt.frag") {
    this->dataFormat_ = DataFloat32::get();
}

QHuntGL::~QHuntGL() {}

void QHuntGL::postProcess() {
	// determine min & max values for output volume->dataMap_
	const VolumeRAM* vol_ram = volume_->getRepresentation<VolumeRAM>();
	const size3_t dims = vol_ram->getDimensions();
	const float* vol_ram_ptr = static_cast<const float*>(vol_ram->getData());
	// find min and max
	float min = std::numeric_limits<float>::max();
	float max = std::numeric_limits<float>::min();	
	const size_t flattened_vol_size = dims.x * dims.y * dims.z;
	float val;
	for(size_t index = 0; index < flattened_vol_size; ++index) {
		val = vol_ram_ptr[index];
		if(val > max) max = val;
		if(val < min) min = val;
	}
	volume_->dataMap_.dataRange = volume_->dataMap_.valueRange = vec2(min, max);
	
}

}  // namespace inviwo
