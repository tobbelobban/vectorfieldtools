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

#include <KTH/vectorfieldtools/processors/volumetesttool.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo VolumeTestTool::processorInfo_{
    "org.inviwo.VolumeTestTool",      // Class identifier
    "Volume Test Tool",                // Display name
    "Undefined",              // Category
    CodeState::Experimental,  // Code state
    Tags::None,               // Tags
};
const ProcessorInfo VolumeTestTool::getProcessorInfo() const { return processorInfo_; }

VolumeTestTool::VolumeTestTool()
    : Processor()
    , volume_inport_("volume_inport")
	, volume_outport_("volume_outport") {

    addPort(volume_inport_);
    addPort(volume_outport_);
}

void VolumeTestTool::process() {
    const auto vol_dims = volume_inport_.getData()->getDimensions();
	const auto vol_basis = volume_inport_.getData()->getBasis();
	const auto vol_offset = volume_inport_.getData()->getOffset();
	const auto vol_spacing = volume_inport_.getData()->getWorldSpaceGradientSpacing();
	
	const auto idx_mat = volume_inport_.getData()->getIndexMatrix();
	const auto wld_mat = volume_inport_.getData()->getWorldMatrix();
	const auto mdl_mat = volume_inport_.getData()->getModelMatrix();
	
	LogProcessorInfo("dimsensions:");
	LogProcessorInfo(vol_dims);

	LogProcessorInfo("basis:");
	LogProcessorInfo(vol_basis);
	
	LogProcessorInfo("offset:");
	LogProcessorInfo(vol_offset);

	LogProcessorInfo("spacing:");
	LogProcessorInfo(vol_spacing);

	LogProcessorInfo("index matrix:");
	LogProcessorInfo(idx_mat);

	LogProcessorInfo("world matrix:");
	LogProcessorInfo(wld_mat);

	LogProcessorInfo("model matrix:");
	LogProcessorInfo(mdl_mat);


	volume_outport_.setData(volume_inport_.getData());
}

}  // namespace inviwo
