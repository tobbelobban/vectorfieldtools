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
#include <KTH/vectorfieldtools/vectorfieldtoolsmodule.h>

#include <modules/opengl/shader/shadermanager.h>
#include <modules/opengl/openglmodule.h>
#include <modules/opengl/openglcapabilities.h>

#include <KTH/vectorfieldtools/processors/okuboweiss.h>
#include <KTH/vectorfieldtools/processors/okuboweissgl.h>
#include <KTH/vectorfieldtools/processors/qhunt.h>
#include <KTH/vectorfieldtools/processors/qhuntgl.h>
#include <KTH/vectorfieldtools/processors/vectorfield3dcurl.h>


namespace inviwo {

VectorFieldToolsModule::VectorFieldToolsModule(InviwoApplication* app) : InviwoModule(app, "VectorFieldTools") {
    // Add a directory to the search path of the Shadermanager
    //ShaderManager::getPtr()->addShaderSearchPath(getPath(ModulePath::GLSL));
    // Register objects that can be shared with the rest of inviwo here:
	ShaderManager::getPtr()->addShaderSearchPath(getPath(ModulePath::GLSL));
	//OpenCL::getPtr()->addCommonIncludeDirectory(getPath(ModulePath::CL));
    // Processors
    registerProcessor<Lambda2>();
    registerProcessor<OkuboWeiss>();
    registerProcessor<OkuboWeissGL>();
    registerProcessor<QHunt>();
    registerProcessor<QHuntGL>();
    registerProcessor<VectorField3DCurl>();    

    // Properties
    // registerProperty<VectorFieldToolsProperty>();

    // Readers and writes
    // registerDataReader(std::make_unique<VectorFieldToolsReader>());
    // registerDataWriter(std::make_unique<VectorFieldToolsWriter>());

    // Data converters
    // registerRepresentationConverter(std::make_unique<VectorFieldToolsDisk2RAMConverter>());

    // Ports
    // registerPort<VectorFieldToolsOutport>();
    // registerPort<VectorFieldToolsInport>();

    // PropertyWidgets
    // registerPropertyWidget<VectorFieldToolsPropertyWidget, VectorFieldToolsProperty>("Default");

    // Dialogs
    // registerDialog<VectorFieldToolsDialog>(VectorFieldToolsOutport);

    // Other things
    // registerCapabilities(std::make_unique<VectorFieldToolsCapabilities>());
    // registerSettings(std::make_unique<VectorFieldToolsSettings>());
    // registerMetaData(std::make_unique<VectorFieldToolsMetaData>());
    // registerPortInspector("VectorFieldToolsOutport", "path/workspace.inv");
    // registerProcessorWidget(std::string processorClassName, std::unique_ptr<ProcessorWidget> processorWidget); 
    // registerDrawer(util::make_unique_ptr<VectorFieldToolsDrawer>());
}

}  // namespace inviwo
