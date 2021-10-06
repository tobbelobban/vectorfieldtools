ivw_module(VectorFieldTools)

set(HEADER_FILES
    include/KTH/vectorfieldtools/algorithm/jacobiancompute.h
    include/KTH/vectorfieldtools/algorithm/omega2s2.h
    include/KTH/vectorfieldtools/processors/okuboweiss.h
    include/KTH/vectorfieldtools/processors/okuboweissgl.h
    include/KTH/vectorfieldtools/processors/qhunt.h
    include/KTH/vectorfieldtools/processors/qhuntgl.h
    include/KTH/vectorfieldtools/processors/vectorfield3dcurl.h
    include/KTH/vectorfieldtools/vectorfieldtoolsmodule.h
    include/KTH/vectorfieldtools/vectorfieldtoolsmoduledefine.h
)
ivw_group("Header Files" ${HEADER_FILES})

set(SOURCE_FILES
    src/algorithm/jacobiancompute.cpp
    src/algorithm/omega2s2.cpp
    src/processors/okuboweiss.cpp
    src/processors/okuboweissgl.cpp
    src/processors/qhunt.cpp
    src/processors/qhuntgl.cpp
    src/processors/vectorfield3dcurl.cpp
    src/vectorfieldtoolsmodule.cpp
)
ivw_group("Source Files" ${SOURCE_FILES})

set(SHADER_FILES
	# Add shaders
)
ivw_group("Shader Files" ${SHADER_FILES})

set(TEST_FILES
    tests/unittests/vectorfieldtools-unittest-main.cpp
)
ivw_add_unittest(${TEST_FILES})

ivw_create_module(${SOURCE_FILES} ${HEADER_FILES} ${SHADER_FILES})

# Add shader directory to install package
#ivw_add_to_module_pack(${CMAKE_CURRENT_SOURCE_DIR}/glsl)
