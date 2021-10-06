
#include "utils/sampler3d.glsl"
#include "utils/gradients.glsl"

uniform sampler3D volume;
uniform VolumeParameters volumeParameters;

in vec4 texCoord_;

void main() {
    vec3 ox = vec3(volumeParameters.reciprocalDimensions.x, 0, 0);
    vec3 oy = vec3(0, volumeParameters.reciprocalDimensions.y, 0);
    vec3 oz = vec3(0, 0, volumeParameters.reciprocalDimensions.z);

    vec3 Fx = (getVoxel(volume, volumeParameters, texCoord_.xyz + ox).xyz -
               getVoxel(volume, volumeParameters, texCoord_.xyz - ox).xyz) /
              (2.0f*volumeParameters.worldSpaceGradientSpacing.x);
    vec3 Fy = (getVoxel(volume, volumeParameters, texCoord_.xyz + oy).xyz -
               getVoxel(volume, volumeParameters, texCoord_.xyz - oy).xyz) /
              (2.0f*volumeParameters.worldSpaceGradientSpacing.y);
    vec3 Fz = (getVoxel(volume, volumeParameters, texCoord_.xyz + oz).xyz -
               getVoxel(volume, volumeParameters, texCoord_.xyz - oz).xyz) /
              (2.0f*volumeParameters.worldSpaceGradientSpacing.z);
	
	mat3 S;
	S[0] = vec3(Fx.x, .5f*(Fx.y+Fy.x), .5f*(Fx.z+Fz.x));
	// column 1
	S[1] = vec3(.5f*(Fx.y+Fy.x), Fy.y, .5f*(Fy.z+Fz.y));
	// column 2
	S[2] = vec3(.5f*(Fx.z+Fz.x), .5f*(Fy.z+Fz.y), Fz.z);
	
	mat3 Omega;
	// column 0
	Omega[0] = vec3(Fx.x, .5f*(Fx.y-Fy.x), .5f*(Fx.z-Fz.x));
	// column 1
	Omega[1] = vec3(.5f*(Fy.x-Fx.y), Fy.y, .5f*(Fy.z-Fz.y));
	// column 2
	Omega[2] = vec3(.5f*(Fz.x-Fx.z), .5f*(Fz.y-Fy.z), Fz.z);

	mat3 Omega2S2 = Omega*Omega + S*S;

	FragData0 = vec4(-Omega2S2[0][0] - Omega2S2[1][1] - Omega2S2[2][2]);
}