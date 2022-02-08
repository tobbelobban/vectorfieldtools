
#include "utils/sampler3d.glsl"
#include "utils/gradients.glsl"

uniform sampler3D volume;
uniform VolumeParameters volumeParameters;

in vec4 texCoord_;

void main() {
	vec3 ox = vec3(volumeParameters.reciprocalDimensions.x, 0, 0);
	vec3 oy = vec3(0, volumeParameters.reciprocalDimensions.y, 0);
	vec3 oz = vec3(0, 0, volumeParameters.reciprocalDimensions.z);

	vec3 Fx = 	(getVoxel(volume, volumeParameters, texCoord_.xyz + ox).xyz -
	       		getVoxel(volume, volumeParameters, texCoord_.xyz - ox).xyz) /
	      		(2.0f*volumeParameters.worldSpaceGradientSpacing.x);
	
	vec3 Fy = 	(getVoxel(volume, volumeParameters, texCoord_.xyz + oy).xyz -
	       		getVoxel(volume, volumeParameters, texCoord_.xyz - oy).xyz) /
	      		(2.0f*volumeParameters.worldSpaceGradientSpacing.y);
	
	vec3 Fz = 	(getVoxel(volume, volumeParameters, texCoord_.xyz + oz).xyz -
	       		getVoxel(volume, volumeParameters, texCoord_.xyz - oz).xyz) /
	      		(2.0f*volumeParameters.worldSpaceGradientSpacing.z);

	mat3 j;
	j[0] = Fx;
	j[1] = Fy;
	j[2] = Fz;

	FragData0 = j;
}
