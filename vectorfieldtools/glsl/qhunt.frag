
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
	
	// S = .5(J + transpose(J))
	mat3 S;
	
	// column 0
	S[0] = vec3(	Fx.x,						// u_x
					.5f * (Fx.y + Fy.x),		// .5(v_x + u_y)
					.5f * (Fx.z + Fz.x)		);	// .5(w_x + u_z)
	// column 1
	S[1] = vec3(	.5f * (Fy.x + Fx.y),		// .5(u_y + v_x)
					Fy.y,						// v_y
					.5f * (Fy.z + Fz.y)		);	// .5(w_y + v_z)
	// column 2
	S[2] = vec3(	.5f * (Fz.x + Fx.z),		// .5(u_z + w_x)
					.5f * (Fz.y + Fy.z),		// .5(v_z + w_y)
					Fz.z					);	// w_z
	
	// Omega = .5(J - transpose(J))
	mat3 Omega;		
	
	// column 0
	Omega[0] = vec3(	0.0f,						// u_x - u_x = 0
						.5f * (Fx.y - Fy.x),		// .5(v_x - u_y)
						.5f * (Fx.z - Fz.x)		);	// .5(w_x - u_z)
	// column 1
	Omega[1] = vec3(	.5f * (Fy.x - Fx.y),		// .5(u_y - v_x)
						0.0f,						// v_y - v_y = 0
						.5f * (Fy.z - Fz.y)		);	// .5(w_y - v_z)
	// column 2
	Omega[2] = vec3(	.5f * (Fz.x - Fx.z),		// .5(u_z - w_x)
						.5f * (Fz.y - Fy.z),		// .5(v_z - w_y)
						0.0f					);	// w_z - w_z = 0
	

	mat3 Omega2S2 = Omega*Omega + S*S;

	FragData0 = vec4(-Omega2S2[0][0] - Omega2S2[1][1] - Omega2S2[2][2]);
}