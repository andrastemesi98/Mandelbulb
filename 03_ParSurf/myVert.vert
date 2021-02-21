#version 130

in vec2 vs_in_pos;

out vec2 vs_out_pos;


void main()
{
	gl_Position = vec4( vs_in_pos, 0, 1);
	vs_out_pos = vs_in_pos;
}