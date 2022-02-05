#version 410 core

layout (location = 0) in vec3 inPos;
out vec3 colour;
uniform mat4 MVP;

void main()
{
    gl_Position = MVP*vec4(inPos.xy, 0, 1);
    colour = vec3(255-inPos.z*100,255,255-inPos.z*100);
}