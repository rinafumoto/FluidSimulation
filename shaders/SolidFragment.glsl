#version 410 core

layout(location=0) out vec4 fragColour;
in vec3 colour;

void main()
{
  fragColour.rgb=vec3(0.7,0.7,0.7);
}