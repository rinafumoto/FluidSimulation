#version 410 core

layout(location=0) in vec3 inPos;
layout(location=1) in vec3 inDir;
uniform mat4 MVP;
out vec4 dir;

void main()
{
  gl_Position=MVP*vec4(inPos,1);  
  dir=MVP*vec4(inDir.x,inDir.y,0,0);
}
