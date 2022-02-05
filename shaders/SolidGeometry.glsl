#version 410 core
// based on https://stackoverflow.com/questions/54686818/glsl-geometry-shader-to-replace-gllinewidth

layout (points) in;
layout (triangle_strip, max_vertices = 4) out;  
uniform vec2 size;
out vec3 colour;

void main()
{
    vec4 p1 = gl_in[0].gl_Position;
    float width = 1.0/size.x;
    float height = 1.0/size.y;

    gl_Position = p1 + vec4(-width,-height, 0.0, 0.0);
    EmitVertex();
    gl_Position = p1 + vec4(-width,height, 0.0, 0.0);
    EmitVertex();
    gl_Position = p1 + vec4(width,-height, 0.0, 0.0);
    EmitVertex();
    gl_Position = p1 + vec4(width,height, 0.0, 0.0);
    EmitVertex();

    EndPrimitive();
}



