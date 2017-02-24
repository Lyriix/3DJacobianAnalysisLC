#version 120

varying vec4 position_3d_original;
varying vec4 position_3d_modelview;
varying vec3 normal;
varying vec4 color;

uniform sampler2D texture;

uniform vec3 light=vec3(0.5f,0.3f,5.0f);


void main (void)
{
    vec3 n=normalize(normal);

    vec3 p=position_3d_modelview.xyz;
    vec3 vertex_to_light=normalize(light-p);
    vec3 reflected_light=reflect(-vertex_to_light,n);
    vec3 user_to_vertex=normalize(-p);

    float diffuse_term=0.8f*clamp(abs(dot(n,vertex_to_light)),0.0f,1.0f);
    float specular_term=0.2f*pow(clamp(dot(reflected_light,user_to_vertex),0.0f,1.0f),128.0f);
    float ambiant_term=0.4f;

    vec4 white=vec4(1.0f,1.0f,1.0f,0.0f);
    vec2 tex_coord=gl_TexCoord[0].xy;
    vec4 color_texture=texture2D(texture,tex_coord);
    vec4 color_final=color*color_texture;

    gl_FragColor = (ambiant_term+diffuse_term)*color_final+specular_term*white;

}
