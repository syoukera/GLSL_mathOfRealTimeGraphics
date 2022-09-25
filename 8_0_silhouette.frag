#version 300 es
precision highp float;
precision highp int;
out vec4 fragColor;
uniform float u_time;
uniform vec2 u_mouse;
uniform vec2 u_resolution;
ivec2 channel;

float PI = 3.1415926;

float text(vec2 st) {
    return mod(floor(st.s) + floor(st.t), 2.0);
}

void main() {
    vec3 cPos = vec3(0.0, 0.0, 0.0);
    vec3 cDir = vec3(0.0, 0.0, -1.0);
    vec3 cUp  = vec3(0.0, 1.0, 0.0);

    vec2 p = (gl_FragCoord.xy * 2.0 - u_resolution) / min(u_resolution.x, u_resolution.y);

    vec3 cSide = cross(cDir, cUp);
    float targetDepth = 1.0;
    vec3 ray = cSide * p.x + cUp * p.y + cDir * targetDepth;

    vec3 groundNormal = vec3(0.0, 1.0, 0.0);
    if (dot(ray, groundNormal) < 0.0) {
        fragColor.rgb = vec3(1.0);
    } else {
        fragColor.rgb = vec3(0.0);
    }
    fragColor.a = 1.0;
}