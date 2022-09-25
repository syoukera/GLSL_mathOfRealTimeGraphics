#version 300 es
precision highp float;
precision highp int;
out vec4 fragColor;
uniform float u_time;
uniform vec2 u_mouse;
uniform vec2 u_resolution;
ivec2 channel;

float PI = 3.1415926;

vec3 hsv2rgb(vec3 c) {
    vec3 rgb = clamp(abs(mod(c.x*6.0 + vec3(0.0, 4.0, 2.0), 6.0) - 3.0) - 1.0, 0.0, 1.0);
    return c.z * mix(vec3(1.0), rgb, c.y);
}

uvec3 k = uvec3(0x456789abu, 0x6789ab45u, 0x89ab4567u);
uvec3 u = uvec3(1, 2, 3);
const uint UINT_MAX = 0xffffffffu;
uvec2 uhash22(uvec2 n) {
    n ^= (n.yx << u.xy);
    n ^= (n.yx >> u.xy);
    n *= k.xy;
    n ^= (n.yx << u.xy);
    return n * k.xy;
}

uvec3 uhash33(uvec3 n) {
    n ^= (n.yzx << u);
    n ^= (n.yzx >> u);
    n *= k;
    n ^= (n.yzx << u);
    return n * k;
}

float hash21(vec2 p) {
    uvec2 n = floatBitsToUint(p);
    return float(uhash22(n).x) / float(UINT_MAX);
}

float hash31(vec3 p) {
    uvec3 n = floatBitsToUint(p);
    return float(uhash33(n).x) / float(UINT_MAX);
}

vec2 hash22(vec2 p) {
    uvec2 n = floatBitsToUint(p);
    return vec2(uhash22(n)) / vec2(UINT_MAX);
}

vec3 hash33(vec3 p) {
    uvec3 n = floatBitsToUint(p);
    return vec3(uhash33(n)) / vec3(UINT_MAX);
}

float vnoise21(vec2 p) {
    vec2 n = floor(p);
    float[4] v;
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 2; i++) {
            v[i+2*j] = hash21(n + vec2(i, j));
        }
    }
    vec2 f = fract(p);
    f = f*f*f*(10.0 - 15.0*f + 6.0*f*f);
    return mix(mix(v[0], v[1], f[0]), mix(v[2], v[3], f[0]), f[1]);
}

float vnoise31(vec3 p) {
    vec3 n = floor(p);
    float[8] v;
    for (int k = 0; k < 2; k++) {
        for (int j = 0; j < 2; j++) {
            for (int i = 0; i < 2; i++) {
                v[i + 2*j + 4*k] = hash31(n + vec3(i, j, k));
            }
        }
    }
    vec3 f = fract(p);
    f = f*f*f*(10.0 - 15.0*f + 6.0*f*f);
    float[2] w;
    for (int i = 0; i < 2; i++) {
        w[i] = mix(mix(v[4*i], v[4*i+1], f[0]), mix(v[4*i+2], v[4*i+3], f[0]), f[1]);
    }

    return mix(w[0], w[1], f[2]);
}


float gnoise21(vec2 p) {
    vec2 n = floor(p);
    vec2 f = fract(p);
    float[4] v;
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 2; i++) {
            vec2 g = normalize(hash22(n + vec2(i, j)) - vec2(0.5));
            v[i + 2*j] = dot(g, f - vec2(i, j));
        }
    }
    f = f*f*f*(10.0 - 15.0*f + 6.0*f*f);
    return 0.5*mix(mix(v[0], v[1], f[0]), mix(v[2], v[3], f[0]), f[1]) + 0.5;
}

float gnoise31(vec3 p) {
    vec3 n = floor(p);
    vec3 f = fract(p);
    float[8] v;
    for (int k = 0; k < 2; k++) {
        for (int j = 0; j < 2; j++) {
            for (int i = 0; i < 2; i++) {
                vec3 g = normalize(hash33(n + vec3(i, j, k)) - vec3(0.5));
                v[i+2*j+4*k] = dot(g, f - vec3(i, j, k));
            }
        }
    }
    f = f*f*f*(10.0 - 15.0*f + 6.0*f*f);
    float[2] w;
    for (int i = 0; i < 2; i++) {
        w[i] = mix(mix(v[4*i], v[4*i+1], f[0]), mix(v[4*i+2], v[4*i+3], f[0]), f[1]);
    }

    return 0.5*mix(w[0], w[1], f[2])+0.5;
}

vec2 grad(vec2 p) {
    float eps = 0.0001;
    return 0.5 * (vec2(
            vnoise21(p + vec2(eps, 0.0)) - vnoise21(p - vec2(eps, 0.0)),
            vnoise21(p + vec2(0.0, eps)) - vnoise21(p - vec2(0.0, eps))
        ))/ eps;
}

vec2 rot2(vec2 p, float t) {
    return vec2(cos(t)*p.x - sin(t) * p.y, sin(t)*p.x + cos(t) * p.y);
}

vec3 rotX(vec3 p, float t) {
    p.yz = rot2(p.yz, t);
    return p;
}

vec3 rotY(vec3 p, float t) {
    p.xz = rot2(p.xz, t);
    return p;
}

vec3 rotZ(vec3 p, float t) {
    p.xy = rot2(p.xy, t);
    return p;
}

float rotNoise21(vec2 p, float ang) {
    vec2 n = floor(p);
    vec2 f = fract(p);
    float[4] v;
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 2; i++) {
            vec2 g = normalize(hash22(n + vec2(i, j)) - vec2(0.5));
            g = rot2(g, ang);
            v[i + 2*j] = dot(g, f - vec2(i, j));
        }
    }
    f = f*f*f*(10.0 - 15.0*f + 6.0*f*f);
    return 0.5*mix(mix(v[0], v[1], f[0]), mix(v[2], v[3], f[0]), f[1]) + 0.5;
}

float gtable2(vec2 lattice, vec2 p) {
    uvec2 n = floatBitsToUint(lattice);
    uint ind = uhash22(n).x >> 29;
    float u = 0.92387953 * (ind < 4u ? p.x : p.y);
    float v = 0.38268343 * (ind < 4u ? p.y : p.x);
    return ((ind &1u) == 0u ? u: -u) + ((ind & 2u) == 0u ? v : -v);
}

float pnoise21(vec2 p) {
    vec2 n = floor(p);
    vec2 f = fract(p);
    float[4] v;
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 2; i++) {
            // vec2 g = normalize(hash22(n + vec2(i, j)) - vec2(0.5));
            v[i+2*j] = gtable2(n + vec2(i, j), f - vec2(i, j));
        }
    }
    f = f*f*f*(10.0 - 15.0*f + 6.0*f*f);
    return 0.5*mix(mix(v[0], v[1], f[0]), mix(v[2], v[3], f[0]), f[1]) + 0.5;
}

float gtable3(vec3 lattice, vec3 p) {
    uvec3 n = floatBitsToUint(lattice);
    uint ind = uhash33(n).x >> 28;
    float u = ind < 8u ? p.x : p.y;
    float v = ind < 4u ? p.y : ind == 12u || ind == 14u ? p.x : p.z;
    return ((ind & 1u) == 0u ? u: - u) + ((ind & 2u) == 0u ? v: -v);
}

float pnoise31(vec3 p) {
    vec3 n = floor(p);
    vec3 f = fract(p);
    float[8] v;
    for (int k = 0; k < 2; k++) {
        for (int j = 0; j < 2; j++) {
            for (int i = 0; i < 2; i++) {
                // vec3 g = normalize(hash33(n + vec3(i, j, k)) - vec3(0.5));
                // v[i+2*j+4*k] = dot(g, f - vec3(i, j, k));
                v[i+2*j+4*k] = gtable3(n + vec3(i, j, k), f - vec3(i, j, k))
                * 0.70710678;
            }
        }
    }
    f = f*f*f*(10.0 - 15.0*f + 6.0*f*f);
    float[2] w;
    for (int i = 0; i < 2; i++) {
        w[i] = mix(mix(v[4*i], v[4*i+1], f[0]), mix(v[4*i+2], v[4*i+3], f[0]), f[1]);
    }

    return 0.5*mix(w[0], w[1], f[2])+0.5;
}

float periodicNoise21(vec2 p, float period) {
    vec2 n = floor(p);
    vec2 f = fract(p);
    float[4] v;
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 2; i++) {
            v[i+2*j] = gtable2(mod(n + vec2(i, j), period), f - vec2(i, j));
        }
    }
    f = f*f*f*(10.0 - 15.0*f + 6.0*f*f);
    return 0.5*mix(mix(v[0], v[1], f[0]), mix(v[2], v[3], f[0]), f[1]) + 0.5;
}

float fbm21(vec2 p, float g) {
    float val = 0.0;
    float amp = 1.0;
    float freq = 1.0;
    for (int i = 0; i < 4; i++) {
        val += amp * (pnoise21(freq*p) - 0.5);
        amp *= g;
        freq *= 2.01;
    }
    return 0.5 * val + 0.5;
}

float warp21(vec2 p, float g) {
    float val = 0.0;
    for (int i = 0; i < 4; i++) {
        val = fbm21(p + g * vec2(cos(2.0 * PI * val), sin(2.0 * PI * val)), 0.5);
    }
    return val;
}

float converter(float v) {
    float time = abs(mod(0.1 * u_time, 2.0) - 1.0);
    float n = floor(8.0 * time);
    return channel == ivec2(1, 0) ? step(time, v) :
        channel == ivec2(2, 0) ? (floor(n * v) + step(0.5, fract(n * v))) / n :
        channel == ivec2(0, 1) ? smoothstep(0.5 * (1.0 - time), 0.5 * (1.0 + time), v) :
        channel == ivec2(1, 1) ? pow(v, 2.0 * time) :
        channel == ivec2(2, 1) ? 0.5 * sin(4.0 * PI * v + u_time) + 0.5 :
        v;
}

vec3 blend(float a, float b) {
    float time = abs(mod(0.1 * u_time, 2.0) - 1.0);
    vec3[2] col2 = vec3[] (
        vec3(a, a, 1), 
        vec3(0, b, b)
    );
    return channel == ivec2(1, 0) ? mix(col2[0], col2[1], time) :
        mix(col2[0], col2[1], smoothstep(0.5 - 0.5 * time, 0.5 + 0.5 * time,
        b / (a + b)));
}

float fdist(vec2 p) {
    vec2 n = floor(p + 0.5);
    float  dist = sqrt(2.0);
    for (float j = -2.0; j <= 2.0; j++) {
        for (float i = -2.0; i <= 2.0; i++) {
            vec2 glid = n + vec2(i, j);
            vec2 jitter = sin(u_time) * (hash22(glid) - 0.5);
            dist = min(dist, distance(glid + jitter, p));
        }
    }
    return dist;
}

float fdist21(vec2 p) {
    vec2 n = floor(p + 0.5);
    float  dist = sqrt(2.0);
    for (float j = -0.0; j <= 2.0; j++) {
        vec2 glid;
        glid.y = n.y + sign(mod(j, 2.0) - 0.5) * ceil(j * 0.5);
        if (abs(glid.y - p.y) - 0.5 > dist) {
            continue;
        }
        for (float i = -1.0; i <= 1.0; i++) {
            glid.x = n.x + i;
            vec2 jitter = (hash22(glid) - 0.5);
            dist = min(dist, distance(glid + jitter, p));
        }
    }
    return dist;
}

float fdist31(vec3 p){
    vec3 n = floor(p + 0.5);
    float dist = sqrt(3.0);
    for(float k = 0.0; k <= 2.0; k ++ ){
            vec3 glid;
            glid.z = n.z + sign(mod(k, 2.0) - 0.5) * ceil(k * 0.5);
            if (abs(glid.z - p.z) - 0.5 > dist){
                continue;
            }
        for(float j = 0.0; j <= 2.0; j ++ ){
            glid.y = n.y + sign(mod(j, 2.0) - 0.5) * ceil(j * 0.5);
            if (abs(glid.y - p.y) - 0.5 > dist){
                continue;
            }
            for(float i = -1.0; i <= 1.0; i ++ ){
                glid.x = n.x + i;
                vec3 jitter = hash33(glid) - 0.5;
                dist = min(dist, length(glid + jitter - p));
            }
        }
    }
    return dist;
}

vec4 sort(vec4 list, float v) {
    bvec4 res = bvec4(step(v, list));
    return res.x ? vec4(v, list.xyz):
        res.y ? vec4(list.x, v, list.yz):
        res.z ? vec4(list.xy, v, list.z):
        res.w ? vec4(list.xyz, v):
        list;
}

vec4 fdist24(vec2 p) {
    vec2 n = floor(p + 0.5);
    vec4  dist4 = vec4(length(1.5 - abs(p - n)));
    for (float j = 0.0; j <= 4.0; j++) {
        vec2 glid;
        glid.y = n.y + sign(mod(j, 2.0) - 0.5) * ceil(j * 0.5);
        // if (abs(glid.y - p.y) - 0.5 > dist4) {
        //     continue;
        // }
        for (float i = -2.0; i <= 2.0; i++) {
            glid.x = n.x + i;
            vec2 jitter = (hash22(glid) - 0.5);
            dist4 = sort(dist4, length(glid + jitter - p));
        }
    }
    return dist4;
}

// float length2(vec2 p) {
//     float t = mod(u_time, 3.0);
//     p = abs(p);
//     return t < 1.0 ? length(p) :
//     t < 2.0 ? dot(p, vec2(1.0)) :
//     max(p.x, p.y);
// }

float length2(vec2 p) {
    p = abs(p);
    float d = 4.0 * sin(0.5 * u_time) + 5.0;
    return pow(pow(p.x, d) + pow(p.y, d), 1.0 / d);
}

float circle(vec2 p, vec2 c, float r) {
    return length2(p - c);
}

float rect(vec2 p, vec2 c, vec2 d) {
    p = abs(p - c);
    return length(max(p - d, vec2(0.0))) + min(max(p.x - d.x, p.y - d.y), 0.0);
}

vec3 contour(float v, float interval) {
    return abs(v) < 0.01 ? vec3(0.0):
        mod(v, interval) < 0.01 ? vec3(1.0) :
        mix(vec3(0, 0, 1), vec3(1, 0, 0), atan(v) / PI + 0.5);
}


vec2 voronoi2(vec2 p) {
    vec2 n = floor(p + 0.5);
    float  dist = sqrt(2.0);
    vec2 id;
    for (float j = -0.0; j <= 2.0; j++) {
        vec2 glid;
        glid.y = n.y + sign(mod(j, 2.0) - 0.5) * ceil(j * 0.5);
        if (abs(glid.y - p.y) - 0.5 > dist) {
            continue;
        }
        for (float i = -1.0; i <= 1.0; i++) {
            glid.x = n.x + i;
            vec2 jitter = sin(u_time)*(hash22(glid) - 0.5);
            if (length2(glid + jitter - p) <= dist) {
                dist = length2(glid + jitter - p);
                id = glid;
            }
        }
    }
    return id;
}

float text(vec2 st) {
    float time = 0.3 * u_time;
    float v0 = warp21(st + time, 1.0);
    float v1 = fdist31(vec3(st + time, time));
    time = abs(mod(time, 2.0) - 1.0);
    return mix(v0, v1, smoothstep(0.25, 0.75, time));
}

void main() {
    vec2 p = (gl_FragCoord.xy * 2.0 - u_resolution) / min(u_resolution.x, u_resolution.y);

    vec3 cPos = vec3(0.0, 0.0, 0.0);
    float t = -0.5 * PI * (u_mouse.y / u_resolution.y);
    vec3 cDir = rotX(vec3(0.0, 0.0, -1.0), t);
    vec3 cUp  = rotX(vec3(0.0, 1.0, 0.0), t);
    vec3 cSide = cross(cDir, cUp);
    float targetDepth = 1.0;
    vec3 ray = cSide * p.x + cUp * p.y + cDir * targetDepth - cPos;
    ray = normalize(ray);
    vec3 groundNormal = vec3(0.0, 1.0, 0.0);
    float groundHeight = 1.0 + (u_mouse.x / u_resolution.x);
    vec3 lPos = vec3(0.0, 0.0, 0.0);
    if (dot(ray, groundNormal) < 0.0) {
        vec3 hit = cPos - ray * groundHeight / dot(ray, groundNormal);
        float diff = max(dot(normalize(lPos - hit), groundNormal), 0.0);
        diff *= 0.5 + u_mouse.y / u_resolution.y;
        diff = pow(diff, 0.5 + u_mouse.x / u_resolution.x);
        fragColor.rgb = vec3(diff*text(hit.zx));
    } else {
        fragColor.rgb = vec3(0.0);
    }
    fragColor.a = 1.0;
}