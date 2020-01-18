/* cd_vector.h - Chemical Development C-style vector library

Everything is contained in this header file, as static inline

It requires `math.h` and `-lm` for the sqrt/etc functions, but that is the only thing

Types:
`vN`: vector of floats, length `N`
`mMxN`: matrix of float, size `M` rows by `N` cols


VERSION 0.1
AUTHOR: Cade Brown <brown.cade@gmail.com>

*/


#pragma once
#ifndef CD_VECTOR_H__
#define CD_VECTOR_H__

// for sqrt*(), hypot*()
#include <math.h>

/* structures */

/* v2 - <f, f>, 2d float */
typedef union {

    // (_), raw data
    float _[2];

    // (x, y) cartesian point
    struct { float x, y; };

    // (u, v) UV/texcoord point
    struct { float u, v; };

} v2;

/* v3- <f, f, f>, 3d float */
typedef union {

    // (_), raw data
    float _[3];

    // (x, y, z) cartesian point
    struct { float x, y, z; };

    // (u, v, w) UV/texcoord point (also baryocentric)
    struct { float u, v, w; };

    // (xy, _z) xy is a swizzle
    struct { v2 xy; float _z; };

    // (_x, yz) ys is a swizzle
    struct { float _x; v2 yz; };

} v3;



/* v4- <f, f, f, f>, 4d float */
typedef union {

    // (_), raw data
    float _[4];

    // (x, y, z, w) cartesian point
    struct { float x, y, z, w; };

    // (xy, zw) xy and zw are swizzles
    struct { v2 xy; v2 zw; };

    // (_x, yz, _w) ys is a swizzle
    struct { float _x; v2 yz; float _w; };

    // (xyz, __w) xyz is a swizzle
    struct { v3 xyz; float __w; };

} v4;


/* m2x2- 2x2d float matrix
    [a, b]
    [c, d]
*/
typedef union {

    // (_), raw data in row major order
    float _[4];

    // [a, b], [c, d] matrix values
    struct { float a, b, c, d; };

    // [row0] [row1] matrix rows
    struct { v2 rows[2]; };

    // [X] [Y] rows corresponding to direction vectors
    struct { v2 X, Y; };

    // row column values
    struct {
        float r0c0, r0c1;
        float r1c0, r1c1;
    };

    // column row values
    struct {
        float c0r0, c1r0;
        float c0r1, c1r1;
    };

} m2x2;


/* m3x3 - 3x3d float matrix

    [a, b, c]
    [d, e, f]
    [g, h, i]

*/
typedef union {

    // (_), raw data in row major order
    float _[9];

    // [a, b, c], [d, e, f], [g, h, i] matrix values
    struct { float a, b, c, d, e, f, g, h, i; };

    // [row0] [row1] [row2] matrix rows
    struct { v3 rows[3]; };

    // [X] [Y] [Z] rows corresponding to direction vectors
    struct { v3 X, Y, Z; };

    // row column values
    struct {
        float r0c0, r0c1, r0c2;
        float r1c0, r1c1, r1c2;
        float r2c0, r2c1, r2c2;
    };

    // column row values
    struct {
        float c0r0, c1r0, c2r0;
        float c0r1, c1r1, c2r1;
        float c0r2, c1r2, c2r2;
    };

} m3x3;



/* m4x4 - 4x4d float matrix

    [a, b, c, d]
    [e, f, g, h]
    [i, j, k, l]
    [m, n, o, p]

*/
typedef union {

    // (_), raw data in row major order
    float _[16];

    // [a, b, c, d], [e, f, g, h], [i, j, k, l], [m, n, o, p] matrix values
    struct { float a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p; };

    // [row0] [row1] matrix rows
    struct { v4 rows[4]; };

    // [X] [Y] [Z] [W] rows corresponding to direction vectors
    struct { v4 X, Y, Z, W; };

    // row/column values
    struct {
        float r0c0, r0c1, r0c2, r0c3;
        float r1c0, r1c1, r1c2, r1c3;
        float r2c0, r2c1, r2c2, r2c3;
        float r3c0, r3c1, r3c2, r3c3;
    };

    // column/row values
    struct {
        float c0r0, c1r0, c2r0, c3r0;
        float c0r1, c1r1, c2r1, c3r1;
        float c0r2, c1r2, c2r2, c3r2;
        float c0r3, c1r3, c2r3, c3r3;
    };


} m4x4;


// aliases for square matrices
/*typedef m2x2 m2;
typedef m3x3 m3;
typedef m4x4 m4;*/


/* constructors */

// constructs a `v2` given x and y components
#define V2(_x, _y) ((v2){ (float)(_x), (float)(_y) })
// constructs a `v3` given x, y, and z components
#define V3(_x, _y, _z) ((v3){ (float)(_x), (float)(_y), (float)(_z) })
// constructs a `v4` given x, y, z, and w components
#define V4(_x, _y, _z, _w) ((v4){ (float)(_x), (float)(_y), (float)(_z), (float)(_w) })

// constructs a `m2x2` given values a, b, c, d such that M=
// [a b]
// [c d]
#define M2X2(_a, _b, _c, _d) ((m2x2){ (float)(_a), (float)(_b), (float)(_c), (float)(_d) })

// constructs a `m3x3` with dense values in row major order
#define M3X3(_a, _b, _c, _d, _e, _f, _g, _h, _i) ((m3x3){ \
    (float)(_a), (float)(_b), (float)(_c), \
    (float)(_d), (float)(_e), (float)(_f), \
    (float)(_g), (float)(_h), (float)(_i) \
})

// constructs a `m4x4` with dense values in row major order
#define M4X4(_a, _b, _c, _d,  _e, _f, _g, _h,  _i, _j, _k, _l,  _m, _n, _o, _p) ((m4x4){ \
    (float)(_a), (float)(_b), (float)(_c), (float)(_d), \
    (float)(_e), (float)(_f), (float)(_g), (float)(_h), \
    (float)(_i), (float)(_j), (float)(_k), (float)(_l), \
    (float)(_m), (float)(_n), (float)(_o), (float)(_p), \
})


/* zeros/ones/unit vectors */
#define V2_00 V2(0, 0)
#define V2_01 V2(0, 1)
#define V2_10 V2(1, 0)
#define V2_11 V2(1, 1)

#define V3_000 V3(0, 0, 0)
#define V3_100 V3(1, 0, 0)
#define V3_010 V3(0, 1, 0)
#define V3_001 V3(0, 0, 1)
#define V3_111 V3(1, 1, 1)

#define V4_0000 V4(0, 0, 0, 0)
#define V4_1000 V4(1, 0, 0, 0)
#define V4_0100 V4(0, 1, 0, 0)
#define V4_0010 V4(0, 0, 1, 0)
#define V4_0001 V4(0, 0, 0, 1)
#define V4_1111 V4(1, 1, 1, 1)

#define M2X2_0 M2X2(0, 0,  0, 0)
#define M2X2_1 M2X2(1, 1,  1, 1)
#define M2X2_I M2X2(1, 0,  0, 1)

#define M3X3_0 M3X3(0, 0, 0,  0, 0, 0,  0, 0, 0)
#define M3X3_1 M3X3(1, 1, 1,  1, 1, 1,  1, 1, 1)
#define M3X3_I M3X3(1, 0, 0,  0, 1, 0,  0, 0, 1)

#define M4X4_0 M4X4(0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0)
#define M4X4_1 M4X4(1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1)
#define M4X4_I M4X4(1, 0, 0, 0,  0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1)


/* common names */
#define V2_0 V2_00
#define V2_1 V2_11
#define V2_X V2_10
#define V2_Y V2_01

#define V3_0 V3_000
#define V3_1 V3_111
#define V3_X V3_100
#define V3_Y V3_010
#define V3_Z V3_001

#define V4_0 V4_0000
#define V4_1 V4_1111
#define V4_X V4_1000
#define V4_Y V4_0100
#define V4_Z V4_0010
#define V4_W V4_0001

/* directions */

#define V2_RIGHT V2(+1, 0)
#define V2_UP    V2(0, +1)
#define V2_LEFT  V2(-1, 0)
#define V2_DOWN  V2(0, -1)

#define V3_RIGHT   V3(+1, 0, 0)
#define V3_UP      V3(0, +1, 0)
#define V3_FORWARD V3(0, 0, +1)
#define V3_LEFT    V3(-1, 0, 0)
#define V3_DOWN    V3(0, -1, 0)
#define V3_BACK    V3(0, 0, -1)


/* printf -style format strings (use the expand macro in the printf macro) */
#define V2_FMT "(%+2.2f, %+2.2f)"
#define V3_FMT "(%+2.2f, %+2.2f, %+2.2f)"
#define V4_FMT "(%+2.2f, %+2.2f, %+2.2f, %+2.2f)"
#define M2X2_FMT "{(%+2.2f, %+2.2f), (%+2.2f, %+2.2f)}"
#define M3X3_FMT "{(%+2.2f, %+2.2f, %+2.2f), (%+2.2f, %+2.2f, %+2.2f), (%+2.2f, %+2.2f, %+2.2f)}"
#define M4X4_FMT "{(%+2.2f, %+2.2f, %+2.2f, %2.2f), (%+2.2f, %+2.2f, %+2.2f, %2.2f), (%+2.2f, %+2.2f, %+2.2f, %2.2f), (%+2.2f, %+2.2f, %+2.2f, %2.2f)}"

/* expansion */
#define V2__(_v2) (_v2).x, (_v2).y
#define V3__(_v3) (_v3).x, (_v3).y, (_v3).z
#define V4__(_v4) (_v4).x, (_v4).y, (_v4).z, (_v4).w
#define M2X2__(_m2x2) (_m2x2).a, (_m2x2).b, (_m2x2).c, (_m2x2).d
#define M3X3__(_m3x3) (_m3x3).a, (_m3x3).b, (_m3x3).c, (_m3x3).d, (_m3x3).e, (_m3x3).f, (_m3x3).g, (_m3x3).h, (_m3x3).i
#define M4X4__(_m4x4) (_m4x4).a, (_m4x4).b, (_m4x4).c, (_m4x4).d, (_m4x4).e, (_m4x4).f, (_m4x4).g, (_m4x4).h, (_m4x4).i, (_m4x4).j, (_m4x4).k, (_m4x4).l, (_m4x4).m, (_m4x4).n, (_m4x4).o, (_m4x4).p

/* operations */
static inline v2 v2_neg(v2 a) { return V2(-a.x, -a.y); }
static inline v3 v3_neg(v3 a) { return V3(-a.x, -a.y, -a.z); }
static inline v4 v4_neg(v4 a) { return V4(-a.x, -a.y, -a.z, -a.w); }


/* statistics/values */
static inline float v2_mag2(v2 a) { return a.x*a.x + a.y*a.y; }
static inline float v3_mag2(v3 a) { return a.x*a.x + a.y*a.y + a.z*a.z; }
static inline float v4_mag2(v4 a) { return a.x*a.x + a.y*a.y + a.z*a.z + a.w*a.w; }

static inline float v2_mag(v2 a) { return hypot(a.x, a.y); }
static inline float v3_mag(v3 a) { return sqrtf(a.x*a.x + a.y*a.y + a.z*a.z); }
static inline float v4_mag(v4 a) { return sqrtf(a.x*a.x + a.y*a.y + a.z*a.z + a.w*a.w); }


/* piecewise operations */

// compute elementwise sum of 'a' and 'b'
static inline v2 v2_add(v2 a, v2 b) { return V2(a.x+b.x, a.y+b.y); }
static inline v3 v3_add(v3 a, v3 b) { return V3(a.x+b.x, a.y+b.y, a.z+b.z); }
static inline v4 v4_add(v4 a, v4 b) { return V4(a.x+b.x, a.y+b.y, a.z+b.z, a.w+b.w); }

// compute elementwise difference of 'a' and 'b'
static inline v2 v2_sub(v2 a, v2 b) { return V2(a.x-b.x, a.y-b.y); }
static inline v3 v3_sub(v3 a, v3 b) { return V3(a.x-b.x, a.y-b.y, a.z-b.z); }
static inline v4 v4_sub(v4 a, v4 b) { return V4(a.x-b.x, a.y-b.y, a.z-b.z, a.w-b.w); }

// compute elementwise product of 'a' and 'b' (this is NOT dot/cross product)
static inline v2 v2_mul(v2 a, v2 b) { return V2(a.x*b.x, a.y*b.y); }
static inline v3 v3_mul(v3 a, v3 b) { return V3(a.x*b.x, a.y*b.y, a.z*b.z); }
static inline v4 v4_mul(v4 a, v4 b) { return V4(a.x*b.x, a.y*b.y, a.z*b.z, a.w*b.w); }

// compute elementwise division of 'a' and 'b'
static inline v2 v2_div(v2 a, v2 b) { return V2(a.x/b.x, a.y/b.y); }
static inline v3 v3_div(v3 a, v3 b) { return V3(a.x/b.x, a.y/b.y, a.z/b.z); }
static inline v4 v4_div(v4 a, v4 b) { return V4(a.x/b.x, a.y/b.y, a.z/b.z, a.w/b.w); }

// compute a vector multiplied by a scalar
static inline v2 v2_scale(v2 a, float b) { return V2(a.x*b, a.y*b); }
static inline v3 v3_scale(v3 a, float b) { return V3(a.x*b, a.y*b, a.z*b); }
static inline v4 v4_scale(v4 a, float b) { return V4(a.x*b, a.y*b, a.z*b, a.w*b); }

// compute dot product, i.e. sum(a[i]*b[i])
static inline float v2_dot(v2 a, v2 b) { return a.x*b.x + a.y*b.y; }
static inline float v3_dot(v3 a, v3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static inline float v4_dot(v4 a, v4 b) { return a.x*b.x + a.y*b.y + a.z*b.z + a.w*b.w; }
#define _MIN(_a, _b) ((_a < _b) ? _a : _b)
#define _MAX(_a, _b) ((_a > _b) ? _a : _b)

// compute elementwise minimum of vectors
static inline v2 v2_min(v2 a, v2 b) {
    return (v2){ _MIN(a.x, b.x), _MIN(a.y, b.y) };
}
static inline v3 v3_min(v3 a, v3 b) {
    return (v3){ _MIN(a.x, b.x), _MIN(a.y, b.y), _MIN(a.z, b.z) };
}
static inline v4 v4_min(v4 a, v4 b) {
    return (v4){ _MIN(a.x, b.x), _MIN(a.y, b.y), _MIN(a.z, b.z), _MIN(a.w, b.w) };
}

// compute elementwise maximum of vectors
static inline v2 v2_max(v2 a, v2 b) {
    return (v2){ _MAX(a.x, b.x), _MAX(a.y, b.y) };
}
static inline v3 v3_max(v3 a, v3 b) {
    return (v3){ _MAX(a.x, b.x), _MAX(a.y, b.y), _MAX(a.z, b.z) };
}
static inline v4 v4_max(v4 a, v4 b) {
    return (v4){ _MAX(a.x, b.x), _MAX(a.y, b.y), _MAX(a.z, b.z), _MAX(a.w, b.w) };
}

// compute cross product
static inline v3 v3_cross(v3 a, v3 b) {
    return V3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

// compute linear interpolation between vectors, i.e. (1-t)a+tb
static inline v2 v2_lerp(v2 a, v2 b, float t) {
    float tp = 1.0f - t;
    return V2(tp*a.x+t*b.x, tp*a.y+t*b.y);
}
static inline v3 v3_lerp(v3 a, v3 b, float t) {
    float tp = 1.0f - t;
    return V3(tp*a.x+t*b.x, tp*a.y+t*b.y, tp*a.z+t*b.z);
}
static inline v4 v4_lerp(v4 a, v4 b, float t) {
    float tp = 1.0f - t;
    return V4(tp*a.x+t*b.x, tp*a.y+t*b.y, tp*a.z+t*b.z, tp*a.w+t*b.w);
}


// compute euclidean distance between 2 vectors
static inline float v2_dist(v2 a, v2 b) {
    return v2_mag(v2_sub(a, b));
}
static inline float v3_dist(v3 a, v3 b) {
    return v3_mag(v3_sub(a, b));
}
static inline float v4_dist(v4 a, v4 b) {
    return v4_mag(v4_sub(a, b));
}


/* transforms */

// returns a unit vector with the same orientation
static inline v2 v2_unit(v2 a) {
    float invmag = 1.0f / v2_mag(a);

    return V2(a.x*invmag, a.y*invmag);
}
static inline v3 v3_unit(v3 a) {
    // inverse magnitude
    float invmag = 1.0f / v3_mag(a);

    return V3(a.x*invmag, a.y*invmag, a.z*invmag);
}
static inline v4 v4_unit(v4 a) {
    // inverse magnitude
    float invmag = 1.0f / v4_mag(a);
    return V4(a.x*invmag, a.y*invmag, a.z*invmag, a.w*invmag);
}


/* matrix operators */

// adds the matrices elementwise
static inline m2x2 m2x2_add(m2x2 A, m2x2 B) {
    return M2X2(
        A.a+B.a, A.b+A.b,
        A.c+B.c, A.d+B.d
    );
}
static inline m3x3 m3x3_add(m3x3 A, m3x3 B) {
    return M3X3(
        A.a+B.a, A.b+A.b, A.c+B.c, 
        A.d+B.d, A.e+B.e, A.f+B.f,
        A.g+B.g, A.h+B.h, A.i+B.i
    );
}
static inline m4x4 m4x4_add(m4x4 A, m4x4 B) {
    return M4X4(
        A.a+B.a, A.b+B.b, A.c+B.c, A.d+B.d, 
        A.e+B.e, A.f+B.f, A.g+B.g, A.h+B.h, 
        A.i+B.i, A.j+B.j, A.k+B.k, A.l+B.l,
        A.m+B.m, A.n+B.n, A.o+B.o, A.p+B.p
    );
}

// subtracts the matrices elementwise
static inline m2x2 m2x2_sub(m2x2 A, m2x2 B) {
    return M2X2(
        A.a-B.a, A.b-A.b,
        A.c-B.c, A.d-B.d
    );
}
static inline m3x3 m3x3_sub(m3x3 A, m3x3 B) {
    return M3X3(
        A.a-B.a, A.b-A.b, A.c-B.c, 
        A.d-B.d, A.e-B.e, A.f-B.f,
        A.g-B.g, A.h-B.h, A.i-B.i
    );
}
static inline m4x4 m4x4_sub(m4x4 A, m4x4 B) {
    return M4X4(
        A.a-B.a, A.b-A.b, A.c-B.c, A.d-B.d, 
        A.e-B.e, A.f-B.f, A.g-B.g, A.h-B.h, 
        A.i-B.i, A.j-B.j, A.k-B.k, A.l-B.l,
        A.m-B.m, A.n-B.n, A.o-B.o, A.p-B.p
    );
}

// matrix transposition, switching rows and columns
static inline m2x2 m2x2_transpose(m2x2 A) {
    return M2X2(
        A.a, A.c,
        A.b, A.d
    );
}
static inline m3x3 m3x3_transpose(m3x3 A) {
    return M3X3(
        A._[0], A._[3], A._[6],
        A._[1], A._[4], A._[7],
        A._[2], A._[5], A._[8]
    );
}
static inline m4x4 m4x4_transpose(m4x4 A) {
    return M4X4(
        A._[0], A._[4], A._[8],  A._[12], 
        A._[1], A._[5], A._[9],  A._[13], 
        A._[2], A._[6], A._[10], A._[14], 
        A._[3], A._[7], A._[11], A._[15]
    );
}


// matrix multiplication, dot products of rows and columns
static inline m2x2 m2x2_mul(m2x2 A, m2x2 B) {
    return (m2x2) {
        // row 1 of result
        A.a*B.a + A.b*B.c, A.a*B.b + A.b*B.d,

        // row 2 of result
        A.c*B.a + A.d*B.c, A.c*B.b + A.d*B.d,

    };
}
static inline m3x3 m3x3_mul(m3x3 A, m3x3 B) {
    return (m3x3) {
        // row 1 of result
        A.r0c0 * B.c0r0 + A.r0c1 * B.c0r1 + A.r0c2 * B.c0r2,
        A.r0c0 * B.c1r0 + A.r0c1 * B.c1r1 + A.r0c2 * B.c1r2,
        A.r0c0 * B.c2r0 + A.r0c1 * B.c2r1 + A.r0c2 * B.c2r2,

        // row 2 of result
        A.r1c0 * B.c0r0 + A.r1c1 * B.c0r1 + A.r1c2 * B.c0r2,
        A.r1c0 * B.c1r0 + A.r1c1 * B.c1r1 + A.r1c2 * B.c1r2,
        A.r1c0 * B.c2r0 + A.r1c1 * B.c2r1 + A.r1c2 * B.c2r2,

        // row 3
        A.r2c0 * B.c0r0 + A.r2c1 * B.c0r1 + A.r2c2 * B.c0r2,
        A.r2c0 * B.c1r0 + A.r2c1 * B.c1r1 + A.r2c2 * B.c1r2,
        A.r2c0 * B.c2r0 + A.r2c1 * B.c2r1 + A.r2c2 * B.c2r2,

    };
}
static inline m4x4 m4x4_mul(m4x4 A, m4x4 B) {
    return (m4x4) {
        // row 1 of result
        A.r0c0 * B.c0r0 + A.r0c1 * B.c0r1 + A.r0c2 * B.c0r2 + A.r0c3 * B.c0r3,
        A.r0c0 * B.c1r0 + A.r0c1 * B.c1r1 + A.r0c2 * B.c1r2 + A.r0c3 * B.c1r3,
        A.r0c0 * B.c2r0 + A.r0c1 * B.c2r1 + A.r0c2 * B.c2r2 + A.r0c3 * B.c2r3,
        A.r0c0 * B.c3r0 + A.r0c1 * B.c3r1 + A.r0c2 * B.c3r2 + A.r0c3 * B.c3r3,

        // row 2 of result
        A.r1c0 * B.c0r0 + A.r1c1 * B.c0r1 + A.r1c2 * B.c0r2 + A.r1c3 * B.c0r3,
        A.r1c0 * B.c1r0 + A.r1c1 * B.c1r1 + A.r1c2 * B.c1r2 + A.r1c3 * B.c1r3,
        A.r1c0 * B.c2r0 + A.r1c1 * B.c2r1 + A.r1c2 * B.c2r2 + A.r1c3 * B.c2r3,
        A.r1c0 * B.c3r0 + A.r1c1 * B.c3r1 + A.r1c2 * B.c3r2 + A.r1c3 * B.c3r3,

        // row 3
        A.r2c0 * B.c0r0 + A.r2c1 * B.c0r1 + A.r2c2 * B.c0r2 + A.r2c3 * B.c0r3,
        A.r2c0 * B.c1r0 + A.r2c1 * B.c1r1 + A.r2c2 * B.c1r2 + A.r2c3 * B.c1r3,
        A.r2c0 * B.c2r0 + A.r2c1 * B.c2r1 + A.r2c2 * B.c2r2 + A.r2c3 * B.c2r3,
        A.r2c0 * B.c3r0 + A.r2c1 * B.c3r1 + A.r2c2 * B.c3r2 + A.r2c3 * B.c3r3,

        // row 4
        A.r3c0 * B.c0r0 + A.r3c1 * B.c0r1 + A.r3c2 * B.c0r2 + A.r3c3 * B.c0r3,
        A.r3c0 * B.c1r0 + A.r3c1 * B.c1r1 + A.r3c2 * B.c1r2 + A.r3c3 * B.c1r3,
        A.r3c0 * B.c2r0 + A.r3c1 * B.c2r1 + A.r3c2 * B.c2r2 + A.r3c3 * B.c2r3,
        A.r3c0 * B.c3r0 + A.r3c1 * B.c3r1 + A.r3c2 * B.c3r2 + A.r3c3 * B.c3r3
    };
}

// matrix determinant, i.e. det(A)
static inline float m2x2_det(m2x2 A) {
    return A.a * A.d - A.b * A.c;
}
static inline float m3x3_det(m3x3 A) {
    return A.a*(A.e*A.i-A.f*A.h)-A.b*(A.d*A.i-A.f*A.g)+A.c*(A.d*A.h-A.e*A.g);
}

// matrix inverse, such that A*A^-1==I. If A is not invertible, return the zero matrix

static inline m2x2 m2x2_inv(m2x2 A) {
    float dA = m2x2_det(A);
    if (fabsf(dA) < 1e-6f) return M2X2_0;
    float scl = 1.0f / dA;
    return M2X2(
        scl*A.d, -scl*A.b,
        -scl*A.c, scl*A.a
    );
}

// inline for easier usage
#define _M2X2_DET(a, b, c, d) ((float)(((a)*(d)-(b)*(c))))

static inline m3x3 m3x3_inv(m3x3 A) {
    float dA = m3x3_det(A);
    if (fabsf(dA) < 1e-6f) return M3X3_0;
    // scaling factor
    float scl = 1.0f / dA;

    return M3X3(
        scl*_M2X2_DET(A.r1c1, A.r1c2, A.r2c1, A.r2c2),
        scl*_M2X2_DET(A.r0c2, A.r0c1, A.r2c2, A.r2c1),
        scl*_M2X2_DET(A.r0c1, A.r0c2, A.r1c1, A.r1c2),

        scl*_M2X2_DET(A.r1c2, A.r1c0, A.r2c2, A.r2c0),
        scl*_M2X2_DET(A.r0c0, A.r0c2, A.r2c0, A.r2c2),
        scl*_M2X2_DET(A.r0c2, A.r0c0, A.r1c2, A.r1c0),

        scl*_M2X2_DET(A.r1c0, A.r1c1, A.r2c0, A.r2c1),
        scl*_M2X2_DET(A.r0c1, A.r0c0, A.r2c1, A.r2c0),
        scl*_M2X2_DET(A.r0c0, A.r0c1, A.r1c0, A.r1c1)
    );
}

static inline m4x4 m4x4_inv(m4x4 A) {
    m4x4 r;
    float det;
    int i;

    r._[0] = A._[5]  * A._[10] * A._[15] - 
             A._[5]  * A._[11] * A._[14] - 
             A._[9]  * A._[6]  * A._[15] + 
             A._[9]  * A._[7]  * A._[14] +
             A._[13] * A._[6]  * A._[11] - 
             A._[13] * A._[7]  * A._[10];

    r._[4] = -A._[4]  * A._[10] * A._[15] + 
              A._[4]  * A._[11] * A._[14] + 
              A._[8]  * A._[6]  * A._[15] - 
              A._[8]  * A._[7]  * A._[14] - 
              A._[12] * A._[6]  * A._[11] + 
              A._[12] * A._[7]  * A._[10];

    r._[8] = A._[4]  * A._[9] * A._[15] - 
             A._[4]  * A._[11] * A._[13] - 
             A._[8]  * A._[5] * A._[15] + 
             A._[8]  * A._[7] * A._[13] + 
             A._[12] * A._[5] * A._[11] - 
             A._[12] * A._[7] * A._[9];

    r._[12] = -A._[4]  * A._[9] * A._[14] + 
               A._[4]  * A._[10] * A._[13] +
               A._[8]  * A._[5] * A._[14] - 
               A._[8]  * A._[6] * A._[13] - 
               A._[12] * A._[5] * A._[10] + 
               A._[12] * A._[6] * A._[9];

    r._[1] = -A._[1]  * A._[10] * A._[15] + 
              A._[1]  * A._[11] * A._[14] + 
              A._[9]  * A._[2] * A._[15] - 
              A._[9]  * A._[3] * A._[14] - 
              A._[13] * A._[2] * A._[11] + 
              A._[13] * A._[3] * A._[10];

    r._[5] = A._[0]  * A._[10] * A._[15] - 
             A._[0]  * A._[11] * A._[14] - 
             A._[8]  * A._[2] * A._[15] + 
             A._[8]  * A._[3] * A._[14] + 
             A._[12] * A._[2] * A._[11] - 
             A._[12] * A._[3] * A._[10];

    r._[9] = -A._[0]  * A._[9] * A._[15] + 
              A._[0]  * A._[11] * A._[13] + 
              A._[8]  * A._[1] * A._[15] - 
              A._[8]  * A._[3] * A._[13] - 
              A._[12] * A._[1] * A._[11] + 
              A._[12] * A._[3] * A._[9];

    r._[13] = A._[0]  * A._[9] * A._[14] - 
              A._[0]  * A._[10] * A._[13] - 
              A._[8]  * A._[1] * A._[14] + 
              A._[8]  * A._[2] * A._[13] + 
              A._[12] * A._[1] * A._[10] - 
              A._[12] * A._[2] * A._[9];

    r._[2] = A._[1]  * A._[6] * A._[15] - 
             A._[1]  * A._[7] * A._[14] - 
             A._[5]  * A._[2] * A._[15] + 
             A._[5]  * A._[3] * A._[14] + 
             A._[13] * A._[2] * A._[7] - 
             A._[13] * A._[3] * A._[6];

    r._[6] = -A._[0]  * A._[6] * A._[15] + 
              A._[0]  * A._[7] * A._[14] + 
              A._[4]  * A._[2] * A._[15] - 
              A._[4]  * A._[3] * A._[14] - 
              A._[12] * A._[2] * A._[7] + 
              A._[12] * A._[3] * A._[6];

    r._[10] = A._[0]  * A._[5] * A._[15] - 
              A._[0]  * A._[7] * A._[13] - 
              A._[4]  * A._[1] * A._[15] + 
              A._[4]  * A._[3] * A._[13] + 
              A._[12] * A._[1] * A._[7] - 
              A._[12] * A._[3] * A._[5];

    r._[14] = -A._[0]  * A._[5] * A._[14] + 
               A._[0]  * A._[6] * A._[13] + 
               A._[4]  * A._[1] * A._[14] - 
               A._[4]  * A._[2] * A._[13] - 
               A._[12] * A._[1] * A._[6] + 
               A._[12] * A._[2] * A._[5];

    r._[3] = -A._[1] * A._[6] * A._[11] + 
              A._[1] * A._[7] * A._[10] + 
              A._[5] * A._[2] * A._[11] - 
              A._[5] * A._[3] * A._[10] - 
              A._[9] * A._[2] * A._[7] + 
              A._[9] * A._[3] * A._[6];

    r._[7] = A._[0] * A._[6] * A._[11] - 
             A._[0] * A._[7] * A._[10] - 
             A._[4] * A._[2] * A._[11] + 
             A._[4] * A._[3] * A._[10] + 
             A._[8] * A._[2] * A._[7] - 
             A._[8] * A._[3] * A._[6];

    r._[11] = -A._[0] * A._[5] * A._[11] + 
               A._[0] * A._[7] * A._[9] + 
               A._[4] * A._[1] * A._[11] - 
               A._[4] * A._[3] * A._[9] - 
               A._[8] * A._[1] * A._[7] + 
               A._[8] * A._[3] * A._[5];

    r._[15] = A._[0] * A._[5] * A._[10] - 
              A._[0] * A._[6] * A._[9] - 
              A._[4] * A._[1] * A._[10] + 
              A._[4] * A._[2] * A._[9] + 
              A._[8] * A._[1] * A._[6] - 
              A._[8] * A._[2] * A._[5];

    det = A._[0] * r._[0] + A._[1] * r._[4] + A._[2] * r._[8] + A._[3] * r._[12];

    if (det == 0.0f) return M4X4_0;

    det = 1.0f / det;

    for (i = 0; i < 16; i++)
        r._[i] = r._[i] * det;

    return r;
}


/* misc. matrix multiplication sub-forms */

// return A * b (b being a column vector)
static inline v3 m3x3_mul_v3(m3x3 A, v3 b) {
    return V3(v3_dot(A.rows[0], b), v3_dot(A.rows[1], b), v3_dot(A.rows[2], b));
}



/* misc. graphics transforms & utilities */

// reflect a vector (I: incidence) over another vector (N: normal)
static inline v3 v3_refl(v3 I, v3 N) {
    return v3_add(I, v3_scale(N, -2.0f * v3_dot(I, N)));
}


/* geometric transforms (m4T_*) 

These functions create matrices `M` such that:

M*(x, y, z, 1) gives a transformation of xyz into M's space

In general it is like:

[xdir.xyz, xoff]
[ydir.xyz, yoff]
[zdir.xyz, zoff]
[0, 0, 0, 1]

*dir are 1x3 vectors that describe `M`'s x, y, or z vector.

*off are the offsets in each direction

Default transform is the identity matrix

*/

// create a transform with an XYZ position. No rotation, no scaling
static inline m4x4 m4T_xyz(float x, float y, float z) {
    return M4X4(
        1, 0, 0, x,
        0, 1, 0, y,
        0, 0, 1, z,
        0, 0, 0, 1
    );
}

// create a transform with XYZ scaling
static inline m4x4 m4T_scale(float x, float y, float z) {
    return M4X4(
        x, 0, 0, 0,
        0, y, 0, 0,
        0, 0, z, 0,
        0, 0, 0, 1
    );
}

// create a transform representing a rotation about the X axis (same as a given pitch)
// obviously, it is in radians
static inline m4x4 m4T_rot_x(float x /*or pitch*/) {
    // compute sin and cosine of the angle
    float sinA = sinf(x), cosA = cosf(x);
    return M4X4(
        1,     0,     0, 0,
        0,  cosA, -sinA, 0,
        0,  sinA,  cosA, 0,
        0,     0,     0, 1
    );
}

// create a transform representing a rotation about the X axis (same as a given pitch)
// obviously, it is in radians
static inline m4x4 m4T_rot_y(float y /*or yaw*/) {
    // compute sin and cosine of the angle
    float sinA = sinf(y), cosA = cosf(y);
    return M4X4(
        cosA,  0, sinA, 0,
        0,     1,    0, 0,
        -sinA, 0, cosA, 0,
        0,     0,    0, 1
    );
}

// create a transform representing a rotation about the Z axis (same as a given roll)
// obviously, it is in radians
static inline m4x4 m4T_rot_z(float z /*or roll*/) {
    // compute sin and cosine of the angle
    float sinA = sinf(z), cosA = cosf(z);
    return M4X4(
        cosA, -sinA, 0, 0,
        sinA,  cosA, 0, 0,
        0,     0,    1, 0,
        0,     0,    0, 1
    );
}


// transforms b into a linear space defined by: A (with the upper-left 3x3 mat being row transforms,
// and and the upper left 3x1 vector is the xyz offset), the W row should be (0, 0, 0, 1)
static inline v3 m4T_transform_direction(m4x4 A, v3 b) {
    return v3_add(v3_add(v3_scale(A.X.xyz, b.x), v3_scale(A.Y.xyz, b.y)), v3_scale(A.Z.xyz, b.z));
}

/* geometric sampling methods */

// generates a unit vector which lies on the unit sphere, given a Z coordinate and an angle
static inline v3 v3_sphere_sample_ZA(float z, float angle) {
    float x = sinf(angle), y = sinf(angle), scl = sqrtf(1.0f - z * z);
    return V3(scl * x, scl * y, z);
}

// creates a view transformation matrix for openGL
static inline m4x4 m4T_view_lookat(v3 pos, v3 target, v3 updir) {
    v3 forward = v3_unit(v3_sub(target, pos));
    v3 right = v3_cross(forward, updir);
    v3 up = v3_cross(right, forward);

    return M4X4(
        right.x, right.y, right.z, -v3_dot(right, pos),
        up.x, up.y, up.z, -v3_dot(up, pos),
        -forward.x, -forward.y, -forward.z, v3_dot(forward, pos),
        0, 0, 0, 1
    );

    m4x4 rT = M4X4_I;
    rT.X.w = -pos.x;
    rT.Y.w = -pos.y;
    rT.Z.w = -pos.z;

    m4x4 rC = M4X4_I;
    rC.X.xyz = right;
    rC.Y.xyz = up;
    rC.Z.xyz = forward;

    return m4x4_mul(rC, rT);

}

// creates a projection matrix for openGL
static inline m4x4 m4T_proj(float FOV, float aspect, float Znear, float Zfar) {
    float tanHalfFOV = tanf(3.1415926535 * FOV / 180.0f / 2.0f);

    m4x4 ret = M4X4_0;
    ret.r0c0 = 1.0f / (aspect * tanHalfFOV);
    ret.r1c1 = 1.0f / (tanHalfFOV);
    ret.r2c2 = -(Zfar+Znear)/(Zfar-Znear);
    ret.r3c2 = -1.0f;
    ret.r2c3 = - 2.0f * Zfar * Znear / (Zfar - Znear);
    return ret;
}


#endif


