//===- drawcontext.cpp -------------------------------------------*- C++ --*-===//
// Copyright 2017  Warren Hunt
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
//===----------------------------------------------------------------------===//

// modified 10.10.2021 (Wilhelm Schnepp): change include path to match adjusted cmake file structure
// modified 10.10.2021 (Wilhelm Schnepp): Extract ShapeBuilder and FontRenderInfo to separate files

#include <thorax_truetype/renderer/build.h>
#include <thorax_truetype/renderer/drawcontext.h>

#include <utility>
#include <cassert>
#include <cmath>
#include <cstdio>

#define DEBUG_NOINLINE
//#define DEBUG_NOINLINE __declspec(noinline)

struct UTMatrix {
    UTMatrix(float x_x, const Vector& y, const Point& w = Point(0, 0)) : x_x(x_x), y(y), w(w) {}
    Vector operator*(const Vector& a) const { return Vector(x_x * a.x + y.x * a.y, y.y * a.y); }
    Point operator*(const Point& a) const { return Point(x_x * a.x + y.x * a.y + w.x, y.y * a.y + w.y); }
    UTMatrix operator*(const UTMatrix& a) const { return UTMatrix(x_x * a.x_x, (*this) * a.y, (*this) * a.w); }
    float x_x;
    Vector y;
    Point w;
};

static inline UTMatrix invert(const UTMatrix& a) {
    auto scale = 1.0f / (a.x_x * a.y.y);
    auto x_x = a.y.y * scale;
    auto y = Vector(-a.y.x, a.x_x) * scale;
    auto w = Point(y.x * -a.w.y - x_x * a.w.x, y.y * -a.w.y);
    return UTMatrix(x_x, y, w);
}

static inline UTMatrix ComputeViewFromSample(Vector x, Vector y) {
    if (std::abs(y.x) > std::abs(x.x)) std::swap(x, y);
    auto xy = x.y;
    auto yy = y.y;
    // TODO: still buggy for smaller distortions
    if (std::abs(xy) * 30 < std::abs(x.x)) xy = 0; // Add dead-zone for highly anisotripic patterns.
    // TODO: guarantee that out.y.y > 0
    auto scale = 1 / std::sqrt(xy * xy + yy * yy);
    xy *= scale;
    yy *= scale;
    return UTMatrix(x.x * yy - y.x * xy, x * xy + y * yy);
}

//==============================================================================
// TileImpl
//==============================================================================

static const size_t STAMP_WIDTH = 4;
static const size_t STAMP_HEIGHT = 2;
static const size_t TILE_WIDTH = 128;
static const size_t TILE_HEIGHT = 64;
static const ptrdiff_t TILE_STRIDE = TILE_WIDTH + STAMP_WIDTH;
static const size_t TILE_TOTAL_HEIGHT = TILE_HEIGHT + STAMP_HEIGHT - 1;
static const size_t TOTAL_TILE_SIZE = TILE_STRIDE * TILE_TOTAL_HEIGHT;

struct Tile {
    Tile(int x, int y, ptrdiff_t renderTargetStride, const Box& viewport, int tileWidth, int tileHeight)
        : pixelport(viewport & Box((float)x, (float)y, (float)x + (float)tileWidth, (float)y + (float)tileHeight)),
          tileAddressOffset(-(y * TILE_STRIDE + x)),
          renderTargetOffset(y * renderTargetStride + x),
          tileWidth(tileWidth),
          tileHeight(tileHeight) {}
    Box pixelport;
    int tileWidth;
    int tileHeight;
    ptrdiff_t tileAddressOffset;
    ptrdiff_t renderTargetOffset;
};

struct ThreadState {
    float8 red;
    float8 green;
    float8 blue;
    short* coverage;
    const Tile* tile;
};

struct PackedMatrix2x3 {
    PackedMatrix2x3(float x_x, float y_x, float w_x, float x_y, float y_y, float w_y)
        : xxyyww00(x_x, x_y, y_x, y_y, w_x, w_y, 0, 0) {}
    explicit PackedMatrix2x3(const Matrix2x3& a) {
        static const __declspec(align(32)) int MASK[] = { -1, -1, -1, -1, -1, -1, 0, 0 };
        xxyyww00 = float8::LoadU(&a.x.x) & *(const float8*)MASK;
    }
    explicit PackedMatrix2x3(float8 xxyyww00) : xxyyww00(xxyyww00) {}
    PackedMatrix2x3 operator*(const PackedMatrix2x3& a) const {
        auto v0v2 = shuffle4<0, 0>(xxyyww00) * shuffle<0, 0, 1, 1>(a.xxyyww00);
        auto v1v3 = shuffle4<0, 0>(xxyyww00) * shuffle<2, 2, 3, 3>(a.xxyyww00);
        return PackedMatrix2x3(shuffle<0, 1, 0, 1>(v0v2, v1v3) + shuffle<2, 3, 2, 3>(v0v2, v1v3) +
            shuffle4<8, 1>(xxyyww00));
    }
    void ComputeBounds(float8& mins, float8& maxs) const {
        auto v0123 = shuffle2<2, 2, 2, 2>(xxyyww00) + shuffle2<3, 0, 3, 0>(xxyyww00) + shuffle2<3, 3, 1, 1>(xxyyww00);
        auto v01 = v0123.v0123();
        auto v23 = v0123.v4567();
        auto vmin = min(v01, v23);
        vmin = min(vmin, shuffle<2, 3, 0, 1>(vmin));
        mins = float8(shuffle<0, 0, 0, 0>(vmin), shuffle<1, 1, 1, 1>(vmin));
        auto vmax = max(v01, v23);
        vmax = max(vmax, shuffle<2, 3, 0, 1>(vmax));
        maxs = float8(shuffle<0, 0, 0, 0>(vmax), shuffle<1, 1, 1, 1>(vmax));
    }
    float8 xxyyww00;
};

// Effectively a Matrix2x3 that can apply to 2 points and 2 vectors simultaneiously (matches the format of Shape)
struct ShapeTransform {
    ShapeTransform() = default;
    explicit ShapeTransform(const UTMatrix& a)
        : x(_mm256_castsi256_ps(_mm256_srli_epi64(_mm256_castps_si256(_mm256_broadcast_ss(&a.x_x)), 32))),
        y(_mm256_castpd_ps(_mm256_broadcast_sd((double*)&a.y.x))),
        w(_mm256_castsi256_ps(_mm256_srli_si256(_mm256_castpd_si256(_mm256_broadcast_sd((double*)&a.w.x)), 8))) {}
    explicit ShapeTransform(const Matrix2x3& a)
        : x(_mm256_castpd_ps(_mm256_broadcast_sd((double*)&a.x.x))),
        y(_mm256_castpd_ps(_mm256_broadcast_sd((double*)&a.y.x))),
        w(_mm256_castsi256_ps(_mm256_srli_si256(_mm256_castpd_si256(_mm256_broadcast_sd((double*)&a.w.x)), 8))) {}
    Shape operator*(const Shape& a) {
        return Shape(madd(shuffle<0, 0, 2, 2>(a.p0v1p2v3), x, madd(shuffle<1, 1, 3, 3>(a.p0v1p2v3), y, w)));
    }
    float8 x;
    float8 y;
    float8 w;
};

//==============================================================================
// DrawContextImpl
//==============================================================================

struct DrawContextImpl final : DrawContext {
    size_t SetRenderTarget(unsigned* renderTarget, int width, int height, size_t stride) override;
    void SetViewport(float lower_x, float lower_y, float upper_x, float upper_y) override;
    void SetFilterKernel(float x0, float y0, float x1, float y1) override;
    void Draw(const Mesh* scene, const Matrix2x3& worldFromScene) override;

    DrawContextImpl(unsigned numThreads);
    bool DrawTrapezoid(short* out, Shape trapezoid, const Box& pixelport);
    bool DrawCurve(short* out, Shape curve, const Box& pixelport);
    void Trace(const ThreadState& state, const Mesh* scene, const Box4Node* node, PackedMatrix2x3 objectFromWorld, PackedMatrix2x3 tileBounds);
    void Trace(const ThreadState& state, const Mesh* scene, const Box4Node* node, PackedMatrix2x3 objectFromWorld, PackedMatrix2x3 tileBounds, float8 mins, float8 maxs);
    void DrawTile(const ThreadState& state, const Mesh* scene, PackedMatrix2x3 layoutFromScreen);

    // RenderTarget related values.
    struct {
        unsigned* data;
        ptrdiff_t stride;
        int width;
        int height;
    } renderTarget;

    // Viewport related values.
    ShapeTransform screenFromWorld;
    Matrix2x3 worldFromScreen;

    // Filter kernel related values.
    ShapeTransform unitFromPixel;
    Box sampleBounds;

    // Local tile storage used during drawing.
    Array<short> coverageBuffer;
    Array<ThreadState> threadStates;
    Array<Tile> tiles;
};

DrawContextImpl::DrawContextImpl(unsigned numThreads) {
    coverageBuffer = Array<short>(TOTAL_TILE_SIZE * numThreads, 32);
    threadStates = Array<ThreadState>(numThreads);
    for (size_t i = 0; i < numThreads; i++)
        threadStates[i].coverage = coverageBuffer + TOTAL_TILE_SIZE * i;
}

DrawContext* DrawContext::Create(unsigned numThreads) {
    auto p = (DrawContextImpl*)_aligned_malloc(sizeof(DrawContextImpl), __alignof(DrawContextImpl));
    new (p) DrawContextImpl(numThreads);
    return p;
}

void DrawContext::Destroy(DrawContext*& context) {
    if (!context) return;
    _aligned_free(context);
    context = nullptr;
}

size_t DrawContextImpl::SetRenderTarget(unsigned* data, int width, int height, size_t stride) {
    // TODO: validate inputs.
    // Set up render-target.
    renderTarget.data = data + (height - 1) * stride;
    renderTarget.width = width;
    renderTarget.height = height;
    renderTarget.stride = -(ptrdiff_t)stride;

    // Initialize the tiles.
    auto viewport = Box(0, 0, (float)width, (float)height);
    auto numTiles = ((width + TILE_WIDTH - 1) / TILE_WIDTH) * ((height + TILE_HEIGHT - 1) / TILE_HEIGHT);
    tiles = Array<Tile>(numTiles);
    auto nextTile = tiles.data;
    for (int y = 0; y < height; y += (int)TILE_HEIGHT) {
        auto tileHeight = height - y < TILE_HEIGHT ? height - y : (int)TILE_HEIGHT;
        for (int x = 0; x < width; x += (int)TILE_WIDTH) {
            auto tileWidth = width - x < TILE_WIDTH ? width - x : (int)TILE_WIDTH;
            *nextTile++ = Tile(x, y, renderTarget.stride, viewport, tileWidth, tileHeight);
        }
    }
    return numTiles;
}

void DrawContextImpl::SetViewport(float lower_x, float lower_y, float upper_x, float upper_y) {
    // TODO: validate inputs
    auto xform =
        UTMatrix((upper_x - lower_x) / (float)renderTarget.width,
                 Vector(0, (upper_y - lower_y) / (float)renderTarget.height), Point(lower_x, lower_y));
    worldFromScreen = Matrix2x3(xform.x_x, xform.y.x, xform.w.x, 0, xform.y.y, xform.w.y);
    screenFromWorld = ShapeTransform(invert(xform));
}

void DrawContextImpl::SetFilterKernel(float x0, float y0, float x1, float y1) {
    // The spaces used during rendering:
    // World: origin is world origin, pixel spacing application defined
    // Screen : origin is viewport.lower, pixel spacing is 1
    // Pixel : pixel spacing is 1, the current pixel center is (0.5, 0.5)
    // View : same as pixel but the current pixel center is (0, 0)
    // Sample : pixel is at the origin, integration window is (-0.5, -0.5) -- (0.5, 0.5)
    // Unit: pixel integration window is (0, 0) -- (1, 1)

    // TODO: validate inputs
    // Compute unitFromPixel
    auto pixelFromView = UTMatrix(1, Vector(0, 1), Point(0.5f, 0.5f));
    auto viewFromSample = ComputeViewFromSample(Vector(x0, y0), Vector(x1, y1));
    auto sampleFromUnit = UTMatrix(1, Vector(0, 1), Point(-0.5f, -0.5f));
    auto pixelFromUnit = pixelFromView * viewFromSample * sampleFromUnit;
    unitFromPixel = ShapeTransform(invert(pixelFromUnit));

    // Compute sampleBounds
    auto v = viewFromSample;
    auto sampleDelta = Vector(std::abs(v.x_x) + std::abs(v.y.x), std::abs(v.y.y)) * 0.5f;
    sampleBounds = Box(-sampleDelta.x, -sampleDelta.y, sampleDelta.x, sampleDelta.y);
}

//==============================================================================
// Wide data structures for drawing.
//==============================================================================

struct Range8 {
    Range8() = default;
    Range8(float8 lower, float8 upper) : lower(lower), upper(upper) {}
    static Range8 Make(float8 a, float8 b) { return Range8(min(a, b), max(a, b)); }
    float8 lower, upper;
};

__forceinline static Range8 Cubic(const Range8& a) {
    auto cubic = [](float8 x) { return nmadd(x, float8(2), float8(3)) * (x * x); };
    return Range8(cubic(a.lower), cubic(a.upper));
}

__forceinline static Range8 Clamp(const Range8& x, const Range8& bounds = Range8(float8::Zero(), float8(1))) {
    return Range8(min(max(x.lower, bounds.lower), bounds.upper), min(max(x.upper, bounds.lower), bounds.upper));
}

struct LinearEqn8 {
    LinearEqn8(float8 a, float8 b, float8 i_a, const Range8& solve) : a(a), b(b), i_a(i_a), solve(solve) {}
    DEBUG_NOINLINE LinearEqn8(float A, float B) {
        a = float8(A);
        b = float8(B);
        // TODO: fix the case where A == 0
        i_a = float8(-1.0f / A);
        auto i_b = i_a * b;
        solve = Range8::Make(i_b, i_b - i_a);
    }
    Range8 operator()(const Range8& t) const { return Range8(madd(a, t.lower, b), madd(a, t.upper, b)); }
    __forceinline LinearEqn8 Apply(float8 x) const {
        return LinearEqn8(a, b + x, i_a, Range8(madd(i_a, x, solve.lower), madd(i_a, x, solve.upper)));
    }
    Range8 Solve() const { return solve; }
    float8 a, b;
    float8 i_a;
    Range8 solve;
};

struct QuadraticEqn8 {
    QuadraticEqn8(float8 a, float8 b, float8 c, float8 i_a, float8 i_b, const Range8& solve)
        : a(a), b(b), c(c), i_a(i_a), i_b(i_b), solve(solve) {}
    DEBUG_NOINLINE QuadraticEqn8(float A, float B, float C) {
        a = float8(A);
        b = float8(B);
        c = float8(C);
        // TODO: fix the case where A == 0
        i_a = float8(-1.0f / A);
        // TODO: write comment about the B == 0 branch.
        i_b = B == 0 ? float8::SignMask() : float8(0.5f * B) * i_a;
        auto i_c = madd(i_b, i_b, i_a * float8(C));
        solve = Range8::Make(i_c, i_c - i_a);
    }
    float8 operator()(float8 t) const { return madd(madd(a, t, b), t, c); }
    Range8 operator()(const Range8& t) const { return Range8(operator()(t.lower), operator()(t.upper)); }
    float8 SolveApex() const { return i_b; }
    Range8 SolveDelta() const {
        return Range8(sqrt(max(float8::Zero(), solve.lower)), sqrt(max(float8::Zero(), solve.upper)));
    }
    __forceinline QuadraticEqn8 Apply(float8 x) const {
        return QuadraticEqn8(a, b, c + x, i_a, i_b, Range8(madd(i_a, x, solve.lower), madd(i_a, x, solve.upper)));
    }
    float8 a, b, c;
    float8 i_a, i_b;
    Range8 solve;
};

struct TrapezoidEqns8 {
    explicit TrapezoidEqns8(const Shape& a)
        : y0(a.p0v1p2v3[3], a.p0v1p2v3[1]),
          x0(a.p0v1p2v3[2], a.p0v1p2v3[0]),
          x1(-a.p0v1p2v3[6], a.p0v1p2v3[4] + a.p0v1p2v3[6]) {}
    LinearEqn8 y0, x0, x1;
};

struct CurveEqns8 {
    explicit CurveEqns8(const Shape& a)
        : x0(-a.p0v1p2v3[6], a.p0v1p2v3[0] + a.p0v1p2v3[6]),
          y0(-a.p0v1p2v3[7], a.p0v1p2v3[1] + a.p0v1p2v3[7]),
          x1(a.p0v1p2v3[6] - 2 * a.p0v1p2v3[2], 2 * a.p0v1p2v3[2], a.p0v1p2v3[0]),
          y1(a.p0v1p2v3[7] - 2 * a.p0v1p2v3[3], 2 * a.p0v1p2v3[3], a.p0v1p2v3[1]) {}
    LinearEqn8 x0, y0;
    QuadraticEqn8 x1, y1;
};

//==============================================================================
// Drawing
//==============================================================================

__forceinline static bool SetupDraw(const Box& box,
                                    int& i_steps,
                                    int& j_steps,
                                    ptrdiff_t& p_offset,
                                    ptrdiff_t& p_stride,
                                    float8& x_start,
                                    float8& y_start) {
    auto origin = round_up(box.data);
    auto rect = int4(origin);
    auto mask = int4::True(rect.data); // create a mask of all 1s
    auto delta = rect + shuffle<2, 3, 0, 1>(rect);
    if (movemask(mask + delta)) return false;
    x_start = float8(origin[0]) - float8(0, 1, 2, 3, 0, 1, 2, 3);
    y_start = float8(origin[1]) - float8(0, 0, 0, 0, 1, 1, 1, 1);
    delta = delta + int4(STAMP_WIDTH - 1, STAMP_HEIGHT - 1, 0, 0);
    i_steps = delta[0] / STAMP_WIDTH;
    j_steps = delta[1] / STAMP_HEIGHT;
    p_offset = -rect[1] * TILE_STRIDE - rect[0];
    p_stride = STAMP_HEIGHT * TILE_STRIDE + STAMP_WIDTH - STAMP_WIDTH * i_steps;
    return true;
}

__forceinline static void WriteAVX(short* out, size_t stride, float8 value) {
    auto x = int8(float8(4096.0f) * value);
    auto y = short8::Pack(x);
    y = y + short8(*(long long*)out, *(long long*)(out + stride));
    *(long long*)out = y.v0123();
    *(long long*)(out + stride) = y.v4567();
}

DEBUG_NOINLINE float8 Integrate(const TrapezoidEqns8& eqns, float8 offset_x, float8 offset_y) {
    auto y0 = eqns.y0.Apply(offset_y);
    auto ty0 = Clamp(y0.Solve());
    auto yy0 = Cubic(y0(ty0));

    float8 area;
    {
        auto x0 = eqns.x0.Apply(offset_x);
        auto tx0 = Clamp(x0.Solve(), ty0);
        auto yx0 = Cubic(y0(tx0));
        auto yx0_mid = yx0.lower + yx0.upper;
        auto xx0 = Cubic(Clamp(x0(tx0)));
        area = xx0.lower * msub(float8(0.5f), yx0_mid, yy0.lower);
        area = madd(xx0.upper, nmadd(float8(0.5f), yx0_mid, yy0.upper), area);
    }
    {
        auto x1 = eqns.x1.Apply(offset_x);
        auto tx1 = Clamp(x1.Solve(), ty0);
        auto yx1 = Cubic(y0(tx1));
        auto yx1_mid = yx1.lower + yx1.upper;
        auto xx1 = Cubic(Clamp(x1(tx1)));
        area = nmadd(xx1.lower, msub(float8(0.5f), yx1_mid, yy0.lower), area);
        area = nmadd(xx1.upper, nmadd(float8(0.5f), yx1_mid, yy0.upper), area);
    }
    return area;
}

bool DrawContextImpl::DrawTrapezoid(short* out, Shape trapezoid, const Box& pixelport) {
    static const float8 x_step = float8(-(float)STAMP_WIDTH);
    static const float8 y_step = float8(-(float)STAMP_HEIGHT);
    int i_steps, j_steps;
    ptrdiff_t offset, stride;
    float8 x_start, y_start;

    trapezoid = screenFromWorld * trapezoid;
    auto box = trapezoid.Bound() + sampleBounds & pixelport;
    if (!SetupDraw(box, i_steps, j_steps, offset, stride, x_start, y_start)) return false;
    trapezoid = unitFromPixel * trapezoid; //< note that pixelFromScreen has been factored into the loop.
    TrapezoidEqns8 eqns(trapezoid);

    auto unitFromPixel_x_x = float8(unitFromPixel.x[0]);
    auto unitFromPixel_y_x = float8(unitFromPixel.y[0]);
    auto unitFromPixel_y_y = float8(unitFromPixel.y[1]);

    auto pixel_y = y_start;
    out += offset;
    auto j = j_steps;
    goto J_ENTRY;
    do {
        pixel_y += y_step;
        out += stride;
    J_ENTRY:
        auto unit_x = unitFromPixel_y_x * pixel_y;
        auto unit_y = unitFromPixel_y_y * pixel_y;
        auto pixel_x = x_start;
        auto i = i_steps;
        goto I_ENTRY;
        do {
            pixel_x += x_step;
            out += STAMP_WIDTH;
        I_ENTRY:
            auto value = Integrate(eqns, madd(unitFromPixel_x_x, pixel_x, unit_x), unit_y);
            WriteAVX(out, TILE_STRIDE, value);
        } while (--i);
    } while (--j);
    return true;
}

DEBUG_NOINLINE float8 Integrate(const CurveEqns8& eqns, float8 offset_x, float8 offset_y) {
    auto y0 = eqns.y0.Apply(offset_y);
    auto ty0 = Clamp(y0.Solve());
    auto yy0 = Cubic(y0(ty0));

    float8 area;
    {
        auto x0 = eqns.x0.Apply(offset_x);
        auto tx0 = Clamp(x0.Solve(), ty0);
        auto yx0 = Cubic(y0(tx0));
        auto yx0_mid = yx0.lower + yx0.upper;
        auto xx0 = Cubic(Clamp(x0(tx0)));
        area = xx0.lower * msub(float8(0.5f), yx0_mid, yy0.lower);
        area = madd(xx0.upper, nmadd(float8(0.5f), yx0_mid, yy0.upper), area);
    }
    {
        auto y1 = eqns.y1.Apply(offset_y);
        auto ty1_apex = y1.SolveApex();
        auto ty1_delta = y1.SolveDelta();
        auto ty1_apex_sign = ty1_apex & float8::SignMask();
        auto ty1 = Clamp(Range8(ty1_apex - (ty1_apex_sign | blend(ty1_delta.upper, ty1_delta.lower, ty1_apex_sign)),
                                ty1_apex - (ty1_apex_sign | blend(ty1_delta.lower, ty1_delta.upper, ty1_apex_sign))));

        auto x1 = eqns.x1.Apply(offset_x);
        auto tx1_apex = x1.SolveApex();
        auto tx1_delta = x1.SolveDelta();
        auto tx1_lower = Clamp(Range8(tx1_apex - tx1_delta.upper, tx1_apex - tx1_delta.lower), ty1);
        auto tx1_upper = Clamp(Range8(tx1_apex + tx1_delta.lower, tx1_apex + tx1_delta.upper), ty1);

        auto yx1_lower = Cubic(y1(tx1_lower));
        auto yx1_lower_mid = yx1_lower.lower + yx1_lower.upper;
        auto xx1_lower = Cubic(Clamp(x1(tx1_lower)));
        area = madd(xx1_lower.lower, msub(float8(0.5f), yx1_lower_mid, yy0.upper), area);
        area = madd(xx1_lower.upper, nmadd(float8(0.5f), yx1_lower_mid, yx1_lower.upper), area);

        auto yx1_upper = Cubic(y1(tx1_upper));
        auto yx1_upper_mid = yx1_upper.lower + yx1_upper.upper;
        auto xx1_upper = Cubic(Clamp(x1(tx1_upper)));
        area = madd(xx1_upper.lower, msub(float8(0.5f), yx1_upper_mid, yx1_lower.upper), area);
        area = madd(xx1_upper.upper, nmadd(float8(0.5f), yx1_upper_mid, yy0.lower), area);
    }
    return area;
}

bool DrawContextImpl::DrawCurve(short* out, Shape curve, const Box& pixelport) {
    static const float8 x_step = float8(-(float)STAMP_WIDTH);
    static const float8 y_step = float8(-(float)STAMP_HEIGHT);
    int i_steps, j_steps;
    ptrdiff_t offset, stride;
    float8 x_start, y_start;

    curve = screenFromWorld * curve;
    auto box = curve.Bound() + sampleBounds & pixelport;
    if (!SetupDraw(box, i_steps, j_steps, offset, stride, x_start, y_start)) return false;
    curve = unitFromPixel * curve; //< note that pixelFromScreen has been factored into the loop.
    CurveEqns8 eqns(curve);

    auto unitFromPixel_x_x = float8(unitFromPixel.x[0]);
    auto unitFromPixel_y_x = float8(unitFromPixel.y[0]);
    auto unitFromPixel_y_y = float8(unitFromPixel.y[1]);

    auto pixel_y = y_start;
    out += offset;
    auto j = j_steps;
    goto J_ENTRY;
    do {
        pixel_y += y_step;
        out += stride;
    J_ENTRY:
        auto unit_x = unitFromPixel_y_x * pixel_y;
        auto unit_y = unitFromPixel_y_y * pixel_y;
        auto pixel_x = x_start;
        auto i = i_steps;
        goto I_ENTRY;
        do {
            pixel_x += x_step;
            out += STAMP_WIDTH;
        I_ENTRY:
            auto value = Integrate(eqns, madd(unitFromPixel_x_x, pixel_x, unit_x), unit_y);
            WriteAVX(out, TILE_STRIDE, value);
        } while (--i);
    } while (--j);
    return true;
}

DEBUG_NOINLINE static void ClearCoverageTile(const short* data) {
    auto count = TILE_STRIDE * TILE_HEIGHT;
    auto p = (char*)(data + count);
    auto i = 0 - sizeof(short) * count;
    auto zero = float8::Zero();
    do {
        *(float8*)(p + i) = zero;
    } while (i += sizeof(float8));
}

DEBUG_NOINLINE static void ResolveCoverageTileAVX(
    unsigned* out, ptrdiff_t stride, const short* local, int width, int height, float8 r, float8 g, float8 b) {
    stride -= width;
    // TODO: fix me to use masked writes and support non-multiples of 16.
    do {
        size_t i = width / 16;
        do {
            auto x = abs(short16::loadu(local));
            x = min(srli<4>(x - srli<7>(x)), short16(0xff));
            auto lo = x.ZeroExtendUnpackLo();
            auto hi = x.ZeroExtendUnpackHi();
            // TODO: remove AVX instructions by using int8
            auto loR = _mm256_cvtps_epi32((float8(_mm256_cvtepi32_ps(lo.data)) * r).data);
            auto loG = _mm256_cvtps_epi32((float8(_mm256_cvtepi32_ps(lo.data)) * g).data);
            auto loB = _mm256_cvtps_epi32((float8(_mm256_cvtepi32_ps(lo.data)) * b).data);
            lo = _mm256_or_si256(loR, _mm256_or_si256(_mm256_slli_epi32(loG, 8), _mm256_slli_epi32(loB, 16)));
            auto hiR = _mm256_cvtps_epi32((float8(_mm256_cvtepi32_ps(hi.data)) * r).data);
            auto hiG = _mm256_cvtps_epi32((float8(_mm256_cvtepi32_ps(hi.data)) * g).data);
            auto hiB = _mm256_cvtps_epi32((float8(_mm256_cvtepi32_ps(hi.data)) * b).data);
            hi = _mm256_or_si256(hiR, _mm256_or_si256(_mm256_slli_epi32(hiG, 8), _mm256_slli_epi32(hiB, 16)));
            unpack4lo(lo, hi).storeu((int*)out);
            unpack4hi(lo, hi).storeu((int*)out + 8);
            local += 16;
            out += 16;
        } while (--i);
        local += TILE_STRIDE - width;
        out += stride;
    } while (--height);
}

DEBUG_NOINLINE static void ResolveCoverageTile(
    unsigned* out, ptrdiff_t stride, const short* local, int width, int height) {
    // TODO: remove me after fixing the AVX version to handle boarder regions
    for (size_t j = 0; j < height; j++)
        for (size_t i = 0; i < width; i++) {
            int x = local[j * TILE_STRIDE + i];
            x = abs(x);
            x = (x - (x >> 7)) >> 4;
            if (x > 0xff) x = 0xff;
            out[j * stride + i] = 0x10101 * x;
        }
}

__forceinline void DrawContextImpl::Trace(const ThreadState& state,
                                          const Mesh* scene,
                                          const Box4Node* node,
                                          PackedMatrix2x3 objectFromWorld,
                                          PackedMatrix2x3 tileBounds) {
    float8 mins, maxs;
    tileBounds.ComputeBounds(mins, maxs);
    Trace(state, scene, node, objectFromWorld, tileBounds, mins, maxs);
}

void DrawContextImpl::Trace(const ThreadState& state,
                            const Mesh* scene,
                            const Box4Node* node,
                            PackedMatrix2x3 objectFromWorld,
                            PackedMatrix2x3 tileBounds,
                            float8 mins,
                            float8 maxs) {
    auto mask = movemask((mins <= float8::Load(node->max_x)) & (float8::Load(node->min_x) <= maxs));
    mask = mask & (mask >> 4);

    for (unsigned i = 0; i < 4; i++) {
        // If we miss the box, continue.
        if (!(mask & (1 << i))) continue;

        // If the box is not a leaf traverse it.
        if (!node->IsLeaf(i)) {
            Trace(state, scene, node->Child(i), objectFromWorld, tileBounds, mins, maxs);
            continue;
        }

        // The box contains leaf geometry or references to other BVHs.
        auto& tile = *state.tile;
        auto worldFromObject = ShapeTransform(invert(*(const Matrix2x3*)&objectFromWorld));

        // Iterate through the leaf geometry/BVH references.
        for (auto p = node->ShapesAt(scene->shapes, i), e = node->ShapesAt(scene->shapes, i + 1); p < e; p++) {
            // If we hit a reference, update our transform and continue traversal.
            if (p->IsBox4NodeRef()) {
                auto ref = node->Box4NodeRefAt(scene->shapes, i);
                auto childFromObject = PackedMatrix2x3(ref->objectFromParent);
                Trace(state, scene + ref->meshIndex, scene[ref->meshIndex].nodes + ref->nodeIndex,
                      childFromObject * objectFromWorld, childFromObject * tileBounds);
                continue;
            }

            // The leaf geometry is a shape, draw it.
            auto shape = worldFromObject * *p;
            if (p->IsTrapazoid())
                DrawTrapezoid(state.coverage + tile.tileAddressOffset, shape, tile.pixelport);
            else
                DrawCurve(state.coverage + tile.tileAddressOffset, shape, tile.pixelport);
        }
    }
}

void DrawContextImpl::DrawTile(const ThreadState& state,
                               const Mesh* scene,
                               PackedMatrix2x3 layoutFromScreen) {
    auto coverage = state.coverage;
    auto& tile = *state.tile;

    // Compute the bounds of the beam to be traversed.
    auto pixelport = tile.pixelport + sampleBounds;
    auto tileBounds = PackedMatrix2x3(worldFromScreen) * PackedMatrix2x3(layoutFromScreen) *
        PackedMatrix2x3(pixelport.data[0] + pixelport.data[2], 0, -pixelport.data[0],
            0, pixelport.data[1] + pixelport.data[3], -pixelport.data[1]);

    ClearCoverageTile(coverage);
    Trace(state, scene, scene->nodes, layoutFromScreen, tileBounds);

    if (tile.tileWidth % 16)
        ResolveCoverageTile(renderTarget.data + tile.renderTargetOffset, renderTarget.stride, coverage,
                            tile.tileWidth, tile.tileHeight);
    else
        ResolveCoverageTileAVX(renderTarget.data + tile.renderTargetOffset, renderTarget.stride, coverage,
                               tile.tileWidth, tile.tileHeight, state.red, state.green, state.blue);
}

void DrawContextImpl::Draw(const Mesh* scene, const Matrix2x3& screenFromLayout) {
    if (!scene) return;
    auto layoutFromScreen = PackedMatrix2x3(invert(screenFromLayout));
    auto& state = threadStates[0];
    state.red = float8(0.4f);
    state.green = float8(1.0f);
    state.blue = float8(0.6f);
    for (size_t i = 0; i < tiles.size; i++) {
        state.tile = tiles + i;
        DrawTile(state, scene, layoutFromScreen);
    }
}