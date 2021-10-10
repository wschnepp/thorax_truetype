//===- ShapeBuilder.cpp -------------------------------------------*- C++ --*-===//
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

#include <thorax_truetype/renderer/ShapeBuilder.h>
#include <thorax_truetype/renderer/types.h>

#include <cassert>
#include <cmath>
#include <cstdio>
#include <utility>

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
// ShapeBuilderImpl
//==============================================================================

struct ShapeBuilderImpl final : ShapeBuilder {
	struct EdgeEqn {
		EdgeEqn(short2 vert, short2 next) : a((float)(next.x - vert.x) / (next.y - vert.y)), b(vert.x - a * vert.y) {}
		float XofY(float Y) const { return a * Y + b; }

		float a;
		float b;
	};

	struct Event {
		Event() = default;
		Event(unsigned short segmentID, bool isBegin, bool isUp, const short2& vert)
		    : segmentID(segmentID), isBegin(isBegin), isUp(isUp), vert(vert) {}
		bool operator<(const Event& a) const { return vert.y < a.vert.y; }
		unsigned short segmentID;
		bool isBegin;
		bool isUp;
		short2 vert;
	};

	struct Edge {
		Edge() = default;
		Edge(float x, unsigned short segmentID, unsigned short slabID) : segmentID(segmentID), slabID(slabID), x(x) {}
		bool operator<(const Edge& a) const { return x < a.x; }
		float x;
		unsigned short segmentID;
		unsigned short slabID;
	};

	struct Slab {
		Slab() = default;
		Slab(float y_min, unsigned short upSegmentID, unsigned short dnSegmentID)
		    : y_min(y_min), upSegmentID(upSegmentID), dnSegmentID(dnSegmentID) {}
		float y_min;
		unsigned short upSegmentID;
		unsigned short dnSegmentID;
	};

	void Clear(size_t reserve) override;
	size_t GenerateShapes(const Segment* segments,
	                      const size_t* contourSizes,
	                      size_t numCountours,
	                      Shape* shapes) override;

	void ProcessSegment(short2 p0, short2 p1, short2 p2);

	void PushTrapezoid(const Slab& slab, float y_max) {
		auto& up = eqns[slab.upSegmentID];
		auto& dn = eqns[slab.dnSegmentID];
		auto y0 = slab.y_min;
		auto y1 = y_max;
		auto x0 = up.XofY(y0);
		auto x1 = up.XofY(y1);
		auto x2 = dn.XofY(y1);
		auto x3 = dn.XofY(y0);
		if (shapes) shapes[numShapes] = Shape(x0, x1, x2, x3, y0, y1);
		numShapes++;
	}

	DynamicArray<EdgeEqn> eqns;
	DynamicArray<Event> events;
	DynamicArray<Edge> upEdges;
	DynamicArray<Edge> dnEdges;
	DynamicArray<Slab> slabs;
	Shape* shapes;
	size_t numShapes;
};

void ShapeBuilderImpl::Clear(size_t reserve) {
	eqns = DynamicArray<EdgeEqn>(reserve);
	events = DynamicArray<Event>(reserve * 2);
	upEdges = DynamicArray<Edge>(reserve);
	dnEdges = DynamicArray<Edge>(reserve);
	slabs = DynamicArray<Slab>(reserve * 2);
}

void ShapeBuilderImpl::ProcessSegment(short2 p0, short2 p1, short2 p2) {
	// Create the curves and generate the trapezoid event list.
	if (p0.y == p2.y) return; // Ignore flat horizontal segments.

	auto v1_x = p1.x - p0.x;
	auto v1_y = p1.y - p0.y;
	auto v2_x = p2.x - p0.x;
	auto v2_y = p2.y - p0.y;
	// If we have a curved segment, push it.
	if (v2_x * v1_y != v1_x * v2_y) {
		if (shapes)
			shapes[numShapes] =
			    Shape(Point(p0.x, p0.y), Vector((float)v1_x, (float)v1_y), Vector((float)v2_x, (float)v2_y));
		numShapes++;
	}

	auto segmentID = (unsigned short)eqns.size;
	if (p0.y < p2.y) {
		events.Push(segmentID, true, true, p0);
		events.Push(segmentID, false, true, p2);
	} else {
		events.Push(segmentID, false, false, p0);
		events.Push(segmentID, true, false, p2);
	}
	eqns.Push(p0, p2);
}

static void InsertEdge(DynamicArray<ShapeBuilderImpl::Edge>& edges,
                       float x,
                       unsigned short segmentID,
                       unsigned short slabID) {
	size_t i = edges.size;
	while (i && edges[i - 1].x < x) i--;
	for (auto j = edges.size++; j > i; j--) edges[j] = edges[j - 1];
	edges[i] = ShapeBuilderImpl::Edge(x, segmentID, slabID);
};

static unsigned short DeleteEdge(DynamicArray<ShapeBuilderImpl::Edge>& edges, unsigned short segmentID) {
	size_t i = 0;
	while (i < edges.size && edges[i].segmentID != segmentID) i++;
	auto slabID = edges[i].slabID;
	for (auto e = edges.size - 1; i < e; i++) edges[i] = edges[i + 1];
	edges.size--;
	return slabID;
};

size_t ShapeBuilderImpl::GenerateShapes(const Segment* segments,
                                        const size_t* contourSizes,
                                        size_t numCountours,
                                        Shape* out) {
	// Set up the outputs.
	shapes = out;
	numShapes = 0;

	// Reset the event and equations lists and repopulate them.  Process contour also
	events.Resize(0);
	eqns.Resize(0);
	for (auto pSize = contourSizes, end = contourSizes + numCountours; pSize < end; pSize++) {
		if (*pSize <= 1) continue;
		for (size_t i = 0, e = *pSize - 1; i < e; i++)
			ProcessSegment(segments[i].vert, segments[i].knot, segments[i + 1].vert);
		ProcessSegment(segments[*pSize - 1].vert, segments[*pSize - 1].knot, segments[0].vert);
		segments += *pSize;
	}
	if (!events.size) return numShapes;

	// Shell Sort the event list.
	static constexpr size_t gapSequence[] = {57, 23, 10, 4, 1};
	for (auto gap : gapSequence)
		for (auto i = gap, e = events.size; i < e; ++i)
			for (auto j = i - gap; j < e && events[j + gap] < events[j]; j -= gap)
				std::swap(events[j], events[j + gap]);

	// Walk the event list and create the trapezoid list.
	upEdges.Resize(0);
	dnEdges.Resize(0);
	slabs.Resize(0);
	auto event = events.begin();
	do {
		auto y = event->vert.y;
		for (auto& edge : upEdges) edge.x = eqns[edge.segmentID].XofY(y);
		for (auto& edge : dnEdges) edge.x = eqns[edge.segmentID].XofY(y);

		// Iterate over all of the events that occur at this y.
		auto baseSlabID = slabs.size;
		do {
			if (event->isUp) {
				if (event->isBegin) {
					auto slabID = slabs.size;
					slabs.Push(y, event->segmentID, (unsigned short)0);
					InsertEdge(upEdges, event->vert.x, event->segmentID, (unsigned short)slabID);
				} else {
					auto slabID = DeleteEdge(upEdges, event->segmentID);
					PushTrapezoid(slabs[slabID], y);
				}
			} else {
				if (event->isBegin)
					InsertEdge(dnEdges, event->vert.x, event->segmentID, (unsigned short)-1);
				else
					DeleteEdge(dnEdges, event->segmentID);
			}
			++event;
		} while (event != events.end() && event->vert.y == y);

		// Re-pair edges and adopt orphaned edges.
		auto numPairs = upEdges.size;
		for (size_t i = 0; i < numPairs; i++) {
			if (upEdges[i].slabID == dnEdges[i].slabID) continue;
			if (upEdges[i].slabID < baseSlabID) {
				PushTrapezoid(slabs[upEdges[i].slabID], y); // We split an old slab.
				auto slabID = slabs.size;
				slabs.Push(y, upEdges[i].segmentID, dnEdges[i].segmentID);
				upEdges[i].slabID = dnEdges[i].slabID = (unsigned short)slabID;
			} else {
				dnEdges[i].slabID = upEdges[i].slabID;
				slabs[dnEdges[i].slabID].dnSegmentID = dnEdges[i].segmentID;
			}
		}
		assert(numPairs == dnEdges.size);
	} while (event != events.end());
	return numShapes;
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
		static const __declspec(align(32)) int MASK[] = {-1, -1, -1, -1, -1, -1, 0, 0};
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