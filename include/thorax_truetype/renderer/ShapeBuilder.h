#pragma once 
//===- ShapeBuilder.h -------------------------------------------*- C++ --*-===//
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

#include "types.h"
#include <thorax_truetype/thorax_truetype.h>

#include <memory>

//==============================================================================
// ShapeBuilder API
//==============================================================================

struct ShapeBuilder {
	virtual ~ShapeBuilder() {}
	virtual void Clear(size_t reserve = 0) = 0;
	virtual size_t GenerateShapes(const Segment* segments,
	                              const size_t* contourSizes,
	                              size_t numCountours,
	                              Shape* shapes) = 0;
};

struct ShapeBuilderDefaultImpl final : ShapeBuilder {

	void Clear(size_t reserve) override;
	size_t GenerateShapes(const Segment* segments,
	                      const size_t* contourSizes,
	                      size_t numCountours,
	                      Shape* shapes) override;

  private:
	struct EdgeEqn {
		EdgeEqn(short2 vert, short2 next);
		float XofY(float Y) const;

		float a;
		float b;
	};

	struct Event {
		Event() = default;
		Event(unsigned short segmentID, bool isBegin, bool isUp, const short2& vert);
		bool operator<(const Event& a) const;
		unsigned short segmentID;
		bool isBegin;
		bool isUp;
		short2 vert;
	};

	struct Edge {
		Edge() = default;
		Edge(float x, unsigned short segmentID, unsigned short slabID);
		bool operator<(const Edge& a) const;
		float x;
		unsigned short segmentID;
		unsigned short slabID;
	};

	struct Slab {
		Slab() = default;
		Slab(float y_min, unsigned short upSegmentID, unsigned short dnSegmentID);
		float y_min;
		unsigned short upSegmentID;
		unsigned short dnSegmentID;
	};


	void ProcessSegment(short2 p0, short2 p1, short2 p2);

	void PushTrapezoid(const Slab& slab, float y_max);

	static void InsertEdge(DynamicArray<ShapeBuilderDefaultImpl::Edge>& edges,
	                       float x,
	                       unsigned short segmentID,
	                       unsigned short slabID);
	static unsigned short DeleteEdge(DynamicArray<ShapeBuilderDefaultImpl::Edge>& edges, unsigned short segmentID);

	DynamicArray<EdgeEqn> eqns;
	DynamicArray<Event> events;
	DynamicArray<Edge> upEdges;
	DynamicArray<Edge> dnEdges;
	DynamicArray<Slab> slabs;
	Shape* shapes;
	size_t numShapes;
};