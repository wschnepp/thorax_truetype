#pragma once
//===- grid_builder.h -------------------------------------------*- C++ --*-===//
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

// This source code has been partially
// modified by Apollo Ellis Jan 2018 and June 2019.

//modified 12.10.2021 (Wilhelm Schnepp): Extract to separate file

#include <thorax_truetype/grid_types.h>
using namespace Grid;

struct GlyphGridBuilder {
	GlyphGridBuilder(const Shape* shapes, size_t numShapes, size_t firstShape);
	struct ShapeEvent {
		ShapeEvent(int id, float start, float end) : id(id), start(start), end(end) {}
		int id;
		float start;
		float end;
	};

	static bool EventCompare(const ShapeEvent& a, const ShapeEvent& b) { return a.start < b.start; }
	bool null = false;
	size_t hres, vres;
	float lowest_x;
	float lowest_y;
	float highest_x;
	float highest_y;
	int cellCount;
	int totalRefs;
	Grid::shape_ptr* shapePtrs;
	Grid::GridCell* shapeOffsets;
};