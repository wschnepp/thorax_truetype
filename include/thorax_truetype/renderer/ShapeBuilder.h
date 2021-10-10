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

#pragma once
#include "types.h"
#include <thorax_truetype/thorax_truetype.h>

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