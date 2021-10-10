//===- drawcontext.h -------------------------------------------*- C++ --*-===//
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
// Drawing Context API
//==============================================================================

struct DrawContext {
	virtual size_t SetRenderTarget(unsigned* renderTarget, int width, int height, size_t stride) = 0;
	virtual void SetViewport(float lower_x, float lower_y, float upper_x, float upper_y) = 0;
	virtual void SetFilterKernel(float x0, float y0, float x1, float y1) = 0; //< in screen space (1 == 1 pixel)
	virtual void Draw(const Mesh* scene, const Matrix2x3& screenFromLayout) = 0;
	virtual ~DrawContext() {}

	static DrawContext* Create(unsigned numThreads = 1);
	static void Destroy(DrawContext*& context);
};