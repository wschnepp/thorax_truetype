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

//==============================================================================
// FontRenderInfo
//==============================================================================

struct FontRenderInfo {
    bool Initialize(const Font* font);
    void Clear();
    size_t LayoutGlyphs(Box4NodeRef* glyphRefs, int meshIndex, const unsigned* codepoints, size_t numCodepoints);

    struct SimpleGlyphInfo {
        unsigned firstShape;
        unsigned entryNode;
        float advanceWidth;
    };

    struct CompoundGlyphInfo {
        unsigned firstElement;
        float advanceWidth;
    };

    struct CompoundElement {
        Matrix2x3 transform;
        unsigned glyphID;
    };

    // codepoint->index buffer
    // Negative codepoints indicate compound glyphs, the index i for i < 0 is compoundGlyphInfos[~i].
    int* codepointIndex = nullptr;                   

    // Glyph information buffers.
    Array<SimpleGlyphInfo> simpleGlyphInfos;     // 
    Array<CompoundGlyphInfo> compoundGlyphInfos; // 
    Array<CompoundElement> compoundElements;     //
    Array<Box4Node> nodes;                       // List of all BVH nodes used by all glyphs.
    Array<Shape> shapes;                         // List of all shapes used by all glyphs.

    // All values in master grids units.
    float emsPerUnit = 0;
    float ascent = 0;         // Distance from baseline of highest ascender.
    float descent = 0;        // Distance from baseline of lowest descender.
    float lineGap = 0;        // Typographic line gap.
    float lineSpacing = 0;    // Total advance height from line to line.
    float caretSlopeRise = 0; // Used to calculate the slope of the caret (rise/run) set to 1 for vertical caret.
    float caretSlopeRun = 0;  // 0 for vertical
    float caretOffset = 0;    // set value to 0 for non-slanted fonts
};
