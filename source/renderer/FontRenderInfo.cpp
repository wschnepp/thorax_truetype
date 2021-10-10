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

#include <thorax_truetype/renderer/FontRenderInfo.h>
#include <thorax_truetype/renderer/ShapeBuilder.h>
#include <thorax_truetype/renderer/build.h>

bool FontRenderInfo::Initialize(const Font* font) {
	Clear();
	FontInfo fontInfo;
	if (!fontInfo.Initialize(font)) return false;
	if (!fontInfo.maximumProfile->numGlyphs) return false;
	if (!fontInfo.maximumProfile->maxComponentElements) return false;

	// Load font metrics.
	emsPerUnit = 1.0f / (int)fontInfo.fontHeader->unitsPerEm;
	ascent = (float)(int)fontInfo.metricsHeader->ascent;
	descent = (float)(int)fontInfo.metricsHeader->descent;
	lineGap = (float)(int)fontInfo.metricsHeader->lineGap;
	caretSlopeRise = (float)(int)fontInfo.metricsHeader->caretSlopeRise;
	caretSlopeRun = (float)(int)fontInfo.metricsHeader->caretSlopeRun;
	caretOffset = (float)(int)fontInfo.metricsHeader->caretOffset;
	lineSpacing = ascent - descent + lineGap;

	size_t numGlyphs = fontInfo.maximumProfile->numGlyphs;

	// Calculate the codes for the sorted glyphs.
	struct GlyphCacheEntry {
		GlyphCacheEntry() = default;
		GlyphCacheEntry(const GlyphData* glyphData, int index) : glyphData(glyphData), index(index) {}
		const GlyphData* glyphData;
		int index;
	};
	DynamicArray<GlyphCacheEntry> cache(numGlyphs);
	size_t numSimpleGlyphs = 0;
	size_t numCompoundGlyphs = 0;
	for (unsigned i = 0, e = (unsigned)numGlyphs; i < e; i++) {
		auto glyphData = fontInfo.GetGlyphDataByIndex(i);
		if (glyphData && glyphData->numberOfContours < 0)
			cache.Push(glyphData, ~(int)numCompoundGlyphs++);
		else
			cache.Push(glyphData, (int)numSimpleGlyphs++);
	}

	// Create the codepointIndex table.
	codepointIndex = new int[0x10000];
	for (size_t i = 0; i < 0x10000; i++)
		codepointIndex[i] = cache[fontInfo.GetGlyphIndexByCodepoint((unsigned)i)].index;

	// Shape building buffers.
	DynamicArray<GlyphPoint> glyphPoints = DynamicArray<GlyphPoint>(fontInfo.maximumProfile->maxPoints);
	DynamicArray<Segment> segments = DynamicArray<Segment>(fontInfo.maximumProfile->maxPoints * 3);
	DynamicArray<size_t> contourSizes = DynamicArray<size_t>(fontInfo.maximumProfile->maxContours);
	auto builder = ShapeBuilder::GetDefaultShapeBuilderInstance();
	builder->Clear(fontInfo.maximumProfile->maxPoints * 3);

	auto ProcessGlyph = [&glyphPoints, &segments, &contourSizes, &builder](const GlyphData* glyphData, Shape* shapes) {
		// Define a helper lambda for processing each contour.
		auto ProcessContour = [&segments, &contourSizes](const GlyphPoint* first, size_t numPoints) {
			auto GetValue = [](const GlyphPoint* a) { return short2{a->value_x, a->value_y}; };
			auto Average = [](const GlyphPoint* a, const GlyphPoint* b) {
				return short2{(short)((a->value_x + b->value_x) / 2), (short)((a->value_y + b->value_y) / 2)};
			};

			auto firstSegmentIndex = segments.size;
			auto last = first + numPoints - 1;
			short2 temp; // An on-curve point we haven't yet used in the output.
			if (first->IsOnCurve())
				temp = short2{first->value_x, first->value_y};
			else if (!last->IsOnCurve())
				segments.Push(Average(last, first), GetValue(first));
			for (auto point = first + 1; point <= last; point++) {
				if (point->IsOnCurve()) {
					if (point[-1].IsOnCurve()) segments.Push(temp, temp);
					temp = GetValue(point);
				} else {
					if (!point[-1].IsOnCurve()) temp = Average(point - 1, point);
					segments.Push(temp, GetValue(point));
				}
			}
			if (last->IsOnCurve()) segments.Push(temp, first->IsOnCurve() ? temp : GetValue(first));

			// Calculate the number of segments and push them onto the contour sizes.
			auto numSegments = segments.size - firstSegmentIndex;
			contourSizes.Push(numSegments);
		};

		// Decode the points into the glyphPoints buffer.
		if (!glyphData->DecodeGlyphPoints(glyphPoints.data)) return (size_t)0;

		// Process the contours to generate segments and contour sizes.
		segments.Resize(0);
		contourSizes.Resize(0);
		ProcessContour(glyphPoints.data, (size_t)glyphData->endPtsOfContours[0] + 1);
		for (size_t i = 1, e = (size_t)(short)glyphData->numberOfContours; i < e; ++i)
			ProcessContour(glyphPoints.data + glyphData->endPtsOfContours[i - 1] + 1,
			               glyphData->endPtsOfContours[i] - glyphData->endPtsOfContours[i - 1]);

		return builder->GenerateShapes(segments.data, contourSizes.data, contourSizes.size, shapes);
	};

	// Count all of the shapes and the compound elements.
	// Allocate the glyph info data structures.
	simpleGlyphInfos = Array<SimpleGlyphInfo>(numSimpleGlyphs + 1);
	compoundGlyphInfos = Array<CompoundGlyphInfo>(numCompoundGlyphs + 1);
	unsigned totalShapes = 0;
	unsigned totalCompoundElements = 0;
	for (size_t i = 0; i < numGlyphs; i++) {
		auto index = cache[i].index;
		auto glyphData = cache[i].glyphData;
		auto advanceWidth = (float)(int)fontInfo.GetHorizontalMetricByIndex((unsigned)i).advanceWidth;
		if (index >= 0) {
			auto count = glyphData ? (unsigned)ProcessGlyph(glyphData, nullptr) : 0;
			simpleGlyphInfos[index].advanceWidth = advanceWidth;
			simpleGlyphInfos[index].firstShape = totalShapes;
			totalShapes += count;
		} else {
			auto count = glyphData->DecodeCompoundGlyph(nullptr);
			compoundGlyphInfos[~index].advanceWidth = advanceWidth;
			compoundGlyphInfos[~index].firstElement = totalCompoundElements;
			totalCompoundElements += count;
		}
	}
	simpleGlyphInfos[numSimpleGlyphs].firstShape = totalShapes;
	compoundGlyphInfos[numCompoundGlyphs].firstElement = totalCompoundElements;

	// Allocate storage
	if (!totalShapes) return Clear(), false;

	// Allocate memory based on the number of elements.
	compoundElements = Array<CompoundElement>(totalCompoundElements);
	Array<Shape> tempShapes(totalShapes);
	Array<GlyphReference> references(fontInfo.maximumProfile->maxComponentElements);
	Array<BuilderMesh*> builders(numSimpleGlyphs);

	// Iterate through the glyphs a second time and fill out the now-allocated glyph data.
	auto totalNodeCount = 0u;
	for (size_t i = 0; i < numGlyphs; i++) {
		auto index = cache[i].index;
		auto glyphData = cache[i].glyphData;
		if (index >= 0) {
			auto& glyphInfo = simpleGlyphInfos[index];
			glyphInfo.entryNode = totalNodeCount;
			if (glyphData) {
				auto firstShape = tempShapes + (size_t)glyphInfo.firstShape;
				auto numShapes = ProcessGlyph(glyphData, firstShape);
				builders[index] = new BuilderMesh(firstShape, numShapes);
				totalNodeCount += (unsigned)builders[index]->TranscribeBVH(nullptr, nullptr, 0);
			} else
				builders[index] = nullptr;
		} else {
			// Decode the compound glyph and initialize the element list.
			auto numCompounds = glyphData->DecodeCompoundGlyph(references.data);
			auto elements = compoundElements + (size_t)compoundGlyphInfos[~index].firstElement;
			for (size_t j = 0; j < numCompounds; j++) {
				auto& ref = references[j];
				elements[j].glyphID = cache[ref.glyphIndex].index;
				elements[j].transform = Matrix2x3(ref.x_x, ref.y_x, ref.offset_x, ref.x_y, ref.y_y, ref.offset_y);
			}
		}
	}

	// Allocate final shape and node array and build the BVHs.
	nodes = Array<Box4Node>(totalNodeCount);
	shapes = Array<Shape>(totalShapes);
	for (size_t i = 0; i < numSimpleGlyphs; i++) {
		if (!builders[i]) continue;
		auto glyphInfo = simpleGlyphInfos + i;
		builders[i]->TranscribeBVH(nodes + (size_t)glyphInfo->entryNode, shapes, glyphInfo->firstShape);
		delete builders[i];
	}

	return true;
}

void FontRenderInfo::Clear() {
	codepointIndex = Array<int>();
	simpleGlyphInfos = Array<SimpleGlyphInfo>();
	compoundGlyphInfos = Array<CompoundGlyphInfo>();
	compoundElements = Array<CompoundElement>();
	shapes = Array<Shape>();
}

size_t FontRenderInfo::LayoutGlyphs(Box4NodeRef* glyphRefs,
                                    int meshIndex,
                                    const unsigned* codepoints,
                                    size_t numCodepoints) {
	if (!glyphRefs) {
		// Count the drawn glyphs.
		size_t drawnGlyphCount = 0;
		for (auto j = 0; j < numCodepoints; j++) {
			if (codepoints[j] == '\n') continue;
			auto index = codepointIndex[codepoints[j]];
			if (index >= 0) {
				auto glyphInfo = simpleGlyphInfos + (size_t)index;
				if (glyphInfo->entryNode != glyphInfo[1].entryNode) drawnGlyphCount++;
				continue;
			}
			auto compoundInfo = compoundGlyphInfos + (size_t)~index;
			for (auto i = compoundInfo->firstElement, e = compoundInfo[1].firstElement; i < e; i++) {
				auto glyphInfo = simpleGlyphInfos + (size_t)compoundElements[i].glyphID;
				if (glyphInfo->entryNode != glyphInfo[1].entryNode) drawnGlyphCount++;
			}
		}
		return drawnGlyphCount;
	}

	auto nextGlyphRef = glyphRefs;
	float pen_x = 0, pen_y = 0;
	for (auto j = 0; j < numCodepoints; j++) {
		// Handle newline.
		if (codepoints[j] == '\n') {
			pen_x = 0;
			pen_y -= lineSpacing;
			continue;
		}
		auto index = codepointIndex[codepoints[j]];
		if (index >= 0) {
			// Layout a simple glyph.
			auto glyphInfo = simpleGlyphInfos + (size_t)index;
			if (glyphInfo->entryNode != glyphInfo[1].entryNode) {
				nextGlyphRef->objectFromParent = invert(Matrix2x3(1, 0, pen_x, 0, 1, pen_y));
				nextGlyphRef->meshIndex = meshIndex;
				nextGlyphRef->nodeIndex = glyphInfo->entryNode;
				nextGlyphRef++;
			}
			pen_x += glyphInfo->advanceWidth;
			continue;
		}
		auto compoundInfo = compoundGlyphInfos + (size_t)~index;
		// Layout a compound glyph.
		for (auto i = compoundInfo->firstElement, e = compoundInfo[1].firstElement; i < e; i++) {
			auto& element = compoundElements[i];
			auto glyphInfo = simpleGlyphInfos + (size_t)element.glyphID;
			if (glyphInfo->entryNode != glyphInfo[1].entryNode) {
				nextGlyphRef->objectFromParent = invert(Matrix2x3(1, 0, pen_x, 0, 1, pen_y) * element.transform);
				nextGlyphRef->meshIndex = meshIndex;
				nextGlyphRef->nodeIndex = glyphInfo->entryNode;
				nextGlyphRef++;
			}
		}
		pen_x += compoundInfo->advanceWidth;
	}
	return nextGlyphRef - glyphRefs;
}