/**
 * Copyright (c) 2018-present, Apollo Ellis
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree. An additional grant
 * of patent rights can be found in the PATENTS file in the same directory.
 */

 // modified 12.10.2021 (Wilhelm Schnepp): Manually merge Appollo Ellis' changes into fresh copy of Thorax True Type


#include <iostream>
#include <thorax_truetype/thorax_bridge.h>
#include <thorax_truetype/renderer/build.h>
#include <thorax_truetype/renderer/FontRenderInfo.h>

Grid::TextGrid* GetTextGridCPU(Array<Shape> &shapes, Array<Grid::GlyphGrid> &glyphGrids, Array<Grid::GlyphGridCell> &glyphGridCells,
	Array<Grid::shape_ptr> &shapePtrs){
	
	Grid::GridShape* cpuShapes = new Grid::GridShape[shapes.size];
	for (int i = 0; i < shapes.size; i++) {
		cpuShapes[i].FromShape(shapes[i].p0v1p2v3.data.m256_f32);
	}
	Grid::TextGrid *textGrid = new Grid::TextGrid;
	textGrid->shapes = cpuShapes;
	textGrid->glyph_grids = glyphGrids.data;
	textGrid->glyph_grid_cells = glyphGridCells.data;
	textGrid->shape_ptrs = shapePtrs.data;
	textGrid->glyph_grid_count = (int)glyphGrids.size;
	textGrid->glyph_grid_cell_count = (int)glyphGridCells.size;
	textGrid->shape_count = (int)shapes.size;
	textGrid->shape_ptr_count = (int)shapePtrs.size;
	return textGrid;
}

Grid::TextGrid* InitializeText(std::string path, std::string filename)
{
	std::string textfilename = "Times\ New\ Roman.ttf";
	const Font *font;
	{
		FILE *f = new FILE();
		fopen_s(&f, (path + textfilename).c_str(), "rb");

		fseek(f, 0, SEEK_END);
		int size = ftell(f);
		rewind(f);
		const char* data = new char[size];
		fread_s((void*)data, size, sizeof(char), size, f);
		font = Font::AsFont(data);
		fclose(f);
	}
	FontRenderInfo *ft_render_info = new FontRenderInfo();
	ft_render_info->Initialize(font);

	std::string words = "To be, or not to be, that is the question :";

	unsigned* points = new unsigned[4096];
	int hres = 0;
	int vresmax = 0;
	int hresmax = 0;
	int counts = 0;
	for (int i = 0; i < words.length(); i++)
	{
		if (words[i] == '\n') {
			vresmax++;
			hres = 0;
		}
		else {
			hres++;
			hresmax = hresmax < hres ? hres : hresmax;
		}
		points[i] = words[i];
		counts++;
	}

	int refs = (int)ft_render_info->LayoutGlyphs((Grid::GridRef*)NULL, 0, (const unsigned *)points, counts);
	Grid::GridRef *glyphRefs = new Grid::GridRef[refs];
	int numRefs = ft_render_info->LayoutGlyphs(glyphRefs, 0, (const unsigned *)points, counts);

	Grid::TextGrid *tgrid = GetTextGridCPU(ft_render_info->shapes, ft_render_info->glyph_grids,
		ft_render_info->glyph_grid_cells, ft_render_info->shape_ptrs);
	GridBuild(tgrid, glyphRefs, numRefs, hresmax, vresmax);
	delete[] points;
	return tgrid;
}