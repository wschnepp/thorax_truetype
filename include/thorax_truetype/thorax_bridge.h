/**
 * Copyright (c) 2018-present, Apollo Ellis
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree. An additional grant
 * of patent rights can be found in the PATENTS file in the same directory.
 */

// modified 12.10.2021 (Wilhelm Schnepp): Manually merge Appollo Ellis' changes into fresh copy of Thorax True Type


#pragma once

#include <string>

#include "grids.h"
#include "renderer/drawcontext.h"

Grid::TextGrid* GetTextGridCPU(Array<Shape> &shapes, Array<Grid::GlyphGrid> &glyphGrids, Array<Grid::GlyphGridCell> &gridCells,
	Array<Grid::shape_ptr> &shapePtrs);

Grid::TextGrid* InitializeText(std::string path, std::string filename);