#pragma once

// Include barsystem.h to get the actual structure definitions
#include "../barsystem.h"

namespace FEM {

// Re-use the global type definitions from barsystem.h
// by creating type aliases within the FEM namespace
using ::Material;
using ::Section;
using ::Node;
using ::Element;
using ::Load;
using ::ConstrainedNode;

} // namespace FEM
