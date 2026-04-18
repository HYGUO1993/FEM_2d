#pragma once

namespace FEM {

// Element types
#define TRUSS 1
#define FRAME 2

// Node types
#define TRUSS_NODE 1
#define FRAME_NODE 2

// Load types
#define FORCE_ON_NODE 1
#define LATERAL_FORCE 2
#define LATERAL_UNIFORM_PRESSURE 3
#define MOMENT_ON_A_POINT 4
#define LATERAL_LINEARLY_PRESSURE 5
#define AXIAL_PRESSURE 6
#define AXIAL_FORCE 7
#define MOMENT_ON_BEAM 8
#define TEMPERATURE 9
#define SUPPORT_MOVE 10

// Direction constants
#define DIRECT_X 0
#define DIRECT_Y 1
#define DIRECT_R 2

/**
 * @brief Material properties
 */
struct Material {
    double dE;       // Young's modulus
    double dMu;      // Poisson's ratio
    double dAlpha;   // Thermal expansion coefficient
};

/**
 * @brief Cross-section properties
 */
struct Section {
    double dA;       // Cross-sectional area
    double dIz;      // Moment of inertia
    double dH;       // Section height
};

/**
 * @brief Node structure
 */
struct Node {
    int iType;           // Node type (TRUSS_NODE or FRAME_NODE)
    double dX, dY;       // Coordinates
    int iaDOFIndex[3];   // DOF indices [Ux, Uy, Rz]
};

/**
 * @brief Element structure
 */
struct Element {
    int iType;              // Element type (TRUSS or FRAME)
    int iaNode[2];          // Node indices
    int iSection;           // Section index
    int iMaterial;          // Material index
    double dLength;         // Element length
    double dSin, dCos;      // Direction cosines
    double daEndInterForce[6];  // End forces
};

/**
 * @brief Load structure
 */
struct Load {
    int iType;          // Load type
    int iDirect;        // Direction (X, Y, or R)
    double dValue;      // Load value
    int iLoadedElem;    // Element index (-1 if node load)
    int iLoadedNode;    // Node index (-1 if element load)
    double dPosition;   // Position along element
    double dT0, dT1;    // Temperature values
};

/**
 * @brief Constrained node structure
 */
struct ConstrainedNode {
    int iNode;              // Node index
    int iaConstrainedDOF[3];  // Constraint flags for [Ux, Uy, Rz]
                              // -1 = constrained, 0 = free, >10000 = coupled
};

} // namespace FEM
