#pragma once

#include <string>
#include <vector>

namespace FEM {

// Forward declarations
struct Node;
struct Element;
struct Material;
struct Section;
struct Load;
struct ConstrainedNode;

/**
 * @brief FEM Model class - encapsulates all model data
 */
class Model {
public:
    Model();
    ~Model();

    // Model building
    int addNode(int type, double x, double y);
    int addElement(int type, int node1, int node2, int section, int material);
    int addMaterial(double E, double mu, double alpha);
    int addSection(double A, double Iz, double H);
    int addLoad(int type, int direction, double value, int elem, int node,
                double position, double T0, double T1);
    int addConstraint(int node, int dofX, int dofY, int dofR);

    // Analysis
    bool solve();

    // Results access
    double* getDisplacements(int& count);
    double* getForces(int& count);
    double* getReactions(int& count);

    // File I/O
    bool loadFromFile(const std::string& filename);
    bool saveToFile(const std::string& filename);
    bool exportResults(const std::string& filename);

    // Getters
    int getNodeCount() const { return nTotalNode; }
    int getElementCount() const { return nTotalElem; }
    int getDOFCount() const { return nTotalDOF; }
    int getFreeDOFCount() const { return nFreeDOF; }

    const Node* getNodes() const { return pNode; }
    const Element* getElements() const { return pElem; }
    const Material* getMaterials() const { return pMate; }
    const Section* getSections() const { return pSect; }
    const Load* getLoads() const { return pLoad; }
    const ConstrainedNode* getConstraints() const { return pConsNode; }

private:
    // Model data
    int nTotalNode;
    int nConstrainedNode;
    int nTotalElem;
    int nMaterialType;
    int nSectionType;
    int nLoad;
    int nTotalDOF;
    int nFreeDOF;

    Node* pNode;
    ConstrainedNode* pConsNode;
    Element* pElem;
    Material* pMate;
    Section* pSect;
    Load* pLoad;

    // Solution data
    double* pDisp;      // Displacement vector
    double* pLoadVect;  // Load vector
    double* pGK;        // Global stiffness matrix (skyline)
    int* pDiag;         // Diagonal addresses
    int** pElemDOF;     // Element DOF mapping

    // Internal methods
    void cleanup();
    bool allocateMemory();
    bool calculateGeometry();
    bool assembleSolver();
};

} // namespace FEM
