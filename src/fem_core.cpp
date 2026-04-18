#include "fem_core.h"
#include "fem_types.h"
#include "../barsystem.h"
#include <cstring>
#include <iostream>

namespace FEM {

Model::Model()
    : nTotalNode(0), nConstrainedNode(0), nTotalElem(0),
      nMaterialType(0), nSectionType(0), nLoad(0),
      nTotalDOF(0), nFreeDOF(0),
      pNode(nullptr), pConsNode(nullptr), pElem(nullptr),
      pMate(nullptr), pSect(nullptr), pLoad(nullptr),
      pDisp(nullptr), pLoadVect(nullptr), pGK(nullptr),
      pDiag(nullptr), pElemDOF(nullptr)
{
}

Model::~Model() {
    cleanup();
}

void Model::cleanup() {
    delete[] pNode;
    delete[] pConsNode;
    delete[] pElem;
    delete[] pMate;
    delete[] pSect;
    delete[] pLoad;
    delete[] pDisp;
    delete[] pLoadVect;
    delete[] pGK;
    delete[] pDiag;

    if (pElemDOF) {
        TwoArrayFree(nTotalElem, pElemDOF);
        pElemDOF = nullptr;
    }

    pNode = nullptr;
    pConsNode = nullptr;
    pElem = nullptr;
    pMate = nullptr;
    pSect = nullptr;
    pLoad = nullptr;
    pDisp = nullptr;
    pLoadVect = nullptr;
    pGK = nullptr;
    pDiag = nullptr;
}

int Model::addNode(int type, double x, double y) {
    // Reallocate node array
    Node* newNodes = new Node[nTotalNode + 1];
    if (pNode) {
        std::memcpy(newNodes, pNode, sizeof(Node) * nTotalNode);
        delete[] pNode;
    }
    pNode = newNodes;

    pNode[nTotalNode].iType = type;
    pNode[nTotalNode].dX = x;
    pNode[nTotalNode].dY = y;
    pNode[nTotalNode].iaDOFIndex[0] = -1;
    pNode[nTotalNode].iaDOFIndex[1] = -1;
    pNode[nTotalNode].iaDOFIndex[2] = -1;

    return nTotalNode++;
}

int Model::addElement(int type, int node1, int node2, int section, int material) {
    if (node1 < 0 || node1 >= nTotalNode || node2 < 0 || node2 >= nTotalNode) {
        return -1;
    }

    Element* newElems = new Element[nTotalElem + 1];
    if (pElem) {
        std::memcpy(newElems, pElem, sizeof(Element) * nTotalElem);
        delete[] pElem;
    }
    pElem = newElems;

    pElem[nTotalElem].iType = type;
    pElem[nTotalElem].iaNode[0] = node1;
    pElem[nTotalElem].iaNode[1] = node2;
    pElem[nTotalElem].iSection = section;
    pElem[nTotalElem].iMaterial = material;
    pElem[nTotalElem].dLength = 0.0;
    pElem[nTotalElem].dSin = 0.0;
    pElem[nTotalElem].dCos = 1.0;
    for (int i = 0; i < 6; i++) {
        pElem[nTotalElem].daEndInterForce[i] = 0.0;
    }

    return nTotalElem++;
}

int Model::addMaterial(double E, double mu, double alpha) {
    Material* newMats = new Material[nMaterialType + 1];
    if (pMate) {
        std::memcpy(newMats, pMate, sizeof(Material) * nMaterialType);
        delete[] pMate;
    }
    pMate = newMats;

    pMate[nMaterialType].dE = E;
    pMate[nMaterialType].dMu = mu;
    pMate[nMaterialType].dAlpha = alpha;

    return nMaterialType++;
}

int Model::addSection(double A, double Iz, double H) {
    Section* newSects = new Section[nSectionType + 1];
    if (pSect) {
        std::memcpy(newSects, pSect, sizeof(Section) * nSectionType);
        delete[] pSect;
    }
    pSect = newSects;

    pSect[nSectionType].dA = A;
    pSect[nSectionType].dIz = Iz;
    pSect[nSectionType].dH = H;

    return nSectionType++;
}

int Model::addLoad(int type, int direction, double value, int elem, int node,
                   double position, double T0, double T1) {
    Load* newLoads = new Load[nLoad + 1];
    if (pLoad) {
        std::memcpy(newLoads, pLoad, sizeof(Load) * nLoad);
        delete[] pLoad;
    }
    pLoad = newLoads;

    pLoad[nLoad].iType = type;
    pLoad[nLoad].iDirect = direction;
    pLoad[nLoad].dValue = value;
    pLoad[nLoad].iLoadedElem = elem;
    pLoad[nLoad].iLoadedNode = node;
    pLoad[nLoad].dPosition = position;
    pLoad[nLoad].dT0 = T0;
    pLoad[nLoad].dT1 = T1;

    return nLoad++;
}

int Model::addConstraint(int node, int dofX, int dofY, int dofR) {
    if (node < 0 || node >= nTotalNode) {
        return -1;
    }

    ConstrainedNode* newCons = new ConstrainedNode[nConstrainedNode + 1];
    if (pConsNode) {
        std::memcpy(newCons, pConsNode, sizeof(ConstrainedNode) * nConstrainedNode);
        delete[] pConsNode;
    }
    pConsNode = newCons;

    pConsNode[nConstrainedNode].iNode = node;
    pConsNode[nConstrainedNode].iaConstrainedDOF[0] = dofX;
    pConsNode[nConstrainedNode].iaConstrainedDOF[1] = dofY;
    pConsNode[nConstrainedNode].iaConstrainedDOF[2] = dofR;

    return nConstrainedNode++;
}

bool Model::loadFromFile(const std::string& filename) {
    // Implementation similar to existing main() function
    std::ifstream fin(filename.c_str());
    if (!fin) return false;

    int nNodes, nCons, nElems, nMats, nSects, nLoads;
    fin >> nNodes >> nCons >> nElems >> nMats >> nSects >> nLoads;

    if (fin.fail() || nNodes <= 0 || nElems <= 0 || nMats <= 0 || nSects <= 0 || nLoads < 0) {
        return false;
    }

    cleanup();

    nTotalNode = nNodes;
    nConstrainedNode = nCons;
    nTotalElem = nElems;
    nMaterialType = nMats;
    nSectionType = nSects;
    nLoad = nLoads;

    pNode = new Node[nTotalNode];
    pConsNode = new ConstrainedNode[nConstrainedNode];
    pElem = new Element[nTotalElem];
    pMate = new Material[nMaterialType];
    pSect = new Section[nSectionType];
    pLoad = new Load[nLoad];

    // Read nodes
    for (int i = 0; i < nTotalNode; i++) {
        fin >> pNode[i].iType >> pNode[i].dX >> pNode[i].dY;
    }

    // Read constraints
    for (int i = 0; i < nConstrainedNode; i++) {
        fin >> pConsNode[i].iNode
            >> pConsNode[i].iaConstrainedDOF[0]
            >> pConsNode[i].iaConstrainedDOF[1]
            >> pConsNode[i].iaConstrainedDOF[2];
    }

    // Read elements
    for (int i = 0; i < nTotalElem; i++) {
        fin >> pElem[i].iType >> pElem[i].iaNode[0] >> pElem[i].iaNode[1] >> pElem[i].iSection;
        std::string rest;
        std::getline(fin, rest);
        std::istringstream iss(rest);
        int matIndex = 0;
        if (!(iss >> matIndex)) matIndex = 0;
        pElem[i].iMaterial = matIndex;
    }

    // Read materials
    for (int i = 0; i < nMaterialType; i++) {
        fin >> pMate[i].dE >> pMate[i].dMu >> pMate[i].dAlpha;
    }

    // Read sections
    for (int i = 0; i < nSectionType; i++) {
        fin >> pSect[i].dA >> pSect[i].dIz >> pSect[i].dH;
    }

    // Read loads
    for (int i = 0; i < nLoad; i++) {
        fin >> pLoad[i].iType >> pLoad[i].iDirect >> pLoad[i].dValue
            >> pLoad[i].iLoadedElem >> pLoad[i].iLoadedNode >> pLoad[i].dPosition
            >> pLoad[i].dT0 >> pLoad[i].dT1;
    }

    return true;
}

bool Model::solve() {
    if (nTotalNode <= 0 || nTotalElem <= 0) {
        return false;
    }

    // Calculate element geometry
    LengthSinCosCalcu(nTotalElem, pElem, pNode);

    // Calculate DOF indices
    int iBuf = 0;
    nFreeDOF = DOFIndexCalcu(iBuf, nTotalNode, nConstrainedNode, pConsNode, pNode);
    nTotalDOF = iBuf;

    if (nFreeDOF <= 0) {
        return false;
    }

    // Allocate solution arrays
    pElemDOF = TwoArrayIntAlloc(nTotalElem, 6);
    pDiag = new int[nTotalDOF];
    pDisp = new double[nTotalDOF];
    pLoadVect = new double[nTotalDOF];

    VectorZeroize(nTotalDOF, pDisp);
    VectorZeroize(nTotalDOF, pLoadVect);

    // Calculate element DOFs
    ElementDOFCalcu(nTotalElem, pNode, pElem, pElemDOF);

    // Calculate bandwidth and diagonal addresses
    BandAndDiagCalcu(nTotalElem, nTotalDOF, pElem, pElemDOF, pDiag);

    // Allocate global stiffness matrix
    pGK = new double[pDiag[nTotalDOF - 1] + 1];
    VectorZeroize(pDiag[nTotalDOF - 1] + 1, pGK);

    // Assemble global stiffness matrix
    GKAssembly(nTotalDOF, nTotalElem, pElem, pNode, pMate, pSect, pDiag, pGK, "");

    // Assemble load vector
    LoadVectorAssembly(nLoad, nTotalDOF, nFreeDOF, pDiag, pGK, pElem, pMate, pSect,
                      pLoad, pNode, pLoadVect, pDisp);

    // Solve system
    if (!LDLTSolve(nFreeDOF, pDiag, pGK, pLoadVect)) {
        return false;
    }

    // Copy solution to displacement vector
    for (int i = 0; i < nFreeDOF; i++) {
        pDisp[i] = pLoadVect[i];
    }

    // Initialize element end forces
    ElementEndForceInit(nTotalElem, pElem);

    // Calculate internal forces
    InternalForceCalcu(nTotalElem, pElem, pNode, pMate, pSect, pDisp);

    // Calculate support reactions
    SupportReactionCalcu(nTotalDOF, nFreeDOF, pDiag, pGK, pDisp, pLoadVect);

    return true;
}

double* Model::getDisplacements(int& count) {
    if (!pDisp) {
        count = 0;
        return nullptr;
    }
    count = nTotalDOF;
    return pDisp;
}

double* Model::getForces(int& count) {
    if (!pElem) {
        count = 0;
        return nullptr;
    }
    count = nTotalElem * 6;
    static std::vector<double> forces;
    forces.resize(count);
    for (int i = 0; i < nTotalElem; i++) {
        for (int j = 0; j < 6; j++) {
            forces[i * 6 + j] = pElem[i].daEndInterForce[j];
        }
    }
    return forces.data();
}

double* Model::getReactions(int& count) {
    if (!pLoadVect || !pConsNode) {
        count = 0;
        return nullptr;
    }
    count = nConstrainedNode * 3;
    static std::vector<double> reactions;
    reactions.resize(count);
    for (int i = 0; i < nConstrainedNode; i++) {
        int node = pConsNode[i].iNode;
        for (int j = 0; j < 3; j++) {
            int dof = pNode[node].iaDOFIndex[j];
            if (dof >= 0) {
                reactions[i * 3 + j] = pLoadVect[dof];
            } else {
                reactions[i * 3 + j] = 0.0;
            }
        }
    }
    return reactions.data();
}

bool Model::exportResults(const std::string& filename) {
    std::ofstream fout(filename.c_str());
    if (!fout) return false;

    // Output model data
    fout << "Total nodes: " << nTotalNode << std::endl;
    fout << "Constrained nodes: " << nConstrainedNode << std::endl;
    fout << "Total elements: " << nTotalElem << std::endl;
    fout << "Material types: " << nMaterialType << std::endl;
    fout << "Section types: " << nSectionType << std::endl;
    fout << "Loads: " << nLoad << std::endl;
    fout << std::endl;

    // Output nodes
    NodeDisplOutput(fout, nTotalNode, pNode, pDisp);

    // Output element forces
    EndInternalForceOutput(fout, nTotalElem, pElem);

    // Output reactions
    SupportReactionOutput(fout, nConstrainedNode, pConsNode, pNode, pLoadVect);

    return true;
}

} // namespace FEM
