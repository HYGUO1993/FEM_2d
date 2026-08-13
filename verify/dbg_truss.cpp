#include <iostream>
#include <cmath>
#include "barsystem.h"
int main() {
    Node nodes[2];
    nodes[0].iType = TRUSS_NODE; nodes[0].dX = 0.0; nodes[0].dY = 0.0;
    nodes[1].iType = TRUSS_NODE; nodes[1].dX = 1.0; nodes[1].dY = 1.0;
    nodes[0].iaDOFIndex[0]=2; nodes[0].iaDOFIndex[1]=3; nodes[0].iaDOFIndex[2]=0;
    nodes[1].iaDOFIndex[0]=0; nodes[1].iaDOFIndex[1]=1; nodes[1].iaDOFIndex[2]=0;
    Element elems[1];
    elems[0].iType = TRUSS; elems[0].iaNode[0]=0; elems[0].iaNode[1]=1;
    elems[0].iSection=0; elems[0].iMaterial=0;
    Material mate[1]; mate[0].dE=210e9; mate[0].dMu=0.3; mate[0].dAlpha=0;
    Section sect[1]; sect[0].dA=0.01; sect[0].dIz=0; sect[0].dH=0;
    LengthSinCosCalcu(1, elems, nodes);
    std::cout << "L=" << elems[0].dLength << " cos=" << elems[0].dCos << " sin=" << elems[0].dSin << "\n";
    int pDiag[4]={0,0,0,0};
    int** edof = TwoArrayIntAlloc(1,6);
    ElementDOFCalcu(1, nodes, elems, edof);
    BandAndDiagCalcu(1, 4, elems, edof, pDiag);
    std::cout << "pDiag: " << pDiag[0] << "," << pDiag[1] << "," << pDiag[2] << "," << pDiag[3] << "\n";
    int ngk = pDiag[3]+1;
    double* pGK = new double[ngk];
    for (int i=0;i<ngk;i++) pGK[i]=0;
    GKAssembly(4, 1, elems, nodes, mate, sect, pDiag, pGK, "");
    double* pLoad = new double[4]{-707.10678, -707.10678, 0, 0};
    double* pDisp = new double[4]{0,0,0,0};
    std::cout << "pGK:";
    for (int i=0;i<ngk;i++) std::cout << " " << pGK[i];
    std::cout << "\n";
    bool ok = LDLTSolve(4, pDiag, pGK, pLoad);
    std::cout << "LDLT ok=" << ok << " u:";
    for (int i=0;i<4;i++) { pDisp[i]=pLoad[i]; std::cout << " " << pDisp[i]; }
    std::cout << "\n";
    InternalForceCalcu(1, elems, nodes, mate, sect, pDisp);
    std::cout << "endf:";
    for (int k=0;k<6;k++) std::cout << " " << elems[0].daEndInterForce[k];
    std::cout << "\n";
    TwoArrayFree(1, edof);
    delete[] pGK; delete[] pLoad; delete[] pDisp;
    return 0;
}
