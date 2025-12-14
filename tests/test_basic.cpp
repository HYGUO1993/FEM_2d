#include <iostream>
#include <cmath>
#include "barsystem.h"

int main() {
    // Test LengthSinCosCalcu
    Node nodes[2];
    nodes[0].dX = 0.0; nodes[0].dY = 0.0; nodes[0].iType = TRUSS_NODE;
    nodes[1].dX = 3.0; nodes[1].dY = 4.0; nodes[1].iType = TRUSS_NODE;
    Element elems[1];
    elems[0].iaNode[0] = 0; elems[0].iaNode[1] = 1; elems[0].iType = TRUSS;
    LengthSinCosCalcu(1, elems, nodes);
    if (fabs(elems[0].dLength - 5.0) > 1e-12) {
        std::cerr << "LengthSinCosCalcu failed" << std::endl; return 1;
    }
    if (fabs(elems[0].dSin - 4.0/5.0) > 1e-12) { std::cerr << "Sin wrong" << std::endl; return 1; }

    // Test MatrixMultiply
    double *aRow0 = new double[2]{1.0,2.0};
    double *aRow1 = new double[2]{3.0,4.0};
    double** A = new double*[2]; A[0]=aRow0; A[1]=aRow1;
    double *bRow0 = new double[2]{2.0,0.0};
    double *bRow1 = new double[2]{1.0,2.0};
    double** B = new double*[2]; B[0]=bRow0; B[1]=bRow1;
    double** C = new double*[2]; C[0]=new double[2]; C[1]=new double[2];
    MatrixMultiply(2,2,2,A,B,C);
    if (fabs(C[0][0] - 4.0) > 1e-12) { std::cerr<<"MatMul failed"<<std::endl; return 1; }

    // Test LDLTSolve for trivial 1x1
    int nRow=1;
    int pDiag[1] = {0};
    double pGK[1] = {10.0};
    double pB[1] = {20.0};
    bool ok = LDLTSolve(1, pDiag, pGK, pB);
    if (!ok) { std::cerr<<"LDLTSolve failed"<<std::endl; return 1; }
    if (fabs(pB[0] - 2.0) > 1e-12) { std::cerr<<"LDLT result wrong"<<std::endl; return 1; }

    // Test LDLTSolve for 2x2 symmetric positive definite matrix
    int pDiag2[2] = {0,2};
    double pGK2[3] = {2.0, 1.0, 2.0};
    double pB2[2] = {3.0, 3.0};
    bool ok2 = LDLTSolve(2, pDiag2, pGK2, pB2);
    if (!ok2) { std::cerr<<"LDLTSolve 2x2 failed"<<std::endl; return 1; }
    if (fabs(pB2[0] - 1.0) > 1e-12 || fabs(pB2[1] - 1.0) > 1e-12) { std::cerr<<"LDLT 2x2 result wrong"<<std::endl; return 1; }

    // Test DOFIndexCalcu
    Node nodes2[2];
    nodes2[0].iType = TRUSS_NODE; nodes2[1].iType = TRUSS_NODE;
    for (int k=0;k<3;k++){ nodes2[0].iaDOFIndex[k]=0; nodes2[1].iaDOFIndex[k]=0; }
    ConstrainedNode cons[1];
    cons[0].iNode = 0; cons[0].iaConstrainedDOF[0] = -1; cons[0].iaConstrainedDOF[1] = 0; cons[0].iaConstrainedDOF[2] = 0;
    int nFreeDOF=0;
    int totalDOF = DOFIndexCalcu(nFreeDOF, 2, 1, cons, nodes2);
    if (nFreeDOF != 3) { std::cerr<<"DOF free count wrong"<<std::endl; return 1; }
    if (totalDOF != 4) { std::cerr<<"DOF total count wrong"<<std::endl; return 1; }
    if (nodes2[0].iaDOFIndex[0] != 3) { std::cerr<<"Constrained index wrong"<<std::endl; return 1; }

    // Test ElementDOFCalcu
    Element elems2[1]; elems2[0].iType = TRUSS; elems2[0].iaNode[0]=0; elems2[0].iaNode[1]=1;
    int** elemDOF = TwoArrayIntAlloc(1,6);
    ElementDOFCalcu(1, nodes2, elems2, elemDOF);
    if (elemDOF[0][0] != nodes2[0].iaDOFIndex[0]) { std::cerr<<"ElemDOF 0 mismatch"<<std::endl; return 1; }

    // Test BandAndDiagCalcu expected pDiag values [0,2,5,9]
    int* pDiagBand = new int[4];
    for (int i=0;i<4;i++) pDiagBand[i]=1;
    BandAndDiagCalcu(1, 4, elems2, elemDOF, pDiagBand);
    if (pDiagBand[0]!=0 || pDiagBand[1]!=2 || pDiagBand[2]!=5 || pDiagBand[3]!=9) { std::cerr<<"BandAndDiagCalcu wrong"<<std::endl; return 1; }

    // Test GKAssembly produces some nonzero entries
    int iBuf = pDiagBand[3] + 1;
    double* pGKBand = new double[iBuf];
    for (int i=0;i<iBuf;i++) pGKBand[i]=0.0;
    Material mates[1]; mates[0].dE = 210e9; mates[0].dMu=0.3; mates[0].dAlpha=0.0;
    Section sects[1]; sects[0].dA = 0.01; sects[0].dIz=0.0; sects[0].dH=0.0;
    GKAssembly(4, 1, elems2, nodes2, mates, sects, pDiagBand, pGKBand);
    bool anyNonZero=false;
    for (int i=0;i<iBuf;i++) if (fabs(pGKBand[i])>1e-12) anyNonZero=true;
    if (!anyNonZero) { std::cerr<<"GKAssembly produced zero matrix"<<std::endl; return 1; }

    // cleanup
    TwoArrayFree(1, elemDOF);
    delete[] pDiagBand; delete[] pGKBand;

    std::cout << "all tests passed" << std::endl;
    return 0;
}
