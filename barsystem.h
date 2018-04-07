#pragma once
#define TRUSS
#define FRAME

#define TRUSS_NODE
#define FRAME_NODE

#define FORCE_ON_NODE
#define LATERAL_FORCE
#define LATERAL_UNIFORM_PRESSURE
#define MOMEN_ON_A_POINT
#define LATERAL_LINEARLY_PRESSURE
#define AXIAL_PRESSURE
#define AXIAL_FORCE
#define MOMENT_ON_BEAM
#define TEMPERATURE
#define SUPPORT_MOVE

#define DIRECT_X
#define DIRECT_Y
#define DIRECT_R

struct Material {
	double dE;					//弹性模量
	double dMu;					//泊松比
	double dAlpha;				//线膨胀系数
};

struct Section {
	double dA;					//横截面面积
	double dIz;					//横截面积惯性矩
	double dH;					//横截面高
};

struct Node {
	int iType;					//节点类型
	double dX, dY;				//节点坐标
	int iaDOFIndex[3];			//节点自由度编号
};

struct Element {
	int iType;					//单元类型号
	int iaNode[2];				//单元两端节点编号
	int iSection;				//单元截面索引号
	int iMaterial;				//单元材料索引号
	double dLength;				//单元长度
	double dSin, dCos;			//单元局部坐标x轴与整体坐标x轴的夹角正余弦

	double daEndInterForce[6];	//单元杆端力向量

};

struct Load {
	int iType;					//载荷类型
	int iDirect;				//荷载作用方向：0-X向，1-Y向，2-R向
	double dValue;				//载荷值
	int iLoadedElem;			//载荷作用的单元号
	int iLoadedNode;			//载荷作用的节点号
	double dPosition;			//载荷作用的位置或分布长度
	double dT0, dT1;			//杆上下表面温度变化值

};

struct ConstrainedNode {
	int iNode;					//受约束节点号
	int iaConstrainedDOF[3];	//10000+节点号 N-与节点N的同类自由度耦合

};

int DOFIndexCalcu(int& iBuf0, int nTotalNode, int nConstrainedNode, ConstrainedNode* pConsNode, Node* pNode);
void ElementDOFCalcu(int nTotalElem, Node* pNode, Element *pElem, int** pElemDOF);
void BandAndDiagCalcu(int nTotalElem, int nTotalDOF, Element *pElem, int** pElemDOF, int* pDiag);
double ** TwoArrayDoubAlloc(int nRow, int nCol);
int** TwoArrayIntAlloc(int nRow, int nCol);
template <class T>
void TwoArrayFree(int nRow, T** pdi);
template <class T>
void MatrixZeroize(int nRow, int nCol, T** pT);
template <class T>
void VectorZeroize(int n, T* pT);
void MatrixMultiply(int nRow, int nCOl1, int nCol, double** pA, double** pB, double ** pC);
void MatrixVectorMultiply(int nRow, int nCol, double** pA, double* pB, double* pC);
void MatrixTrans(int nRow, int nCol, double** pA, double** pAT);
void LengthSinCosCalcu(int nTotalElem, Element* pElem, Node* pNode);
void TrussElemStiffCalcu(std::ofstream& fout1, Element* pElem, Material* pMate, Section* pSect, double** pKe);
void FrameElemStiffCalcu(std::ofstream& fout1, Element* pElem, Material* pMate, Section* pSect, double** pKe);
void GKAssembly(int nTotalDOF, int nTotalElem, Element* pElem, Node* pNode, Material* pMate, Section* pSect, int* pDiag, double* pGK);
bool LDLTSolve(int nRow, int* pDiag, double* pGK, double* pB);
void FixedEndForceCalcu(Element* pElem, Material* pMate, Section* pSect, Load* pLoad, double* pFixedEndF, int i);
void LoadVectorAssembly(int nLoad, int nTotalDOF, int nFreeDOF, int* pDiag, double* pGK, Element* pElem, Material* pMate, Section* pSect, Load* pLoad, Node* pNode, double* pLoadVect, double* pDisp);
void ElementEndForceInit(int nTotalElem, Element* pElem);
double GetElementInGK(int nRow, int iRow, int iCol, int* pDiag, double* pGK);
void LoadVectorModify(int nTotalDOF, int nFreeDOF, int* pDiag, double* pGK, double*pDisp, double* pLoadVect);
void InternalForceCalcu(int nTotalElem, Element* pElem, Node* pNode, double* pDisp);
void TrussInternalForceCalcu(Element* pElem, Node* pNode, double* pDisp, double** pKe);
void FrameInternalForceCalcu(Element* pElem, Node* pNode, double* pDisp, double** pKe);
void SupportReactionCalcu(int nTotalDOF, int nFreeDOF, int* pDiag, double* pGK, double* pDisp, double* pLoadVect);
void NodeDisplOutput(std::ofstream& fout, int nTotalNode, Node* pNode, double* pDisp);
void EndInternalForceOutput(std::ofstream& fout0, int nTotalElem, Element* pElem);
void SupportReactionOutput(std::ofstream& fout0, int nConstrainedNode, ConstrainedNode* pConsNode, Node* pNode, double* pLoadVect);
