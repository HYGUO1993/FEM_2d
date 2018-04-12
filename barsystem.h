#pragma once
#define TRUSS
#define FRAME

#define TRUSS_NODE
#define FRAME_NODE

#define FORCE_ON_NODE
#define LATERAL_FORCE
#define LATERAL_UNIFORM_PRESSURE
#define MOMENT_ON_A_POINT
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
	double dE;					//����ģ��
	double dMu;					//���ɱ�
	double dAlpha;				//������ϵ��
};

struct Section {
	double dA;					//��������
	double dIz;					//���������Ծ�
	double dH;					//������
};

struct Node {
	int iType;					//�ڵ�����
	double dX, dY;				//�ڵ�����
	int iaDOFIndex[3];			//�ڵ����ɶȱ��
};

struct Element {
	int iType;					//��Ԫ���ͺ�
	int iaNode[2];				//��Ԫ���˽ڵ���
	int iSection;				//��Ԫ����������
	int iMaterial;				//��Ԫ����������
	double dLength;				//��Ԫ����
	double dSin, dCos;			//��Ԫ�ֲ�����x������������x��ļн�������

	double daEndInterForce[6];	//��Ԫ�˶�������

};

struct Load {
	int iType;					//�غ�����
	int iDirect;				//�������÷���0-X��1-Y��2-R��
	double dValue;				//�غ�ֵ
	int iLoadedElem;			//�غ����õĵ�Ԫ��
	int iLoadedNode;			//�غ����õĽڵ��
	double dPosition;			//�غ����õ�λ�û�ֲ�����
	double dT0, dT1;			//�����±����¶ȱ仯ֵ

};

struct ConstrainedNode {
	int iNode;					//��Լ���ڵ��
	int iaConstrainedDOF[3];	//10000+�ڵ�� N-��ڵ�N��ͬ�����ɶ����

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
void GKAssembly(int nTotalDOF, int nTotalElem, Element* pElem, Node* pNode, Material* pMate,
	Section* pSect, int* pDiag, double* pGK);
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
