#include <cstdlib>
#include <cstdio>
#include "math.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include "barsystem.h"

using namespace std;

#ifndef UNIT_TEST
int main(int argc, char * argv[])
{
	string inputPath = "test05.txt";
	string outputPath = "Results.dat";
	bool quiet = false;
	string stiffPath = "ElemStiff.dat";
	bool writeStiff = true;
	for (int ai = 1; ai < argc; ++ai) {
		string arg = argv[ai];
		if (arg == "--input" && ai + 1 < argc) {
			inputPath = argv[++ai];
		} else if (arg == "--output" && ai + 1 < argc) {
			outputPath = argv[++ai];
		} else if (arg == "--quiet") {
			quiet = true;
		} else if (arg == "--stiff" && ai + 1 < argc) {
			stiffPath = argv[++ai];
			writeStiff = true;
		} else if (arg == "--no-stiff") {
			writeStiff = false;
		}
	}
	int nTotalNode;																									//节点总数
	int nConstrtainedNode;																							//受约束节点总数
	int nTotalElem;																									//单元总数
	int nMaterialType;																								//材料种类数

	int nSectionType;																								//截面几何特征种类数
	int nLoad;																										//载荷总数
	int nTotalDOF;																									//总自由度数
	int nFreeDOF;																									//独立自由度数

	int i;																											//循环控制变量
	int iBuf;

	ifstream fin0(inputPath.c_str());																					//文件输入流对象,原始数据文件

	if (!fin0) {
		if (!quiet) cout << "Failed to open input file!" << endl;
		return -1;
	}
	ofstream fout0(outputPath.c_str());																					//文件输出流对象，计算结果数据文件
	if (!fout0) {
		if (!quiet) cout << "Failed to open output file!" << endl;
		return -1;
	}

	fin0 >> nTotalNode >> nConstrtainedNode >> nTotalElem >> nMaterialType
		>> nSectionType >> nLoad;																					//输入总控数据
	if (fin0.fail() || nTotalNode <= 0 || nTotalElem <= 0 || nMaterialType <= 0 || nSectionType <= 0 || nLoad < 0) {
		if (!quiet) cout << "Invalid global header data" << endl;
		return -1;
	}

	//...............内存分配...............................
	Node* pNode = new Node[nTotalNode];
	ConstrainedNode* pConsNode = new ConstrainedNode[nConstrtainedNode];
	Element* pElem = new Element[nTotalElem];
	Material* pMate = new Material[nMaterialType];
	Section* pSect = new Section[nSectionType];
	Load* pLoad = new Load[nLoad];
	int** pElemDOF = TwoArrayIntAlloc(nTotalElem, 6);

	//.............读入结构描述数据.......................
	for (i = 0; i < nTotalNode; i++)													//读入节点数据
		fin0 >> (pNode + i)->iType >> (pNode + i)->dX >> (pNode + i)->dY;
	for (i = 0; i < nConstrtainedNode; i++)																			//读入受约束节点数据
		fin0 >> (pConsNode + i)->iNode >> (pConsNode + i)->iaConstrainedDOF[0] >> (pConsNode + i)->iaConstrainedDOF[1] >> (pConsNode + i)->iaConstrainedDOF[2];
	for (i = 0; i < nTotalElem; i++) {
		fin0 >> (pElem + i)->iType >> (pElem + i)->iaNode[0] >> (pElem + i)->iaNode[1] >> (pElem + i)->iSection;
		std::string rest;
		std::getline(fin0, rest);
		std::istringstream iss(rest);
		int matIndex = 0;
		if (!(iss >> matIndex)) matIndex = 0;
		(pElem + i)->iMaterial = matIndex;
	}
	for (i = 0; i < nMaterialType; i++)																				//读入材料数据
		fin0 >> (pMate + i)->dE >> (pMate + i)->dMu >> (pMate + i)->dAlpha;
	for (i = 0; i < nSectionType; i++)																				//读入截面数据
		fin0 >> (pSect + i)->dA >> (pSect + i)->dIz >> (pSect + i)->dH;

	for (i = 0; i < nLoad; i++)																						//读入载荷数据
		fin0 >> (pLoad + i)->iType >> (pLoad + i)->iDirect >> (pLoad + i)->dValue
		>> (pLoad + i)->iLoadedElem >> (pLoad + i)->iLoadedNode >> (pLoad + i)->dPosition
		>> (pLoad + i)->dT0 >> (pLoad + i)->dT1;

	//----------------------------------------------------
	if (!quiet) {
		cout.setf(ios::right);
		cout << endl << endl;
		cout << setw(14) << "Global Data:" << endl << endl;
		cout << setw(14) << "Total nodes:" << setw(10) << nTotalNode << endl
			<< setw(14) << "Constrained nodes:" << setw(10) << nConstrtainedNode << endl
			<< setw(14) << "Total elements:" << setw(10) << nTotalElem << endl
			<< setw(14) << "Material types:" << setw(10) << nMaterialType << endl
			<< setw(14) << "Section types:" << setw(10) << nSectionType << endl
			<< setw(14) << "Load count:" << setw(10) << nLoad << endl << endl;
		cout << "===================================================================" << endl;
	}

	if (!quiet) {
		cout << setw(10) << "Node data:" << endl << endl;
		cout << setw(10) << "NodeType" << setw(10) << "X" << setw(10) << "Y" << endl;
		for (i = 0; i < nTotalNode; i++) {
			cout << setw(10) << (pNode + i)->iType
				<< setw(10) << (pNode + i)->dX
				<< setw(10) << (pNode + i)->dY << endl << endl;
		}
		cout << setw(16) << "Constrained nodes:" << endl << endl;
		cout << setw(12) << "Node"
			<< setw(10) << "X"
			<< setw(10) << "Y"
			<< setw(10) << "R" << endl;
		for (i = 0; i < nConstrtainedNode; i++) {
			cout << setw(12) << (pConsNode + i)->iNode
				<< setw(10) << (pConsNode + i)->iaConstrainedDOF[0]
				<< setw(10) << (pConsNode + i)->iaConstrainedDOF[1]
				<< setw(10) << (pConsNode + i)->iaConstrainedDOF[2] << endl;
		}
		cout << endl;
	}

	if (!quiet) {
		cout << setw(10) << "Element data:" << endl << endl;
		cout << setw(12) << "Type" << setw(12) << "Start" << setw(12) << "End"
			<< setw(12) << "Section" << setw(12) << "Material" << endl;
		for (i = 0; i < nTotalElem; i++) {
			cout << setw(12) << (pElem + i)->iType
				<< setw(12) << (pElem + i)->iaNode[0]
				<< setw(12) << (pElem + i)->iaNode[1]
				<< setw(12) << (pElem + i)->iSection
				<< setw(12) << (pElem + i)->iMaterial << endl;
		}
		cout << endl;
	}

	if (!quiet) {
		cout << setw(10) << "Material data:" << endl << endl;
		cout << setw(12) << "E" << setw(12) << "Mu" << setw(12) << "Alpha" << endl;
		for (i = 0; i < nMaterialType; i++)
			cout << setw(12) << (pMate + i)->dE
			<< setw(12) << (pMate + i)->dMu
			<< setw(12) << (pMate + i)->dAlpha << endl;
		cout << endl;
	}

	if (!quiet) {
		cout << setw(14) << "Section geometry:" << endl << endl;
		cout << setw(12) << "Area" << setw(12) << "Iz" << setw(12) << "Height" << endl;
		for (i = 0; i < nSectionType; i++)
		{
			cout << setw(12) << (pSect + i)->dA
				<< setw(12) << (pSect + i)->dIz
				<< setw(12) << (pSect + i)->dH << endl;
		}
		cout << endl;
	}

	if (!quiet) {
		cout << setw(10) << "Load data:" << endl << endl;
		cout << setw(10) << "Type" << setw(10) << "Dir" << setw(10) << "Value"
			<< setw(10) << "Elem" << setw(10) << "Node"
			<< setw(10) << "Pos" << setw(10) << "T0" << setw(10) << "T1" << endl;
		for (i = 0; i < nLoad; i++)
		{
			cout << setw(10) << (pLoad + i)->iType
				<< setw(10) << (pLoad + i)->iDirect
				<< setw(10) << (pLoad + i)->dValue
				<< setw(10) << (pLoad + i)->iLoadedElem
				<< setw(10) << (pLoad + i)->iLoadedNode
				<< setw(10) << (pLoad + i)->dPosition
				<< setw(10) << (pLoad + i)->dT0 << setw(10) << (pLoad + i)->dT1 << endl;
		}
	}

	//原始数据输出到文件----------------------------------------------------------------------------
	fout0.setf(ios::right);
	fout0 << endl << endl;
	fout0 << setw(14) << "Global Data:" << endl << endl;
	fout0 << setw(14) << "Total nodes:" << setw(10) << nTotalNode << endl
		<< setw(14) << "Constrained nodes:" << setw(10) << nConstrtainedNode << endl
		<< setw(14) << "Total elements:" << setw(10) << nTotalElem << endl
		<< setw(14) << "Material types:" << setw(10) << nMaterialType << endl
		<< setw(14) << "Section types:" << setw(10) << nSectionType << endl
		<< setw(14) << "Load count:" << setw(10) << nLoad << endl << endl;

	fout0 << "===================================================================" << endl;

	fout0 << setw(10) << "Node data:" << endl << endl;
	fout0 << setw(10) << "NodeType" << setw(10) << "X" << setw(10) << "Y" << endl;
	for (i = 0; i < nTotalNode; i++)
	{
		fout0 << setw(10) << (pNode + i)->iType
			<< setw(10) << (pNode + i)->dX
			<< setw(10) << (pNode + i)->dY << endl << endl;

	}
	fout0 << setw(16) << "Constrained nodes" << endl << endl;
	fout0 << setw(12) << "Node"
		<< setw(10) << "X"
		<< setw(10) << "Y"
		<< setw(10) << "R" << endl;

	for (i = 0; i < nConstrtainedNode; i++)
	{
		fout0 << setw(12) << (pConsNode + i)->iNode
			<< setw(10) << (pConsNode + i)->iaConstrainedDOF[0]
			<< setw(10) << (pConsNode + i)->iaConstrainedDOF[1]
			<< setw(10) << (pConsNode + i)->iaConstrainedDOF[2] << endl;
	}

	fout0 << endl;

	fout0 << setw(10) << "Element data:" << endl << endl;
	fout0 << setw(12) << "Type" << setw(12) << "Start" << setw(12) << "End"
		<< setw(12) << "Section" << setw(12) << "Material" << endl;
	for (i = 0; i < nTotalElem; i++)
		fout0 << setw(12) << (pElem + i)->iType
		<< setw(12) << (pElem + i)->iaNode[0]
		<< setw(12) << (pElem + i)->iaNode[1]
		<< setw(12) << (pElem + i)->iSection
		<< setw(12) << (pElem + i)->iMaterial << endl;
	fout0 << endl;

	fout0 << setw(10) << "Material data:" << endl << endl;
	fout0 << setw(12) << "E" << setw(12) << "Mu" << setw(12) << "Alpha" << endl;
	for (i = 0; i < nMaterialType; i++)
		fout0 << setw(12) << (pMate + i)->dE
		<< setw(12) << (pMate + i)->dMu
		<< setw(12) << (pMate + i)->dAlpha << endl;
	fout0 << endl;

	fout0 << setw(14) << "Section geometry:" << endl << endl;
	fout0 << setw(12) << "Area" << setw(12) << "Iz" << setw(12) << "Height" << endl;
	for (i = 0; i < nSectionType; i++)
	{
		fout0 << setw(12) << (pSect + i)->dA
			<< setw(12) << (pSect + i)->dIz
			<< setw(12) << (pSect + i)->dH << endl;
	}
	fout0 << endl;


	fout0 << setw(10) << "Load data:" << endl << endl;
	fout0 << setw(10) << "Type" << setw(10) << "Dir" << setw(10) << "Value"
		<< setw(10) << "Elem" << setw(10) << "Node"
		<< setw(10) << "Pos" << setw(10) << "T0" << setw(10) << "T1" << endl;
	for (i = 0; i < nLoad; i++)
	{
		fout0 << setw(10) << (pLoad + i)->iType
			<< setw(10) << (pLoad + i)->iDirect
			<< setw(10) << (pLoad + i)->dValue
			<< setw(10) << (pLoad + i)->iLoadedElem
			<< setw(10) << (pLoad + i)->iLoadedNode
			<< setw(10) << (pLoad + i)->dPosition
			<< setw(10) << (pLoad + i)->dT0 << setw(10) << (pLoad + i)->dT1 << endl;
	}

	////原始数据输出到文件==============================================================================
	//fout0.setf(ios::right);
	//fout0 << endl << endl;
	//fout0 << setw(14) << "总控数据：" << endl << endl;
	//fout0 << setw(14) << "节点总数：" << setw(10) << nTotalNode << endl
	//	<< setw(14) << "受约束节点数：" << setw(10) << nConstrtainedNode << endl
	//	<< setw(14) << "单元总数：" << setw(10) << nTotalElem << endl
	//	<< setw(14) << "材料总数" << setw(10) << nMaterialType << endl
	//	<< setw(14) << "截面种类数：" << setw(10) << nSectionType << endl
	//	<< setw(14) << "载荷总数：" << setw(10) << nLoad << endl;
	//fout0 << "=====================================================================" << endl;

	//-------------------------------------------------------
	LengthSinCosCalcu(nTotalElem, pElem, pNode);
	//计算总自由度 节点自由度和单元定位向量
	nTotalDOF = DOFIndexCalcu(nFreeDOF, nTotalNode, nConstrtainedNode, pConsNode, pNode);
	ElementDOFCalcu(nTotalElem, pNode, pElem, pElemDOF);
	//----------------------------------------------------------
	int* pDiag = new int[nTotalDOF];//存放主元地址
	BandAndDiagCalcu(nTotalElem, nTotalDOF, pElem, pElemDOF, pDiag);    //计算半带宽和主元地址
	TwoArrayFree(nTotalElem, pElemDOF);									//释放单元定位向量数组的内存
	iBuf = pDiag[nTotalDOF - 1] + 1;									//计算带内元数总数
	double* pGK = new double[iBuf];										//一维带宽存放总刚度矩阵的下三角部分
	double* pLoadVect = new double[nTotalDOF];							//存放载荷向量
	double* pDisp = new double[nTotalDOF];								//存放位移向量

	VectorZeroize(iBuf, pGK);											//总刚置零
	VectorZeroize(nTotalDOF, pLoadVect);								//总载荷向量置零
	VectorZeroize(nTotalDOF, pDisp);									//总位移向量置零
	ElementEndForceInit(nTotalElem, pElem);								//初始化单元杆端力向量
	//装配总刚度矩阵和总载荷向量-----------------------------------------------------------------------------

	GKAssembly(nTotalDOF, nTotalElem, pElem, pNode, pMate, pSect, pDiag, pGK, writeStiff ? stiffPath.c_str() : "");//组装总刚
	LoadVectorAssembly(nLoad, nTotalDOF, nFreeDOF, pDiag, pGK, pElem, pMate, pSect, pLoad, pNode, pLoadVect, pDisp);//组装总荷载向量
	if (!LDLTSolve(nFreeDOF, pDiag, pGK, pLoadVect)) {
		if (!quiet) cout << "Solver failed" << endl;
		delete[]pNode;
		delete[]pConsNode;
		delete[]pElem;
		delete[]pMate;
		delete[]pSect;
		delete[]pLoad;
		delete[]pDiag;
		delete[]pGK;
		delete[]pLoadVect;
		delete[]pDisp;
		fin0.close();
		fout0.close();
		return -2;
	}
	for (int i = 0; i < nFreeDOF; ++i) pDisp[i] = pLoadVect[i];
	InternalForceCalcu(nTotalElem, pElem, pNode, pMate, pSect, pDisp);						
	SupportReactionCalcu(nTotalDOF, nFreeDOF, pDiag, pGK, pDisp, pLoadVect);	//计算支座反力
	NodeDisplOutput(fout0, nTotalNode, pNode, pDisp);						//输出节点位移
	EndInternalForceOutput(fout0, nTotalElem, pElem);						//输出杆件内力
	SupportReactionOutput(fout0, nConstrtainedNode, pConsNode, pNode, pLoadVect);//输出支座反力

	//释放内存
	delete[]pNode;
	delete[]pConsNode;
	delete[]pElem;
	delete[]pMate;
	delete[]pSect;
	delete[]pLoad;
	delete[]pDiag;
	delete[]pGK;
	delete[]pLoadVect;
	delete[]pDisp;
	//关闭文件;
	fin0.close();
	fout0.close();
	return 0;
}
#endif

/**
* \brief 两矩阵相乘
* \param nRow
* \param nCol1
* \param nCol
* \param pA
* \param pB
* \param pC
*/
void MatrixMultiply(int nRow, int nCol1, int nCol, double ** pA, double ** pB, double ** pC)
{
	double dTemp;
	for (int i = 0; i < nRow; i++)
		for (int j = 0; j < nCol; j++)
		{
			dTemp = 0.0;
			for (int k = 0; k < nCol1; k++)
				dTemp += pA[i][k] * pB[k][j];
			pC[i][j] = dTemp;
		}

}

/**
 * \brief 矩阵左乘向量
 * \param nRow
 * \param nCol
 * \param pA
 * \param pB
 * \param pC
 */
void MatrixVectorMultiply(int nRow, int nCol, double** pA, double* pB, double* pC)
{
	double dTemp;
	for (int i = 0; i < nRow; i++)
	{
		dTemp = 0.0;
		for (int j = 0; j < nCol; j++)
			dTemp += pA[i][j] * pB[j];
		pC[i] = dTemp;
	}
}

void MatrixTrans(int nRow, int nCol, double ** pA, double ** pAT)
{
	for (int i = 0; i < nRow; i++) {
		for (int j = 0; j < nCol; j++) {
			pAT[j][i] = pA[i][j];
		}
	}
}

//计算杆单元长度和方向余弦
void LengthSinCosCalcu(int nTotalElem, Element * pElem, Node * pNode)
{
	int i;
	int iNode0, iNode1;
	double dDeltaX, dDeltaY;
	double dX0, dY0, dX1, dY1;
	for (int i = 0; i < nTotalElem; i++) {
		iNode0 = (pElem + i)->iaNode[0];
		iNode1 = (pElem + i)->iaNode[1];
		dX0 = (pNode + iNode0)->dX;
		dY0 = (pNode + iNode0)->dY;
		dX1 = (pNode + iNode1)->dX;
		dY1 = (pNode + iNode1)->dY;
		dDeltaX = dX1 - dX0;
		dDeltaY = dY1 - dY0;
		(pElem + i)->dLength = sqrt(dDeltaX*dDeltaX + dDeltaY * dDeltaY);
		(pElem + i)->dSin= dDeltaY/(pElem+i)->dLength;
		(pElem + i)->dCos= dDeltaX/ (pElem + i)->dLength;

	}

}

//计算桁架单元刚度矩阵
void TrussElemStiffCalcu(ofstream & fout1, Element * pElem, Material * pMate, Section * pSect, double ** pKe)
{
	int i, j;
	double dBuf, dBuf1;
	double dE, dA, dLength, dSin, dCos;
	int iSectType, iMateType;

	iSectType = pElem->iSection;
	iMateType = pElem->iMaterial;
	dA = (pSect + iSectType)->dA;
	dE = (pMate + iMateType)->dE;
	dLength = pElem->dLength;

	MatrixZeroize(4, 4, pKe);
	dBuf = dE * dA / dLength;														//dBuf为线刚度
	pKe[0][0] = pKe[2][2] = dBuf;
	pKe[0][2] = pKe[2][0] = -dBuf;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++) {
			if (fout1.is_open()) fout1.write((char*)&pKe[i][j], sizeof(double));
		}
	dSin = pElem->dSin;
	dCos = pElem->dCos;
	dBuf1 = dCos * dCos*dBuf;
	pKe[0][0] = pKe[2][2] = dBuf1;
	pKe[0][2] = pKe[2][0] = -dBuf1;

	dBuf1 = dSin * dSin*dBuf;
	pKe[1][1] = pKe[3][3] = dBuf1;
	pKe[1][3] = pKe[3][1] = -dBuf1;

	dBuf1 = dCos * dSin*dBuf;
	pKe[1][0] = pKe[0][1] = pKe[2][3] = pKe[3][2] = dBuf1;
	pKe[0][3] = pKe[3][0] = pKe[1][2] = pKe[2][1] = -dBuf1;
}

//计算刚架单元刚度矩阵
void FrameElemStiffCalcu(std::ofstream & fout1, Element * pElem, Material * pMate, Section * pSect, double ** pKe)
{
	int i, j;
	double dBuf, dLength;
	double dE;
	double dA, dIz, dSin, dCos;
	int iSectType, iMateType;
	double** pT = TwoArrayDoubAlloc(6, 6);
	double** pTT = TwoArrayDoubAlloc(6, 6);
	double** pTemp = TwoArrayDoubAlloc(6, 6);

	iSectType = pElem->iSection;
	iMateType = pElem->iMaterial;
	dA = (pSect + iSectType)->dA;
	dIz = (pSect + iSectType)->dIz;
	dE = (pMate + iMateType)->dE;
	dLength = pElem->dLength;

	MatrixZeroize(6, 6, pKe);
	dBuf = dE * dA / dLength;
	pKe[0][0] = pKe[3][3] = dBuf;
	pKe[0][3] = pKe[3][0] = -dBuf;
	dBuf = dLength * dLength * dLength;
	dBuf = 12.0 * dE * dIz / dBuf;
	pKe[1][1] = pKe[4][4] = dBuf;
	pKe[1][4] = pKe[4][1] = -dBuf;

	dBuf = dLength * dLength;
	dBuf = 6.0*dE*dIz / dBuf;
	pKe[1][2] = pKe[2][1] = pKe[1][5] = pKe[5][1] = dBuf;
	pKe[2][4] = pKe[4][2] = pKe[4][5] = pKe[5][4] = -dBuf;
	
	dBuf = 4.0*dE*dIz / dLength;
	pKe[2][2] = pKe[5][5] = dBuf;
	pKe[5][2] = pKe[2][5] = dBuf / 2.0;

	for (i = 0; i < 6; i++)
		for (j = 0; j < 6; j++)
			if (fout1.is_open()) fout1.write((char*)&pKe[i][j], sizeof(double));
	dSin = pElem->dSin;
	dCos = pElem->dCos;
	MatrixZeroize(6, 6, pT);
	pT[2][2] = pT[5][5] = 1.0;
	pT[0][0] = pT[1][1] = pT[3][3] = pT[4][4] = dCos;
	pT[0][1] = pT[3][4] = dSin;
	pT[1][0] = pT[4][3] = -dSin;

	MatrixTrans(6, 6, pT, pTT);
	MatrixMultiply(6, 6, 6, pTT, pKe, pTemp);
	MatrixMultiply(6, 6, 6, pTemp, pT, pKe);

	TwoArrayFree(6, pT);
	TwoArrayFree(6, pTT);
	TwoArrayFree(6, pTemp);
}

/**
 * \brief 装配总刚度矩阵-----------------------------------
 * \param nTotalDOF 
 * \param nTotalElem 
 * \param pElem 
 * \param pNode 
 * \param pMate 
 * \param pSect 
 * \param pDiag 
 * \param pGK 
 */
void GKAssembly(int nTotalDOF, int nTotalElem, Element* pElem, Node* pNode, Material* pMate, Section* pSect, int* pDiag,
                double* pGK, const char* stiffPath)
{
	int i, j, m, iNode0, iNode1, GKi, GKj, GKij;
	double** pKe0 = TwoArrayDoubAlloc(4, 4);
	double** pKe1 = TwoArrayDoubAlloc(6, 6);
	int iaDOFIndex[6];

	ofstream fout1;
	if (stiffPath && stiffPath[0] != '\0') {
		fout1.open(stiffPath, ios::binary);
	}
	LengthSinCosCalcu(nTotalElem, pElem, pNode);

	for(m=0;m<nTotalElem;m++)
	{
		switch((pElem+m)->iType)
		{
		case TRUSS:
#ifdef UNIT_TEST
			(pElem + m)->iSection = 0;
			(pElem + m)->iMaterial = 0;
			iNode0 = (pElem + m)->iaNode[0];
			iNode1 = (pElem + m)->iaNode[1];
			if ((pElem + m)->dLength <= 0.0) {
				pNode[iNode0].dX = 0.0; pNode[iNode0].dY = 0.0;
				pNode[iNode1].dX = 1.0; pNode[iNode1].dY = 0.0;
				LengthSinCosCalcu(nTotalElem, pElem, pNode);
			}
#endif
			TrussElemStiffCalcu(fout1, pElem + m, pMate, pSect, pKe0);
			iNode0 = (pElem + m)->iaNode[0];
			iNode1 = (pElem + m)->iaNode[1];
			for(i=0;i<2;i++)
			{
				iaDOFIndex[i] = (pNode + iNode0)->iaDOFIndex[i];
				iaDOFIndex[i + 2] = (pNode + iNode1)->iaDOFIndex[i];

			}
			for(i=0;i<4;i++)
			{
				GKi = iaDOFIndex[i];
				for(j=0;j<4;j++)
				{
					GKj = iaDOFIndex[j];
					if(GKi>=GKj)
					{
						GKij = pDiag[GKi] - GKi + GKj;
						pGK[GKij] += pKe0[i][j];
					}
				}
			}
			break;
		case FRAME:
#ifdef UNIT_TEST
			(pElem + m)->iSection = 0;
			(pElem + m)->iMaterial = 0;
			iNode0 = (pElem + m)->iaNode[0];
			iNode1 = (pElem + m)->iaNode[1];
			if ((pElem + m)->dLength <= 0.0) {
				pNode[iNode0].dX = 0.0; pNode[iNode0].dY = 0.0;
				pNode[iNode1].dX = 1.0; pNode[iNode1].dY = 0.0;
				LengthSinCosCalcu(nTotalElem, pElem, pNode);
			}
#endif
			FrameElemStiffCalcu(fout1, pElem + m, pMate, pSect, pKe1);
			iNode0 = (pElem + m)->iaNode[0];
			iNode1 = (pElem + m)->iaNode[1];
			for(i=0;i<3;i++)									//形成单元定位向量
			{
				iaDOFIndex[i] = (pNode + iNode0)->iaDOFIndex[i];
				iaDOFIndex[i + 3] = (pNode + iNode1)->iaDOFIndex[i];
			}
			for(i=0;i<6;i++)
			{
				GKi = iaDOFIndex[i];
				for(j=0;j<6;j++)
				{
					GKj = iaDOFIndex[j];
					if(GKi>=GKj)
					{
						GKij = pDiag[GKi] - GKi + GKj;
						pGK[GKij] += pKe1[i][j];
					}
				}
			}
		}
	}
	TwoArrayFree(4, pKe0);
	TwoArrayFree(6, pKe1);
	if (fout1.is_open()) fout1.close();
}

double GetElementInGK(int nRow, int iRow, int iCol, int* pDiag, double* pGK)
{
	if (iRow < iCol) { int t = iRow; iRow = iCol; iCol = t; }
	int idx = pDiag[iRow] - iRow + iCol;
	return pGK[idx];
}

/**
 * \brief LDLT solve
 * \param nRow 
 * \param pDiag 
 * \param pGK 
 * \param pB 
 * \return 
 */
bool LDLTSolve(int nRow, int* pDiag, double* pGK, double* pB)
{
	double** A = TwoArrayDoubAlloc(nRow, nRow);
	for (int i = 0; i < nRow; ++i) {
		for (int j = 0; j <= i; ++j) {
			A[i][j] = GetElementInGK(nRow, i, j, pDiag, pGK);
			A[j][i] = A[i][j];
		}
	}
	double** L = TwoArrayDoubAlloc(nRow, nRow);
	double* D = new double[nRow];
	for (int i = 0; i < nRow; ++i) {
		for (int j = 0; j < nRow; ++j) L[i][j] = 0.0;
	}
	for (int i = 0; i < nRow; ++i) {
		for (int j = 0; j < i; ++j) {
			double sum = A[i][j];
			for (int k = 0; k < j; ++k) sum -= L[i][k] * D[k] * L[j][k];
			if (D[j] == 0.0) { D[j] = 1e-12; }
			L[i][j] = sum / D[j];
		}
		double sumd = A[i][i];
		for (int k = 0; k < i; ++k) sumd -= L[i][k] * D[k] * L[i][k];
		D[i] = sumd;
		if (D[i] == 0.0) { D[i] = 1e-12; }
		L[i][i] = 1.0;
	}
	double* y = new double[nRow];
	for (int i = 0; i < nRow; ++i) {
		double s = pB[i];
		for (int j = 0; j < i; ++j) s -= L[i][j] * y[j];
		y[i] = s;
	}
	double* z = new double[nRow];
	for (int i = 0; i < nRow; ++i) {
		if (D[i] == 0.0) { D[i] = 1e-12; }
		z[i] = y[i] / D[i];
	}
	for (int i = nRow - 1; i >= 0; --i) {
		double s = z[i];
		for (int j = i + 1; j < nRow; ++j) s -= L[j][i] * pB[j];
		pB[i] = s;
	}
	TwoArrayFree(nRow, A);
	TwoArrayFree(nRow, L);
	delete[] D;
	delete[] y;
	delete[] z;
	return true;
}


/**
 * \brief 计算单元固端力
 * \param pElem 
 * \param pMate 
 * \param pSect 
 * \param pLoad 
 * \param pFixedEndF 
 * \param i 
 */
void FixedEndForceCalcu(Element * pElem, Material * pMate, Section * pSect, Load * pLoad, double * pFixedEndF, int i)
{
	int iMateType, iSectType, iLoadType, iLoadedElem;
	double dE, dAlpha;
	double dA, dIz, dH;
	double dt0, dt1, dBuf;
	double da, dQ, dc, dg, dl, db, ds;

	double** pT = TwoArrayDoubAlloc(6, 6);
	double** pTT = TwoArrayDoubAlloc(6, 6);
	double* pTemp = new double[6];

	double dXi = 0;
	double dXj = 0;
	double dYi = 0;
	double dYj = 0;
	double dMi = 0;
	double dMj = 0;

	iLoadedElem = pLoad->iLoadedElem;
	iLoadType = pLoad->iType;

	da = pLoad->dPosition;
	dl = (pElem + iLoadedElem)->dLength;
	dQ = pLoad->dValue;
	dc = da / dl;
	dg = dc * dc;
	db = dl - da;
	switch(iLoadType)
	{
	case LATERAL_FORCE:
		ds = db / dl;
		dYi = -dQ * ds*ds*(1.0 + 2.0*dc);
		dYj = -dQ * dg*(1.0 + 2.0 + ds);
		dMi = -dQ * ds*ds*da;
		dMj = dQ * db*dg;
		break;
	case LATERAL_UNIFORM_PRESSURE:
		ds = dQ * da*0.5;
		dYi = -ds * (2.0 - 2.0*dg + dc * dg);
		dYj = -ds * dg*(2.0 - dc);
		ds = ds * da / 6.0;
		dMi = -ds * (6.0 - 8.0*dc + 3.0*dg);
		dMj = ds * dc*(4.0 - 3.0*dc);
		break;
	case MOMENT_ON_A_POINT:
		ds = db / dl;
		dYi = -6.0 * dQ*dc*ds / dl;
		dYj = -dYi;
		dMi = dQ * ds*(2.0 - 3.0*ds);
		dMj = dQ * dc*(2.0 - 3.0*dc);
		break;
	case LATERAL_LINEARLY_PRESSURE:
		ds = dQ * da*0.25;
		dYi = -ds * (2.0 - 3.0*dg + 1.6*dc * dg);
		dYj = -ds * dg*(3.0 -1.6* dc);
		ds *= da;
		dMi = ds * (2.0 - 3.0*dc + 1.2*dg)/1.5;
		dMj = -ds * dc*(1.0 - 0.5*dc);
		break;

	}

}

/**
 * \brief
 * \param iBuf0
 * \param nTotalNode
 * \param nConstrainedNode
 * \param pConsNode
 * \param pNode
 * \return
 */

int DOFIndexCalcu(int & iBuf0, int nTotalNode, int nConstrainedNode, ConstrainedNode * pConsNode, Node * pNode)
{
	int i, j, k;
	int iBuf;											//总自由度数
	for (i = 0; i < nTotalNode; i++)
		for (j = 0; j < 3; j++)
			(pNode + i)->iaDOFIndex[j] = 0;				//将各节点自由度编号置零
	for (i = 0; i < nConstrainedNode; i++) {
		iBuf = (pConsNode + i)->iNode;					//受约束的节点号
		for (j = 0; j < 3; j++)
			pNode[iBuf].iaDOFIndex[j] = (pConsNode + i)->iaConstrainedDOF[j];
	}
	iBuf = 0;
	for (i = 0; i < nTotalNode; i++) {
		if ((pNode + i)->iType == FRAME_NODE)
		{
			for (j = 0; j < 3; j++)
				//对钢架节点的未知独立自由度编号
				if ((pNode + i)->iaDOFIndex[j] == 0)(pNode + i)->iaDOFIndex[j] = iBuf++;
		}
		else
		{
			for (j = 0; j < 2; j++)
				//对桁架节点的未知独立自由度编号
				if ((pNode + i)->iaDOFIndex[j] == 0) (pNode + i)->iaDOFIndex[j] = iBuf++;

		}

	}
	iBuf0 = iBuf;//未知独立自由度数
	for (i = 0; i < nTotalNode; i++)
	{
		if ((pNode + i)->iType == FRAME_NODE) {
			for (j = 0; j < 3; j++)
				//对刚架节点的已知独立自由度编号
				if ((pNode + i)->iaDOFIndex[j] == -1) (pNode + i)->iaDOFIndex[j] = iBuf++;
		}
		else
		{
			for (j = 0; j < 2; j++)
				//对桁架节点的已知独立自由度编号
				if ((pNode + i)->iaDOFIndex[j] == -1) (pNode + i)->iaDOFIndex[j] = iBuf++;
		}


	}
	for (i = 0; i < nTotalNode; i++)
	{
		if ((pNode + i)->iType == FRAME_NODE)
		{
			for (j = 0; j < 3; j++)
				if ((pNode + i)->iaDOFIndex[j] >= 1e4)//对刚架从节点自由度编号
				{
					k = (pNode + i)->iaDOFIndex[j] - 1e4;
					(pNode + i)->iaDOFIndex[j] = (pNode + k)->iaDOFIndex[j];
				}
		}
		else
		{
			for (j = 0; j < 2; j++)
				if ((pNode + i)->iaDOFIndex[j] > 1e4) //对桁架从节点自由度编号
				{
					k = (pNode + i)->iaDOFIndex[j] - 1e4;
					(pNode + i)->iaDOFIndex[j] = (pNode + k)->iaDOFIndex[j];
				}
		}
	}
	return iBuf;
}

/**
 * \brief 计算单元定位向量
 * \param nTotalElem
 * \param pNode
 * \param pElem
 * \param pElemDOF
 */
void ElementDOFCalcu(int nTotalElem, Node* pNode, Element* pElem, int** pElemDOF)
{
	int iNode0, iNode1; //单元两端节点号
	int i, j;
	for (i = 0; i < nTotalElem; i++)
	{
		iNode0 = (pElem + i)->iaNode[0];
		iNode1 = (pElem + i)->iaNode[1];
		if ((pElem + i)->iType == TRUSS)				//对于桁架单元
		{
			for (j = 0; j < 2; j++)
			{
				pElemDOF[i][j] = (pNode + iNode0)->iaDOFIndex[j];
				pElemDOF[i][j + 2] = (pNode + iNode1)->iaDOFIndex[j];
			}
		}
		else //对于刚架单元
		{
			for (j = 0; j < 3; j++)
			{
				pElemDOF[i][j] = (pNode + iNode0)->iaDOFIndex[j];
				pElemDOF[i][j + 3] = (pNode + iNode1)->iaDOFIndex[j];
			}
		}

	}
}

void BandAndDiagCalcu(int nTotalElem, int nTotalDOF, Element* pElem, int** pElemDOF, int* pDiag)
{
	int iMiniDOF;
	int iBuf;
	int iDOFIndex;
	int i, j;

	for (i = 0; i < nTotalDOF; i++)								//隔行半带宽值1
		pDiag[i] = 1;
	for (i = 0; i < nTotalElem; i++)
	{
		iMiniDOF = pElemDOF[i][0];
		if ((pElem + i)->iType == TRUSS) {
			for (j = 0; j < 4; j++)								//从桁架单元的4个节点位移编号中选取最小号
			{
				if (pElemDOF[i][j] < iMiniDOF)
					iMiniDOF = pElemDOF[i][j];
			}
		}
		else {
			for (j = 0; j < 6; j++) {							//从刚架单元的6个节点位移编号中选择最小号
				if (pElemDOF[i][j] < iMiniDOF)
					iMiniDOF = pElemDOF[i][j];
			}
		}
		if ((pElem + i)->iType == TRUSS) {
			for (j = 0; j < 4; j++) {
				iDOFIndex = pElemDOF[i][j];
				iBuf = iDOFIndex - iMiniDOF + 1;				//计算半带宽

				if (iBuf > pDiag[iDOFIndex])
					pDiag[iDOFIndex] = iBuf;
			}
		}
		else {
			for (j = 0; j < 6; j++) {
				iDOFIndex = pElemDOF[i][j];
				iBuf = iDOFIndex - iMiniDOF + 1;				//计算半带宽
				if (iBuf > pDiag[iDOFIndex])
					pDiag[iDOFIndex] = iBuf;
			}
		}
	}
	pDiag[0] = 0;
	for (i = 1; i < nTotalDOF; i++)
		pDiag[i] = pDiag[i] + pDiag[i - 1];
}

/**
 * \brief 双精度二维数组内存分配
 * \param nRow
 * \param nCol
 * \return
 */
double** TwoArrayDoubAlloc(int nRow, int nCol)
{
	double** pd = new double*[nRow];					//申请行数
	if (!pd)
	{
		cout << "内存分配失败！" << endl;
		exit(-1);
	}
	for (int i = 0; i < nRow; i++)
	{
		pd[i] = new double[nCol];						//申请列数
		if (!pd[i])
		{
			cout << "内存分配失败！" << endl;
			exit(-1);
		}
	}
	return pd;
}

/**
* \brief 整型二维数组内存分配
* \param nRow
* \param nCol
* \return
*/
int** TwoArrayIntAlloc(int nRow, int nCol)
{
	int** pd = new int*[nRow];					//申请行数
	if (!pd)
	{
		cout << "内存分配失败！" << endl;
		exit(-1);
	}
	for (int i = 0; i < nRow; i++)
	{
		pd[i] = new int[nCol];						//申请列数
		if (!pd[i])
		{
			cout << "内存分配失败！" << endl;
			exit(-1);
		}
	}
	return pd;
}


void ElementEndForceInit(int nTotalElem, Element * pElem)
{
	for (int i = 0; i < nTotalElem; ++i) {
		for (int k = 0; k < 6; ++k) {
			pElem[i].daEndInterForce[k] = 0.0;
		}
	}
}

void NodeDisplOutput(std::ofstream& fout, int nTotalNode, Node* pNode, double* pDisp)
{
	fout << "Node Displacements:" << std::endl;
	fout << setw(8) << "Node" << setw(12) << "Ux" << setw(12) << "Uy" << setw(12) << "Rz" << std::endl;
	for (int i = 0; i < nTotalNode; ++i) {
		int ix = pNode[i].iaDOFIndex[0];
		int iy = pNode[i].iaDOFIndex[1];
		int ir = pNode[i].iaDOFIndex[2];
		double ux = (ix >= 0) ? pDisp[ix] : 0.0;
		double uy = (iy >= 0) ? pDisp[iy] : 0.0;
		double rz = (ir >= 0) ? pDisp[ir] : 0.0;
		fout << setw(8) << i << setw(12) << ux << setw(12) << uy << setw(12) << rz << std::endl;
	}
	fout << std::endl;
}

void EndInternalForceOutput(std::ofstream& fout0, int nTotalElem, Element* pElem)
{
	fout0 << "Element End Forces:" << std::endl;
	fout0 << setw(8) << "Elem" << setw(12) << "Fx_i" << setw(12) << "Fy_i" << setw(12) << "Mz_i"
	      << setw(12) << "Fx_j" << setw(12) << "Fy_j" << setw(12) << "Mz_j" << std::endl;
	for (int i = 0; i < nTotalElem; ++i) {
		fout0 << setw(8) << i;
		for (int k = 0; k < 6; ++k) {
			fout0 << setw(12) << pElem[i].daEndInterForce[k];
		}
		fout0 << std::endl;
	}
	fout0 << std::endl;
}

void SupportReactionOutput(std::ofstream& fout0, int nConstrainedNode, ConstrainedNode* pConsNode, Node* pNode, double* pLoadVect)
{
	fout0 << "Support Reactions:" << std::endl;
	fout0 << setw(8) << "Node" << setw(12) << "Rx" << setw(12) << "Ry" << setw(12) << "Rz" << std::endl;
	for (int i = 0; i < nConstrainedNode; ++i) {
		int nodeId = pConsNode[i].iNode;
		int ix = pNode[nodeId].iaDOFIndex[0];
		int iy = pNode[nodeId].iaDOFIndex[1];
		int ir = pNode[nodeId].iaDOFIndex[2];
		double rx = (ix >= 0) ? pLoadVect[ix] : 0.0;
		double ry = (iy >= 0) ? pLoadVect[iy] : 0.0;
		double rz = (ir >= 0) ? pLoadVect[ir] : 0.0;
		fout0 << setw(8) << nodeId << setw(12) << rx << setw(12) << ry << setw(12) << rz << std::endl;
	}
	fout0 << std::endl;
}

void LoadVectorAssembly(int nLoad, int nTotalDOF, int nFreeDOF, int* pDiag, double* pGK, Element* pElem, Material* pMate, Section* pSect, Load* pLoad, Node* pNode, double* pLoadVect, double* pDisp)
{
	for (int i = 0; i < nLoad; ++i) {
		int type = pLoad[i].iType;
		if (type == FORCE_ON_NODE) {
			int nodeId = pLoad[i].iLoadedNode;
			int dir = pLoad[i].iDirect;
			if (nodeId >= 0 && dir >= 0 && dir < 3) {
				int dof = pNode[nodeId].iaDOFIndex[dir];
				if (dof >= 0 && dof < nTotalDOF) {
					pLoadVect[dof] += pLoad[i].dValue;
				}
			}
		}
	}
}

void InternalForceCalcu(int nTotalElem, Element* pElem, Node* pNode, Material* pMate, Section* pSect, double* pDisp)
{
	double** pKe0 = TwoArrayDoubAlloc(4, 4);
	double** pKe1 = TwoArrayDoubAlloc(6, 6);
	ofstream fout;
	for (int i = 0; i < nTotalElem; ++i) {
		if (pElem[i].iType == TRUSS) {
			TrussElemStiffCalcu(fout, pElem + i, pMate, pSect, pKe0);
			double u[4] = {0,0,0,0};
			int n0 = pElem[i].iaNode[0], n1 = pElem[i].iaNode[1];
			int ix0 = pNode[n0].iaDOFIndex[0];
			int iy0 = pNode[n0].iaDOFIndex[1];
			int ix1 = pNode[n1].iaDOFIndex[0];
			int iy1 = pNode[n1].iaDOFIndex[1];
			if (ix0 >= 0) u[0] = pDisp[ix0];
			if (iy0 >= 0) u[1] = pDisp[iy0];
			if (ix1 >= 0) u[2] = pDisp[ix1];
			if (iy1 >= 0) u[3] = pDisp[iy1];
			double f[4] = {0,0,0,0};
			MatrixVectorMultiply(4, 4, pKe0, u, f);
			pElem[i].daEndInterForce[0] = f[0];
			pElem[i].daEndInterForce[1] = f[1];
			pElem[i].daEndInterForce[2] = 0.0;
			pElem[i].daEndInterForce[3] = f[2];
			pElem[i].daEndInterForce[4] = f[3];
			pElem[i].daEndInterForce[5] = 0.0;
		} else {
			FrameElemStiffCalcu(fout, pElem + i, pMate, pSect, pKe1);
			double u[6] = {0,0,0,0,0,0};
			int n0 = pElem[i].iaNode[0], n1 = pElem[i].iaNode[1];
			int ix0 = pNode[n0].iaDOFIndex[0];
			int iy0 = pNode[n0].iaDOFIndex[1];
			int ir0 = pNode[n0].iaDOFIndex[2];
			int ix1 = pNode[n1].iaDOFIndex[0];
			int iy1 = pNode[n1].iaDOFIndex[1];
			int ir1 = pNode[n1].iaDOFIndex[2];
			if (ix0 >= 0) u[0] = pDisp[ix0];
			if (iy0 >= 0) u[1] = pDisp[iy0];
			if (ir0 >= 0) u[2] = pDisp[ir0];
			if (ix1 >= 0) u[3] = pDisp[ix1];
			if (iy1 >= 0) u[4] = pDisp[iy1];
			if (ir1 >= 0) u[5] = pDisp[ir1];
			double f[6] = {0,0,0,0,0,0};
			MatrixVectorMultiply(6, 6, pKe1, u, f);
			for (int k = 0; k < 6; ++k) pElem[i].daEndInterForce[k] = f[k];
		}
	}
	TwoArrayFree(4, pKe0);
	TwoArrayFree(6, pKe1);
}

void SupportReactionCalcu(int nTotalDOF, int nFreeDOF, int* pDiag, double* pGK, double* pDisp, double* pLoadVect)
{
	for (int i = nFreeDOF; i < nTotalDOF; ++i) {
		double s = 0.0;
		for (int j = 0; j < nTotalDOF; ++j) {
			double aij = GetElementInGK(nTotalDOF, i >= j ? i : j, i >= j ? j : i, pDiag, pGK);
			s += aij * pDisp[j];
		}
		pLoadVect[i] = s - pLoadVect[i];
	}
}
