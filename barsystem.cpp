#include "stdafx.h"
#include <cstdlib>
#include <cstdio>
#include "math.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "barsystem.h"

using namespace std;

int main(int argc, char * argv[])
{
	int nTotalNode;																									//�ڵ�����
	int nConstrtainedNode;																							//��Լ���ڵ�����
	int nTotalElem;																									//��Ԫ����
	int nMaterialType;																								//����������

	int nSectionType;																								//���漸������������
	int nLoad;																										//�غ�����
	int nTotalDOF;																									//�����ɶ���
	int nFreeDOF;																									//�������ɶ���

	int i;																											//ѭ�����Ʊ���
	int iBuf;

	ifstream fin0("test05.txt");																					//�ļ�����������,ԭʼ�����ļ�

	if (!fin0) {
		cout << "ԭʼ�����ļ���ʧ�ܣ�" << endl;
		exit(-1);
	}
	ofstream fout0("Results.dat");																					//�ļ���������󣬼����������ļ�
	if (!fout0) {
		cout << "����������ļ���ʧ��" << endl;
		exit(-1);
	}

	fin0 >> nTotalNode >> nConstrtainedNode >> nTotalElem >> nMaterialType
		>> nSectionType >> nLoad;																					//�����ܿ�����

	//...............�ڴ����...............................
	Node* pNode = new Node[nTotalNode];
	ConstrainedNode* pConsNode = new ConstrainedNode[nConstrtainedNode];
	Element* pElem = new Element[nTotalElem];
	Material* pMate = new Material[nMaterialType];
	Section* pSect = new Section[nSectionType];
	Load* pLoad = new Load[nLoad];
	int** pElemDOF = TwoArrayIntAlloc(nTotalElem, 6);

	//.............����ṹ��������.......................
	for (i = 0; i < nTotalElem; i++)																				//����ڵ�����
		fin0 >> (pNode + i)->iType >> (pNode + i)->dX >> (pNode + i)->dY;
	for (i = 0; i < nConstrtainedNode; i++)																			//������Լ���ڵ�����
		fin0 >> (pConsNode + i)->iNode >> (pConsNode + i)->iaConstrainedDOF[0] >> (pConsNode + i)->iaConstrainedDOF[1] >> (pConsNode + i)->iaConstrainedDOF[2];
	for (i = 0; i < nTotalElem; i++)																				//���뵥Ԫ����
		fin0 >> (pElem + i)->iType >> (pElem + i)->iaNode[0] >> (pElem + i)->iaNode[1] >> (pElem + i)->iSection >> (pElem + i)->iMaterial;
	for (i = 0; i < nMaterialType; i++)																				//�����������
		fin0 >> (pMate + i)->dE >> (pMate + i)->dMu >> (pMate + i)->dAlpha;
	for (i = 0; i < nSectionType; i++)																				//�����������
		fin0 >> (pSect + i)->dA >> (pSect + i)->dIz >> (pSect + i)->dH;

	for (i = 0; i < nLoad; i++)																						//�����غ�����
		fin0 >> (pLoad + i)->iType >> (pLoad + i)->iDirect >> (pLoad + i)->dValue
		>> (pLoad + i)->iLoadedElem >> (pLoad + i)->iLoadedNode >> (pLoad + i)->dPosition
		>> (pLoad + i)->dT0 >> (pLoad + i)->dT1;

	//----------------------------------------------------
	cout.setf(ios::right);
	cout << endl << endl;
	cout << setw(14) << "�ܿ����ݣ�" << endl << endl;
	cout << setw(14) << "�ڵ�������" << setw(10) << nTotalNode << endl
		<< setw(14) << "��Լ���Ľڵ�����" << setw(10) << nConstrtainedNode << endl
		<< setw(14) << "��Ԫ������" << setw(10) << nTotalElem << endl
		<< setw(14) << "����������" << setw(10) << nMaterialType << endl
		<< setw(14) << "������������" << setw(10) << nSectionType << endl
		<< setw(14) << "������������" << setw(10) << nLoad << endl << endl;
	cout << "===================================================================" << endl;

	cout << setw(10) << "�ڵ����ݣ�" << endl << endl;
	cout << setw(10) << "�ڵ����ͣ�" << setw(10) << "X" << setw(10) << "Y" << endl;
	for (i = 0; i < nTotalNode; i++) {
		cout << setw(10) << (pNode + i)->iType
			<< setw(10) << (pNode + i)->dX
			<< setw(10) << (pNode + i)->dY << endl << endl;

	}
	cout << setw(16) << "��Լ���ڵ����ݣ�" << endl << endl;
	cout << setw(12) << "��Լ���ڵ��"
		<< setw(10) << "X��������"
		<< setw(10) << "Y��������"
		<< setw(10) << "R��������" << endl;

	for (i = 0; i < nConstrtainedNode; i++) {
		cout << setw(12) << (pConsNode + i)->iNode
			<< setw(10) << (pConsNode + i)->iaConstrainedDOF[0]
			<< setw(10) << (pConsNode + i)->iaConstrainedDOF[1]
			<< setw(10) << (pConsNode + i)->iaConstrainedDOF[2] << endl;
	}
	cout << endl;

	cout << setw(10) << "��Ԫ���ݣ�" << endl << endl;
	cout << setw(12) << "��Ԫ����" << setw(12) << "ʼ�˽ڵ��" << setw(12) << "�ն˽ڵ��"
		<< setw(12) << "�����ͺ�" << setw(12) << "����������" << endl;

	for (i = 0; i < nTotalElem; i++) {
		cout << setw(12) << (pElem + i)->iType
			<< setw(12) << (pElem + i)->iaNode[0]
			<< setw(12) << (pElem + i)->iaNode[1]
			<< setw(12) << (pElem + i)->iSection
			<< setw(12) << (pElem + i)->iMaterial << endl;

	}
	cout << endl;

	cout << setw(10) << "�������ݣ�" << endl << endl;
	cout << setw(12) << "����ģ��" << setw(12) << "���ɱ�" << setw(12) << "������ϵ��" << endl;
	for (i = 0; i < nMaterialType; i++)
		cout << setw(12) << (pMate + i)->dE
		<< setw(12) << (pMate + i)->dMu
		<< setw(12) << (pMate + i)->dAlpha << endl;
	cout << endl;

	cout << setw(14) << "���漸��������" << endl << endl;
	cout << setw(12) << "���" << setw(12) << "���Ծ�" << setw(12) << "�����" << endl;
	for (i = 0; i < nSectionType; i++)
	{
		cout << setw(12) << (pSect + i)->dA
			<< setw(12) << (pSect + i)->dIz
			<< setw(12) << (pSect + i)->dH << endl;
	}
	cout << endl;

	cout << setw(10) << "�������ݣ�" << endl << endl;
	cout << setw(10) << "���ͺ�" << setw(10) << "����" << setw(10) << "�غ�ֵ"
		<< setw(10) << "��Ԫ��" << setw(10) << "�ڵ��"
		<< setw(10) << "λ��" << setw(10) << "�¶�T0" << setw(10) << "�¶�T1" << endl;
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

	//ԭʼ����������ļ�----------------------------------------------------------------------------
	fout0.setf(ios::right);
	fout0 << endl << endl;
	fout0 << setw(14) << "�ܿ����ݣ�" << endl << endl;
	fout0 << setw(14) << "�ڵ�������" << setw(10) << nTotalNode << endl
		<< setw(14) << "��Լ���ڵ�����" << setw(10) << nConstrtainedNode << endl
		<< setw(14) << "��Ԫ������" << setw(10) << nTotalElem << endl
		<< setw(14) << "����������" << setw(10) << nMaterialType << endl
		<< setw(14) << "������������" << setw(10) << nSectionType << endl
		<< setw(14) << "�غ���������" << setw(10) << nLoad << endl << endl;

	fout0 << "===================================================================" << endl;

	fout0 << setw(10) << "�ڵ����ݣ�" << endl << endl;
	fout0 << setw(10) << "�ڵ�����" << setw(10) << "X" << setw(10) << "Y" << endl;
	for (i = 0; i < nTotalNode; i++)
	{
		fout0 << setw(10) << (pNode + i)->iType
			<< setw(10) << (pNode + i)->dX
			<< setw(10) << (pNode + i)->dY << endl << endl;

	}
	fout0 << setw(16) << "��Լ���ڵ�����" << endl << endl;
	fout0 << setw(12) << "��Լ���ڵ��"
		<< setw(10) << "X��������"
		<< setw(10) << "Y��������"
		<< setw(10) << "R��������" << endl;

	for (i = 0; i < nConstrtainedNode; i++)
	{
		fout0 << setw(12) << (pConsNode + i)->iNode
			<< setw(10) << (pConsNode + i)->iaConstrainedDOF[0]
			<< setw(10) << (pConsNode + i)->iaConstrainedDOF[1]
			<< setw(10) << (pConsNode + i)->iaConstrainedDOF[2] << endl;
	}

	fout0 << endl;

	fout0 << setw(10) << "��Ԫ���ݣ�" << endl << endl;
	fout0 << setw(12) << "��Ԫ�ͺ�" << setw(12) << "ʼ�˽ڵ��" << setw(12) << "�ն˽ڵ��"
		<< setw(12) << "�����ͺ�" << setw(12) << "����������" << endl;
	for (i = 0; i < nTotalElem; i++)
		fout0 << setw(12) << (pElem + i)->iType
		<< setw(12) << (pElem + i)->iaNode[0]
		<< setw(12) << (pElem + i)->iaNode[1]
		<< setw(12) << (pElem + i)->iSection
		<< setw(12) << (pElem + i)->iMaterial << endl;
	fout0 << endl;

	fout0 << setw(10) << "�������ݣ�" << endl << endl;
	fout0 << setw(12) << "����ģ��" << setw(12) << "���ɱ�" << setw(12) << "������ϵ��" << endl;
	for (i = 0; i < nMaterialType; i++)
		fout0 << setw(12) << (pMate + i)->dE
		<< setw(12) << (pMate + i)->dMu
		<< setw(12) << (pMate + i)->dAlpha << endl;
	fout0 << endl;

	cout << setw(14) << "���漸��������" << endl << endl;
	cout << setw(12) << "���" << setw(12) << "���Ծ�" << setw(12) << "�����" << endl;
	for (i = 0; i < nSectionType; i++)
		cout << setw(12) << (pSect + i)->dA
		<< setw(12) << (pSect + i)->dIz
		<< setw(12) << (pSect + i)->dH << endl;
	cout << endl;


	cout << setw(10) << "�غ����ݣ�" << endl << endl;
	cout << setw(10) << "���ͺ�" << setw(10) << "����" << setw(10) << "�غ�ֵ"
		<< setw(10) << "��Ԫ��" << setw(10) << "�ڵ��"
		<< setw(10) << "λ��" << setw(10) << "�¶�T0" << setw(10) << "�¶�T1" << endl;
	for (i = 0; i < nLoad; i++)
		cout << setw(10) << (pLoad + i)->iType
		<< setw(10) << (pLoad + i)->iDirect
		<< setw(10) << (pLoad + i)->dValue
		<< setw(10) << (pLoad + i)->iLoadedElem
		<< setw(10) << (pLoad + i)->iLoadedNode
		<< setw(10) << (pLoad + i)->dPosition
		<< setw(10) << (pLoad + i)->dT0 << setw(10) << (pLoad + i)->dT1 << endl;

	////ԭʼ����������ļ�==============================================================================
	//fout0.setf(ios::right);
	//fout0 << endl << endl;
	//fout0 << setw(14) << "�ܿ����ݣ�" << endl << endl;
	//fout0 << setw(14) << "�ڵ�������" << setw(10) << nTotalNode << endl
	//	<< setw(14) << "��Լ���ڵ�����" << setw(10) << nConstrtainedNode << endl
	//	<< setw(14) << "��Ԫ������" << setw(10) << nTotalElem << endl
	//	<< setw(14) << "��������" << setw(10) << nMaterialType << endl
	//	<< setw(14) << "������������" << setw(10) << nSectionType << endl
	//	<< setw(14) << "�غ�������" << setw(10) << nLoad << endl;
	//fout0 << "=====================================================================" << endl;

	//-------------------------------------------------------
	LengthSinCosCalcu(nTotalElem, pElem, pNode);
	//���������ɶ� �ڵ����ɶȺ͵�Ԫ��λ����
	nTotalDOF = DOFIndexCalcu(nFreeDOF, nTotalNode, nConstrtainedNode, pConsNode, pNode);
	ElementDOFCalcu(nTotalElem, pNode, pElem, pElemDOF);
	//----------------------------------------------------------
	int* pDiag = new int[nTotalDOF];//�����Ԫ��ַ
	BandAndDiagCalcu(nTotalElem, nTotalDOF, pElem, pElemDOF, pDiag);    //�����������Ԫ��ַ
	TwoArrayFree(nTotalElem, pElemDOF);									//�ͷŵ�Ԫ��λ����������ڴ�
	iBuf = pDiag[nTotalDOF - 1] + 1;									//�������Ԫ������
	double* pGK = new double[iBuf];										//һά�������ܸնȾ���������ǲ���
	double* pLoadVect = new double[nTotalDOF];							//����غ�����
	double* pDisp = new double[nTotalDOF];								//���λ������

	VectorZeroize(iBuf, pGK);											//�ܸ�����
	VectorZeroize(nTotalDOF, pLoadVect);								//���غ���������
	VectorZeroize(nTotalDOF, pDisp);									//��λ����������
	ElementEndForceInit(nTotalElem, pElem);								//��ʼ����Ԫ�˶�������
	//װ���ܸնȾ�������غ�����-----------------------------------------------------------------------------

	GKAssembly(nTotalDOF, nTotalElem, pElem, pNode, pMate, pSect, pDiag, pGK);//��װ�ܸ�
	LoadVectorAssembly(nLoad, nTotalDOF, nFreeDOF, pDiag, pGK, pElem, pMate, pSect, pLoad, pNode, pLoadVect, pDisp);//��װ�ܺ�������
	LDLTSolve(nFreeDOF, pDiag, pGK, pDisp);									//��ƽ�ⷽ����λ��
	InternalForceCalcu(nTotalDOF, pElem, pNode, pDisp);						//����˼�����
	SupportReactionCalcu(nTotalDOF, nFreeDOF, pDiag, pGK, pDisp, pLoadVect);	//����֧������
	NodeDisplOutput(fout0, nTotalNode, pNode, pDisp);						//����ڵ�λ��
	EndInternalForceOutput(fout0, nTotalElem, pElem);						//����˼�����
	SupportReactionOutput(fout0, nConstrtainedNode, pConsNode, pNode, pLoadVect);//���֧������

	//�ͷ��ڴ�
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
	//�ر��ļ�;
	fin0.close();
	fout0.close();
	return 0;
}

/**
* \brief ���������
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
 * \brief �����������
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

//����˵�Ԫ���Ⱥͷ�������
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
		dDeltaY = dY1 - dX0;
		(pElem + i)->dLength = sqrt(dDeltaX*dDeltaX + dDeltaY * dDeltaY);
		(pElem + i)->dSin= dDeltaY/(pElem+i)->dLength;
		(pElem + i)->dCos= dDeltaX/ (pElem + i)->dLength;

	}

}

//������ܵ�Ԫ�նȾ���
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
	dBuf = dE * dA / dLength;														//dBufΪ�߸ն�
	pKe[0][0] = pKe[2][2] = dBuf;
	pKe[0][2] = pKe[2][2] = -dBuf;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++) {
			fout1.write((char*)&pKe[i][j], sizeof(double));							//���ֲ�����ϵ���ռ���ת������������ϵ
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

//����ռܵ�Ԫ�նȾ���
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
			fout1.write((char*)&pKe[i][j], sizeof(double));				//���ֲ�����ϵ���ռ���ת������������ϵ
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
 * \brief װ���ܸնȾ���-----------------------------------
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
                double* pGK)
{
	int i, j, m, iNode0, iNode1, GKi, GKj, GKij;
	double** pKe0 = TwoArrayDoubAlloc(4, 4);
	double** pKe1 = TwoArrayDoubAlloc(6, 6);
	int iaDOFIndex[6];

	ofstream fout1("ElemStiff.dat", ios::binary);
	if(!fout1)
	{
		cout << "��Ԫ�նȾ�������ļ���ʧ��!" << endl;
		exit(-1);
	}

	for(m=0;m<nTotalElem;m++)
	{
		switch((pElem+m)->iType)
		{
		case TRUSS:
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
					if(GKi>GKj)
					{
						GKij = pDiag[GKi] - GKi + GKj;
						pGK[GKij] += pKe0[i][j];
					}
				}
			}
			break;
		case FRAME:
			FrameElemStiffCalcu(fout1, pElem + m, pMate, pSect, pKe1);
			iNode0 = (pElem + m)->iaNode[0];
			iNode1 = (pElem + m)->iaNode[1];
			for(i=0;i<3;i++)									//�γɵ�Ԫ��λ����
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
	fout1.close();
}

/**
 * \brief LDLT solve
 * \param nRow 
 * \param pDiag 
 * \param pGK 
 * \param pB 
 * \return 
 */
bool LDLTSolve(int nRow, int * pDiag, double * pGK, double * pB)
{
	int i;
	int j, k;
	int iBegin, jBegin, iAdd_ii, iAdd_jj, iAdd_kk;
	int iRow, nHalfBand;
	double dBuf;

	for(i=1;i<nRow;i++)
	{
		iBegin = i - pDiag[i] + pDiag[i - 1] + 1;
		iAdd_ii = pDiag[i];
		for(j=iBegin+1;j<i+1;j++)
		{
			jBegin = j - pDiag[j] + pDiag[j - 1] + i;
			if (jBegin < iBegin)
				jBegin = iBegin;
			iAdd_jj = pDiag[j];
			dBuf = 0.0;
			for(k=jBegin;k<j;k++)
			{
				if(pGK[pDiag[k]]==0.0)
				{
					cout << "�ܸ�ϵ���д�!" << endl;
					return false;

				}
				dBuf += pGK[iAdd_ii - i + k] * pGK[iAdd_jj - j + k]
					/ pGK[pDiag[k]];

			}
			pGK[iAdd_ii - i + j] -= dBuf;
		}
	}
	if(pGK[0]=0.0)
	{
		cout << "�ܸ�ϵ���д�" << endl;
		return false;

	}
	pB[0] = pB[0] / pGK[0];
	for(i=1;i<nRow;i++)
	{
		dBuf = 0.0;
		iBegin = i - pDiag[i] + pDiag[i - 1] + 1;
		iAdd_ii = pDiag[i];
		for(j=iBegin;j<i;j++)
		{
			dBuf += pB[j] * pGK[iAdd_ii - i + j];
		}
		pB[i] -= dBuf;
		if(pGK[iAdd_ii]==0.0)
		{
			cout << "�ܸ�ϵ���д�" << endl;
			return false;
		}
		pB[i] /= pGK[iAdd_ii];

	}
	nHalfBand = 1;
	iRow = nRow - 1;
	for(i=nRow-2;i>=0;i--)
	{
		nHalfBand++;
		while(nHalfBand>(pDiag[iRow]-pDiag[iRow-1]))
		{
			nHalfBand--;
			iRow--;
		}
		dBuf = 0.0;
		iAdd_ii = pDiag[i];
		for(j=1;j<nHalfBand;j++)
		{
			if(j<(pDiag[j+i]-pDiag[j+i-1]))
			{
				if(pGK[iAdd_ii]==0.0)
				{
					cout << "�ܸ�ϵ���д�!" << endl;
					return false;
				}
				iAdd_kk = pDiag[j + 1];
				dBuf += pGK[iAdd_kk - j] * pB[j + i] / pGK[iAdd_ii];
			}
		}
		pB[i] -= dBuf;
	}
	return true;
}


/**
 * \brief ���㵥Ԫ�̶���
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
	int iBuf;											//�����ɶ���
	for (i = 0; i < nTotalNode; i++)
		for (j = 0; j < 3; j++)
			(pNode + i)->iaDOFIndex[j] = 0;				//�����ڵ����ɶȱ������
	for (i = 0; i < nConstrainedNode; i++) {
		iBuf = (pConsNode + i)->iNode;					//��Լ���Ľڵ��
		for (j = 0; j < 3; j++)
			pNode[iBuf].iaDOFIndex[j] = (pConsNode + i)->iaConstrainedDOF[j];
	}
	iBuf = 0;
	for (i = 0; i < nTotalNode; i++) {
		if ((pNode + i)->iType == FRAME_NODE)
		{
			for (j = 0; j < 3; j++)
				//�Ըּܽڵ��δ֪�������ɶȱ��
				if ((pNode + i)->iaDOFIndex[j] == 0)(pNode + i)->iaDOFIndex[j] = iBuf++;
		}
		else
		{
			for (j = 0; j < 2; j++)
				//����ܽڵ��δ֪�������ɶȱ��
				if ((pNode)->iaDOFIndex[j] == 0) (pNode + i)->iaDOFIndex[j] = iBuf++;

		}

	}
	iBuf0 = iBuf;//δ֪�������ɶ���
	for (i = 0; i < nTotalNode; i++)
	{
		if ((pNode + i)->iType == FRAME_NODE) {
			for (j = 0; j < 3; j++)
				//�Ըռܽڵ����֪�������ɶȱ��
				if ((pNode + i)->iaDOFIndex[j] == -1) (pNode + i)->iaDOFIndex[j] = iBuf++;
		}
		else
		{
			for (j = 0; j < 2; j++)
				//����ܽڵ����֪�������ɶȱ��
				if ((pNode + i)->iaDOFIndex[j] == -1) (pNode + i)->iaDOFIndex[j] = iBuf++;
		}


	}
	for (i = 0; i < nTotalNode; i++)
	{
		if ((pNode + i)->iType == FRAME_NODE)
		{
			for (j = 0; j < 3; j++)
				if ((pNode + i)->iaDOFIndex[j] >= 1e4)//�Ըռܴӽڵ����ɶȱ��
				{
					k = (pNode + i)->iaDOFIndex[j] - 1e4;
					(pNode + i)->iaDOFIndex[j] = (pNode + k)->iaDOFIndex[j];
				}
		}
		else
		{
			for (j = 0; j < 2; j++)
				if ((pNode + i)->iaDOFIndex[j] > 1e4) //����ܴӽڵ����ɶȱ��
				{
					k = (pNode + i)->iaDOFIndex[j] - 1e4;
					(pNode + i)->iaDOFIndex[j] = (pNode + k)->iaDOFIndex[j];
				}
		}
	}
	return iBuf;
}

/**
 * \brief ���㵥Ԫ��λ����
 * \param nTotalElem
 * \param pNode
 * \param pElem
 * \param pElemDOF
 */
void ElementDOFCalcu(int nTotalElem, Node* pNode, Element* pElem, int** pElemDOF)
{
	int iNode0, iNode1; //��Ԫ���˽ڵ��
	int i, j;
	for (i = 0; i < nTotalElem; i++)
	{
		iNode0 = (pElem + i)->iaNode[0];
		iNode1 = (pElem + i)->iaNode[1];
		if ((pElem + i)->iType == TRUSS)				//������ܵ�Ԫ
		{
			for (i = 0; j < 2; j++)
			{
				pElemDOF[i][j] = (pNode + iNode0)->iaDOFIndex[j];
				pElemDOF[i][j + 2] = (pNode + iNode1)->iaDOFIndex[j];
			}
		}
		else //���ڸռܵ�Ԫ
		{
			for (j = 0; j < 3; j++)
			{
				pElemDOF[i][j] = (pNode + iNode0)->iaDOFIndex[i];
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

	for (i = 0; i < nTotalDOF; i++)								//���а����ֵ1
		pDiag[i] = 1;
	for (i = 0; i < nTotalElem; i++)
	{
		iMiniDOF = pElemDOF[i][0];
		if ((pElem + i)->iType == TRUSS) {
			for (j = 0; j < 4; j++)								//����ܵ�Ԫ��4���ڵ�λ�Ʊ����ѡȡ��С��
			{
				if (pElemDOF[i][j] < iMiniDOF)
					iMiniDOF = pElemDOF[i][j];
			}
		}
		else {
			for (j = 0; j < 6; j++) {							//�Ӹռܵ�Ԫ��6���ڵ�λ�Ʊ����ѡ����С��
				if (pElemDOF[i][j] < iMiniDOF)
					iMiniDOF = pElemDOF[i][j];
			}
		}
		if ((pElem + i)->iType == TRUSS) {
			for (j = 0; j < 4; j++) {
				iDOFIndex = pElemDOF[i][j];
				iBuf = iDOFIndex - iMiniDOF + 1;				//��������

				if (iBuf > pDiag[iDOFIndex])
					pDiag[iDOFIndex] = iBuf;
			}
		}
		else {
			for (j = 0; j < 6; j++) {
				iDOFIndex = pElemDOF[i][j];
				iBuf = iDOFIndex - iMiniDOF + 1;				//��������
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
 * \brief ˫���ȶ�ά�����ڴ����
 * \param nRow
 * \param nCol
 * \return
 */
double** TwoArrayDoubAlloc(int nRow, int nCol)
{
	double** pd = new double*[nRow];					//��������
	if (!pd)
	{
		cout << "�ڴ����ʧ�ܣ�" << endl;
		exit(-1);
	}
	for (int i = 0; i < nRow; i++)
	{
		pd[i] = new double[nCol];						//��������
		if (!pd[i])
		{
			cout << "�ڴ����ʧ�ܣ�" << endl;
			exit(-1);
		}
	}
	return pd;
}

/**
* \brief ���Ͷ�ά�����ڴ����
* \param nRow
* \param nCol
* \return
*/
int** TwoArrayIntAlloc(int nRow, int nCol)
{
	int** pd = new int*[nRow];					//��������
	if (!pd)
	{
		cout << "�ڴ����ʧ�ܣ�" << endl;
		exit(-1);
	}
	for (int i = 0; i < nRow; i++)
	{
		pd[i] = new int[nCol];						//��������
		if (!pd[i])
		{
			cout << "�ڴ����ʧ�ܣ�" << endl;
			exit(-1);
		}
	}
	return pd;
}

//��ά�����ڴ��ͷ�
template <class T>
void TwoArrayFree(int nRow, T** pdi)
{
	for (int i = 0; i < nRow; i++)					//�����пռ�
		delete[]pdi[i];
	delete[]pdi;									//�����пռ�
}

//��������
template <class T>
void MatrixZeroize(int nRow, int nCol, T** pT)
{
	for (int i = 0; i < nRow; i++)
		for (int j = 0; j < nCol; j++)
			pT[i][j] = 0;
}

//��������
template <class T>
void VectorZeroize(int n, T* pT)
{
	for (int i = 0; i < n; i++)
		pT[i] = 0;
}

void ElementEndForceInit(int nTotalElem, Element * pElem)
{

}
