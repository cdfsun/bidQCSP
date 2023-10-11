#include <ilcplex/ilocplex.h>
#include  <fstream>
#include  <iostream>

//#define TEST
#define ReduConst// use reduced set in constraints
//#define UserActive//activate the usercallback
#define _AFXDLL
//#define OUTBRANCH//输出分支节点处的解
//#define SubTour

#include "afxtempl.h"
#include  <afx.h>
#include <afxdb.h>
#include <time.h>
#include "stdafx.h"
#include <iomanip>
#include <vector>


//#ifdef _DEBUG
//#define new DEBUG_NEW
//#endif



// 唯一的应用程序对象

CWinApp theApp;

using namespace std;

#undef UNICODE
#undef _UNICODE

const   IloInt nbTask = 50;
const   IloInt nbBay=15;
const   IloInt nbCrane=4;
const   IloInt nbTime=100;

#define instance_no 10//选择测试的算例数量
#define start_instance 1//开始测试的算例


#define QCmove_time 1
#define safe_margin 1
#define multiple_1 1

const   IloInt crane_move_time = 1;

const IloBool ORIGINAL=1;

IloInt SubCutNum;
IloNum ObjVal;
IloNum LB;
IloNum LB0;
IloNum UB;
IloNum GapF;

double  duration_SP = 0;

int cut_count = 0;
int cb_cut_count = 0;
int ST_cb_cut_count = 0;
int user_cut_count = 0;
int prec_cut_count = 0;


int IC_48_cut_count = 0;
int IC_49_cut_count = 0;
int inf_cut_count = 0;
int opt_cut_count = 0;
int opt_cut_count_split = 0;
int SP_count = 0;
int SP_count_split = 0;
int BC_node_count = 0;

int nbLocation_BD[nbTask];
int nbb_BD[nbCrane];



IloNum GapF2;

#define ModNo 0//选择要调用的模型编号, 0：第一种移动方式（一次直，一次回）；1： 第二种

#define CB_cut_Activate 1  //0 不采用组合cut，1 采用
//#define precedence_ineq
//#define re_ineq

//#define split_cut2
//#define strong_prec_cut

#define Normal_detour
#define MPobj

//#define prec_cut
////#define prec_cut2



#define valid_ineq
#define pre_process

#define no_prec_cross

//#define split_cut

//#define LB_Comput


#define waiting_ineq//计算等待时间的有效不等式
#define retracing_ineq
#define retrace_constraints//模型里是否包含等待约束


//#define ceshiUB
//#define ceshi_UB 756

#define lp_tolerance 0.9
#define lp_tolerance2 0.85



typedef IloArray<IloNumArray>    NumMatrix;
typedef IloArray<IloBoolArray>    BoolMatrix;
typedef IloArray<IloArray<IloBoolArray> > BoolMatrix2;
typedef IloArray<IloArray<IloNumArray> > NumMatrix2;
typedef IloArray<IloArray<IloArray<IloBoolArray> > > BoolMatrix3;

typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<IloBoolVarArray>    BoolVarMatrix;
typedef IloArray<IloIntVarArray>    IntVarMatrix;
typedef IloArray<IloArray<IloBoolVarArray> > BoolVarMatrix2;
typedef IloArray<IloArray<IloNumVarArray> > NumVarMatrix2;
typedef IloArray<IloArray<IloArray<IloBoolVarArray> > > BoolVarMatrix3;

ILOSTLBEGIN


//BoolMatrix B_SI;
//IloNumArray  xFsol[nbTime][nbSystem][nbBay];
//IloNumVarArray   uF[nbTime][nbSystem];
//IloNumArray   uFsol[nbTime][nbSystem];
//IloBoolVarArray   zF[nbTime][nbSystem];
//IloBoolVarArray   xF[nbTime][nbSystem];
//IloNumArray   zFsol[nbTime][nbSystem];
//IloBoolArray   B_P_KI[nbMaterial];
//IloBoolArray   B_Q_SK[nbSystem];

BOOL SP_1(IloModel model, IloNum nbs, IloIntArray  nbb, IloNumArray nbQ, IloIntArray  nbLocation, BoolMatrix nbprecR, IloIntArray  nbreadyT,
BoolVarMatrix xF, BoolVarMatrix yF, BoolVarMatrix zF, IloIntVar CF, IloNumVarArray QC_CF, IloIntVarArray thetaF, IloBoolVarArray vF,
BoolMatrix yF_best, BoolMatrix CxF_best, IloNum UB, IloNum *ObjVal, NumMatrix wF_best, NumMatrix CwF_best);

BOOL Chen_Same_Direction(IloModel model, IloNum nbs, IloIntArray  nbb, IloNumArray nbQ, IloIntArray  nbreadyT,
	BoolVarMatrix yF, IloIntVar CF, IloInt *CF_best, BoolMatrix yF_best,
	IloNum *ObjVal, IloNum *UB, IloIntArray  nbLocation, BoolMatrix nbprecR, int my);

BOOL LB_Chen_Same_Direction(IloModel model, IloNum nbs, IloIntArray  nbb, IloNumArray nbQ, IloIntArray  nbreadyT,
	BoolVarMatrix yF, IloIntVar CF, IloInt *CF_best, BoolMatrix yF_best,
	IloNum *ObjVal, IloNum *UB, IloIntArray  nbLocation, BoolMatrix nbprecR, int my);

BOOL LB_Same_Direction(IloModel model, IloNum nbs, IloIntArray  nbb, IloNumArray nbQ, IloIntArray  nbreadyT,
	BoolVarMatrix yF, IloIntVar CF, IloInt *CF_best, BoolMatrix yF_best,
	IloNum *ObjVal, IloNum *UB, IloIntArray  nbLocation, BoolMatrix nbprecR, int my);

BOOL Same_Direction(IloModel model, IloNum nbs, IloIntArray  nbb, IloNumArray nbQ, IloIntArray  nbreadyT,
	BoolVarMatrix yF, IloIntVar CF, IloInt *CF_best, BoolMatrix yF_best,
	IloNum *ObjVal, IloNum *UB, IloIntArray  nbLocation, BoolMatrix nbprecR, int my);

BOOL MP_Same_Direction(IloModel model, IloNum nbs, IloIntArray  nbb, IloNumArray nbQ,
	BoolVarMatrix yF, IloIntVar CF, IloInt *CF_best, BoolMatrix yF_best, NumVarMatrix mF, BoolVarMatrix endCzF, IloNumVarArray betaF,
	IloNum *ObjVal, IloNum *UB, IloIntArray  nbLocation, BoolMatrix nbprecR, NumVarMatrix xF1, NumVarMatrix xF2);

BOOL Split_SP_1(IloModel model, IloNum nbs, IloIntArray  nbb, IloNumArray nbQ, IloIntArray  nbLocation, BoolMatrix nbprecR, IloIntArray  nbreadyT,
	BoolVarMatrix xF, BoolVarMatrix yF, BoolVarMatrix zF, IloIntVar CF, IloNumVarArray QC_CF, IloIntVarArray thetaF, IloBoolVarArray vF,
	BoolMatrix yF_best, BoolMatrix CxF_best, IloNum UB, IloNum *ObjVal, NumMatrix wF_best, NumMatrix CwF_best, int QC);

IloInt Split_SP_Callback(IloInt solmode, IloCplex subcplex, IloNumArray nbQ, int* nbLocation, BoolMatrix nbprecR, IloNumArray2 yF_best, IloInt QC);

ILOUSERCUTCALLBACK7(BendersUserCallback, BoolVarMatrix, yF, NumVarMatrix, xF1, NumVarMatrix, xF2, IloNumArray, nbQ, IloIntArray, nbLocation, BoolMatrix, nbprecR, IloCplex, subcplex)
{
	IloInt i;
	IloEnv masterEnv = getEnv();
	IloInt numNodes = yF.getSize();
	IloInt numNodes_x2 = xF1.getSize();
	IloInt numNodes_y2 = xF2.getSize();

	int split_ind = 0;

	int QC_integer[nbCrane];

	IloNumArray2 yF_best(masterEnv, numNodes);
	for (i = 0; i < numNodes; ++i) {
		//cout << "xF[]: "<<xF[i].getSize() << endl;
		yF_best[i] = IloNumArray(masterEnv, yF[i].getSize());
		getValues(yF_best[i], yF[i]);
	}

	for (int k = 0; k < nbCrane; k++)
	{//检查：但凡取正的yF，是否都逼近整数解
		int che = 0;
		for (i = 0; i < nbTask; i++)
			if (yF_best[k][i] > 0.1)
			{
			if (yF_best[k][i] < 0.9) che++;

			if (yF_best[k][i] >= 0.9)
				yF_best[k][i] = 1;

			}
			else
				yF_best[k][i] = 0;

		if (che==0)
			QC_integer[k] = 1;
		else
			QC_integer[k] = 0;
	}
		
	IloIntArray  QC_ind0(masterEnv, nbCrane);
	IloIntArray  QC_ind(masterEnv, nbCrane);
	for (i = 0; i < nbCrane; ++i)
	{
		QC_ind0[i] = 0;//是否有违背情况，0-否
		QC_ind[i] = 0;//是否求过子问题，0-否
	}
				

	IloNumArray2 xF1_best(masterEnv, numNodes_x2);
	for (i = 0; i < numNodes; ++i) {
		xF1_best[i] = IloNumArray(masterEnv, xF1[i].getSize());
		getValues(xF1_best[i], xF1[i]);
	}

	IloNumArray2 xF2_best(masterEnv, numNodes_y2);
	for (i = 0; i < numNodes; ++i) {
		xF2_best[i] = IloNumArray(masterEnv, xF2[i].getSize());
		getValues(xF2_best[i], xF2[i]);
	}

	
	for (int k = 1; k < nbCrane; k++)
		if (QC_integer[k - 1] == 1 & QC_integer[k] == 1)
			//if (QC_integer[k - 1] == 1)
		{
		int re_ind = 0;
		for (int i2 = 0; i2 < nbTask; i2++)
			if (xF2_best[k - 1][i2] > 0.2)
			{
			re_ind++;
			}
		if (re_ind>=1)
		{
			for (i = 0; i < nbTask; i++)
			{
				if (xF2_best[k][i] > lp_tolerance)//有detour了
				{
					IloInt R0_sum, R1_sum, R2_sum, RT_sum;
					R0_sum = 0;
					RT_sum = 0;
					R1_sum = 0;
					R2_sum = 0;
					IloInt R0_count, RT_count, R1_count, R2_count;
					R0_count = 0;
					RT_count = 0;
					R1_count = 0;
					R2_count = 0;

					IloNum sumx = 0;//计算一下左端项的目前的值的和

					sumx += xF2_best[k][i];

					vector<int> setS1;
					//vector<int> setS2;
					vector<int> setS3;

					IloExpr cbcut0(masterEnv);
					IloExpr cbcut1(masterEnv);
					IloExpr cbcut2(masterEnv);
					IloExpr cbcutT(masterEnv);

					for (int j = 0; j < nbTask; j++)
						if (j != i)
						{
						if (nbLocation[j] <= nbLocation[i] - 2)
						{
							if (xF1_best[k - 1][j] > lp_tolerance)
							{
								R0_sum += nbQ[j];
								R0_count++;
								cbcut0 += xF1[k - 1][j];

								//if (xF1_best[k - 1][j] < lp_tolerance2)
								setS1.push_back(j);

								sumx += xF1_best[k - 1][j];
							}

						}
						if (nbLocation[j] >= nbLocation[i] - 1 && yF_best[k - 1][j] > lp_tolerance && i != j)
						{
							R1_sum += nbQ[j];
							R1_count++;
							cbcut1 += yF[k - 1][j];

							//setS2.push_back(j);

							sumx += yF_best[k - 1][j];
						}
						if (nbLocation[j] <= nbLocation[i] - 1)
						{
							if (xF2_best[k][j] > lp_tolerance)
							{
								RT_sum += nbQ[j];
								RT_count++;
								cbcutT += xF2[k][j];

								setS3.push_back(j);

								sumx += xF2_best[k][j];
							}
						}
						if (nbprecR[i][j] == 1)
						{
							if (yF_best[k][j] > lp_tolerance)
							{
								RT_sum += nbQ[j];
								RT_count++;
								cbcutT += yF[k][j];

								//setS3.push_back(j);

								sumx += yF_best[k][j];
							}
						}
						}

					if (R0_sum + RT_sum + R1_sum + nbQ[i] >= UB)
					{
						QC_ind0[k - 1] = 1;
						add(cbcut0 + xF2[k][i] + cbcut1 + cbcutT <= R0_count + RT_count + R1_count);
						IC_48_cut_count++;
						if (sumx > R0_count + RT_count + R1_count)
						{
							//add(cbcut0 + xF2[k][i] + cbcut1 + cbcutT <= R0_count + RT_count + R1_count);
							QC_ind[k - 1] = 1;
							cb_cut_count++;



							split_ind = 1;
						}

					}

#ifdef strong_prec_cut
					if (setS1.size() > 0.1)
					{

						if (setS1.size() >= 3)
						{
							//排序
							for (int i22 = 0; i22 < setS1.size() - 1; i22++) {
								for (int j22 = 0; j22 < setS1.size() - 1 - i22; j22++) {
									if (xF1_best[k - 1][setS1[j22]] < xF1_best[k - 1][setS1[j22 + 1]]) {        // 相邻元素两两对比
										int temp = setS1[j22 + 1];        // 元素交换
										setS1[j22 + 1] = setS1[j22];
										setS1[j22] = temp;
									}
									else if (xF1_best[k - 1][setS1[j22]] == xF1_best[k - 1][setS1[j22 + 1]])
									{
										if (nbQ[setS1[j22]] < nbQ[setS1[j22 + 1]]) {        // 相邻元素两两对比
											int temp = setS1[j22 + 1];        // 元素交换
											setS1[j22 + 1] = setS1[j22];
											setS1[j22] = temp;
										}
									}
								}
							}

							//删除多个元素

							while (setS1.size() >= 1)
							{
								IloExpr cbcut0_new(masterEnv);
								IloInt R0_sum_new = 0;
								for (int j2 = 0; j2 < setS1.size(); j2++)
								{
									cbcut0_new += xF1[k - 1][setS1[j2]];
									R0_sum_new += nbQ[setS1[j2]];
								}

								if (R0_sum_new + RT_sum + R1_sum + nbQ[i] - nbQ[setS1[setS1.size() - 1]] >= UB)
								{
									setS1.pop_back();
									//cout << "set S1 no: " << setS1.size()<< endl;
								}
								else
								{
									//cout << "strong cut jiale!:: " << endl;
									int count_r0;
									count_r0 = setS1.size();
									add(cbcut0_new + xF2[k][i] + cbcut1 + cbcutT <= count_r0 + RT_count + R1_count);
									IC_48_cut_count++;
									cbcut0_new.end();

									break;
								}
								cbcut0_new.end();
							}

						}
						else
						{
							for (int j = 0; j < setS1.size(); j++)
							{
								if (R0_sum + RT_sum + R1_sum + nbQ[i] - nbQ[setS1[j]] >= UB)
								{
									add(cbcut0 + xF2[k][i] + cbcut1 + cbcutT - xF1[k - 1][setS1[j]] <= R0_count + RT_count + R1_count - 1);
									ST_cb_cut_count++;
									IC_48_cut_count++;

									//if (sumx - xF1_best[k - 1][setS1[j]]> R0_count + RT_count + R1_count - 1)
									//{
									//	split_ind = 1;
									//}
								}
							}
						}
					}

					if (setS3.size() > 0.1)
					{
						if (setS3.size() >= 3)
						{
							//排序
							for (int i22 = 0; i22 < setS3.size() - 1; i22++) {
								for (int j22 = 0; j22 < setS3.size() - 1 - i22; j22++) {
									if (xF2_best[k][setS3[j22]] < xF2_best[k][setS3[j22 + 1]]) {        // 相邻元素两两对比
										int temp = setS3[j22 + 1];        // 元素交换
										setS3[j22 + 1] = setS3[j22];
										setS3[j22] = temp;
									}
									else if (xF2_best[k][setS3[j22]] == xF2_best[k][setS3[j22 + 1]])
									{
										if (nbQ[setS3[j22]] < nbQ[setS3[j22 + 1]]) {        // 相邻元素两两对比
											int temp = setS3[j22 + 1];        // 元素交换
											setS3[j22 + 1] = setS3[j22];
											setS3[j22] = temp;
										}
									}
								}
							}

							//删除多个元素

							while (setS3.size() >= 1)
							{
								IloExpr cbcutT_new(masterEnv);
								IloInt RT_sum_new = 0;
								for (int j2 = 0; j2 < setS3.size(); j2++)
								{
									cbcutT_new += xF2[k][setS3[j2]];
									RT_sum_new += nbQ[setS3[j2]];
								}

								if (R0_sum + RT_sum_new + R1_sum + nbQ[i] - nbQ[setS3[setS3.size() - 1]] >= UB)
								{
									setS3.pop_back();
									//cout << "set S1 no: " << setS1.size()<< endl;
								}
								else
								{
									//cout << "strong cut jiale!:: " << endl;
									int count_rT;
									count_rT = setS3.size();
									add(cbcut0 + xF2[k][i] + cbcut1 + cbcutT_new <= count_rT + R0_count + R1_count);
									IC_49_cut_count++;
									cbcutT_new.end();

									break;
								}
								cbcutT_new.end();
							}

						}
						else
						{
							for (int j = 0; j < setS3.size(); j++)
							{
								if (R0_sum + RT_sum + R1_sum + nbQ[i] - nbQ[setS3[j]] >= UB)
								{
									add(cbcut0 + xF2[k][i] + cbcut1 + cbcutT - xF2[k][setS3[j]] <= R0_count + RT_count + R1_count - 1);
									ST_cb_cut_count++;
									IC_48_cut_count++;

									if (sumx - xF2_best[k][setS3[j]] > R0_count + RT_count + R1_count - 1)
									{
										split_ind = 1;
									}
								}
							}
						}
					}
#endif
					cbcut0.end();
					cbcut1.end();
					cbcut2.end();
					cbcutT.end();


					vector<int>().swap(setS1);
					//vector<int>().swap(setS2);
					vector<int>().swap(setS3);


				}

			else if (xF1_best[k][i] > lp_tolerance)
			{
				IloInt R0_sum, R1_sum, R2_sum, RT_sum;
				R0_sum = 0;
				RT_sum = 0;
				R1_sum = 0;
				R2_sum = 0;
				IloInt R0_count, RT_count, R1_count, R2_count;
				R0_count = 0;
				RT_count = 0;
				R1_count = 0;
				R2_count = 0;

				vector<int> setS1;
				//vector<int> setS2;
				vector<int> setS3;

				IloExpr cbcut0(masterEnv);
				IloExpr cbcut1(masterEnv);
				IloExpr cbcut2(masterEnv);
				IloExpr cbcutT(masterEnv);

				IloNum sumx = 0;//计算一下左端项的目前的值的和

				sumx += xF1_best[k][i];

				for (int j = 0; j < nbTask; j++)
					if (j != i)
					{
						if (nbLocation[j] <= nbLocation[i] - 1)
						{
							if (xF1_best[k][j] > lp_tolerance)
							{
								R0_sum += nbQ[j];
								R0_count++;
								cbcut0 += xF1[k][j];

								//if (xF1_best[k][j] < lp_tolerance2)
									setS1.push_back(j);

								sumx += xF1_best[k][j];
							}

						}
						if (nbprecR[j][i] == 1)
						{
							if (yF_best[k][j] > lp_tolerance)
							{
								R2_sum += nbQ[j];
								R2_count++;
								cbcut2 += yF[k][j];

								//setS3.push_back(j);

								sumx += yF_best[k][j];
							}
						}
						if (nbLocation[j] <= nbLocation[i] - 2)
						{

							if (xF2_best[k - 1][j] > lp_tolerance)
							{
								RT_sum += nbQ[j];
								RT_count++;
								cbcutT += xF2[k - 1][j];

								//if (xF2_best[k - 1][j] < lp_tolerance2)
									setS3.push_back(j);

								sumx += xF2_best[k - 1][j];
							}

						}
						if (nbLocation[j] >= nbLocation[i] - 1 && yF_best[k - 1][j] > lp_tolerance)
						{
							R1_sum += nbQ[j];
							R1_count++;
							cbcut1 += yF[k - 1][j];

							//setS2.push_back(j);

							sumx += yF_best[k - 1][j];
						}

					}

				if (R0_sum + RT_sum + R1_sum + R2_sum + nbQ[i] >= UB)
				{
					QC_ind0[k - 1] = 1;
					add(cbcut0 + xF1[k][i] + cbcut2 + cbcutT <= R0_count + RT_count + R1_count + R2_count);
					IC_49_cut_count++;

					//if (sumx>R0_count + RT_count + R1_count + R2_count)
					//{
					//	//add(cbcut0 + xF1[k][i] + cbcut2 + cbcutT <= R0_count + RT_count + R2_count);
					//	QC_ind[k - 1] =1;
					//	split_ind = 1;

					//	//cb_cut_count++;
					//}

				}

#ifdef strong_prec_cut
				if (setS1.size() > 0.1)
				{

					if (setS1.size() >= 3)
					{
						//排序
						for (int i22 = 0; i22 < setS1.size() - 1; i22++) {
							for (int j22 = 0; j22 < setS1.size() - 1 - i22; j22++) {
								if (xF1_best[k - 1][setS1[j22]] < xF1_best[k - 1][setS1[j22 + 1]]) {        // 相邻元素两两对比
									int temp = setS1[j22 + 1];        // 元素交换
									setS1[j22 + 1] = setS1[j22];
									setS1[j22] = temp;
								}
								else if (xF1_best[k - 1][setS1[j22]] == xF1_best[k - 1][setS1[j22 + 1]])
								{
									if (nbQ[setS1[j22]] < nbQ[setS1[j22 + 1]]) {        // 相邻元素两两对比
										int temp = setS1[j22 + 1];        // 元素交换
										setS1[j22 + 1] = setS1[j22];
										setS1[j22] = temp;
									}
								}
							}
						}

						//删除多个元素

						while (setS1.size() >= 1)
						{
							IloExpr cbcut0_new(masterEnv);
							IloInt R0_sum_new = 0;
							for (int j2 = 0; j2 < setS1.size(); j2++)
							{
								cbcut0_new += xF1[k - 1][setS1[j2]];
								R0_sum_new += nbQ[setS1[j2]];
							}

							if (R0_sum_new + RT_sum + R1_sum + R2_sum + nbQ[i] - nbQ[setS1[setS1.size() - 1]] >= UB)
							{
								setS1.pop_back();
								//cout << "set S1 no: " << setS1.size()<< endl;
							}
							else
							{
								//cout << "strong cut jiale!:: " << endl;
								int count_r0;
								count_r0 = setS1.size();
								add(cbcut0_new + xF1[k][i] + cbcut1 + cbcut2 + cbcutT <= count_r0 + RT_count + R1_count + R2_count);
								IC_49_cut_count++;
								cbcut0_new.end();

								break;
							}
							cbcut0_new.end();
						}

					}
					else
					{
						for (int j = 0; j < setS1.size(); j++)
						{
							if (R0_sum + RT_sum + R1_sum + nbQ[i] - nbQ[setS1[j]] >= UB)
							{
								add(cbcut0 + xF1[k][i] + cbcut1 + cbcut2 + cbcutT - xF1[k - 1][setS1[j]] <= R0_count + RT_count + R1_count + R2_count - 1);
								ST_cb_cut_count++;
								IC_49_cut_count++;

								//if (sumx - xF1_best[k - 1][setS1[j]]> R0_count + RT_count + R1_count - 1)
								//{
								//	split_ind = 1;
								//}
							}
						}
					}
				}

				if (setS3.size() > 0.1)
				{
					if (setS3.size() >= 3)
					{
						//排序
						for (int i22 = 0; i22 < setS3.size() - 1; i22++) {
							for (int j22 = 0; j22 < setS3.size() - 1 - i22; j22++) {
								if (xF2_best[k][setS3[j22]] < xF2_best[k][setS3[j22 + 1]]) {        // 相邻元素两两对比
									int temp = setS3[j22 + 1];        // 元素交换
									setS3[j22 + 1] = setS3[j22];
									setS3[j22] = temp;
								}
								else if (xF2_best[k][setS3[j22]] == xF2_best[k][setS3[j22 + 1]])
								{
									if (nbQ[setS3[j22]] < nbQ[setS3[j22 + 1]]) {        // 相邻元素两两对比
										int temp = setS3[j22 + 1];        // 元素交换
										setS3[j22 + 1] = setS3[j22];
										setS3[j22] = temp;
									}
								}
							}
						}

						//删除多个元素

						while (setS3.size() >= 1)
						{
							IloExpr cbcutT_new(masterEnv);
							IloInt RT_sum_new = 0;
							for (int j2 = 0; j2 < setS3.size(); j2++)
							{
								cbcutT_new += xF2[k][setS3[j2]];
								RT_sum_new += nbQ[setS3[j2]];
							}

							if (R0_sum + RT_sum_new + R1_sum + R2_sum + nbQ[i] - nbQ[setS3[setS3.size() - 1]] >= UB)
							{
								setS3.pop_back();
								//cout << "set S1 no: " << setS1.size()<< endl;
							}
							else
							{
								//cout << "strong cut jiale!:: " << endl;
								int count_rT;
								count_rT = setS3.size();
								add(cbcut0 + xF1[k][i] + cbcut1 + cbcut2 + cbcutT_new <= count_rT + R0_count + R1_count + R2_count);
								IC_49_cut_count++;
								cbcutT_new.end();

								break;
							}
							cbcutT_new.end();
						}

					}
					else
					{
						for (int j = 0; j < setS3.size(); j++)
						{
							if (R0_sum + RT_sum + R1_sum + R2_sum + nbQ[i] - nbQ[setS3[j]] >= UB)
							{
								add(cbcut0 + xF1[k][i] + cbcut1 + cbcut2 + cbcutT - xF2[k][setS3[j]] <= R0_count + RT_count + R1_count + R2_count - 1);
								ST_cb_cut_count++;
								IC_49_cut_count++;

								//if (sumx - xF2_best[k][setS3[j]] > R0_count + RT_count + R1_count + R2_count - 1)
								//{
								//	split_ind = 1;
								//}
							}
						}
					}
				}
#endif	
				cbcut0.end();
				cbcut1.end();
				cbcut2.end();
				cbcutT.end();

				vector<int>().swap(setS1);
				//vector<int>().swap(setS2);
				vector<int>().swap(setS3);


			}
			}




			//if (QC_ind[k - 1] == 0 && QC_ind0[k - 1] == 1 && split_ind == 0)
			//{
			//	clock_t start_SP = 0, finish_SP = 0;
			//	start_SP = clock();
			//	IloInt obj;
			//	obj = Split_SP_Callback(1, subcplex, nbQ, nbLocation, nbprecR, yF_best, k - 1);
			//	finish_SP = clock();
			//	duration_SP += (double)(finish_SP - start_SP) / CLOCKS_PER_SEC;

			//	if (obj > UB || obj == 0)
			//	{
			//		IloExpr BDcut(masterEnv);

			//		for (int k2 = k - 1; k2 <= k; k2++)
			//		{
			//			for (int i2 = 0; i2 < nbTask; i2++)
			//			{
			//				if (yF_best[k2][i2] > 0.9)
			//					BDcut += 1 - yF[k2][i2];
			//			}
			//		}
			//		add(BDcut >= 1);
			//		BDcut.end();

			//		cut_count++;
			//	}
			//}

		}
		}

}

ILOLAZYCONSTRAINTCALLBACK7(BendersLazyCallback, BoolVarMatrix, yF, NumVarMatrix, xF1, NumVarMatrix, xF2, IloNumArray, nbQ, IloIntVar, CF, BoolMatrix, nbprecR,
	IloCplex, subcplex)
{
	IloInt i;
	IloEnv masterEnv = getEnv();
	IloInt numNodes = yF.getSize();
	IloInt numNodes_x2 = xF1.getSize();
	IloInt numNodes_y2 = xF2.getSize();

	IloIntArray  QC_ind0(masterEnv, nbCrane);
	IloIntArray  QC_ind(masterEnv, nbCrane);

	for (i = 0; i < nbCrane; ++i)
	{
		QC_ind0[i] = 0;//是否有违背情况，0-否
		QC_ind[i] = 0;//是否求过子问题，0-否
	}

	IloNumArray2 yF_best(masterEnv, numNodes);
	for (i = 0; i < numNodes; ++i) {
		//cout << "xF[]: "<<xF[i].getSize() << endl;
		yF_best[i] = IloNumArray(masterEnv, yF[i].getSize());
		getValues(yF_best[i], yF[i]);
	}

	IloNumArray2 xF1_best(masterEnv, numNodes_x2);
	for (i = 0; i < numNodes; ++i) {
		xF1_best[i] = IloNumArray(masterEnv, xF1[i].getSize());
		getValues(xF1_best[i], xF1[i]);
	}

	IloNumArray2 xF2_best(masterEnv, numNodes_y2);
	for (i = 0; i < numNodes; ++i) {
		xF2_best[i] = IloNumArray(masterEnv, xF2[i].getSize());
		getValues(xF2_best[i], xF2[i]);
	}

//	cout << "lazy: " << endl;
	//for (int k = 0; k < nbCrane; k++)
	//{
	//	cout << k << ":  ";

	//	for (i = 0; i < nbTask; i++)
	//	{

	//		if (xF1_best[k][i]> 0.01)
	//			cout << i << "-f-" << xF1_best[k][i] << "  ";
	//		if (xF2_best[k][i]> 0.01)
	//			cout << i << "-d-" << xF2_best[k][i] << "  ";

	//	}
	//	cout << endl;
	//}
	//cout << endl << endl;

	////////////////////////////////////////////////////////////////////////
	////构造combinatorial cut
	////////////////////////////////////////////////////////////////////////

	int split_ind = 0;


	SP_count++;

	for (int k = 0; k < nbCrane; k++)
		for (i = 0; i < nbTask; i++)
			if (yF_best[k][i]>0.9)
				yF_best[k][i] = 1;

	clock_t start_SP = 0, finish_SP = 0;
	start_SP = clock();

	//if (split_ind == 0)
	{
		IloInt obj;
		obj = Split_SP_Callback(0, subcplex, nbQ, nbLocation_BD, nbprecR, yF_best, 0);
		if (obj == 0)
		{
			IloExpr BDcut(masterEnv);

			for (int k2 = 0; k2 < nbCrane; k2++)
			{
				for (int i2 = 0; i2 < nbTask; i2++)
				{
					if (yF_best[k2][i2] > 0.9)
						BDcut += 1 - yF[k2][i2];
				}
			}
			add(BDcut >= 1);
			BDcut.end();

			inf_cut_count++;

			cut_count++;
		}
		else
		{
			if (obj < UB)
			{
				split_ind = 1;
				UB = obj;
			}
			//else
			{
				IloExpr BDcut(masterEnv);

				for (int k2 = 0; k2 < nbCrane; k2++)
				{
					for (int i2 = 0; i2 < nbTask; i2++)
					{
						if (yF_best[k2][i2] > 0.9)
							BDcut += (obj - LB0)*(1 - yF[k2][i2]);
					}
				}
				add(CF - obj + BDcut >= 0);
				BDcut.end();

				opt_cut_count++;

				cut_count++;
			}




		}



	}

	finish_SP = clock();
	duration_SP += (double)(finish_SP - start_SP) / CLOCKS_PER_SEC;


	////////////////////////////////////////////////////////////////////////
	////构造combinatorial cut
	////////////////////////////////////////////////////////////////////////

	//split_ind = 1;
	if (CB_cut_Activate)
	{

	for (int k = 1; k < nbCrane; k++)
	{
		for (i = 0; i < nbTask; i++)
		{
			if (xF2_best[k][i] > lp_tolerance)//有detour了
			{
				IloInt R0_sum, R1_sum, R2_sum, RT_sum;
				R0_sum = 0;
				RT_sum = 0;
				R1_sum = 0;
				R2_sum = 0;
				IloInt R0_count, RT_count, R1_count, R2_count;
				R0_count = 0;
				RT_count = 0;
				R1_count = 0;
				R2_count = 0;

				IloNum sumx = 0;//计算一下左端项的目前的值的和

				sumx += xF2_best[k][i];

				vector<int> setS1;
				//vector<int> setS2;
				vector<int> setS3;

				IloExpr cbcut0(masterEnv);
				IloExpr cbcut1(masterEnv);
				IloExpr cbcut2(masterEnv);
				IloExpr cbcutT(masterEnv);

				for (int j = 0; j < nbTask; j++)
					if (j != i)
				{
					if (nbLocation_BD[j] <= nbLocation_BD[i] - 2)
					{
						if (xF1_best[k - 1][j] > lp_tolerance)
						{
							R0_sum += nbQ[j];
							R0_count++;
							cbcut0 += xF1[k - 1][j];

							//if (xF1_best[k - 1][j] < lp_tolerance2)
							setS1.push_back(j);

							sumx += xF1_best[k - 1][j];
						}

					}
					if (nbLocation_BD[j] >= nbLocation_BD[i] - 1 && yF_best[k - 1][j] > lp_tolerance && i != j)
					{
						R1_sum += nbQ[j];
						R1_count++;
						cbcut1 += yF[k - 1][j];

						//setS2.push_back(j);

						sumx += yF_best[k - 1][j];
					}
					if (nbLocation_BD[j] <= nbLocation_BD[i] - 1)
					{
						if (xF2_best[k][j] > lp_tolerance)
						{
							RT_sum += nbQ[j];
							RT_count++;
							cbcutT += xF2[k][j];

							setS3.push_back(j);

							sumx += xF2_best[k][j];
						}
					}
					if (nbprecR[i][j] == 1)
					{
						if (yF_best[k][j] > lp_tolerance)
						{
							R2_sum += nbQ[j];
							R2_count++;
							cbcut2 += yF[k][j];

							//setS3.push_back(j);

							sumx += yF_best[k][j];
						}
					}
				}

				if (R0_sum + RT_sum + R1_sum + R2_sum + nbQ[i] >= UB)
				{
					QC_ind0[k - 1] = 1;
					add(cbcut0 + xF2[k][i] + cbcut1 + cbcut2 + cbcutT <= R0_count + RT_count + R1_count + R2_count);
					IC_48_cut_count++;
					//cout << "haha1" << endl;
					//if (sumx > R0_count + RT_count + R1_count)
					//{
					//	//add(cbcut0 + xF2[k][i] + cbcut1 + cbcutT <= R0_count + RT_count + R1_count);
					//	QC_ind[k - 1] =1;
					//	split_ind = 1;
					//	//cb_cut_count++;
					//	//IC_48_cut_count++;
					//}				

				}

#ifdef strong_prec_cut
				if (setS1.size() > 0.1)
				{

					if (setS1.size() >= 3)
					{
						//排序
						for (int i22 = 0; i22 < setS1.size() - 1; i22++) {
							for (int j22 = 0; j22 < setS1.size() - 1 - i22; j22++) {
								if (xF1_best[k - 1][setS1[j22]] < xF1_best[k - 1][setS1[j22+1]]) {        // 相邻元素两两对比
									int temp = setS1[j22+1];        // 元素交换
									setS1[j22 + 1] = setS1[j22];
									setS1[j22] = temp;
								}
								else if (xF1_best[k - 1][setS1[j22]] == xF1_best[k - 1][setS1[j22 + 1]])
								{
									if (nbQ[setS1[j22]] < nbQ[setS1[j22 + 1]]) {        // 相邻元素两两对比
										int temp = setS1[j22 + 1];        // 元素交换
										setS1[j22 + 1] = setS1[j22];
										setS1[j22] = temp;
									}
								}
							}
						}

						//删除多个元素
						
						while (setS1.size()>=1)
						{
							IloExpr cbcut0_new(masterEnv);
							IloInt R0_sum_new=0;
							for (int j2 = 0; j2 < setS1.size(); j2++)
							{
								cbcut0_new += xF1[k - 1][setS1[j2]];
								R0_sum_new += nbQ[setS1[j2]];
							}				
							
							if (R0_sum_new + RT_sum + R1_sum + R2_sum + nbQ[i] - nbQ[setS1[setS1.size() - 1]] >= UB)
							{
								setS1.pop_back();
								//setS1.erase(setS1.end());
								//cout << "set S1 no: " << setS1.size()<< endl;
							}
							else
							{
								//cout << "strong cut jiale!:: " << endl;
								int count_r0;
								count_r0 = setS1.size();
								add(cbcut0_new + xF2[k][i] + cbcut1 + cbcut2 + cbcutT <= count_r0 + RT_count + R1_count + R2_count);
								//add(cbcut0_new + xF2[k][i] + cbcut1 + cbcut2 + cbcutT <= R0_count + RT_count + R1_count + R2_count);
								IC_48_cut_count++;
								cbcut0_new.end();

								break;
							}
							cbcut0_new.end();
						}

					}
					else
					{
						for (int j = 0; j < setS1.size(); j++)
						{
							if (R0_sum + RT_sum + R1_sum + nbQ[i] - nbQ[setS1[j]] >= UB)
							{
								add(cbcut0 + xF2[k][i] + cbcut1 + cbcut2 + cbcutT - xF1[k - 1][setS1[j]] <= R0_count + RT_count + R1_count + R2_count - 1);
								ST_cb_cut_count++;
								IC_48_cut_count++;

								//if (sumx - xF1_best[k - 1][setS1[j]]> R0_count + RT_count + R1_count - 1)
								//{
								//	split_ind = 1;
								//}
							}
						}
					}
				}

					if (setS3.size() > 0.1)
					{
						if (setS3.size() >= 3)
						{
							//排序
							for (int i22 = 0; i22 < setS3.size() - 1; i22++) {
								for (int j22 = 0; j22 < setS3.size() - 1 - i22; j22++) {
									if (xF2_best[k][setS3[j22]] < xF2_best[k][setS3[j22 + 1]]) {        // 相邻元素两两对比
										int temp = setS3[j22 + 1];        // 元素交换
										setS3[j22 + 1] = setS3[j22];
										setS3[j22] = temp;
									}
									else if (xF2_best[k][setS3[j22]] == xF2_best[k][setS3[j22 + 1]])
									{
										if (nbQ[setS3[j22]] < nbQ[setS3[j22 + 1]]) {        // 相邻元素两两对比
											int temp = setS3[j22 + 1];        // 元素交换
											setS3[j22 + 1] = setS3[j22];
											setS3[j22] = temp;
										}
									}
								}
							}

							//删除多个元素

							while (setS3.size() >= 1)
							{
								IloExpr cbcutT_new(masterEnv);
								IloInt RT_sum_new = 0;
								for (int j2 = 0; j2 < setS3.size(); j2++)
								{
									cbcutT_new += xF2[k][setS3[j2]];
									RT_sum_new += nbQ[setS3[j2]];
								}

								if (R0_sum + RT_sum_new + R1_sum + R2_sum + nbQ[i] - nbQ[setS3[setS3.size() - 1]] >= UB)
								{
									setS3.pop_back();
									//cout << "set S1 no: " << setS1.size()<< endl;
								}
								else
								{
									//cout << "strong cut jiale!:: " << endl;
									int count_rT;
									count_rT = setS3.size();
									add(cbcut0 + xF2[k][i] + cbcut1 + cbcut2 + cbcutT_new <= count_rT + R0_count + R1_count + R2_count);
									IC_48_cut_count++;
									cbcutT_new.end();

									break;
								}
								cbcutT_new.end();
							}

						}
						else
						{
							for (int j = 0; j < setS3.size(); j++)
							{
								if (R0_sum + RT_sum + R1_sum + R2_sum + nbQ[i] - nbQ[setS3[j]] >= UB)
								{
									add(cbcut0 + xF2[k][i] + cbcut1 + cbcut2 + cbcutT - xF2[k][setS3[j]] <= R0_count + RT_count + R1_count + R2_count - 1);
									ST_cb_cut_count++;
									IC_48_cut_count++;

									//if (sumx - xF2_best[k][setS3[j]] > R0_count + RT_count + R1_count + R2_count - 1)
									//{
									//	split_ind = 1;
									//}
								}
							}
						}
					}
#endif
				cbcut0.end();
				cbcut1.end();
				cbcut2.end();
				cbcutT.end();


				vector<int>().swap(setS1);
				//vector<int>().swap(setS2);
				vector<int>().swap(setS3);


			}

			else if (xF1_best[k][i] > lp_tolerance)
			{
				IloInt R0_sum, R1_sum, R2_sum, RT_sum;
				R0_sum = 0;
				RT_sum = 0;
				R1_sum = 0;
				R2_sum = 0;
				IloInt R0_count, RT_count, R1_count, R2_count;
				R0_count = 0;
				RT_count = 0;
				R1_count = 0;
				R2_count = 0;

				vector<int> setS1;
				//vector<int> setS2;
				vector<int> setS3;

				IloExpr cbcut0(masterEnv);
				IloExpr cbcut1(masterEnv);
				IloExpr cbcut2(masterEnv);
				IloExpr cbcutT(masterEnv);

				IloNum sumx = 0;//计算一下左端项的目前的值的和

				sumx += xF1_best[k][i];

				for (int j = 0; j < nbTask; j++)
					if (j!=i)
				{
					if (nbLocation_BD[j] <= nbLocation_BD[i] - 1)
					{
						if (xF1_best[k][j] > lp_tolerance)
						{
							R0_sum += nbQ[j];
							R0_count++;
							cbcut0 += xF1[k][j];

							//if (xF1_best[k][j] < lp_tolerance2)
							setS1.push_back(j);

							sumx += xF1_best[k][j];
						}

					}
					if (nbprecR[j][i] == 1)
					{
						if (yF_best[k][j] > lp_tolerance)
						{
							R2_sum += nbQ[j];
							R2_count++;
							cbcut2 += yF[k][j];

							//setS3.push_back(j);

							sumx += yF_best[k][j];
						}
					}
					if (nbLocation_BD[j] <= nbLocation_BD[i] - 2)
					{

						if (xF2_best[k - 1][j] > lp_tolerance)
						{
							RT_sum += nbQ[j];
							RT_count++;
							cbcutT += xF2[k - 1][j];

							//if (xF2_best[k - 1][j] < lp_tolerance2)
							setS3.push_back(j);

							sumx += xF2_best[k - 1][j];
						}

					}
					if (nbLocation_BD[j] >= nbLocation_BD[i] - 1 && yF_best[k - 1][j] > lp_tolerance)
					{
						R1_sum += nbQ[j];
						R1_count++;
						cbcut1 += yF[k-1][j];

						//setS2.push_back(j);

						sumx += yF_best[k-1][j];
					}
	
				}

				if (R0_sum + RT_sum + R1_sum + R2_sum + nbQ[i] >= UB)
				{
					QC_ind0[k-1] = 1;
					add(cbcut0 + xF1[k][i] + cbcut2 + cbcutT <= R0_count + RT_count + R1_count + R2_count);
					IC_49_cut_count++; 

					//cout << "haha2" << endl;

					//if (sumx>R0_count + RT_count + R1_count + R2_count)
					//{
					//	//add(cbcut0 + xF1[k][i] + cbcut2 + cbcutT <= R0_count + RT_count + R2_count);
					//	QC_ind[k - 1] =1;
					//	split_ind = 1;

					//	//cb_cut_count++;
					//}

				}

#ifdef strong_prec_cut
				if (setS1.size() > 0.1)
				{

					if (setS1.size() >= 3)
					{
						//排序
						for (int i22 = 0; i22 < setS1.size() - 1; i22++) {
							for (int j22 = 0; j22 < setS1.size() - 1 - i22; j22++) {
								if (xF1_best[k - 1][setS1[j22]] < xF1_best[k - 1][setS1[j22 + 1]]) {        // 相邻元素两两对比
									int temp = setS1[j22 + 1];        // 元素交换
									setS1[j22 + 1] = setS1[j22];
									setS1[j22] = temp;
								}
								else if (xF1_best[k - 1][setS1[j22]] == xF1_best[k - 1][setS1[j22 + 1]])
								{
									if (nbQ[setS1[j22]] < nbQ[setS1[j22 + 1]]) {        // 相邻元素两两对比
										int temp = setS1[j22 + 1];        // 元素交换
										setS1[j22 + 1] = setS1[j22];
										setS1[j22] = temp;
									}
								}
							}
						}

						//删除多个元素

						while (setS1.size() >= 1)
						{
							IloExpr cbcut0_new(masterEnv);
							IloInt R0_sum_new = 0;
							for (int j2 = 0; j2 < setS1.size(); j2++)
							{
								cbcut0_new += xF1[k - 1][setS1[j2]];
								R0_sum_new += nbQ[setS1[j2]];
							}

							if (R0_sum_new + RT_sum + R1_sum + R2_sum + nbQ[i] - nbQ[setS1[setS1.size() - 1]] >= UB)
							{
								setS1.pop_back();
								//cout << "set S1 no: " << setS1.size()<< endl;
							}
							else
							{
								//cout << "strong cut jiale!:: " << endl;
								int count_r0;
								count_r0 = setS1.size();
								add(cbcut0_new + xF1[k][i] + cbcut1 + cbcut2 + cbcutT <= count_r0 + RT_count + R1_count + R2_count);
								IC_49_cut_count++;
								cbcut0_new.end();

								break;
							}
							cbcut0_new.end();
						}

					}
					else
					{
						for (int j = 0; j < setS1.size(); j++)
						{
							if (R0_sum + RT_sum + R1_sum + nbQ[i] - nbQ[setS1[j]] >= UB)
							{
								add(cbcut0 + xF1[k][i] + cbcut1 + cbcut2 + cbcutT - xF1[k - 1][setS1[j]] <= R0_count + RT_count + R1_count + R2_count - 1);
								ST_cb_cut_count++;
								IC_49_cut_count++;

								//if (sumx - xF1_best[k - 1][setS1[j]]> R0_count + RT_count + R1_count - 1)
								//{
								//	split_ind = 1;
								//}
							}
						}
					}
				}

				if (setS3.size() > 0.1)
				{
					if (setS3.size() >= 3)
					{
						//排序
						for (int i22 = 0; i22 < setS3.size() - 1; i22++) {
							for (int j22 = 0; j22 < setS3.size() - 1 - i22; j22++) {
								if (xF2_best[k][setS3[j22]] < xF2_best[k][setS3[j22 + 1]]) {        // 相邻元素两两对比
									int temp = setS3[j22 + 1];        // 元素交换
									setS3[j22 + 1] = setS3[j22];
									setS3[j22] = temp;
								}
								else if (xF2_best[k][setS3[j22]] == xF2_best[k][setS3[j22 + 1]])
								{
									if (nbQ[setS3[j22]] < nbQ[setS3[j22 + 1]]) {        // 相邻元素两两对比
										int temp = setS3[j22 + 1];        // 元素交换
										setS3[j22 + 1] = setS3[j22];
										setS3[j22] = temp;
									}
								}
							}
						}

						//删除多个元素

						while (setS3.size() >= 1)
						{
							IloExpr cbcutT_new(masterEnv);
							IloInt RT_sum_new = 0;
							for (int j2 = 0; j2 < setS3.size(); j2++)
							{
								cbcutT_new += xF2[k][setS3[j2]];
								RT_sum_new += nbQ[setS3[j2]];
							}

							if (R0_sum + RT_sum_new + R1_sum + R2_sum + nbQ[i] - nbQ[setS3[setS3.size() - 1]] >= UB)
							{
								setS3.pop_back();
								//cout << "set S1 no: " << setS1.size()<< endl;
							}
							else
							{
								//cout << "strong cut jiale!:: " << endl;
								int count_rT;
								count_rT = setS3.size();
								add(cbcut0 + xF1[k][i] + cbcut1 + cbcut2 + cbcutT_new <= count_rT + R0_count + R1_count + R2_count);
								IC_49_cut_count++;
								cbcutT_new.end();

								break;
							}
							cbcutT_new.end();
						}

					}
					else
					{
						for (int j = 0; j < setS3.size(); j++)
						{
							if (R0_sum + RT_sum + R1_sum + R2_sum + nbQ[i] - nbQ[setS3[j]] >= UB)
							{
								add(cbcut0 + xF1[k][i] + cbcut1 + cbcut2 + cbcutT - xF2[k][setS3[j]] <= R0_count + RT_count + R1_count + R2_count - 1);
								ST_cb_cut_count++;
								IC_49_cut_count++;

								//if (sumx - xF2_best[k][setS3[j]] > R0_count + RT_count + R1_count + R2_count - 1)
								//{
								//	split_ind = 1;
								//}
							}
						}
					}
				}
#endif
	
				cbcut0.end();
				cbcut1.end();
				cbcut2.end();
				cbcutT.end();

				vector<int>().swap(setS1);
				//vector<int>().swap(setS2);
				vector<int>().swap(setS3);


			}
		}
	}


#ifdef split_cut2

	start_SP = clock();

	/////****加上split速度慢***/////////
	for (int k = 1; k < nbCrane; k++)
		if (QC_ind0[k - 1] == 1)
	{
		SP_count_split++;

		IloInt obj;
		obj = Split_SP_Callback(1, subcplex, nbQ, nbLocation_BD, nbprecR, yF_best, k - 1);

		if (obj == 0)
		{
			IloExpr BDcut(masterEnv);

			for (int k2 = k - 1; k2 <= k; k2++)
			{
				for (int i2 = 0; i2 < nbTask; i2++)
				{
					if (yF_best[k2][i2] > 0.9)
						BDcut += 1 - yF[k2][i2];
				}
			}
			add(BDcut >= 1);
			BDcut.end();
			//split_ind = 1;

			opt_cut_count_split++;
			cut_count++;
		}
		else if (obj > UB)
		{

			IloExpr BDcut(masterEnv);

			for (int k2 = k - 1; k2 <= k; k2++)
			{
				for (int i2 = 0; i2 < nbTask; i2++)
				{
					if (yF_best[k2][i2] > 0.9)
						BDcut += (obj - LB0)*(1 - yF[k2][i2]);
				}
			}
			add(CF - obj + BDcut >= 0);
			BDcut.end();

			opt_cut_count_split++;

			cut_count++;



		}
	}

	finish_SP = clock();
	duration_SP += (double)(finish_SP - start_SP) / CLOCKS_PER_SEC;
#endif

		}
	





	////cout << "2--" << endl;

}

ILOINCUMBENTCALLBACK7(IncumbentCallback, BoolVarMatrix, yF, NumVarMatrix, xF1, NumVarMatrix, xF2, IloNumArray, nbQ, IloIntArray, nbLocation, BoolMatrix, nbprecR,
	IloCplex, subcplex)
{
	IloInt i;
	IloEnv masterEnv = getEnv();
	IloInt numNodes = yF.getSize();
	IloInt numNodes_x2 = xF1.getSize();
	IloInt numNodes_y2 = xF2.getSize();

	int split_ind = 0;

	IloIntArray  QC_ind(masterEnv, nbCrane);


	IloNumArray2 yF_best(masterEnv, numNodes);
	for (i = 0; i < numNodes; ++i) {
		//cout << "xF[]: "<<xF[i].getSize() << endl;
		yF_best[i] = IloNumArray(masterEnv, yF[i].getSize());
		getValues(yF_best[i], yF[i]);
	}

		for (int k = 0; k < nbCrane; k++)
			for (i = 0; i < nbTask; i++)
				if (yF_best[k][i]>0.98)
					yF_best[k][i] = 1;

		clock_t start_SP = 0, finish_SP = 0;
		start_SP = clock();

#ifdef split_cut2
		/////****加上split速度慢***/////////
		for (int k = 1; k < nbCrane; k++)
			{
			IloInt obj;
			obj = Split_SP_Callback(1, subcplex, nbQ, nbLocation_BD, nbprecR, yF_best, k - 1);

			if ( obj == 0)
			{
				reject();
				break;
			}
			}
#endif

		{
			IloInt obj;
			obj = Split_SP_Callback(0, subcplex, nbQ, nbLocation_BD, nbprecR, yF_best, 0);
			if (obj == 0)
			{
				inf_cut_count++;
				reject();
			}
			//else if (obj > UB)
			//{
			//	opt_cut_count++;
			//	reject();
			//}

		}

		finish_SP = clock();
		duration_SP += (double)(finish_SP - start_SP) / CLOCKS_PER_SEC;





	////cout << "2--" << endl;

}


int _tmain(int argc, char* argv[], char* envp[])
{
	double GapAve[instance_no - start_instance+1];
	double GapAve2[instance_no - start_instance + 1];
	double IterAve[instance_no - start_instance + 1];

	double ObjAve[instance_no - start_instance + 1];
	int lowerBound[instance_no - start_instance + 1];
	int upperBound[instance_no - start_instance + 1];
	int cut37[instance_no - start_instance + 1];
	int cut38[instance_no - start_instance + 1];
	int cut43[instance_no - start_instance + 1];
	int cut44[instance_no - start_instance + 1];
	int cut45[instance_no - start_instance + 1];
	int cut48[instance_no - start_instance + 1];
	int cut49[instance_no - start_instance + 1];
	int nodesNO[instance_no - start_instance + 1];
	int spNO[instance_no - start_instance + 1];
	double DuraAve_MP[instance_no - start_instance + 1];
	double DuraAve_SP[instance_no - start_instance + 1];
	double DuraAve[instance_no - start_instance + 1];

	double DuraAve_SF[instance_no - start_instance + 1];
	double DuraAve_LB[instance_no - start_instance + 1];


	for (int my = start_instance; my <= instance_no; my++)
	{
		cout << "data-" << my << "数据运行中..." << endl;
		char* filename;
		char dream[100] = "test/data";
		filename = dream;
		char C1[3];
		char C2[3];
		char C3[3];
		char C4[3];
		itoa(my, C1, 10);

		//此处可编辑规模，以输入
		itoa(nbBay, C2, 10);
		itoa(nbCrane, C3, 10);
		itoa(nbTask, C4, 10);

		strcpy(filename, "test/");
		//此处可编辑规模，以输入
		strcat(filename, C4);
		strcat(filename, "-");
		strcat(filename, C2);
		strcat(filename, "-");
		strcat(filename, C3);

		strcat(filename, "/data");
		strcat(filename, "-");
		strcat(filename, C1);
		strcat(filename, ".txt");
		clock_t start = 0, finish = 0;
		clock_t start_MP = 0, finish_MP = 0;
		clock_t start_CB = 0, finish_CB = 0;
		clock_t start_LB = 0, finish_LB = 0;
		//clock_t start_SP = 0, finish_SP = 0;
		//start = clock();
		double  duration=0;
		double  duration_MP = 0;


		if (argc <= 1)
		{
			cerr << "Usage: " << argv[0] << " <model>" << endl;
			cerr << "  model = 0 -> convex piecewise linear model, " << endl;
			cerr << "  model = 1 -> concave piecewise linear model. [default]" << endl;
		}

		IloBool convex;
		if (argc <= 1)
			convex = IloFalse;
		else
			convex = atoi(argv[1]) == 0 ? IloTrue : IloFalse;

		//*************************//
		//     打开文件data.txt    //
		//*************************//
		ifstream fin(filename);
		IloEnv env;
		IloEnv subEnv;
		try 
		{

			//	定义原问题模型
			IloModel model(env);

			IloNum gap;

			//	定义临时变量
			IloInt i,k,j,t,kk; 

			//	定义连续决策变量CF
			IloIntVar   CF(env,0,1800);

			//IloNumArray tem_C(env, nbCrane);// 子问题计算的completion time
			IloNumVarArray QC_CF(env, nbCrane,0, IloInfinity);// each QC's completion time
			IloNumVarArray Task_CF(env, nbTask, 0, IloInfinity);// each task's completion time

			//	定义决策变量xF,CxF,CzF
			BoolVarMatrix xF(env, nbCrane);
			for(k = 0; k < nbCrane; k++)
			{				
				xF[k]=IloBoolVarArray(env,nbTask);
			}

			BoolVarMatrix yF(env, nbCrane);
			for (k = 0; k < nbCrane; k++)
			{
				yF[k] = IloBoolVarArray(env, nbTask);
			}


			IloBoolVarArray vF(env, nbTask);
			IloIntVarArray thetaF(env, nbCrane, 0, nbBay);

			BoolVarMatrix zF(env, nbCrane);
			for (i = 0; i < nbCrane; i++)
			{
				zF[i] = IloBoolVarArray(env, nbBay);
			}

			BoolVarMatrix xF_end(env, nbCrane);
			for (i = 0; i < nbCrane; i++)
			{
				xF_end[i] = IloBoolVarArray(env, nbTask);
			}

			NumVarMatrix mF(env, nbCrane);
			for (k = 0; k < nbCrane; k++)
				mF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);
			//**********************************//
			//            新变量           //
			//**********************************//
			//	定义保存决策变量xF,CxF,CzF的最优值
			BoolMatrix xF_best(env, nbCrane);
			for(k = 0; k < nbCrane; k++)
			{				
				xF_best[k]=IloBoolArray(env,nbTask);				
			}
			BoolMatrix yF_best(env, nbCrane);
			for (k = 0; k < nbCrane; k++)
			{
				yF_best[k] = IloBoolArray(env, nbTask);
			}

			BoolMatrix CxF_best(env, nbCrane);
			for (k = 0; k < nbCrane; k++)
			{
				CxF_best[k] = IloBoolArray(env, nbTask);
			}
			IloBoolArray CyF_best(env, nbTask);


			BoolMatrix ubF_best(env, nbCrane);
			for (k = 0; k < nbCrane; k++)
			{
				ubF_best[k] = IloBoolArray(env, nbTask);
			}			
			BoolMatrix vbF_best(env, nbCrane);
			for (k = 0; k < nbCrane; k++)
			{
				vbF_best[k] = IloBoolArray(env, nbTask);
			}

			BoolMatrix zF_best(env, nbCrane);
			for (k = 0; k < nbCrane; k++)
			{
				zF_best[k] = IloBoolArray(env, nbBay);
			}


			IloInt   CF_best;
			IloNumArray QC_CF_best(env, nbCrane);
			IloNumArray Task_CF_best(env, nbTask);

			IloIntArray left_bay(env, nbCrane);
			IloIntArray right_bay(env, nbCrane);

			///////////测试用
			BoolMatrix2 FxF_best(env);
			for (k = 0; k < nbCrane; k++)
			{
				FxF_best.add(BoolMatrix(env));
				for (i = 0; i <nbTask; i++)
					FxF_best[k].add(IloBoolArray(env, nbTask));
			}
			BoolMatrix FzF_best(env, nbTask);
			for (i = 0; i <nbTask; i++)
			{
				FzF_best[i] = IloBoolArray(env, nbTask);
			}

			BoolMatrix xF_begin_best(env, nbCrane);
			for (k = 0; k < nbCrane; k++)
			{
				xF_begin_best[k] = IloBoolArray(env, nbTask);
			}

			BoolMatrix xF_end_best(env, nbCrane);
			for (k = 0; k < nbCrane; k++)
			{
				xF_end_best[k] = IloBoolArray(env, nbTask);
			}


			NumMatrix wF_best(env, nbCrane);
			for (k = 0; k < nbCrane; k++)
				wF_best[k] = IloNumArray(env, nbBay);

			NumMatrix CwF_best(env, nbCrane);
			for (k = 0; k < nbCrane; k++)
				CwF_best[k] = IloNumArray(env, nbBay);

			//*******************************************************//
			//                 定义主问题输出的参数                         //
			//*******************************************************// 	 

			//约束中参数
			//	定义需代入子问题计算的work zone数目
			IloInt   nbn;
			//定义 work zone 起始 section
			IloIntArray nbn_start(env, nbBay);
			//	定义 work zone 终止 section
			IloIntArray  nbn_end(env,nbBay);
			//	定义每个work zone的crane个数
			IloIntArray  nbn_crane_no(env,nbBay);
			//	定义每个work zone的最后一个crane的编号
			IloIntArray  nbn_crane_index(env,nbBay);
			//	定义每个work zone的第一个crane的编号
			IloIntArray  nbn_crane_start(env,nbBay);
			


			IloInt cut_1_no=0;
			IloInt cut_2_no = 0;


			IC_48_cut_count = 0;
			IC_49_cut_count = 0;
			inf_cut_count = 0;
			opt_cut_count = 0;
			SP_count = 0;
			BC_node_count = 0;

			//*******************************************************//
			//                 定义输入参数                          //
			//*******************************************************// 	 
			//约束中参数
			//	定义能力常量参数s
			IloNum   nbs;
			//定义吊机位置初始状态参数
			IloBoolArray nbb(env, nbCrane);
			//	定义任务量参数nbQ
			IloNumArray  nbQ(env, nbTask);

			IloNum  aveQ=0;

			IloIntArray  nbreadyT(env, nbCrane);

			//	定义任务所在贝位参数nbQ
			IloIntArray  nbLocation(env, nbTask);

			BoolMatrix nbprecR(env,nbTask);//优先级关系
			for (i = 0; i < nbTask; i++)  nbprecR[i] = IloBoolArray(env, nbTask);
			for (int ai = 0; ai < nbTask; ai++)
				for (int bi = 0; bi < nbTask; bi++)
					nbprecR[ai][bi] = 0;

			//*******************************************************//
			//                 定义callback参数                      //
			//*******************************************************// 	 


	 
			//*******************************************************//
			//                 读入参数数据                          //
			//*******************************************************// 	 

			IloIntArray ls1(env, 8);

			fin >> ls1;
			//cout << ls1[2] << "  ";
			fin >> nbQ;
			fin >> nbLocation;
			fin >> nbreadyT;
			fin >> nbb;	

			for (i = 0; i < nbTask; i++)
				nbLocation_BD[i] = nbLocation[i];

			for (i = 0; i < nbCrane; i++)
				nbb_BD[i] = nbb[i];


			////读入数据nbLocation[i]
			//for (i = 0; i < nbTask; i++)
			//{
			//	//fin >> nbLocation[i];
			//	cout<<nbLocation[i]<<"  ";
			//}
			////cout<<endl;

			////	读入数据nbb[i]
			//for (i = 0; i < nbCrane; i++)
			//{
			//	//fin >> nbb[i];
			//	cout<<nbb[i]<<"  ";
			//}
			//cout<<endl;

			IloIntArray prec(env, 2);

			//	读入数据nbprecR
			for (i = 0; i < ls1[2]; i++)
			{
				fin >> prec;

				for (int ai = 0; ai < nbTask; ai++)
					if (ai == prec[0] - 1)
						for (int bi = 0; bi < nbTask; bi++)
							if (bi == prec[1] - 1)
								nbprecR[ai][bi] = 1;


			}


			//for (int ai = 0; ai < nbTask; ai++)
			//	for (int bi = 0; bi < nbTask; bi++)
			//		if (nbprecR[ai][bi] == 1)
			//		{
			//	cout << "(" << ai + 1 << "," << bi + 1 << "), ";
			//		}cout << endl;


			nbs = 1;

			NumMatrix travel_time(env, nbTask);
			for (i = 0; i < nbTask; i++)
			{
				travel_time[i] = IloNumArray(env, nbTask);
			}
			for (i = 0; i < nbTask; i++)
			{
				for (j = 0; j < nbTask; j++)
				{

					if (nbLocation[i] == nbLocation[j]) travel_time[i][j] = 0;
					else if (nbLocation[i] > nbLocation[j])
					{
						travel_time[i][j] = (nbLocation[i] - nbLocation[j])*crane_move_time;
					}
					else
					{
						travel_time[i][j] = (nbLocation[j] - nbLocation[i])*crane_move_time;
					}

				}

			}

			for (i = 0; i < nbBay; i++)
			{
				cout << i + 1 << ":  ";
				for (j = 0; j < nbTask; j++)
					if (nbLocation[j] == i + 1)
					{
						cout << j + 1 << "  ";
					}
				cout << endl;
			}cout << "    "; cout << endl;

			int prec_cut_no = 0;

			IloIntArray botom_task_index(env,nbBay);
			for (j = 0; j < nbTask; j++)
				botom_task_index[nbLocation[j]-1] = j;

			////为加cut定义的辅助数组，避免重复加cut
			int botom_violate_task[nbCrane][nbBay];
			for (i = 0; i < nbBay; i++)
				for (k = 0; k < nbCrane; k++)
					botom_violate_task[k][i] = 1;


			//linear relaxation
			NumVarMatrix xF2(env, nbCrane);
			for (k = 0; k < nbCrane; k++)
				xF2[k] = IloNumVarArray(env, nbTask, 0, 1);

			NumVarMatrix xF1(env, nbCrane);
			for (k = 0; k < nbCrane; k++)
				xF1[k] = IloNumVarArray(env, nbTask, 0, 1);


			//************************************************************************************//
			//      至此已知数据输入完毕                                                          //
			//************************************************************************************//
			char* filename1;
			char dream1[100] = "result/C_data";
			filename1 = dream1;

			itoa(my, C1, 10);

			//此处可编辑规模，以输入
			itoa(nbBay, C2, 10);
			itoa(nbCrane, C3, 10);
			itoa(nbTask, C4, 10);

			strcpy(filename1, "result/");
			//此处可编辑规模，以输入
			strcat(filename1, C4);
			strcat(filename1, "-");
			strcat(filename1, C2);
			strcat(filename1, "-");
			strcat(filename1, C3);

			strcat(filename1, "/Cdata");
			strcat(filename1, "-");
			strcat(filename1, C1);
			strcat(filename1, ".txt");

			ofstream fout(filename1);

			cout<<filename1<<endl;


			for (int b = 0; b < nbBay; b++)
			{
				cout << "bay " << b << ": ";
				fout << "bay " << b << ": ";
				for (i = 0; i < nbTask; i++)
				{
					if (nbLocation[i] == b + 1)
					{
						cout << i << "  ";
						fout << i << "  ";
					}

				}
				cout << endl;
				fout << endl;
			}

			//*******************************************************//
			//                 建立主问题 CBMP                       //
			//*******************************************************// 
			int sumiter = 0;//统计迭代数目
			IloNum UB0;

			UB = 2000;

			start = clock();
			start_LB = clock();

			IloNum LB0;
			//LB_Same_Direction(model, nbs, nbb, nbQ, nbreadyT,
			//	yF, CF, &CF_best, yF_best,
			//	&ObjVal, &LB, nbLocation, nbprecR, my);
			//LB0 = LB;

			LB_Chen_Same_Direction(model, nbs, nbb, nbQ, nbreadyT,
				yF, CF, &CF_best, yF_best,
				&ObjVal, &LB, nbLocation, nbprecR, my);
			LB0 = LB;

			finish_LB = clock();
			start_CB = clock();


			//Same_Direction(model, nbs, nbb, nbQ, nbreadyT,
			//	yF, CF, &CF_best, yF_best,
			//	&ObjVal, &UB0, nbLocation, nbprecR, my);
			Chen_Same_Direction(model, nbs, nbb, nbQ, nbreadyT,
				yF, CF, &CF_best, yF_best,
				&ObjVal, &UB0, nbLocation, nbprecR, my);


			IloBoolArray nbb2(env, nbCrane);
			IloNumArray  nbQ2(env, nbTask);
			//IloNumArray  nbs2(env, nbCrane);
			IloIntArray  nbLocation2(env, nbTask);
			IloIntArray  nbreadyT2(env, nbCrane);

			//re-numbering the bays
			for (j = 0; j < nbTask; j++)
			{
				nbLocation2[j] = nbBay + 1 - nbLocation[j];
			}
			for (k = 0; k < nbCrane; k++)
				nbb2[nbCrane - 1 - k] = nbBay + 1 - nbb[k];
			//for (k = 0; k < nbCrane; k++)
			//	nbs2[nbCrane - 1 - k] = nbs[k];
			for (k = 0; k < nbCrane; k++)
				nbreadyT2[nbCrane - 1 - k] = nbreadyT[k];

			//Same_Direction(model, nbs, nbb2, nbQ, nbreadyT2,
			//	yF, CF, &CF_best, xF_best,
			//	&ObjVal, &UB, nbLocation2, nbprecR, my);
			Chen_Same_Direction(model, nbs, nbb2, nbQ, nbreadyT2,
				yF, CF, &CF_best, xF_best,
				&ObjVal, &UB, nbLocation2, nbprecR, my);


			if (UB<UB0)
			{
				UB0 = UB;

				for (k = 0; k < nbCrane; k++)
					for (j = 0; j < nbTask; j++)
						yF_best[k][j] = xF_best[k][j];
			}

			finish_CB = clock();

			UB = UB0;

			fout << "\LB0 = " << LB0 << ", time=" << (double)(finish_LB - start_LB) / CLOCKS_PER_SEC << endl << endl;;
			fout << "\Unidirectional Makespan = " << UB0 << ", time="<<(double)(finish_CB - start_CB) / CLOCKS_PER_SEC<<endl << endl;

			int cut_ind_0 = 0;
			 
			if (LB - UB<1 && UB - LB<1)//解至最优最优
			{
				model.end();
				cut_ind_0 = 1;

				finish = clock();

				sumiter++;
				goto out_end;
			}

			IloNumVarArray betaF(env, nbCrane, 0, 100);

			BoolVarMatrix CzF(env, nbCrane);
			for (k = 0; k < nbCrane; k++)
			{
				CzF[k] = IloBoolVarArray(env, nbBay);
			}

			BoolVarMatrix CzF2(env, nbCrane);
			for (k = 0; k < nbCrane; k++)
			{
				CzF2[k] = IloBoolVarArray(env, nbBay);
			}

			MP_Same_Direction(model, nbs, nbb, nbQ, 
				yF, CF, &CF_best, yF_best, mF, CzF2, betaF,
				&ObjVal, &UB, nbLocation, nbprecR, xF1, xF2);

			IloCplex cplex(env);
			cplex.extract(model);





#ifdef UserActive
			//cplex.use(CandyUserCallback(env));
			//cplex.setParam(IloCplex::MIPEmphasis, 2);
			//cplex.setParam(IloCplex::HeurFreq, 1);
			//cplex.setParam(IloCplex::ParallelMode, 1);
			//cplex.setParam(IloCplex::Threads, 4);
#endif



			//cplex.setParam(IloCplex::MIPEmphasis, CPX_MIPEMPHASIS_FEASIBILITY);
			//cplex.setParam(IloCplex::EpGap, 0.001);
			//cplex.setOut(env.getNullStream());
			cplex.setParam(IloCplex::TiLim, 1800);
			//cplex.setParam(IloCplex::Threads, 4);

			//cplex.setParam(IloCplex::VarSel, 3);
			//cplex.setParam(IloCplex::MIPEmphasis, 3);
			//cplex.setParam(IloCplex::Probe, 3);
			
	//cout<<endl<<"number of binary variables = "<<cplex.getNbinVars()<<endl;
	//cout<<endl<<"number of constraints = "<<cplex.getNrows()<<endl;



			//**********************************//
			//             Benders 主算法           //
			//**********************************//	

			gap = 100;
			int countNO = 0;

			//int cbcut_count = 0;

			int update_no = 0;// if bidirectional solution is updated
			int sp_no;//which sp is the optimal, 0 or 1



			//**********************************//
			//             Benders 子问题           //
			//**********************************//	

			IloCplex subcplex(subEnv);
			//	定义决策变量xF,CxF,CzF
			//BoolVarMatrix xF1_sub(subEnv, nbCrane);
			//for (k = 0; k < nbCrane; k++)
			//{
			//	xF1_sub[k] = IloBoolVarArray(subEnv, nbTask);
			//}
			//BoolVarMatrix xF2_sub(subEnv, nbCrane);
			//for (k = 0; k < nbCrane; k++)
			//{
			//	xF2_sub[k] = IloBoolVarArray(subEnv, nbTask);
			//}





			//cplex.use(BendersUserCallback(env, yF, xF1, xF2, nbQ, nbLocation, nbprecR, subcplex));
			cplex.use(BendersLazyCallback(env, yF, xF1, xF2, nbQ, CF, nbprecR, subcplex));
			//cplex.use(IncumbentCallback(env, yF, xF1, xF2, nbQ, nbLocation, nbprecR, subcplex));

			
				//**********************************//
				//            开始求解		        //
				//**********************************//

			//for (k = 0; k < nbCrane; k++)
			//	cplex.setObject()

			IloNumVarArray startVar(env);
			IloNumArray startVal(env);
			for (int i = 0; i < nbCrane; ++i)
				for (int j = 0; j < nbTask; ++j) {
				startVar.add(yF[i][j]);
				startVal.add(yF_best[i][j]);
				}
			cplex.addMIPStart(startVar, startVal);
			startVal.end();
			startVar.end();


				BOOL h1;

				start_MP = clock();

				h1 = cplex.solve();

				finish_MP = clock();
				finish = clock();


				//cout<<"h0 "<<h0<<endl;
				if (!h1)
				{
					cout << "\nno feasible solution has been found，algorithm terminate" << endl;
					fout << "\nno feasible solution has been found，algorithm terminate" << endl;
					cplex.clearModel();
					cplex.clear();
					cplex.end();
					model.end();

					//sumiter++;
					goto out_end;//此LB下没有可行解，因此算法结束

					//break;

					//return FALSE;
				}
				//**********************************//
				//             记录CBMP最好解       //
				//**********************************//
				//ObjVal = cplex.getValue(CF);
				//LB=cplex.getValue(CF);
				
				//if (sumiter == 0) 
					LB0 = LB;

				LB = cplex.getBestObjValue();

				BC_node_count =cplex.getNnodes();

				

				if (finish_MP - start_MP > 1000)
					GapAve2[my - 1] = cplex.getMIPRelativeGap();
				else
					GapAve2[my - 1] = 0;


				////finish = clock();

				//int l_bay_k[nbCrane];
				////cout<<"xF_best"<<endl;
				//for (k = 0; k < nbCrane; k++)
				//{
				//	cout << k << ":  ";
				//	l_bay_k[k] = nbBay;

				//	for (i = 0; i < nbTask; i++)
				//	{

				//		if (cplex.isExtracted(yF[k][i]))
				//		{
				//			if (cplex.getValue(yF[k][i]) > 0.1)
				//			{
				//				yF_best[k][i] = 1;
				//				if (nbLocation[i] - 1 < l_bay_k[k])
				//					l_bay_k[k] = nbLocation[i] - 1;

				//			}

				//			else
				//				yF_best[k][i] = 0;
				//		}
				//		else
				//			yF_best[k][i] = 0;

				//		if (yF_best[k][i] == 1) cout << i << "  ";

				//	}
				//	cout << "    "; cout << endl;
				//}
				//cout<<endl<<endl;

				////NumMatrix xF1_best(env, nbCrane);
				////for (k = 0; k < nbCrane; k++)
				////{
				////	xF1_best[k] = IloNumArray(env, nbTask);
				////}
				////NumMatrix xF2_best(env, nbCrane);
				////for (k = 0; k < nbCrane; k++)
				////{
				////	xF2_best[k] = IloNumArray(env, nbTask);
				////}

				////for (k = 0; k < nbCrane; k++)
				////{
				////	cout << k << ":  ";
				////	l_bay_k[k] = nbBay;

				////	for (i = 0; i < nbTask; i++)
				////	{

				////		if (cplex.isExtracted(xF1[k][i]))
				////		{
				////			if (cplex.getValue(xF1[k][i]) > 0.01)
				////			{
				////				xF1_best[k][i] = cplex.getValue(xF1[k][i]);
				////				cout << i << "-f-" << xF1_best[k][i] << "  ";
				////			}
				////			else
				////				xF1_best[k][i] = 0;
				////		}
				////		else
				////			xF1_best[k][i] = 0;

				////		if (cplex.isExtracted(xF2[k][i]))
				////		{
				////			if (cplex.getValue(xF2[k][i]) > 0.01)
				////			{ 
				////				cout << i << "-d-" << cplex.getValue(xF2[k][i]) << "  ";
				////				xF2_best[k][i] = cplex.getValue(xF2[k][i]);
				////			}
				////			else
				////				xF2_best[k][i] = 0;
				////				
				////		}
				////		else
				////			xF2_best[k][i] = 0;
				////	}
				////	 cout << endl;
				////}
				////cout<<endl<<endl;

				////	输出下界
				//fout << "LB" << sumiter + 1 << "=" << LB << "-time-" << (double)(finish - start) / CLOCKS_PER_SEC << "\t";
				////fout << "LB" << sumiter + 1 << "(makespan)=" << CF_iter_LB << "\t";
				//duration_MP = (double)(finish_MP - start_MP) / CLOCKS_PER_SEC;

				////if ((double)(finish - start) / CLOCKS_PER_SEC > 1800)
				////	break;

				////**********************************//
				////             求解子问题		       //
				////**********************************//

				//BOOL h0;

				////子问题
				//h0 = SP_1(model, nbs, nbb, nbQ, nbLocation, nbprecR, nbreadyT,
				//	xF, yF, zF, CF, QC_CF, thetaF, vF, yF_best, CxF_best, UB, &ObjVal, wF_best, CwF_best);

				//if (! h0)
				//{
				//	cout << "h1= " << h0 << endl;
				//}
				//

				//if (h0 == 1)
				//{
				//	/////更新上界
				//	if (ObjVal < UB)
				//	{
				//		UB = ObjVal;
				//		//cplex.addCut(CF<=UB-1);//更新主问题模

				//		update_no++;
				//		sp_no = 1;

				//	}

				//	//IloExpr BDcut(env);

				//	//for (k = 0; k < nbCrane; k++)
				//	//{
				//	//	for (i = 0; i < nbTask; i++)
				//	//	{
				//	//		if (yF_best[k][i] == 1)
				//	//			BDcut += (ObjVal - LB0)*(1 - yF[k][i]);
				//	//	}
				//	//}
				//	////cplex.addCut(BDcut >= 1);
				//	//cplex.addCut(BDcut + CF >= ObjVal);
				//	//BDcut.end();

				//}
				////else
				////{
				////	//IloExpr BDcut(env);

				////	//for (k = 0; k < nbCrane; k++)
				////	//{
				////	//	for (i = 0; i < nbTask; i++)
				////	//	{
				////	//		if (yF_best[k][i] == 1)
				////	//			BDcut += 1 - yF[k][i];
				////	//	}
				////	//}
				////	//cplex.addCut(BDcut >= 1);
				////	//BDcut.end();
				////}
				//	if (h0)
				//		fout << "UB_current" << sumiter + 1 << "=" << ObjVal << "-time-" << (double)(finish - start) / CLOCKS_PER_SEC << "      ";
				//	else
				//		fout << "UB_current" << sumiter + 1 << " inf" << "-time-" << (double)(finish - start) / CLOCKS_PER_SEC << "      ";				
				//	fout << "UB" << sumiter + 1 << "=" << UB << "-time-" << (double)(finish - start) / CLOCKS_PER_SEC << endl;
				//	//duration_SP += (double)(finish_SP - start_SP) / CLOCKS_PER_SEC;

				//	gap = (UB - LB) / UB;
				//	sumiter++;

			
				cplex.clearModel();
				cplex.clear();
				cplex.end();



			//**********************************//
			//             输出最好解           //
			//**********************************//	


				//**********************************//
				//             计算下界           //
				//**********************************//	
out_end:		cout << "update_no: " << update_no << endl;
				



			//**********************************//
			//             计算目标值           //
			//**********************************//	

			//double a,b; 
			//a=0,b=0;
			////	计算OBJ1
			////for(i = 0; i < nbMaterial; i++)for(j = 0; j < nbBay; j++) for(k = 0; k < nbBay; k++) a+=nbcF[i][j][k]*xF_best[i][j][k];			
			////for(j = 0; j < nbBay; j++) for(k = 0; k < nbBay; k++) a+=nNs*wF_best[k][k];


			////计算OBJ2
			//b=CF_best;

			//GapAve[my-1]=gap;


			//fout.precision(10);
			////fout<<"OBJ1 = "<<a<<endl;
			//fout<<"OBJ2 = "<<b<<endl;
			//fout.precision(10);
			//fout<< "\nobj1+obj2 = " << ObjVal << endl;





			fout.precision(10);
			fout << endl;
			//fout << "\detour = " << LB << endl;
			fout << "LB0 = " << LB0 << endl;
			fout << "Unidirectional Makespan = " << UB0 << endl;
			fout << "LB = " << LB << endl;
			fout<< "UB = " << UB << endl;
			//if (sumiter>0) fout << "real makespan = " << CF_best << endl;
			//fout << "\LB0 = " << LB0 << endl; 
			
			//fout<< "\H_UB = " << H_UB << endl;

			//fout<< "\niteration = " << sumiter << endl;
			//fout << "\Logic_Cut = " << countNO << endl;
			//fout << "comb_Cut = " << cb_cut_count << endl;
			//fout << "Strong comb_Cut = " << ST_cb_cut_count<< endl;
			//fout << "prec_Cut= " << prec_cut_count << endl;
			//fout << "bd_Cut= " <<cut_count << endl;
			//fout << "user_Cut= " << user_cut_count << endl;


			fout << "IC_48_cut_count = " << IC_48_cut_count << endl;
			fout << "IC_49_cut_count = " << IC_49_cut_count << endl;
			fout << "inf_cut_count= " << inf_cut_count << endl;
			fout << "opt_cut_count= " << opt_cut_count << endl;
			fout << "opt_cut_count_split= " << opt_cut_count_split << endl;
			fout << "SP_Count= " << SP_count << endl;
			fout << "SP_Count_split= " << SP_count_split << endl;
			fout << "BC_node_count= " << BC_node_count << endl;



			ObjAve[my - start_instance] = UB;
			lowerBound[my - start_instance] = LB0;
			upperBound[my - start_instance] = UB0;
			if (cut_ind_0)
			{
				cut37[my - start_instance] = 0;
				cut38[my - start_instance] = 0;
				cut43[my - start_instance] = 0;
				cut44[my - start_instance] = 0;
				cut45[my - start_instance] = 0;
				cut48[my - start_instance] = 0;
				cut49[my - start_instance] = 0;
				nodesNO[my - start_instance] = 0;
				spNO[my - start_instance] = 0;

			} 
			else
			{
				cut37[my - start_instance] = inf_cut_count;
				cut38[my - start_instance] = opt_cut_count;
				cut43[my - start_instance] = ls1[2] * (nbCrane - 1);
				cut44[my - start_instance] = ls1[2] * nbCrane;
				cut45[my - start_instance] = nbTask*nbCrane;
				cut48[my - start_instance] = IC_48_cut_count;
				cut49[my - start_instance] = IC_49_cut_count;
				nodesNO[my - start_instance] = BC_node_count;
				spNO[my - start_instance] = SP_count;

			}

			double duaCB, duaLB;
			duration = (double)(finish - start) / CLOCKS_PER_SEC;
			duration_MP = (double)(finish_MP - start_MP) / CLOCKS_PER_SEC;
			duaCB=(double)(finish_CB - start_CB) / CLOCKS_PER_SEC;
			duaLB = (double)(finish_LB - start_LB) / CLOCKS_PER_SEC;

			//fout<< "\nGAP = " << GapF << endl;

			fout << "MP_Time=" << duration_MP - duration_SP << endl;
			fout << "SP_Time=" << duration_SP << endl;

			fout << "CB_Time=" << duaCB << endl;
			fout << "LB_Time=" << duaLB << endl;

			fout << "MP_gap=" << GapAve2[my - start_instance] << endl;

			

			fout<<"Time="<<duration<<endl;
			cout<<"Time="<<duration<<endl;

			DuraAve[my- start_instance] = duration;
			DuraAve_MP[my - start_instance] = duration_MP-duration_SP;
			DuraAve_SP[my - start_instance] = duration_SP;
			ObjAve[my - start_instance] = UB;
			IterAve[my - start_instance] = sumiter;

			


			DuraAve_SF[my - start_instance] = duaCB;
			DuraAve_LB[my - start_instance] = duaLB;



			duration_SP = 0;
			cb_cut_count = 0;
			cut_count = 0;
			user_cut_count = 0;
			prec_cut_count = 0;
			ST_cb_cut_count = 0;
			SP_count_split = 0;
			opt_cut_count_split = 0;
			//cout<<SubCutNum<<endl;
			//cout<<"lala"<<endl;

			model.end();

		}
		catch (IloException& e) 
		{
			cerr << "ERROR: " << e.getMessage() << endl;
		}
		catch (...) 
		{
			cerr << "Error" << endl;
		}
TERMINATE : env.end();
	}
	double gapa=0;
	double gapa2 = 0;
	double durationa=0;
	double durationa_mp = 0;
	double durationa_sp = 0;
	double obja = 0;
	double itera = 0;
	for (int i=0;i<instance_no-start_instance+1;i++)
	{
		gapa+=GapAve[i];
		gapa2 += GapAve2[i];
		durationa += DuraAve[i];
		durationa_mp += DuraAve_MP[i];
		durationa_sp += DuraAve_SP[i];
		obja += ObjAve[i];
		itera += IterAve[i];
	}


	durationa = durationa / (instance_no - start_instance + 1);
	durationa_mp = durationa_mp / (instance_no - start_instance + 1);
	durationa_sp = durationa_sp / (instance_no - start_instance + 1);
	gapa = gapa / (instance_no - start_instance + 1);
	gapa2 = gapa2 / (instance_no - start_instance + 1);
	obja = obja / (instance_no - start_instance + 1);
	itera = itera / (instance_no - start_instance + 1);
	cout<<"average gap: "<<gapa<<endl;
	cout << "average gap: " << gapa2 << endl;
	cout<<"average CPU: "<<durationa<<endl;
	cout << "average CPU of MP: " << durationa_mp << endl;
	cout << "average CPU of SP: " << durationa_sp << endl;
	cout << "average obj: " << obja << endl;
	cout << "average iter: " << itera << endl;

	ofstream oFile, oFile1;
	string thefilename1 = "haha.csv";

	for (int i = start_instance; i <= instance_no; i++) {
		oFile1.open(thefilename1, ios::app);//| ios::trunc // 这样就很容易的输出一个需要的excel 文件  

		oFile1 << i << "," << nbTask << "," << multiple_1*ObjAve[i - start_instance] << "," << lowerBound[i - start_instance] << "," << upperBound[i - start_instance] << ","
			<< cut37[i - start_instance] << "," << cut38[i - start_instance] << "," << cut43[i - start_instance] << "," << cut44[i - start_instance] << "," << cut45[i - start_instance] << "," 
			<< cut48[i - start_instance] << "," << cut49[i - start_instance] << ","
			<< nodesNO[i - start_instance] << "," << spNO[i - start_instance] << ","
			<< fixed << setprecision(2) << DuraAve_MP[i - start_instance] << "," 
			<< fixed << setprecision(2) << DuraAve_SP[i - start_instance] << ","
			<< fixed << setprecision(2) << DuraAve[i - start_instance] << "," << ","
			<< upperBound[i - start_instance] << "," << fixed << setprecision(2) << DuraAve_SF[i - start_instance] << "," << ","
			<< lowerBound[i - start_instance] << "," << fixed << setprecision(2) << DuraAve_LB[i - start_instance] << endl;

		oFile1.close();
	}


	system("pause");
	return 0;
}




BOOL MP_Same_Direction(IloModel model, IloNum nbs, IloIntArray  nbb, IloNumArray nbQ,
	BoolVarMatrix yF, IloIntVar CF, IloInt *CF_best, BoolMatrix yF_best, NumVarMatrix mF, BoolVarMatrix endCzF, IloNumVarArray betaF,
	IloNum *ObjVal, IloNum *UB, IloIntArray  nbLocation, BoolMatrix nbprecR, NumVarMatrix xF1, NumVarMatrix xF2)
{
	IloEnv env = model.getEnv();
	IloInt i, j, k;

	//问题模型
	//IloModel submodel(env);
	IloNumVarArray t0kF(env, nbCrane, 0, 100);
	
	//leaving time at each bay
	NumVarMatrix TF2(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		TF2[k] = IloNumVarArray(env, nbBay, -100, IloInfinity);

	BoolVarMatrix CKF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		CKF[k] = IloBoolVarArray(env, nbBay);

	BoolVarMatrix start_zF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		start_zF[k] = IloBoolVarArray(env, nbBay);
	}

	//NumVarMatrix DzF(env, nbCrane);
	//for (k = 0; k < nbCrane; k++)
	//{
	//	DzF[k] = IloNumVarArray(env, nbBay, 0,1);
	//}

	BoolVarMatrix DzF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		DzF[k] = IloBoolVarArray(env, nbBay);
	}

	//**********************************//
	//            松弛问题变量       //
	//**********************************//


	//问题模型
	//IloNumVarArray t0kF(env, nbCrane, 0, 100);
	IloNumVarArray gammaF(env, nbCrane, 0, 100);//

	NumVarMatrix wF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		wF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);

	NumVarMatrix CwF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		CwF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);

	NumVarMatrix CzF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		CzF[k] = IloNumVarArray(env, nbBay, 0, 1);
	}


	//BoolVarMatrix CzF(env, nbCrane);
	//for (k = 0; k < nbCrane; k++)
	//{
	//	CzF[k] = IloBoolVarArray(env, nbBay);
	//}


	//NumVarMatrix zF(env, nbCrane);
	//for (k = 0; k < nbCrane; k++)
	//{
	//	zF[k] = IloNumVarArray(env, nbBay, 0, 1);
	//}


	//BoolVarMatrix zF(env, nbCrane);
	//for (k = 0; k < nbCrane; k++)
	//{
	//	zF[k] = IloBoolVarArray(env, nbBay);
	//}


	//leaving time at each bay
	NumVarMatrix TF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		TF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);
	NumVarMatrix CTF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		CTF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);

	//IloBoolVarArray vF(env, nbCrane);
	IloNumVarArray vF(env, nbCrane,0,1);


	//**********************************//
	//            原问题目标函数        //
	//**********************************//

	//IloExpr obj1(env); 


	IloExpr  obj2(env);
	//  建立子问题目标函数表达式 

	obj2 += CF;

	//	将目标函数加入到原问题模型
	model.add(IloMinimize(env, obj2));//
	//obj1.end();
	obj2.end();


	//model.add(CF<= *UB);

	//**********************************//
	//            MP问题 约束           //
	//**********************************//


	//定义TF
	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbBay; i++)
		{
		IloExpr  epa(env);
		epa += t0kF[k];

		for (j = 0; j <= i; j++)
		{
			epa -= QCmove_time*j * start_zF[k][j];
		}

		for (j = 0; j < nbTask; j++)
			if (nbLocation[j] - 1 <= i)
				epa += (nbQ[j] / nbs)*yF[k][j];
		for (j = 0; j <= i; j++)
		{

			epa += mF[k][j];
		}
		model.add(epa + QCmove_time*i - TF2[k][i] == 0);
		epa.end();

		}

	IloRangeArray  c2(env);
	for (k = 0; k < nbCrane; k++)
		//for (i = 0; i < nbTask; i++)
	{
		IloExpr  epa(env);
		epa += CF - t0kF[k];
		//epa -= betaF[k];
		epa -= gammaF[k];
		for (int ii = 0; ii < nbTask; ii++)
			epa -= (nbQ[ii] / nbs)*yF[k][ii];
		//for (j = 0; j < nbLocation[i]-1; j++)	
		for (j = 0; j < nbBay; j++)
			epa -= mF[k][j];
			//epa -= wF[k][j] + CwF[k][j];
		//epa -= (nbLocation[i] - 1)*yF[k][i];

		//epa -= betaF[k];

		for (j = 0; j < nbBay; j++)
			epa -= QCmove_time*j*endCzF[k][j];

		//for (j = 0; j < nbLocation[i]; j++)
		for (j = 0; j < nbBay; j++)
			epa += QCmove_time*j*start_zF[k][j];

		c2.add(epa >= 0);
		epa.end();
	}
	model.add(c2);
	c2.end();


	////endCzF取值
	for (j = 0; j < nbTask; j++)
		if (nbLocation[j] > 1)
			for (k = 0; k < nbCrane; k++)
			{
		IloExpr  epa(env);
		for (i = 0; i < nbLocation[j] - 1; i++)
			epa += endCzF[k][i];
		epa += yF[k][j];

		model.add(epa <= 1);
		epa.end();
			}
	//////endCzF
	IloRangeArray  c00(env);
	for (k = 0; k < nbCrane; k++)
		for (j = 0; j < nbTask; j++)
		{
		IloExpr  epa(env);

		for (i = 0; i < nbBay; i++)
			epa += i*endCzF[k][i];
		//epa += betaF[k];
		epa -= (nbLocation[j] - 1)*yF[k][j];

		c00.add(epa >= 0);
		epa.end();
		}
	model.add(c00);
	c00.end();
	//for (j = 0; j <nbTask; j++)
	//	for (k = 0; k < nbCrane; k++)					
	//	{
	//	IloExpr  epa(env);

	//	 epa += yF[k][j];
	//	 for (i = 0; i < nbBay; i++)
	//		 if (nbLocation[j]<= i+1)
	//			 epa -= endCzF[k][i];

	//	model.add(epa <= 0);
	//	epa.end();
	//	}
	//////endCzF进一步限定
	//for (i = 0; i < nbBay; i++)
	//	for (k = 0; k < nbCrane; k++)
	//	{
	//	IloExpr  epa(env);

	//	epa += endCzF[k][i];
	//	for (j = 0; j < nbTask; j++)
	//		if (nbLocation[j] == i + 1)
	//			epa -= yF[k][j];

	//	model.add(epa <= 0);
	//	epa.end();
	//	}





	//约束（3）
	IloRangeArray  c3(env);
	for (i = 0; i < nbBay - 1 - safe_margin; i++)
		for (k = 0; k < nbCrane - 1; k++)
		{
		IloExpr  epa(env);

		epa += TF2[k][i];

		//for (j = 0; j < i; j++)
		//{
		//	epa -= 1000 * CzF[k][j];
		//}

		epa -= TF2[k + 1][i + 1 + safe_margin];

		for (j = 0; j <= i + 1 + safe_margin; j++)
		{
			epa -= 2000 * start_zF[k + 1][j];
		}
		for (j = 0; j < i; j++)
		{
			epa += 2000 * endCzF[k][j];
		}

		epa += 2000;

		c3.add(epa >= 0);
		epa.end();
		}
	model.add(c3);
	c3.end();




	//约束（4）
	IloRangeArray  c4(env);
	for (i = 0; i < nbTask; i++)
	{
		IloExpr  epa(env);
		for (k = 0; k < nbCrane; k++)
			epa += yF[k][i];
		c4.add(epa == 1);
		epa.end();
	}
	model.add(c4);
	c4.end();

	//约束（40）
	IloRangeArray  c40(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbTask; i++)
			epa += yF[k][i];
		c40.add(epa >= 1);
		epa.end();
	}
	model.add(c40);
	c40.end();


	////约束（5）///endCzF进一步限定
	//IloRangeArray  c5(env);
	//for (i = 0; i < nbBay; i++)
	//	for (k = 0; k < nbCrane; k++)
	//	{
	//	IloExpr  epa(env);
	//	for (j = 0; j < nbTask; j++)
	//		if (nbLocation[j] - 1 == i)
	//			epa += yF[k][j];
	//	epa -= endCzF[k][i];
	//	c5.add(epa >= 0);
	//	epa.end();
	//	}
	//submodel.add(c5);
	//c5.end();

	////////约束（6）// start_zF取值

	for (j = 0; j < nbTask; j++)
		if (nbLocation[j] < nbBay)
			for (k = 0; k < nbCrane; k++)
			{
		IloExpr  epa(env);
		for (i = nbLocation[j]; i < nbBay; i++)
			epa += start_zF[k][i];
		epa += yF[k][j];

		model.add(epa <= 1);
		epa.end();
			}
	//IloRangeArray  c6(env);
	//for (i = 1; i < nbBay; i++)
	//	for (k = 0; k < nbCrane; k++)
	//	{
	//	IloExpr  epa(env);
	//	for (j = 0; j < nbTask; j++)
	//		if (nbLocation[j] - 1 < i)
	//			epa += yF[k][j];
	//	epa += nbTask * start_zF[k][i];
	//	c6.add(epa <= nbTask);
	//	epa.end();
	//	}
	////for (j = 0; j < nbTask; j++)
	////	for (k = 0; k < nbCrane; k++)
	////	{
	////	IloExpr  epa(env);
	////	for (i = 0; i < nbBay; i++)
	////		if (nbLocation[j] - 1 >= i)
	////			epa += CzF[k][i];
	////	epa -= yF[k][j];
	////	c6.add(epa >= 0);
	////	epa.end();
	////	}
	//model.add(c6);
	//c6.end();

	//for (k = 0; k < nbCrane; k++)
	//	for (i = 0; i < nbBay; i++)
	//		if (i + 1 > nbb[k])
	//			model.add(start_zF[k][i] == 0);






	//约束（7）
	IloRangeArray  c7(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbBay; i++)
			epa += start_zF[k][i];
		c7.add(epa == 1);
		epa.end();

		IloExpr  epa2(env);
		for (i = 0; i < nbBay; i++)
			epa2 += endCzF[k][i];
		c7.add(epa2 == 1);
		epa2.end();
	}
	model.add(c7);
	c7.end();


	//	建立约束(8)
	IloRangeArray  c8(env);
	for (k = 0; k < nbCrane - 1; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbBay; i++)
			epa += i*start_zF[k + 1][i] - i*start_zF[k][i];
		c8.add(epa >= 1 + safe_margin);
		epa.end();

		IloExpr  epa2(env);
		for (i = 0; i < nbBay; i++)
			epa2 += i*endCzF[k + 1][i] - i*endCzF[k][i];
		c8.add(epa2 >= 1 + safe_margin);
		epa2.end();

	}
	model.add(c8);
	c8.end();


	////	建立约束(9)
	//IloRangeArray  c9(env);
	//for (i = 0; i < nbTask - 1; i++)
	//	for (j = i + 1; j < nbTask; j++)
	//		if (nbprecR[i][j] == 1)
	//		{
	//	for (int l = 0; l < nbCrane - 1; l++)
	//		for (k = l + 1; k < nbCrane; k++)
	//		{
	//		IloExpr  epa(env);
	//		epa += yF[l][i] + yF[k][j];
	//		c9.add(epa <= 1);
	//		}

	//		}
	//model.add(c9);
	//c9.end();

	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa1(env);
		IloExpr  epa2(env);
		epa1 += t0kF[k];
		for (i = 0; i < nbBay; i++)
			epa2 += QCmove_time*i*start_zF[k][i];
		epa2 -= QCmove_time*nbb[k] - QCmove_time;
		model.add(epa1 - epa2 >= 0);
		model.add(epa1 + epa2 >= 0);
		epa1.end();
		epa2.end();
	}

	////no overlap
	//for (i = 0; i < nbTask - 1; i++)
	//	for (j = i + 1; j < nbTask; j++)
	//	{
	//	if (nbprecR[i][j] == 1 || nbLocation[i] < nbLocation[j])
	//	{
	//		IloExpr  epa(env);
	//		for (k = 0; k < nbCrane; k++)
	//			epa += yF[k][i] - yF[k][j];
	//		submodel.add(epa <= 0);
	//	}

	//	}



	//IloRangeArray  c00(env);
	//for (k = 0; k < nbCrane; k++)		
	//	{
	//	IloExpr  epa(env);
	//	epa += betaF[k];
	//	for (j = 0; j < nbBay; j++)
	//		epa -= j*endCzF[k][j];

	//	c00.add(epa >= 0);
	//	epa.end();
	//	}
	//model.add(c00);
	//c00.end();

	for (k = 0; k < nbCrane; k++)
	{
		if ((1 + safe_margin)*k + 1 > nbBay)
		{
			for (i = 0; i < nbTask; i++)
				model.add(xF1[k][i] + xF2[k][i] == 0);
		}

		//else if ((1 + safe_margin)*(k+1) + 1 > nbBay && k+1< nbCrane)
		//{
		//	for (i = 0; i < nbTask; i++)
		//		if (nbLocation[i]<=(1 + safe_margin)*k)
		//			model.add(xF[k+1][i] + yF[k+1][i] == 0);
		//}
	}



	//////////////prec 有效不等式
	//for (i = 0; i < nbTask - 2; i++)
	//	for (int j2 = i + 1; j2 < nbTask - 1; j2++)
	//		for (int j3 = j2 + 1; j3 < nbTask; j3++)
	//			if (nbprecR[i][j2] == 1 && nbprecR[j2][j3] == 1)
	//			{
	//	for (k = 1; k < nbCrane; k++)
	//	{
	//		IloExpr  epa(env);
	//		epa += yF[k][j2];
	//		for (int kk = 0; kk < k; kk++)
	//			epa += yF[kk][i] + yF[kk][j3];

	//		model.add(epa <= 2);
	//		epa.end();
	//	}
	//			}


				//**********************************//
				//            有效不等式          //
				//**********************************//

				//变量固定
				for (k = 0; k < nbCrane; k++)
					for (i = 0; i < nbTask; i++)
						model.add(xF1[k][i] + xF2[k][i] - yF[k][i] == 0);

#ifdef precedence_ineq
				for (i = 0; i < nbTask - 1; i++)
					for (j = i + 1; j < nbTask; j++)
						if (nbprecR[i][j] == 1)
						{

					for (k = 1; k < nbCrane; k++)
					{
						IloExpr  epa(env);
						epa += xF2[k][j] - yF[k][j];
						for (int k2 = 0; k2 < k; k2++)
							epa -= yF[k2][i];
						model.add(epa + 1 >= 0);
						epa.end();
					}

					for (k = 1; k < nbCrane; k++)
					{
						IloExpr  epa(env);
						epa += xF1[k][i] - yF[k][i];
						for (int k2 = 0; k2 < k; k2++)
							epa -= yF[k2][j];
						model.add(epa + 1 >= 0);
						epa.end();
					}


					for (k = 0; k < nbCrane; k++)
					{
						IloExpr  epa(env);
						epa += xF1[k][i] - yF[k][i];
						for (int k2 = 0; k2 <= k; k2++)
							epa -= xF1[k2][j];
						model.add(epa + 1 >= 0);
						epa.end();
					}


					for (k = 0; k < nbCrane; k++)
					{
						IloExpr  epa(env);
						epa += xF2[k][j] - yF[k][j];
						for (int k2 = 0; k2 <= k; k2++)
							epa -= xF2[k2][i];
						model.add(epa + 1 >= 0);
						epa.end();
					}



						}


#endif
#ifdef	re_ineq
				//for (k = 0; k < nbCrane - 1; k++)
				//{
				//	for (int bb = 0; bb < nbBay; bb++)
				//	{
				//		IloExpr  epa(env);
				//		epa += nbTask* vF[k];
				//		for (int bb2 = 0; bb2 <= bb; bb2++)
				//			epa += nbTask*yF[k][bb2] ;
				//		for (i = 0; i < nbTask; i++)
				//			if (nbLocation[i] - 1 <= bb)
				//				epa -= xF2[k + 1][i];
				//		model.add(epa >= 0);
				//		epa.end();
				//	}
				//}

				for (k = 0; k < nbCrane - 1; k++)
				{
					for (i = 0; i < nbTask; i++)
					{
						IloExpr  epa(env);
						epa += vF[k] - xF2[k + 1][i];
						for (int bb = 0; bb < nbLocation[i]; bb++)
							epa += yF[k][bb];
						model.add(epa >= 0);
						epa.end();
					}
				}


#endif

	////**********************************//
	////            松弛约束          //
	////**********************************//
	////model.add(vF[nbCrane - 1] - vF[nbCrane - 2] <= 0);


	//for (k = 0; k < nbCrane; k++)
	//	for (i = 0; i < nbBay; i++)
	//		if (i + 1 > nbb[k])
	//			model.add(zF[k][i] == 0);


	IloRangeArray  c20(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		epa += CF - t0kF[k] - gammaF[k];
		for (i = 0; i < nbTask; i++)
			epa -= (nbQ[i] / nbs)*yF[k][i];
		for (j = 0; j < nbBay; j++)
		{
			epa -= QCmove_time*j*endCzF[k][j] - QCmove_time*j*start_zF[k][j];
			epa -= wF[k][j] + CwF[k][j];
		}
		//epa -= thetaF[k];

		////for (j = 0; j < nbLocation[i]; j++)
		//for (j = 0; j < nbBay; j++)
		//	epa += j*zF[k][j];

		c20.add(epa >= 0);
		epa.end();
	}
	model.add(c20);
	c20.end();

	////// thetaF 取值
	////for (k = 0; k < nbCrane; k++)
	////	for (j = 0; j < nbTask; j++)
	////	{
	////	//model.add(thetaF[k] - (nbLocation[j] - 1) * xF[k][j] >= 0);
	////	IloExpr  epa(env);
	////	for (i = 0; i < nbBay; i++)
	////		epa += i*endCzF[k][i];
	////	model.add(epa - (nbLocation[j] - 1) * xF1[k][j] >= 0);
	////	model.add(epa - nbLocation[j] * xF2[k][j] >= 0);
	////	epa.end();
	////	}
	//////for (k = 0; k < nbCrane; k++)
	//////{
	//////	IloExpr  epa(env);
	//////	for (i = 0; i < nbBay; i++)
	//////		epa += i*endCzF[k][i];
	//////	model.add(epa - thetaF[k] == 0);
	//////	epa.end();
	//////}

	//for (k = 0; k < nbCrane; k++)
	//{
	//	IloExpr  epa(env);
	//	for (i = 0; i < nbBay; i++)
	//		epa += i*zF[k][i];
	//	model.add(epa - nbb[k] + 1 <= 0);
	//	epa.end();
	//}


	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbTask; i++)
		{
		model.add(xF2[k][i] - vF[k] <= 0);
		//model.add(nbLocation[i] * yF[k][i] - thetaF[k] <= 0);
		}



	//////约束（4）// 所有任务被分配到QC上
	////IloRangeArray  c4(env);
	////for (i = 0; i < nbTask; i++)
	////{
	////	IloExpr  epa(env);
	////	for (k = 0; k < nbCrane; k++)
	////		epa += xF1[k][i] + xF2[k][i];
	////	c4.add(epa == 1);
	////	epa.end();
	////}
	////submodel.add(c4);
	////c4.end();

	//////约束（5）//endCzF 取在最后一个xF处
	////IloRangeArray  c5(env);
	////for (i = 0; i < nbBay; i++)
	////	for (k = 0; k < nbCrane; k++)
	////	{
	////	IloExpr  epa(env);
	////	for (j = 0; j < nbTask; j++)
	////		if (nbLocation[j] - 1 == i)
	////			epa += xF[k][j];
	////	epa -= endCzF[k][i];
	////	c5.add(epa >= 0);
	////	epa.end();
	////	}
	////model.add(c5);
	////c5.end();

	////约束（6）// zF 小于最小的xF的bay
	//IloRangeArray  c60(env);
	//for (i = 0; i < nbBay; i++)
	//	for (k = 0; k < nbCrane; k++)
	//	{
	//	IloExpr  epa(env);
	//	IloExpr  epa2(env);
	//	for (j = 0; j < nbTask; j++)
	//		if (nbLocation[j] - 1 < i)
	//		{
	//		epa += xF2[k][j];
	//		epa2 += xF1[k][j];
	//		}
	//	for (j = i; j < nbBay; j++)
	//	{
	//		epa += 100 * CzF[k][j];
	//		epa2 += 100 * zF[k][j];
	//	}
	//	c60.add(epa <= 100);
	//	c60.add(epa2 <= 100);
	//	epa.end();
	//	epa2.end();
	//	}
	//model.add(c60);
	//c60.end();
	for (j = 0; j < nbTask; j++)
		if (nbLocation[j] < nbBay)
			for (k = 0; k < nbCrane; k++)
			{
		//IloExpr  epa(env);
		IloExpr  epa2(env);

		//for (i = nbLocation[j]; i < nbBay; i++)
		//	epa += zF[k][i];
		//epa += xF1[k][j];

		for (i = nbLocation[j]; i < nbBay; i++)
			epa2 += CzF[k][i];
		epa2 += xF2[k][j];

		//model.add(epa <= 1);
		model.add(epa2 <= 1);
		//epa.end();
		epa2.end();
			}



	//约束（7） zF unique
	IloRangeArray  c70(env);
	for (k = 0; k < nbCrane; k++)
	{
		//IloExpr  epa(env);
		IloExpr  epa2(env);
		//IloExpr  epa3(env);
		for (i = 0; i < nbBay; i++)
		{
			//epa += zF[k][i];
			epa2 += CzF[k][i];
			//epa3 += endCzF[k][i];
		}
		//epa2 -= vF[k];
		c70.add(epa2 <= 1);
		//c70.add(epa2 == 0);
		//c7.add(epa3 == 1);
		//epa.end();
		epa2.end();
		//epa3.end();
	}
	model.add(c70);
	c70.end();


	//	建立约束(8)  z 和 z之间隔 /delta +1
	IloRangeArray  c80(env);
	for (k = 0; k < nbCrane - 1; k++)
	{
		//IloExpr  epa(env);
		//for (i = 0; i < nbBay; i++)
		//	epa += i*zF[k + 1][i] - i*zF[k][i];
		//c80.add(epa >= 2);
		//epa.end();

		IloExpr  epa2(env);
		for (i = 0; i < nbBay; i++)
			epa2 += i*CzF[k + 1][i] - i*CzF[k][i];
		//epa2 -= nbBay * (vF[k] + vF[k + 1] - 2);
		c80.add(epa2 >= 1 + safe_margin);
		epa2.end();

	}
	model.add(c80);
	c80.end();

	//for (k = 0; k < nbCrane; k++)
	//{
	//	IloExpr  epa1(env);
	//	IloExpr  epa2(env);
	//	epa1 += t0kF[k];
	//	for (i = 0; i < nbBay; i++)
	//		epa2 += i*zF[k][i];
	//	model.add(epa1 + epa2 - nbb[k] + 1 >= 0);
	//	model.add(epa1 - epa2 + nbb[k] - 1 >= 0);
	//	epa1.end();
	//	epa2.end();
	//}


	////QC travel limits
	//IloRangeArray  v000(env);
	//for (k = 0; k < nbCrane; k++)
	//{
	//	IloExpr  epa(env);
	//	for (i = 0; i < nbTask; i++)
	//	{
	//		if (nbLocation[i] < 2 * k + 1)
	//			epa += xF1[k][i] + xF2[k][i];
	//		if (nbLocation[i] > nbBay - 2 * (nbCrane - k - 1))
	//			epa += xF1[k][i] + xF2[k][i];
	//	}
	//	v000.add(epa <= 0);
	//	epa.end();
	//}
	//model.add(v000);
	//v000.end();



	//////////////////sub-problem

	// vF 取值(2)
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbTask; i++)
			epa += xF2[k][i];
		epa -= vF[k];
		model.add(epa >= 0);
		epa.end();
	}


	//for (k = 0; k < nbCrane; k++)
	//{
	//	IloExpr  epa(env);
	//	for (i = 0; i < nbTask; i++)
	//		epa += xF1[k][i];
	//	model.add(epa >= 1);
	//	epa.end();
	//}
	//// vF 取值(3) vF与violation关系
	//for (i = 0; i < nbTask - 1; i++)
	//	for (j = i + 1; j < nbTask; j++)
	//		if (nbprecR[i][j] == 1)
	//		{
	//	for (k = 0; k < nbCrane; k++)
	//	{
	//		IloExpr epa(env);
	//		IloExpr epa2(env);

	//		epa += vF[k] + 1;
	//		epa2 += vF[k] + 1;

	//		if (k < nbCrane - 1)
	//		{
	//			epa -= xF1[k][i] + xF2[k][i];
	//			for (int kk = k + 1; kk < nbCrane; kk++)
	//				epa -= xF1[kk][j] + xF2[kk][j];
	//			model.add(epa >= 0);
	//		}

	//		if (k > 0)
	//		{
	//			epa2 -= xF1[k][j] + xF2[k][j];
	//			for (int kk = 0; kk < k; kk++)
	//				epa2 -= xF1[kk][i] + xF2[kk][i];
	//			model.add(epa2 >= 0);
	//		}
	//		epa.end();
	//		epa2.end();
	//	}
	//		}


	///////precedence
	for (i = 0; i < nbTask - 1; i++)
		for (j = i + 1; j < nbTask; j++)
			if (nbprecR[i][j] == 1)
			{
		for (k = 0; k < nbCrane; k++)
		{
			//1
			if (k < nbCrane - 1)
			{
				IloExpr epa3(env);
				epa3 += xF1[k][i];
				for (int kk = k + 1; kk < nbCrane; kk++)
					epa3 += xF1[kk][j];
				//epa3 += xF[kk][j] + yF[kk][j];
				model.add(epa3 <= 1);
				epa3.end();
			}

			//3
			{
				IloExpr epa(env);
				epa += xF2[k][i];
				for (int kk = k; kk < nbCrane; kk++)
					epa -= xF2[kk][j];
				model.add(epa <= 0);
				epa.end();
			}
			//2
			{
				IloExpr epa2(env);
				epa2 += xF1[k][j];
				for (int kk = k; kk < nbCrane; kk++)
					epa2 -= xF1[kk][i];
				model.add(epa2 <= 0);
				epa2.end();
			}

			//4
			if (k < nbCrane - 1)
			{
				IloExpr epa4(env);
				epa4 += xF2[k][j];

				//for (int kk = k + 1; kk < nbCrane; kk++)
				//	epa4 += yF[kk][i];
				//epa4 -= 1;


				for (int kk = k + 1; kk < nbCrane; kk++)
					epa4 += xF2[kk][i];
				//for (int kk = 0; kk <= k; kk++)
				//	epa4 -= xF1[kk][i] + xF2[kk][i];
				model.add(epa4 <= 1);
				epa4.end();
			}

			//很有效果
			for (int kk = 0; kk <= k; kk++)
			{
				model.add(xF2[k][j] + xF1[kk][i] - vF[kk] <= 1);
			}

		}
			}

	////back time gammaF
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		epa += gammaF[k];
		for (i = 0; i < nbBay; i++)
			epa += QCmove_time*i*CzF[k][i] - QCmove_time*i*endCzF[k][i];
		//epa += nbBay - nbBay*vF[k];
		model.add(epa >= 0);
		epa.end();
	}

	////////////////////////////////////////
	////////waiting time ///////
	///////////////////////////////////////

	///////forward trip
	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbBay; i++)
		{
		IloExpr  epa(env);
		epa += t0kF[k];
		for (j = 0; j < nbTask; j++)
			if (nbLocation[j] - 1 <= i)
				epa += (nbQ[j] / nbs)*xF1[k][j];
		for (j = 0; j <= i; j++)
		{

			epa += wF[k][j];
			epa -= QCmove_time*j*start_zF[k][j];
			//epa -= j*zF[k][j];
			//epa -= 1000 * CzF[k][j];
		}
		model.add(epa + QCmove_time*i - TF[k][i] == 0);
		epa.end();

		}

	/////约束（3）
	IloRangeArray  c30(env);
	for (i = 0; i < nbBay - 1 - safe_margin; i++)
		for (k = 0; k < nbCrane - 1; k++)
		{
		IloExpr  epa(env);

		epa += TF[k][i] - TF[k + 1][i + 1 + safe_margin];
		epa += 1600;
		for (j = 0; j <= i + 1 + safe_margin; j++)
			epa -= 1600 * start_zF[k + 1][j];
		for (j = 0; j < i; j++)
			epa += 1600 * endCzF[k][j];

		c30.add(epa >= 0);
		epa.end();
		}
	model.add(c30);
	c30.end();

	//////retracing trip
	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbBay; i++)
		{
		IloExpr  epa(env);
		epa += t0kF[k];
		for (j = 0; j < nbTask; j++)
			epa += (nbQ[j] / nbs)*xF1[k][j];
		for (j = 0; j < nbBay; j++)
		{
			epa += 2 * QCmove_time*j*endCzF[k][j] - QCmove_time*j*start_zF[k][j];
			epa += wF[k][j];
		}


		for (j = 0; j < nbTask; j++)
			if (nbLocation[j] - 1 >= i)
				epa += (nbQ[j] / nbs)*xF2[k][j];
		for (j = i; j < nbBay; j++)
		{
			epa += CwF[k][j];
		}

		model.add(epa - QCmove_time*i - CTF[k][i] == 0);
		epa.end();

		}

	for (k = 0; k < nbCrane - 1; k++)
		for (i = 0; i < nbBay - 1 - safe_margin; i++)
		{
		IloExpr  epa(env);
		epa += CTF[k + 1][i + 1 + safe_margin] - CTF[k][i];

		epa += 4500 - 1500 * vF[k] - 1500 * vF[k + 1];
		for (j = 0; j < i; j++)
			epa += 1500 * endCzF[k][j];
		for (j = 0; j <= i + 1 + safe_margin; j++)
			epa -= 1500 * CzF[k + 1][j];

		model.add(epa >= 0);
		epa.end();
		}

	// forward and retrace
	for (k = 0; k < nbCrane - 1; k++)
	{
		IloExpr  epa(env);
		for (j = 0; j < nbBay; j++)
			epa += j * CzF[k + 1][j] - j * endCzF[k][j];
		epa += 1500 + 1500 * vF[k] - 1500 * vF[k + 1];
		model.add(epa >= 1 + safe_margin);
		epa.end();

	}



	return TRUE;
}

BOOL SP_1(IloModel model, IloNum nbs, IloIntArray  nbb, IloNumArray nbQ, IloIntArray  nbLocation, BoolMatrix nbprecR, IloIntArray  nbreadyT,
	BoolVarMatrix xF, BoolVarMatrix yF, BoolVarMatrix zF, IloIntVar CF, IloNumVarArray QC_CF, IloIntVarArray thetaF, IloBoolVarArray vF,
	BoolMatrix yF_best, BoolMatrix CxF_best, IloNum UB, IloNum *ObjVal, NumMatrix wF_best, NumMatrix CwF_best)
{

	//xF: forward trip, yF: retracing trip


	IloEnv env = model.getEnv();
	IloInt i, j, k;

	//问题模型
	IloModel submodel(env);

	//问题模型
	IloNumVarArray t0kF(env, nbCrane, 0, 100);
	IloNumVarArray gammaF(env, nbCrane, 0, 100);//

	IloNumVarArray ckF(env, nbCrane, 0, 3000);//
	IloNumVarArray chaF(env, nbCrane, 0, 3000);//

	NumVarMatrix wF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		wF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);

	NumVarMatrix CwF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		CwF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);

	BoolVarMatrix CzF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		CzF[k] = IloBoolVarArray(env, nbBay);
	}

	BoolVarMatrix endCzF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		endCzF[k] = IloBoolVarArray(env, nbBay);
	}

	//leaving time at each bay
	NumVarMatrix TF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		TF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);
	NumVarMatrix CTF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		CTF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);


	//**********************************//
	//           objective 原问题目标函数        //
	//**********************************//

	//IloExpr obj1(env); 


	IloExpr  obj2(env);
	//  建立子问题目标函数表达式 

	obj2 += CF;

	//for (k = 0; k < nbCrane; k++)
	//	obj2 += chaF[k];

	//	将目标函数加入到原问题模型
	submodel.add(IloMinimize(env, obj2));//
	//obj1.end();
	obj2.end();


	//**********************************//
	//            variable fixing           //
	//**********************************//
	for (i = 0; i < nbTask - 1; i++)
		for (j = i + 1; j < nbTask; j++)
			if (nbprecR[i][j] == 1)
			{
		int i_k, j_k;
		for (k = 0; k < nbCrane; k++)
		{
			if (yF_best[k][i] == 1)
				i_k = k;
			if (yF_best[k][j] == 1)
				j_k = k;
		}

		if (i_k>j_k)
		{
			submodel.add(xF[i_k][i] == 1);
			submodel.add(yF[i_k][i] == 0);
		}
		else if (i_k < j_k)
		{
			submodel.add(xF[j_k][j] == 0);
			submodel.add(yF[j_k][j] == 1);
		}
			}


	//**********************************//
	//            Constraints           //
	//**********************************//
	//model.add(vF[nbCrane - 1] - vF[nbCrane - 2] <= 0);

	//变量固定
	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbTask; i++)
			submodel.add(xF[k][i] + yF[k][i] == yF_best[k][i]);

	for (k = 0; k < nbCrane; k++)
		submodel.add(CF - ckF[k] >= 0);

	//for (k = 0; k < nbCrane; k++)
	//	submodel.add(chaF[k] - ckF[k] + UB>= 0);

	IloRangeArray  c2(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		//epa += CF - t0kF[k] - gammaF[k];

		epa += ckF[k] - t0kF[k] - gammaF[k];

		for (i = 0; i < nbTask; i++)
			epa -= (nbQ[i] / nbs)*(xF[k][i] + yF[k][i]);
		for (j = 0; j < nbBay; j++)
		{
			epa -= j*endCzF[k][j];
			epa -= wF[k][j] + CwF[k][j];
		}
		//epa -= thetaF[k];

		//for (j = 0; j < nbLocation[i]; j++)
		for (j = 0; j < nbBay; j++)
			epa += j*zF[k][j];

		c2.add(epa >= 0);
		epa.end();
	}
	submodel.add(c2);
	c2.end();

	// thetaF 取值
	for (k = 0; k < nbCrane; k++)
		for (j = 0; j < nbTask; j++)
		{
		//model.add(thetaF[k] - (nbLocation[j] - 1) * xF[k][j] >= 0);
		IloExpr  epa(env);
		for (i = 0; i < nbBay; i++)
			epa += i*endCzF[k][i];
		submodel.add(epa - (nbLocation[j] - 1) * xF[k][j] >= 0);
		submodel.add(epa - nbLocation[j] * yF[k][j] >= 0);
		epa.end();
		}
	//for (k = 0; k < nbCrane; k++)
	//{
	//	IloExpr  epa(env);
	//	for (i = 0; i < nbBay; i++)
	//		epa += i*endCzF[k][i];
	//	model.add(epa - thetaF[k] == 0);
	//	epa.end();
	//}

	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbBay; i++)
			epa += i*zF[k][i];
		submodel.add(epa - nbb[k] + 1 <= 0);
		epa.end();
	}


	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbTask; i++)
		{
		submodel.add(yF[k][i] - vF[k] <= 0);
		//model.add(nbLocation[i] * yF[k][i] - thetaF[k] <= 0);
		}
	for (k = 0; k < nbCrane - 1; k++)
	{
		//model.add(thetaF[k + 1] - thetaF[k] >= 2);
		IloExpr  epa(env);
		for (i = 0; i < nbBay; i++)
			epa += i*endCzF[k + 1][i] - i*endCzF[k][i];
		submodel.add(epa >= 2);
		epa.end();
	}


	//约束（4）// 所有任务被分配到QC上
	IloRangeArray  c4(env);
	for (i = 0; i < nbTask; i++)
	{
		IloExpr  epa(env);
		for (k = 0; k < nbCrane; k++)
			epa += xF[k][i] + yF[k][i];
		c4.add(epa == 1);
		epa.end();
	}
	submodel.add(c4);
	c4.end();

	////约束（5）//endCzF 取在最后一个xF处
	//IloRangeArray  c5(env);
	//for (i = 0; i < nbBay; i++)
	//	for (k = 0; k < nbCrane; k++)
	//	{
	//	IloExpr  epa(env);
	//	for (j = 0; j < nbTask; j++)
	//		if (nbLocation[j] - 1 == i)
	//			epa += xF[k][j];
	//	epa -= endCzF[k][i];
	//	c5.add(epa >= 0);
	//	epa.end();
	//	}
	//model.add(c5);
	//c5.end();

	//约束（6）// zF 小于最小的xF的bay
	IloRangeArray  c6(env);
	for (i = 0; i < nbBay; i++)
		for (k = 0; k < nbCrane; k++)
		{
		IloExpr  epa(env);
		IloExpr  epa2(env);
		for (j = 0; j < nbTask; j++)
			if (nbLocation[j] - 1 < i)
			{
			epa += yF[k][j];
			epa2 += xF[k][j];
			}
		for (j = i; j < nbBay; j++)
		{
			epa += 100 * CzF[k][j];
			epa2 += 100 * zF[k][j];
		}


		c6.add(epa <= 100);
		c6.add(epa2 <= 100);
		epa.end();
		epa2.end();
		}
	submodel.add(c6);
	c6.end();



	//约束（7） zF unique
	IloRangeArray  c7(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		IloExpr  epa2(env);
		IloExpr  epa3(env);
		for (i = 0; i < nbBay; i++)
		{
			epa += zF[k][i];
			epa2 += CzF[k][i];
			epa3 += endCzF[k][i];
		}
		epa2 -= vF[k];
		c7.add(epa == 1);
		c7.add(epa2 == 0);
		c7.add(epa3 == 1);
		epa.end();
		epa2.end();
		epa3.end();
	}
	submodel.add(c7);
	c7.end();


	//	建立约束(8)  z 和 z之间隔 /delta +1
	IloRangeArray  c8(env);
	for (k = 0; k < nbCrane - 1; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbBay; i++)
			epa += i*zF[k + 1][i] - i*zF[k][i];
		c8.add(epa >= 2);
		epa.end();

		IloExpr  epa2(env);
		for (i = 0; i < nbBay; i++)
			epa2 += i*CzF[k + 1][i] - i*CzF[k][i];
		epa2 -= nbBay * (vF[k] + vF[k + 1] - 2);
		c8.add(epa2 >= 2);
		epa2.end();

	}
	submodel.add(c8);
	c8.end();

	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa1(env);
		IloExpr  epa2(env);
		epa1 += t0kF[k];
		for (i = 0; i < nbBay; i++)
			epa2 += i*zF[k][i];
		submodel.add(epa1 + epa2 - nbb[k] + 1 >= 0);
		submodel.add(epa1 - epa2 + nbb[k] - 1 >= 0);
		epa1.end();
		epa2.end();
	}


	//QC travel limits
	IloRangeArray  v00(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbTask; i++)
		{
			if (nbLocation[i] < 2 * k + 1)
				epa += xF[k][i] + yF[k][i];
			if (nbLocation[i] > nbBay - 2 * (nbCrane - k - 1))
				epa += xF[k][i] + yF[k][i];
		}
		v00.add(epa <= 0);
		epa.end();
	}
	submodel.add(v00);
	v00.end();



	////////////////sub-problem

	// vF 取值(2)
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbTask; i++)
			epa += yF[k][i];
		epa -= vF[k];
		submodel.add(epa >= 0);
		epa.end();
	}
	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbTask; i++)
			submodel.add(yF[k][i] - vF[k] <= 0);

	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbTask; i++)
			epa += xF[k][i];
		submodel.add(epa >= 1);
		epa.end();
	}
	// vF 取值(3) vF与violation关系
	for (i = 0; i < nbTask - 1; i++)
		for (j = i + 1; j < nbTask; j++)
			if (nbprecR[i][j] == 1)
			{
		for (k = 0; k < nbCrane; k++)
		{
			IloExpr epa(env);
			IloExpr epa2(env);

			epa += vF[k] + 1;
			epa2 += vF[k] + 1;

			if (k < nbCrane - 1)
			{
				epa -= xF[k][i] + yF[k][i];
				for (int kk = k + 1; kk < nbCrane; kk++)
					epa -= xF[kk][j] + yF[kk][j];
				submodel.add(epa >= 0);
			}

			if (k > 0)
			{
				epa2 -= xF[k][j] + yF[k][j];
				for (int kk = 0; kk < k; kk++)
					epa2 -= xF[kk][i] + yF[kk][i];
				submodel.add(epa2 >= 0);
			}
			epa.end();
			epa2.end();
		}
			}


	//precedence
	for (i = 0; i < nbTask - 1; i++)
		for (j = i + 1; j < nbTask; j++)
			if (nbprecR[i][j] == 1)
			{
		for (k = 0; k < nbCrane; k++)
		{

			{
				IloExpr epa(env);
				epa += yF[k][i];
				for (int kk = k; kk < nbCrane; kk++)
					epa -= yF[kk][j];
				submodel.add(epa <= 0);
				epa.end();
			}

			{
				IloExpr epa2(env);
				epa2 += xF[k][j];
				for (int kk = k; kk < nbCrane; kk++)
					epa2 -= xF[kk][i];
				submodel.add(epa2 <= 0);
				epa2.end();
			}

			if (k < nbCrane - 1)
			{
				IloExpr epa3(env);
				epa3 += xF[k][i];
				for (int kk = k + 1; kk < nbCrane; kk++)
					epa3 += xF[kk][j];
				//epa3 += xF[kk][j] + yF[kk][j];
				submodel.add(epa3 <= 1);
				epa3.end();
			}


			if (k < nbCrane - 1)
			{
				IloExpr epa4(env);
				epa4 += yF[k][j];

				//for (int kk = k + 1; kk < nbCrane; kk++)
				//	epa4 += yF[kk][i];
				//epa4 -= 1;


				for (int kk = k + 1; kk < nbCrane; kk++)
					epa4 -= xF[kk][i];
				for (int kk = 0; kk <= k; kk++)
					epa4 -= xF[kk][i] + yF[kk][i];
				submodel.add(epa4 <= 0);
				epa4.end();
			}

			for (int kk = 0; kk <= k; kk++)
			{
				submodel.add(yF[k][j] + xF[kk][i] - vF[kk] <= 1);
			}

		}
			}

	////back time gammaF
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		epa += gammaF[k];
		for (i = 0; i < nbBay; i++)
			epa += i*CzF[k][i] - i*endCzF[k][i];
		epa += nbBay - nbBay*vF[k];
		submodel.add(epa >= 0);
		epa.end();
	}

	//////////////////////////////////////
	//////waiting time ///////
	/////////////////////////////////////

	///////forward trip
	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbBay; i++)
		{
		IloExpr  epa(env);
		epa += t0kF[k];
		for (j = 0; j < nbTask; j++)
			if (nbLocation[j] - 1 <= i)
				epa += (nbQ[j] / nbs)*xF[k][j];
		for (j = 0; j <= i; j++)
		{

			epa += wF[k][j];
			epa -= j*zF[k][j];
			//epa -= 1000 * CzF[k][j];
		}
		submodel.add(epa + i - TF[k][i] == 0);
		epa.end();

		}

	/////约束（3）
	IloRangeArray  c3(env);
	for (i = 0; i < nbBay - 2; i++)
		for (k = 0; k < nbCrane - 1; k++)
		{
		IloExpr  epa(env);

		epa += TF[k][i] - TF[k + 1][i + 2];
		epa += 1600;
		for (j = 0; j <= i + 2; j++)
			epa -= 1600 * zF[k + 1][j];
		for (j = 0; j < i; j++)
			epa += 1600 * endCzF[k][j];

		c3.add(epa >= 0);
		epa.end();
		}
	submodel.add(c3);
	c3.end();

	//////retracing trip
	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbBay; i++)
		{
		IloExpr  epa(env);
		epa += t0kF[k];
		for (j = 0; j < nbTask; j++)
			epa += (nbQ[j] / nbs)*xF[k][j];
		for (j = 0; j < nbBay; j++)
		{
			epa += 2 * j*endCzF[k][j] - j*zF[k][j];
			epa += wF[k][j];
		}


		for (j = 0; j < nbTask; j++)
			if (nbLocation[j] - 1 >= i)
				epa += (nbQ[j] / nbs)*yF[k][j];
		for (j = i; j < nbBay; j++)
		{
			epa += CwF[k][j];
		}

		submodel.add(epa - i - CTF[k][i] == 0);
		epa.end();

		}

	for (k = 0; k < nbCrane - 1; k++)
		for (i = 0; i < nbBay - 2; i++)
		{
		IloExpr  epa(env);
		epa += CTF[k + 1][i + 2] - CTF[k][i];

		epa += 6000 - 2000 * vF[k] - 2000 * vF[k + 1];
		for (j = 0; j < i; j++)
			epa += 2000 * endCzF[k][j];
		for (j = 0; j <= i + 2; j++)
			epa -= 2000 * CzF[k + 1][j];

		submodel.add(epa >= 0);
		epa.end();
		}

	// forward and retrace
	for (k = 0; k < nbCrane - 1; k++)
	{
		IloExpr  epa(env);
		for (j = 0; j < nbBay; j++)
			epa += j * CzF[k + 1][j] - j * endCzF[k][j];
		epa += 2000 + 2000 * vF[k] - 2000 * vF[k + 1];
		submodel.add(epa >= 2);
		epa.end();

	}

	//**********************************//
	//            开始求解		        //
	//**********************************//
	IloCplex cplex(env);
	cplex.extract(submodel);

#ifdef UserActive
	//cplex.use(CandyUserCallback(env));
	//cplex.setParam(IloCplex::MIPEmphasis, 2);
	//cplex.setParam(IloCplex::HeurFreq, 1);
	//cplex.setParam(IloCplex::ParallelMode, 1);
	cplex.setParam(IloCplex::Threads, 4);
#endif



	//cplex.setParam(IloCplex::MIPEmphasis ,CPX_MIPEMPHASIS_HIDDENFEAS);
	//cplex.setParam(IloCplex::EpGap,0.01);
	//cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::TiLim, 7200);


	//cout<<endl<<"number of binary variables = "<<cplex.getNbinVars()<<endl;
	//cout<<endl<<"number of constraints = "<<cplex.getNrows()<<endl;

	BOOL h0;
	h0 = cplex.solve();
	//cout<<"h0 "<<h0<<endl;
	if (!h0)
	{
		cout << "\nno feasible solution has been found" << endl;
		cplex.clearModel();
		cplex.clear();
		cplex.end();
		submodel.end();

		return FALSE;
	}
	//**********************************//
	//             记录最好解           //
	//**********************************//

	//*ObjVal = cplex.getBestObjValue();
	*ObjVal = cplex.getValue(CF);

	for (k = 0; k < nbCrane; k++)
	{
		for (i = 0; i < nbTask; i++)
		{
			if (cplex.isExtracted(xF[k][i]))
			{
				CxF_best[k][i] = cplex.getValue(xF[k][i]);
				//cout<< wF_best[k][i] << "; ";
			}
			else
				CxF_best[k][i] = 0;

		}//cout << endl;
	}
	//cout << endl << endl;

	cout << "final3:  " << *ObjVal<<endl;
	for (k = 0; k < nbCrane; k++)
	{
		cout << k << ":  ";
		for (i = 0; i < nbTask; i++)
		{
			if (yF_best[k][i]==1)
			{
				if (CxF_best[k][i] == 1)
					cout<< i << " ";
				else
					cout << i << "r ";
			}


		}cout << endl;
	}
	cout  << endl;

	int recb = 1;
	cout << recb << ":  ";
	for (i = 0; i < nbTask; i++)
	{
		if (nbLocation[i] != recb)
		{
			recb = nbLocation[i];
			cout << endl << recb << ":  ";			
		}
		cout << i << "-";
		for (k = 0; k < nbCrane; k++)
			if (yF_best[k][i] == 1)
				cout << k << " ";
	}
	cout << endl;

	//cout << "QC_CF_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	cout << cplex.getValue(QC_CF[k]) << "  ";
	//}
	//cout << endl << endl;

	//cout << "Task_CF_best" << endl;
	//for (i = 0; i < nbTask; i++)
	//{
	//	cout << cplex.getValue(Task_CF[i]) << "  ";
	//}
	//cout << endl << endl;

	for (k = 0; k < nbCrane; k++)
	{
		cout << k << ":  ";
		for (i = 0; i < nbBay; i++)
		{
			if (cplex.isExtracted(wF[k][i]))
			{
				
				wF_best[k][i] = cplex.getValue(wF[k][i]);
				cout<< wF_best[k][i] << "; ";
			}
			else
				wF_best[k][i] = 0;

		}cout << endl;
	}
	cout << endl;

	for (k = 0; k < nbCrane; k++)
	{
		cout << k << ":  ";
		for (i = 0; i < nbBay; i++)
		{
			if (cplex.isExtracted(CwF[k][i]))
			{
				CwF_best[k][i] = cplex.getValue(CwF[k][i]);
				cout<< CwF_best[k][i] << "; ";
			}
			else
				CwF_best[k][i] = 0;

		}cout << endl;
	}
	cout << endl;


	for (k = 0; k < nbCrane; k++)
	{
		cout << k << ":  "<<cplex.getValue(ckF[k]) << endl;
	}
	cout << endl;


	//cout << "CF_best: " << *CF_best << endl;
	cout << "ObjVal: " << *ObjVal << endl;

	cplex.clearModel();
	cplex.clear();
	cplex.end();
	submodel.end();

	return TRUE;

}

BOOL Chen_Same_Direction(IloModel model, IloNum nbs, IloIntArray  nbb, IloNumArray nbQ, IloIntArray  nbreadyT,
	BoolVarMatrix yF, IloIntVar CF, IloInt *CF_best, BoolMatrix yF_best,
	IloNum *ObjVal, IloNum *UB, IloIntArray  nbLocation, BoolMatrix nbprecR, int my)

{
	IloEnv env = model.getEnv();
	IloInt i, j, k;


	//问题模型
	IloModel submodel(env);
	BoolVarMatrix CuF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		CuF[k] = IloBoolVarArray(env, nbBay);
	}
	BoolVarMatrix CvF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		CvF[k] = IloBoolVarArray(env, nbBay);
	}
	NumVarMatrix mF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		mF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);

	IloIntVarArray startF(env, nbCrane, 0, 100);
	IloIntVarArray endF(env, nbCrane, 0, 100);

	IloNumVarArray t0kF(env, nbCrane, 0, 100);

	//leaving time at each bay
	NumVarMatrix TF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		TF[k] = IloNumVarArray(env, nbBay, -100, IloInfinity);

	//**********************************//
	//            原问题目标函数        //
	//**********************************//

	//IloExpr obj1(env); 


	IloExpr  obj2(env);
	//  建立子问题目标函数表达式 

	obj2 += CF;

	//	将目标函数加入到原问题模型
	submodel.add(IloMinimize(env, obj2));//
	//obj1.end();
	obj2.end();


	//**********************************//
	//            MP问题 约束           //
	//**********************************//

	//IloRangeArray  c0(env);
	//for (k = 0; k < nbCrane; k++)
	//{
	//	IloExpr  epa(env);
	//	for (j = 0; j < nbTask; j++)
	//		if (nbLocation[j] >(nbBay - 2 * (nbCrane - k - 1)))
	//			epa += yF[k][j];
	//	c0.add(epa == 0);
	//	epa.end();
	//}
	//submodel.add(c0);
	//c0.end();


	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbBay; i++)
		{
		IloExpr  epa(env);
		epa += nbreadyT[k] + t0kF[k] - QCmove_time*startF[k];
		for (j = 0; j < nbTask; j++)
			if (nbLocation[j] - 1 <= i)
				epa += (nbQ[j] / nbs)*yF[k][j];
		for (j = 0; j <= i; j++)
		{

			epa += mF[k][j];
		}
		submodel.add(epa + QCmove_time*i - TF[k][i] == 0);
		epa.end();

		}

	IloRangeArray  c2(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		epa += CF - t0kF[k] - nbreadyT[k];
		for (i = 0; i < nbTask; i++)
			epa -= (nbQ[i] / nbs)*yF[k][i];
		for (i = 0; i < nbBay; i++)
			epa -= mF[k][i];


		epa += QCmove_time*startF[k];
		epa -= QCmove_time*endF[k];

		c2.add(epa >= 0);
		epa.end();
	}
	submodel.add(c2);
	c2.end();

	//约束（3）
	IloRangeArray  c3(env);
	for (i = 0; i < nbBay - 2; i++)
		for (k = 0; k < nbCrane - 1; k++)
		{
		IloExpr  epa(env);
		epa += TF[k][i] - TF[k + 1][i + 2];

		epa += 2000 * (CvF[k][i] + CuF[k + 1][i + 2]);
		//epa += 2000 * (CuF[k][i] + CuF[k + 1][i + 2]);

		c3.add(epa >= 0);
		epa.end();
		}
	submodel.add(c3);
	c3.end();

	//约束（4）
	IloRangeArray  c4(env);
	for (i = 0; i < nbTask; i++)
	{
		IloExpr  epa(env);
		for (k = 0; k < nbCrane; k++)
			epa += yF[k][i];
		c4.add(epa == 1);
		epa.end();
	}
	submodel.add(c4);
	c4.end();

	////约束（40）
	//IloRangeArray  c40(env);
	//for (k = 0; k < nbCrane; k++)
	//{
	//	IloExpr  epa(env);
	//	for (i = 0; i < nbTask; i++)
	//		epa += yF[k][i];
	//	c40.add(epa >= 1);
	//	epa.end();
	//}
	//submodel.add(c40);
	//c40.end();

	//约束（5）
	IloRangeArray  c5(env);
	for (i = 0; i < nbTask; i++)
		for (k = 0; k < nbCrane; k++)
		{
		IloExpr  epa(env);
		epa += (nbLocation[i] - 1) * yF[k][i];
		epa -= startF[k];
		epa += 2 * nbBay * (1 - yF[k][i]);
		c5.add(epa >= 0);
		epa.end();
		}
	submodel.add(c5);
	c5.end();

	//约束（6）
	IloRangeArray  c6(env);
	for (i = 0; i < nbTask; i++)
		for (k = 0; k < nbCrane; k++)
		{
		IloExpr  epa(env);
		epa -= (nbLocation[i] - 1) * yF[k][i];
		epa += endF[k];
		c6.add(epa >= 0);
		epa.end();
		}
	submodel.add(c6);
	c6.end();



	//	建立约束(9)
	IloRangeArray  c9(env);
	for (i = 0; i < nbTask; i++)
		for (j = 0; j < nbTask; j++)
			if (nbprecR[i][j] == 1)
			{
		for (int l = 0; l < nbCrane - 1; l++)
			for (k = l + 1; k < nbCrane; k++)
			{
			IloExpr  epa(env);
			epa += yF[l][i] + yF[k][j];
			c9.add(epa <= 1);
			}

			}
	submodel.add(c9);
	c9.end();


	//约束（5）
	IloRangeArray  c10a(env);
	for (i = 0; i < nbBay; i++)
		for (k = 0; k < nbCrane; k++)
		{
		IloExpr  epa1(env);
		IloExpr  epa2(env);
		IloExpr  epa3(env);
		epa1 += i - startF[k];
		epa2 += 2 * nbBay * (1 - CuF[k][i]);
		epa3 += 2 * nbBay * CuF[k][i];
		c10a.add(epa1 - epa2 + 0.1 <= 0);
		c10a.add(epa3 + epa1 >= 0);
		epa1.end();
		epa2.end();
		epa3.end();
		}
	submodel.add(c10a);
	c10a.end();


	//约束（5）
	IloRangeArray  c10a2(env);
	for (i = 0; i < nbBay; i++)
		for (k = 0; k < nbCrane; k++)
		{
		IloExpr  epa1(env);
		IloExpr  epa2(env);
		IloExpr  epa3(env);
		epa1 += endF[k] - i;
		epa2 += 2 * nbBay * (1 - CvF[k][i]);
		epa3 += 2 * nbBay * CvF[k][i];
		c10a2.add(epa1 - epa2 + 0.1 <= 0);
		c10a2.add(epa3 + epa1 >= 0);
		epa1.end();
		epa2.end();
		epa3.end();
		}
	submodel.add(c10a2);
	c10a2.end();

	for (k = 0; k < nbCrane - 1; k++)
	{
		IloExpr  epa1(env);
		IloExpr  epa2(env);
		epa1 += startF[k + 1] - startF[k];
		epa2 += endF[k + 1] - endF[k];
		submodel.add(epa1 >= 1 + safe_margin);
		submodel.add(epa2 >= 1 + safe_margin);
		epa1.end();
		epa2.end();
	}

	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa1(env);
		IloExpr  epa2(env);
		epa1 += t0kF[k];
		epa2 += QCmove_time*startF[k] - QCmove_time*nbb[k] + QCmove_time;
		submodel.add(epa1 - epa2 >= 0);
		submodel.add(epa1 + epa2 >= 0);
		//submodel.add(startF[k] <= nbb[k] - 1);
		//submodel.add(mF[k][0] == 0);
		epa1.end();
		epa2.end();
	}

	//for (k = 0; k < nbCrane; k++)
	//{
	//	if ((1 + safe_margin)*k + 1 > nbBay)
	//	{
	//		for (i = 0; i < nbTask; i++)
	//			model.add( yF[k][i] == 0);
	//	}
	//}


	//**********************************//
	//            开始求解		        //
	//**********************************//
	IloCplex cplex(env);
	cplex.extract(submodel);
	//cplex.exportModel("LP_format_of_CB.LP");
	//************************************************************************************//
	//      output LP format                                                         //
	//************************************************************************************//
	//char* filename1;
	//char dream1[100] = "LP_file/CB_LtoR";
	//filename1 = dream1;
	//char C1[3];
	//char C2[3];
	//char C3[3];
	//char C4[3];
	//itoa(my, C1, 10);

	////此处可编辑规模，以输入
	//itoa(nbBay, C2, 10);
	//itoa(nbCrane, C3, 10);
	//itoa(nbTask, C4, 10);

	//strcpy(filename1, "LP_file/");
	////此处可编辑规模，以输入
	//strcat(filename1, C4);
	//strcat(filename1, "-");
	//strcat(filename1, C2);
	//strcat(filename1, "-");
	//strcat(filename1, C3);

	//strcat(filename1, "/CB_LtoR");
	//strcat(filename1, "-");
	//strcat(filename1, C1);
	//strcat(filename1, ".LP");

	//cplex.exportModel(filename1);
#ifdef UserActive
	cplex.use(CandyUserCallback(env));
	//cplex.setParam(IloCplex::MIPEmphasis, 2);
	//cplex.setParam(IloCplex::HeurFreq, 1);
	//cplex.setParam(IloCplex::ParallelMode, 1);
	cplex.setParam(IloCplex::Threads, 4);
#endif



	//cplex.setParam(IloCplex::MIPEmphasis ,CPX_MIPEMPHASIS_HIDDENFEAS);
	//cplex.setParam(IloCplex::EpGap,0.001);
	//cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::TiLim, 1800);


	//cout<<endl<<"number of binary variables = "<<cplex.getNbinVars()<<endl;
	//cout<<endl<<"number of constraints = "<<cplex.getNrows()<<endl;



	BOOL h0;
	h0 = cplex.solve();
	//cout<<"h0 "<<h0<<endl;
	if (!h0)
	{
		cout << "\nno feasible solution has been found" << endl; return FALSE;
		cplex.clearModel();
		cplex.clear();
		cplex.end();
		submodel.end();

		return FALSE;
	}
	//**********************************//
	//             记录最好解           //
	//**********************************//

	//*GapF=cplex.getMIPRelativeGap();
	*ObjVal = cplex.getBestObjValue();

	*UB = cplex.getValue(CF);
	*CF_best = cplex.getValue(CF);
	//*UB=(int)(*UB)+1; 

	//cout<<"yF_best"<<endl;
	for (k = 0; k < nbCrane; k++)
	{
		//int sumQK = 0;
		for (i = 0; i < nbTask; i++)
		{
			if (cplex.isExtracted(yF[k][i]))
			{
				if (cplex.getValue(yF[k][i]) > 0.1)
				{
					yF_best[k][i] = 1;
					//cout << i + 1 << "  ";
				}

				else
					yF_best[k][i] = 0;
				//cout << yF_best[k][i] << "  ";

				//sumQK += nbQ[i] * yF_best[k][i];
			}
			else
				yF_best[k][i] = 0;
		}//cout << sumQK<<"    "<<endl;
	}
	//cout << endl << endl;

	//cout << "mF_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	for (i = 0; i < nbBay; i++)
	//	{
	//		cout << cplex.getValue(mF[k][i]) << "  ";
	//	}cout << "    "; cout << endl;
	//}
	//cout << endl << endl;

	//cout << "uF_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	for (i = 0; i < nbBay; i++)
	//	{
	//		cout << cplex.getValue(CuF[k][i]) << "  ";
	//	}cout << "    "; cout << endl;
	//}
	//cout << endl << endl;

	//cout << "vF_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	for (i = 0; i < nbBay; i++)
	//	{
	//		cout << cplex.getValue(CvF[k][i]) << "  ";
	//	}cout << "    "; cout << endl;
	//}
	//cout << endl << endl;

	//cout << "aF_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	cout << cplex.getValue(startF[k]) << "    "<<cplex.getValue(endF[k]); cout << endl;
	//}
	//cout << endl << endl;

	//	cout << "t0F_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	cout << k << "QC    " << cplex.getValue(t0kF[k]); cout << endl;
	//}
	//cout << endl << endl;

	//for (i = 0; i < nbBay; i++)
	//{
	//	cout << i+1 << ":  ";
	//	for (j = 0; j < nbTask; j++)
	//		if (nbLocation[j]==i+1)
	//		{
	//			cout <<j+1 << "  ";
	//		}
	//	cout << endl;
	//}cout << "    "; cout << endl;


	cplex.clearModel();
	cplex.clear();
	cplex.end();
	submodel.end();

	return TRUE;



}

BOOL LB_Chen_Same_Direction(IloModel model, IloNum nbs, IloIntArray  nbb, IloNumArray nbQ, IloIntArray  nbreadyT,
	BoolVarMatrix yF, IloIntVar CF, IloInt *CF_best, BoolMatrix yF_best,
	IloNum *ObjVal, IloNum *UB, IloIntArray  nbLocation, BoolMatrix nbprecR, int my)
{
	IloEnv env = model.getEnv();
	IloInt i, j, k;


	//问题模型
	IloModel submodel(env);
	BoolVarMatrix CuF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		CuF[k] = IloBoolVarArray(env, nbBay);
	}
	BoolVarMatrix CvF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		CvF[k] = IloBoolVarArray(env, nbBay);
	}
	NumVarMatrix mF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		mF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);

	IloIntVarArray startF(env, nbCrane, 0, 100);
	IloIntVarArray endF(env, nbCrane, 0, 100);

	IloNumVarArray t0kF(env, nbCrane, 0, 100);

	//leaving time at each bay
	NumVarMatrix TF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		TF[k] = IloNumVarArray(env, nbBay, -100, IloInfinity);

	//**********************************//
	//            原问题目标函数        //
	//**********************************//

	//IloExpr obj1(env); 


	IloExpr  obj2(env);
	//  建立子问题目标函数表达式 

	obj2 += CF;

	//	将目标函数加入到原问题模型
	submodel.add(IloMinimize(env, obj2));//
	//obj1.end();
	obj2.end();


	//**********************************//
	//            MP问题 约束           //
	//**********************************//

	//IloRangeArray  c0(env);
	//for (k = 0; k < nbCrane; k++)
	//{
	//	IloExpr  epa(env);
	//	for (j = 0; j < nbTask; j++)
	//		if (nbLocation[j] >(nbBay - 2 * (nbCrane - k - 1)))
	//			epa += yF[k][j];
	//	c0.add(epa == 0);
	//	epa.end();
	//}
	//submodel.add(c0);
	//c0.end();


	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbBay; i++)
		{
		IloExpr  epa(env);
		epa += nbreadyT[k] + t0kF[k] - QCmove_time*startF[k];
		for (j = 0; j < nbTask; j++)
			if (nbLocation[j] - 1 <= i)
				epa += (nbQ[j] / nbs)*yF[k][j];
		for (j = 0; j <= i; j++)
		{

			epa += mF[k][j];
		}
		submodel.add(epa + QCmove_time*i - TF[k][i] == 0);
		epa.end();

		}

	IloRangeArray  c2(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		epa += CF - t0kF[k] - nbreadyT[k];
		for (i = 0; i < nbTask; i++)
			epa -= (nbQ[i] / nbs)*yF[k][i];
		for (i = 0; i < nbBay; i++)
			epa -= mF[k][i];


		epa += QCmove_time*startF[k];
		epa -= QCmove_time*endF[k];

		c2.add(epa >= 0);
		epa.end();
	}
	submodel.add(c2);
	c2.end();

	//约束（3）
	IloRangeArray  c3(env);
	for (i = 0; i < nbBay - 2; i++)
		for (k = 0; k < nbCrane - 1; k++)
		{
		IloExpr  epa(env);
		epa += TF[k][i] - TF[k + 1][i + 2];

		epa += 2000 * (CvF[k][i] + CuF[k + 1][i + 2]);
		//epa += 2000 * (CuF[k][i] + CuF[k + 1][i + 2]);

		c3.add(epa >= 0);
		epa.end();
		}
	submodel.add(c3);
	c3.end();

	//约束（4）
	IloRangeArray  c4(env);
	for (i = 0; i < nbTask; i++)
	{
		IloExpr  epa(env);
		for (k = 0; k < nbCrane; k++)
			epa += yF[k][i];
		c4.add(epa == 1);
		epa.end();
	}
	submodel.add(c4);
	c4.end();

	////约束（40）
	//IloRangeArray  c40(env);
	//for (k = 0; k < nbCrane; k++)
	//{
	//	IloExpr  epa(env);
	//	for (i = 0; i < nbTask; i++)
	//		epa += yF[k][i];
	//	c40.add(epa >= 1);
	//	epa.end();
	//}
	//submodel.add(c40);
	//c40.end();

	//约束（5）
	IloRangeArray  c5(env);
	for (i = 0; i < nbTask; i++)
		for (k = 0; k < nbCrane; k++)
		{
		IloExpr  epa(env);
		epa += (nbLocation[i] - 1) * yF[k][i];
		epa -= startF[k];
		epa += 2 * nbBay * (1 - yF[k][i]);
		c5.add(epa >= 0);
		epa.end();
		}
	submodel.add(c5);
	c5.end();

	//约束（6）
	IloRangeArray  c6(env);
	for (i = 0; i < nbTask; i++)
		for (k = 0; k < nbCrane; k++)
		{
		IloExpr  epa(env);
		epa -= (nbLocation[i] - 1) * yF[k][i];
		epa += endF[k];
		c6.add(epa >= 0);
		epa.end();
		}
	submodel.add(c6);
	c6.end();



	////	建立约束(9)
	//IloRangeArray  c9(env);
	//for (i = 0; i < nbTask; i++)
	//	for (j = 0; j < nbTask; j++)
	//		if (nbprecR[i][j] == 1)
	//		{
	//	for (int l = 0; l < nbCrane - 1; l++)
	//		for (k = l + 1; k < nbCrane; k++)
	//		{
	//		IloExpr  epa(env);
	//		epa += yF[l][i] + yF[k][j];
	//		c9.add(epa <= 1);
	//		}

	//		}
	//submodel.add(c9);
	//c9.end();


	//约束（5）
	IloRangeArray  c10a(env);
	for (i = 0; i < nbBay; i++)
		for (k = 0; k < nbCrane; k++)
		{
		IloExpr  epa1(env);
		IloExpr  epa2(env);
		IloExpr  epa3(env);
		epa1 += i - startF[k];
		epa2 += 2 * nbBay * (1 - CuF[k][i]);
		epa3 += 2 * nbBay * CuF[k][i];
		c10a.add(epa1 - epa2 + 0.1 <= 0);
		c10a.add(epa3 + epa1 >= 0);
		epa1.end();
		epa2.end();
		epa3.end();
		}
	submodel.add(c10a);
	c10a.end();


	//约束（5）
	IloRangeArray  c10a2(env);
	for (i = 0; i < nbBay; i++)
		for (k = 0; k < nbCrane; k++)
		{
		IloExpr  epa1(env);
		IloExpr  epa2(env);
		IloExpr  epa3(env);
		epa1 += endF[k] - i;
		epa2 += 2 * nbBay * (1 - CvF[k][i]);
		epa3 += 2 * nbBay * CvF[k][i];
		c10a2.add(epa1 - epa2 + 0.1 <= 0);
		c10a2.add(epa3 + epa1 >= 0);
		epa1.end();
		epa2.end();
		epa3.end();
		}
	submodel.add(c10a2);
	c10a2.end();

	for (k = 0; k < nbCrane - 1; k++)
	{
		IloExpr  epa1(env);
		IloExpr  epa2(env);
		epa1 += startF[k + 1] - startF[k];
		epa2 += endF[k + 1] - endF[k];
		submodel.add(epa1 >= 1 + safe_margin);
		submodel.add(epa2 >= 1 + safe_margin);
		epa1.end();
		epa2.end();
	}

	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa1(env);
		IloExpr  epa2(env);
		epa1 += t0kF[k];
		epa2 += QCmove_time*startF[k] - QCmove_time*nbb[k] + QCmove_time;
		submodel.add(epa1 - epa2 >= 0);
		submodel.add(epa1 + epa2 >= 0);
		//submodel.add(startF[k] <= nbb[k] - 1);
		//submodel.add(mF[k][0] == 0);
		epa1.end();
		epa2.end();
	}

	//for (k = 0; k < nbCrane; k++)
	//{
	//	if ((1 + safe_margin)*k + 1 > nbBay)
	//	{
	//		for (i = 0; i < nbTask; i++)
	//			model.add( yF[k][i] == 0);
	//	}
	//}


	//**********************************//
	//            开始求解		        //
	//**********************************//
	IloCplex cplex(env);
	cplex.extract(submodel);
	//cplex.exportModel("LP_format_of_CB.LP");
	//************************************************************************************//
	//      output LP format                                                         //
	//************************************************************************************//
	//char* filename1;
	//char dream1[100] = "LP_file/CB_LtoR";
	//filename1 = dream1;
	//char C1[3];
	//char C2[3];
	//char C3[3];
	//char C4[3];
	//itoa(my, C1, 10);

	////此处可编辑规模，以输入
	//itoa(nbBay, C2, 10);
	//itoa(nbCrane, C3, 10);
	//itoa(nbTask, C4, 10);

	//strcpy(filename1, "LP_file/");
	////此处可编辑规模，以输入
	//strcat(filename1, C4);
	//strcat(filename1, "-");
	//strcat(filename1, C2);
	//strcat(filename1, "-");
	//strcat(filename1, C3);

	//strcat(filename1, "/CB_LtoR");
	//strcat(filename1, "-");
	//strcat(filename1, C1);
	//strcat(filename1, ".LP");

	//cplex.exportModel(filename1);
#ifdef UserActive
	cplex.use(CandyUserCallback(env));
	//cplex.setParam(IloCplex::MIPEmphasis, 2);
	//cplex.setParam(IloCplex::HeurFreq, 1);
	//cplex.setParam(IloCplex::ParallelMode, 1);
	cplex.setParam(IloCplex::Threads, 4);
#endif



	//cplex.setParam(IloCplex::MIPEmphasis ,CPX_MIPEMPHASIS_HIDDENFEAS);
	//cplex.setParam(IloCplex::EpGap,0.001);
	//cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::TiLim, 1800);


	//cout<<endl<<"number of binary variables = "<<cplex.getNbinVars()<<endl;
	//cout<<endl<<"number of constraints = "<<cplex.getNrows()<<endl;



	BOOL h0;
	h0 = cplex.solve();
	//cout<<"h0 "<<h0<<endl;
	if (!h0)
	{
		cout << "\nno feasible solution has been found" << endl; return FALSE;
		cplex.clearModel();
		cplex.clear();
		cplex.end();
		submodel.end();

		return FALSE;
	}
	//**********************************//
	//             记录最好解           //
	//**********************************//

	//*GapF=cplex.getMIPRelativeGap();
	*ObjVal = cplex.getBestObjValue();

	*UB = cplex.getValue(CF);
	*CF_best = cplex.getValue(CF);
	//*UB=(int)(*UB)+1; 

	//cout<<"yF_best"<<endl;
	for (k = 0; k < nbCrane; k++)
	{
		//int sumQK = 0;
		for (i = 0; i < nbTask; i++)
		{
			if (cplex.isExtracted(yF[k][i]))
			{
				if (cplex.getValue(yF[k][i]) > 0.1)
				{
					yF_best[k][i] = 1;
					//cout << i + 1 << "  ";
				}

				else
					yF_best[k][i] = 0;
				//cout << yF_best[k][i] << "  ";

				//sumQK += nbQ[i] * yF_best[k][i];
			}
			else
				yF_best[k][i] = 0;
		}//cout << sumQK<<"    "<<endl;
	}
	//cout << endl << endl;

	//cout << "mF_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	for (i = 0; i < nbBay; i++)
	//	{
	//		cout << cplex.getValue(mF[k][i]) << "  ";
	//	}cout << "    "; cout << endl;
	//}
	//cout << endl << endl;

	//cout << "uF_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	for (i = 0; i < nbBay; i++)
	//	{
	//		cout << cplex.getValue(CuF[k][i]) << "  ";
	//	}cout << "    "; cout << endl;
	//}
	//cout << endl << endl;

	//cout << "vF_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	for (i = 0; i < nbBay; i++)
	//	{
	//		cout << cplex.getValue(CvF[k][i]) << "  ";
	//	}cout << "    "; cout << endl;
	//}
	//cout << endl << endl;

	//cout << "aF_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	cout << cplex.getValue(startF[k]) << "    "<<cplex.getValue(endF[k]); cout << endl;
	//}
	//cout << endl << endl;

	//	cout << "t0F_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	cout << k << "QC    " << cplex.getValue(t0kF[k]); cout << endl;
	//}
	//cout << endl << endl;

	//for (i = 0; i < nbBay; i++)
	//{
	//	cout << i+1 << ":  ";
	//	for (j = 0; j < nbTask; j++)
	//		if (nbLocation[j]==i+1)
	//		{
	//			cout <<j+1 << "  ";
	//		}
	//	cout << endl;
	//}cout << "    "; cout << endl;


	cplex.clearModel();
	cplex.clear();
	cplex.end();
	submodel.end();

	return TRUE;



}

BOOL Same_Direction(IloModel model, IloNum nbs, IloIntArray  nbb, IloNumArray nbQ, IloIntArray  nbreadyT,
	BoolVarMatrix yF, IloIntVar CF, IloInt *CF_best, BoolMatrix yF_best,
	IloNum *ObjVal, IloNum *UB, IloIntArray  nbLocation, BoolMatrix nbprecR, int my)
{
	IloEnv env = model.getEnv();
	IloInt i, j, k;

	//问题模型
	IloModel submodel(env);
	IloNumVarArray t0kF(env, nbCrane, 0, 100);
	//IloNumVarArray betaF(env, nbCrane, 0, 100);

	NumVarMatrix mF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		mF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);

	BoolVarMatrix CzF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		CzF[k] = IloBoolVarArray(env, nbBay);
	}

	BoolVarMatrix CzF2(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		CzF2[k] = IloBoolVarArray(env, nbBay);
	}

	//leaving time at each bay
	NumVarMatrix TF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		TF[k] = IloNumVarArray(env, nbBay, -100, IloInfinity);


	//**********************************//
	//            原问题目标函数        //
	//**********************************//

	//IloExpr obj1(env); 


	IloExpr  obj2(env);
	//  建立子问题目标函数表达式 

	obj2 += CF;

	//	将目标函数加入到原问题模型
	submodel.add(IloMinimize(env, obj2));//
	//obj1.end();
	obj2.end();


	//**********************************//
	//            MP问题 约束           //
	//**********************************//




	IloRangeArray  c2(env);
	for (k = 0; k < nbCrane; k++)
		//for (i = 0; i < nbTask; i++)
	{
		IloExpr  epa(env);
		epa += CF - t0kF[k] - nbreadyT[k];
		for (int ii = 0; ii < nbTask; ii++)
			epa -= (nbQ[ii] / nbs)*yF[k][ii];
		//for (j = 0; j < nbLocation[i]-1; j++)	
		for (j = 0; j < nbBay; j++)
			epa -= mF[k][j];
		//epa -= (nbLocation[i] - 1)*yF[k][i];

		//epa -= betaF[k];

		for (j = 0; j < nbBay; j++)//for (j = 2 * k; j < nbBay - 2 * (nbCrane - k - 1); j++)
			epa -= QCmove_time*j*CzF2[k][j];

		//for (j = 0; j < nbLocation[i]; j++)
		for (j = 0; j < nbBay; j++)//for (j = 2 * k; j < nbBay - 2 * (nbCrane - k - 1); j++)
			epa += QCmove_time*j*CzF[k][j];

		c2.add(epa >= 0);
		epa.end();
	}
	submodel.add(c2);
	c2.end();

	//约束（4）
	IloRangeArray  c4(env);
	for (i = 0; i < nbTask; i++)
	{
		IloExpr  epa(env);
		for (k = 0; k < nbCrane; k++)
			epa += yF[k][i];
		c4.add(epa == 1);
		epa.end();
	}
	submodel.add(c4);
	c4.end();

	for (k = 0; k < nbCrane; k++)
		for (j = 0; j < nbTask; j++)
			if (nbLocation[j] < nbBay)
			{
		IloExpr  epa(env);

		for (i = nbLocation[j]; i < nbBay; i++)
			epa += CzF[k][i];
		//epa += betaF[k];
		epa += yF[k][j];

		submodel.add(epa <= 1);
		epa.end();
			}


	for (k = 0; k < nbCrane; k++)
		for (j = 0; j < nbTask; j++)
		{
		//model.add(thetaF[k] - (nbLocation[j] - 1) * xF[k][j] >= 0);
		IloExpr  epa(env);
		for (i = 0; i < nbBay; i++)
			epa += i*CzF2[k][i];
		submodel.add(epa - (nbLocation[j] - 1) *  yF[k][j] >= 0);
		epa.end();
		}

	//约束（7）
	IloRangeArray  c7(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbBay; i++)
			epa += CzF[k][i];
		c7.add(epa == 1);
		epa.end();

		IloExpr  epa2(env);
		for (i = 0; i < nbBay; i++)
			epa2 += CzF2[k][i];
		c7.add(epa2 == 1);
		epa2.end();
	}
	submodel.add(c7);
	c7.end();

	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa1(env);
		IloExpr  epa2(env);
		epa1 += t0kF[k];
		for (i = 0; i < nbBay; i++)
			epa2 += QCmove_time*i*CzF[k][i];
		epa2 -= QCmove_time*nbb[k] - QCmove_time;
		submodel.add(epa1 - epa2 >= 0);
		submodel.add(epa1 + epa2 >= 0);
		epa1.end();
		epa2.end();
	}

	//13o
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbBay; i++)
			epa += i*CzF[k][i] - i*CzF2[k][i];
		submodel.add(epa <= 0);
		epa.end();
	}


	//	建立约束(8)
	IloRangeArray  c8(env);
	for (k = 0; k < nbCrane - 1; k++)
	{
		if ((1 + safe_margin)*(k + 1) + 1 <= nbBay)
		{
			IloExpr  epa(env);
			for (i = 0; i < nbBay; i++)
				epa += i*CzF[k + 1][i] - i*CzF[k][i];
			c8.add(epa >= 1 + safe_margin);
			epa.end();
		}

		if ((1 + safe_margin)*(k + 2) <= nbBay)
		{
			IloExpr  epa2(env);
			for (i = 0; i < nbBay; i++)
				epa2 += i*CzF2[k + 1][i] - i*CzF2[k][i];
			c8.add(epa2 >= 1 + safe_margin);
			epa2.end();
		}
		//IloExpr  epa3(env);
		//for (i = 0; i < nbBay; i++)
		//	epa3 += i*CzF2[k][i] - i*CzF[k][i];
		//c8.add(epa3 >= 0);
		//epa3.end();

	}
	submodel.add(c8);
	c8.end();

	//for (k = 0; k < nbCrane; k++)
	//{
	//	if ((1 + safe_margin)*k + 1 > nbBay)
	//	{
	//		for (i = 0; i < nbTask; i++)
	//			submodel.add(yF[k][i] == 0);
	//	}
	//}

	IloRangeArray  c010(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		for (j = 0; j < nbTask; j++)
			if (nbLocation[j] >(nbBay - 2 * (nbCrane - k - 1)))
				epa += yF[k][j];
		c010.add(epa == 0);
		epa.end();
	}
	model.add(c010);
	c010.end();


	//	建立约束(9)
	IloRangeArray  c9(env);
	for (i = 0; i < nbTask; i++)
		for (j = 0; j < nbTask; j++)
			if (nbprecR[i][j] == 1)
			{
		for (int l = 0; l < nbCrane - 1; l++)
			for (k = l + 1; k < nbCrane; k++)
			{
			IloExpr  epa(env);
			epa += yF[l][i] + yF[k][j];
			c9.add(epa <= 1);
			}

			}
	submodel.add(c9);
	c9.end();



	//定义TF
	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbBay; i++)//for (i = (1 + safe_margin) * k; i < nbBay - (1 + safe_margin) * (nbCrane - k - 1); i++)
		{
		IloExpr  epa(env);
		epa += t0kF[k] + nbreadyT[k];

		for (j = 0; j < nbTask; j++)
			if (nbLocation[j] - 1 <= i)//if (nbLocation[j] - 1 <= i &&nbLocation[j] - 1 >= 2 * k)
				epa += (nbQ[j] / nbs)*yF[k][j];

		for (j = 0; j <= i; j++)//for (j = 2 * k; j <= i; j++)
		{
			epa += mF[k][j];
			epa -= QCmove_time*j * CzF[k][j];
		}

		submodel.add(epa + QCmove_time*i - TF[k][i] == 0);
		epa.end();

		}

	//约束（3）
	IloRangeArray  c3(env);
	for (k = 0; k < nbCrane - 1; k++)
		for (i = 0; i < nbBay - 1 - safe_margin; i++)//for (i = 2 * k; i < nbBay - 2 * (nbCrane - k - 1); i++)
		{
		IloExpr  epa(env);

		epa += TF[k][i] - TF[k + 1][i + 1 + safe_margin];

		for (j = 0; j <= i + 1 + safe_margin; j++)
		{
			epa -= 2000 * CzF[k + 1][j];
		}
		for (j = 0; j < i; j++)
		{
			epa += 2000 * CzF2[k][j];
		}

		epa += 2000;

		c3.add(epa >= 0);
		epa.end();
		}
	submodel.add(c3);
	c3.end();


	//**********************************//
	//            开始求解		        //
	//**********************************//
	IloCplex cplex(env);
	cplex.extract(submodel);

	//cplex.exportModel("LP_F.LP");

	//************************************************************************************//
	//      output LP format                                                         //
	//************************************************************************************//
	//char* filename1;
	//char dream1[100] = "LP_file/F_LtoR";
	//filename1 = dream1;
	//char C1[3];
	//char C2[3];
	//char C3[3];
	//char C4[3];
	//itoa(my, C1, 10);

	////此处可编辑规模，以输入
	//itoa(nbBay, C2, 10);
	//itoa(nbCrane, C3, 10);
	//itoa(nbTask, C4, 10);

	//strcpy(filename1, "LP_file/");
	////此处可编辑规模，以输入
	//strcat(filename1, C4);
	//strcat(filename1, "-");
	//strcat(filename1, C2);
	//strcat(filename1, "-");
	//strcat(filename1, C3);

	//strcat(filename1, "/F_LtoR");
	//strcat(filename1, "-");
	//strcat(filename1, C1);
	//strcat(filename1, ".LP");

	//cplex.exportModel(filename1);

#ifdef UserActive
	cplex.use(CandyUserCallback(env));
	//cplex.setParam(IloCplex::MIPEmphasis, 2);
	//cplex.setParam(IloCplex::HeurFreq, 1);
	//cplex.setParam(IloCplex::ParallelMode, 1);
	cplex.setParam(IloCplex::Threads, 4);
#endif



	//cplex.setParam(IloCplex::MIPEmphasis ,CPX_MIPEMPHASIS_HIDDENFEAS);
	//cplex.setParam(IloCplex::EpGap,0.001);
	//cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::TiLim, 1800);


	//cout<<endl<<"number of binary variables = "<<cplex.getNbinVars()<<endl;
	//cout<<endl<<"number of constraints = "<<cplex.getNrows()<<endl;



	BOOL h0;
	h0 = cplex.solve();
	//cout<<"h0 "<<h0<<endl;
	if (!h0)
	{
		cout << "\nno feasible solution has been found" << endl; return FALSE;
		cplex.clearModel();
		cplex.clear();
		cplex.end();
		submodel.end();

		return FALSE;
	}
	//**********************************//
	//             记录最好解           //
	//**********************************//

	//*GapF=cplex.getMIPRelativeGap();
	*ObjVal = cplex.getBestObjValue();

	*UB = cplex.getValue(CF);
	*CF_best = cplex.getValue(CF);
	//*UB=(int)(*UB)+1; 


	//cout<<"yF_best"<<endl;
	for (k = 0; k < nbCrane; k++)
	{
		//int sumQK = 0;
		for (i = 0; i < nbTask; i++)
		{
			if (cplex.isExtracted(yF[k][i]))
			{
				if (cplex.getValue(yF[k][i]) > 0.1)
				{
					yF_best[k][i] = 1;
					//cout << i + 1 << "  ";
				}

				else
					yF_best[k][i] = 0;
				//cout << yF_best[k][i] << "  ";

				//sumQK += nbQ[i] * yF_best[k][i];
			}
			else
				yF_best[k][i] = 0;
		}//cout << sumQK<<"    "<<endl;
	}
	//cout << endl << endl;

	//cout << "mF_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	for (i = 0; i < nbBay; i++)
	//	{
	//		cout << cplex.getValue(mF[k][i]) << "  ";
	//	}cout << "    "; cout << endl;
	//}
	//cout << endl << endl;

	//cout << "aF_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	for (i = 0; i < nbBay; i++)
	//		if (cplex.getValue(CzF[k][i])>0.1)
	//			cout << i << "    " << k; 
	//	cout << endl;
	//}
	//cout << endl << endl;

	//cout << "betaF_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	for (i = 0; i < nbBay; i++)
	//		if (cplex.getValue(CzF2[k][i])>0.1)
	//			cout << i << "    " << k; 
	//	cout << endl;
	//}
	//cout << endl << endl;

	//cout << "t0F_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	cout << k << "QC    " << cplex.getValue(t0kF[k]); cout << endl;
	//}
	//cout << endl << endl;

	cplex.clearModel();
	cplex.clear();
	cplex.end();
	submodel.end();

	return TRUE;




}

BOOL LB_Same_Direction(IloModel model, IloNum nbs, IloIntArray  nbb, IloNumArray nbQ, IloIntArray  nbreadyT,
	BoolVarMatrix yF, IloIntVar CF, IloInt *CF_best, BoolMatrix yF_best,
	IloNum *ObjVal, IloNum *UB, IloIntArray  nbLocation, BoolMatrix nbprecR, int my)
{
	IloEnv env = model.getEnv();
	IloInt i, j, k;

	//问题模型
	IloModel submodel(env);
	IloNumVarArray t0kF(env, nbCrane, 0, 100);
	//IloNumVarArray betaF(env, nbCrane, 0, 100);

	NumVarMatrix mF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		mF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);

	BoolVarMatrix CzF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		CzF[k] = IloBoolVarArray(env, nbBay);
	}

	BoolVarMatrix CzF2(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		CzF2[k] = IloBoolVarArray(env, nbBay);
	}

	//leaving time at each bay
	NumVarMatrix TF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		TF[k] = IloNumVarArray(env, nbBay, -100, IloInfinity);


	//**********************************//
	//            原问题目标函数        //
	//**********************************//

	//IloExpr obj1(env); 


	IloExpr  obj2(env);
	//  建立子问题目标函数表达式 

	obj2 += CF;

	//	将目标函数加入到原问题模型
	submodel.add(IloMinimize(env, obj2));//
	//obj1.end();
	obj2.end();


	//**********************************//
	//            MP问题 约束           //
	//**********************************//




	IloRangeArray  c2(env);
	for (k = 0; k < nbCrane; k++)
		//for (i = 0; i < nbTask; i++)
	{
		IloExpr  epa(env);
		epa += CF - t0kF[k] - nbreadyT[k];
		for (int ii = 0; ii < nbTask; ii++)
			epa -= (nbQ[ii] / nbs)*yF[k][ii];
		//for (j = 0; j < nbLocation[i]-1; j++)	
		for (j = 0; j < nbBay; j++)
			epa -= mF[k][j];
		//epa -= (nbLocation[i] - 1)*yF[k][i];

		//epa -= betaF[k];

		for (j = 0; j < nbBay; j++)//for (j = 2 * k; j < nbBay - 2 * (nbCrane - k - 1); j++)
			epa -= QCmove_time*j*CzF2[k][j];

		//for (j = 0; j < nbLocation[i]; j++)
		for (j = 0; j < nbBay; j++)//for (j = 2 * k; j < nbBay - 2 * (nbCrane - k - 1); j++)
			epa += QCmove_time*j*CzF[k][j];

		c2.add(epa >= 0);
		epa.end();
	}
	submodel.add(c2);
	c2.end();

	//约束（4）
	IloRangeArray  c4(env);
	for (i = 0; i < nbTask; i++)
	{
		IloExpr  epa(env);
		for (k = 0; k < nbCrane; k++)
			epa += yF[k][i];
		c4.add(epa == 1);
		epa.end();
	}
	submodel.add(c4);
	c4.end();

	for (k = 0; k < nbCrane; k++)
		for (j = 0; j < nbTask; j++)
			if (nbLocation[j] < nbBay)
			{
		IloExpr  epa(env);

		for (i = nbLocation[j]; i < nbBay; i++)
			epa += CzF[k][i];
		//epa += betaF[k];
		epa += yF[k][j];

		submodel.add(epa <= 1);
		epa.end();
			}


	for (k = 0; k < nbCrane; k++)
		for (j = 0; j < nbTask; j++)
		{
		//model.add(thetaF[k] - (nbLocation[j] - 1) * xF[k][j] >= 0);
		IloExpr  epa(env);
		for (i = 0; i < nbBay; i++)
			epa += i*CzF2[k][i];
		submodel.add(epa - (nbLocation[j] - 1) *  yF[k][j] >= 0);
		epa.end();
		}

	//约束（7）
	IloRangeArray  c7(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbBay; i++)
			epa += CzF[k][i];
		c7.add(epa == 1);
		epa.end();

		IloExpr  epa2(env);
		for (i = 0; i < nbBay; i++)
			epa2 += CzF2[k][i];
		c7.add(epa2 == 1);
		epa2.end();
	}
	submodel.add(c7);
	c7.end();

	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa1(env);
		IloExpr  epa2(env);
		epa1 += t0kF[k];
		for (i = 0; i < nbBay; i++)
			epa2 += QCmove_time*i*CzF[k][i];
		epa2 -= QCmove_time*nbb[k] - QCmove_time;
		submodel.add(epa1 - epa2 >= 0);
		submodel.add(epa1 + epa2 >= 0);
		epa1.end();
		epa2.end();
	}

	//13o
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbBay; i++)
			epa += i*CzF[k][i] - i*CzF2[k][i];
		submodel.add(epa <= 0);
		epa.end();
	}


	//	建立约束(8)
	IloRangeArray  c8(env);
	for (k = 0; k < nbCrane - 1; k++)
	{
		if ((1 + safe_margin)*(k + 1) + 1 <= nbBay)
		{
			IloExpr  epa(env);
			for (i = 0; i < nbBay; i++)
				epa += i*CzF[k + 1][i] - i*CzF[k][i];
			c8.add(epa >= 1 + safe_margin);
			epa.end();
		}

		if ((1 + safe_margin)*(k + 2) <= nbBay)
		{
			IloExpr  epa2(env);
			for (i = 0; i < nbBay; i++)
				epa2 += i*CzF2[k + 1][i] - i*CzF2[k][i];
			c8.add(epa2 >= 1 + safe_margin);
			epa2.end();
		}
		//IloExpr  epa3(env);
		//for (i = 0; i < nbBay; i++)
		//	epa3 += i*CzF2[k][i] - i*CzF[k][i];
		//c8.add(epa3 >= 0);
		//epa3.end();

	}
	submodel.add(c8);
	c8.end();

	//for (k = 0; k < nbCrane; k++)
	//{
	//	if ((1 + safe_margin)*k + 1 > nbBay)
	//	{
	//		for (i = 0; i < nbTask; i++)
	//			submodel.add(yF[k][i] == 0);
	//	}
	//}

	IloRangeArray  c010(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		for (j = 0; j < nbTask; j++)
			if (nbLocation[j] >(nbBay - 2 * (nbCrane - k - 1)))
				epa += yF[k][j];
		c010.add(epa == 0);
		epa.end();
	}
	model.add(c010);
	c010.end();


	////	建立约束(9)
	//IloRangeArray  c9(env);
	//for (i = 0; i < nbTask; i++)
	//	for (j = 0; j < nbTask; j++)
	//		if (nbprecR[i][j] == 1)
	//		{
	//	for (int l = 0; l < nbCrane - 1; l++)
	//		for (k = l + 1; k < nbCrane; k++)
	//		{
	//		IloExpr  epa(env);
	//		epa += yF[l][i] + yF[k][j];
	//		c9.add(epa <= 1);
	//		}

	//		}
	//submodel.add(c9);
	//c9.end();




	//定义TF
	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbBay; i++)//for (i = (1 + safe_margin) * k; i < nbBay - (1 + safe_margin) * (nbCrane - k - 1); i++)
		{
		IloExpr  epa(env);
		epa += t0kF[k] + nbreadyT[k];

		for (j = 0; j < nbTask; j++)
			if (nbLocation[j] - 1 <= i)//if (nbLocation[j] - 1 <= i &&nbLocation[j] - 1 >= 2 * k)
				epa += (nbQ[j] / nbs)*yF[k][j];

		for (j = 0; j <= i; j++)//for (j = 2 * k; j <= i; j++)
		{
			epa += mF[k][j];
			epa -= QCmove_time*j * CzF[k][j];
		}

		submodel.add(epa + QCmove_time*i - TF[k][i] == 0);
		epa.end();

		}

	//约束（3）
	IloRangeArray  c3(env);
	for (k = 0; k < nbCrane - 1; k++)
		for (i = 0; i < nbBay - 1 - safe_margin; i++)//for (i = 2 * k; i < nbBay - 2 * (nbCrane - k - 1); i++)
		{
		IloExpr  epa(env);

		epa += TF[k][i] - TF[k + 1][i + 1 + safe_margin];

		for (j = 0; j <= i + 1 + safe_margin; j++)
		{
			epa -= 2000 * CzF[k + 1][j];
		}
		for (j = 0; j < i; j++)
		{
			epa += 2000 * CzF2[k][j];
		}

		epa += 2000;

		c3.add(epa >= 0);
		epa.end();
		}
	submodel.add(c3);
	c3.end();


	//**********************************//
	//            开始求解		        //
	//**********************************//
	IloCplex cplex(env);
	cplex.extract(submodel);

	//cplex.exportModel("LP_F.LP");

	//************************************************************************************//
	//      output LP format                                                         //
	//************************************************************************************//
	//char* filename1;
	//char dream1[100] = "LP_file/F_LtoR";
	//filename1 = dream1;
	//char C1[3];
	//char C2[3];
	//char C3[3];
	//char C4[3];
	//itoa(my, C1, 10);

	////此处可编辑规模，以输入
	//itoa(nbBay, C2, 10);
	//itoa(nbCrane, C3, 10);
	//itoa(nbTask, C4, 10);

	//strcpy(filename1, "LP_file/");
	////此处可编辑规模，以输入
	//strcat(filename1, C4);
	//strcat(filename1, "-");
	//strcat(filename1, C2);
	//strcat(filename1, "-");
	//strcat(filename1, C3);

	//strcat(filename1, "/F_LtoR");
	//strcat(filename1, "-");
	//strcat(filename1, C1);
	//strcat(filename1, ".LP");

	//cplex.exportModel(filename1);

#ifdef UserActive
	cplex.use(CandyUserCallback(env));
	//cplex.setParam(IloCplex::MIPEmphasis, 2);
	//cplex.setParam(IloCplex::HeurFreq, 1);
	//cplex.setParam(IloCplex::ParallelMode, 1);
	cplex.setParam(IloCplex::Threads, 4);
#endif



	//cplex.setParam(IloCplex::MIPEmphasis ,CPX_MIPEMPHASIS_HIDDENFEAS);
	//cplex.setParam(IloCplex::EpGap,0.001);
	//cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::TiLim, 1800);


	//cout<<endl<<"number of binary variables = "<<cplex.getNbinVars()<<endl;
	//cout<<endl<<"number of constraints = "<<cplex.getNrows()<<endl;



	BOOL h0;
	h0 = cplex.solve();
	//cout<<"h0 "<<h0<<endl;
	if (!h0)
	{
		cout << "\nno feasible solution has been found" << endl; return FALSE;
		cplex.clearModel();
		cplex.clear();
		cplex.end();
		submodel.end();

		return FALSE;
	}
	//**********************************//
	//             记录最好解           //
	//**********************************//

	//*GapF=cplex.getMIPRelativeGap();
	*ObjVal = cplex.getBestObjValue();

	*UB = cplex.getValue(CF);
	*CF_best = cplex.getValue(CF);
	//*UB=(int)(*UB)+1; 


	//cout<<"yF_best"<<endl;
	for (k = 0; k < nbCrane; k++)
	{
		//int sumQK = 0;
		for (i = 0; i < nbTask; i++)
		{
			if (cplex.isExtracted(yF[k][i]))
			{
				if (cplex.getValue(yF[k][i]) > 0.1)
				{
					yF_best[k][i] = 1;
					//cout << i + 1 << "  ";
				}

				else
					yF_best[k][i] = 0;
				//cout << yF_best[k][i] << "  ";

				//sumQK += nbQ[i] * yF_best[k][i];
			}
			else
				yF_best[k][i] = 0;
		}//cout << sumQK<<"    "<<endl;
	}
	//cout << endl << endl;

	//cout << "mF_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	for (i = 0; i < nbBay; i++)
	//	{
	//		cout << cplex.getValue(mF[k][i]) << "  ";
	//	}cout << "    "; cout << endl;
	//}
	//cout << endl << endl;

	//cout << "aF_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	for (i = 0; i < nbBay; i++)
	//		if (cplex.getValue(CzF[k][i])>0.1)
	//			cout << i << "    " << k; 
	//	cout << endl;
	//}
	//cout << endl << endl;

	//cout << "betaF_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	for (i = 0; i < nbBay; i++)
	//		if (cplex.getValue(CzF2[k][i])>0.1)
	//			cout << i << "    " << k; 
	//	cout << endl;
	//}
	//cout << endl << endl;

	//cout << "t0F_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	cout << k << "QC    " << cplex.getValue(t0kF[k]); cout << endl;
	//}
	//cout << endl << endl;

	cplex.clearModel();
	cplex.clear();
	cplex.end();
	submodel.end();

	return TRUE;




}




BOOL Split_SP_1(IloModel model, IloNum nbs, IloIntArray  nbb, IloNumArray nbQ, IloIntArray  nbLocation, BoolMatrix nbprecR, IloIntArray  nbreadyT,
	BoolVarMatrix xF, BoolVarMatrix yF, BoolVarMatrix zF, IloIntVar CF, IloNumVarArray QC_CF, IloIntVarArray thetaF, IloBoolVarArray vF,
	BoolMatrix yF_best, BoolMatrix CxF_best, IloNum UB, IloNum *ObjVal, NumMatrix wF_best, NumMatrix CwF_best, int QC)
{

	//xF: forward trip, yF: retracing trip


	IloEnv env = model.getEnv();
	IloInt i, j, k;

	//问题模型
	IloModel submodel(env);

	//问题模型
	IloNumVarArray t0kF(env, nbCrane, 0, 100);
	IloNumVarArray gammaF(env, nbCrane, 0, 100);//

	IloNumVarArray ckF(env, nbCrane, 0, 3000);//
	IloNumVarArray chaF(env, nbCrane, 0, 3000);//

	NumVarMatrix wF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		wF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);

	NumVarMatrix CwF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		CwF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);

	BoolVarMatrix CzF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		CzF[k] = IloBoolVarArray(env, nbBay);
	}

	BoolVarMatrix endCzF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		endCzF[k] = IloBoolVarArray(env, nbBay);
	}

	//leaving time at each bay
	NumVarMatrix TF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		TF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);
	NumVarMatrix CTF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		CTF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);


	//**********************************//
	//           objective 原问题目标函数        //
	//**********************************//

	//IloExpr obj1(env); 


	IloExpr  obj2(env);
	//  建立子问题目标函数表达式 

	obj2 += CF;

	//for (k = 0; k < nbCrane; k++)
	//	obj2 += chaF[k];

	//	将目标函数加入到原问题模型
	submodel.add(IloMinimize(env, obj2));//
	//obj1.end();
	obj2.end();


	//**********************************//
	//            variable fixing           //
	//**********************************//
	//for (i = 0; i < nbTask - 1; i++)
	//	for (j = i + 1; j < nbTask; j++)
	//		if (nbprecR[i][j] == 1)
	//		{
	//	int i_k, j_k;
	//	i_k = -1;
	//	j_k = -1;
	//	for (k = QC; k <= QC + 1; k++)
	//	{
	//		if (yF_best[k][i] == 1)
	//			i_k = k;
	//		if (yF_best[k][j] == 1)
	//			j_k = k;
	//	}

	//	if (i_k >= 0 && j_k >= 0)
	//	{
	//		if (i_k > j_k)
	//		{
	//			submodel.add(xF[i_k][i] == 1);
	//			submodel.add(yF[i_k][i] == 0);
	//		}
	//		else if (i_k < j_k)
	//		{
	//			submodel.add(xF[j_k][j] == 0);
	//			submodel.add(yF[j_k][j] == 1);
	//		}
	//	}

	//		}


	//**********************************//
	//            Constraints           //
	//**********************************//
	//model.add(vF[nbCrane - 1] - vF[nbCrane - 2] <= 0);

	////约束（4）// 所有任务被分配到QC上
	//IloRangeArray  c4(env);
	//for (i = 0; i < nbTask; i++)
	//{
	//	IloExpr  epa(env);
	//	for (k = 0; k < nbCrane; k++)
	//		epa += xF[k][i] + yF[k][i];
	//	c4.add(epa == 1);
	//	epa.end();
	//}
	//submodel.add(c4);
	//c4.end();

	//变量固定
	for (k = QC; k <= QC + 1; k++)
		for (i = 0; i < nbTask; i++)
			if (yF_best[k][i]==1)
			{	submodel.add(xF[k][i] + yF[k][i] == 1);
			}


	for (k = QC; k <= QC + 1; k++)
		submodel.add(CF - ckF[k] >= 0);

	//for (k = 0; k < nbCrane; k++)
	//	submodel.add(chaF[k] - ckF[k] + UB>= 0);

	IloRangeArray  c2(env);
	for (k = QC; k <= QC + 1; k++)
	{
		IloExpr  epa(env);
		//epa += CF - t0kF[k] - gammaF[k];

		epa += ckF[k] - t0kF[k] - gammaF[k];

		for (i = 0; i < nbTask; i++)
			epa -= (nbQ[i] / nbs)*(xF[k][i] + yF[k][i]);
		for (j = 0; j < nbBay; j++)
		{
			epa -= j*endCzF[k][j];
			epa -= wF[k][j] + CwF[k][j];
		}
		//epa -= thetaF[k];

		//for (j = 0; j < nbLocation[i]; j++)
		for (j = 0; j < nbBay; j++)
			epa += j*zF[k][j];

		c2.add(epa >= 0);
		epa.end();
	}
	submodel.add(c2);
	c2.end();

	// thetaF 取值
	for (k = QC; k <= QC + 1; k++)
		for (j = 0; j < nbTask; j++)
		{
		//model.add(thetaF[k] - (nbLocation[j] - 1) * xF[k][j] >= 0);
		IloExpr  epa(env);
		for (i = 0; i < nbBay; i++)
			epa += i*endCzF[k][i];
		submodel.add(epa - (nbLocation[j] - 1) * xF[k][j] >= 0);
		submodel.add(epa - nbLocation[j] * yF[k][j] >= 0);
		epa.end();
		}


	for (k = QC; k <= QC + 1; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbBay; i++)
			epa += i*zF[k][i];
		submodel.add(epa - nbb[k] + 1 <= 0);
		epa.end();
	}


	for (k = QC; k <= QC + 1; k++)
		for (i = 0; i < nbTask; i++)
		{
		submodel.add(yF[k][i] - vF[k] <= 0);
		//model.add(nbLocation[i] * yF[k][i] - thetaF[k] <= 0);
		}
	for (k = QC; k <= QC; k++)
	{
		//model.add(thetaF[k + 1] - thetaF[k] >= 2);
		IloExpr  epa(env);
		for (i = 0; i < nbBay; i++)
			epa += i*endCzF[k + 1][i] - i*endCzF[k][i];
		submodel.add(epa >= 2);
		epa.end();
	}

	//约束（6）// zF 小于最小的xF的bay
	IloRangeArray  c6(env);
	for (i = 0; i < nbBay; i++)
		for (k = QC; k <= QC + 1; k++)
		{
		IloExpr  epa(env);
		IloExpr  epa2(env);
		for (j = 0; j < nbTask; j++)
			if (nbLocation[j] - 1 < i)
			{
			epa += yF[k][j];
			epa2 += xF[k][j];
			}
		for (j = i; j < nbBay; j++)
		{
			epa += 100 * CzF[k][j];
			epa2 += 100 * zF[k][j];
		}


		c6.add(epa <= 100);
		c6.add(epa2 <= 100);
		epa.end();
		epa2.end();
		}
	submodel.add(c6);
	c6.end();



	//约束（7） zF unique
	IloRangeArray  c7(env);
	for (k = QC; k <= QC + 1; k++)
	{
		IloExpr  epa(env);
		IloExpr  epa2(env);
		IloExpr  epa3(env);
		for (i = 0; i < nbBay; i++)
		{
			epa += zF[k][i];
			epa2 += CzF[k][i];
			epa3 += endCzF[k][i];
		}
		epa2 -= vF[k];
		c7.add(epa == 1);
		c7.add(epa2 == 0);
		c7.add(epa3 == 1);
		epa.end();
		epa2.end();
		epa3.end();
	}
	submodel.add(c7);
	c7.end();


	//	建立约束(8)  z 和 z之间隔 /delta +1
	IloRangeArray  c8(env);
	for (k = QC; k <= QC; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbBay; i++)
			epa += i*zF[k + 1][i] - i*zF[k][i];
		c8.add(epa >= 2);
		epa.end();

		IloExpr  epa2(env);
		for (i = 0; i < nbBay; i++)
			epa2 += i*CzF[k + 1][i] - i*CzF[k][i];
		epa2 -= nbBay * (vF[k] + vF[k + 1] - 2);
		c8.add(epa2 >= 2);
		epa2.end();

	}
	submodel.add(c8);
	c8.end();

	for (k = QC; k <= QC + 1; k++)
	{
		IloExpr  epa1(env);
		IloExpr  epa2(env);
		epa1 += t0kF[k];
		for (i = 0; i < nbBay; i++)
			epa2 += i*zF[k][i];
		submodel.add(epa1 + epa2 - nbb[k] + 1 >= 0);
		submodel.add(epa1 - epa2 + nbb[k] - 1 >= 0);
		epa1.end();
		epa2.end();
	}


	//QC travel limits
	IloRangeArray  v00(env);
	for (k = QC; k <= QC + 1; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbTask; i++)
		{
			if (nbLocation[i] < 2 * k + 1)
				epa += xF[k][i] + yF[k][i];
			if (nbLocation[i] > nbBay - 2 * (nbCrane - k - 1))
				epa += xF[k][i] + yF[k][i];
		}
		v00.add(epa <= 0);
		epa.end();
	}
	submodel.add(v00);
	v00.end();



	////////////////sub-problem

	// vF 取值(2)
	for (k = QC; k <= QC + 1; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbTask; i++)
			epa += yF[k][i];
		epa -= vF[k];
		submodel.add(epa >= 0);
		epa.end();
	}
	for (k = QC; k <= QC + 1; k++)
		for (i = 0; i < nbTask; i++)
			submodel.add(yF[k][i] - vF[k] <= 0);

	for (k = QC; k <= QC + 1; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbTask; i++)
			epa += xF[k][i];
		submodel.add(epa >= 1);
		epa.end();
	}
	//// vF 取值(3) vF与violation关系
	//for (i = 0; i < nbTask - 1; i++)
	//	for (j = i + 1; j < nbTask; j++)
	//		if (nbprecR[i][j] == 1)
	//		{
	//	for (k = 0; k < nbCrane; k++)
	//	{
	//		IloExpr epa(env);
	//		IloExpr epa2(env);

	//		epa += vF[k] + 1;
	//		epa2 += vF[k] + 1;

	//		if (k < nbCrane - 1)
	//		{
	//			epa -= xF[k][i] + yF[k][i];
	//			for (int kk = k + 1; kk < nbCrane; kk++)
	//				epa -= xF[kk][j] + yF[kk][j];
	//			submodel.add(epa >= 0);
	//		}

	//		if (k > 0)
	//		{
	//			epa2 -= xF[k][j] + yF[k][j];
	//			for (int kk = 0; kk < k; kk++)
	//				epa2 -= xF[kk][i] + yF[kk][i];
	//			submodel.add(epa2 >= 0);
	//		}
	//		epa.end();
	//		epa2.end();
	//	}
	//		}


	//precedence
	for (i = 0; i < nbTask - 1; i++)
		for (j = i + 1; j < nbTask; j++)
			if (nbprecR[i][j] == 1)
			{
		for (k = QC; k <= QC + 1; k++)
		{

			{
				IloExpr epa(env);
				epa += yF[k][i];
				for (int kk = k; kk <= QC+1; kk++)
					epa -= yF[kk][j];
				submodel.add(epa <= 0);
				epa.end();
			}

			{
				IloExpr epa2(env);
				epa2 += xF[k][j];
				for (int kk = k; kk <= QC+1; kk++)
					epa2 -= xF[kk][i];
				submodel.add(epa2 <= 0);
				epa2.end();
			}

			if (k == QC)
			{
				IloExpr epa3(env);
				epa3 += xF[k][i];
				for (int kk = k + 1; kk <= QC+1; kk++)
					epa3 += xF[kk][j];
				//epa3 += xF[kk][j] + yF[kk][j];
				submodel.add(epa3 <= 1);
				epa3.end();
			}


			if (k == QC)
			{
				IloExpr epa4(env);
				epa4 += yF[k][j];

				//for (int kk = k + 1; kk < nbCrane; kk++)
				//	epa4 += yF[kk][i];
				//epa4 -= 1;


				for (int kk = k + 1; kk <= QC+1; kk++)
					epa4 -= xF[kk][i];
				for (int kk = QC; kk <= k; kk++)
					epa4 -= xF[kk][i] + yF[kk][i];
				submodel.add(epa4 <= 0);
				epa4.end();
			}

			for (int kk = QC; kk <= k; kk++)
			{
				submodel.add(yF[k][j] + xF[kk][i] - vF[kk] <= 1);
			}

		}
			}

	////back time gammaF
	for (k = QC; k <= QC + 1; k++)
	{
		IloExpr  epa(env);
		epa += gammaF[k];
		for (i = 0; i < nbBay; i++)
			epa += i*CzF[k][i] - i*endCzF[k][i];
		epa += nbBay - nbBay*vF[k];
		submodel.add(epa >= 0);
		epa.end();
	}

	//////////////////////////////////////
	//////waiting time ///////
	/////////////////////////////////////

	///////forward trip
	for (k = QC; k <= QC + 1; k++)
		for (i = 0; i < nbBay; i++)
		{
		IloExpr  epa(env);
		epa += t0kF[k];
		for (j = 0; j < nbTask; j++)
			if (nbLocation[j] - 1 <= i)
				epa += (nbQ[j] / nbs)*xF[k][j];
		for (j = 0; j <= i; j++)
		{

			epa += wF[k][j];
			epa -= j*zF[k][j];
			//epa -= 1000 * CzF[k][j];
		}
		submodel.add(epa + i - TF[k][i] == 0);
		epa.end();

		}

	/////约束（3）
	IloRangeArray  c3(env);
	for (i = 0; i < nbBay - 2; i++)
		for (k = QC; k <= QC; k++)
		{
		IloExpr  epa(env);

		epa += TF[k][i] - TF[k + 1][i + 2];
		epa += 1600;
		for (j = 0; j <= i + 2; j++)
			epa -= 1600 * zF[k + 1][j];
		for (j = 0; j < i; j++)
			epa += 1600 * endCzF[k][j];

		c3.add(epa >= 0);
		epa.end();
		}
	submodel.add(c3);
	c3.end();

	//////retracing trip
	for (k = QC; k <= QC + 1; k++)
		for (i = 0; i < nbBay; i++)
		{
		IloExpr  epa(env);
		epa += t0kF[k];
		for (j = 0; j < nbTask; j++)
			epa += (nbQ[j] / nbs)*xF[k][j];
		for (j = 0; j < nbBay; j++)
		{
			epa += 2 * j*endCzF[k][j] - j*zF[k][j];
			epa += wF[k][j];
		}


		for (j = 0; j < nbTask; j++)
			if (nbLocation[j] - 1 >= i)
				epa += (nbQ[j] / nbs)*yF[k][j];
		for (j = i; j < nbBay; j++)
		{
			epa += CwF[k][j];
		}

		submodel.add(epa - i - CTF[k][i] == 0);
		epa.end();

		}

	for (k = QC; k <= QC; k++)
		for (i = 0; i < nbBay - 2; i++)
		{
		IloExpr  epa(env);
		epa += CTF[k + 1][i + 2] - CTF[k][i];

		epa += 6000 - 2000 * vF[k] - 2000 * vF[k + 1];
		for (j = 0; j < i; j++)
			epa += 2000 * endCzF[k][j];
		for (j = 0; j <= i + 2; j++)
			epa -= 2000 * CzF[k + 1][j];

		submodel.add(epa >= 0);
		epa.end();
		}

	// forward and retrace
	for (k = QC; k <= QC; k++)
	{
		IloExpr  epa(env);
		for (j = 0; j < nbBay; j++)
			epa += j * CzF[k + 1][j] - j * endCzF[k][j];
		epa += 2000 + 2000 * vF[k] - 2000 * vF[k + 1];
		submodel.add(epa >= 2);
		epa.end();

	}

	//**********************************//
	//            开始求解		        //
	//**********************************//
	IloCplex cplex(env);
	cplex.extract(submodel);

#ifdef UserActive
	cplex.use(CandyUserCallback(env));
	//cplex.setParam(IloCplex::MIPEmphasis, 2);
	//cplex.setParam(IloCplex::HeurFreq, 1);
	//cplex.setParam(IloCplex::ParallelMode, 1);
	cplex.setParam(IloCplex::Threads, 4);
#endif



	//cplex.setParam(IloCplex::MIPEmphasis ,CPX_MIPEMPHASIS_HIDDENFEAS);
	//cplex.setParam(IloCplex::EpGap,0.01);
	//cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::TiLim, 7200);


	//cout<<endl<<"number of binary variables = "<<cplex.getNbinVars()<<endl;
	//cout<<endl<<"number of constraints = "<<cplex.getNrows()<<endl;

	BOOL h0;
	h0 = cplex.solve();
	//cout<<"h0 "<<h0<<endl;
	if (!h0)
	{
		cout << "\nno feasible solution has been found" << endl;
		cplex.clearModel();
		cplex.clear();
		cplex.end();
		submodel.end();

		return FALSE;
	}
	//**********************************//
	//             记录最好解           //
	//**********************************//

	//*ObjVal = cplex.getBestObjValue();
	*ObjVal = cplex.getValue(CF);

	for (k = QC; k <= QC + 1; k++)
	{
		for (i = 0; i < nbTask; i++)
		{
			if (cplex.isExtracted(xF[k][i]))
			{
				CxF_best[k][i] = cplex.getValue(xF[k][i]);
				//cout<< wF_best[k][i] << "; ";
			}
			else
				CxF_best[k][i] = 0;

		}//cout << endl;
	}
	//cout << endl << endl;

	cout << "final1:  " << *ObjVal << endl;
	for (k = QC; k <= QC + 1; k++)
	{
		cout << k << ":  ";
		for (i = 0; i < nbTask; i++)
		{
			if (yF_best[k][i] == 1)
			{
				if (CxF_best[k][i] == 1)
					cout << i << " ";
				else
					cout << i << "r ";
			}


		}cout << endl;
	}
	cout << endl;

	//int recb = 1;
	//cout << recb << ":  ";
	//for (i = 0; i < nbTask; i++)
	//{
	//	if (nbLocation[i] != recb)
	//	{
	//		recb = nbLocation[i];
	//		cout << endl << recb << ":  ";
	//	}
	//	cout << i << "-";
	//	for (k = 0; k < nbCrane; k++)
	//		if (yF_best[k][i] == 1)
	//			cout << k << " ";
	//}
	//cout << endl;

	//cout << "QC_CF_best" << endl;
	//for (k = 0; k < nbCrane; k++)
	//{
	//	cout << cplex.getValue(QC_CF[k]) << "  ";
	//}
	//cout << endl << endl;

	//cout << "Task_CF_best" << endl;
	//for (i = 0; i < nbTask; i++)
	//{
	//	cout << cplex.getValue(Task_CF[i]) << "  ";
	//}
	//cout << endl << endl;

	//for (k = 0; k < nbCrane; k++)
	//{
	//	cout << k << ":  ";
	//	for (i = 0; i < nbBay; i++)
	//	{
	//		if (cplex.isExtracted(wF[k][i]))
	//		{

	//			wF_best[k][i] = cplex.getValue(wF[k][i]);
	//			cout << wF_best[k][i] << "; ";
	//		}
	//		else
	//			wF_best[k][i] = 0;

	//	}cout << endl;
	//}
	//cout << endl;

	//for (k = 0; k < nbCrane; k++)
	//{
	//	cout << k << ":  ";
	//	for (i = 0; i < nbBay; i++)
	//	{
	//		if (cplex.isExtracted(CwF[k][i]))
	//		{
	//			CwF_best[k][i] = cplex.getValue(CwF[k][i]);
	//			cout << CwF_best[k][i] << "; ";
	//		}
	//		else
	//			CwF_best[k][i] = 0;

	//	}cout << endl;
	//}
	//cout << endl;


	//for (k = 0; k < nbCrane; k++)
	//{
	//	cout << k << ":  " << cplex.getValue(ckF[k]) << endl;
	//}
	//cout << endl;


	//cout << "CF_best: " << *CF_best << endl;
	//cout << "ObjVal: " << *ObjVal << endl;

	cplex.clearModel();
	cplex.clear();
	cplex.end();
	submodel.end();

	return TRUE;

}

IloInt Split_SP_Callback(IloInt solmode, IloCplex cplex, IloNumArray nbQ, int*  nbLocation, BoolMatrix nbprecR, IloNumArray2 yF_best, IloInt QC)
{

	//xF: forward trip, yF: retracing trip


	IloEnv env = cplex.getEnv();
	IloInt i, j, k;
	IloInt ObjVal;

	//问题模型
	IloModel submodel(env);

	IloIntVar CF(env, 0, 1800);

	BoolVarMatrix xF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		xF[k] = IloBoolVarArray(env, nbTask);
	}
	BoolVarMatrix yF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		yF[k] = IloBoolVarArray(env, nbTask);
	}
	IloBoolVarArray vF(env,nbCrane);

	//问题模型
	IloNumVarArray t0kF(env, nbCrane, 0, 100);
	IloNumVarArray gammaF(env, nbCrane, 0, 100);//

	IloIntVarArray thetaF(env, nbCrane, 0, 3000);

	IloNumVarArray ckF(env, nbCrane, 0, 3000);//
	//IloNumVarArray chaF(env, nbCrane, 0, 3000);//

	NumVarMatrix wF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		wF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);

	NumVarMatrix CwF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		CwF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);

	BoolVarMatrix zF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		zF[k] = IloBoolVarArray(env, nbBay);
	}

	BoolVarMatrix CzF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		CzF[k] = IloBoolVarArray(env, nbBay);
	}

	BoolVarMatrix endCzF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		endCzF[k] = IloBoolVarArray(env, nbBay);
	}

	//leaving time at each bay
	NumVarMatrix TF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		TF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);
	NumVarMatrix CTF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		CTF[k] = IloNumVarArray(env, nbBay, 0, IloInfinity);


	//**********************************//
	//           objective 原问题目标函数        //
	//**********************************//

	//IloExpr obj1(env); 


	IloExpr  obj2(env);
	//  建立子问题目标函数表达式 

	obj2 += CF;

	//for (k = 0; k < nbCrane; k++)
	//	obj2 += chaF[k];

	//	将目标函数加入到原问题模型
	submodel.add(IloMinimize(env, obj2));//
	//obj1.end();
	obj2.end();




	if (solmode==1)
	{
		//**********************************//
		//            variable fixing           //
		//**********************************//
		//for (i = 0; i < nbTask - 1; i++)
		//	for (j = i + 1; j < nbTask; j++)
		//		if (nbprecR[i][j] == 1)
		//		{
		//	int i_k, j_k;
		//	i_k = -1;
		//	j_k = -1;
		//	for (k = QC; k <= QC + 1; k++)
		//	{
		//		if (yF_best[k][i] == 1)
		//			i_k = k;
		//		if (yF_best[k][j] == 1)
		//			j_k = k;
		//	}

		//	if (i_k >= 0 && j_k >= 0)
		//	{
		//		if (i_k > j_k)
		//		{
		//			submodel.add(xF[i_k][i] == 1);
		//			submodel.add(yF[i_k][i] == 0);
		//		}
		//		else if (i_k < j_k)
		//		{
		//			submodel.add(xF[j_k][j] == 0);
		//			submodel.add(yF[j_k][j] == 1);
		//		}
		//	}

		//		}

		//**********************************//
		//            Constraints           //
		//**********************************//
		//model.add(vF[nbCrane - 1] - vF[nbCrane - 2] <= 0);

		////约束（4）// 所有任务被分配到QC上
		//IloRangeArray  c4(env);
		//for (i = 0; i < nbTask; i++)
		//{
		//	IloExpr  epa(env);
		//	for (k = 0; k < nbCrane; k++)
		//		epa += xF[k][i] + yF[k][i];
		//	c4.add(epa == 1);
		//	epa.end();
		//}
		//submodel.add(c4);
		//c4.end();

		//变量固定
		for (k = QC; k <= QC + 1; k++)
			for (i = 0; i < nbTask; i++)

			{
			if (yF_best[k][i] == 1)
			{
				submodel.add(xF[k][i] + yF[k][i] == 1);
			}
			else
			{
				submodel.add(xF[k][i] + yF[k][i] == 0);
			}
			}


		for (k = QC; k <= QC + 1; k++)
			submodel.add(CF - ckF[k] >= 0);

		//for (k = 0; k < nbCrane; k++)
		//	submodel.add(chaF[k] - ckF[k] + UB>= 0);

		IloRangeArray  c2(env);
		for (k = QC; k <= QC + 1; k++)
		{
			IloExpr  epa(env);
			//epa += CF - t0kF[k] - gammaF[k];

			epa += ckF[k] - t0kF[k] - gammaF[k];

			for (i = 0; i < nbTask; i++)
				epa -= nbQ[i] * (xF[k][i] + yF[k][i]);
			for (j = 0; j < nbBay; j++)
			{
				epa -= j*endCzF[k][j];
				epa -= wF[k][j] + CwF[k][j];
			}
			//epa -= thetaF[k];

			//for (j = 0; j < nbLocation[i]; j++)
			for (j = 0; j < nbBay; j++)
				epa += j*zF[k][j];

			c2.add(epa >= 0);
			epa.end();
		}
		submodel.add(c2);
		c2.end();

		// thetaF 取值
		for (k = QC; k <= QC + 1; k++)
			for (j = 0; j < nbTask; j++)
			{
			//model.add(thetaF[k] - (nbLocation[j] - 1) * xF[k][j] >= 0);
			IloExpr  epa(env);
			for (i = 0; i < nbBay; i++)
				epa += i*endCzF[k][i];
			submodel.add(epa - (nbLocation[j] - 1) * xF[k][j] >= 0);
			submodel.add(epa - nbLocation[j] * yF[k][j] >= 0);
			epa.end();
			}


		//for (k = QC; k <= QC + 1; k++)
		//{
		//	IloExpr  epa(env);
		//	for (i = 0; i < nbBay; i++)
		//		epa += i*zF[k][i];
		//	submodel.add(epa - 2 * k <= 0);
		//	epa.end();
		//}


		for (k = QC; k <= QC + 1; k++)
			for (i = 0; i < nbTask; i++)
			{
			submodel.add(yF[k][i] - vF[k] <= 0);
			//model.add(nbLocation[i] * yF[k][i] - thetaF[k] <= 0);
			}
		for (k = QC; k <= QC; k++)
			if ((1 + safe_margin)*(k + 1) + 1 <= nbBay)
		{
			//model.add(thetaF[k + 1] - thetaF[k] >= 2);
			IloExpr  epa(env);
			for (i = 0; i < nbBay; i++)
				epa += i*endCzF[k + 1][i] - i*endCzF[k][i];
			submodel.add(epa >= 1 + safe_margin);
			epa.end();
		}

		//约束（6）// zF 小于最小的xF的bay
		IloRangeArray  c6(env);
		for (i = 0; i < nbBay; i++)
			for (k = QC; k <= QC + 1; k++)
			{
			IloExpr  epa(env);
			IloExpr  epa2(env);
			for (j = 0; j < nbTask; j++)
				if (nbLocation[j] - 1 < i)
				{
				epa += yF[k][j];
				epa2 += xF[k][j];
				}
			for (j = i; j < nbBay; j++)
			{
				epa += 100 * CzF[k][j];
				epa2 += 100 * zF[k][j];
			}


			c6.add(epa <= 100);
			c6.add(epa2 <= 100);
			epa.end();
			epa2.end();
			}
		submodel.add(c6);
		c6.end();



		//约束（7） zF unique
		IloRangeArray  c7(env);
		for (k = QC; k <= QC + 1; k++)
		{
			IloExpr  epa(env);
			IloExpr  epa2(env);
			IloExpr  epa3(env);
			for (i = 0; i < nbBay; i++)
			{
				epa += zF[k][i];
				epa2 += CzF[k][i];
				epa3 += endCzF[k][i];
			}
			epa2 -= vF[k];
			c7.add(epa == 1);
			c7.add(epa2 == 0);
			c7.add(epa3 == 1);
			epa.end();
			epa2.end();
			epa3.end();
		}
		submodel.add(c7);
		c7.end();


		//	建立约束(8)  z 和 z之间隔 /delta +1
		IloRangeArray  c8(env);
		for (k = QC; k <= QC; k++)
		{
			IloExpr  epa(env);
			for (i = 0; i < nbBay; i++)
				epa += i*zF[k + 1][i] - i*zF[k][i];
			c8.add(epa >= 2);
			epa.end();

			IloExpr  epa2(env);
			for (i = 0; i < nbBay; i++)
				epa2 += i*CzF[k + 1][i] - i*CzF[k][i];
			epa2 -= nbBay * (vF[k] + vF[k + 1] - 2);
			c8.add(epa2 >= 2);
			epa2.end();

		}
		submodel.add(c8);
		c8.end();

		for (k = QC; k <= QC + 1; k++)
		{
			IloExpr  epa1(env);
			IloExpr  epa2(env);
			epa1 += t0kF[k];
			for (i = 0; i < nbBay; i++)
				epa2 += i*zF[k][i];
			submodel.add(epa1 + epa2 - 2 * k >= 0);
			submodel.add(epa1 - epa2 + 2 * k >= 0);
			epa1.end();
			epa2.end();
		}


		//QC travel limits
		IloRangeArray  v00(env);
		for (k = QC; k <= QC + 1; k++)
		{
			IloExpr  epa(env);
			for (i = 0; i < nbTask; i++)
			{
				if (nbLocation[i] < 2 * k + 1)
					epa += xF[k][i] + yF[k][i];
				if (nbLocation[i] > nbBay - 2 * (nbCrane - k - 1))
					epa += xF[k][i] + yF[k][i];
			}
			v00.add(epa <= 0);
			epa.end();
		}
		submodel.add(v00);
		v00.end();



		////////////////sub-problem

		// vF 取值(2)
		for (k = QC; k <= QC + 1; k++)
		{
			IloExpr  epa(env);
			for (i = 0; i < nbTask; i++)
				epa += yF[k][i];
			epa -= vF[k];
			submodel.add(epa >= 0);
			epa.end();
		}
		for (k = QC; k <= QC + 1; k++)
			for (i = 0; i < nbTask; i++)
				submodel.add(yF[k][i] - vF[k] <= 0);

		for (k = QC; k <= QC + 1; k++)
		{
			IloExpr  epa(env);
			for (i = 0; i < nbTask; i++)
				epa += xF[k][i];
			submodel.add(epa >= 1);
			epa.end();
		}
		//// vF 取值(3) vF与violation关系
		//for (i = 0; i < nbTask - 1; i++)
		//	for (j = i + 1; j < nbTask; j++)
		//		if (nbprecR[i][j] == 1)
		//		{
		//	for (k = 0; k < nbCrane; k++)
		//	{
		//		IloExpr epa(env);
		//		IloExpr epa2(env);

		//		epa += vF[k] + 1;
		//		epa2 += vF[k] + 1;

		//		if (k < nbCrane - 1)
		//		{
		//			epa -= xF[k][i] + yF[k][i];
		//			for (int kk = k + 1; kk < nbCrane; kk++)
		//				epa -= xF[kk][j] + yF[kk][j];
		//			submodel.add(epa >= 0);
		//		}

		//		if (k > 0)
		//		{
		//			epa2 -= xF[k][j] + yF[k][j];
		//			for (int kk = 0; kk < k; kk++)
		//				epa2 -= xF[kk][i] + yF[kk][i];
		//			submodel.add(epa2 >= 0);
		//		}
		//		epa.end();
		//		epa2.end();
		//	}
		//		}


		//precedence
		for (i = 0; i < nbTask - 1; i++)
			for (j = i + 1; j < nbTask; j++)
				if (nbprecR[i][j] == 1)
				{
			for (k = QC; k <= QC + 1; k++)
			{

				{
					IloExpr epa(env);
					epa += yF[k][i];
					for (int kk = k; kk <= QC + 1; kk++)
						epa -= yF[kk][j];
					submodel.add(epa <= 0);
					epa.end();
				}

			{
				IloExpr epa2(env);
				epa2 += xF[k][j];
				for (int kk = k; kk <= QC + 1; kk++)
					epa2 -= xF[kk][i];
				submodel.add(epa2 <= 0);
				epa2.end();
			}

				if (k == QC)
				{
					IloExpr epa3(env);
					epa3 += xF[k][i];
					for (int kk = k + 1; kk <= QC + 1; kk++)
						epa3 += xF[kk][j];
					//epa3 += xF[kk][j] + yF[kk][j];
					submodel.add(epa3 <= 1);
					epa3.end();
				}


				if (k == QC)
				{
					IloExpr epa4(env);
					epa4 += yF[k][j];

					//for (int kk = k + 1; kk < nbCrane; kk++)
					//	epa4 += yF[kk][i];
					//epa4 -= 1;


					for (int kk = k + 1; kk <= QC + 1; kk++)
						epa4 -= xF[kk][i];
					for (int kk = QC; kk <= k; kk++)
						epa4 -= xF[kk][i] + yF[kk][i];
					submodel.add(epa4 <= 0);
					epa4.end();
				}

				for (int kk = QC; kk <= k; kk++)
				{
					submodel.add(yF[k][j] + xF[kk][i] - vF[kk] <= 1);
				}

			}
				}

		////back time gammaF
		for (k = QC; k <= QC + 1; k++)
		{
			IloExpr  epa(env);
			epa += gammaF[k];
			for (i = 0; i < nbBay; i++)
				epa += i*CzF[k][i] - i*endCzF[k][i];
			epa += nbBay - nbBay*vF[k];
			submodel.add(epa >= 0);
			epa.end();
		}

		//////////////////////////////////////
		//////waiting time ///////
		/////////////////////////////////////

		///////forward trip
		for (k = QC; k <= QC + 1; k++)
			for (i = 0; i < nbBay; i++)
			{
			IloExpr  epa(env);
			epa += t0kF[k];
			for (j = 0; j < nbTask; j++)
				if (nbLocation[j] - 1 <= i)
					epa += nbQ[j] * xF[k][j];
			for (j = 0; j <= i; j++)
			{

				epa += wF[k][j];
				epa -= j*zF[k][j];
				//epa -= 1000 * CzF[k][j];
			}
			submodel.add(epa + i - TF[k][i] == 0);
			epa.end();

			}

		/////约束（3）
		IloRangeArray  c3(env);
		for (i = 0; i < nbBay - 2; i++)
			for (k = QC; k <= QC; k++)
			{
			IloExpr  epa(env);

			epa += TF[k][i] - TF[k + 1][i + 2];
			epa += 1600;
			for (j = 0; j <= i + 2; j++)
				epa -= 1600 * zF[k + 1][j];
			for (j = 0; j < i; j++)
				epa += 1600 * endCzF[k][j];

			c3.add(epa >= 0);
			epa.end();
			}
		submodel.add(c3);
		c3.end();

		//////retracing trip
		for (k = QC; k <= QC + 1; k++)
			for (i = 0; i < nbBay; i++)
			{
			IloExpr  epa(env);
			epa += t0kF[k];
			for (j = 0; j < nbTask; j++)
				epa += nbQ[j] * xF[k][j];
			for (j = 0; j < nbBay; j++)
			{
				epa += 2 * j*endCzF[k][j] - j*zF[k][j];
				epa += wF[k][j];
			}


			for (j = 0; j < nbTask; j++)
				if (nbLocation[j] - 1 >= i)
					epa += nbQ[j] * yF[k][j];
			for (j = i; j < nbBay; j++)
			{
				epa += CwF[k][j];
			}

			submodel.add(epa - i - CTF[k][i] == 0);
			epa.end();

			}

		for (k = QC; k <= QC; k++)
			for (i = 0; i < nbBay - 2; i++)
			{
			IloExpr  epa(env);
			epa += CTF[k + 1][i + 2] - CTF[k][i];

			epa += 6000 - 2000 * vF[k] - 2000 * vF[k + 1];
			for (j = 0; j < i; j++)
				epa += 2000 * endCzF[k][j];
			for (j = 0; j <= i + 2; j++)
				epa -= 2000 * CzF[k + 1][j];

			submodel.add(epa >= 0);
			epa.end();
			}

		// forward and retrace
		for (k = QC; k <= QC; k++)
		{
			IloExpr  epa(env);
			for (j = 0; j < nbBay; j++)
				epa += j * CzF[k + 1][j] - j * endCzF[k][j];
			epa += 2000 + 2000 * vF[k] - 2000 * vF[k + 1];
			submodel.add(epa >= 2);
			epa.end();

		}
	}
	else
	{
		//**********************************//
		//            variable fixing           //
		//**********************************//
		for (i = 0; i < nbTask - 1; i++)
			for (j = i + 1; j < nbTask; j++)
				if (nbprecR[i][j] == 1)
				{
			int i_k, j_k;
			for (k = 0; k < nbCrane; k++)
			{
				if (yF_best[k][i] == 1)
					i_k = k;
				if (yF_best[k][j] == 1)
					j_k = k;
			}

			if (i_k > j_k)
			{
				submodel.add(xF[i_k][i] == 1);
				submodel.add(yF[i_k][i] == 0);
			}
			else if (i_k < j_k)
			{
				submodel.add(xF[j_k][j] == 0);
				submodel.add(yF[j_k][j] == 1);
			}
				}


		//**********************************//
		//            Constraints           //
		//**********************************//
		//model.add(vF[nbCrane - 1] - vF[nbCrane - 2] <= 0);

		//变量固定
		for (k = 0; k < nbCrane; k++)
			for (i = 0; i < nbTask; i++)
				submodel.add(xF[k][i] + yF[k][i] == yF_best[k][i]);

		for (k = 0; k < nbCrane; k++)
			submodel.add(CF - ckF[k] >= 0);

		//for (k = 0; k < nbCrane; k++)
		//	submodel.add(chaF[k] - ckF[k] + UB>= 0);

		IloRangeArray  c2(env);
		for (k = 0; k < nbCrane; k++)
		{
			IloExpr  epa(env);
			//epa += CF - t0kF[k] - gammaF[k];

			epa += ckF[k] - t0kF[k] - gammaF[k];

			for (i = 0; i < nbTask; i++)
				epa -= nbQ[i]*(xF[k][i] + yF[k][i]);
			for (j = 0; j < nbBay; j++)
			{
				epa -= QCmove_time*j*endCzF[k][j];
				epa -= wF[k][j] + CwF[k][j];
			}
			//epa -= thetaF[k];

			//for (j = 0; j < nbLocation[i]; j++)
			for (j = 0; j < nbBay; j++)
				epa += QCmove_time*j*zF[k][j];

			c2.add(epa >= 0);
			epa.end();
		}
		submodel.add(c2);
		c2.end();

		// thetaF 取值
		for (k = 0; k < nbCrane; k++)
			for (j = 0; j < nbTask; j++)
			{
			//model.add(thetaF[k] - (nbLocation[j] - 1) * xF[k][j] >= 0);
			IloExpr  epa(env);
			for (i = 0; i < nbBay; i++)
				epa += i*endCzF[k][i];
			submodel.add(epa - (nbLocation[j] - 1) * xF[k][j] >= 0);
			submodel.add(epa - nbLocation[j] * yF[k][j] >= 0);
			epa.end();
			}
		//for (k = 0; k < nbCrane; k++)
		//{
		//	IloExpr  epa(env);
		//	for (i = 0; i < nbBay; i++)
		//		epa += i*endCzF[k][i];
		//	model.add(epa - thetaF[k] == 0);
		//	epa.end();
		//}

		//for (k = 0; k < nbCrane; k++)
		//{
		//	IloExpr  epa(env);
		//	for (i = 0; i < nbBay; i++)
		//		epa += i*zF[k][i];
		//	submodel.add(epa - 2*k <= 0);
		//	epa.end();
		//}


		for (k = 0; k < nbCrane; k++)
			for (i = 0; i < nbTask; i++)
			{
			submodel.add(yF[k][i] - vF[k] <= 0);
			//model.add(nbLocation[i] * yF[k][i] - thetaF[k] <= 0);
			}
		for (k = 0; k < nbCrane - 1; k++)
			if ((1 + safe_margin)*(k + 2) <= nbBay)
		{
			//model.add(thetaF[k + 1] - thetaF[k] >= 2);
			IloExpr  epa(env);
			for (i = 0; i < nbBay; i++)
				epa += i*endCzF[k + 1][i] - i*endCzF[k][i];
			submodel.add(epa >= 1 + safe_margin);
			epa.end();
		}


		//约束（4）// 所有任务被分配到QC上
		IloRangeArray  c4(env);
		for (i = 0; i < nbTask; i++)
		{
			IloExpr  epa(env);
			for (k = 0; k < nbCrane; k++)
				epa += xF[k][i] + yF[k][i];
			c4.add(epa == 1);
			epa.end();
		}
		submodel.add(c4);
		c4.end();

		////约束（5）//endCzF 取在最后一个xF处
		//IloRangeArray  c5(env);
		//for (i = 0; i < nbBay; i++)
		//	for (k = 0; k < nbCrane; k++)
		//	{
		//	IloExpr  epa(env);
		//	for (j = 0; j < nbTask; j++)
		//		if (nbLocation[j] - 1 == i)
		//			epa += xF[k][j];
		//	epa -= endCzF[k][i];
		//	c5.add(epa >= 0);
		//	epa.end();
		//	}
		//model.add(c5);
		//c5.end();

		//约束（6）// zF 小于最小的xF的bay
		IloRangeArray  c6(env);
		for (i = 0; i < nbBay; i++)
			for (k = 0; k < nbCrane; k++)
			{
			IloExpr  epa(env);
			IloExpr  epa2(env);
			for (j = 0; j < nbTask; j++)
				if (nbLocation[j] - 1 < i)
				{
				epa += yF[k][j];
				epa2 += xF[k][j];
				}
			for (j = i; j < nbBay; j++)
			{
				epa += 100 * CzF[k][j];
				epa2 += 100 * zF[k][j];
			}


			c6.add(epa <= 100);
			c6.add(epa2 <= 100);
			epa.end();
			epa2.end();
			}
		submodel.add(c6);
		c6.end();



		//约束（7） zF unique
		IloRangeArray  c7(env);
		for (k = 0; k < nbCrane; k++)
		{
			IloExpr  epa(env);
			IloExpr  epa2(env);
			IloExpr  epa3(env);
			for (i = 0; i < nbBay; i++)
			{
				epa += zF[k][i];
				epa2 += CzF[k][i];
				epa3 += endCzF[k][i];
			}
			epa2 -= vF[k];
			c7.add(epa == 1);
			c7.add(epa2 == 0);
			c7.add(epa3 == 1);
			epa.end();
			epa2.end();
			epa3.end();
		}
		submodel.add(c7);
		c7.end();


		//	建立约束(8)  z 和 z之间隔 /delta +1
		IloRangeArray  c8(env);
		for (k = 0; k < nbCrane - 1; k++)
			if ((1 + safe_margin)*(k + 1) + 1 <= nbBay)
		{
			IloExpr  epa(env);
			for (i = 0; i < nbBay; i++)
				epa += i*zF[k + 1][i] - i*zF[k][i];
			c8.add(epa >= 1 + safe_margin);
			epa.end();

			IloExpr  epa2(env);
			for (i = 0; i < nbBay; i++)
				epa2 += i*CzF[k + 1][i] - i*CzF[k][i];
			epa2 -= nbBay * (vF[k] + vF[k + 1] - 2);
			c8.add(epa2 >= 1 + safe_margin);
			epa2.end();

		}
		submodel.add(c8);
		c8.end();

		for (k = 0; k < nbCrane; k++)
		{
			IloExpr  epa1(env);
			IloExpr  epa2(env);
			epa1 += t0kF[k];
			for (i = 0; i < nbBay; i++)
				epa2 += QCmove_time*i*zF[k][i];
			epa2 -= QCmove_time*nbb_BD[k] - QCmove_time;
			submodel.add(epa1 - epa2 >= 0);
			submodel.add(epa1 + epa2 >= 0);
			epa1.end();
			epa2.end();
		}


		////QC travel limits
		//IloRangeArray  v00(env);
		//for (k = 0; k < nbCrane; k++)
		//{
		//	IloExpr  epa(env);
		//	for (i = 0; i < nbTask; i++)
		//	{
		//		if (nbLocation[i] < 2 * k + 1)
		//			epa += xF[k][i] + yF[k][i];
		//		if (nbLocation[i] > nbBay - 2 * (nbCrane - k - 1))
		//			epa += xF[k][i] + yF[k][i];
		//	}
		//	v00.add(epa <= 0);
		//	epa.end();
		//}
		//submodel.add(v00);
		//v00.end();



		////////////////sub-problem

		// vF 取值(2)
		for (k = 0; k < nbCrane; k++)
		{
			IloExpr  epa(env);
			for (i = 0; i < nbTask; i++)
				epa += yF[k][i];
			epa -= vF[k];
			submodel.add(epa >= 0);
			epa.end();
		}
		for (k = 0; k < nbCrane; k++)
			for (i = 0; i < nbTask; i++)
				submodel.add(yF[k][i] - vF[k] <= 0);

		for (k = 0; k < nbCrane; k++)
		{
			IloExpr  epa(env);
			for (i = 0; i < nbTask; i++)
				epa += xF[k][i];
			submodel.add(epa >= 1);
			epa.end();
		}
		// vF 取值(3) vF与violation关系
		for (i = 0; i < nbTask - 1; i++)
			for (j = i + 1; j < nbTask; j++)
				if (nbprecR[i][j] == 1)
				{
			for (k = 0; k < nbCrane; k++)
			{
				IloExpr epa(env);
				IloExpr epa2(env);

				epa += vF[k] + 1;
				epa2 += vF[k] + 1;

				if (k < nbCrane - 1)
				{
					epa -= xF[k][i] + yF[k][i];
					for (int kk = k + 1; kk < nbCrane; kk++)
						epa -= xF[kk][j] + yF[kk][j];
					submodel.add(epa >= 0);
				}

				if (k > 0)
				{
					epa2 -= xF[k][j] + yF[k][j];
					for (int kk = 0; kk < k; kk++)
						epa2 -= xF[kk][i] + yF[kk][i];
					submodel.add(epa2 >= 0);
				}
				epa.end();
				epa2.end();
			}
				}


		//precedence
		for (i = 0; i < nbTask - 1; i++)
			for (j = i + 1; j < nbTask; j++)
				if (nbprecR[i][j] == 1)
				{
			for (k = 0; k < nbCrane; k++)
			{

				{
					IloExpr epa(env);
					epa += yF[k][i];
					for (int kk = k; kk < nbCrane; kk++)
						epa -= yF[kk][j];
					submodel.add(epa <= 0);
					epa.end();
				}

			{
				IloExpr epa2(env);
				epa2 += xF[k][j];
				for (int kk = k; kk < nbCrane; kk++)
					epa2 -= xF[kk][i];
				submodel.add(epa2 <= 0);
				epa2.end();
			}

				if (k < nbCrane - 1)
				{
					IloExpr epa3(env);
					epa3 += xF[k][i];
					for (int kk = k + 1; kk < nbCrane; kk++)
						epa3 += xF[kk][j];
					//epa3 += xF[kk][j] + yF[kk][j];
					submodel.add(epa3 <= 1);
					epa3.end();
				}


				if (k < nbCrane - 1)
				{
					IloExpr epa4(env);
					epa4 += yF[k][j];

					//for (int kk = k + 1; kk < nbCrane; kk++)
					//	epa4 += yF[kk][i];
					//epa4 -= 1;


					for (int kk = k + 1; kk < nbCrane; kk++)
						epa4 -= xF[kk][i];
					for (int kk = 0; kk <= k; kk++)
						epa4 -= xF[kk][i] + yF[kk][i];
					submodel.add(epa4 <= 0);
					epa4.end();
				}

				for (int kk = 0; kk <= k; kk++)
				{
					submodel.add(yF[k][j] + xF[kk][i] - vF[kk] <= 1);
				}

			}
				}

		////back time gammaF
		for (k = 0; k < nbCrane; k++)
		{
			IloExpr  epa(env);
			epa += gammaF[k];
			for (i = 0; i < nbBay; i++)
				epa += i*CzF[k][i] - i*endCzF[k][i];
			epa += nbBay - nbBay*vF[k];
			submodel.add(epa >= 0);
			epa.end();
		}

		//////////////////////////////////////
		//////waiting time ///////
		/////////////////////////////////////

		///////forward trip
		for (k = 0; k < nbCrane; k++)
			for (i = 0; i < nbBay; i++)
			{
			IloExpr  epa(env);
			epa += t0kF[k];
			for (j = 0; j < nbTask; j++)
				if (nbLocation[j] - 1 <= i)
					epa += nbQ[j] *xF[k][j];
			for (j = 0; j <= i; j++)
			{

				epa += wF[k][j];
				epa -= QCmove_time*j*zF[k][j];
				//epa -= 1000 * CzF[k][j];
			}
			submodel.add(epa + QCmove_time*i - TF[k][i] == 0);
			epa.end();

			}

		/////约束（3）
		IloRangeArray  c3(env);
		for (i = 0; i < nbBay - 1 - safe_margin; i++)
			for (k = 0; k < nbCrane - 1; k++)
			{
			IloExpr  epa(env);

			epa += TF[k][i] - TF[k + 1][i + 1 + safe_margin];
			epa += 1600;
			for (j = 0; j <= i + 1 + safe_margin; j++)
				epa -= 1600 * zF[k + 1][j];
			for (j = 0; j < i; j++)
				epa += 1600 * endCzF[k][j];

			c3.add(epa >= 0);
			epa.end();
			}
		submodel.add(c3);
		c3.end();

		//////retracing trip
		for (k = 0; k < nbCrane; k++)
			for (i = 0; i < nbBay; i++)
			{
			IloExpr  epa(env);
			epa += t0kF[k];
			for (j = 0; j < nbTask; j++)
				epa += nbQ[j]*xF[k][j];
			for (j = 0; j < nbBay; j++)
			{
				epa += 2 * QCmove_time* j*endCzF[k][j] - QCmove_time*j*zF[k][j];
				epa += wF[k][j];
			}


			for (j = 0; j < nbTask; j++)
				if (nbLocation[j] - 1 >= i)
					epa += nbQ[j]*yF[k][j];
			for (j = i; j < nbBay; j++)
			{
				epa += CwF[k][j];
			}

			submodel.add(epa - QCmove_time*i - CTF[k][i] == 0);
			epa.end();

			}

		for (k = 0; k < nbCrane - 1; k++)
			for (i = 0; i < nbBay - 1 - safe_margin; i++)
			{
			IloExpr  epa(env);
			epa += CTF[k + 1][i + 1 + safe_margin] - CTF[k][i];

			epa += 6000 - 2000 * vF[k] - 2000 * vF[k + 1];
			for (j = 0; j < i; j++)
				epa += 2000 * endCzF[k][j];
			for (j = 0; j <= i + 1 + safe_margin; j++)
				epa -= 2000 * CzF[k + 1][j];

			submodel.add(epa >= 0);
			epa.end();
			}

		// forward and retrace
		for (k = 0; k < nbCrane - 1; k++)
		{
			IloExpr  epa(env);
			for (j = 0; j < nbBay; j++)
				epa += j * CzF[k + 1][j] - j * endCzF[k][j];
			epa += 2000 + 2000 * vF[k] - 2000 * vF[k + 1];
			submodel.add(epa >= 1 + safe_margin);
			epa.end();

		}
	}
	

	//**********************************//
	//            开始求解		        //
	//**********************************//
	IloCplex cplex(env);
	cplex.extract(submodel);

#ifdef UserActive
	cplex.use(CandyUserCallback(env));
	//cplex.setParam(IloCplex::MIPEmphasis, 2);
	//cplex.setParam(IloCplex::HeurFreq, 1);
	//cplex.setParam(IloCplex::ParallelMode, 1);
	cplex.setParam(IloCplex::Threads, 4);
#endif



	//cplex.setParam(IloCplex::MIPEmphasis ,CPX_MIPEMPHASIS_HIDDENFEAS);
	//cplex.setParam(IloCplex::EpGap,0.01);
	//cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::TiLim, 7200);


	//cout<<endl<<"number of binary variables = "<<cplex.getNbinVars()<<endl;
	//cout<<endl<<"number of constraints = "<<cplex.getNrows()<<endl;

	BOOL h0;
	h0 = cplex.solve();
	//cout<<"h0 "<<h0<<endl;
	if (!h0)
	{
		cout << "\nno feasible solution has been found" << endl;
		cplex.clearModel();
		cplex.clear();
		cplex.end();
		submodel.end();

		return 0;
	}
	//**********************************//
	//             记录最好解           //
	//**********************************//

	//*ObjVal = cplex.getBestObjValue();
	ObjVal = cplex.getValue(CF);

	//for (k = QC; k <= QC + 1; k++)
	//{
	//	for (i = 0; i < nbTask; i++)
	//	{
	//		if (cplex.isExtracted(xF[k][i]))
	//		{
	//			CxF_best[k][i] = cplex.getValue(xF[k][i]);
	//			//cout<< wF_best[k][i] << "; ";
	//		}
	//		else
	//			CxF_best[k][i] = 0;

	//	}//cout << endl;
	//}
	//cout << endl << endl;


	//for (k = QC; k <= QC + 1; k++)
	//{
	//	cout << k << ":  ";
	//	for (i = 0; i < nbTask; i++)
	//	{
	//		if (yF_best[k][i] == 1)
	//		{
	//			if (CxF_best[k][i] == 1)
	//				cout << i << " ";
	//			else
	//				cout << i << "r ";
	//		}


	//	}cout << endl;
	//}
	//cout << endl;

	cout << "final2:  " << ObjVal << endl;

	cplex.clearModel();
	cplex.clear();
	cplex.end();
	submodel.end();

	return ObjVal;

}
