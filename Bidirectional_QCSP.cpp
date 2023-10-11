#include <ilcplex/ilocplex.h>
#include  <fstream>
#include  <iostream>

//#define TEST
#define ReduConst// use reduced set in constraints
//#define UserActive//activate the usercallback
#define _AFXDLL
//#define OUTBRANCH//输出分支节点处的解
//#define SubTour

//#include "afxtempl.h"
//#include  <afx.h>
//#include <afxdb.h>
#include <time.h>
#include "stdafx.h"
#include <iomanip>


//#ifdef _DEBUG
//#define new DEBUG_NEW
//#endif



// 唯一的应用程序对象

//CWinApp theApp;

using namespace std;

const   IloInt nbTask = 100;
const   IloInt nbBay=24;
const   IloInt nbCrane=6;
const   IloInt nbTime=100;


const   IloInt crane_move_time = 1;

const IloBool ORIGINAL=1;

IloInt SubCutNum;
IloNum ObjVal;
IloNum LB;
IloNum LB0;
IloNum UB;
IloNum GapF;

IloNum GapF2;

#define ModNo  1//选择要调用的模型编号 1 or 2
#define QCmove_time 1
#define safe_margin 1

#define waiting_constraints

#define instance_no 10//选择测试的算例数量
#define start_instance 4//开始测试的算例

#define Normal_detour

#define retracing_ineq
#define precedence_ineq


#define multiple_para 1
//#define ceshiUB
//#define ceshi_UB 756



typedef IloArray<IloNumArray>    NumMatrix;
typedef IloArray<IloBoolArray>    BoolMatrix;
typedef IloArray<IloArray<IloBoolArray> > BoolMatrix2;
typedef IloArray<IloArray<IloNumArray> > NumMatrix2;
typedef IloArray<IloArray<IloArray<IloBoolArray> > > BoolMatrix3;

typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<IloBoolVarArray>    BoolVarMatrix;
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


bool CBMP_1(IloModel model, IloNum nbs, IloIntArray  nbb, IloNumArray nbQ, IloIntArray  nbLocation, IloBoolArray  *nbprecR, IloIntArray  nbreadyT,
	BoolVarMatrix xF, BoolVarMatrix yF, BoolVarMatrix zF, IloIntVar CF, IloNumVarArray QC_CF, BoolVarMatrix uF, IloIntVarArray thetaF, IloBoolVarArray vF,
	IloNum UB);// chun y

bool CSP(IloModel model, IloNum nbs, IloBoolArray  nbb, IloNumArray nbQ, IloBoolArray  *nbprecR, IloIntArray  nbLocation, IloIntArray  nbreadyT,
	BoolMatrix2  xF_best, BoolMatrix yF_best, IloInt *CF_best, BoolMatrix xF_begin_best, BoolMatrix xF_end_best,
	IloNum *ObjVal);//验证是否对

int _tmain(int argc, char* argv[], char* envp[])
{
	double GapAve[instance_no - start_instance + 1];
	double GapAve2[instance_no - start_instance + 1];
	double DuraAve[instance_no - start_instance + 1];
	double ObjAve[instance_no - start_instance + 1];
	double IterAve[instance_no - start_instance + 1];
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

		//cout << filename1 << endl;




		clock_t start = 0, finish = 0;
		clock_t start_MP = 0, finish_MP = 0;
		clock_t start_SP = 0, finish_SP = 0;
		//start = clock();
		double  duration;
		double  duration_MP = 0;
		double  duration_SP = 0;

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
		try 
		{

			//	定义原问题模型
			IloModel model(env);

			IloNum gap;

			//	定义临时变量
			IloInt i,k,j,t,kk; 

			//	定义连续决策变量CF
			IloIntVar   CF(env,0,2000);

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

			BoolVarMatrix uF(env, nbTask);
			for (i = 0; i < nbTask; i++)
			{				
				uF[i]=IloBoolVarArray(env,nbTask);
			}

			IloBoolVarArray vF(env, nbTask);
			IloIntVarArray thetaF(env, nbCrane, 0, safe_margin*nbBay);

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
			//**********************************//
			//            新变量                //
			//**********************************//

			BoolVarMatrix vbF(env, nbCrane);//如果k在i之后要detour
			for (k = 0; k < nbCrane; k++)
			{
				vbF[k] = IloBoolVarArray(env, nbTask);
			}

			BoolVarMatrix ubF(env, nbCrane); //如果i是k的detour的第一个任务
			for (k = 0; k < nbCrane; k++)
			{
				ubF[k] = IloBoolVarArray(env, nbTask);
			}

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
			BoolMatrix FyF_best(env, nbCrane);
			for (k = 0; k < nbCrane; k++)
			{
				FyF_best[k] = IloBoolArray(env, nbTask);
			}

			BoolMatrix CxF_best(env, nbCrane);
			for (k = 0; k < nbCrane; k++)
			{
				CxF_best[k] = IloBoolArray(env, nbTask);
			}
			BoolMatrix CyF_best(env, nbCrane);
			for (k = 0; k < nbCrane; k++)
			{
				CyF_best[k] = IloBoolArray(env, nbTask);
			}


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
			BoolMatrix uF_best(env, nbTask);
			for (i = 0; i <nbTask; i++)
			{				
				uF_best[i]=IloBoolArray(env,nbTask);
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
			IloIntArray  nbreadyT2(env, nbCrane);

			//	定义任务所在贝位参数nbQ
			IloIntArray  nbLocation(env, nbTask);

			IloBoolArray  nbprecR[nbTask];//优先级关系
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

			////展示task位置
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
			//cout<<endl;

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




			//************************************************************************************//
			//      至此已知数据输入完毕                                                          //
			//************************************************************************************//


			//*******************************************************//
			//                 建立主问题 CBMP                       //
			//*******************************************************// 
			int retracingMod=1;

			IloIntArray nbb2(env, nbCrane);// reverse the numbering of the QCs' initial bays

			switch (ModNo)
			{
			case 1: CBMP_1(model, nbs, nbb, nbQ, nbLocation, nbprecR, nbreadyT,
				xF, yF, zF, CF, QC_CF, uF, thetaF, vF, UB);
				break;

			case 2: 
				
				//re-numbering the bays
				for (j = 0; j < nbTask; j++)
				{
					nbLocation[j] = nbBay + 1 - nbLocation[j];
				}
				for (k = 0; k < nbCrane; k++)
					nbb2[nbCrane - 1 - k] = nbBay + 1 - nbb[k];
				for (k = 0; k < nbCrane; k++)
				{
					nbb[k] = nbb2[k];
					//cout << k << " QC: " << nbb[k] << endl;
				}
				for (k = 0; k < nbCrane; k++)
					nbreadyT2[nbCrane - 1 - k] = nbreadyT[k];
				for (k = 0; k < nbCrane; k++)
				{
					nbreadyT[k] = nbreadyT2[k];
					//cout << k << " QC: " << nbb[k] << endl;
				}

				//solve reverse model
				CBMP_1(model, nbs, nbb, nbQ, nbLocation, nbprecR, nbreadyT,
				xF, yF, zF, CF, QC_CF, uF, thetaF, vF, UB);
				retracingMod = 2;
				break;
			};

			IloCplex cplex(env);
			cplex.extract(model);
#ifdef UserActive
			cplex.use(CandyUserCallback(env));
			//cplex.setParam(IloCplex::MIPEmphasis, 2);
			//cplex.setParam(IloCplex::HeurFreq, 1);
			//cplex.setParam(IloCplex::ParallelMode, 1);
			cplex.setParam(IloCplex::Threads, 4);
#endif



			//cplex.setParam(IloCplex::MIPEmphasis, CPX_MIPEMPHASIS_FEASIBILITY);
			//cplex.setParam(IloCplex::EpGap, 0.0006);
			//cplex.setOut(env.getNullStream());
			cplex.setParam(IloCplex::TiLim, 7200);
			//cplex.setParam(IloCplex::Threads, 4);

			//cplex.setParam(IloCplex::VarSel, 3);
			//cplex.setParam(IloCplex::MIPEmphasis, 3);
			//cplex.setParam(IloCplex::Probe, 3);
			
	//cout<<endl<<"number of binary variables = "<<cplex.getNbinVars()<<endl;
	//cout<<endl<<"number of constraints = "<<cplex.getNrows()<<endl;

			cplex.exportModel("model.lp");

			//**********************************//
			//             Benders 主算法           //
			//**********************************//	

			gap = 100;
			int sumiter = 0;//统计迭代数目
			int countNO = 0;

			start = clock();

			int update_no = 0;// if bidirectional solution is updated
			int sp_no;//which sp is the optimal, 0 or 1

			LB = 0;
			UB = 100000;

			//while (UB - LB > 0.1 && (double)(finish - start) / CLOCKS_PER_SEC<600)
			//while ((double)(finish - start) / CLOCKS_PER_SEC <60)
			//while (gap>0.0008)
			//{

				//**********************************//
				//            开始求解		        //
				//**********************************//
				bool h1;

				//start_MP = clock();				
				
				h1 = cplex.solve();

				//finish_MP = clock();

				finish = clock();

				//cout<<"h0 "<<h0<<endl;
				if (!h1)
				{
					//LB = cplex.getBestObjValue();
					//UB = cplex.getValue(CF);

					cout << "\nno feasible solution has been found，algorithm terminate" << endl;
					fout << "\nno feasible solution has been found，algorithm terminate" << endl;
					//cplex.clearModel();
					//cplex.clear();
					//cplex.end();
					//model.end();

					//sumiter++;
					//goto out_end;//此LB下没有可行解，因此算法结束

					//break;

					//return false;
				}
				//**********************************//
				//             记录CBMP最好解       //
				//**********************************//
				//ObjVal = cplex.getValue(CF);
				//LB=cplex.getValue(CF);
				LB = cplex.getBestObjValue();
				UB = cplex.getValue(CF);

				//if (sumiter == 0) LB0 = LB;

				GapAve2[my - 1] = cplex.getMIPRelativeGap();


				fout << "xF_best" << endl;
				for (k = 0; k < nbCrane; k++)
				{
					fout << k << ":  ";
					for (i = 0; i < nbTask; i++)
					{

						if (cplex.isExtracted(xF[k][i]))
						{
							if (cplex.getValue(xF[k][i])>0.1)
								xF_best[k][i] = 1;
							else
								xF_best[k][i] = 0;
						}
						else
							xF_best[k][i] = 0;

						if (xF_best[k][i] == 1) fout << i << "  ";

						if (cplex.isExtracted(yF[k][i]))
						{
							if (cplex.getValue(yF[k][i])>0.1)
								yF_best[k][i] = 1;
							else
								yF_best[k][i] = 0;
						}
						else
							yF_best[k][i] = 0;

						if (yF_best[k][i] == 1) fout << i << "r  ";

					}
					fout << "    "; cout << endl;
				}
				fout << endl << endl;


				//////转化解
				//for (k = 0; k < nbCrane; k++)
				//{
				//	int current_i = 1080;
				//	for (i = 0; i < nbTask; i++)
				//	{
				//		//第一轮
				//		if (xF_best[k][i] == 1)
				//		{
				//			if (current_i == 1080)
				//				xF_begin_best[k][i] = 1;
				//			else
				//				FxF_best[k][current_i][i] = 1;

				//			if (current_i < nbTask && FxF_best[k][current_i][i] == 1)
				//			{
				//				//cout << " QC " << k << " move from " << current_i << " to " << i << endl;
				//				fout << " QC " << k << " move from " << current_i << " to " << i << endl;
				//			}

				//			current_i = i;
				//		}
				//	}
				//	//retracing
				//	for (j = nbBay; j >= 0; j--)
				//	{
				//		for (i = 0; i < nbTask; i++)
				//			if ((nbLocation[i] - 1) == j)
				//			{
				//			if (yF_best[k][i] == 1)
				//			{
				//				FxF_best[k][current_i][i] = 1;
				//				//cout << " QC " << k << " move from " << current_i << " to " << i << endl;
				//				fout << " QC " << k << " move from " << current_i << " to " << i << endl;
				//				current_i = i;
				//			}
				//			}
				//	}
				//	xF_end_best[k][current_i] = 1;

				//}


			
				cplex.clearModel();
				cplex.clear();
				cplex.end();

				//**********************************//
				//             计算下界           //
				//**********************************//	
				//cout << "update_no: " << update_no << endl;
				//finish = clock();


//out_end:	

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





			//UB = LB;

			fout.precision(10);
			fout << endl;

			fout << "\LB = " << multiple_para*LB << endl;
			fout << "\OPT = " << multiple_para*UB << endl;


			fout << "gap=" << GapAve2[my - 1] << endl;

			
			duration = (double)(finish - start) / CLOCKS_PER_SEC;
			fout<<"Time="<<duration<<endl;
			cout<<"Time="<<duration<<endl;

			DuraAve[my - start_instance] = duration;
			ObjAve[my - start_instance] = UB*multiple_para;
			//IterAve[my - 1] = sumiter;

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
	double gapa = 0;
	double gapa2 = 0;
	double durationa = 0;
	double obja = 0;
	double itera = 0;
	for (int i = start_instance; i <= instance_no; i++)
	{
		gapa += GapAve[i - start_instance];
	}
	for (int i = start_instance; i <= instance_no; i++)
	{
		gapa2 += GapAve2[i - start_instance];
	}

	for (int i = start_instance; i <= instance_no; i++)
	{
		durationa += DuraAve[i - start_instance];
	}

	for (int i = start_instance; i <= instance_no; i++)
	{
		obja += ObjAve[i - start_instance];
	}
	for (int i = start_instance; i <= instance_no; i++)
	{
		itera += IterAve[i - start_instance];
	}


	durationa = durationa / (instance_no - start_instance + 1);
	gapa = gapa / (instance_no - start_instance + 1);
	gapa2 = gapa2 / (instance_no - start_instance + 1);
	obja = obja / (instance_no - start_instance + 1);
	cout<<"average gap: "<<gapa<<endl;
	cout << "average gap: " << gapa2 << endl;
	cout<<"average CPU: "<<durationa<<endl;
	cout << "average obj: " << obja << endl;
	cout << "average iter: " << itera << endl;

	ofstream oFile, oFile1;
	string thefilename1 = "haha.csv";

	for (int i = start_instance; i <= instance_no; i++) {
		oFile1.open(thefilename1, ios::app);//| ios::trunc // 这样就很容易的输出一个需要的excel 文件  

		oFile1 << i - start_instance + 1 << "," << ObjAve[i - start_instance] << "," << fixed << setprecision(2) << DuraAve[i - start_instance] << endl;

		oFile1.close();
	}

	system("pause");
	return 0;
}

bool CBMP_1(IloModel model, IloNum nbs, IloIntArray  nbb, IloNumArray nbQ, IloIntArray  nbLocation, IloBoolArray  *nbprecR, IloIntArray  nbreadyT,
	BoolVarMatrix xF, BoolVarMatrix yF, BoolVarMatrix zF, IloIntVar CF, IloNumVarArray QC_CF, BoolVarMatrix uF, IloIntVarArray thetaF, IloBoolVarArray vF,
	IloNum UB)
{

	//xF: forward trip, yF: retracing trip


	IloEnv env = model.getEnv();
	IloInt i, j, k;

	//问题模型
	IloNumVarArray t0kF(env, nbCrane, 0, 100);
	IloNumVarArray gammaF(env, nbCrane, 0, 100);//retracing time
	IloNumVarArray thetaF2(env, nbCrane, 0, 100);


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
		TF[k] = IloNumVarArray(env, nbBay, -100, IloInfinity);
	NumVarMatrix CTF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
		CTF[k] = IloNumVarArray(env, nbBay, -100, IloInfinity);
	

	//**********************************//
	//            原问题目标函数        //
	//**********************************//

	//IloExpr obj1(env); 


	IloExpr  obj2(env);
	//  建立子问题目标函数表达式 

	obj2 += CF;

	//for (i = 0; i < nbTask - 1; i++)
	//	for (j = i + 1; j < nbTask; j++)
	//		if (nbprecR[i][j] == 1)
	//			obj2 += uF[i][j];

	//	将目标函数加入到原问题模型
	model.add(IloMinimize(env, obj2));//
	//obj1.end();
	obj2.end();


	//**********************************//
	//            MP问题 约束           //
	//**********************************//
	//model.add(vF[nbCrane - 1] - vF[nbCrane - 2] <= 0);

	////////////////测试////////////////
	//for (k = 0; k < nbCrane; k++)
	//	for (i = 5*(k+1); i < nbBay; i++)
	//		if (i < nbBay)
	//		{
	//			model.add(zF[k][i] == 0);
	//			for (j = 0; j < nbTask; j++)
	//				if (nbLocation[j]==i+1)
	//					model.add(xF[k][j] + yF[k][j] == 0);
	//		}
				

	////////////////问题约束////////////////
	//13c
	IloRangeArray  c2(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		epa += CF - t0kF[k] - gammaF[k]-nbreadyT[k]-thetaF2[k];
		for (i = 0; i < nbTask; i++)
			epa -= (nbQ[i] / nbs)*(xF[k][i] + yF[k][i]);
		for (j = 0; j < nbBay; j++)
		{
			//epa -= QCmove_time*j*endCzF[k][j];
			epa -= wF[k][j] + CwF[k][j];
		}

		//for (j = 0; j < nbBay; j++)
		//	epa += QCmove_time*j*zF[k][j];

		c2.add(epa >= 0);
		epa.end();
	}
	model.add(c2);
	c2.end();

	//约束 13d // 所有任务被分配到QC上
	IloRangeArray  c4(env);
	for (i = 0; i < nbTask; i++)
	{
		IloExpr  epa(env);
		for (k = 0; k < nbCrane; k++)
			epa += xF[k][i] + yF[k][i];
		c4.add(epa == 1);
		epa.end();
	}
	model.add(c4);
	c4.end();

	//13e
	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbTask; i++)
			model.add(yF[k][i] - vF[k] <= 0);


	////13f 13g
	for (k = 0; k < nbCrane; k++)
		for (j = 0; j < nbTask; j++)
			if (nbLocation[j] < nbBay)
			{
		IloExpr  epa(env);

		for (i = nbLocation[j]; i < nbBay; i++)
			epa += zF[k][i];
		epa += xF[k][j];

		model.add(epa <= 1);
		epa.end();

		IloExpr  epa2(env);

		for (i = nbLocation[j]; i < nbBay; i++)
			epa2 += CzF[k][i];
		epa2 += yF[k][j];

		model.add(epa2 <= 1);
		epa2.end();

		//IloExpr  epa3(env);
		//for (i = 0; i < nbLocation[j] - 1; i++)
		//	epa3 += endCzF[k][i];
		//epa3 += xF[k][j];

		//model.add(epa3 <= 1);
		//epa3.end();

			}



	//13h, 13i,ending bay
	for (k = 0; k < nbCrane; k++)
		for (j = 0; j < nbTask; j++)
		{			
			//model.add(thetaF[k] - (nbLocation[j] - 1) * xF[k][j] >= 0);
			IloExpr  epa(env); 
			for (i = 0; i < nbBay; i++)
				epa += i*endCzF[k][i];
			model.add(epa - (nbLocation[j] - 1) * (xF[k][j] + yF[k][j]) >= 0);
			model.add(epa - nbLocation[j] * yF[k][j] >= 0);
			epa.end();
		}


	//约束 13j,13k : zF unique
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
	model.add(c7);
	c7.end();

	//13l
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa1(env);
		IloExpr  epa2(env);
		epa1 += t0kF[k];
		for (i = 0; i < nbBay; i++)
			epa2 += QCmove_time*i*zF[k][i];
		model.add(epa1 + epa2 - QCmove_time*nbb[k] + QCmove_time >= 0);
		model.add(epa1 - epa2 + QCmove_time*nbb[k] - QCmove_time >= 0);
		epa1.end();
		epa2.end();
	}


	//13m
	IloRangeArray  c13m(env);
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		epa += thetaF2[k];
		for (j = 0; j < nbBay; j++)
		{
			epa -= QCmove_time*j*endCzF[k][j] - QCmove_time*j*zF[k][j];
		}
		c13m.add(epa >= 0);
		epa.end();
	}
	model.add(c13m);
	c13m.end();

	////13n : back time gammaF
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		epa += gammaF[k];
		for (i = 0; i < nbBay; i++)
			epa += QCmove_time*i*CzF[k][i] - QCmove_time*i*endCzF[k][i];
		epa += QCmove_time*nbBay - QCmove_time*nbBay*vF[k];
		model.add(epa >= 0);
		epa.end();
	}

	//13o
	for (k = 0; k < nbCrane; k++)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbBay; i++)
			epa += i*zF[k][i] - i*endCzF[k][i];
		model.add(epa <= 0);
		epa.end();
	}

	//	建立约束 13p,13q:  z 和 z之间隔 /delta +1
	IloRangeArray  c8(env);
	for (k = 0; k < nbCrane - 1; k++)
		if ((1 + safe_margin)*(k+1) + 1 <= nbBay)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbBay; i++)
			epa += i*zF[k + 1][i] - i*zF[k][i];
		c8.add(epa >= 1 + safe_margin);
		epa.end();

		IloExpr  epa2(env);
		for (i = 0; i < nbBay; i++)
			epa2 += i*CzF[k + 1][i] - i*CzF[k][i];
		epa2 += nbBay * (1 - vF[k + 1]);
		c8.add(epa2 >= 1 + safe_margin);
		epa2.end();

	}
	model.add(c8);
	c8.end();

	//13r
	for (k = 0; k < nbCrane - 1; k++)
		if ((1 + safe_margin)*(k+2)  <= nbBay)
	{
		IloExpr  epa(env);
		for (i = 0; i < nbBay; i++)
			epa += i*endCzF[k+1][i] - i*endCzF[k][i];
		model.add(epa >= 1+safe_margin);
		epa.end();
	}

	for (k = 0; k < nbCrane; k++)
	{ 
		if ((1 + safe_margin)*k + 1 > nbBay)
		{
			for (i = 0; i < nbTask; i++)
				model.add(xF[k][i] + yF[k][i] == 0);
		}

		//else if ((1 + safe_margin)*(k+1) + 1 > nbBay && k+1< nbCrane)
		//{
		//	for (i = 0; i < nbTask; i++)
		//		if (nbLocation[i]<=(1 + safe_margin)*k)
		//			model.add(xF[k+1][i] + yF[k+1][i] == 0);
		//}
	}

	//////1+delta 版本有错误
	////	建立约束 13p,13q:  z 和 z之间隔 /delta +1
	//IloRangeArray  c8(env);
	//for (k = 0; k < nbCrane - 1; k++)
	//{
	//	IloExpr  epa(env);
	//	for (i = 0; i < nbBay; i++)
	//		epa += i*zF[k + 1][i] - i*zF[k][i];
	//	c8.add(epa >= 1 );
	//	epa.end();

	//	IloExpr  epa2(env);
	//	for (i = 0; i < nbBay; i++)
	//		epa2 += i*CzF[k + 1][i] - i*CzF[k][i];
	//	epa2 += nbBay * (1 - vF[k + 1]);
	//	c8.add(epa2 >= 1);
	//	epa2.end();

	//}
	//model.add(c8);
	//c8.end();

	////13r
	//for (k = 0; k < nbCrane - 1; k++)
	//	{
	//	IloExpr  epa(env);
	//	for (i = 0; i < nbBay; i++)
	//		epa += i*endCzF[k + 1][i] - i*endCzF[k][i];
	//	model.add(epa >= 1 );
	//	epa.end();
	//	}
	



	//////约束（5）//endCzF 取在最后一个xF处之后
	//////约束（6）// zF 小于最小的xF的bay
	//for (j = 0; j < nbTask; j++)
	//	for (k = 0; k < nbCrane; k++)
	//	{
	//	IloExpr  epa(env);

	//	epa += yF[k][j];
	//	for (i = 0; i < nbLocation[j]; i++)
	//		epa -= CzF[k][i];

	//	model.add(epa <= 0);
	//	epa.end();

	//	IloExpr  epa2(env);
	//	epa2 += xF[k][j];
	//	for (i = 0; i < nbLocation[j]; i++)
	//		epa2 -= zF[k][i];
	//	model.add(epa2 <= 0);
	//	epa2.end();

	//	IloExpr  epa3(env);
	//	epa3 += xF[k][j];
	//	for (i = nbLocation[j] - 1; i < nbBay; i++)
	//		epa3 -= endCzF[k][i];
	//	model.add(epa3 <= 0);
	//	epa3.end();

	//	}


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

	//	if (k < nbCrane)
	//		for (j = nbBay - 2 * (nbCrane - k - 1); j < nbBay; j++)
	//			model.add(endCzF[k][j] == 0);
	//	if (k > 0)
	//		for (j = 0; j < 2 * k; j++)
	//			model.add(zF[k][j]+CzF[k][j] == 0);

	//}
	//model.add(v00);
	//v00.end();



////////////////sub-problem

	//// vF 取值(2)
	//for (k = 0; k < nbCrane; k++)
	//{
	//	IloExpr  epa(env);
	//	for (i = 0; i < nbTask; i++)
	//			epa += yF[k][i];
	//	epa -= vF[k];
	//	model.add(epa >= 0);
	//	epa.end();
	//}


	//for (k = 0; k < nbCrane; k++)
	//{
	//	IloExpr  epa(env);
	//	for (i = 0; i < nbTask; i++)
	//		epa += xF[k][i];
	//	model.add(epa >= 1);
	//	epa.end();
	//}

	//// vF 取值(3) vF与violation关系
	//for (i = 0; i < nbTask; i++)
	//	for (j = 0; j < nbTask; j++)
	//		if (nbprecR[i][j] == 1)
	//		{ 
	//		for (k = 0; k < nbCrane; k++)
	//		{
	//			IloExpr epa(env);
	//			IloExpr epa2(env);
	//			
	//			epa += vF[k] +1;
	//			epa2 += vF[k] + 1;

	//			if (k < nbCrane - 1)
	//			{
	//				epa -= xF[k][i] + yF[k][i];
	//				for (int kk = k+1; kk < nbCrane; kk++)
	//					epa -= xF[kk][j] + yF[kk][j];
	//				model.add(epa >= 0);
	//			}
	//			
	//			if (k > 0)
	//			{
	//				epa2 -= xF[k][j] + yF[k][j];
	//				for (int kk = 0; kk < k; kk++)
	//					epa2 -= xF[kk][i]+yF[kk][i];
	//				model.add(epa2 >= 0);
	//			}
	//			epa.end();
	//			epa2.end();
	//		}
	//		}


	//precedence
for (i = 0; i < nbTask; i++)
	for (j = 0; j < nbTask; j++)
		if (nbprecR[i][j] == 1)
		{
		for (k = 0; k < nbCrane; k++)
		{

			{
				IloExpr epa(env);
				epa += yF[k][i];
				for (int kk = k; kk < nbCrane; kk++)
					epa -=  yF[kk][j];
				model.add(epa <= 0);
				epa.end();
			}

			{
				IloExpr epa2(env);
				epa2 += xF[k][j];
				for (int kk = k; kk < nbCrane; kk++)
					epa2 -= xF[kk][i];
				model.add(epa2 <= 0);
				epa2.end();
			}

			if (k < nbCrane - 1)
			{
				IloExpr epa3(env);
				epa3 += xF[k][i];
				for (int kk = k + 1; kk < nbCrane; kk++)
					epa3 += xF[kk][j];
					//epa3 += xF[kk][j] + yF[kk][j];
				model.add(epa3 <= 1);
				epa3.end();
			}

			if (k < nbCrane - 1)
			{
				/*for (int kk = k + 1; kk < nbCrane; kk++)
				{
				IloExpr epa3(env);
				IloExpr epa4(env);
				epa3 += xF[k][i] + xF[kk][j];
				epa4 += yF[k][j] + yF[kk][i];
				model.add(epa3 <= 1);
				model.add(epa4 <= 1);
				epa3.end();
				epa4.end();
				}*/
				IloExpr epa3(env);
				IloExpr epa4(env);
				epa3 += xF[k][i];
				epa4 += yF[k][j];
				for (int kk = k + 1; kk < nbCrane; kk++)
				{
					epa3 += xF[kk][j];
					epa4 += yF[kk][i];
				}
				model.add(epa3 <= 1);
				model.add(epa4 <= 1);
				epa3.end();
				epa4.end();

			}

			//if (k < nbCrane - 1)
			//{
			//	IloExpr epa4(env);
			//	epa4 += yF[k][j];

			//	//for (int kk = k + 1; kk < nbCrane; kk++)
			//	//	epa4 += yF[kk][i];
			//	//epa4 -= 1;
			//	

			//	for (int kk = k + 1; kk < nbCrane; kk++)
			//		epa4 -= xF[kk][i];
			//	for (int kk = 0; kk <= k; kk++)
			//		epa4 -= xF[kk][i] + yF[kk][i];
			//	model.add(epa4 <= 0);
			//	epa4.end();
			//}

			//for (int kk = 0; kk <= k; kk++)
			//{
			//	model.add(yF[k][j] + xF[kk][i]-vF[kk]<=1);
			//}
			
		}
	}



//////////////////////////////////////
	//////waiting time ///////
/////////////////////////////////////

#ifdef waiting_constraints
	///////forward trip
	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbBay; i++)
		{
			IloExpr  epa(env);
			epa += t0kF[k] + nbreadyT[k];
			for (j = 0; j < nbTask; j++)
				if (nbLocation[j] - 1 <= i)
					epa += (nbQ[j] / nbs)*xF[k][j];
			for (j = 0; j <= i; j++)
			{

				epa += wF[k][j];
				epa -= QCmove_time*j*zF[k][j];
				//epa -= 1000 * CzF[k][j];
			}
			model.add(epa + QCmove_time*i - TF[k][i] == 0);
			epa.end();

		}

	/////约束（3）
	IloRangeArray  c3(env);
	for (i = 0; i < nbBay - 1 - safe_margin; i++)
		for (k = 0; k < nbCrane - 1; k++)
		{
			IloExpr  epa(env);			

			epa += TF[k][i] - TF[k + 1][i + 1 + safe_margin];
			epa += 2000;
			for (j = 0; j <= i + 1+safe_margin; j++)
				epa -= 2000 * zF[k + 1][j];
			for (j = 0; j < i; j++)
				epa += 2000 * endCzF[k][j];

			c3.add(epa >= 0);
			epa.end();
		}
	model.add(c3);
	c3.end();

	//////retracing trip
	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbBay; i++)
		{
			IloExpr  epa(env);
			epa += t0kF[k] + nbreadyT[k];
			for (j = 0; j < nbTask; j++)
				epa += (nbQ[j] / nbs)*xF[k][j];
			for (j = 0; j < nbBay; j++)
			{
				epa += 2 * QCmove_time*j*endCzF[k][j] - QCmove_time*j*zF[k][j];
				epa += wF[k][j];
			}
			

			for (j = 0; j < nbTask; j++)
				if (nbLocation[j] - 1 >= i)
					epa += (nbQ[j] / nbs)*yF[k][j];
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
			epa += CTF[k + 1][i + 1+safe_margin] - CTF[k][i];

			epa += 6000 - 2000 * vF[k] - 2000 * vF[k + 1];
			for (j = 0; j < i; j++)
				epa += 2000 * endCzF[k][j];
			for (j = 0; j <= i+1+safe_margin; j++)
				epa -= 2000 * CzF[k+1][j];

			model.add(epa>=0);
			epa.end();
		}

	// forward and retrace
	for (k = 0; k < nbCrane - 1; k++)
	{
		IloExpr  epa(env);
		for (j = 0; j <nbBay; j++)
			epa += j * CzF[k + 1][j] - j * endCzF[k][j];
		epa += 2000 + 2000 * vF[k] - 2000 * vF[k + 1];
		model.add(epa >= 1 + safe_margin);
		epa.end();

	}

#endif
	return true;

}

bool CSP(IloModel model, IloNum nbs, IloBoolArray  nbb, IloNumArray nbQ, IloBoolArray  *nbprecR, IloIntArray  nbLocation, IloIntArray  nbreadyT,
	BoolMatrix2  xF_best, BoolMatrix yF_best, IloInt *CF_best, BoolMatrix xF_begin_best, BoolMatrix xF_end_best,
	IloNum *ObjVal)
{
	IloEnv env = model.getEnv();
	IloInt i, j, k;
	IloInt nbM = 2000;

	//	定义连续决策变量CF
	IloIntVar   CF(env, 0, 1800);

	//IloNumArray tem_C(env, nbCrane);// 子问题计算的completion time
	IloNumVarArray QC_CF(env, nbCrane, 0, IloInfinity);// each QC's completion time
	IloNumVarArray Task_CF(env, nbTask, 0, IloInfinity);// each task's completion time

	//	定义决策变量xF,CyF,CzF
	BoolVarMatrix yF(env, nbCrane);
	for (k = 0; k < nbCrane; k++)
	{
		yF[k] = IloBoolVarArray(env, nbTask);
	}

	BoolVarMatrix zF(env, nbTask);
	for (i = 0; i < nbTask; i++)
	{
		zF[i] = IloBoolVarArray(env, nbTask);
	}

	//问题模型
	IloModel submodel(env);

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


	IloNumArray ctime_t[nbTask];// the travelling time of QC moving from the bay i to the bay j
	for (i = 0; i < nbTask; i++)
		ctime_t[i] = IloNumArray(env, nbTask);


	for (i = 0; i < nbTask; i++)
	{
		for (j = 0; j < nbTask; j++)
		{
			if (nbLocation[i] == nbLocation[j]) ctime_t[i][j] = 0;
			else if (nbLocation[i] > nbLocation[j])
			{
				ctime_t[i][j] = (nbLocation[i] - nbLocation[j])*crane_move_time;
			}
			else
			{
				ctime_t[i][j] = (nbLocation[j] - nbLocation[i])*crane_move_time;
			}
		}

	}

	//**********************************//
	//            子问题 约束           //
	//**********************************//


	IloRangeArray  c2(env);
	for (k = 0; k < nbCrane; k++)
	{

		IloExpr  epa(env);
		epa += CF;
		epa -= QC_CF[k];

		c2.add(epa >= 0);
		epa.end();

	}
	submodel.add(c2);
	c2.end();

	IloRangeArray  c7(env);
	for (k = 0; k < nbCrane; k++)
	{
		for (i = 0; i < nbTask; i++)
		{
			for (j = 0; j < nbTask; j++)
				if (xF_best[k][i][j] == 1)
				{
				IloExpr  epa(env);
				epa += Task_CF[i] + ctime_t[i][j] + nbQ[j] - Task_CF[j];

				c7.add(epa <= 0);
				epa.end();
				}
		}

	}
	submodel.add(c7);
	c7.end();

	IloRangeArray  c8(env);
	for (i = 0; i < nbTask; i++)
	{
		for (j = 0; j < nbTask; j++)
		{
			if (nbLocation[i] != nbLocation[j])
			{
				IloExpr  epa(env);
				epa += Task_CF[i] + nbQ[j] - Task_CF[j];
				epa -= nbM*(1 - zF[i][j]);
				c8.add(epa <= 0);
				epa.end();
			}
		}
	}
	submodel.add(c8);
	c8.end();

	IloRangeArray  c9(env);
	for (i = 0; i < nbTask; i++)
	{
		for (j = 0; j < nbTask; j++)
			if (nbLocation[i] != nbLocation[j])
			{
			IloExpr  epa(env);
			epa += Task_CF[j] - nbQ[j] - Task_CF[i];
			epa -= nbM*zF[i][j];

			c9.add(epa <= 0);
			epa.end();
			}
	}

	submodel.add(c9);
	c9.end();

	IloRangeArray  c10(env);
	for (k = 0; k < nbCrane; k++)
		for (i = 0; i < nbTask; i++)
		{
		for (j = 0; j < nbTask; j++)
		{
			if (nbLocation[i]<nbLocation[j])
			{
				IloExpr  epa(env);
				for (int kk = 0; kk <= k; kk++)
					epa += yF_best[kk][j];
				for (int kk = k; kk < nbCrane; kk++)
					epa += yF_best[kk][i];
				epa -= zF[i][j] + zF[j][i];
				c10.add(epa <= 1);
				epa.end();
			}
		}
		}
	submodel.add(c10);
	c10.end();

	IloRangeArray  c11(env);
	for (i = 0; i < nbTask; i++)
	{
		for (j = 0; j < nbTask; j++)
			if (nbLocation[i] == nbLocation[j])
			{
			IloExpr  epa(env);
			epa += Task_CF[i] + nbQ[j] - Task_CF[j];
			epa -= nbM*(1 - zF[i][j]);
			for (k = 0; k < nbCrane; k++)
			{
				epa += crane_move_time*xF_begin_best[k][j];
				for (int ii = 0; ii < nbTask; ii++)
					if (nbLocation[ii] != nbLocation[i])
						epa += crane_move_time*xF_best[k][ii][j];
			}
			c11.add(epa <= 0);
			epa.end();
			}
	}

	submodel.add(c11);
	c11.end();


	IloRangeArray  c12(env);

	for (i = 0; i < nbTask; i++)
	{
		for (j = 0; j < nbTask; j++)
			if (nbLocation[i] == nbLocation[j])
			{
			IloExpr  epa(env);
			epa += Task_CF[j] - nbQ[j] - Task_CF[i];
			epa -= nbM*zF[i][j];
			for (k = 0; k < nbCrane; k++)
			{
				epa -= xF_begin_best[k][j];
				for (int ii = 0; ii < nbTask; ii++)
					if (nbLocation[ii] != nbLocation[i])
						epa -= crane_move_time*xF_best[k][ii][j];
			}
			c12.add(epa <= 0);
			epa.end();
			}
	}

	submodel.add(c12);
	c12.end();


	IloRangeArray  c13(env);
	for (i = 0; i < nbTask; i++)
	{
		for (j = 0; j < nbTask; j++)
			if (i != j)
			{
			//if (nbLocation[i] == nbLocation[j] || nbLocation[i] == nbLocation[j] + 1 || nbLocation[j] == nbLocation[i] + 1)
			if (nbLocation[i] == nbLocation[j] + 1 || nbLocation[j] == nbLocation[i] + 1 || nbLocation[j] == nbLocation[i])
				//if (nbLocation[i] == nbLocation[j]) 
			{
				IloExpr  epa(env);
				epa += zF[i][j] + zF[j][i];

				c13.add(epa == 1);
				epa.end();
			}
			}
	}
	submodel.add(c13);
	c13.end();

	IloRangeArray  c14(env);//i<j, and both located in the same bay, then i should be processed before j
	for (i = 0; i < nbTask; i++)
	{
		for (j = 0; j < nbTask; j++)
		{
			if (nbprecR[i][j] == 1)
			{
				c14.add(zF[i][j] == 1);
				c14.add(zF[j][i] == 0);
			}
		}
	}
	submodel.add(c14);
	c14.end();

	IloRangeArray  c15(env);
	for (k = 0; k < nbCrane; k++)
	{
		for (j = 0; j < nbTask; j++)
			if (xF_begin_best[k][j] == 1)
			{
			IloExpr  epa(env);
			epa += nbQ[j] - Task_CF[j];
			epa += nbreadyT[k];//QC ready time, already defined in the main function
			if (nbb[k] >= nbLocation[j])
				epa += nbb[k] - nbLocation[j];
			else
				epa += nbLocation[j] - nbb[k];

			c15.add(epa <= 0);
			epa.end();
			}

	}
	submodel.add(c15);
	c15.end();

	IloRangeArray  c16(env);
	for (k = 0; k < nbCrane; k++)
	{
		for (j = 0; j < nbTask; j++)
			if (xF_end_best[k][j] == 1)
			{
			IloExpr  epa(env);
			epa += Task_CF[j] - QC_CF[k];

			c16.add(epa <= 0);
			epa.end();
			}

	}
	submodel.add(c16);
	c16.end();




	IloRangeArray  c0(env);
	for (k = 0; k < nbCrane; k++)
	{
		for (i = 0; i < nbTask; i++)
		{
			for (j = 0; j < nbTask; j++)
				if (xF_best[k][i][j] == 1)
					c0.add(zF[i][j] == 1);
		}

	}
	submodel.add(c0);
	c0.end();


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



	bool h0;
	h0 = cplex.solve();
	//cout<<"h0 "<<h0<<endl;
	if (!h0)
	{
		cout << "\nno feasible solution has been found" << endl; 
		cplex.clearModel();
		cplex.clear();
		cplex.end();
		submodel.end();

		return false;
	}
	//**********************************//
	//             记录最好解           //
	//**********************************//

	*ObjVal = cplex.getBestObjValue();




	*CF_best = cplex.getValue(CF);


	cout << "QC_CF_best" << endl;
	for (k = 0; k < nbCrane; k++)
	{
		cout << cplex.getValue(QC_CF[k]) << "  ";
	}
	cout << endl << endl;

	cout << "Task_CF_best" << endl;
	for (i = 0; i < nbTask; i++)
	{
		cout << i << ": " <<cplex.getValue(Task_CF[i]) << ";  ";
	}
	cout << endl << endl;

	////cout << "CF_best: " << *CF_best << endl;
	//cout << "ObjVal: " << *ObjVal << endl;

	cplex.clearModel();
	cplex.clear();
	cplex.end();
	submodel.end();

	return true;

}

