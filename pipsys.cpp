#include <algorithm>
#include<iostream>
#include<QDebug>
#include<QMessageBox>

#include "pipsys.h"
#include"my_math.h"

using namespace std;

Pipsys::Pipsys(QWidget *parrent):parrent(parrent)
{

    //Creat elements name list
     elementsList<<"reservoir"<<"pipe"<<"node"<<"branch";
	//计算时间步数
	steps = int(totalTime / timeStep + 0.5) + 1; //int向下取整，+0.5四舍五入
	//时间序列
	timeSeries = regspace<vec>(0, timeStep, totalTime);
	//读取管道

	//初始化节点
	//读取边界条件
	for_each(reservoirs.begin(), reservoirs.end(), [this](Reservoirs& r)
	{
		r.Connect(nodes);
		if (r.type == RE::Upper)
		{
			nUpperReservoirs++; // 上游水库数
		}
		else
			nLowerReservoirs++; // 下游水库数
	});
	//水轮机

	for_each(turbines.begin(), turbines.end(), [this](Turbines& tur)
	{
		tur.GetGuide(timeSeries);
		tur.IniGuideLine();
		tur.Connect(nodes);
		iniQ = tur.Q0;
		tur.nTPoint = steps;
		tur.printDensity = printDensity;
	});
}

void Pipsys::addReservoir(){
    nReservoirs++;
reservoirs.push_back(Reservoirs(nReservoirs));
}
void Pipsys::addPipe(){
    nPipes++;
pipes.push_back(Pipes(nPipes));
}
void Pipsys::addNode(){
    nNodes++;
nodes.push_back(Nodes(nNodes));
}

void Pipsys::addElement(QString eleName)
{
    if(eleName=="reservoir")
        addReservoir();
    else if(eleName=="pipe")
        addPipe();
    else if(eleName=="node")
        addNode();
    else
        error("Unknown element!");
};

double direction(int n1, int n2, mat const& QS, mat const& S);

void Pqr(int, int, double, double, colvec&, mat&, mat&, mat&);

void Pipsys::Steady()
{
	//恒定流计算管道头尾节点的流量和水头
	const int IMAX = 200;
	int dim = nNodes;
	mat QS(dim, dim, fill::zeros);
	mat SY(dim, dim, fill::zeros);
	mat SJ(dim, dim, fill::zeros);
	mat S(dim, dim, fill::zeros);
	colvec E(dim, fill::zeros);
	//初始化
	SteadyIni(SY, SJ, S, QS, E);
	//E.print();
#pragma region Iteration
	int iter;
	for (iter = 0; iter < IMAX; iter++)
	{
		mat p(dim, dim, fill::zeros);
		mat q(dim, dim, fill::zeros);
		mat r(dim, dim, fill::zeros);
		mat A(dim, dim, fill::zeros);
		colvec B(dim, fill::zeros);
		colvec DE(dim, fill::zeros);
		// 初始化pqr
		for (int i = 0; i < nPipes; i++)
		{
			int je = pipes[i].je;
			int js = pipes[i].js;
			double s = direction(js, je, QS, S);
			double fQ, dfQ;
			//尾节点是水轮机的管道
			if (nodes[je].type == BD::Turbine)
				dfQ = turbines[nodes[je].typeID].fQdfQ(QS(js, je), fQ);
			else
				fQ = 0, dfQ = 0;
			//计算能量损失（流量的函数）fQ及其导数
			fQ = fQ + s * fabs(QS(js, je)) * QS(js, je);
			dfQ = dfQ + 2 * s * fabs(QS(js, je));
			Pqr(js, je, dfQ, fQ, E, p, q, r);
		}
		////计算总体单元矩阵A和B
		////A(i,i) = 求和 pij ，A(i,j)=q(i, j)
		////B(i)=ci- 求和 Q - 求和 r	
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				if (i == j)
				{
					for (int k2 = 0; k2 < dim; k2++)
						A(i, i) = A(i, i) + p(i, k2); // total pij of Ei
				}
				else
					A(i, j) = q(i, j); //qif of Ej
				B[i] = B[i] - QS(i, j) - r(i, j); //
			}
			//水库单元方程 ΔE(i,i)=0
			if (nodes[i].type == BD::UpperReservoir || nodes[i].type == BD::LowerReservoir) //reservoir
			{
				//dEi=0
				for (int j = 0; j < dim; j++)
				{
					if (i == j)
						A(i, j) = 1;
					else
						A(i, j) = 0;
				}
				B[i] = 0;
			}
		}
		//解单元方程组
		DE = solve(A, B, solve_opts::refine);
		//更新流量
		for (const auto& pipe : pipes)
		{
			int js = pipe.js;
			int je = pipe.je;
			QS(js, je) = QS(js, je) + p(js, je) * DE[js] + q(js, je) * DE[je] + r(js, je); //flow of every pipe 
			QS(je, js) = -QS(js, je);
		}
		//更新水头
		//初始化访问变量
		vector<bool> visit(dim, false);
		for (auto r : reservoirs)
		{
			int id = r.nodeID;
			E[id] = reservoirs[nodes[id].typeID].waterLevel;
			visit[id] = true;
			Q2E(nodes[id], E, S, QS, visit);
		}
		//判断结果精度
		if (abs(DE).max() < 0.00001)
			break;
	}
	if (iter == IMAX)
        error("恒定流结果不收敛!");
	cout << "恒定流计算完成" << endl;
	E.print("E");
	//利用恒定流计算结果，计算管道的水头、流量
	for_each(pipes.begin(), pipes.end(), [&](auto& pipe)
	{
		int js = pipe.js;
		int je = pipe.je;
		//管道流量
		pipe.Q[0] = QS(js, je);
		//管道头结点的水头，能量减去速度头
		pipe.H[0] = E(js) - pipe.Q[0] * pipe.Q[0] / 2 / g / pipe.A / pipe.A;
		//从头结点到尾结点，依次计算管道其余节点的水头
		for (int j = 1; j <= pipe.NN; j++)
		{
			double s;
			//管道流量不变
			pipe.Q[j] = pipe.Q[0];
			//判断流量方向
			if (pipe.Q[j] > 0)
				s = SY(js, je);
			else s = SY(je, js);	
			pipe.H[j] = pipe.H[j - 1] - s * fabs(pipe.Q[j]) * pipe.Q[j] / pipe.NN;
		}
	});
	
	//耦合20 s的特征线法，消除恒定流和非恒定流之间的跳跃
	for (double t = 0;t <15;t=t+timeStep)
		MOC();
	//保存恒定流状态
	for_each(pipes.begin(), pipes.end(), [&](auto& pipe)
	{
		pipe.H0 = pipe.H;
		pipe.Q0 = pipe.Q;
	});
}

void Pipsys::SteadyIni(mat& SY, mat& SJ, mat& S, mat& QS, vec& E)
{
	//初始化管道流量和流量系数，dH=sQ^2
	for (auto pipe : pipes)
	{
		int js = pipe.js;
		int je = pipe.je;
		//流量
		QS(js, je) = iniQ;
		QS(je, js) = -QS(js, je);
		//沿程损失流量系数
		SY(js, je) = 8.0 * pipe.f * pipe.length / PI / PI / pow(pipe.D, 5.0) / g; //on-way loss factor;
		SY(je, js) = SY(js, je);
		//局部损失流量系数
		SJ(js, je) = pipe.zeta[0] / 2 / g / pipe.A / pipe.A;
		SJ(je, js) = pipe.zeta[1] / 2 / g / pipe.A / pipe.A;
		//流量系数
		S(js, je) = SY(js, je) + SJ(js, je);
		S(je, js) = SY(je, js) + SJ(je, js);
	}
	//初始化访问变量
	vector<bool> visit(nNodes, false);
	//初始化水头，从水库开始遍历
	for (auto r : reservoirs)
	{
		int id = r.nodeID;
		E[id] = reservoirs[nodes[id].typeID].waterLevel;
		visit[id] = true;
		Q2E(nodes[id], E, S, QS, visit);
	}
}

void Pipsys::Q2E(Nodes& node, mat& E, mat& S, mat& QS, vector<bool>& visit)
{
	////目的：通过流量计算各节点的能量
	////输入：流量，流量系数，访问标记序列
	////输出：各节点的能量
	//当前节点ID
	int n1 = node.id;
	//if (search[n1])
	//	return;
	double DH; //能量损失
	//遍历节点处的管道
	for (auto p2 : node.pID)
	{
		//判断下个节点
		int n2;
		if (pipes[p2].js == node.id)
			n2 = pipes[p2].je;
		else n2 = pipes[p2].js;
		//判断是否已初始化
		if (visit[n2])
			continue;
		int typeID = nodes[n2].typeID;
		//判断流量方向
		double s = direction(n1, n2, QS, S);
		//计算单元能量损失
		switch (nodes[n2].type)
		{
		case BD::UpperReservoir:
			//上游节点
			E(n2) = reservoirs[typeID].waterLevel;
			visit[n2] = true;
			continue;
		case BD::LowerReservoir:
			//下游节点
			E(n2) = reservoirs[typeID].waterLevel;
			visit[n2] = true;
			continue;
		case BD::Turbine:
			//水轮机
			if (p2 == nodes[n2].spID[0])
				//正向,损失计入前管道
				turbines[typeID].fQdfQ(QS(n1, n2), DH);
			else
				return;
			break;
		default:
			DH = 0;
			break;
		}
		E(n2) = E(n1) - s * fabs(QS(n1, n2)) * QS(n1, n2) - DH;
		visit[n2] = true;
		//遍历下一个节点
		Q2E(nodes[n2], E, S, QS, visit);
	}
}

double direction(int n1, int n2, mat const& QS, mat const& S)
{
	////目的：判断流量方向，返回对应方向的流量系数
	if (QS(n1, n2) > 0)
		return S(n1, n2);
	return S(n2, n1);
}

void Pqr(int sta, int end, double dfQ, double fQ, colvec& E, mat& p, mat& q, mat& r)
{
	////输入：单元首末节点，流量函数的导数，流量函数（ΔE=f(Q)），
	////输出：p,q,r
	////目的：计算单元方程的系数p,q,r,为构成单元方程组的系数矩阵A,B做准备
	//dfQ不能为0
	if (dfQ == 0)
		dfQ = 1E-11;
	p(sta, end) = 1 / dfQ;
	q(sta, end) = -p(sta, end);
	r(sta, end) = p(sta, end) * (E(sta) - E(end) - fQ);
	p(end, sta) = p(sta, end);
	q(end, sta) = q(sta, end);
	r(end, sta) = -r(sta, end);
}

void Pipsys::MOC()
{
	////目的：单步非恒定流计算————特征线法
	//管道
	for_each(pipes.begin(), pipes.end(), [](Pipes & pipe) {
		pipe.MOC();});
	//节点
	for_each(nodes.begin(), nodes.end(), [this](Nodes & node) {
		int typeID = node.typeID;
		if (node.type == BD::UpperReservoir || node.type == BD::LowerReservoir)
			reservoirs[typeID].MOC(pipes);
		else if (node.type == BD::Series)
			node.SeriesMOC(pipes);
		else if (node.type == BD::Branch)
			node.BranchMOC(pipes);
		else if (node.type == BD::Turbine)
			turbines[typeID].MOC(pipes, T);
		else
			error("未知节点类型！");
	});
	for_each(pipes.begin(), pipes.end(), [](Pipes & pipe) {
		pipe.H = pipe.HP;
		pipe.Q= pipe.QP;
		//pipe.QP.print("QP=");
		//pipe.HP.print("L=");
		});
}

void Pipsys::Run()
{
    clock_t time_begin = clock();
    qDebug()<<"Calculation start!";
	////目的：开始管网计算，包括恒定流和非恒定流
	//恒定流计算
	if (!isSteady)
	{
		T = 0;
		Steady();
		isSteady = true;
	}
//        qDebug<<"Steady flow calculation has been done successfully!";
    //非恒定流计算c ff
    for (T = 0;T < totalTime+ timeStep;T+= timeStep)
		MOC();
	//输出机组非恒定流信息
	for(int i=0;i<nTurbines;i++)
        turbines[i].SaveTurbine();
//     qDebug<<"Transient calculation has been done successfully!";
     QMessageBox *mocFinish;
     mocFinish=new QMessageBox();
     mocFinish->setWindowTitle("Message");
      mocFinish->setText("Transient calculation has been done successfully!");
     mocFinish->setAttribute(Qt::WA_DeleteOnClose);
     mocFinish->exec();
    //结束计时并打印
    clock_t time_end = clock();
    double time_cost = (double)(time_end - time_begin) / CLOCKS_PER_SEC;
//    QMessageBox::information(nullptr,"Time","time cost is"+QString::number(time_cost));
    qDebug()<<"Cost time: "<<time_cost;
}

