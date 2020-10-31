#include<iostream>
#include"elements.h"
#include"my_math.h"


Pipes::Pipes(int id, int js, int je, double timeStep, double a, double D, double len, double roughness,
             double zeta1,
             double zeta2)
	: id(id), js(js), je(je), timeStep(timeStep), a(a), D(D), length(len), roughness(roughness), zeta{zeta1, zeta2}
{
	//分段
	NN = (int)(length / timeStep / a + 1);
	//计算空间步长
	dx = length / NN;
	//调整波速
	a = length / timeStep / NN;
	//计算管道面积
	A = PI * D * D / 4;
	//计算沿程阻力系数
	/*曼宁公式 C=R^(1/6)， R=D/4
	谢才公式 v=C(RJ)^0.5, J=hf/l
	DW公式   hf=flv^2/(2gD)
			f=8gn^2(D/4)^(1/3)
	*/
	f = 8 * g * roughness * roughness * pow((D / 4), -1.0 / 3.0);

	//状态变量
	H0.resize(NN + 1); //初始水头                  
	Q0.resize(NN + 1); //初始流量
	H.resize(NN + 1); //上一时刻水头                      
	Q.resize(NN + 1); //上一时刻流量
	HP.resize(NN + 1); //当前时刻水头
	QP.resize(NN + 1); //当前时刻流量


	m1 = a / g / A; // a / gA
	m2 = timeStep * f * a / 2.0 / g / D / A / A; // fadt / (2gDA2)
	m3 = g * A * timeStep / dx; // gAdt / dx
	m4 = 1.0 / 2 / g / A / A; //  1/(2gA^2)
}

Pipes::Pipes()
{
};

double Pipes::QCP(int J)
{
	double qcp = ((Q[J - 1]) * m1 + H[J - 1]) / (m1 + m2 * abs(Q[J - 1]));
	return qcp;
};

double Pipes::CQP(int J)
{
	double cqp = 1 / (m1 + m2 * abs(Q[J - 1]));
	return cqp;
};

double Pipes::QCM(int J)
{
	double qcm = (Q[J + 1] * m1 - H[J + 1]) / (m1 + m2 * abs(Q[J + 1]));
	return qcm;
};

double Pipes::CQM(int J)
{
	double cqm = 1 / (m1 + m2 * abs(Q[J + 1]));
	return cqm;
};

void Pipes::MOC()
{
	for (int j = 1;j < NN;j++)
	{
		HP[j] = (QCP(j) - QCM(j)) / (CQP(j) + CQM(j));
		QP[j] = QCP(j) - CQP(j) * HP[j];
	}		
}

Nodes::Nodes(int id, int typeID, BD type) :
	id(id), typeID(typeID), type(type)
{
}

void Nodes::AddPipe(Pipes& pipe, bool start)
{
	//更新管道信息
	if (start)
	{
		// 前管道数+1
		n_sp++;
		// 新增前管道id
		spID.push_back(pipe.id);
	}
	else
	{
		// 后管道数+1
		n_ep++;
		// 新增后管道id
		epID.push_back(pipe.id);
	}
}

void Nodes::SeriesMOC(vector<Pipes>& pipes)
{
	Pipes & sp = pipes[spID[0]];
	Pipes & ep = pipes[epID[0]];
	double a1 =  sp.m4 - ep.m4 + sp.zeta[0] * sp.m4;
	double a2 = sp.m4 - ep.m4 - sp.zeta[1] *  sp.m4;
	int & NN = sp.NN;
	////Init end////
	int a0;
	double b0 = 1 / sp.CQP(NN) + 1 / ep.CQM(0);
	double c0 = -sp.QCP(NN) / sp.CQP(NN) - ep.QCM(0) / ep.CQM(0);
	if (sp.Q[NN] >= 0)
		a0=a1;
	else a0=a2;
	if (fabs(a0) < 0.000001)
		sp.QP[NN] = -c0 / b0;
	else
	{
		double y1 = (-b0 - sqrt(b0*b0 - 4.0*a0 * c0)) / (2.0*a0);
		double y2 = (-b0 + sqrt(b0*b0 - 4.0*a0 * c0)) / (2.0*a0);
		if (fabs(y1 - sp.Q[NN]) < fabs(y2 - sp.Q[NN]))
			sp.QP[NN] = y1;
		else
			sp.QP[NN] = y2;
	}
	ep.QP[0] = sp.QP[NN];
	sp.HP[NN] = (sp.QCP(NN) - sp.QP[NN]) / sp.CQP(NN);
	ep.HP[0] = (ep.QP[0] - ep.QCM(0)) / ep.CQM(0);
}

void Nodes::BranchMOC(vector<Pipes>& L)
{
	double l = 0.0, m = 0.0, B = 0.0;
	//		int n = JP->StartPipeNumber + JP->EndPipeNumber;
	for (auto sp:spID)
	{
		l = l + L[sp].QCP(L[sp].NN);
		B = B + L[sp].CQP(L[sp].NN);
	}
	for (auto ep : epID)
	{
		m = m - L[ep].QCM(0);
		B = B + L[ep].CQM(0);
	}
	double hp = (l + m) / B;
	for (auto sp : spID)
	{
		L[sp].HP[L[sp].NN] = hp;
		L[sp].QP[L[sp].NN] = L[sp].QCP(L[sp].NN) - L[sp].CQP(L[sp].NN)*L[sp].HP[L[sp].NN];
	}
	for (auto ep : epID)
	{
		L[ep].HP[0] = hp;
		L[ep].QP[0] = L[ep].QCM(0) + L[ep].CQM(0)*L[ep].HP[0];
	}
};

Reservoirs::Reservoirs(int id, int nodeID, double waterLevel) :
	id(id), nodeID(nodeID), waterLevel(waterLevel)
{
};

void Reservoirs::Connect(vector<Nodes>& nodes)
{
	//获取水库所在节点
	Nodes& node = nodes[nodeID];
	// 判断水库类型
	//上游水库
	if (node.n_sp == 0 && node.n_ep == 1)
	{
		type = RE::Upper;
		epID = node.epID[0];
		node.type = BD::UpperReservoir;
	}
		//下游水库
	else if (node.n_ep == 0 && node.n_sp == 1)
	{
		type = RE::Lower;
		spID = node.spID[0];
		node.type = BD::LowerReservoir;
	}
	else
	{
		cout << "水库-管道定义错误！" << endl;
		system("pause");
		exit(0);
	}
	node.typeID = id;
}

void Reservoirs::MOC(vector<Pipes>& pipes)

{
	if (type == RE::Upper)
	{
		Pipes& pipe = pipes[epID];
		/*pipe.HP[0] = waterLevel;
		pipe.QP[0] = pipe.QCM(0) + pipe.CQM(0) * pipe.HP[0];*/
		double a2 =  pipe.m4 * pipe.CQM(0);
		pipe.QP[0] = (sqrt(1 + 4 * a2 * (pipe.QCM(0) + pipe.CQM(0)*waterLevel)) - 1) / 2 / a2;
		pipe.HP[0] = waterLevel - fabs(pipe.QP[0]) * pipe.QP[0] * pipe.m4;
	}
	else if (type == RE::Lower)
	{
		Pipes& pipe = pipes[spID];
		int NN = pipe.NN;
		pipe.HP[NN] = waterLevel;
		pipe.QP[NN] = pipe.QCP(NN) - pipe.CQP(NN) * pipe.HP[NN];
	}
}
 
Turbines::Turbines(int id, string name, int nodeID, double n0, double H0, double Q0, double M0, double D1, double D2,
                   double GD2, double elevation,
                   double convertTime, double guideHeight, bool spine, bool writeSuter) :
	id(id), name(name), nodeID(nodeID), n0(n0), H0(H0), Q0(Q0), M0(M0), D1(D1), D2(D2), GD2(GD2), elevation(elevation),
	convertTime(convertTime),
	guideHeight(guideHeight), spine(spine), writeSuter(writeSuter)
{
	if (n0 > 0)
		isTurMode = true;
	else isTurMode = false;
	//Suter变换
	SuterTransform();
	//保存suter曲线
	if (writeSuter)
		SaveSuter();
	//构建三次样条
	if (spine)
		SpineSuter();
	//读取导叶关闭规律
	ReadGuide();
}

void Turbines::Connect(vector<Nodes>& nodes)
{
	Nodes& node = nodes[nodeID];
	node.type = BD::Turbine;
	node.typeID = id;
	spID = node.spID[0];
	epID = node.epID[0];
};

void Turbines::SuterTransform()
{
	/*
	 *从txt文件中读取特性曲线
	 *对特性曲线进行suter变换
	*/

	//打开特性曲线文件
	string ccName = name + "_cc.txt";
	//读取数据
	mat cc;
	cc.load(ccName);
	//曲线根数
	LCC = cc.n_rows / 3;
	//每根曲线上的点数
	NCC = cc.n_cols - 1;
	//截取出开度序列
	guide = cc.submat(0, 0, LCC - 1, 0);
	//去掉开度列
	cc = cc.cols(1, cc.n_cols - 1);
	//读取单位参数
	mat n = cc.rows(0, LCC - 1);
	mat Q = cc.rows(LCC, LCC * 2 - 1);
	mat M = cc.rows(LCC * 2, LCC * 3 - 1);


	//H0 = 30, Q0 = 0.15, M0 =941;
	double n01 = n0 * D1 / sqrt(H0);
	double Q01 = Q0 / D1 / D1 / sqrt(H0);
	double M01 = M0 / pow(D1, 3) / H0;
	//相对单位参数
	mat nr = n / n01, Qr = Q / Q01, Mr = M / M01;

	////Suter变换
	mat xTemp, WHTemp, WBTemp;
	// x=arctan[(q+k1)/n]
	xTemp = atan((Qr + k1) / nr);
	// x=arctan[(q+k1)/n]+pi (n<=0)
	xTemp.for_each([](auto& val)
	{
		if (val < 0)
			val += 3.1415926;
	});
	// WH=(y+cy)^2/(q^2+a^2+ch)
	WHTemp = 1 / (Qr % Qr + nr % nr + Ch);
	for (int i = 0; i < LCC; i++)
		WHTemp.row(i) = pow(guide(i) + Cy, 2) * WHTemp.row(i);
	// Wb=(m+k2)*WH
	WBTemp = (Mr + k2) % WHTemp;

	//排序
	x.resize(LCC, NCC);
	WH.resize(LCC, NCC);
	WB.resize(LCC, NCC);
	for (int i = 0; i < LCC; i++)
	{
		uvec indices = sort_index(xTemp.row(i));
		for (int j = 0; j < NCC; j++)
		{
			int k = indices(j);
			x(i, j) = xTemp(i, k);
			WH(i, j) = WHTemp(i, k);
			WB(i, j) = WBTemp(i, k);
		}
	}
}

void Turbines::SaveSuter()
{
	//导出suter变换曲线
	FILE* suterFile;
	fopen_s(&suterFile, "sut1erCC.xls", "w");
	fprintf(suterFile, "ANG%t	WH%t	WB\n");
	for (int m2 = 0; m2 < LCC; m2++)
	{
		for (int m1 = 0; m1 < NCC; m1++)
			fprintf(suterFile, "%lf	%lf	%lf\n", x(m2, m1), WH(m2, m1), WB(m2, m1));
	}
	fclose(suterFile);
}

void Turbines::SpineSuter()
{
	//三次样条建模，获取系数h,MH,MB，参考数值分析--颜庆津。
	//初始化变量
	mat hyh = zeros<mat>(LCC, NCC), hyb = zeros<mat>(LCC, NCC);
	mat ar = zeros<mat>(LCC, NCC), ga = zeros<mat>(LCC, NCC), bh = zeros<mat>(LCC, NCC), bb = zeros<mat>(LCC, NCC), hd =
		    zeros<mat>(LCC, NCC);
	h.zeros(LCC, NCC);
	MH.reshape(LCC, NCC);
	MB.reshape(LCC, NCC);
	//A
	for (int j = 1; j < NCC; j++)
	{
		h.col(j) = x.col(j) - x.col(j - 1);
		hyh.col(j) = WH.col(j) - WH.col(j - 1);
		hyb.col(j) = WB.col(j) - WB.col(j - 1);
	}
	for (int j = 1; j < NCC - 1; j++)
	{
		hd.col(j) = h.col(j) + h.col(j + 1);
		ar.col(j) = h.col(j + 1) / hd.col(j);
		ga.col(j) = 1 - ar.col(j);
		bh.col(j) = 6 / hd.col(j) % (hyh.col(j + 1) / h.col(j + 1) - hyh.col(j) / h.col(j));
		bb.col(j) = 6 / hd.col(j) % (hyb.col(j + 1) / h.col(j + 1) - hyb.col(j) / h.col(j));
	}
	for (int i = 0; i < LCC; i++)
	{
		mat A(NCC, NCC, fill::zeros);
		for (int j = 0; j < NCC; j++)
		{
			if (j == 0)
			{
				A(0, 0) = -h(i, 2);
				A(0, 1) = h(i, 1) + h(i, 2);
				A(0, 2) = -h(i, 1);
			}
			else if (j == NCC - 1)
			{
				A(NCC - 1, NCC - 1) = -h(i, NCC - 1);
				A(NCC - 1, NCC - 2) = h(i, NCC - 2) + h(i, NCC - 1);
				A(NCC - 1, NCC - 3) = -h(i, NCC - 1);
			}
			else
			{
				A(j, j - 1) = ga(i, j);
				A(j, j) = 2;
				A(j, j + 1) = ar(i, j);
			}
		}
		MH.row(i) = trans(solve(A, trans(bh.row(i)), solve_opts::refine));
		MB.row(i) = trans(solve(A, trans(bb.row(i)), solve_opts::refine));
	}
}

void Turbines::ReadGuide()
{
	//读取导叶关闭规律文件
	string gName = name + "_guide.txt";
	mat guideData;
	//读取文件
	guideData.load(gName);
	//解析数据
	ServoGuide = guideData.submat(1, 0, guideData(0, 0), 1);
	TimeServo = guideData.submat(guideData(0, 0) + 1, 0, guideData.n_rows - 1, 1);
};

void Turbines::GetGuide(vec& timeSeries)
{
	//按时间序列插值得到相应的导叶开度序列
	vec servo;
	interp1(TimeServo.col(0), TimeServo.col(1), timeSeries, servo, "*linear");
	interp1(ServoGuide.col(0), ServoGuide.col(1), servo, guideSeries, "*linear");
};

double Turbines::GetGuide(double T)
{
	//某时间对应的导叶开度
	double servo=interp1(TimeServo.col(0), TimeServo.col(1), T);
	return interp1(ServoGuide.col(0), ServoGuide.col(1), servo);
};

double Turbines::GetWH(double gui, double xi)
{
	int m;
	double g1, g2;
	double y[2];
	//查找相邻开度线,第一第二最小值对应开度线
	vec delta = abs(guide - gui);
	uvec index = sort_index(delta);
	g1 = guide(index(0)), g2 = guide(index(1));
	//在两根开度线上分别插值
	for (int i = 0; i < 2; i++)
	{
		int k = index(i);
		//样条插值
		if (spine)
		{
			for (m = 1; m < NCC; m++)
			{
				if ((xi >= x(k, m) && xi <= x(k, m - 1)) || (xi <= x(k, m) && xi >= x(k, m - 1)))
				{
					y[i] = MH(k, m - 1) / 6 / h(k, m) * pow(x(k, m) - xi, 3) + MH(k, m) / 6 / h(k, m) *
						pow(xi - x(k, m - 1), 3) + (WH(k, m - 1) / h(k, m) - MH(k, m - 1) / 6 * h(k, m)) * (x(k, m) - xi
						) + (WH(k, m) / h(k, m) - MH(k, m) / 6 * h(k, m)) * (xi - x(k, m - 1));
					break;
				}
			}
			//xi超过区间
			if(m==NCC)
			error("特性曲线插值超出范围");
		}
		else
			//线性插值
			y[i] = interp1(x.row(k), WH.row(k), xi);
	}
	return Linear(gui, g1, g2, y[0], y[1]);
}

double Turbines::GetWB(double gui, double xi)
{
	int m;
	double g1, g2;
	double y[2];
	//查找相邻开度线,第一第二最小值对应开度线
	vec delta = abs(guide - gui);
	uvec index = sort_index(delta);
	g1 = guide(index(0)), g2 = guide(index(1));
	//在两根开度线上分别插值
	for (int i = 0; i < 2; i++)
	{
		int k = index(i);
		//样条插值
		if (spine)
		{
			for (m = 1; m < NCC; m++)
			{
				if ((xi >= x(k, m) && xi <= x(k, m - 1)) || (xi <= x(k, m) && xi >= x(k, m - 1)))
				{
					y[i] = MB(k, m - 1) / 6 / h(k, m) * pow(x(k, m) - xi, 3) + MB(k, m) / 6 / h(k, m)
						* pow(xi - x(k, m - 1), 3) + (WB(k, m - 1) / h(k, m) - MB(k, m - 1) / 6 * h(k, m))
						* (x(k, m) - xi) + (WB(k, m) / h(k, m) - MB(k, m) / 6 * h(k, m)) * (xi - x(k, m - 1));
					break;
				}
			}
			//xi超过区间
			if (m == NCC)
			error("特性曲线插值超出范围");
		}
		else
			y[i] = interp1(x.row(k), WB.row(k), xi);
	}
	return Linear(gui, g1, g2, y[0], y[1]);
}

double Turbines::fQdfQ(double Q, double& fQ)
{
	double q = Q / Q0;
	double xi = atan((q + k1) / 1); //x
	double A0, A1;
	int k;
	for (k = 1; k < NCC; k++)
	{
		if ((xi >= x0[k] && xi <= x0[k - 1]) || (xi <= x0[k] && xi >= x0[k - 1]))
		{
			A0 = WH0[k - 1];
			A1 = (WH0[k] - WH0[k - 1]) / (x0[k] - x0[k - 1]);
			double WHi = A0 + A1 * (xi - x0[k - 1]);
			fQ = H0 * (1 + pow(q, 2)) * WHi / (pow(guide0 + Cy, 2) - WHi * Ch);
			break;
		}
	}
	if (k == NCC)
		error("恒定流计算：水轮机边界溢出！");
	return H0 * (2 * (A0 + A1 * xi) * q / Q0 + A1 / Q0);
};

void Turbines::IniGuideLine()
{
	int m1, m2;
	//初始开度
	guide0 = guideSeries(0);
	//查找相邻两根开度线
	interp1_index(guide, guide0, m1, m2);
	x0.resize(NCC), WH0.resize(NCC);
	//在两根开度线之间插值，得到初始开度线
	for (int i = 0; i < NCC; i++)
	{
		x0(i) = Linear(guide0, guide[m1], guide[m2], x(m1, i), x(m2, i));
		WH0(i) = Linear(guide0, guide[m1], guide[m2], WH(m1, i), WH(m2, i));
	}
}

void Turbines::MOC(vector<Pipes>& pipes, double T)
{
	//用于MOC
	static vector<double> H_t, n_t, Q_t;
	Pipes& sp = pipes[spID];
	Pipes& ep = pipes[epID];
	int & NN = sp.NN;
	static int num;
	double Qi, ni, Hi; //To compare with estimated q,n
	double guidei;
	//初始化
	if (T < 0.000001)
	{
		//分配内存
		H_t.resize(3);
		n_t.resize(3);
		Q_t.resize(3);
		timeSeries.resize(nTPoint);
		nlist.resize(nTPoint);
		Qlist.resize(nTPoint);
		n11.resize(nTPoint);
		Q11.resize(nTPoint);
		DH.resize(nTPoint);
		Pdra.resize(nTPoint);
		Pspi.resize(nTPoint);
		num = 0;
		for (int j = 0; j < n_t.size(); j++)
			n_t[j] = n0, Q_t[j] = ep.Q[0], H_t[j] = sp.H[NN] - ep.H[0];
	}
	//当前导叶开度
	guidei = GetGuide(T);
	//预估转速，流量，水头，等于上一时刻物理量
	n_t[2] = n_t[1];
	Q_t[2] = Q_t[1];
	///////开始迭代
	int iter;
	for (iter = 0; iter < 100; iter++)
	{
		double n_R, Q_R; //相对参数
		n_R = n_t[2] / n0, Q_R = Q_t[2] / Q0, Hi= H_t[2]; //Relative n,Q
		double xi = atan((Q_R + k1 * sqrt(Hi / H0)) / n_R); //Horizontal ordinate
		if (n_R < 0)
			xi += PI;
		double WHi = GetWH(guidei, xi);
		double WBi = GetWB(guidei, xi);
		double H_R = WHi * (n_R * n_R + Q_R * Q_R) / (pow(guidei + Cy, 2) - WHi * Ch);
		double M_R = (WBi / WHi - k2) * H_R;
		Hi = H_R * H0;
	/*	if (guidei)*/
			Qi = (sp.QCP(NN) * ep.CQM(0) + ep.QCM(0) * sp.CQP(NN) - ep.CQM(0) * sp.CQP(NN) * Hi) / (ep.CQM(0) + sp.
				CQP(NN));
		//else
		//	Qi = 0;
		/////////
		if (T > convertTime)
			//甩负荷
			ni = n_t[1] + M0 / GD2 * ep.timeStep * M_R * 0.375;
		else
			//正常运行
			ni = n_t[1];
		if (((fabs(Qi - Q_t[2]) / Q_t[2]) < 0.00001) && (fabs(ni - n_t[2]) < 0.00001))
			break;
		n_t[2] = ni;
		Q_t[2] = Qi;
		H_t[2] = Hi;
	}

	if (iter == 100)
		cout << "水轮机瞬态过程迭代超出范围！" << endl;
	else
	{
	}
	for (int j = 0; j < 2; j++)
	{
		n_t[j] = n_t[j + 1];
		Q_t[j] = Q_t[j + 1];
		H_t[j] = H_t[j + 1];
	}
	sp.QP[NN] = Qi;
	sp.HP[NN] = sp.QCP(NN) / sp.CQP(NN) - 1.0 / sp.CQP(NN) * sp.QP[NN];
	ep.QP[0] = Qi;
	ep.HP[0] = ep.QP[0] / ep.CQM(0) - ep.QCM(0) / ep.CQM(0);
	timeSeries[num] = T;
	guideSeries[num] = guidei;
	nlist[num] = ni;
	Qlist[num] = Qi;
	n11[num] = ni / sqrt(Hi) * D1;
	Q11[num] = Qi / sqrt(Hi) / D1 / D1 * 100; //Scale
	DH[num] = Hi;
	Pspi[num] = sp.HP[NN] - elevation- Qi*Qi/2/g/(PI*D1*guideHeight)/ (PI*D1*guideHeight);
	Pdra[num] = ep.HP[0] - elevation - Qi*Qi / 2 / g / ep.A / ep.A;
	num++;
}

void Turbines::SaveTurbine()
{
	ofstream ft;
	mat data1 = join_rows(timeSeries, guideSeries,nlist,Qlist);
	mat data2 = join_rows( Pspi, Pdra, n11, Q11);
	mat data = join_rows(data1, data2);
	field<string> header(data.n_cols);

	header(0) = "时间";
	header(1) = "开度";
	header(2) = "转速";
	header(3) = "流量";
	header(4) = "蜗壳压力";
	header(5) = "尾水管压力";
	header(6) = "单位转速";
	header(7) = "单位流量";
	saveMat(name+"_机组J" + to_string(id+1) + "非恒定流信息.xls", header,data);
}
