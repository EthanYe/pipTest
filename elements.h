#ifndef COMPONENTS_H
#define	COMPONENTS_H

#include<vector>
#include<armadillo>
#include<QAction>
#include<QMetaEnum>
using namespace std;
using namespace arma;

const double g = 9.81; // 重力加速度
const double PI = 3.1415926; //圆周率



enum class BD
{
    Series = 0,
    UpperReservoir = 1,
    LowerReservoir = 2,
    Branch = 3,
    Turbine = 4,
    BallValve = 5,
    EndBallValve = 6,
    Pump = 7,
    EndValve = 8,
    number = 9
};

enum class RE
{
    Upper = 0,
    Lower=1
};

enum class POS
{
    Upper = 0,
    Internal=1,
    Lower = 2
};

class Elements
{
    static int counter;
    int id;
    string text;
};

class Pipes
{
public:
    static const int MAX_PIPE = 100;
    //需要输入的管道参数
    int id; //编号
    double length; //长度
    double D; //直径
    int js; //头节点
    int je; //尾节点
    double timeStep; //时间步长
    double a; //波速
    double roughness; //糙率
    double zeta[2]; //局部损失系数
    //计算得到的参数
    int NN; //分段数
    double A; //面积
    double dx; //空间步长
    double f; //沿程阻力系数
    //状态变量
    vec H0; //初始水头
    vec Q0; //初始流量
    vec H; //上一时刻水头
    vec Q; //上一时刻流量
    vec HP; //当前时刻水头
    vec QP; //当前时刻流量

public:
    Pipes();
    Pipes(int id, int js=0, int je=0, double timeStep=0.01, double a=1000, double D=0.1, double len=1, double roughness=0, double zeta1=0,
          double zeta2=0);
    /*
    QP= QCP - CQP * HP
    QP= QCM + CQM * HP
    */
    double QCP(int J);
    double CQP(int J);
    double QCM(int J);
    double CQM(int J);
    void MOC();
    //1 / (2gA ^ 2)
    double m4;
private:
    double m1;
    double m2;
    double m3;

};

class Nodes
{
public:
    int id;
    BD type;
    int typeID;
    int n_sp=0; //起始（流入）管道数
    int n_ep=0; //终止（流出管道数）
    POS position; //上游节点，中间节点，下游节点
    vector<int> spID;
    vector<int> epID;

    int np = 0; //总管道数
    vector<int> pID;//全部管道ID
    Nodes(){};
    Nodes(int id, int typeID = 0, BD type = BD::Series);

    void AddPipe(Pipes& pipe, bool start = true);
    void SeriesMOC(vector<Pipes>& pipes);
    void BranchMOC(vector<Pipes>& pipes);
};

class Reservoirs
{
public:
    int id;
    int nodeID;
    double waterLevel;
    RE type;
    int epID;
    int spID;
public:
    Reservoirs(int id, int nodeID=0, double waterLevel=0);
    void Connect(vector<Nodes>& nodes);
    void MOC(vector<Pipes>& pipes);
};

class Turbines
{
public:
    int id;
    int nodeID;
    string name; //名称
    int nTPoint; //工况轨迹点数
    int NCC; //每根特性曲线点数
    int LCC; //特性曲线条数
    double guide0; //初始导叶开度
    vec x0, WH0;
    bool isTurMode;
    double k1 = 1.2, k2 = 4.2, Ch = 0.4, Cy = 0.2;
    double n0, Q0, M0, H0; //额定参数
    double printDensity; //打印密度
private:
    int spID; //前管道
    int epID; //后管道
    double GD2; //GD2
    double D1; //水轮机直径
    double D2; //尾水管进口直径
    double elevation; //安装高程
    double convertTime; //工况转换时间
    double guideHeight; //导叶高度
    mat h, MH, MB, WH, WB, x; //x[i]-x[i-1],len(NCC+1)
    vec guide; //特性点应的导叶开度
    mat ServoGuide; //接力器行程和导叶相对开度关系
    mat TimeServo; //时间和接力器行程关
    vec guideSeries; //对应时间序列下的导叶开度序列
    vec timeSeries; //对应时间序列下的导叶开度序列
    vec nlist,Qlist,n11, Q11, DH, Pspi, Pdra;
    bool spine; //特性曲线是否样条插值,默认否
    bool writeSuter; //是否导出suter曲线，默认否


public:
    Turbines(int id, string name, int node, double n0,double H0, double Q0,double M0, double D1, double D2, double GD2, double elevation,
             double convertTime,
             double guideHeight, bool spine = false, bool writerSuter = false);
    void Connect(vector<Nodes>& nodes);
    double GetWH(double gui, double xi);
    double GetWB(double gui, double xi);
    void GetGuide(vec&);
    double GetGuide(double T);
    double fQdfQ(double Q, double &fQ);
    void MOC(vector<Pipes>& pipes,double T);
    void IniGuideLine();//初始导叶开度线
    void IniMOC();//初始导叶开度线
    void SaveTurbine();

private:
    void SuterTransform();
    void ReadGuide();
    void SpineSuter();
    void SaveSuter();
};
#endif // !COMPONENTS
