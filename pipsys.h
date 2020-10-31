#ifndef PIPSYS_H
#define PIPSYS_H
#include<QObject>
#include<string>
#include<QString>
#include<QStringList>
#include"elements.h"
using namespace std;


class Pipsys :public QObject {
    Q_OBJECT
public:
	const double g = 9.81; // 重力加速度
	const double PI = 3.1415926; //圆周率
	double timeStep; //时间步长
	double totalTime; //计算总时间
	double steps; //计算步数
	int printDensity; //打印密度
	string name; //系统名称
	vec timeSeries; //时间序列
	double T=0; //当前时间

    QStringList elementsList;
	//定义对象
	int nPipes; // 管道数
	vector<Pipes> pipes;   //管道
	int nNodes; //节点数
	vector<Nodes> nodes; //节点
	//水库
    int nReservoirs=0; // 水库数
	int nUpperReservoirs=0;// 上游水库数
	int nLowerReservoirs=0;// 下游水库数 
	vector<Reservoirs> reservoirs; //水库
	//水轮机
	int nTurbines; //水轮机数
	vector<Turbines> turbines; //水轮机
	bool isSteady = false;
	double iniQ = 0.1;
    QWidget *parrent;
	
public:
    Pipsys(QWidget *);
    void Run();
    void addElement(QString);
void addReservoir();
void addPipe();
void addNode();
private:
	void Steady();
	void MOC();
	void SteadyIni(mat& SY,mat& SJ ,mat& S,mat& QS, vec& E);
	void Q2E(Nodes& node, mat &E, mat &S, mat &QS, vector<bool>& visit);
    //ui

};


#endif // !PIPSYS_H

