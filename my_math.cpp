//#include<iostream>
#include"my_math.h"
double Linear(double x, double x1, double x2, double y1, double y2)
{
	double y;
	y = y1 + (x - x1) / (x2 - x1)*(y2 - y1);
	return y;
}

void saveMat(string filename, field<string> const& header,mat const& data)
{
	ofstream pfile;
	pfile.open(filename, ios::out);
	for (int j = 0;j<data.n_cols;j++)
	pfile << header(j) << '\t' ;
	pfile << '\n';
	for (int i = 0;i < data.n_rows;i++)
	{
		for (int j = 0;j < data.n_cols;j++)
			pfile << data(i, j) << '\t';
		pfile << '\n';
	}
}
void interp1_index(vec x,  double xi, int& index1,int& index2)
{
	vec delta = abs(x - xi);
	uvec index = sort_index(delta);
	index1 = index(0);
	index2= index(1);
}

void error(string message)
{
	cout << message << endl;
	system("pause");
	exit(0);
}