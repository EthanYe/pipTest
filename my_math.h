#pragma once
#include<armadillo>
//#include <algorithm>

using namespace std;
using namespace arma;
double Linear(double, double, double, double, double);

void saveMat(string filename, field<string> const& header, mat const& data);

template <class T>
double interp1(T x, T y, double xi)
{
        vec XI(1), YI(1);
        XI = xi;
        YI = xi;
        interp1(x, y, XI, YI,"*linear");
        return YI(0);
}

void interp1_index(vec x, double xi, int& index1, int& index2);

void error(string message);

template <class T>
void merge(vector<T>& vec1, vector<T>& vec2, vector<T>& vec3)
{
        vec3.insert(vec3.end(), vec1.begin(), vec1.end());
        vec3.insert(vec3.end(), vec2.begin(), vec2.end());
}

template <class T,class T2>
void sort(vector<T>& vec, T2 func)
{
        //��id��С����
        sort(vec.begin(),vec.end(), func);
}

template <class T>
void AutoID(vector<T>& vec)
{
        //��id��С����
        sort(vec, [](const T& p1, const T& p2) {return p1.id< p2.id; });
        //���id������
        for (int i = 0;i < vec.size();i++)
                vec[i].id=i;
}


template <class T, class T2,class T3>
bool find(vector<T> const& vec, T2 member, T3 func, int &index)
{//���ң��Զ��庯������T �� T2
        for (auto it= vec.begin();it!=vec.end();it++)
                if (func(*it, member))
                {
                        index = it-vec.begin();
                        return true;
                }
        return false;
}

template <class T, class T2>
T find(vector<T>& vec,  T2 func,int &index)
{//���ң�����true���������ֵ��p
        //�Զ��庯������T �� T

        T obj = vec[0];
        for (auto it = vec.begin();it != vec.end();it++)
                if (func(it, obj))
                {
                        obj = *it;
                        index = it - vec.begin();
                        return obj;
                }
}

