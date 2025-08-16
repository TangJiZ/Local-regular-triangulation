#pragma once
#include <unordered_set>

#include "Point_d.h"
#include "Local_Regular.h"

constexpr double epsilon = 1e-15;

/// @brief ���������power distance��point2��point1�ľ��룬��ȥ����point1��Ȩֵ
/// @param point1 �����1
/// @param point2 �����2
/// @return power distance
double compute_weight_distance(Point_d point1, Point_d point2);

/// @brief ��ⷽ�̼�����������
/// @param point1 �����1
/// @param point2 �����2
/// @param point3 �����3
/// @param center ��������
void compute_circle_alg(Point_d point1, Point_d point2, Point_d point3, Point_d& center);

/// @brief ��ʱ��������жϣ������c��p2��p1�ķ���
/// @param p1 
/// @param p2 
/// @param c ��Ե�
/// @return ��ʱ���򷵻�true��˳ʱ���򷵻�false
bool ccw_compute(Point_d p1, Point_d p2, Point_d c);

/// @brief ���������powerƽ����
/// @param p1 
/// @param p2 
/// @param a ƽ����ͨ�÷��̲���A��B��C
/// @param b 
/// @param c 
void weightedBisector(Point_d p1, Point_d p2, double& a, double& b, double& c);

/// @brief ��ֱ�߽������
/// @param A1 ֱ��1ͨ�÷��̲���
/// @param B1 
/// @param C1 
/// @param A2 ֱ��2ͨ�÷��̲���
/// @param B2 
/// @param C2 
/// @param x ����x
/// @param y ����y
void intersection(double A1, double B1, double C1, double A2, double B2, double C2, double& x, double& y);

/// @brief �����Ƿ���ֱ��ͬ���ж�
/// @param point1 
/// @param point2 
/// @param A ֱ��ͨ�÷���
/// @param B 
/// @param C 
/// @return ������ֱ��ͬ��true������ͬ��false
bool same_side_judge(Point_d point1, Point_d point2, double A, double B, double C);

/// @brief p1��c��Ȩֵƽ���ߵķ����������㣬p1c��p2c��Ȩֵƽ���߹�ͬ����c��power diagram���˷���Ϊ��������
/// @param p1 �����
/// @param p2 �ο���
/// @param c ���ĵ�
/// @param x ������������
/// @param y 
void dv_compute(Point_d p1, Point_d p2, Point_d c, double& x, double& y);

/// @brief ����ֱ�߹�ʽ����
/// @param p1 �����
/// @param p2 
/// @param A ��������ֱ��ͨ�÷��̲���
/// @param B 
/// @param C 
void line(Point_d p1, Point_d p2, double& A, double& B, double& C);

/// @brief �жϵ�b�����c��a��˳ʱ�뷽������ʱ�뷽�򣬲������
/// @param a 
/// @param b 
/// @param c ��Ե�
/// @return ˳ʱ��false����ʱ��true
bool cross_product(Point_d a, Point_d b, Point_d c);

/// @brief compare_RT�������ͣ������ȽϽ�������������Ķ��в���
struct CompareResult
{
	bool equal = true;
	std::unordered_set<RTriangle, RTriangleHash> unSet1;
	std::unordered_set<RTriangle, RTriangleHash> unSet2;
};

/// @brief �Ƚ�����Local_Regular���������Ƿ���ͬ
/// @param rt1 ���Ƚ�������1
/// @param rt2 ���Ƚ�������2
/// @return �ȽϽ�������в���
CompareResult compare_RT(Local_Regular* rt1, Local_Regular* rt2);


void compareRTSaveAndPlot(CompareResult& compareRt, std::string fileName, std::string path = "./", bool plot = false);