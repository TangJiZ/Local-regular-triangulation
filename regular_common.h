#pragma once
#include "Point_d.h"

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

/// @brief 
/// @param p1 
/// @param p2 
/// @param c 
/// @param x 
/// @param y 
void dv_compute(Point_d p1, Point_d p2, Point_d c, double& x, double& y);