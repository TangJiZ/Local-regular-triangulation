#pragma once
#include <unordered_set>

#include "Point_d.h"
#include "Local_Regular.h"

constexpr double epsilon = 1e-15;

/// @brief 计算两点的power distance，point2到point1的距离，减去的是point1的权值
/// @param point1 计算点1
/// @param point2 计算点2
/// @return power distance
double compute_weight_distance(Point_d point1, Point_d point2);

/// @brief 求解方程计算正交中心
/// @param point1 计算点1
/// @param point2 计算点2
/// @param point3 计算点3
/// @param center 正交中心
void compute_circle_alg(Point_d point1, Point_d point2, Point_d point3, Point_d& center);

/// @brief 逆时针序计算判断，相对于c，p2在p1的方向
/// @param p1 
/// @param p2 
/// @param c 相对点
/// @return 逆时针序返回true，顺时针序返回false
bool ccw_compute(Point_d p1, Point_d p2, Point_d c);

/// @brief 计算两点的power平分线
/// @param p1 
/// @param p2 
/// @param a 平分线通用方程参数A，B，C
/// @param b 
/// @param c 
void weightedBisector(Point_d p1, Point_d p2, double& a, double& b, double& c);

/// @brief 两直线交点计算
/// @param A1 直线1通用方程参数
/// @param B1 
/// @param C1 
/// @param A2 直线2通用方程参数
/// @param B2 
/// @param C2 
/// @param x 交点x
/// @param y 交点y
void intersection(double A1, double B1, double C1, double A2, double B2, double C2, double& x, double& y);

/// @brief 两点是否在直线同侧判断
/// @param point1 
/// @param point2 
/// @param A 直线通用方程
/// @param B 
/// @param C 
/// @return 两点在直线同侧true，不在同侧false
bool same_side_judge(Point_d point1, Point_d point2, double A, double B, double C);

/// @brief p1和c的权值平分线的方向向量计算，p1c和p2c的权值平分线共同构成c的power diagram，此方向为方向向量
/// @param p1 计算点
/// @param p2 参考点
/// @param c 中心点
/// @param x 方向向量坐标
/// @param y 
void dv_compute(Point_d p1, Point_d p2, Point_d c, double& x, double& y);

/// @brief 两点直线公式计算
/// @param p1 计算点
/// @param p2 
/// @param A 两点所在直线通用方程参数
/// @param B 
/// @param C 
void line(Point_d p1, Point_d p2, double& A, double& B, double& C);

/// @brief 判断点b相对于c在a的顺时针方向还是逆时针方向，叉积计算
/// @param a 
/// @param b 
/// @param c 相对点
/// @return 顺时针false，逆时针true
bool cross_product(Point_d a, Point_d b, Point_d c);

/// @brief compare_RT返回类型，包含比较结果，各三角网的独有部分
struct CompareResult
{
	bool equal = true;
	std::unordered_set<RTriangle, RTriangleHash> unSet1;
	std::unordered_set<RTriangle, RTriangleHash> unSet2;
};

/// @brief 比较两个Local_Regular的三角网是否相同
/// @param rt1 待比较三角网1
/// @param rt2 待比较三角网2
/// @return 比较结果及独有部分
CompareResult compare_RT(Local_Regular* rt1, Local_Regular* rt2);


void compareRTSaveAndPlot(CompareResult& compareRt, std::string fileName, std::string path = "./", bool plot = false);