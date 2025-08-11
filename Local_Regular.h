#pragma once
#include <vector>

#include "Point_d.h"

/// @brief power图中边界点存储结构，链表单元结构
struct PCell
{
	Point_d q;				// 边界交点，power外接圆圆心；p1,p2中存在infinate点时，用做方向向量
	double r;				// 半径
	Point_d* p1, * p2;		// 形成此交点的两邻居节点
	PCell* next;				// power单元边界逆时针序的下一交点
};

/// @brief 局部方法计算regular三角剖分类
class Local_Regular
{
/// @brief value
private:
	/// @brief 无穷点
	Point_d* infinate;

	/// @brief 待计算点集和点集中节点数量
	Point_d* pts;	// 待计算点集
	int n;			// 节点数量

	/// @brief 点集中xy最大最小值，w的最大值
	double xmin, xmax, ymin, ymax, wmax;

	/// @brief 点集分区相关参数
	int xnum, ynum;				// x，y方向分区数量
	double xcell, ycell;		// 分区单元大小
	int* bucket;				// 存储每个分区中第一个节点
	int* next;					// 存储每个节点在分区中的下一个节点，多节点同分区时发挥作用

public:


/// @brief function
private:
	/// @brief 获取待计算点集xy最大最小值，w最大值
	void getMinMax_ptr();

	/// @brief 待计算点集区域划分
	void partition();

	/// @brief 外扩展计算点p的邻近节点，保证节点数量大于2，逆螺旋至少为周围8格
	/// @param p 待查找节点
	/// @param vp 保存结果，p的邻近节点
	void near_find(Point_d* p, std::vector<int>& vp);

	/// @brief power diagram初始构造，要求交点逆时针序排列，只加入前两个节点
	/// @param nf 邻近节点
	/// @param first power diagram交点链表头
	/// @param last power diagram交点链表尾
	void Power_diagram_build_initial(std::vector<int>& nf, PCell*& first, PCell*& last);

public:
	/// @brief 构造函数同时内部执行剖分
	/// @param p 待计算点集
	/// @param num 点集中节点数量
	Local_Regular(Point_d* p, int num);
};

