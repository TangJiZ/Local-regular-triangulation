#pragma once
#include <vector>

#include "Point_d.h"

/// @brief 在PCell扩展权值圆内的节点，链表存储
struct inPCell
{
	Point_d* p;
	inPCell* next;
};

/// @brief power图中边界点存储结构，链表单元结构
struct PCell
{
	Point_d q;				// 边界交点，power外接圆圆心；p1,p2中存在infinate点时，用做方向向量
	Point_d* p1, * p2;		// 形成此交点的两邻居节点
	PCell* next;				// power单元边界逆时针序的下一交点
	inPCell* ipc = nullptr;
	double r;				// 外扩展最大权值后的半径
};

/// @brief 圆外网格包围盒
struct Box
{
	int fy_min, fy_max, fx_min, fx_max;
};

/// @brief 节点在regular triangulaton中的相邻连接节点链表
struct Nei
{
	Point_d* p;
	Nei* next;
};

/// @brief regular triangle存储
struct RTriangle
{
	Point_d a, b, c;

	/// @brief 重载==运算符
	/// @param other 待比较对象
	/// @return 是否相同
	bool operator==(const RTriangle& other) const;
};

/// @brief RTriangleHash结构体，用于包含RTriangle的unorder_set等集合类型创建时的hash模板输入
struct RTriangleHash
{
	size_t operator()(const RTriangle& rt) const
	{
		double all = 0;
		all += rt.a.x + rt.a.y;
		all += rt.b.x + rt.b.y;
		all += rt.c.x + rt.c.y;
		/// 截断double到小数点后15位，double的误差会导致hash的计算不同，从而使得比较时产生错误
		all = floor(all * 1e15) / 1e15;
		return std::hash<double>{}(all);
	}
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

	std::vector<Nei*> nei;		// 存储每个节点计算完成后的相邻邻居节点

public:

	std::vector<RTriangle> rtri;		// 存储最终三角形集合
	int valid_num = 0;					// 有效节点数量
	int invalid_num = 0;				// 无效节点数量


/// @brief function
private:
	/// @brief 获取待计算点集xy最大最小值，w最大值
	void getMinMax_ptr();

	/// @brief 待计算点集区域划分
	void partition();

	/// @brief 外扩展计算点p的邻近节点，保证节点数量大于2，逆螺旋至少为周围8格，同时获取搜索的矩形网格范围
	/// @param p 待查找节点
	/// @param vp 保存结果，p的邻近节点
	/// @param box 最终搜索的矩形网格包围盒 
	void near_find(Point_d* p, std::vector<int>& vp, Box& box);

	/// @brief power diagram初始构造，要求交点逆时针序排列，只加入前两个节点
	/// @param nf 邻近节点
	/// @param c 中心节点
	/// @param first power diagram交点链表头
	/// @param last power diagram交点链表尾
	void power_diagram_build_initial(std::vector<int>& nf, Point_d c, PCell*& first, PCell*& last);

	/// @brief 显示power diagram链
	/// @param iter power diagram链开始指针
	void show_pcell(PCell* iter);

	/// @brief 向power diagram图中插入新节点，并更新图结构，判断中心点是否为多余，可继续计算
	/// @param p 插入节点
	/// @param c 中心节点
	/// @param first power diagram交点链表头，插入完成后更新
	/// @param last power diagram交点链表尾，插入完成后更新
	/// @return 中心点为多余节点无后续计算false，中心点非多余继续计算true 
	bool power_diagram_insert(Point_d* p, Point_d c, PCell* &first, PCell* &last);

	/// @brief 包含无穷点的PCell是否在直线判断点j侧的判断
	/// @param cell 包含无穷点的PCell
	/// @param A 直线统一方程参数
	/// @param B 
	/// @param C 
	/// @param j 判断点j，判断cell是否和j在直线的同一侧
	/// @param c 中心点c
	/// @return 在同一侧true，不在同一侧false
	bool infinate_same_side_judge(PCell* cell, double A, double B, double C, Point_d j, Point_d c);

	/// @brief 显示包含在cell扩展外接圆并在夹角范围中的节点
	/// @param iter power diagram交点
	void show_inPCell(PCell* iter);

/// @brief 计算power diagram链中交点扩展权值圆的包含范围，分为交点形成包含无穷点和不包含无穷点两类
/// ------不包含无穷点------
	
	/// @brief 计算pcell的包围盒
	/// @param pcell 待计算power diagram交点
	/// @param box 矩形包围盒
	void pcell_box(PCell* pcell, Box& box);

	/// @brief 判断点p是否在cell扩展外接圆内，同时在中心点c和其两构成点形成的夹角中
	/// @param p 判断点
	/// @param cell power diagram交点
	/// @param c 中心点
	/// @return 在外扩展圆和夹角中true，其余false
	bool in_circle_and_angle(Point_d p, PCell* cell, Point_d c);

	/// @brief 获取在同时cell扩展外接圆和夹角中的节点，存储在cell中的ipc中，链表
	/// @param cell power diagram交点
	/// @param exbox 中心点初始构造时已查询的包围盒范围
	/// @param c 中心点
	void findPoints_in_PCell(PCell* cell, Box exbox, Point_d c);

/// ------包含无穷点------
	
	/// @brief x方向上，直线和长方形相交部分，每个xid对应的yid范围
	struct Bound
	{
		int x_min, x_max;	// 通常只使用x_min，当直线垂直于x轴时，单个Bound边可表示所有范围
		int y_min, y_max;
	};

	/// @brief 计算外扩展直线和节点外长方形框所围成的网格范围，方向为无穷点对应的方向向量
	/// @param cell 待计算power diagram交点，包含infinate点
	/// @param c 中心点，直线为中心点和非infinate点形成
	/// @param vb 围成的网格范围，x方向逐网格存储
	void pcell_bound_infinate(PCell* cell, Point_d c, std::vector<Bound>& vb);

	/// @brief 获取在外扩展线内的节点，存储在cell中的ipc中，链表
	/// @param cell power diagram交点，包含infinate节点
	/// @param exbox 中心点初始构造时已查询的包围盒范围
	/// @param c 中心点
	void findPoints_in_PCell_infinate(PCell* cell, Box exbox, Point_d c);

/// ------******------

	/// @brief 插入节点后，更新需改动交点对应的inPCell（圆内节点），被替代的交点圆内节点要么在新交点圆外，要么在新交点圆内
	/// @param start 被替代链的链头
	/// @param end 被替代链的链尾
	/// @param c 中心节点
	/// @param new1 新产生交点1
	/// @param new2 新产生交点2
	void update_inPCell(PCell* start, PCell* end, Point_d c, PCell* new1, PCell* new2);

	/// @brief 计算包含无穷点的cell外扩直线，并生成一个在外扩直线内部方向的节点
	/// @param cell	power diagram交点，包含infinate
	/// @param c 中心节点
	/// @param A 外扩直线通用参数
	/// @param B 
	/// @param C 
	/// @param p 在外扩直线内部方向一点
	void compute_ext_line_and_point(PCell* cell, Point_d c, double& A, double& B, double& C, Point_d& p);

	/// @brief 获取最终三角形集合，每个节点只和编号比自身大的节点形成三角形
	void get_triangle();

public:
	/// @brief 构造函数同时内部执行剖分
	/// @param p 待计算点集
	/// @param num 点集中节点数量
	Local_Regular(Point_d* p, int num);

	/// @brief 默认构造函数
	Local_Regular() {}

	/// @brief 读取三角形
	/// @param path 读取文件路径
	void read_triangle(std::string path);

	/// @brief 保存，绘制三角网
	/// @param name 文件名，不含后缀；数据文件自动为文件名+txt，图像文件自动为文件名+jpg
	/// @param path 文件保存路径
	/// @param plot 选择是否绘制图像
	void saveAndPlot(std::string name, std::string path, bool plot);
};

