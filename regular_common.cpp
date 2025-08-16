#include <string>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

#include "regular_common.h"

double compute_weight_distance(Point_d point1, Point_d point2)
{
	return (point1.x - point2.x) * (point1.x - point2.x) + (point1.y - point2.y) * (point1.y - point2.y) - point1.w;
}

void compute_circle_alg(Point_d point1, Point_d point2, Point_d point3, Point_d& center)
{
	double t32 = point3.x * point3.x - point2.x * point2.x + point3.y * point3.y - point2.y * point2.y + point2.w - point3.w;
	double t21 = point2.x * point2.x - point1.x * point1.x + point2.y * point2.y - point1.y * point1.y + point1.w - point2.w;
	double D = (point3.x - point2.x) * (point2.y - point1.y) - (point2.x - point1.x) * (point3.y - point2.y);
	center.x = (t32 * (point2.y - point1.y) - t21 * (point3.y - point2.y)) / (2 * D);
	center.y = -(t32 * (point2.x - point1.x) - t21 * (point3.x - point2.x)) / (2 * D);
	center.w = (center.x - point1.x) * (center.x - point1.x) + (center.y - point1.y) * (center.y - point1.y) - point1.w;
}

bool ccw_compute(Point_d p1, Point_d p2, Point_d c)
{
	double xy1[2]{ p1.x - c.x, p1.y - c.y };
	double xy2[2]{ p2.x - c.x, p2.y - c.y };
	double cross = xy1[0] * xy2[1] - xy2[0] * xy1[1];
	if (cross > 0)
		return true;
	else
		return false;
}

void weightedBisector(Point_d p1, Point_d p2, double& a, double& b, double& c)
{
	a = p2.x - p1.x;
	b = p2.y - p1.y;
	c = ((p1.x * p1.x + p1.y * p1.y) - (p2.x * p2.x + p2.y * p2.y) + (p2.w - p1.w)) / 2;
}

void intersection(double A1, double B1, double C1, double A2, double B2, double C2, double& x, double& y)
{
	x = (B1 * C2 - B2 * C1) / (A1 * B2 - A2 * B1);
	y = (A2 * C1 - A1 * C2) / (A1 * B2 - A2 * B1);
}

bool same_side_judge(Point_d point1, Point_d point2, double A, double B, double C)
{
	double D1 = A * point1.x + B * point1.y + C;
	double D2 = A * point2.x + B * point2.y + C;
	double D = D1 * D2;
	if (abs(D) < epsilon || D >= 0)
		return true;
	else
		return false;
}

void dv_compute(Point_d p1, Point_d p2, Point_d c, double& x, double& y)
{
	double A1, B1, C1;
	weightedBisector(p1, c, A1, B1, C1);
	double A2, B2, C2;
	weightedBisector(p2, c, A2, B2, C2);
	Point_d q;
	intersection(A1, B1, C1, A2, B2, C2, q.x, q.y);
	Point_d m;
	if (B1 == 0)
	{
		m.x = q.x;
		m.y = q.y + 1;
	}
	else if (A1 == 0)
	{
		m.x = q.x + 1;
		m.y = q.y;
	}
	else
	{
		m.x = q.x + 1;
		m.y = -(A1 * m.x + C1) / B1;
	}

	bool side = same_side_judge(m, c, A2, B2, C2);
	// 逆序的下一节点p2，如果其权值超出其和c的欧氏距离的平方
	// 则其power平分线会超过c，使得整体cell脱离c，此时方向向量的判断条件和cell未脱离c时相反
	if ((p2.w - c.w) > (p2.x - c.x) * (p2.x - c.x) + (p2.y - c.y) * (p2.y - c.y))
	{
		if (side)
		{
			x = q.x - m.x;
			y = q.y - m.y;
		}
		else
		{
			x = m.x - q.x;
			y = m.y - q.y;
		}
	}
	else
	{
		if (side)
		{
			x = m.x - q.x;
			y = m.y - q.y;
		}
		else
		{
			x = q.x - m.x;
			y = q.y - m.y;
		}
	}	
}

void line(Point_d p1, Point_d p2, double& A, double& B, double& C)
{
	A = p2.y - p1.y;
	B = -(p2.x - p1.x);
	C = p2.x * p1.y - p1.x * p2.y;
}

bool cross_product(Point_d a, Point_d b, Point_d c)
{
	Point_d ca = { a.x - c.x, a.y - c.y };
	Point_d cb = { b.x - c.x, b.y - c.y };
	double cp = ca.x * cb.y - cb.x * ca.y;
	if (cp > 0)
		return true;
	else
		return false;
}

CompareResult compare_RT(Local_Regular* rt1, Local_Regular* rt2)
{
	std::unordered_set<RTriangle, RTriangleHash> unSet1;
	for (int i = 0; i < rt1->rtri.size(); i++)
	{
		unSet1.insert(rt1->rtri[i]);
	}

	std::unordered_set<RTriangle, RTriangleHash> unSet2;
	for (int i = 0; i < rt2->rtri.size(); i++)
	{
		unSet2.insert(rt2->rtri[i]);
	}

	CompareResult compareResult;
	if (unSet1 == unSet2)
	{
		compareResult.equal = true;
		return compareResult;
	}
	else
	{
		/// 计算unSet1和unSet2的差集
		std::unordered_set<RTriangle, RTriangleHash> differenceSet1;
		for (RTriangle rt : unSet1)
		{
			if (unSet2.find(rt) == unSet2.end())
			{
				differenceSet1.insert(rt);
			}
		}
		std::unordered_set<RTriangle, RTriangleHash> differenceSet2;
		for (RTriangle rt : unSet2)
		{
			if (unSet1.find(rt) == unSet1.end())
			{
				differenceSet2.insert(rt);
			}
		}
		compareResult.equal = false;
		compareResult.unSet1 = differenceSet1;
		compareResult.unSet2 = differenceSet2;
		return compareResult;
	}
}

void compareRTSaveAndPlot(CompareResult& compareRt, std::string fileName, std::string path, bool plot)
{
	using std::string;
	using std::ofstream;
	using std::endl;
	using std::setprecision;

	/// 创建path下的name文件夹
	string filePath = path + '/' + fileName;
	if (!std::filesystem::create_directory(filePath))
	{
		//cout << "文件夹创建失败" << endl;
	}

	/// 保存比较结果总结
	string summaryPath = filePath + "/summary.txt";
	ofstream file(summaryPath);
	file << std::fixed << setprecision(15);
	file << "结果比较是否相同：" << compareRt.equal << endl;
	file << "比较对象1独有三角形个数：" << compareRt.unSet1.size() << endl;
	int i = 0;
	file << "比较对象2独有三角形个数：" << compareRt.unSet2.size() << endl;
	file.close();

	/// 单独保存各对象的独有三角形文件
	string set1Path = filePath + "/set1.txt";
	ofstream set1file(set1Path);
	set1file << std::fixed << setprecision(12);
	i = 0;
	for (const RTriangle& rt : compareRt.unSet1)
	{
		set1file << i << ' ' << rt.a.x << ' ' << rt.a.y << ' ' << rt.a.w << endl;
		set1file << i << ' ' << rt.b.x << ' ' << rt.b.y << ' ' << rt.b.w << endl;
		set1file << i << ' ' << rt.c.x << ' ' << rt.c.y << ' ' << rt.c.w << endl;
		i++;
	}
	set1file.close();

	string set2Path = filePath + "/set2.txt";
	ofstream set2file(set2Path);
	set2file << std::fixed << setprecision(12);
	i = 0;
	for (const RTriangle& rt : compareRt.unSet2)
	{
		set2file << i << ' ' << rt.a.x << ' ' << rt.a.y << ' ' << rt.a.w << endl;
		set2file << i << ' ' << rt.b.x << ' ' << rt.b.y << ' ' << rt.b.w << endl;
		set2file << i << ' ' << rt.c.x << ' ' << rt.c.y << ' ' << rt.c.w << endl;
		i++;
	}
	set2file.close();

	/// 绘图
	if (plot)
	{
		/// set1绘图
		string scriptPath = "D:/Microsoft Visual Studio/code/Regular_Triangulation/regular_2d/code/plotScript.py";
		string name1 = "set1";
		// 构造完整的命令字符串
		std::string command1 = "python \"" + scriptPath + "\" \"" + set1Path + "\" \"" + name1 + "\"";
		// 使用system()函数调用Python脚本
		system(command1.c_str());

		/// set2绘图
		string name2 = "set2";
		// 构造完整的命令字符串
		std::string command2 = "python \"" + scriptPath + "\" \"" + set2Path + "\" \"" + name2 + "\"";
		// 使用system()函数调用Python脚本
		system(command2.c_str());
	}
}
