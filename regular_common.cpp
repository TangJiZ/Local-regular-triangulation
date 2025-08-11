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
	c = ((p2.x * p2.x + p2.y * p2.y) - (p1.x * p1.x + p1.y * p1.y) - (p1.w - p2.w)) / 2;
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
	if (D1 * D2 >= 0)
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
