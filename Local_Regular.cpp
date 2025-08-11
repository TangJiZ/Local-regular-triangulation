#include <iostream>

#include "Local_Regular.h"
#include "regular_common.h"

void Local_Regular::getMinMax_ptr()
{
	xmin = ymin = std::numeric_limits<double>::max();
	xmax = ymax = -std::numeric_limits<double>::max();
	wmax = -std::numeric_limits<double>::max();

	for (size_t i = 0; i < n; ++i) {
		xmin = std::min(xmin, pts[i].x);
		xmax = std::max(xmax, pts[i].x);

		ymin = std::min(ymin, pts[i].y);
		ymax = std::max(ymax, pts[i].y);

		wmax = std::max(wmax, pts[i].w);
	}
}

void Local_Regular::partition()
{
	double x = xmax - xmin;
	double y = ymax - ymin;

	xnum = int(sqrt(n * x / y)) + 1;
	ynum = int(sqrt(n * y / x)) + 1;

	xcell = x / xnum;
	ycell = y / ynum;

	x = x + xcell * 0.01;
	y = y + ycell * 0.01;

	xcell = x / xnum;
	ycell = y / ynum;

	bucket = new int[xnum * ynum];
	next = new int[n];

	for (int i = 0; i < n; i++)
	{
		next[i] = -1;
	}

	for (int i = 0; i < xnum * ynum; i++)
	{
		bucket[i] = -1;
	}

	for (int i = 0; i < n; i++)
	{
		int ix = int((pts[i].x - xmin) / xcell);
		int iy = int((pts[i].y - ymin) / ycell);
		int b = ix * ynum + iy;
		next[i] = bucket[b];
		bucket[b] = i;
	}
}

void Local_Regular::near_find(Point_d* p, std::vector<int>& vp)
{
	int x = int((p->x - xmin) / xcell);
	int y = int((p->y - ymin) / ycell);

	// ��ѯ�ڵ�����cell
	int now = bucket[x * ynum + y];
	while (now != -1)
	{
		Point_d* pn = pts + now;
		if (pts + now != p)
		{
			vp.push_back(now);
		}
		now = next[now];
	}

	int step = 1;
	do
	{
		int fy_min, fy_max, fx_min, fx_max;
		bool hy_min, hy_max, hx_min, hx_max;
		if (x + step < xnum)
		{
			fx_max = x + step;
			hx_max = true;
		}
		else
		{
			fx_max = xnum - 1;
			hx_max = false;
		}

		if (x - step < 0)
		{
			fx_min = 0;
			hx_min = false;
		}
		else
		{
			fx_min = x - step;
			hx_min = true;
		}

		if (y + step < ynum)
		{
			fy_max = y + step - 1;
			hy_max = true;
		}
		else
		{
			fy_max = ynum - 1;
			hy_max = false;
		}

		if (y - step < 0)
		{
			fy_min = 0;
			hy_min = false;
		}
		else
		{
			fy_min = y - step + 1;
			hy_min = true;
		}

		if (hx_min)
		{
			for (int i = fy_min; i <= fy_max; i++)
			{
				now = bucket[fx_min * ynum + i];
				while (now != -1)
				{
					Point_d* pn = pts + now;
					vp.push_back(now);
					now = next[now];
				}
			}
		}

		if (hx_max)
		{
			for (int i = fy_min; i <= fy_max; i++)
			{
				now = bucket[fx_max * ynum + i];
				while (now != -1)
				{
					Point_d* pn = pts + now;
					vp.push_back(now);
					now = next[now];
				}
			}
		}

		if (hy_min)
		{
			for (int i = fx_min; i <= fx_max; i++)
			{
				now = bucket[i * ynum + fy_min - 1];
				while (now != -1)
				{
					Point_d* pn = pts + now;
					vp.push_back(now);
					now = next[now];
				}
			}
		}

		if (hy_max)
		{
			for (int i = fx_min; i <= fx_max; i++)
			{
				now = bucket[i * ynum + fy_max + 1];
				while (now != -1)
				{
					Point_d* pn = pts + now;
					vp.push_back(now);
					now = next[now];
				}
			}
		}
		
		step++;
	} while (vp.size() < 2);
}

void Local_Regular::Power_diagram_build_initial(std::vector<int>& nf, PCell*& first, PCell*& last)
{
	first = new PCell;
	first->p1 = pts;
	first->p2 = infinate;

}

Local_Regular::Local_Regular(Point_d* p, int num)
{
	pts = p;
	n = num;

	// ��ȡ�㼯�����Сֵ
	getMinMax_ptr();

	// ���������Сֵ�γɵľ�������Խڵ���з���
	partition();

	// �Ե㼯��ÿ���ڵ���оֲ����㣬�õ��ڵ��ȫ���ھ�����
	for (int i = 0; i < n; i++)
	{
		// ����չ�����ѯ�ڽ��ڵ㣬��֤�ڵ���������Ϊ2
		std::vector<int> nf;
		near_find(pts + i, nf);

		//std::cout << "point: " << i << std::endl;
		//for (int i = 0; i < nf.size(); i++)
		//{
		//	std::cout << "near " << i << ": " << nf[i] << std::endl;
		//}
		//std::cout << std::endl;

		// power diagram������������
		PCell* first, last;
	}
}
