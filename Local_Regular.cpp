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

void Local_Regular::near_find(Point_d* p, std::vector<int>& vp, Box& box)
{
	int x = int((p->x - xmin) / xcell);
	int y = int((p->y - ymin) / ycell);

	// 查询节点所在cell
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
	bool hy_min, hy_max, hx_min, hx_max;
	do
	{
		if (x + step < xnum)
		{
			box.fx_max = x + step;
			hx_max = true;
		}
		else
		{
			box.fx_max = xnum - 1;
			hx_max = false;
		}

		if (x - step < 0)
		{
			box.fx_min = 0;
			hx_min = false;
		}
		else
		{
			box.fx_min = x - step;
			hx_min = true;
		}

		if (y + step < ynum)
		{
			box.fy_max = y + step - 1;
			hy_max = true;
		}
		else
		{
			box.fy_max = ynum - 1;
			hy_max = false;
		}

		if (y - step < 0)
		{
			box.fy_min = 0;
			hy_min = false;
		}
		else
		{
			box.fy_min = y - step + 1;
			hy_min = true;
		}

		if (hx_min)
		{
			for (int i = box.fy_min; i <= box.fy_max; i++)
			{
				now = bucket[box.fx_min * ynum + i];
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
			for (int i = box.fy_min; i <= box.fy_max; i++)
			{
				now = bucket[box.fx_max * ynum + i];
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
			for (int i = box.fx_min; i <= box.fx_max; i++)
			{
				now = bucket[i * ynum + box.fy_min - 1];
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
			for (int i = box.fx_min; i <= box.fx_max; i++)
			{
				now = bucket[i * ynum + box.fy_max + 1];
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

	// 恢复fy_max，fy_min，用于后续的节点查找
	if (hy_max)
		box.fy_max++;
	if (hy_min)
		box.fy_min--;
}

void Local_Regular::power_diagram_build_initial(std::vector<int>& nf, Point_d c, PCell*& first, PCell*& last)
{
	first = new PCell;
	first->p1 = pts + nf[0];
	first->p2 = infinate;

	PCell* second = new PCell;
	compute_circle_alg(pts[nf[0]], pts[nf[1]], c, second->q);
	bool ccw = ccw_compute(pts[nf[0]], pts[nf[1]], c);
	if (ccw)
	{
		second->p1 = first->p1;
		second->p2 = pts + nf[1];
	}
	else
	{
		second->p1 = pts + nf[1];
		second->p2 = first->p1;
	}

	first->p1 = infinate;
	first->p2 = second->p1;
	dv_compute(*second->p1, *second->p2, c, first->q.x, first->q.y);
	first->next = second;

	PCell* third = new PCell;
	third->p1 = second->p2;
	third->p2 = infinate;
	dv_compute(*second->p2, *second->p1, c, third->q.x, third->q.y);
	second->next = third;
	third->next = first;
	last = third;
}

void Local_Regular::show_pcell(PCell* iter)
{
	if (iter == nullptr)
		return;

	PCell* start = iter;
	int i = 0;
	std::cout << "Power Diagram Cell:" << std::endl;
	do
	{
		std::cout << "num: " << i << std::endl;
		if (iter->p1 == infinate || iter->p2 == infinate)
		{
			std::cout << "方向向量：" << iter->q.x << " " << iter->q.y << std::endl;
		}
		else
		{
			std::cout << "q" << i << ": " << iter->q.x << " " << iter->q.y << std::endl;
		}
		if (iter->p1 == infinate)
			std::cout << "p1: " << "infinate" << std::endl;
		else
			std::cout << "p1: " << iter->p1->x << " " << iter->p1->y << std::endl;

		if (iter->p2 == infinate)
			std::cout << "p2: " << "infinate" << std::endl;
		else
			std::cout << "p2: " << iter->p2->x << " " << iter->p2->y << std::endl;
		iter = iter->next;
		i++;
	} while (iter != start);
	std::cout << std::endl;
}

bool Local_Regular::power_diagram_insert(Point_d* p, Point_d c, PCell* &first, PCell* &last)
{
	// p和c的权值平分线
	double A, B, C;
	weightedBisector(*p, c, A, B, C);

	// 计算权值平分线有效半区内的方向向量，用于后续的点在半区内判断
	Point_d j;
	double Aj, Bj, Cj;
	line(*p, c, Aj, Bj, Cj);
	intersection(Aj, Bj, Cj, A, B, C, j.x, j.y);
	Point_d pc{ c.x - p->x, c.y - p->y };	// 计算pc向量，此向量一定指向power diagram中c的区域
	j.x = j.x + pc.x;
	j.y = j.y + pc.y;

	PCell* iter = first;
	PCell* start = nullptr;		// 不在区域中的子链链头
	PCell* startUp = nullptr;	// 不在区域中的子链链头的上一链单元
	PCell* end = nullptr;		// 不在区域中的子链链尾
	PCell* up = last;
	int i = 0;

	// 计算power diagram交点链中不在pc权值平分线指向区域中的子链
	do
	{
		bool side;
		if (iter->p1 == infinate || iter->p2 == infinate)
			side = infinate_same_side_judge(iter, A, B, C, j, c);
		else
			side = same_side_judge(iter->q, j, A, B, C);

		if (!side && i == 1)
		{
			end = iter;
		}
		else if (!side && (i == 0 || i == 2))
		{
			start = iter;
			startUp = up;
			i++;

			if (i == 3)
				break;

			end = iter;
		}
		else if (side && i == 1)
		{
			i++;
		}

		iter = iter->next;
		up = up->next;
	} while (iter != first);

	// 中心点是否多余，是否需要继续计算
	bool no_redundant = true;

	if (start != nullptr && end != nullptr)	// 存在子链不在区域中，替换，并删除子链
	{
		if (start == end->next)		// 整个链均不在指向区域内，则此中心节点c不存在power cell，为多余点
		{
			first = last = nullptr;
			no_redundant = false;
		}
		else
		{
			no_redundant = true;

			Point_d* p1 = start->p1;
			Point_d* p2 = end->p2;

			PCell* new1 = new PCell;
			new1->p1 = p1;
			new1->p2 = p;
			if (p1 == infinate)
			{
				dv_compute(*p, *p2, c, new1->q.x, new1->q.y);
			}
			else
			{
				compute_circle_alg(*(new1->p1), *(new1->p2), c, new1->q);
			}

			PCell* new2 = new PCell;
			new2->p1 = p;
			new2->p2 = p2;
			if (p2 == infinate)
			{
				dv_compute(*p, *end->p1, c, new2->q.x, new2->q.y);
			}
			else
			{
				compute_circle_alg(*(new2->p1), *(new2->p2), c, new2->q);
			}

			new1->next = new2;
			startUp->next = new1;
			new2->next = end->next;

			/// 更新first
			first = new1;
			last = startUp;
		}

		/// 删除start到end
		end->next = nullptr;
		PCell* d = start;
		PCell* now = start;
		while (now != nullptr)
		{
			d = now->next;
			delete now;
			now = d;
		}
	}

	return no_redundant;
}

bool Local_Regular::infinate_same_side_judge(PCell* cell, double A, double B, double C, Point_d j, Point_d c)
{
	Point_d* p;
	if (cell->p1 == infinate)
		p = cell->p2;
	else
		p = cell->p1;

	double Ap, Bp, Cp;
	weightedBisector(*p, c, Ap, Bp, Cp);

	double x, y;
	intersection(A, B, C, Ap, Bp, Cp, x, y);

	Point_d k{ x + cell->q.x, y + cell->q.y };
	return same_side_judge(k, j, A, B, C);
}

void Local_Regular::show_inPCell(PCell* iter)
{
	if (iter == nullptr)
		return;

	PCell* start = iter;
	int i = 0;
	std::cout << "Power Diagram Cell ------ inPCell:" << std::endl;

	do
	{
		std::cout << "num: " << i << std::endl;
		if (iter->p1 == infinate || iter->p2 == infinate)
		{
			std::cout << "方向向量：" << iter->q.x << " " << iter->q.y << std::endl;
		}
		else
		{
			std::cout << "q" << i << ": " << iter->q.x << " " << iter->q.y << "  r: " << iter->r << std::endl;
		}

		inPCell* ipcStart = iter->ipc;
		std::cout << "inPCell points:" << std::endl;
		if(ipcStart == nullptr)
			std::cout << "NULL" << std::endl;
		else
		{
			int pn = 0;
			while (ipcStart != nullptr)
			{
				std::cout << "ipcNum: " << pn << "  x: " << ipcStart->p->x << "  y: " << ipcStart->p->y << std::endl;
				pn++;
				ipcStart = ipcStart->next;
			}
		}
		iter = iter->next;
		i++;
	} while (iter != start);
	std::cout << std::endl;
}

void Local_Regular::pcell_box(PCell* pcell, Box& box)
{
	// 计算外扩展最大权值后的半径
	double r = sqrt(pcell->q.w + wmax);
	pcell->r = r;

	double cymin = (pcell->q.y - r - ymin) / ycell;
	double cymax = (pcell->q.y + r - ymin) / ycell;
	double cxmin = (pcell->q.x - r - xmin) / xcell;
	double cxmax = (pcell->q.x + r - xmin) / xcell;

	if (cymin < 0)
		box.fy_min = 0;
	else
		box.fy_min = int(cymin);

	if (cymax >= ynum)
		box.fy_max = ynum - 1;
	else
		box.fy_max = int(cymax);

	if (cxmin < 0)
		box.fx_min = 0;
	else
		box.fx_min = int(cxmin);

	if (cxmax >= xnum)
		box.fx_max = xnum - 1;
	else
		box.fx_max = int(cxmax);
}

bool Local_Regular::in_circle_and_angle(Point_d p, PCell* cell, Point_d c)
{
	double d = sqrt(pow(p.x - cell->q.x, 2) + pow(p.y - cell->q.y, 2));
	if (d < cell->r)
	{
		bool cp1 = cross_product(*cell->p1, p, c);
		bool cp2 = cross_product(*cell->p2, p, c);
		if (cp1 != cp2)
			return true;
		else
			return false;
	}
	else
	{
		return false;
	}
}

void Local_Regular::findPoints_in_PCell(PCell* cell, Box exbox, Point_d c)
{
	Box box;
	pcell_box(cell, box);

	for (int i = box.fx_min; i <= box.fx_max; i++)
	{
		for (int j = box.fy_min; j <= box.fy_max; j++)
		{
			if (i >= exbox.fx_min && i <= exbox.fx_max && j >= exbox.fy_min && j <= exbox.fy_max)
				continue;
			else
			{
				int n = bucket[i * ynum + j];
				while (n != -1)
				{
					if (in_circle_and_angle(pts[n], cell, c))
					{
						inPCell* inpcell = new inPCell;
						inpcell->p = pts + n;
						inpcell->next = cell->ipc;
						cell->ipc = inpcell;
					}
					n = next[n];
				}
			}
		}
	}
}

Local_Regular::Local_Regular(Point_d* p, int num)
{
	pts = p;
	n = num;

	// 获取点集最大最小值
	getMinMax_ptr();

	// 依据最大最小值形成的矩形区域对节点进行分区
	partition();

	// 对点集中每个节点进行局部计算，得到节点的全局邻居连接
	for (int i = 4; i < 5; i++)
	{
		// 外扩展计算查询邻近节点，保证节点数量至少为2
		Box box;		// 外扩展查找过的网格包围盒
		std::vector<int> nf;
		near_find(pts + i, nf, box);

		std::cout << "point: " << i << std::endl;
		for (int i = 0; i < nf.size(); i++)
		{
			std::cout << "near " << i << ": " << nf[i] << std::endl;
		}
		std::cout << std::endl;

		// power diagram交点链表，环形
		PCell* first, * last;

		// 计算前两邻居点和当前中心计算节点形成的初始power diagram
		power_diagram_build_initial(nf, pts[i], first, last);
		show_pcell(first);

		// 逐个加入邻近节点，构造完整的初始power diagram
		for (int j = 2; j < nf.size(); j++)
		{
			power_diagram_insert(pts + nf[j], pts[i], first, last);
		}
		//show_pcell(first);

		// 计算初始power diagram交点链表单元PCell的扩展权值圆内节点链表
		PCell* iter = first;
		do
		{
			if (iter->p1 != infinate && iter->p2 != infinate)
			{
				findPoints_in_PCell(iter, box, pts[i]);
			}
			else
			{

			}
			iter = iter->next;
		} while (iter != first);
		show_inPCell(first);
	}
}
