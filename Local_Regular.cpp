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

			// 更新链中的inPCell
			end->next = nullptr;
			update_inPCell(start, end, c, new1, new2);
		}

		/// 删除start到end
		end->next = nullptr;
		PCell* d = start;
		PCell* now = start;
		while (now != nullptr)
		{
			// 删除ipc
			inPCell* inp = now->ipc;
			while (inp != nullptr)
			{
				inPCell* t = inp->next;
				delete inp;
				inp = t;
			}

			// 删除cell
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

void Local_Regular::pcell_bound_infinate(PCell* cell, Point_d c, std::vector<Bound>& vb)
{
	// 查找非infinate点
	Point_d t;
	if (cell->p1 == infinate)
		t = *(cell->p2);
	else
		t = *(cell->p1);

	// 计算外扩展直线
	double A, B, C;
	line(t, c, A, B, C);
	double n = sqrt(wmax);	// 移动距离
	if (A * cell->q.x > 0)	// 当前直线的法向量和power diagram方向向量同向，则外扩方向为直线法向量向后
	{
		C = C + n * sqrt(A * A + B * B);
	}
	else  // 当前直线的法向量和power diagram方向向量反向，则外扩方向为直线法向量
	{
		C = C - n * sqrt(A * A + B * B);
	}

	// 计算包围范围
	if (B == 0)		// 外扩展直线垂直于x轴，单个Bound即可表示所有范围
	{
		double x = -C / A;
		int x_id = int((x - xmin) / xcell);
		if (x_id >= 0 && x_id < xnum)
		{
			Bound b;
			if (cell->q.x > 0)
			{
				b.x_min = x_id;
				b.x_max = xnum - 1;
			}
			else
			{
				b.x_min = 0;
				b.x_max = x_id;
			}
			b.y_min = 0;
			b.y_max = ynum - 1;
			vb.push_back(b);
		}
	}
	else     // 外扩展直线非垂直于x轴
	{
		double k = -A / B;	// 直线斜率
		double change;	// 查询方向y的变化数值

		// 确保有效区域向上时沿着直线y值增大的方向进行，有效区域向下时，沿着y值减小的方向进行
		int for_start, for_end, step;
		if (cell->q.y > 0)		// 有效区域向上，沿着直线y值增大的方向
		{
			if (k >= 0)
			{
				for_start = 0;
				for_end = xnum;
				step = 1;
			}
			else
			{
				for_start = xnum - 1;
				for_end = 1;	// 循环遍历时采用i * step与for_end进行判断，正常应该是-1，乘以step后使用1
				step = -1;
			}
			change = abs(k) * xcell;
		}
		else		// 有效区域向下，沿着直线y值减小的方向
		{
			if (k >= 0)
			{
				for_start = xnum - 1;
				for_end = 1;
				step = -1;
			}
			else
			{
				for_start = 0;
				for_end = xnum;
				step = 1;
			}
			change = -abs(k) * xcell;
		}

		double y_1d;
		int y_1;
		
		if (for_start == 0)
		{
			y_1d = -(A * xmin + C) / B;
			y_1 = int((y_1d - ymin) / ycell);
		}
		else
		{
			y_1d = -(A * (xmin + xnum * xcell) + C) / B;
			y_1 = int((y_1d - ymin) / ycell);
		}
		
		for (int i = for_start; i * step < for_end; i = i + step)
		{
			double y_2d = y_1d + change;
			int y_2 = int((y_2d - ymin) / ycell);

			Bound b;
			int y_min, y_max;
			if (cell->q.y > 0)	// 有效区域向上，y_2 > y_1
			{
				if (y_1 >= ynum)
					break;
				else if (y_1 < 0)
					y_min = 0;
				else
					y_min = y_1;
				y_max = ynum - 1;
			}
			else  // 有效区域向下，y_2 < y_1
			{
				if (y_1 < 0)
					break;
				else if (y_1 >= ynum)
					y_max = ynum - 1;
				else
					y_max = y_1;
				y_min = 0;
			}

			b.x_min = i;
			b.x_max = i;
			b.y_min = y_min;
			b.y_max = y_max;
			vb.push_back(b);

			y_1 = y_2;
			y_1d = y_2d;
		}
	}
}

void Local_Regular::findPoints_in_PCell_infinate(PCell* cell, Box exbox, Point_d c)
{
	std::vector<Bound> vb;
	pcell_bound_infinate(cell, c, vb);

	// vb显示代码
	//std::cout << "infinate Bound show: " << std::endl;
	//if (cell->p1 == infinate)
	//{
	//	std::cout << "p1: infinate" << std::endl;
	//	std::cout << "p2: " << cell->p2->x << " " << cell->p2->y << std::endl;
	//}
	//else
	//{
	//	std::cout << "p1: " << cell->p1->x << " " << cell->p1->y << std::endl;
	//	std::cout << "p2: infinate" << std::endl;
	//}
	//for (int i = 0; i < vb.size(); i++)
	//{
	//	std::cout << "Bound xmin: " << vb[i].x_min << "  xmax: " << vb[i].x_max << "  ymin: " << vb[i].y_min << "  ymax: " << vb[i].y_max << std::endl;
	//}
	//std::cout << std::endl;

	for (int k = 0; k < vb.size(); k++)
	{
		for (int i = vb[k].x_min; i <= vb[k].x_max; i++)
		{
			for (int j = vb[k].y_min; j <= vb[k].y_max; j++)
			{
				if (i >= exbox.fx_min && i <= exbox.fx_max && j >= exbox.fy_min && j <= exbox.fy_max)
					continue;
				else
				{
					int n = bucket[i * ynum + j];
					while (n != -1)
					{
						inPCell* inpcell = new inPCell;
						inpcell->p = pts + n;
						inpcell->next = cell->ipc;
						cell->ipc = inpcell;
						n = next[n];
					}
				}
			}
		}
	}
}

void Local_Regular::update_inPCell(PCell* start, PCell* end, Point_d c, PCell* new1, PCell* new2)
{
	double A1, B1, C1, A2, B2, C2;
	Point_d p1, p2;
	bool new1_infinate, new2_infinate;

	// 判断new节点中是否有infinate，有则需进行半平面判断，无则进行圆内判断
	if (new1->p1 != infinate && new1->p2 != infinate)
	{
		new1->r = sqrt(new1->q.w + wmax);
		new1_infinate = false;
	}
	else
	{
		compute_ext_line_and_point(new1, c, A1, B1, C1, p1);
		new1_infinate = true;
	}

	if (new2->p1 != infinate && new2->p2 != infinate)
	{
		new2->r = sqrt(new2->q.w + wmax);
		new2_infinate = false;
	}
	else
	{
		compute_ext_line_and_point(new2, c, A2, B2, C2, p2);
		new2_infinate = true;
	}

	PCell* s = start;
	while (s != nullptr)	// start end链判断
	{
		inPCell* ip = s->ipc;
		while (ip != nullptr)	// 链上节点内部ipc节点链判断
		{
			inPCell* next = ip->next;
			bool in1, in2;
			if (new1_infinate)	// 有无穷点，进行半平面判断
			{
				in1 = same_side_judge(p1, *(ip->p), A1, B1, C1);
			}
			else    // 无无穷点，进行园内及夹角内判断
			{
				in1 = in_circle_and_angle(*(ip->p), new1, c);
			}

			if (in1)
			{
				inPCell* t;
				t = new1->ipc;
				new1->ipc = ip;
				ip->next = t;
			}
			else
			{
				if (new2_infinate)	// 有无穷点，进行半平面判断
				{
					in2 = same_side_judge(p2, *(ip->p), A2, B2, C2);
				}
				else    // 无无穷点，进行园内及夹角内判断
				{
					in2 = in_circle_and_angle(*(ip->p), new2, c);
				}

				if (in2)
				{
					inPCell* t;
					t = new2->ipc;
					new2->ipc = ip;
					ip->next = t;
				}
				else   // 不在new1或new2中，则删除当前inPCell
				{
					delete ip;
				}
			}

			ip = next;
		}

		s = s->next;
	}
}

void Local_Regular::compute_ext_line_and_point(PCell* cell, Point_d c, double& A, double& B, double& C, Point_d& p)
{
	// 查找非infinate点
	Point_d t;
	if (cell->p1 == infinate)
		t = *(cell->p2);
	else
		t = *(cell->p1);

	// 计算外扩展直线
	A, B, C;
	line(t, c, A, B, C);
	double n = sqrt(wmax);	// 移动距离
	if (A * cell->q.x > 0)	// 当前直线的法向量和power diagram方向向量同向，则外扩方向为直线法向量向后
	{
		C = C + n * sqrt(A * A + B * B);
	}
	else  // 当前直线的法向量和power diagram方向向量反向，则外扩方向为直线法向量
	{
		C = C - n * sqrt(A * A + B * B);
	}

	p.x = t.x + cell->q.x;
	p.y = t.y + cell->q.y;
}

Local_Regular::Local_Regular(Point_d* p, int num)
{
	pts = p;
	n = num;

	nei.resize(n);

	// 获取点集最大最小值
	getMinMax_ptr();

	// 依据最大最小值形成的矩形区域对节点进行分区
	partition();

	// 对点集中每个节点进行局部计算，得到节点的全局邻居连接
	for (int i = 0; i < n; i++)
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

		// 当前中心点是否为多余点
		bool valid = true;
		// 逐个加入邻近节点，构造完整的初始power diagram
		for (int j = 2; j < nf.size(); j++)
		{
			valid = power_diagram_insert(pts + nf[j], pts[i], first, last);
			if (!valid)
				break;
		}
		//show_pcell(first);

		// 当前中心点非有效节点，当前中心节点的邻域连接点为NULL，无需后续计算步骤
		if (!valid)
			continue;

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
				findPoints_in_PCell_infinate(iter, box, pts[i]);
			}
			iter = iter->next;
		} while (iter != first);
		//show_inPCell(first);

		// 循环遍历power diagram链，不断加入节点并更新范围内节点链表，直到无节点包含
		bool has_point = true;		// 是否存在power diagram交点包含圆内节点，存在则需进行下一轮循环，无则计算完成
		while (has_point)		// 每次循环遍历完整链一次
		{
			has_point = false;
			iter = first;
			do
			{
				Point_d* p;
				// 存在圆内点，取出首点，修改ipc链表，修改has_point为true，需进行下次循环判断
				if (iter->ipc != nullptr)
				{
					has_point = true;
					p = iter->ipc->p;
					inPCell* del = iter->ipc;
					iter->ipc = iter->ipc->next;
					delete del;
				}
				else
				{
					iter = iter->next;
					continue;
				}

				valid = power_diagram_insert(p, pts[i], first, last);
				if (!valid)
				{
					has_point = false;
					break;
				}

				iter = iter->next;
			} while (iter != first);
		}

		if (!valid)		// 当前中心点为多余点，无需操作
			continue;
		else			// 当前中心点非多余，统一保存邻居连接节点
		{
			last->next = nullptr;
			PCell* pc = first;
			while (pc != nullptr)
			{
				Nei* e = new Nei;
				e->p = pc->p1;
				e->next = nei[i];
				nei[i] = e;
				
				pc = pc->next;
			}
		}

		std::cout << "邻居连接点：" << std::endl;
		//std::cout << "num: " << i << std::endl;
		int m = 0;
		Nei* ni = nei[i];
		while (ni != nullptr)
		{
			std::cout << "seq: " << m << ":  x: " << ni->p->x << "  y: " << ni->p->y << std::endl;
			ni = ni->next;
		}
		std::cout << std::endl;
	}
}
