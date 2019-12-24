# 二维非稳态导热问题的有限体积数值解法#
#from scipy import linalg
#import os

def solve(N_x_grid, N_y_grid ,t,eff1,eff2,delta_t):
	# 物性参数：发动机长度，高温辐射区域长度，两层材料的宽度、热传导系数、密度、比热容，燃气温度，初始温度,玻尔兹曼常数
	L = 0.05
	L1 = 0.02
	HH = 0.004
	hh = 0.002
	k1 = 53.6
	k2 = 0.2093
	k0 = 0.41697
	rou1 = 7830
	rou2 = 1500
	Cp1 = 465
	Cp2 = 1465
	Tair= 293
	Tgas = 923
	T0 = 293  # 初始温度
	sigma = 5.67 * 10 ** -8
	
	# 网格尺寸：
	delta_x = L / N_x_grid  # 20个网格为0.0025
	delta_y = (2 * HH + hh) / N_y_grid  # 20个网格为0.0005
	
	# 时间步长及迭代次数：
	#delta_t = 0.002#0.004*(20/N_x_grid)*(20/N_y_grid)
	cal_num = int(t / delta_t)
	
	# 系数
	aw1 = ae1 = k1 / delta_x * delta_y
	as1 = an1 = k1 / delta_y * delta_x
	aw2 = ae2 = k2 / delta_x * delta_y
	as2 = an2 = k2 / delta_y * delta_x
	as0 = an0 = k0 / delta_y * delta_x
	ap1 = rou1 * Cp1 * delta_x * delta_y / delta_t
	ap2 = rou2 * Cp1 * delta_x * delta_y / delta_t
	qgas=eff1 * sigma * Tgas**4
	#h=10
	qcon=h*Tair
	
	# 定义求解线性方程组AA*T=CC的动态数组AA、CC：
	AA = np.zeros((N_x_grid * N_y_grid, N_x_grid * N_y_grid))
	CC = np.zeros((N_x_grid * N_y_grid))
	
	# 定义初始0矩阵:
	a_w = np.zeros((N_x_grid, N_y_grid))
	a_e = np.zeros((N_x_grid, N_y_grid))
	a_n = np.zeros((N_x_grid, N_y_grid))
	a_s = np.zeros((N_x_grid, N_y_grid))
	sp = np.zeros((N_x_grid, N_y_grid))
	su = np.zeros((N_x_grid, N_y_grid))
	a_p = np.zeros((N_x_grid, N_y_grid))
	a_p0 = np.zeros((N_x_grid, N_y_grid))
	T = np.zeros((N_x_grid, N_y_grid))
	X = np.zeros((N_x_grid, N_y_grid))
	Y = np.zeros((N_x_grid, N_y_grid))
	
	# 离散方程系数计算
	for k in range(0, cal_num):
		# 给定初始温度场
		if k== 0:
			for i in range(0, N_x_grid):
				for j in range(0, N_y_grid):
					T[i][j] = T0
					Y[i][j] = delta_y * (j+1)
					X[i][j] = delta_x * (i+1)
					
		else:
			for i in range(0, N_x_grid):
				for j in range(0, N_y_grid):
					# 结点位置：
					Y[i][j] = delta_y * (j + 1)
					X[i][j] = delta_x * (i + 1)
					
					'''#对流换热系数计算：
					tm=(T[i][j]-253)/2
					beta=2/(T[i][j]+293)
					Pr=4*10**(-14)*tm**6-3*10**(-11)*tm**5+6*10**(-9)*tm**4-7*10**(-7)*tm**3+4*10**(-5)*tm**2-0.0012*tm+0.7163
					v=9*10**(-11)*tm**2+9.13*10**(-8)*tm+1.3175*10**(-5)
					Gr=9.8*beta*(T[i][j]-293)*0.03**3/v/v
					if Gr<0:
						print(i,j,Gr,T[i][j],Tair)
						break
					Nu=0.54*(Gr*Pr)**0.25
					h=2.59/3*Nu'''
					#h = 10
					
					# 1.判断结点属于哪个材料，给出a_e,a_w,a_s,a_n,a_p相应系数
					if j <= 0.4 * N_y_grid - 1 or j >= 0.6 * N_y_grid:  # 20个网格对应0——7，12——19
						a_w[i][j] = aw1  # 10.72
						a_e[i][j] = ae1  # 10.72
						a_s[i][j] = as1  # 268
						a_n[i][j] = an1  # 268
						a_p[i][j] = ap1
					else:  # 20个网格对应8——11
						a_w[i][j] = aw2
						a_e[i][j] = ae2
						a_s[i][j] = as2
						a_n[i][j] = an2
						a_p[i][j] = ap2
						
					# 2.交界面处系数
					if j == 0.6 * N_y_grid - 1 or j == 0.4 * N_y_grid - 1:  # 7，11
						a_n[i][j] = an0
					if j == 0.6 * N_y_grid or j == 0.4 * N_y_grid:  # 8，12
						a_s[i][j] = as0
					
					# 3.边界条件
					if i == 0:  # 左边界为绝热边界
						a_w[i][j] = 0.0
						su[i][j] = 0.0
						sp[i][j] = 0.0
					if i == N_x_grid - 1:  # 右边界为辐射换热边界
						a_e[i][j] = 0.0
						if j <= 0.4 * N_y_grid - 1 or j >= 0.6 * N_y_grid: #材料1辐射换热系数eff1
							su[i][j] = qgas*delta_y-eff1 * sigma * (T[i][j]) ** 4 * delta_y
						else : #材料2辐射换热系数eff2
							su[i][j] = qgas*delta_y-eff2 * sigma * (T[i][j]) ** 4 * delta_y
						sp[i][j] = 0.0
					if j == N_y_grid - 1:  # 上边界为混合边界
						if i <= N_x_grid * 0.6 - 1:  # 上边界左侧为对流换热
							a_n[i][j] = 0.0
							su[i][j] = qcon*delta_x-h * T[i][j] * delta_x
							sp[i][j] = 0
						else:  #改动1： 上边界右侧为辐射换热
							a_n[i][j] = 0.0
							su[i][j] = qgas*delta_x-eff1 * sigma * (T[i][j]) ** 4 * delta_x
							sp[i][j] = 0.0
					if j == 0:  # 下边界为绝热边界
						a_s[i][j] = 0.0
						su[i][j] = 0.0
						sp[i][j] = 0.0
						
					# 4.四个角点系数
					if i == 0 and j == 0:  # 左下角点
						su[i][j] = 0.0
					if i == 0 and j == N_y_grid - 1:  # 左上角点
						su[i][j] = qcon*delta_x-h * T[i][j] * delta_x
					if i == N_x_grid - 1 and j == N_y_grid - 1:  # 右上角点
						su[i][j] = qgas*delta_y-eff1 * sigma * (T[i][j]) ** 4 * delta_y+qgas*delta_x-eff1 * sigma * (T[i][j]) ** 4 * delta_x
					if i == N_x_grid - 1 and j == 0:  # 右下角点
						su[i][j] = qgas*delta_y-eff1 * sigma * (T[i][j]) ** 4 * delta_y
					
					# 5.最终a_p0系数
					a_p0[i][j] = a_p[i][j] - (a_w[i][j] + a_e[i][j] + a_s[i][j] + a_n[i][j] - sp[i][j])
			
			# 将系数给入系数矩阵AA,CC
					AA[i * N_y_grid + j][i * N_y_grid + j] = a_p[i][j]
					CC[i * N_y_grid + j] = su[i][j] + a_p0[i][j] * T[i][j]
					if j != 0:
						CC[i * N_y_grid + j] = CC[i * N_y_grid + j] + a_s[i][j] * T[i][j - 1]
					if j != N_y_grid - 1:
						CC[i * N_y_grid + j] = CC[i * N_y_grid + j] + a_n[i][j] * T[i][j + 1]
					if i != 0:
						CC[i * N_y_grid + j] = CC[i * N_y_grid + j] + a_w[i][j] * T[i - 1][j]
					if i != N_x_grid - 1:
						CC[i * N_y_grid + j] = CC[i * N_y_grid + j] + a_e[i][j] * T[i + 1][j]
			
			#调用scipy库求解线性方程组
			# result=np.linalg.solve(AA, CC)
			# for i in range(0, N_x_grid):
			# 	for j in range(0, N_y_grid):
			# 		if result[i * N_y_grid + j]>293:
			# 			T[i][j] = result[i * N_y_grid + j]
			result = np.linalg.solve(AA, CC)  # 求解温度场
			T = result.reshape(N_x_grid, N_x_grid)

		#if (k + 1) * delta_t % 5 == 0:
		数据写入文件(eff1, eff2, h, N_x_grid, N_y_grid, k, X, Y, T, a_p, a_n, a_s, a_w, a_e, sp, su,delta_t)
		if k % 20 == 0:
			print('已完成:{:.2f}% '.format(k / cal_num * 100))

# 调用函数的主程序
import numpy as np
from gaoti import 数据写入文件,打开tecplot_linux
#N_x_grid,N_y_grid,t,eff1,eff2=input("请输入x方向网格数量，y方向网格数量，计算时间，材料1辐射率，材料2辐射率").split()
#solve(int(N_x_grid),int(N_y_grid),int(t),float(eff1),float(eff2))#N_xgrid=20,N_ygrid=20,t=50,eff1=0.6949,eff2=0.7515
if __name__== '__main__':
	N_xgrid = 20;N_ygrid = 20;t = 30;eff1 = 0.8;eff2 = 0.3;h = 0.3;delta_t = 0.004
	#N_xgrid = 20;N_ygrid = 20; t = 30; eff1 = 0.6949; eff2 = 0.7515;h=10;delta_t=0.002
	solve (N_xgrid, N_ygrid, t, eff1, eff2,delta_t)  #
	打开tecplot_linux(N_xgrid, N_ygrid, delta_t, eff1, h, eff2)
