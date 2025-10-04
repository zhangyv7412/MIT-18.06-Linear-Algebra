#Math/LinearAlgebra 
# 微分方程
$$\large \frac{\mathrm{d}\boldsymbol{u}}{\mathrm{d}t} = A \boldsymbol{u}$$

## 一阶线性微分方程组

### 从例子开始
$$
\begin{cases}
\frac{\mathrm{d}u_1}{\mathrm{d}t} = -u_1 + 2u_2 \\ 
\frac{\mathrm{d}u_2}{\mathrm{d}t} = u_1 - 2u_2 \\
u_1(0) = 1 \\ u_2(0) = 0
\end{cases}
$$
改写成矩阵形式
$$ 
\begin{bmatrix}
	\frac{\mathrm{d}u_1}{\mathrm{d}t} \\ \frac{\mathrm{d}u_2}{\mathrm{d}t}
	\end{bmatrix}
=
\begin{bmatrix}
	-1 & 2 \\ 1 & -2 \\
	\end{bmatrix}
\begin{bmatrix}
	u_1 \\ u_2
	\end{bmatrix}
$$
$$
u(0) = 
\begin{bmatrix}
	1 \\ 0
	\end{bmatrix}
$$

因矩阵$A$是奇异矩阵,故一个特征值为$0$,又因为迹为$-3$,所以另一个特征值为$-3$.
将其记为$${\lambda}_1 = 0,{\lambda}_2 = -3$$
容易求得特征向量分别为
$$
x_1 =
\begin{bmatrix}
	2 \\ 1
	\end{bmatrix}
,x_2 =
	\begin{bmatrix}
	1 \\ -1
	\end{bmatrix}
$$
**\*我们知道解的结构应该形如**
$$\Large u = c_1u_1 + c_2u_2= c_1 e^{\lambda_1 t} x_1 +c_2 e^{\lambda_2 t} x_2$$
代入特征值和特征向量
$$u = c_1 \begin{bmatrix} 2 \\ 1 \end{bmatrix} +c_2 e^{-3 t} \begin{bmatrix} 1 \\ -1 \end{bmatrix}$$
代入初始条件
$$
\begin{equation}
u(0) = c_1 \begin{bmatrix} 2 \\ 1 \end{bmatrix} +c_2 \begin{bmatrix} 1 \\ -1 \end{bmatrix} = \begin{bmatrix} 1 \\ 0 \end{bmatrix}
\iff u(0) = \begin{bmatrix} 2 & 1 \\ 1 & -2 \end{bmatrix} \begin{bmatrix} c_1 \\ c_2 \end{bmatrix} = \begin{bmatrix} 1 \\ 0 \end{bmatrix}
\end{equation}
$$
容易解得$c_1 = c_2 = \frac{1}{3}$,最终答案为$$u = \frac{1}{3} \begin{bmatrix} 2 \\ 1 \end{bmatrix} + \frac{1}{3}e^{-3t} \begin{bmatrix} 1 \\ -1 \end{bmatrix}$$

### 一般性的推导
让我们沿着例子[[微分方程和exp(At)#一阶线性微分方程组]]中的过程,得到更一般的结果.
下面我们来求解
$$\large \frac{\mathrm{d}\boldsymbol{u}}{\mathrm{d}t} = A \boldsymbol{u}$$
其中,矩阵$A$说明$u$是耦合的,关键在于如何解耦.解耦在此处就是对角化.
令
$$u = Sv$$
代入得
$$S \frac{\mathrm{d}v}{\mathrm{d}t} = ASv$$
$$\frac{\mathrm{d}v}{\mathrm{d}t} = S^{-1}ASv = \Lambda v$$
>$S^{-1}AS = \Lambda$,详见[[对角化]]]

得到解耦后的结果
$$\frac{\mathrm{d}v}{\mathrm{d}t} = \Lambda v$$
解得
$$v(t) = e^{\Lambda t}v(0)$$

>根据[[微分方程和exp(At)#指数矩阵#定义]],只需要把每个矩阵加起来,对角线各元素位置就是$e^x$的级数展开式,所以有
>$$e^{\Lambda t} = 
>\begin{bmatrix}
>	e^{\lambda_1 t} & 0 & \cdots & 0 \\
>	0 & e^{\lambda_2 t} & \cdots & 0 \\
>	\vdots & \vdots & \ddots & \vdots \\
>	0 & 0 & \cdots & e^{\lambda_n t} \\
>	\end{bmatrix}
>$$

代入$u = Sv$,得
$$u = Se^{\Lambda t}v(0) = SS^{-1}e^{At}Sv(0) = e^{At}u(0)$$

>此处使用等式 $e^{At}= Se^{\Lambda t}S^{-1}$
>$$
>\begin{align}
>e^{At} &= I + At + \frac{(At)^2}{2} + \cdots \\
>&= SS^{-1} + S \Lambda S^{-1}t + \frac{(S \Lambda S^{-1}t)^2}{2} + \cdots \\
>&= SS^{-1} + S \Lambda S^{-1}t + \frac{S (\Lambda)^2 S^{-1}t^2}{2} + \cdots \\
>&= S(I + \Lambda t + \frac{\Lambda^2 t^2}{2} + \cdots)S^{-1} \\ 
>&= Se^{\Lambda t}S^{-1}
>\end{align}
>$$
>>$A = S \Lambda S^{-1}$,详见[[对角化]]

## 二阶线性微分方程
>[!note] 关于二阶线性微分方程的求解
>$$y'' + by'  + ky = 0$$
>>[!tip] 主要的想法
>>**把二阶方程转换为一阶方程组**.

只需要给原方程添加一个等式,得到的方程组就可以用矩阵表示.
$$
\begin{align} &\begin{cases} y'' + by' + ky &= 0 \\ y' &= y' \end{cases} \iff \begin{cases} y'' &= -by' - ky \\ y' &= y' \end{cases} \end{align}
$$
令
$$u = \begin{bmatrix}
y' \\ y 
\end{bmatrix}$$
方程组可以改写成
$$\begin{bmatrix}
y'' \\ y'
\end{bmatrix}
= 
\begin{bmatrix}
-b & -k \\ 1 & 0
\end{bmatrix}
\begin{bmatrix}
y' \\ y
\end{bmatrix}
\iff
u' = Au
$$
下面的计算就和[[微分方程和exp(At)#一阶线性微分方程组]]中所提到的一样
>[!summary] 总结
>大多数情况下,对高阶/高维形式的问题处理总是想办法将其降阶/降维,然后从简单的情形入手解决问题.






# 指数矩阵
## 定义
正如定义指数函数是通过幂级数展开的方式,我们在此处也通过**幂级数**来定义指数矩阵.
>$eg1.$指数函数
>$$e^x = \sum_{n = 0}^{\infty} \frac{x^n}{n!}$$

仿照$e^x$的幂级数展开,$1$替换为单位阵$I$,$x$替换为$At$,得到如下定义式:
$$\large e^{At} = I + At + \frac{(At)^2}{2} + \cdots = \sum_{n = 0}^{\infty} \frac{(At)^n}{n!}$$
>$eg2.$几何级数
>$$\frac{1}{1-x} = \sum_{n = 0}^{\infty} x^n$$

同样的,可以得到:
$$\large (I - At)^{-1} = I + At + (At)^2 + \cdots = \sum_{n = 0}^{\infty} (At)^n$$
