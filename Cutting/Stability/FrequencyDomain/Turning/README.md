## Calculation of the phase angle $\psi$ in polar plots

注意：这里将车削条件设为正交插入车削，即重叠系数视为$1$.

<div align = "center">

<img src = "Orthogonal plunge turning.png"  width = "500" height = "400" alt = "Orthogonal plunge turning" title = "Orthogonal plunge turning">

</div>

<p align = "center"><b>Figure 1.  Orthogonal plunge turning</b> </p>

切削厚度为：
$$h\left( t \right) = {h _0} - \left[ {y\left( t \right) - y\left( {t - T} \right)} \right]$$

系统的运动方程可以写为：
$$\begin{gathered}
  {m _y}\ddot y\left( t \right) + {c _y}\dot y\left( t \right) + {k _y}y\left( t \right) = {F _f}\left( t \right) = {K_f}ah\left( t \right)\\\\
   = {K_j}a\left[ {{h_0} + y\left( {t - T} \right) - y\left( t \right)} \right]
\end{gathered} $$
将时变的切削厚度转换到拉式域：
$$h\left( s \right) = {h_0} - y\left( s \right) + {e^{ - sT}}y\left( s \right) = {h_0} + \left( {{e^{ - sT}} - 1} \right)y\left( s \right)$$
时变的切削力为：
$${F_f}\left( s \right) = {K_f}ah\left( s \right)$$
时变切削力激励结构产生振动：
$$y\left( s \right) = {F _f}\left( s \right)\Phi \left( s \right) = {K _f}ah\left( s \right)\Phi \left( s \right)$$
代入到时变切削厚度中有：
$$h\left( s \right) = {h _0} + \left( {{e^{ - sT}} - 1} \right){K _f}ah\left( s \right)\Phi \left( s \right)$$
那么得到实际切削厚度与名义切削厚度之间的传递函数：
$$\frac{{h\left( s \right)}}{{{h_0}}} = \frac{1}{{1 + \left( {1 - {e^{ - sT}}} \right){K_f}a\Phi \left( s \right)}}$$
如图$1$所示闭环系统。那么，该闭环系统的稳定性由传递函数特征多项式（即传递函数的分母）的根决定。
令特征多项式为$0$：
$$1 + \left( {1 - {e^{ - j{\omega _c}T}}} \right){K _f}{a _{\lim }}\Phi \left( {j{\omega _c}} \right) = 0$$

由于传递函数为复值函数，令  。
上式可展开为：
$$ \lbrace {1 + {K _f}{a _{\lim }}\left[ {G\left( {1 - \cos {\omega _c}T} \right) - H\sin {\omega _c}T} \right]}  \rbrace + j \lbrace {{K _f}{a _{\lim }}\left[ {G\sin {\omega _c}T + H\left( {1 - \cos {\omega _c}T} \right)} \right]}  \rbrace = 0$$
实部虚部分别都得等于零，对于虚部有：
$$G\sin {\omega _c}T + H\left( {1 - \cos {\omega _c}T} \right) = 0$$
那么频响函数相位满足：
$$\tan \psi  = \frac{{H\left( {{\omega _c}} \right)}}{{G\left( {{\omega _c}} \right)}} = \frac{{\sin {\omega _c}T}}{{\cos {\omega _c}T - 1}} =  - \frac{{\cos \left( {{\omega _c}\frac{T}{2}} \right)}}{{\sin \left( {{\omega _c}\frac{T}{2}} \right)}} = \tan \left[ {\frac{{\left( {{\omega _c}T} \right)}}{2} - \frac{{\left( {3\pi } \right)}}{2}} \right]$$
这里需要注意的一点是这个$ - \frac{{3\pi }}{2}$ ，一般来说，都是取$ - \cot \alpha  = \tan \left( {\alpha  - \frac{\pi }{2}} \right)$ ，这样$\alpha$可以取$(0,90^\circ)$锐角。但试验发现，如果这样取，每阶包含负实部的模态都会引入一个多余的负频率解。
由上面的公式可以看出：
$$\begin{gathered}
  \psi  + k\pi  = \frac{{\left( {{\omega _c}T} \right)}}{2} - \frac{{\left( {3\pi } \right)}}{2}  \\\\
  {\omega _c}T = 3\pi  + 2\psi  + 2k\pi {\kern 1pt} \left( {k = 0,1,2...} \right)  \\
\end{gathered} $$
即，$\psi $  值的大小会影响${\omega _c}$的正负，需要保证频率为正。
$$\begin{gathered}
  {f_c}T = k + \frac{\varepsilon }{{2\pi }}\left( {k = 0,1,2...} \right) \\\\
  2\pi {f_c}T = 2\pi k + \varepsilon  \\\\
  {\omega _c}T = 2\pi k + \varepsilon
\end{gathered} $$
由此可知：
$$\varepsilon  = 3\pi  + 2\psi $$
如果上面不取$ - \frac{{3\pi }}{2}$而取其它值，则会造成负频率（大于$ - \frac{{3\pi }}{2}$时），或低阶$Lobe$缺失（小于$ - \frac{{3\pi }}{2}$时）。
这是因为，$\psi $  值只能在$(-180^{\circ} ,-90^\circ)$之间变化。

<div align = "center">

<img src = "Polar Plots.png"  width = "400" height = "400" alt = "phase angle $\psi$ in polar plots" title = "phase angle $\psi$ in polar plots">

</div>

<p align = "center"><b>Figure 2.  phase angle $\psi$ in polar plots</b> </p>

一方面可以从特征方程的实部看出：

$$1 + {K_f}{a_{\lim }}\left[ {G\left( {1 - \cos {\omega _c}T} \right) - H\sin {\omega _c}T} \right] = 0$$

$$\begin{gathered}
  {a_{\lim }} = \frac{{ - 1}}{{{K_f}{a_{\lim }}\left[ {G\left( {1 - \cos {\omega _c}T} \right) - H\sin {\omega _c}T} \right]}} = \frac{{ - 1}}{{{K_f}{a_{\lim }}G\left[ {\left( {1 - \cos {\omega _c}T} \right) - \frac{H}{G}\sin {\omega _c}T} \right]}} \\\\
   = \frac{{ - 1}}{{2{K_f}G\left( {{\omega _c}} \right)}} \\\\
\end{gathered} $$

由于，轴向切深为正值，所以，如果颤振发生，那么发生颤振的频率 ${\omega _c}$ 对应的频响的实部必须<b>为负</b>。
$\psi $  角只能在二、三象限。
$\psi \in (-180^{\circ} ,-90^\circ)$，所以有$\psi + {\pi } \in (0^{\circ} ,90^\circ)$
把这个 $\pi$ 移到等号另一边，这样就能解释上面的 $ - \frac{{3\pi }}{2}$ 了。

另一个角度也可以说明，为什么颤振频率比固有频率要高，因为固有频率处实部为零，而颤振频率为实部为负的地方，显然此处的频率要大于该阶模态对应的固有频率。

从机械这个二阶系统来说，对于实部为负时，其相位只能在第二、三象限， 也能解释得通。