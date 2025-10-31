下面用HEPS参数演示用法。

首先，为了方便，先将加速器参数打包进一个结构体中：
```julia
par = ParDBRF(3.6759e6, 0.1908, 2.0391, 5.6879795185567140, 756, 3, 1360.4, 6e9, 1.56e-5, 0, 106.27, 1.06e-3)
```
这里的各个数值参考函数`ParDBRF`的变量含义——在左边搜索框搜索即可。

这里参数实际有的时候是不必要的，比如你只想要纵向密度分布，$\nu_y$是没有必要的，随便填点就可以。

---

# 双高频纵向密度分布

```julia
zvec, ρvec = getNormedDensity(par)
```
程序会先计算出相稳定区的左右边界，然后在其中均分`10_0001`个点，并计算相应归一化密度分布。

该程序无法考虑纵向阻抗。


!!! note "需要时刻保证纵向密度分布是否光滑"
    这里选用`10_0001`个点，有点多。因为存在一些极端情况，粒子只占相稳定区中间极少部分。这种情况下，如果粒子数仅取1000个，则相应曲线就不光滑，甚至都看不到峰值。

    这`10_0001`个点是我写死的。不敢保证你的参数下，这10w个点是否足够。所以，最好画出来观察。

---

# 计算束长

```julia
σz = getRMSWidth(zvec, ρvec)
```
计算$\rho(z)$的RMS宽度。


---

# 双高频系统纵向工作点分布
```julia
zvec, νsvec = getSynchrotronTunes(numPoints, par)
```
根据加速器参数`par`，以及采样点数`numPoints`计算纵向位置`zvec`以及相应位置的工作点`νsvec`。

!!! note "注意画出来观察一下"
    不知道什么原因，有时程序在个别的点，会算出两倍的频率结果。
    
    我知道这个结果是错的，所以我通常会手动降低采样点数，这通常能解决这个问题。

!!! note "这一段对于分析不稳定性是不必要的"
    如果你只是想分析双高频系统不稳定性，没必要运行这一段。这一段是用来单独观察纵向工作点分布$\nu_s(z)$的。双高频不稳定性的代码中，有单独的通过正则变换计算方法得到$\nu_s(J)$。

    而且这里`zvec`和前面纵向密度分布变量重名，你不注意就会覆盖掉前面的`zvec`并产生冲突。如果你非要在同一个代码中运行，请注意更改变量名。

---

# 双高频系统不稳定性分析

计算不稳定性，需要先对相空间采样，再根据采样点分析不稳定性。


计算采样点有两种方法：
1. 一种是我常用的，纯双高频系统、没有势阱畸变的。
2. 一种我没怎么用过，原本是为了纵向势阱畸变留的方案：根据具体potential进行计算。

## 纯双高频系统
如果是纯双高频系统，则可直接用下面的方式获得相空间采样：
```julia
EJ, JJ, ΔJJ, zJθ, δJθ, ΦJθ, νsJ, ψJ, nJs = getDatas(200, 199, 3σz, zvec, ρvec, par)
```
画出来的采样最好再观察一下，特别是过拉伸，哪部分是头部采样点、哪部分是尾部、有没有漏掉的，最好看看。

得到采样之后，可以利用采样数据分析不稳定性。下面的程序展示的`0~0.5`nC的不稳定性结果。
```julia
rvec, gvec, tvec=scanNb(collect(0:0.02:1), 0.5 * 1e-9/elementcharge, JJ, ΔJJ, ψJ, νsJ, zJθ, ΦJθ, ωvec, imp_trans, 3, par)
```
这里`rvec`是归一化的电荷量，表示的是`N/Nb`；`gvec`表示每秒增长率；`tvec`表示频移；`elementcharge`表示元电荷量，需要自己定义`const elementcharge = 1.6e-19`。

!!! note "阻抗的处理"
    这里采用的阻抗是全频域的阻抗，频率单位是`rad/s`。通常ELEGANT导出的阻抗是正频域的，频率单位是`GHz`。所以需要统一下格式，直接算肯定是错误的。

    如果是ELEGANT导出的阻抗，可以用下面的代码处理一下：
    ```julia
    using CSV, DataFrames

    function read_impedance(path::String=raw"./ZTransverse_f_Re_Im.csv")
        df=CSV.read(path, DataFrame, header=["f[GHz]","Re","Im"])
        tmpωvec = 1e9 * 2π * df[!, "f[GHz]"]
        tmpreal = df[!, "Re"]
        tmpimag = df[!,"Im"]
        ωvec = vcat(-tmpωvec[end:-1:1], tmpωvec)
        reZ = vcat(-tmpreal[end:-1:1], tmpreal)
        imZ = vcat(tmpimag[end:-1:1], tmpimag)
        ωvec, reZ .+ 1.0im .* imZ
    end
    ```

    这里没有`CSV`和`DataFrames`包的话，需要先安装一下。这两个一个是读取csv文件的包，一个是处理DataFrame的包。



## 双高频系统+势阱畸变
如果包含势阱畸变，先要求出准确的`zvec`, `ρvec`，根据$\rho(z)$，计算直角坐标系中$(z,\delta)$相空间的密度分布$\psi(z,\delta)$：
```julia
funψzδ(z, δ)=LinearInterpolation(zvec, ρvec)(z)*1/(sqrt(2π)*par.σδ)*exp(-0.5*(δ/par.σδ)^2)
```
还要计算势能的函数`func_potential`：
```julia
pvec = getPotenCavity.(zvec, par) + 畸变势阱
func_potential = Spline1D(zvec, pvec)

EJ, JJ, ΔJJ, zJθ, δJθ, ΦJθ, νsJ, ψJ, nJs = getDatas(200, 199, 3σz, funψzδ, func_potential, -0.0075, 0.0075, [-0.002740890061896845, 0.0027408900618968233], [0.0], 1360.4, par.η)
```
得到采样之后，和上面一样的方法分析不稳定性：
```julia
rvec, gvec, tvec=scanNb(collect(0:0.02:1), 0.5 * 1e-9/elementcharge, JJ, ΔJJ, ψJ, νsJ, zJθ, ΦJθ, ωvec, imp_trans, 3, par)
```

!!! note "要检查的地方"
    一个是要检查纵向密度分布$\rho(z)$是否正确、是否光滑。不光滑会导致径向密度分布$\psi(J)$的错误。

    另一个要检查的是，径向密度分布$\psi(J)$的积分：$2 \pi \int_0^\infty \psi(J) dJ$是否接近1。
