# Library

## 目录
```@contents
Pages = ["api.md"]
```


## 导出函数
```@docs
getDatas
```

```@docs
scanNb
```



## 参数结构体
```@docs
ParDBRF
```



## Basic Equations
### 电压
The total voltage is:
```math
V_{tot}(z) =V_1 \mathcal{V}(z) - 4 \pi \epsilon_0 N_b e \int_z^\infty dz' \rho(z') W_0'(z-z'),
```
where
```math
\mathcal{V}(z)=\sin{\left( \phi_s-\dfrac{2h_1\pi}{C}z \right)}-\sin{\phi_s}+r\sin{\left[\phi_{2s}-\dfrac{2h_1 h \pi}{C}z \right]}-r\sin{\phi_{2s}}
```
, called `normed voltage`. The first term of RHS of ``V_{tot}(z)`` is cavity voltage and second term is wake voltage.


```@docs
getCavityVoltage
```

### 电压导数
```@docs
DoubleRFs.∇normedVoltage
```



### Hamiltonian
我们这里采取的哈密顿量形式为：
```math
H(z,\delta;s)=-\dfrac{\eta_p}{2}\delta^2-\dfrac{eV_1}{2\pi h_1 E}\mathcal{P}(z)+\dfrac{4 \pi \epsilon_0 N_0 r_e}{\gamma C}\int_{0}^{z}dz'' \int_{z''}^\infty dz' \rho(z') W_0'(z''-z'),
```
这里：
```math
\mathcal{P}(z)= \cos{\left( \phi_s - \dfrac{2 \pi h_1}{C}z \right)}+\dfrac{r}{h} \cos{\left( \phi_{2s}-\dfrac{2\pi h_2}{C}z \right)}+\left(-\dfrac{2\pi h_1}{C}z\right)\left(\sin{\phi_s}+r\sin{\phi_{2s}}\right)-\cos{\phi_s}-\dfrac{r}{h}\cos{\phi_{2s}},
```

the second and the last term in RHS of ``H(z,\delta;s)`` is our "cavity potential" and "wake potential".(see `getPotenCavity` and `getPotenWake`)

The Hamiltonian with this form ``H(z, \delta; s)`` looks like an inverted bowl.

```@docs
hamiltonian
```


### Density Distribution
容易发现，
```math
\exp{\left( \frac{H(z, \delta; s)}{\eta_p \sigma_\delta^2} \right)}
```
关于``\delta``是纯高斯的。这很符合一种称为Haissinski分布。分解去``\delta``剩下的关于``z``的就是纵向分布``\rho(z)``。






## Potential Term in Hamiltonian


## Hamiltonian

## Density Distribution


## Synchrotron Tunes
```@docs
getBucketBounds
```
