# 文档编写
如果你觉得文档中说得不对、有歧义或想补充自己的想法。就需要更改这里的文档。

下面介绍更改文档的做法。

## 函数注释
参考`src/`目录下的函数注释方式，采用下面的格式装饰在函数前面：
```julia
@doc raw"""
    getCavityVoltage(z, v1, r, ϕs, ϕ2s, h1, h, circum)::Float64

这里添加函数说明。
"""
```
上面的缩进，等同于代码块：
````markdown
```julia
getCavityVoltage(z, v1, r, ϕs, ϕ2s, h1, h, circum)::Float64
```
````
编译后，会先显示代码块，再显示函数说明。

同时，`docs/src/`的目录下会有不少`.md`文件，只要在其中加上
````markdown
```@docs
getCavityVoltage
```
或者：
```@docs
DoubleRFs.getCavityVoltage
```
````
即可编译。两者区别为：当`DoubleRFs`包没有`export`该函数时，需要用`DoubleRFs.函数名`.

注意编译时，会自动将各种多重分派的函数全部收集起来一起显示。所以，如果你存在函数进行了多重分派，那么只要在`.md`文件中声明一次，就会把所有同名的代码块注释编译过来。


如何在注释中使用``\LaTeX``语法参考[Documenter的链接](https://documenter.juliadocs.org/stable/man/latex/)，语法案例参考[这里](https://documenter.juliadocs.org/stable/showcase/)。


## 编译
当完成之后，切换进`docs/`目录，在终端中：
```cmd
julia --project make.jl
```
即可编译文档。

编译完成后，在`docs/build/`目录下即可找到生成的`.html`文档。