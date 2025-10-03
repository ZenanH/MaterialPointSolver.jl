# 粒子更新流程（每时间步）

---

## 预备输入（步 $n$ 已知）

- $\mathbf F_n$ （存储为 9 元展开，使用时还原为 $3\times3$）
- 若走 **速率式压力**：$p_n$
- 物性与数值参数：
  - $A(T)$  
  - 幂指数 $n$  
  - 正则化 $\dot\varepsilon_0$  
  - 体模量 $K_s$  
  - 时间步 $\Delta t$  
  - 初始体积 $V_{p0}$ 或初始密度 $\rho_0$
- 从网格插值得到的粒子速度梯度：  

  $$
  \mathbf L = \nabla\mathbf v = \sum_i \mathbf v_i \otimes \nabla N_i(\mathbf x_p)
  $$

---

## 计算顺序（统一到步 $n \to n{+}1$）

### 1. 分解速度梯度
$$
\mathbf D = \tfrac12(\mathbf L+\mathbf L^\mathsf T), \qquad
\mathbf W = \tfrac12(\mathbf L-\mathbf L^\mathsf T)
$$

（后续只用 $\mathbf D$）

---

### 2. 更新变形梯度与体积比
$$
\mathbf F_{n+1}=(\mathbf I+\Delta t\,\mathbf L)\,\mathbf F_n, \qquad
J_{n+1}=\det\mathbf F_{n+1}
$$

并得：
$$
V_{p}^{\,n+1}=J_{n+1}V_{p0}, \qquad
\rho_{n+1}=\rho_0/J_{n+1}
$$

---

### 3. 体变与偏变形率
$$
\operatorname{tr}\mathbf D, \qquad
\mathbf D'=\mathbf D-\tfrac13(\operatorname{tr}\mathbf D)\mathbf I
$$

---

### 4. 等效应变率（带正则）
$$
\dot\varepsilon_\mathrm{e}=\sqrt{\tfrac12\,\mathbf D':\mathbf D'}, \qquad
\tilde{\dot\varepsilon}_\mathrm{e}=\sqrt{\dot\varepsilon_0^2+\dot\varepsilon_\mathrm{e}^2}
$$

---

### 5. 等效黏度（Glen–Nye，应变率式）
$$
\eta=\tfrac12\,A(T)^{-1/n}\,\tilde{\dot\varepsilon}_\mathrm{e}^{\,(1-n)/n}
$$

（必要时对 $\eta$ 设上下限以稳态）

---

### 6. 偏应力（Cauchy 的偏部分）
$$
\boldsymbol{\tau}=2\,\eta\,\mathbf D'
$$

---

### 7. 压力闭合（两选一）

**A. 速率式（推荐，易与显式步进配合）**

$$
p_{n+1}=p_n- K_s\,\Delta t\,\operatorname{tr}\mathbf D
$$

**B. 由体积比（能量一致、无漂移）**

$$
p_{n+1}=-K_s\,\ln J_{n+1}
$$

---

### 8. 总 Cauchy 应力
$$
\boldsymbol{\sigma}_{n+1}=-p_{n+1}\mathbf I+\boldsymbol{\tau}
$$

六分量展开：

$$
(\sigma_{xx},\sigma_{yy},\sigma_{zz},\sigma_{xy},\sigma_{yz},\sigma_{zx})
=(\sigma_{11},\sigma_{22},\sigma_{33},\sigma_{12},\sigma_{23},\sigma_{31})
$$

---

### 9. （用于内力装配）粒子体积与内力
$$
V_p^{\,n+1}=J_{n+1}V_{p0}, \qquad
\mathbf f_{\text{int}}^{\,i}=-\,V_p^{\,n+1}\,\boldsymbol{\sigma}_{n+1}\,\nabla N_i(\mathbf x_p^{\,n+1})
$$

---

## 小结

$$
\boxed{
\mathbf L \;\Rightarrow\; \mathbf D,\mathbf W \;\Rightarrow\; \mathbf F_{n+1},J_{n+1}
\;\Rightarrow\; \mathbf D',\dot\varepsilon_\mathrm{e} \;\Rightarrow\; \eta
\;\Rightarrow\; \boldsymbol{\tau} \;\Rightarrow\; p_{n+1}
\;\Rightarrow\; \boldsymbol{\sigma}_{n+1}
}
$$

并同步更新 $V_p^{\,n+1}$（或 $\rho_{n+1}$）以用于内力与后续步。

> 采用此流程，纯刚体转动时 $\mathbf D=0 \Rightarrow \boldsymbol{\tau}=0$，天然客观；  
> 体积变化主要由 $K_s$ 控制。
