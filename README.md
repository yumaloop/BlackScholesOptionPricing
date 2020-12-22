# BlackScholesOptionPricing

This is the `R` implementation with [YUIMA](https://yuimaproject.com/) library
for A comparison between the Monte Carlo siumlation and the Black-Scholes option pricing model.

### Option Pricing

**ブラックショールズ式によるオプション価格の導出**

満期$T$で原資産価格(株式価格)が連続時間確率過程$S = {(S_t)}_{t \in [0,
T]}$に従うオプションの時刻$t \in [0, T]$における価格$C(t, S_t)$について考える．$S$が確率微分方程式
\begin{align}
    d S_t = \sigma S_t dt + \mu S_t d W_t
\end{align}
の解で与えられているとする($S$は幾何ブラウン運動に従う)．時刻$t \in [0, T]$において，オプションの原資産価格$S_t$とオプションの行使価格$K$が与えられたとき，ブラック・ショールズ式によってオプションの理論価格$C(t, S_t)$は

<a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{150}&space;\fn_cm&space;\begin{align}&space;C(t,&space;S_t)&space;&=&space;S_t&space;\Phi(d_1)&space;-&space;K&space;e^{-r(T-t)}&space;\Phi(d_2)&space;\nonumber&space;\\&space;where&space;~~&space;d_1&space;&=&space;\frac{\log&space;\left(&space;\frac{S_t}{K}&space;\right)&space;&plus;&space;\left(&space;r&space;&plus;&space;\frac{\sigma^2}{2}&space;\right)&space;T}{\sigma&space;\sqrt{T}}&space;\nonumber&space;\\&space;d_2&space;&=&space;\frac{\log&space;\left(&space;\frac{S_t}{K}&space;\right)&space;&plus;&space;\left(&space;r&space;-&space;\frac{\sigma^2}{2}&space;\right)&space;T}{\sigma&space;\sqrt{T}}&space;\nonumber&space;\end{align}" target="_blank"><img src="https://latex.codecogs.com/png.latex?\dpi{150}&space;\fn_cm&space;\begin{align}&space;C(t,&space;S_t)&space;&=&space;S_t&space;\Phi(d_1)&space;-&space;K&space;e^{-r(T-t)}&space;\Phi(d_2)&space;\nonumber&space;\\&space;where&space;~~&space;d_1&space;&=&space;\frac{\log&space;\left(&space;\frac{S_t}{K}&space;\right)&space;&plus;&space;\left(&space;r&space;&plus;&space;\frac{\sigma^2}{2}&space;\right)&space;T}{\sigma&space;\sqrt{T}}&space;\nonumber&space;\\&space;d_2&space;&=&space;\frac{\log&space;\left(&space;\frac{S_t}{K}&space;\right)&space;&plus;&space;\left(&space;r&space;-&space;\frac{\sigma^2}{2}&space;\right)&space;T}{\sigma&space;\sqrt{T}}&space;\nonumber&space;\end{align}" title="\begin{align} C(t, S_t) &= S_t \Phi(d_1) - K e^{-r(T-t)} \Phi(d_2) \nonumber \\ where ~~ d_1 &= \frac{\log \left( \frac{S_t}{K} \right) + \left( r + \frac{\sigma^2}{2} \right) T}{\sigma \sqrt{T}} \nonumber \\ d_2 &= \frac{\log \left( \frac{S_t}{K} \right) + \left( r - \frac{\sigma^2}{2} \right) T}{\sigma \sqrt{T}} \nonumber \end{align}" /></a>

となる．ただし，$r$は無リスク資産の利率とする．

ここで，簡単のために満期$T$を$1$とすると，現在($t=0$)のオプション価格の理論値は，

<a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{150}&space;\fn_cm&space;\begin{align}&space;C(0,&space;S_0)&space;&=&space;S_0&space;\Phi(d_1)&space;-&space;K&space;e^{-r}&space;\Phi(d_2)&space;\nonumber&space;\\&space;where&space;~~&space;d_1&space;&=&space;\frac{\log&space;\left(&space;\frac{S_0}{K}&space;\right)&space;&plus;&space;\left(&space;r&space;&plus;&space;\frac{\sigma^2}{2}&space;\right)}{\sigma&space;}&space;\nonumber&space;\\&space;d_2&space;&=&space;\frac{\log&space;\left(&space;\frac{S_0}{K}&space;\right)&space;&plus;&space;\left(&space;r&space;-&space;\frac{\sigma^2}{2}&space;\right)}{\sigma}&space;\nonumber&space;\end{align}" target="_blank"><img src="https://latex.codecogs.com/png.latex?\dpi{150}&space;\fn_cm&space;\begin{align}&space;C(0,&space;S_0)&space;&=&space;S_0&space;\Phi(d_1)&space;-&space;K&space;e^{-r}&space;\Phi(d_2)&space;\nonumber&space;\\&space;where&space;~~&space;d_1&space;&=&space;\frac{\log&space;\left(&space;\frac{S_0}{K}&space;\right)&space;&plus;&space;\left(&space;r&space;&plus;&space;\frac{\sigma^2}{2}&space;\right)}{\sigma&space;}&space;\nonumber&space;\\&space;d_2&space;&=&space;\frac{\log&space;\left(&space;\frac{S_0}{K}&space;\right)&space;&plus;&space;\left(&space;r&space;-&space;\frac{\sigma^2}{2}&space;\right)}{\sigma}&space;\nonumber&space;\end{align}" title="\begin{align} C(0, S_0) &= S_0 \Phi(d_1) - K e^{-r} \Phi(d_2) \nonumber \\ where ~~ d_1 &= \frac{\log \left( \frac{S_0}{K} \right) + \left( r + \frac{\sigma^2}{2} \right)}{\sigma } \nonumber \\ d_2 &= \frac{\log \left( \frac{S_0}{K} \right) + \left( r - \frac{\sigma^2}{2} \right)}{\sigma} \nonumber \end{align}" /></a>

となる．

**モンテカルロ法によるオプション価格の導出**

ブラックショールズ式に含まれる$S_T$の期待値計算をモンテカルロ法によって近似することを考える．すなわち，幾何ブラウン運動に従うサンプルパス$S = {(S_t)}_{t \in [0,
T]}$を大量に生成することで，$S_T$の期待値を求める．

モンテカルロシミュレーションによって生成された，満期$T$における原資産価格$S_T$の$n$個のサンプルを$(s^{(1)}_{T}, \cdots, s^{(n)}_{T})$とすると，
時刻$t \in [0, T]$におけるオプションの価格の推定値$\hat{C}(t, S_t)$は

<a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{150}&space;\fn_cm&space;\begin{align}&space;\hat{C}(t,&space;S_t)&space;&=&space;\frac{1}{n}&space;\sum_{i=1}^{n}&space;e^{-r(T-t)}&space;\cdot&space;max(s^{(i)}_{T}&space;-&space;K,&space;0)&space;\nonumber&space;\end{align}" target="_blank"><img src="https://latex.codecogs.com/png.latex?\dpi{150}&space;\fn_cm&space;\begin{align}&space;\hat{C}(t,&space;S_t)&space;&=&space;\frac{1}{n}&space;\sum_{i=1}^{n}&space;e^{-r(T-t)}&space;\cdot&space;max(s^{(i)}_{T}&space;-&space;K,&space;0)&space;\nonumber&space;\end{align}" title="\begin{align} \hat{C}(t, S_t) &= \frac{1}{n} \sum_{i=1}^{n} e^{-r(T-t)} \cdot max(s^{(i)}_{T} - K, 0) \nonumber \end{align}" /></a>

と求められる．
ここで，簡単のために満期$T$を$1$とすると，現在($t=0$)のオプション価格の推定値は，

<a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{150}&space;\fn_cm&space;\begin{align}&space;\hat{C}(0,&space;S_0)&space;&=&space;\frac{1}{n}&space;\sum_{i=1}^{n}&space;e^{-r}&space;\cdot&space;max(s^{(i)}_{1}&space;-&space;K,&space;0)&space;\nonumber&space;\end{align}" target="_blank"><img src="https://latex.codecogs.com/png.latex?\dpi{150}&space;\fn_cm&space;\begin{align}&space;\hat{C}(0,&space;S_0)&space;&=&space;\frac{1}{n}&space;\sum_{i=1}^{n}&space;e^{-r}&space;\cdot&space;max(s^{(i)}_{1}&space;-&space;K,&space;0)&space;\nonumber&space;\end{align}" title="\begin{align} \hat{C}(0, S_0) &= \frac{1}{n} \sum_{i=1}^{n} e^{-r} \cdot max(s^{(i)}_{1} - K, 0) \nonumber \end{align}" /></a>

となる．

<img src="https://github.com/yumaloop/BlackScholesOptionPricing/blob/main/figures/BSCallOptionMC1.png">

<img src="https://github.com/yumaloop/BlackScholesOptionPricing/blob/main/figures/BSCallOptionMC2.png">
