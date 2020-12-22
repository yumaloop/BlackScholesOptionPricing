# BlackScholesOptionPricing

This is the `R` implementation with [YUIMA](https://yuimaproject.com/) library
for A comparison between the Monte Carlo siumlation and the Black-Scholes option pricing model.

See [qitta post](https://qiita.com/yumaloop/items/58e134f03a452eaed335).


<img src="https://github.com/yumaloop/BlackScholesOptionPricing/blob/main/figures/BSCallOptionMC1.png">

- パラメータ設定
    - n: 1000 # モンテカルロシミュレーションの回数
    - K: 900 # コールオプションの権利行使価格 (t=T)
    - S: 1000 # 株式の現在価格 (t=0)
    - r: 0.005 # 幾何ブラウン運動のdrift
    - sigma: 0.3 # 幾何ブラウン運動のdiffusion
    - T: 1 # オプションの満期

<img src="https://github.com/yumaloop/BlackScholesOptionPricing/blob/main/figures/BSCallOptionMC2.png">

- パラメータ設定
    - n: 1000 # モンテカルロシミュレーションの回数
    - K: 1100 # コールオプションの権利行使価格 (t=T)
    - S: 1000 # 株式の現在価格 (t=0)
    - r: 0.005 # 幾何ブラウン運動のdrift
    - sigma: 0.3 # 幾何ブラウン運動のdiffusion
    - T: 1 # オプションの満期
