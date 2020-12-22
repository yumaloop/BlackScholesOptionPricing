BlackScholesCallPrice = function(S, K, r, sigma, T=1)
{
    d1 <- ( log(S/K) + (r + sigma^2/2) * T)/( sigma * sqrt(T))
    d2 <- ( log(S/K) + (r - sigma^2/2) * T)/( sigma * sqrt(T))
    C0 <- S * pnorm(d1) - K * exp(-r * T) * pnorm(d2)
    return(C0)
}

MonteCarloCallPrice = function(S, K, r, sigma, n, T=1)
{
    n_sample <- 1000
    c0_list <- list()
    c0 <- 0
    for (i in 1:n) {
        resultGBM <- GBM_sample(S, r, sigma, n_sample)
        sT <- resultGBM@data@original.data[n_sample]
        c0 <- (1/i) * exp(-1*r*T) * max(sT - K, 0) + ((i-1)/i) * c0
        c0_list <- append(c0_list, list(c0))
    }
    return(c0_list)
}

GBM_sample = function(x0, alpha, beta, n_sample, T=1)
{
    # Step1:確率微分方程式モデルの定義
    # dS_t = alpha * S_t * dt + beta * S_t * dW_t
    mod <-setModel(drift="alpha*x", diffusion="beta*x")
    
    # Step2:サンプル時点の定義
    samp <-setSampling(Initial=0, Terminal=T, n=n_sample)
    
    # Step3:統計モデルの定義
    smod <-setYuima(model=mod, sampling=samp)
    
    # Step4:シミュレーションの実行
    # set.seed(9999)
    xinit <- x0
    param <- list(alpha=alpha, beta=beta) #パラメータの値
    resultGBM <- simulate(smod, xinit=xinit, true.parameter=param)
    return(resultGBM)
}

# Call-Option Pricing
n <- 1000 # Num. of MonteCalro simulation
K <- 900 # Option Exercise Price (t=T)
S <- 1000 # Curent Stock Price (t=0)
r <- 0.005 # Drift for Geometric Brownian motion
sigma <- 0.3 # Diffusion for Geometric Brownian motion
T <- 1 # Optional Term

bs_price <- BlackScholesCallPrice(S, K, r, sigma, T=1)
mc_price <- MonteCarloCallPrice(S, K, r, sigma, n, T=1)
print(bs_price)
plot(1:n, mc_price, main="Monte Carlo Simulation:\nBlack-Scholes Option Pricing Model", xlab="Number of sample paths: # of ST", ylab="Option Price: C0", cex=0.5)
abline(h=bs_price, col='red', lwd=1, lty=2)
