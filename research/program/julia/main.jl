using DifferentialEquations
using LinearAlgebra
using FFTW


# フーリエ級数・チェビシェフ級数を扱う関数が記述されている（区間演算なし）。
include("FourierChebyshev.jl")

# 丸め向きを変更しない区間行列積（int_mul）や区間演算を使ったFFTのverifyfft、FFTを使った離散畳み込みの計算（powerconvfourier）が実装されている。入力のタイプによって区間行列積の場合分けをしている。
include("IntervalFunctions.jl")


# van der Pol方程式
function vanderpol(du, u, μ, t)
  x, y = u
  du[1] = y
  du[2] = μ * (1 - x^2) * y - x
end

# 初期値
u₀ = [0.0; 2.0]
tspan = (0.0, 300)
μ = 1.0
prob = ODEProblem(vanderpol, u₀, tspan, μ)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

#println("prob: $prob")
#println("sol: $sol")

#おおよその周期
a = 30
app_period = 6.55
timestep = 0.1

f_tmp = sol(a+app_period/2:timestep:a+3*app_period/2)
find_period = abs.(f_tmp .- sol(a))
(~, ind) = findmin(find_period[1, :])
b = a + app_period / 2 + timestep * (ind - 1)
# abs.(sol(b) .- sol(a))

#a function of  fourier coeffs
function odefouriercoeffs(f, N, I, n=1)
  a = I[1]
  b = I[2]
  # x_j: equidistance node points
  h = (b - a) / (2N - 1)
  j = 0:2N-2
  xⱼ = a .+ j * h
  # f_j: function values on node points
  fⱼ = f(xⱼ)[n, :]
  return (fftshift(fft(fⱼ))) / (2 * N - 1)
end

# compute Fourier fouriercoeffs
N = 61 # size of Fourier
a₀ = odefouriercoeffs(sol, N, [a, b])

println("$a₀")


