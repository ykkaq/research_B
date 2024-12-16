using LinearAlgebra, DifferentialEquations, FFTW
include("FourierChebyshev.jl")

# van der Pol方程式
function vanderpol(du, u, μ, t)
  x, y = u
  du[1] = y
  du[2] = μ * (1 - x^2) * y - x
end

# F^(N)(x_n)
function F_fourier(x, μ, η₀)
  N = length(x) / 2
  ω = x[1]
  a = x[2:end]
  (a³, ~) = powerconvfourier(a, 3)
  eta = sum(a) - η₀

  k = -(N - 1):(N-1)
  f = (-k .^ 2 * ω^2 - μ * im * k * ω .+ 1) .* a + μ * im * k * ω .* a³ / 3

  return [eta; f]
end

function DF_fourier(x, μ)
  N = Int((length(x)) / 2)
  ω = x[1]
  a = x[2:end]
  k = (-N+1):(N-1)
  (a³, ~) = powerconvfourier(a, 3)

  DF = zeros(ComplexF64, 2N, 2N)

  DF[1, 2:end] .= 1
  DF[2:end, 1] = (-2 * ω * k .^ 2 - μ * im * k) .* a + μ * im * k .* a³ / 3

  (~, a2) = powerconvfourier(a, 2)

  M = zeros(ComplexF64, 2 * N - 1, 2 * N - 1)

  for j = (-N+1):(N-1)
    M[k.+N, j+N] = μ * im * k * ω .* a2[k.-j.+(2*N-1)]
  end

  L = diagm(-k .^ 2 * ω^2 - μ * im * k * ω .+ 1)

  DF[2:end, 2:end] = L + M
  return DF
end

# a function of  fourier coeffs
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

# たぶん，初期値設定
u_0 = [0.0; 2.0]
tspan = (0.0, 300)
μ = 1.0
prob = ODEProblem(vanderpol, u_0, tspan, μ)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
u = hcat(sol.u...)
ind = floor(Int, length(sol.t) / 2)

# おおよその周期
#a = 30
#b = 36.55
a = 30
app_period = 6.55
timestep = 0.1

f_tmp = sol(a+app_period/2:timestep:a+3*app_period/2)
find_period = abs.(f_tmp .- sol(a))
(~, ind) = findmin(find_period[1, :])
b = a + app_period / 2 + timestep * (ind - 1)
#calc fouriercoeffs

N = 50 # size of Fourier

println("size of Fourier = $N")
a_0 = odefouriercoeffs(sol, N, [a, b])

include("IntervalFunctions.jl")
# Initial value of Newton method
η_0 = 0.0
x = [2 * pi / (b - a); a_0]

# Newton iteration
tol = 5e-12
F = F_fourier(x, μ, η_0)
println("Before step #1, ||F||_1 = $(norm(F,1))")
num_itr = 0

while num_itr ≤ 100
  global x = x - DF_fourier(x, μ) \ F
  global num_itr += 1
  global F = F_fourier(x, μ, η_0)
  # println("After step #$(num_itr), ||F||_1 = $(norm(F,1))")
  if norm(F, 1) < tol
    break
  end
end

# A^(N)
ix = map(interval, x)
iω̄ = map(interval, real(x[1]))
iā = map(interval, x[2:end])

ix = x
iω̄ = real(x[1])
iā = x[2:end]

function DF_fourier(x::Vector{Complex{Interval{T}}}, μ) where {T}
  N = Int((length(x)) / 2)
  ω = x[1]
  a = x[2:end]
  k = (-N+1):(N-1)

  println("pcf Input: ", a, size(a), typeof(a))

  (a³, ~) = powerconvfourier(a, 3)
  DF = zeros(Complex{Interval{T}}, 2N, 2N)
  DF[1, 2:end] .= 1
  DF[2:end, 1] = (-2 * ω * k .^ 2 - μ * im * k) .* a + μ * im * k .* a³ / 3
  (~, a2) = powerconvfourier(a, 2)
  M = zeros(Complex{Interval{T}}, 2 * N - 1, 2 * N - 1)
  for j = (-N+1):(N-1)
    M[k.+N, j+N] = μ * im * k * ω .* a2[k.-j.+(2*N-1)]
  end
  L = diagm(-k .^ 2 * ω^2 - μ * im * k * ω .+ 1)
  DF[2:end, 2:end] = L + M
  return DF
end

iDF = DF_fourier(ix, μ);
iA = inv(iDF) # map(Interval,inv(mid.(iDF)))

## =======================
## I get a and omega by x.
## =======================
"""
println("x = $x")
println("μ = $μ")
println("N = $N")

μ
"""

omega = x[1]
a = x[2:end]
mu = μ

extend_size = 3 * N + 1
topleft_size = 2 * N

# define DF[]
zero_padding = zeros(ComplexF64, Int(extend_size))
extend_x = vcat(omega, zero_padding, a, zero_padding)
extend_DF = DF_fourier(extend_x, mu)

println(size(x), size(extend_x), size(extend_DF))

# define A
DF = DF_fourier(x, mu)
T_inv = inv(DF)
## lambda 行列
lambda_array = zeros(ComplexF64, size(extend_DF) .- size(T_inv))
for k in 1:size(lambda_array)[1]
  kk = k + N - 1
  lambda_array[k, k] = 1 / (-1 * kk * omega^2 - mu * im * kk * omega + 1)
end
topleft_size = size(T_inv)[1]

# norm_D
DF_bottomright = extend_DF[topleft_size+1:end, topleft_size+1:end]
D = lambda_array * DF_bottomright
D = 1.0I(size(D)[1]) - D
normD = maximum(sum(abs.(D), dims=1))

# norm_C
DF_bottomleft = extend_DF[topleft_size+1:end, 1:topleft_size]
C = lambda_array * DF_bottomleft
normC = maximum(sum(abs.(C), dims=1))

# norm_B
DF_topright = extend_DF[1:topleft_size, topleft_size+1:end]
B = T_inv * DF_topright
normB = maximum(sum(abs.(B), dims=1))
normDF_topright = maximum(sum(abs.(DF_topright), dims=1))
normT_inv = maximum(sum(abs.(T_inv), dims=1))

println("||I - D|| = ", normD)
println("||C|| = ", normC)
println("||B|| = ", normB)
println("||T^-1|| = ", normT_inv)
println("||M|| = ", normDF_topright)


norm_result = normD + normC * normB
norm_result2 = normD + normC * normT_inv * normDF_topright

println("||I - D|| + ||C|| ||B|| = ", norm_result)
println("||D|| + ||C|| ||T^-1|| ||M|| = ", norm_result2)

