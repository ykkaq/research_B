#van der Pol方程式
include("FourierChebyshev.jl")

function vanderpol(du, u, μ, t)
  x, y = u
  du[1] = y
  du[2] = μ * (1 - x^2) * y - x
end

function F_fourier(x, μ, η_0)
  N = length(x) / 2
  ω = x[1]
  a = x[2:end]
  (a³, ~) = powerconvfourier(a, 3)
  eta = sum(a) - η_0

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


using DifferentialEquations
u_0 = [0.0; 2.0]
tspan = (0.0, 300)
μ = 1.0
prob = ODEProblem(vanderpol, u_0, tspan, μ)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
u = hcat(sol.u...)
ind = floor(Int, length(sol.t) / 2)


#おおよその周期
# a = 30
# b = 36.55
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


using LinearAlgebra
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
include("IntervalFunctions.jl")
setformat(:full)

ix = map(Interval, x)
iω̄ = map(Interval, real(x[1]))
iā = map(Interval, x[2:end])
ν = 1.05

function DF_fourier(x::Vector{Complex{Interval{T}}}, μ) where {T}
  N = Int((length(x)) / 2)
  ω = x[1]
  a = x[2:end]
  k = (-N+1):(N-1)
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
iA = map(Interval, inv(mid.(iDF)))
A_a0 = iA[1, 2:end]
A_a1 = iA[2:end, 2:end]
A_01 = iA[2:end, 1];

function F_fourier_ext(x::Vector{Complex{Interval{T}}}, μ, η_0) where {T}
  N = length(x) / 2
  ω = x[1]
  a = [zeros(Complex{Interval{T}}, 2 * (Int(N) - 1)); x[2:end]; zeros(Complex{Interval{T}}, 2 * (Int(N) - 1))]
  (~, a³) = powerconvfourier(x[2:end], 3)
  eta = sum(a) - η_0

  k = -3*(N-1):3*(N-1)
  f = (-k .^ 2 * ω^2 - μ * im * k * ω .+ 1) .* a + μ * im * k * ω .* a³ / 3

  return [eta; f]
end

function wnorm(a, ν)
  N = (length(a) + 1) / 2 # length(a) = 2*N-1
  k = (-N+1):(N-1)
  w = ν .^ abs.(k)
  return sum(abs.(a) .* w)
end

δ = F_fourier_ext(ix, μ, η_0)
δ_0 = δ[1]
δ_1 = δ[2:end]
δ_1_N = δ_1[2*(N-1)+1:end-2*(N-1)] #N-1 ,1 , N-1 = 2N-1
δ_1[2*(N-1)+1:end-2*(N-1)] .= 0
δ_1_tail = δ_1

λ_k(k, ω) = -k .^ 2 * ω^2 - μ * im * k * ω .+ 1

k_tail = -3*(N-1):3*(N-1)
Y0 = sup(max(abs(iA[1, 1] * δ_0 + dot(A_a0, δ_1_N)), wnorm(A_01 * δ_0 + A_a1 * δ_1_N, ν) + wnorm(δ_1_tail ./ (abs.(λ_k(map(Interval, Vector(k_tail)), iω̄))), ν)))

@show Y0;


# Z0 bounds
function wnorm_mat(A, ν)
  m = size(A, 1) # m = 2*N-1
  N = (m + 1) / 2
  k = -N+1:N-1
  w = ν .^ abs.(k)
  return maximum(sum(w .* abs.(A), dims=1) ./ w')
end

function wsnorm(a, ν) # the input should be vector
  # the norm of dual space of the weighted ell^1
  m = length(a) # m = 2*N-1
  N = (m + 1) / 2
  k = -N+1:N-1
  w = ν .^ abs.(k)
  return maximum(abs.(a) ./ w)
end

B = I - iA * iDF #2N × 2N
Z0_0 = abs(B[1, 1]) + wsnorm(B[1, 2:end], ν)
Z0_1 = wnorm(B[2:end, 1], ν) + wnorm_mat(B[2:end, 2:end], ν)
Z0 = sup(max(Z0_0, Z0_1))
println("Z0 = $Z0")


# Z1 bounds
(~, ia²) = powerconvfourier(iā, 2)
(~, ia³) = powerconvfourier(iā, 3)

ζ = map(Interval, zeros(2 * N - 1))
for ell = -N+1:N-1
  j = ell-2*(N-1):-N
  if isempty(j)
    ζ_1 = -1
  else
    w_j = ν .^ abs.(j)
    ζ_1 = abs(μ * im * ell * iω̄) * maximum(abs.(ia²[ell.-j.+2*N.-1]) ./ w_j)
  end
  j = N:ell+2*(N-1)
  if isempty(j)
    ζ_2 = -1
  else
    w_j = ν .^ abs.(j)
    ζ_2 = abs(μ * im * ell * iω̄) * maximum(abs.(ia²[ell.-j.+2*N.-1]) ./ w_j)
  end
  ζ[ell+N] = max(ζ_1, ζ_2)
end

conv = map(Interval, 0)
for k = N:2*(N-1)
  #positive
  global conv += abs(μ * im * k * ia³[k+4(N-1)+1]) * ν^(k) / (3 * abs(λ_k(k, iω̄)))
  #negative
  global conv += abs(-μ * im * k * ia³[-k+2*(N-1)+1] * ν^(k)) / (3 * abs(λ_k(-k, iω̄)))
end

w_n = ν^(N)
iā_norm = wnorm(iā, ν)
Z1_0 = abs(iA[1, 1]) / w_n + dot(abs.(A_a0), ζ)
Z1_1 = wnorm(A_01, ν) / w_n + wnorm(abs.(A_a1) * ζ, ν) + conv + abs(μ * im * iω̄) * iā_norm^2 / (N * iω̄^2 - 1 / N)
Z1 = sup(max(Z1_0, Z1_1))
println("Z1 = $Z1")


#Z2 bound
function bopnorm(A, tail_es, ν) # the operator norm of bounded operators with tail
  return max(wnorm_mat(A, ν), tail_es)
end

k = -N+1:N-1
Ã = abs.(k) .* abs.(A_a0)
B̃ = (k .^ 2) .* abs.(A_a0)
Ã_norm = wsnorm(Ã, ν)
B̃_norm = wsnorm(B̃, ν)
A_norm = wsnorm(A_a0, ν)

Z2_20 = B̃_norm * (4 * iω̄ + 2 * iā_norm) + 2 * μ * Ã_norm * (1 + iω̄ * iā_norm + iā_norm^2)
Z2_30 = 3 * B̃_norm + μ * Ã_norm * (3 * iā_norm + iω̄)
Z2_40 = 4 * μ * Ã_norm / 3

tA = transpose(abs.(k)) .* abs.(A_a1)
tB = transpose(k .^ 2) .* abs.(A_a1)
tA_bopnorm = bopnorm(tA, (N + 1) / abs(λ_k(N + 1, iω̄)), ν)
tB_bopnorm = bopnorm(tB, 1 / (iω̄^2 - 1 / (N^2)), ν)

Z2_21 = tB_bopnorm * (4 * iω̄ + 2 * iā_norm) + 2 * μ * tA_bopnorm * (1 + iω̄ * iā_norm + iā_norm^2)
Z2_31 = 3 * tB_bopnorm + μ * tA_bopnorm * (3 * iā_norm + iω̄)
Z2_41 = 4 * μ * tA_bopnorm / 3

Z2_2 = sup(max(Z2_20, Z2_21))
Z2_3 = sup(max(Z2_30, Z2_31))
Z2_4 = sup(max(Z2_40, Z2_41))

@show Z2_2
@show Z2_3
@show Z2_4

Z2(r) = Z2_4 * r .^ 3 + Z2_3 * r .^ 2 + Z2_2 * r


using ForwardDiff
#ニュートン法で近似解を計算する
function newton(F, x0)
  tol = 5e-10
  count = 0
  x = x0
  Fx = F(x)
  while maximum(abs, Fx) ≥ tol && count ≤ 20
    DF = ForwardDiff.derivative(F, x)
    x -= DF \ Fx
    Fx = F(x)
    count += 1
  end
  return x
end

#クラフチック写像
function krawczyk(F, X)
  iDF = ForwardDiff.derivative(F, X)
  c = mid.(X)
  ic = map(Interval, c)
  DF = ForwardDiff.derivative(F, c)
  R = inv(DF)
  M = I - R * iDF
  return c - R * F(ic) + M * (X - c)
end

function verifynlss_krawczyk(F, c)
  DF = ForwardDiff.derivative(F, c)
  R = inv(DF)
  r = abs.(R * F(c))
  u = r .+ (sum(r) / length(r))
  X = c .± u
  K = krawczyk(F, X)
  if all(K .⊂ X)
    tol = 5e-10
    count = 0
    while maximum(radius, K) >= tol && count ≤ 100
      K = krawczyk(F, K)
      count += 1
    end
    success = 1
    return success, K
  end
  println("Oh my way, verification is failed...return a improved approximatesolution")
  success = 0
  return success, newton(F, c)
end


rp(r) = Z2(r) .* r - (1 - Z1 - Z0) * r + Y0 # radii-polynomial
r0_mid = newton(rp, 1e-10)
success, r0 = verifynlss_krawczyk(rp, r0_mid)
@show success
@show r0
rp(interval(sup(r0))) < 0