using LinearAlgebra, DifferentialEquations, FFTW

# van der Pol方程式
function vanderpol(du, u, μ, t)
  x, y = u
  du[1] = y
  du[2] = μ * (1 - x^2) * y - x
end

# 畳み込み
function powerconvfourier(a::Vector{Complex{T}}, p) where {T}
  M = Int((length(a) + 1) / 2)
  N = (p - 1) * M
  ta = [zeros(N, 1); a; zeros(N, 1)] # 1. Padding zeros: size(ta) = 2pM-1
  tb = ifft(ifftshift(ta)) # 2. IFFT of ta
  tbᵖ = tb .^ p # 3. tb*^tb
  cᵖ = fftshift(fft(tbᵖ)) * (2.0 * p * M - 1)^(p - 1)
  return cᵖ[N+1:end-N], cᵖ[p:end-(p-1)]# return (truncated, full) version
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

# DF^(N)(x_n)
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


## =======================
## I get a and omega by x.
## =======================
"""
println("x = $x")
println("μ = $μ")
println("N = $N")
"""

# norm_D
function norm_D(x, mu, N)
  ## get necessary value
  omega = x[1]
  a = x[2:end]
  finite_size = 3 * N

  ## Declare value
  (~, a_pcf2) = powerconvfourier(a, 2)
  D_finite = zeros(ComplexF64, finite_size, finite_size)
  ret = zeros(ComplexF64, length(a_pcf2) + 2)

  ## calc sup
  ret[1] = 1
  for k in N+3:N+3+length(a_pcf2)-1
    top = mu * im * k * omega * a_pcf2[k-(N+2)]
    lambda = -1 * k^2 * omega^2 - mu * im * k * omega + 1
    ret[k-N] = top / lambda
  end

  """
  ### diagonal condition
  for j in N+1:finite_size+N
    D_finite[j-N, j-N] = ComplexF64(1)
  end

  ### declare value
  for j in N+1:N+length(a_pcf2)+2 ##memo: fix range
    for k in N+1:N+length(a_pcf2)+2
      #memo: tate++ -> k++, yoko++ -> j++
      if 2 <= k - j && k - j <= length(a_pcf2) + 2   ## pcf condition
        top = mu * im * k * omega * a_pcf2[k-j+1]
        lambda = -1 * k^2 * omega^2 - mu * im * k * omega + 1
        D_finite[j-N+1, k-N+1] = top / lambda
      end
    end
  end
  """

  #println(size(ret))
  return norm(ret, 1)
end


# norm_C
function norm_C(x, mu, N)
  ## get necessary value
  omega = x[1]
  a = x[2:end]

  ## Declare value
  (~, a_pcf2) = powerconvfourier(a, 2)
  (~, a_pcf3) = powerconvfourier(a, 3)
  C_finite = zeros(ComplexF64, N, length(a_pcf3) + 3)  # 二次元配列として初期化


  ## calc sup
  ### differentiation omega
  for k in N+1:length(a_pcf3)+2
    lambda = -1 * k^2 * omega^2 - mu * im * k * omega + 1
    C_finite[1, k-N] = mu * im * k * omega * a_pcf3[k-2] / 3 / lambda
  end

  ### differentiation a_j
  for j in 2:N
    for k in N+1:N+length(a_pcf3)
      #memo: tate++ -> k++, yoko++ -> j++

      if 2 <= k - j && k - j < length(a_pcf2) + 2   ## pcf condition
        ## declare value
        top = mu * im * k * omega * a_pcf2[k-j-1]
        lambda = -1 * k^2 * omega^2 - mu * im * k * omega + 1

        C_finite[j, k-N] = top / lambda
      end
    end
  end

  #println(size(C_finite))
  return norm(C_finite, 1)
end

# norm_T_inv
function norm_T_inv(x, mu)
  T_inv = DF_fourier(x, mu)
  return norm(T_inv)
end

"""
# old_norm_M_01
function old_norm_M_01()

end
"""

# norm_B
function norm_B(x, mu)
  T_inv = DF_fourier(x, mu)
  return norm(T_inv[1, 1], 1)
end

println("norm_D = ", norm_D(x, μ, N))
println("norm_C = ", norm_C(x, μ, N))
println("norm_B = ", norm_B(x, μ))
#println("norm_T_inv = ", norm_T_inv(x, μ))
#println("old_norm_M_01 = ", old_norm_T_inv(x, μ,N))

