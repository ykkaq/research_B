using DifferentialEquations
using Plots
using FFTW
using LinearAlgebra


function vanderpol(du, u, μ, t)
  x, y = u
  du[1] = y
  du[2] = μ * (1 - x^2) * y - x
end

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

function plot_solution(u, index) # u = [ω, a_{-N+1}, ..., a_0, ..., a_{N-1}], length(u) = 2N
  # index = 1: profile of solution
  #         2: Fourier mode
  #         3: phase profile
  ω = real(u[1])
  L = 2π / ω
  a = u[2:end]
  N = length(u) / 2 # N: size of Fourier
  n_pad = 1000
  a_pad = [zeros(n_pad); a; zeros(n_pad)]
  N_pad = N + n_pad
  dx = L / (2 * N_pad - 1)
  x = dx * (0:2*N_pad-2)
  if index == 1
    # Plot profile:
    plot(x, real((2N_pad - 1) * ifft(ifftshift(a_pad))),
      xlabel="\$t\$",
      ylabel="\$x\\,(t)\$",
      line=1.6,
      title="Profile of the solution",
      size=(800, 400),
      legend=false,
    )
  elseif index == 2
    # Plot Fourier coefficients:
    plot((-N+1):(N-1), abs.(a), yscale=:log10,
      xlabel="\$k\$",
      ylabel="\$|a_k\\,|\$",
      line=1.6,
      title="Fourier coefficients of the solution",
      size=(800, 400),
      legend=false,
    )
  elseif index == 3
    # Plot phase:
    k = (-N_pad+1):(N_pad-1)
    plot(real((2N_pad - 1) * ifft(ifftshift(a_pad))), real((2N_pad - 1) * ifft(ifftshift(a_pad .* (im * k * ω)))),
      xlabel="\$x(t)\$",
      ylabel="\$\\dot{x}\\,(t)\$",
      line=1.6,
      title="Phase plot of a numerical solution",
      size=(800, 400),
      legend=false,
    )
  end
end

function plot_solution!(u)
  L = 2π / real(u[1])
  a = u[2:end]
  N = length(u) / 2
  n_pad = 1000
  a_pad = [zeros(n_pad); a; zeros(n_pad)]
  N_pad = N + n_pad
  k = (-N_pad+1):(N_pad-1)
  dx = L / (2 * N_pad - 1)
  x = dx * (0:2*N_pad-2)
  plot!(real((2N_pad - 1) * ifft(ifftshift(a_pad))), real((2N_pad - 1) * ifft(ifftshift(a_pad .* (im * k)))),
    line=1.6,
    label="Fourier interpolation"
  )
end

function powerconvfourier(a::Vector{Complex{T}}, p) where {T}
  M = Int((length(a) + 1) / 2)
  N = (p - 1) * M
  ta = [zeros(N, 1); a; zeros(N, 1)] # 1. Padding zeros: size(ta) = 2pM-1
  tb = ifft(ifftshift(ta)) # 2. IFFT of ta
  tbᵖ = tb .^ p # 3. tb*^tb
  cᵖ = fftshift(fft(tbᵖ)) * (2.0 * p * M - 1)^(p - 1)
  return cᵖ[N+1:end-N], cᵖ[p:end-(p-1)]# return (truncated, full) version
end

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