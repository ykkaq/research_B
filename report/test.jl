F(x) = x^2 - 2
Df(x) = 2 * x

num_itr = 0
tol = 1e-12
x_0 = 1.0
x = x_0
Fx = F(x)
DF = Df(x)
println("Before iteration: $(Fx)")


#ニュートン法を計算する
while num_itr ≤ 10
  dx = Fx / Df(x)
  x = x - dx
  num_itr += 1
  Fx = F(x)
  println("After $(num_itr) th iteration: $(Fx)")
  if abs(Fx) < tol
    break
  end
end

println(x)