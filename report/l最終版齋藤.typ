// --- settings ---

#set page(
  paper: "a4",
  margin: (
    x:15mm,
    y:10mm
  )
)

#set par(
  first-line-indent: 1em,
  //linebreaks: "optimized",
  justify: false,
  leading: 0.75em
)

#show par: set block(
  spacing: 0.65em
)

// font
#set text(
  lang:"ja",
  font: "Harano Aji Mincho",
  size: 12pt
)

// heading
#set heading(
  //  headingに1.をつける
  numbering: "1.",
)
#show heading: set text(
  //  headingのフォントを変更
  font: "Harano Aji Gothic"
  // size: 15pt,
)
#show heading: it => {
  //  最初の行をインデントする．
  it
  par(text(size: 0pt, ""))
}

// math numbering
#set math.equation(
  numbering: "(1)",
  number-align: right
)

// shortcut
#let fc() = "Fr"+str.from-unicode(233)+"chet"
#let nk() = "Newton-Kantorovich"
#let rp() = "radii-polynomial approach"

#show ref: it => {
  let eq = math.equation
  let el = it.element
  if el != none and el.func() == eq {
    // Override equation references.
    numbering(
      el.numbering,
      ..counter(eq).at(el.location())
    )
  } else {
    // Other references as usual.
    it
  }
}

//---------------------//1d
// --- main content ---

#align(center, text(
  20pt, font: "Harano Aji Gothic"
  )[
  無限次元ガウスの消去法を用いた\
  //#nk()型定理の改良
  #rp()の改良
])

#align(center)[
\
関根研究室　2131701　齋藤 悠希
]

#line(length: 100%)

#show: rest => columns(
  2,rest
)

= 背景と目的
非線形微分方程式は，さまざまな現象を数学モデル化することができる．この非線形微分方程式を解くことで，数学モデルの現象の解析を行うことができる．

#rp()は，非線形方程式の解の精度保証付き数値計算に関する定理である．従来の#rp()は，難しい計算部分を簡略化している．そのため，精度保証の性能が低くなる．そこで，この簡略化をせずに、無限次元ガウスの消去法を用いて，精度保証の性能を向上させる．

= #rp()

- $X,Y$はBanach空間
- $cal(L) paren.l X,Y paren.r $は$X$から$Y$への有界線形作用素集合
- 有界線形作用素$A^dagger in cal(L)(X,Y), A in cal(L)(Y,X)$
- 作用素$F:X arrow.r Y$が$C^1$-#fc()微分可能

$tilde(x) in X$に対して，正定数$Y_0, Z_0, Z_1$および非減少関数$Z_2(r)(r>0)$が存在して，次の不等式を満たすとする．
$
||A F (tilde(x))||_X &lt.eq Y_0
$<y0>
$
||I-A A^dagger||_(cal(L)(X)) &lt.eq Z_0 \
$
$
||A (D F(tilde(x))-A^dagger)||_(cal(L)(X)) &lt.eq Z_1 \
$
$
||A (D F(b)-D F (tilde(x)))||_(cal(L)(X)) lt.eq Z_2(r)& \
forall b in overline(B(tilde(x),r))&
$

このとき，radii polynomialを以下で定義する．
$
p(r) := Z_2(r)r^2 - (1-Z_1-Z_0)r + Y_0
$

これに対し，$p(r_0)<0$となる$r_0>0$が存在するならば，$F(x)=0$を満たす解$tilde(x)$が$overline(B(x,r))$内に一意に存在する．

ここで，$D F (macron(x))$を$F$の $macron(x)$における#fc()微分，$A^dagger$を$D F (macron(x))$の近似，$A$を$A^dagger$の近似左逆作用素($A A^dagger approx I$)とする．


= 提案手法
従来の#rp()では，$A^dagger,A$は，簡略化のために近似した作用素としていた．

提案手法では，$A^dagger,A$を以下のように，正確なニュートン法である作用素として計算をする．
$
A^dagger = D F(macron(x)), #h(10pt) A = D F(macron(x))^(-1)
$<cor0>

この式@cor0 を，#rp()の式@y0 に代入する．

$
||D F(macron(x))^(-1) F (tilde(x))||_X &lt.eq Y_0
$<tf0>


ここで，$phi.alt := D F (macron(x))^(-1) F (tilde(x))$とする．この式に，無限次元ガウス消去法を用いて変形する．ただし、$Pi_N$は射影作用素とする．

$
  D F (macron(x)) phi.alt = F(tilde(x))\
$

#set text(size: 9.5pt)
$
  &mat(
    Pi_N D F Pi_N, Pi_N D F (I-Pi_N) phi.alt;
    (I-Pi_N) D F Pi_N, (I-Pi_N) D F Pi_N;
  )
  mat(
    Pi_N phi.alt;
    (I-Pi_N) phi.alt;
  )\
  &=mat(
    Pi_N F(tilde(x));
    (I-Pi_N) F(tilde(x));
  )
$<tf1>


#set text(size: 12pt)
= 今後の課題
- 提案手法で提示した式@tf1 のガウスの消去法による展開
- Juliaを用いたプログラムの実証

// 参考文献
#set heading(numbering: none)
#set enum(numbering: "[1]")
= 参考文献
+ 高安亮紀，Julia言語を使った精度保証付き数値計算のチュートリアル


//初めにの一貫性　全体の話を書く
//具体性？かな

//式６がおも
//式1を式６に <- 1.に書く
//引用はあと
//式８は行列
//
