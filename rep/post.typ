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

// --- main content ---

#align(center, text(
  20pt, font: "Harano Aji Gothic"
  )[
  無限次元ガウスの消去法を用いた\
  //#nk()型定理の改良
  #rp()の改良
])

#align(center)[
    （指導教員　関根 晃汰　准教授）\
    関根研究室　2131701 齋藤 悠希
]

#line(length: 100%)

#show: rest => columns(
  2,rest
)

= 背景と目的
非線形微分方程式は，さまざまな現象を数学モデル化することができる．この非線形微分方程式を解くことで，数学モデルの現象の解析を行うことができる．非線形微分方程式の解の導出には計算機が用いられるが，/*計算機容量の有限性のために，*/有限次元として問題を解くと，導出された解と実際の解には誤差が生じる．そのため，解の導出の精度を上げるために，無限次元で考えられる精度保証付き数値計算が必要となる．

非線形微分方程式の精度保証付き数値計算の手法の一つに，/*#nk()型定理*/#rp()がある．この定理は，非線形方程式を有限次元の問題として解を導出している．

本研究では，解の導出精度向上のために，#rp()無限次元ガウスの消去法を用いた改良手法を提案する．

//= #nk()型定理
= #rp()

$X,Y$をBanach空間，$cal(L) paren.l X,Y paren.r $を$X$から$Y$への有界線形作用素の集合とする．有界線形作用素$A^dagger in cal(L)(X,Y), A in cal(L)(Y,X)$を考え，作用素$F:X arrow.r Y$が$C^1$-#fc()微分可能とする．いま，$tilde(x) in X$に対して，正定数$Y_0, Z_0, Z_1$および非減少関数$Z_2(r)(r>0)$が存在して，次に不等式を満たすとする．
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

\

このとき，radii polynomialを以下で定義する．
$
p(r) := Z_2(r)r^2 - (1-Z_1-Z_0)r + Y_0
$

これに対し，$p(r_0)<0$となる$r_0>0$が存在するならば，$F(tilde(x))=0$を満たす解$tilde(x)$が$b in overline(B(tilde(x),r))$内に一意に存在する．

ここで，$D F (macron(x))$をFの $macron(x)$における#fc()微分，$A^dagger$を$D F (macron(x))$の近似，$A$を$A^dagger$の近似左逆作用素($A A^dagger approx I$)とする．


= 提案手法
#rp()で定義された作用素$A$は，有限次元上では作用素$A^dagger$の近似逆作用素である．無限次元上の場合，この$A^dagger$は真の逆作用素$A$となる．$A=D F^(-1)$より，式@y0 は以下になる．

$
||D F^(-1) F (tilde(x))||_X &lt.eq Y_0
$<tf0>


ここで，$phi.alt := D F^(-1) F (tilde(x))$とし，式@tf0 を以下のように変形して，ガウスの消去法を適用する．
$Pi_N$は射影作用素とする．

$
  D F phi.alt = F(tilde(x))\
$

$
  mat(
    Pi_N D F Pi_N, Pi_N D F (I-Pi_N) phi.alt;
    (I-Pi_N) D F Pi_N, (I-Pi_N) D F Pi_N;
  )&\
  mat(
    Pi_N phi.alt;
    (I-Pi_N) phi.alt;
  )
  =mat(
    Pi_N F(tilde(x));
    (I-Pi_N) F(tilde(x));
  )&
$<tf1>

= 今後の課題
提案手法で提示した式@tf1 のガウスの消去法による展開や，Juliaを用いたプログラムの実証を行う．

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
