// Get Polylux from the official package repository
#import "@preview/polylux:0.3.1": *
//https://typst.app/universe/package/polylux/
#import "@preview/cetz:0.2.2"
//https://github.com/cetz-package/cetz
#import table: cell, header

// Make the paper dimensions fit for a presentation and the text larger
#set page(
  paper: "presentation-16-9"
)
#import themes.simple: *

#set text(
  lang:"ja",
  font: "Harano Aji Mincho",
  size: 26pt,
)

#show heading: set text(
  //  headingのフォントを変更
  font: "Harano Aji Gothic"
)

// shortcut
#let fc() = "Fr"+str.from-unicode(233)+"chet"
#let nk() = "Newton-Kantorovich"
#let rp() = "radii-polynomial approach"
#set underline(
  offset: 3pt
)

// --- maintain ---

// タイトルスライド
#slide[
  #set align(horizon + center)
  #set text(size:32pt)
  == 無限次元ガウスの消去法を用いた\ #rp()の改良
  \
  #set text(size:25pt)
  関根研究室　2131701　齋藤 悠希
]

// 2枚目
#slide[
  ==  背景と目的
  //#set align(left + horizon)
  #set align(center + horizon)

  #rp()・・・非線形方程式の解の精度保証に使われる定理

  // table
  #set align(center + horizon)
  既存手法と提案手法の比較#table(
    columns: 3,
    align: center,
    inset:10pt,
    header(
      [],
      [*計算*],
      [*精度*],
    ),
    [既存手法],
    [簡略化],
    cell(
      align: center,
      fill: blue.lighten(80%),
      [悪い],//[#sym.arrow.b],
    ),
    [提案手法],
    [無限次元ガウスの消去法],
    cell(
      align: center,
      fill: red.lighten(80%),
      [良い],//[#sym.arrow.t],
    ),
  )

  \
  簡略化部分を、無限次元ガウスの消去法で計算する\ #sym.arrow.r 精度を改善する


/*
  - 既存の#rp()では，計算が簡略化されている \ #sym.arrow.r 精度に難がある

  /*
  #set align(center + horizon)
  #sym.arrow.b
  */

  #set align(left + horizon)
  - 無限次元ガウスの消去法を用いて，#rp()の\ 精度を改善する．
*/

]


// 3枚目
#slide[
  == #rp()
  #set align(horizon)
  $tilde(x) in X$に対して、正定数$Y_0,Z_0,Z_1$および、非減少関数$Z_2(r)(r>0)$が存在して、次の式を満たすとする．

  $
  ||A F (tilde(x))||_X &lt.eq Y_0 \
  ||I-A A^dagger||_(cal(L)(X)) &lt.eq Z_0 \
  ||A (D F(tilde(x))-A^dagger)||_(cal(L)(X)) &lt.eq Z_1 \
  ||A (D F(b)-D F (tilde(x)))||_(cal(L)(X)) &lt.eq Z_2(r), #h(10pt) forall b &in overline(B(tilde(x),r))
  $
]

#slide[
  == #rp()
  #set align(horizon)
  このとき，radii polynomialを以下で定義する．
  $
  p(r) := Z_2(r)r^2 - (1-Z_1-Z_0)r + Y_0
  $

  これに対し、$p(r_0)<0$となる$r_0>0$が存在するならば、$F(tilde(x))=0$を満たす解 $tilde(x)$ が $b in overline(B(tilde(x),r))$ 内に一意に存在する．
]


// --- tmp --- //
/*
// 3枚目
#slide[
  == #rp()
  #set align(horizon)

  - $X,Y$を Banach空間

  - $cal(L) paren.l X,Y paren.r $を $X$から$Y$への有界線形作用素の集合

  - 有界線形作用素$A^dagger in cal(L)(X,Y), A in cal(L)(Y,X)$

  - 作用素$F:X arrow.r Y$が$C^1$-#fc()微分可能とする．
]

#slide[
  == #rp()
  #set align(horizon)

  $tilde(x) in X$に対して，正定数$Y_0, Z_0, Z_1$および非減少関数$Z_2(r)(r>0)$が存在して，次に不等式を満たすとする．

  $
  ||A F (tilde(x))||_X &lt.eq Y_0 \
  ||I-A A^dagger||_(cal(L)(X)) &lt.eq Z_0 \
  ||A (D F(tilde(x))-A^dagger)||_(cal(L)(X)) &lt.eq Z_1 \
  ||A (D F(b)-D F (tilde(x)))||_(cal(L)(X)) &lt.eq Z_2(r) \
  forall b &in overline(B(tilde(x),r))
  $

]

#slide[
  == #rp()
  #set align(horizon)

  このとき，radii polynomialを以下で定義する．
  $
  p(r) := Z_2(r)r^2 - (1-Z_1-Z_0)r + Y_0
  $

  これに対し，$p(r_0)<0$となる$r_0>0$が存在するならば，$F(tilde(x))=0$を満たす解$tilde(x)$が$b in overline(B(tilde(x),r))$内に一意に存在する．
]
*/
// --- tmp-end --- //

#slide[
  == 提案手法

  #set align(horizon+left)
  #rp()の一部， $Y_0$の評価式 \

  #set align(center)
  $||A F (tilde(x))||_X $\

  #set align(horizon+left)
  に対して，無限次元ガウスの消去法を適用する．

  $A = D F (macron(x)) ^(-1), #h(30pt) phi.alt := D F (macron(x))^(-1) F(tilde(x))$より

  $
  D F (macron(x)) phi.alt = F(tilde(x))\
  $

  #set text(size:20pt)
  $
  mat(
    Pi_N D F (macron(x)) Pi_N, Pi_N D F (macron(x)) (I-Pi_N);
    (I - Pi_N) D F (macron(x)) Pi_N, (I - Pi_N) D F (macron(x)) (I-Pi_N);
  )
  mat(
    Pi_N phi.alt;
    (I-Pi_N) phi.alt
  )=mat(
    Pi_N F(tilde(x)) ;
    (I-Pi_N) F(tilde(x))
  )
  $

]

#slide[
  #set text(size:22pt)
  == 今後の課題
  #set align(horizon + center)

  - 無限次元ガウスの消去法を用いた$Y_0$の展開

  \

  - Juliaを用いたプログラムの実証
]