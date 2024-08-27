// Get Polylux from the official package repository
#import "@preview/polylux:0.3.1": *
//https://typst.app/universe/package/polylux/
#import "@preview/cetz:0.2.2"
//https://github.com/cetz-package/cetz

// Make the paper dimensions fit for a presentation and the text larger
#set page(
  paper: "presentation-16-9"
)
#import themes.simple: *

#set text(
  lang:"ja",
  font: "Harano Aji Mincho",
  size: 22pt,
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
  == 無限次元ガウスの消去法を用いた\
  #rp()の改良

  関根研究室　2131701　齋藤 悠希
]

// 2枚目
#slide[
  == 背景と目的
  #set align(left + horizon)

  - 非線形微分方程式の解の精度保証の保証#rp()の\
    解の保証精度を向上させる

  #set align(center + horizon)
  $+$

  #set align(left + horizon)
  - 無限次元上で定理を進めることで精度が向上する

  /*
  #set align(center + horizon)
  #sym.arrow.b

  #set align(left + horizon)
  - #rp()を無限次元
  */
]



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

#slide[
  == 提案手法
  #set align(horizon)

  $Y_0$の評価式
  $||A F (tilde(x))||_X &lt.eq Y_0$
  に対して，ガウスの無限次元消去法を適用する．$Pi_N$は射影作用素とする．

  $
  A = D F^(-1),  phi.alt := D F^(-1) F(tilde(x))\
  D F phi.alt = F(tilde(x))\
  cases(
    Pi_N D F (Pi_N phi.alt + (I-Pi_N) phi.alt) &= Pi_N F(tilde(x)) ,
    (I-Pi_N) D F (Pi_N phi.alt + (I-Pi_N) phi.alt) &= (I-Pi_N) F(tilde(x)) ,
  )
  $

　/*
  $
  mat(
    Pi_N D F Pi_N, Pi_N D F (I-Pi_N);
    (I - Pi_N) D F Pi_N, (I - Pi_N) D F (I-Pi_N);
  )
  mat(
    Pi_N phi.alt;
    (I-Pi_N) phi.alt
  )=mat(
    Pi_N F(tilde(x)) ;
    (I-Pi_N) F(tilde(x))
  )
  $
  */
]

#slide[
  == 今後の課題
  #set align(horizon + center)

  - ガウスの消去法を用いた$Y_0$の展開

  - Juliaを用いた，プログラムの実証
]