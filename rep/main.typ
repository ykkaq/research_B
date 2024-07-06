// --- settings ---

#set page(
  paper: "a4",
  margin: (
    x:20mm,
    y:15mm
  )
)

#set par(
  first-line-indent: 1em,
  //linebreaks: "optimized",
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
  numbering: "1."
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
  number-align: 
)


// --- main content ---

#align(center, text(
  20pt, font: "Harano Aji Gothic"
  )[
  radii polynomial approachにおける\
  無限次元ガウスの消去法
])

#align(center)[
    （指導教員　関根 晃汰　准教授）\
    関根研究室　2131701 齋藤 悠希
]

#line(length: 100%)

#show: rest => columns(
  2,rest
)

= 研究背景と目的
微分方程式の解の精度保証付き数値計算では，Newton法の半局所的収束定理を精度保証付き数値計算に利用したNewton-Kantorovich定理がある．この定理は，連立一次方程式や非線形方程式，偏微分方程式などほとんどの微分方程式に用いることができる．しかし，十分条件が厳しく非線形作用素の解を求める必要があり，数値計算に多くの手間がかかる．

本研究では，Newton-Kantorovich型定理の一部を変更し，数値計算を減らすとともに，既存手法との精度比較を行うことを目的とする．

= Newton-Kantorovichの定理

　$X,Y$をBanach空間，$cal(L) paren.l X,Y paren.r $を$X$から$Y$への有界線形作用素の集合とする．有界線形作用素$A^dagger in cal(L)(X,Y), A in cal(L)(X,Y)$を考え，作用素$F:X arrow.r Y$が$C^1 - "Fr"acute(e)"chet"$微分可能とする．いま，$tilde(x) in X$に対して，正定数$Y_0, Z_0, Z_1$および非減少関数$Z_2(r)(r>0)$が存在して，次に不等式を満たすとする．

$
||A F (tilde(x))||_X &lt.eq Y_0 \
||I-A A^dagger||_(cal(L)(X)) &lt.eq Z_0 \
||A (D F(tilde(x))-A^dagger)||_(cal(L)(X)) &lt.eq Z_1 \
||A (D F(b)-D F (tilde(x)))||_(cal(L)(X)) &lt.eq Z_2(r) \
$