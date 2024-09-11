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
  size: 22pt,
)

#show heading: set text(
  //  headingのフォントを変更
  font: "Harano Aji Gothic",
  size: 30pt,
)


#show figure.where(
  kind: table
): set figure.caption(position: top)

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

#slide[
  == 背景
  #set align(horizon)
  #underline[*radii polynomial approach* 非線形方程式の解の精度保証に使われる定理]

  \

  $tilde(x) in X$に対して、正定数$Y_0,Z_0,Z_1$および、非減少関数$Z_2(r)(r>0)$が存在して、次の式を満たすとする．

  $
  ||A F (tilde(x))||_X &lt.eq Y_0 \
  ||I-A A^dagger||_(cal(L)(X)) &lt.eq Z_0 \
  ||A (D F(tilde(x))-A^dagger)||_(cal(L)(X)) &lt.eq Z_1 \
  ||A (D F(b)-D F (tilde(x)))||_(cal(L)(X)) &lt.eq Z_2(r), #h(10pt) forall b &in overline(B(tilde(x),r))
  $
]

#slide[
  == 背景
  #set align(horizon)
  このとき，radii polynomialを以下で定義する．
  $
  p(r) := Z_2(r)r^2 - (1-Z_1-Z_0)r + Y_0
  $

  これに対し、$p(r_0)<0$となる$r_0>0$が存在するならば、$F(x)=0$を満たす解 $tilde(x)$ が $overline(B(x,r))$ 内に一意に存在する．
]

#slide[
  ==  目的
  // table
  #set align(center + horizon)

  #rect(
    radius: 7pt,
    inset: 10pt,
    [
    $A$を$D F (macron(x))^(-1)$の近似した作用素で代替 \ #sym.arrow.r 無限次元ガウスの消去法で $D F (macron(x))^(-1)$ を計算
    ],
  )

  #figure(
    supplement: [],
    numbering: none,
    caption: [従来手法と提案手法の簡単な比較],
    table(
      columns: 3,
      align: center,
      inset:9pt,
      header(
        [],
        [*計算*],
        [*精度*],
      ),
      [従来手法],
      [簡略化],
      cell(
        align: center,
        fill: blue.lighten(80%),
        [悪い],//[#sym.arrow.b],
      ),
      [提案手法],
      [無限次元ガウスの\ 消去法],
      cell(
        align: center,
        fill: red.lighten(80%),
        [良い],//[#sym.arrow.t],
      ),
    )
  )
]

#slide[
  == 目的
  #set align(center + horizon)
  #set text(size: 28pt)

  #underline[*簡略化部分を無限次元ガウスの消去法で計算*]\ #sym.arrow.r 精度の改善
]


#slide[
  == 提案手法

  #set align(horizon + left)
  #rp()の一部， $Y_0$の評価式 \

  #set align(center)
  $||A F (tilde(x))||_X <= Y_0$\

  #set align(horizon+left)
  に対して，無限次元ガウスの消去法を適用する．

  $A = D F (macron(x)) ^(-1), #h(30pt) phi.alt := D F (macron(x))^(-1) F(tilde(x))$として

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
  == 今後の課題
  #set align(horizon + left)

  - 無限次元ガウスの消去法を用いた$Y_0$の展開

  \

  - Juliaを用いたプログラムの実証
]