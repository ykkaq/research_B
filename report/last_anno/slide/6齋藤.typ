// Get Polylux from the official package repository
#import "@preview/polylux:0.3.1": *
//https://typst.app/universe/package/polylux/
#import "@preview/cetz:0.3.1"
//https://github.com/cetz-package/cetz
#import "@preview/showybox:2.0.3":*

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
#let fre() = "Fr"+str.from-unicode(233)+"chet"
#let nk() = "Newton-Kantorovich"
#let rad() = "radii-polynomial approach"
#let vdp()= "van der Pol方程式"
#let infg()= "無限次元ガウスの消去法"

#set underline(
  offset: 3pt
)

#set math.mat(delim: "[")

// --- maintain ---

// タイトルスライド
#slide[
  #set align(horizon + center)
  #set text(size:32pt)
  == 無限次元ガウスの消去法を用いた\ #rad()の改良
  \
  #set text(size:25pt)
  関根研究室　2131701　齋藤 悠希
]


#slide[
  == はじめに
  #set align(horizon+center)

  #showybox(
    body-style: (align: center)
  )[
    微分方程式を計算機で解くとき，\ \
    計算機の資源が #underline[有限] という特徴のために\ \ 
    方程式の解に誤差が発生する．
  ]

  #sym.arrow.b
  
   解の誤差を評価し，精度を保証する．\

  #sym.arrow.r #underline[精度保証付き数値計算]
]

#slide[
  == 背景 - #vdp()
  #set align(horizon)
  #side-by-side[
    #showybox(
      frame: (
        title-color: blue.darken(30%),
        border-color: blue.darken(30%),
        body-color: aqua.lighten(80%)
        ),
      title: [#vdp()],
      title-style: (weight: 600)
    )[
      $
        frac(d^2x, d t^2) - mu (1-x^2)+x=0
      $
      - 未知関数　 : $x(t)$
      - パラメータ : $mu>0$
    ]
  ][
    #figure(
      image("vdp.png"),
      caption:[van der Pol方程式 \ 初期値 $(0,2), mu = 1.0$]
    )
  ]

]

#slide[
  == 背景 - 先行研究
  #set align(horizon)

  #set table(
    stroke: (x,y)=>(
      x: {0pt},
      y: if y!=0 {1pt},
      right: none,
      bottom: 1pt
    )
  )

  #showybox(
    frame: (title-color: blue.darken(30%), border-color: blue.darken(30%), body-color: aqua.lighten(80%)),
    title: [#rad() [1]],
    title-style: (weight: 600)
    )[
      #set text(size:21pt)
      #set align(horizon)
      #side-by-side()[
        #table(
          columns:(auto, auto),
          table.header(
          [],[]
          ),
          align: center,
          inset:9pt,
        )[
          $X,Y$
        ][
          Banach空間
        ][
          $cal(L)(X,Y)$
        ][
          $X arrow.r Y$への\ 有界線形作用素の集合
        ][
          $A^dagger$
        ][
          $cal(L)(X,Y)$の要素
        ][
          $ A $
        ][
          $cal(L)(Y,X)$の要素
        ][
          $F$
        ][
          /*$X arrow.r Y$で*/$C^1$-#fre() 微分\ 可能な作用素
        ]
      ][
        $
        norm(A F (macron(x)))_X &lt.eq Y_0 \
        norm(I-A A^dagger)_(cal(L)(X)) &lt.eq Z_0 \
        norm(A (D F(macron(x))-A^dagger))_(cal(L)(X)) &lt.eq Z_1 \
        norm(A (D F(b)-D F (macron(x))))_(cal(L)(X)) &lt.eq Z_2(r), \
        forall b &in overline(B(tilde(x),r))
        $
      ]
  ]
]


#slide[
  == 背景 - 先行研究
  #set align(horizon + center)

   #showybox(
    frame: (title-color: blue.darken(30%), border-color: blue.darken(30%), body-color: aqua.lighten(80%)),
    title: [#rad() [1] （続き）],
    title-style: (weight: 600)
    )[
      #side-by-side()[
      #set text(size:21pt)
      #set align(horizon)

      radii polynomialを以下で定義する．
      $
      p(r) := Z_2(r)r^2 - (1-Z_1-Z_0)r + Y_0
      $
      $r_0>0$ かつ $p(r_0)<0$ なら，\
      $F(tilde(x))=0$ となる解 $tilde(x)$ が \
      $overline(B(macron(x),r))$ 内に一意に存在する．
      ][
        #showybox()[
          #cetz.canvas(length: 3cm, {
            import cetz.draw: *

            let rate = 1

            circle((0,0), radius: (rate*1.5,rate*1))
            content((rate*1.5*calc.cos(135deg)+0.1, 0*calc.sin(135deg)+0.2),text()[$overline(B(macron(x),r))$])

            circle((0,0), radius:2pt, fill:gray)
            content((-0.1,-0.2),[近似解 $macron(x)$],)

            //circle((0,-rate/2-0.1), radius: (rate/2,rate/3))

            circle((0,-1*rate/2), radius:2pt, fill:gray)
            content((-0.15,-1*rate/2-0.2),[真の解 $ tilde(x)$],)

            line((0,0),(1.5*calc.cos(30deg), calc.sin(30deg)), name:"r")
            content(("r.start", 50%, "r.end"), text()[$ r $],anchor: "south",padding:.1)
          })
        ]
      ]
    ]
    #align(bottom+right)[#text(size:18pt, gray)[参考：[1]高安亮紀, Julia言語を使った精度保証付き数値計算のチュートリアル]]
]

/*
#slide[
  == 背景 - 先行研究
  #set align(horizon)

  西窪の研究[2]では，#rad()において作用素 $A$ を，\
  「$A^dagger$の#underline()[近似逆作用素] $arrow.r$ #underline()[真の逆作用素]」とおき，
  $
    A approx D F(x)^(-1) arrow.r A = D F(x)^(-1)\
    A A^dagger approx I arrow.r A A^dagger = I
  $

  ノルムの計算を簡略化．

  $
  text("例）") quad norm(A(D F(macron(x))-A^dagger)) lt.eq Z_1  arrow.double.r norm(A D F(macron(x)) - I) lt.eq Z_1
  $

  #sym.arrow.r 精度は大きく低下しない，計算時間は短縮

  #align(bottom+right)[#text(size:18pt, gray)[参考：[2]西窪壱華, #rad()における零点探索手順の削除]]
]
*/

#slide[
  == 既存手法と問題点
  #set align(horizon)

  ノルムの計算に，重み付き $l^1$ ノルムを定義．

  #side-by-side()[
    #showybox(frame: (title-color: gray.darken(30%), border-color: black.darken(30%), ),
    title: [重み付き$l^1$ノルム（既存手法）],
    title-style: (weight: 600),
    body-style: (align: center),
    )[
      $
        norm(a)_omega := sum_(k in bb(Z)) abs(a_k) omega_k < oo, (omega_k > 1)
      $
    ]
  ][
    #showybox(frame: (title-color: gray.darken(30%), border-color: black.darken(30%), ),
    title: [$l^1$ノルム（$l^1$空間）],
    title-style: (weight: 600),
    body-style: (align: center),
    )[
      $
        norm(a) := sum_(k in bb(Z)) abs(a_k) < oo
      $
    ]
  ]

  重み付き $l^1$ ノルムでは，$omega_k$があるため，$a$の条件が厳しくなる．\

  #sym.arrow.r 精度保証できる条件が限られる．


  #sym.arrow.r $l_1$空間を使うことで，条件を緩和できる．
]


#slide[
  == 目的
  #set align(center + horizon)
  #set text(size:24pt)

  重みを外した $l_1$ 空間上で計算することで，\ \
  #rad()の適用できる問題の範囲を大きくする
]

#slide[
  == 提案手法 - 概要
  #set align(horizon + center)

  精度保証するために，$norm(D F (macron(x)) F (tilde(x)))$を計算しなければならない．

  #sym.arrow.b

  $l^1$空間上で，$D F(macron(x))$は全単射でなければならない．

  #sym.arrow.b

  無限次元ガウスの消去法[3]を用いて，$D F(macron(x))$が全単射であるか確かめる．

  #align(bottom+right)[#text(size:16pt, gray)[参考：[3]Kouta Sekine, Mitsuhiro T. Nakao, and Shin’ichi Oishi:, "Numerical verification methods for a system of elliptic PDEs, and their software library"]]
]

#slide[
  == 提案手法 - 無限次元ガウスの消去法
  #set align(horizon)

  #showybox()[
    $phi.alt := D F (macron(x))^(-1) F(tilde(x))$とおくと，
    $
    D F(macron(x)) phi.alt = F(tilde(x))
    $

    両辺に作用素$A$を掛け，

    $
    A D F(macron(x)) phi.alt = A F(tilde(x))
    $
  ]
]

#let zero_padding = $0, dots.h.c,0$

#slide[
  == 提案手法 - 無限次元ガウスの消去法
  #set align(horizon)

  射影演算子$Pi_N$より，以下の作用素を定義する．

  #showybox(
    frame: (
      title-color: blue.darken(30%),
      border-color: blue.darken(30%),
      body-color: aqua.lighten(80%)
      ),
    //title: [],
    //title-style: (weight: 600)
  )[
    #set text(size:20pt)
    $
      &T:= Pi_N A D F(macron(x))|_(X_1):X_1 arrow.r X_1,quad &&B:= Pi_N A D F(macron(x))|_(X_2):X_2 arrow.r X_1,\
      &C:= (I-Pi_N) A D F(macron(x))|_(X_1):X_1 arrow.r X_2,quad &&E:= (I-Pi_N) A D F(macron(x))|_(X_2):X_2 arrow.r X_2
    $
  ]

  $D F(macron(x)) phi.alt = F(tilde(x))$は，作用素の定義より，以下に変形できる．
  #showybox()[
    #set math.mat(delim: "(")
    $
    mat(
      T,B;
      C,E;
    )
    mat(
      Pi_N phi.alt;
      (I-Pi_N) phi.alt
    ) = mat(
      Pi_N A F(tilde(x)) ;
      (I-Pi_N) A F(tilde(x))
    )
    $
  ]
]

#slide[
  == 提案手法 - 無限次元ガウスの消去法
  #set align(horizon)

  $S$ を以下のように定義し，$A$ と $D F(macron(x))$ から求められる．
  #showybox(
    body-style: (align: center)
  )[
    $
    S :=& E - C T^(-1) B \
    /*
    =& (I-Pi_N) A_M D F(macron(x)) - ((I-Pi_N) A_M D F(macron(x)))\ & (Pi_N A_M D F(macron(x)))^(-1)(Pi_N A_M D F(macron(x)))
    */
    =& (I-Pi_N) A D F(macron(x)) - ((I-Pi_N) A D F(macron(x))) & (Pi_N A D F(macron(x)))^(-1)(Pi_N A D F(macron(x)))
    $
  ]

  #showybox(
    //body-style: (align: center)
  )[
    #set math.mat(delim: "(")
    #align(center)[$norm(I_(X_2) - S) < 1$]
    となれば，$S$ は全単射となる．($I_(X_2) = I-Pi_N$)
  ]

  #set align(center)
  $T^(-1)$が存在することを確認し，\ $underline(S #text[が全単射であれば，]D F(macron(x)) #text[が全単射となる])$
]

// -------------------------


#slide[
  == 実験手法
  #set align(horizon + center)

  van der Pol方程式 \
  
  #align(center, block[
    #sym.arrow.b 
  ])フーリエ・スペクトル法\

  近似解 $macron(x)$ を導出\

  $
    f(t) = a_0/2 + sum_(n=1)^N (a_n cos n x + b_n sin n x)
  $

  フーリエ係数の次数$N$を変化

]


#slide[
  == 実験手法
  #set align(horizon)

  #figure(
    caption: [実験環境],
    table(
      columns: 2,
      align: center,
      inset:12pt,
      header(
        [*環境*],[*詳細*]
      ),
      [CPU],[12th Gen Intel(R) Core(TM) i7-12700],
      [OS],[Ubuntu 24.04.1 LTS],
      [コンパイラ],[Julia 1.11.2],
      [数値計算ライブラリ],[IntervalArithmetic v0.20.9],
    )
  )
]

#slide[
  == 実験結果
  #set align(horizon+center)
  #show table.cell.where(y: 0): set text(weight: "bold")

  #set table(
    stroke: (x,y)=>(
      x: if x == 1 {1pt} else {0pt},
      y: if y== 0 or y==1 {1pt} else {0pt},
      right: none,
      bottom: 1pt
    )
  )

  #figure(
    //supplement: [],
    //numbering: none,
    caption: [フーリエ係数の次数の変更による$norm(I_(X_2) - S)$の比較],
    table(
      columns: 2,
      align: center + horizon,
      inset:12pt,
      header(
        [次数], [$norm(I_(X_2) - S)$],
      ),
      [#text(size:20pt)[50]],[#text(size:20pt)[0.22815114629236252]],
      [#text(size:20pt)[100]],[#text(size:20pt)[0.11455533660051737]],
      [#text(size:20pt)[150]],[#text(size:20pt)[0.07655718822651922]],
      [#text(size:20pt)[200]],[#text(size:20pt)[0.05749210273025131]],
    )
  )

  #set align(horizon+left)

  - すべての次数条件において，$norm(I_(X_2) - S)<1$ を満たした．\
  - 次数が上がるにつれ，ノルム値が減少．

]

#slide[
  == まとめ
  #set align(horizon)

  - 無限次元ガウスの消去法を用いた#rad()の\ 改良方法を提案した

  - 数値実験での検証により，$l_1$空間上で$D F(macron(x))$が全単射であることが\ わかった

  #sym.arrow.r 提案方法で改良可能であることがわかった
]