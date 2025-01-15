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
  #set align(horizon)

  微分方程式を計算機で解くとき，\
  計算機の資源が有限という特徴のために\
  方程式の解に誤差が発生する．

  #sym.arrow.r 解の誤差を評価し，精度を保証する．

  #underline[精度保証付き数値計算]
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
      - $x(t)$ : 未知関数
      - $mu>0$ : 非線形の減衰の強さを表すパラメータ
    ]
  ][
    #figure(
      image("vdp.png"),
      caption:[van der Pol方程式 \ 初期値 $(x,y)=(0,2), mu = 1.0$]
    )
  ]

  //フーリエ・スペクトル法 #sym.arrow.r 近似周期解 #sym.arrow.r 方程式の精度保証
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
      #showybox()[
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
            content((rate*1.5*calc.cos(135deg)+0.1, 0*calc.sin(135deg)),text()[$overline(B(macron(x),r))$])

            circle((0,0), radius:2pt, fill:gray)
            content((-0.15,-0.1),[$ macron(x)$],)

            //circle((0,-rate/2-0.1), radius: (rate/2,rate/3))

            circle((0,-1*rate/2), radius:2pt, fill:gray)
            content((-0.15,-1*rate/2),[$ tilde(x)$],)

            line((0,0),(1.5*calc.cos(30deg), calc.sin(30deg)), name:"r")
            content(("r.start", 50%, "r.end"), text()[$ r $],anchor: "south",padding:.1)
          })
        ]
      ]
    ]
    #align(bottom+right)[#text(size:18pt, gray)[参考：[1]高安亮紀, Julia言語を使った精度保証付き数値計算のチュートリアル]]
]

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

#slide[
  == 既存手法と問題点
  #set align(horizon)

  ノルムの計算に，L1ノルムに重み $omega_k$ を付ける．

  #side-by-side()[
    #showybox(frame: (title-color: gray.darken(30%), border-color: black.darken(30%), ),
    title: [重み付きノルム],
    title-style: (weight: 600),
    body-style: (align: center),
    )[
      $
        norm(a) := sum_(k in bb(Z)) abs(a_k) omega_k < oo
      $
    ]
  ][
    #showybox(
    frame: (title-color: gray.darken(30%), border-color: black.darken(30%), ),
    title: [L1ノルム],
    title-style: (weight: 600),
    body-style: (align: center),
    )[
      $
        norm(a) := sum_(k in bb(Z)) abs(a_k)
      $
    ]
  ]

  ノルムの値は，「重み付きノルム」 $>$ 「L1ノルム」

  #sym.arrow.r $a$ に制限がかかり，条件に合わない問題が出てくる．

  #sym.arrow.r L1ノルムで計算すればよい．
]

#slide[
  == 目的
  #set align(center + horizon)
  //#set text(size:24pt)

  #rad()の適用できる問題の範囲を増やす

  #sym.arrow.b

  無限次元ガウスの消去法を用いて，重みを削除する

  #sym.arrow.b

  無限次元ガウスの消去法を用いた計算が可能か検証する
]

#let zero_padding = $0, dots.h.c,0$
/*
#slide[
  == 提案手法
  $D F(x)$と$A_M$を定義する．
  #set align(horizon)

  #showybox(
    frame: (
      title-color: blue.darken(30%),
      border-color: blue.darken(30%),
      body-color: aqua.lighten(80%)
      ),
    title: [ヤコビ行列$D F(x)$],
    title-style: (weight: 600)
  )[
    /*周期とフーリエ係数列*/近似周期解$(omega, a)$より，$ x = (omega, underbrace(#zero_padding, "M"), a, underbrace(#zero_padding, "M"))$と定め，
    $
      D F(x) = mat(
        0, dots.h.c, 1, dots.h.c;
        dots.v, , dots.v;
        partial_omega f_k, dots.h.c, partial_(a_j) f_k, dots.h.c;
        dots.v, , dots.v;
        augment: #(hline: 1, vline: 1)
      )
    $
  ]
]

#slide[
  == 提案手法
  #set align(horizon)
  #showybox(
    frame: (
      title-color: blue.darken(30%),
      border-color: blue.darken(30%),
      body-color: aqua.lighten(80%)
      ),
    title: [作用素$A_M$],
    title-style: (weight: 600)
  )[
    $macron(x) = (omega, underbrace(#zero_padding, "M"), a, underbrace(#zero_padding, "M"))$と定め，
    $
      A_M = mat(
        D F(macron(x))^(-1), 0, dots.h.c,  dots.h.c;
        0, lambda_N^(-1),,0;
        dots.v, , lambda_(N+1)^(-1), ;
        dots.v, 0, , dots.down;
        augment: #(hline: 1, vline: 1)
      ),quad
      (lambda_k := -k^2 omega^2 - i mu k omega + 1)
    $
  ]
]
*/

#slide[
  == 提案手法
  #set align(horizon)
  //#set text(size:20pt)
  //$norm(A F(macron(x))) lt.eq Y_0$をもとに，無限次元ガウスの消去法を用いて求める

  #showybox()[
  /*$A = D F(x)^(-1)$とおき，*/$phi.alt := D F (macron(x))^(-1) F(tilde(x))$とおくと，
  $
  D F(macron(x)) phi.alt = F(tilde(x))
  $
  ]
  無限次元ガウスの消去法[3]を用いて，$D F(macron(x))^(-1)$の全単射性を確かめる．

  #set align(center)

  //$D F(macron(x))^(-1)$の全単射性 $arrow.double.r$ $phi.alt$が計算可能

  #align(bottom+right)[#text(size:16pt, gray)[参考：[3]関根晃太, 中尾充宏, 大石進一, "Numerical verification methods for a system of elliptic
PDEs, and their software library"]]
]

#slide[
  == 提案手法
  #set align(horizon)

  射影演算子$Pi_N$と作用素$A$より，以下の作用素を定義する．

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
    /*
    $
      &T:= Pi_N A_M D F(macron(x))|_(X_1):X_1 arrow.r X_1,quad &&B:= Pi_N A_M D F(macron(x))|_(X_2):X_2 arrow.r X_1,\
      &C:= (I-Pi_N) A_M D F(macron(x))|_(X_1):X_1 arrow.r X_2,quad &&E:= (I-Pi_N) A_M D F(macron(x))|_(X_2):X_2 arrow.r X_2
    $
    */
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
      /*
      Pi_N A_M F(tilde(x)) ;
      (I-Pi_N) A_M F(tilde(x))
      */
    )
    $
  ]

  
]

/*
#slide[
  == 提案手法

  #set align(horizon)
  以下のように変形する．
  #showybox(
    body-style: (
      align: center
    )
  )[
    $
      norm(I_(X_2) - S) &= norm(I_(X_2) - (D - C T^(-1) B)) \
      &lt.eq norm(I_(X_2) + D) + norm(C) norm(T^(-1)) norm(B)\
      &< 1
    $
  ]
]
*/

#slide[
  == 提案手法
  #set align(horizon)
  #set text(size:20pt)

  #showybox(
    //body-style: (align: center)
  )[
    #set math.mat(delim: "(")
    $S := D - C T^(-1) B$ としたとき，\
    #align(center)[$norm(I_(X_2) - S) < 1$]
    となれば，$S$は全単射となる．
  ]

  $S$は，$A$と$D F(macron(x))$から求められる．
  #showybox(
    body-style: (align: center)
  )[
    $
    S :=& D - C T^(-1) B \
    /*
    =& (I-Pi_N) A_M D F(macron(x)) - ((I-Pi_N) A_M D F(macron(x)))\ & (Pi_N A_M D F(macron(x)))^(-1)(Pi_N A_M D F(macron(x)))
    */
    =& (I-Pi_N) A D F(macron(x)) - ((I-Pi_N) A D F(macron(x)))\ & (Pi_N A D F(macron(x)))^(-1)(Pi_N A D F(macron(x)))
    $
  ]

  /*
  #set align(horizon)
  #showybox(
    body-style: (
      align: center
    )
  )[
    $
      A_M D F(macron(x))
      & = mat(
        T^(-1), 0;
        0, lambda;
        augment: #(hline: 1, vline: 1)
      )
      mat(
        T, M_(12);
        M_(21), M_(22);
        augment: #(hline: 1, vline: 1)
      )=mat(
        T^(-1)T, T^(-1)M_(12);
        lambda M_(21), lambda M_(22);
        //augment: #(hline: 1, vline: 1)
      )\
      & = mat(
        I, B;
        C, E;
      )
    $
  ]
  */
]

// -------------------------

#slide[
  == 実行環境
  #set text(size:21pt)
  #set align(horizon + center)

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
      [OS],[Ubuntu 24.04.1],
      [コンパイラ],[Julia 1.11.2]
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
      [50],[0.22815114629236252],
      [100],[0.11455533660051737],
      [150],[0.07655718822651922],
      [200],[0.05749210273025131],
    )
  )
]

#slide[
  == まとめ
  #set align(horizon)

  - 無限次元ガウスの消去法を用いた#rad()の改良手法を提案した

  - 数値実験での検証により，提案手法で改良可能であることが確かめられた．
]