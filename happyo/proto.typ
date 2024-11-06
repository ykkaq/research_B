// Get Polylux from the official package repository
#import "@preview/polylux:0.3.1": *
//https://typst.app/universe/package/polylux/
#import "@preview/cetz:0.3.1"
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
  size: 24pt,
)

#show heading: set text(
  //  headingのフォントを変更
  font: "Harano Aji Gothic",
  size: 34pt,
)


#show figure.where(
  kind: table
): set figure.caption(position: top)
#show figure: set block(breakable: false)

// shortcut
#let fc() = "Fr"+str.from-unicode(233)+"chet"
#let nk() = "Newton-Kantorovich"
#let rp() = "radii-polynomial approach"
#set underline(
  offset: 3pt
)
#set par(leading: 0.8em)


// --- maintain ---

// タイトルスライド
#slide[
  #set align(horizon + center)
  #set text(size:32pt)
  == 3年第n回ゼミ
  \
  #set text(size:25pt)
  関根研究室　2131701　齋藤 悠希
]

#slide[
  == おしながき
  #set align(horizon)

  #grid(
    columns: 2,
    [
      - 自己紹介
      - 知識
        - 有界線形作用素
        - #fc()微分
      - 研究テーマの紹介
    ],[
      #align(right + bottom)[
        #image("img/challenge_mokuhyou_businessman.png", width: 50%)
      ]
    ]
  )

]

#slide[
  == 自己紹介
  #set align(horizon)

  #grid(
    columns: 2,
    rows: 2,
    [
      - 名前
        - 齋藤 悠希
    ],
    (
      image("img/x.png", width: 50%)
    ),
    [
      - 趣味
        - ゲーム (Valorantとか)
        - プログラミングとか
        - ドライブ
    ],
    [
      #grid(
        columns: 3,
        image("img/valo.png", width: 70%),
        image("img/linux_penguin.jpg", width: 50%),
        image("img/rust.png", width: 90%),
      )
    ]
  )
]


#slide[
  == おしながき
  #set align(horizon)

  #grid(
    columns: 2,
    [
      - #emph(text(gray)[自己紹介])
      - 知識
        - 有界線形作用素
        - #emph(text(gray)[#fc()微分])
      - #emph(text(gray)[研究テーマの紹介])
    ],[
      #align(right + bottom)[
        #image("img/mathsha.png", width: 50%)
      ]
    ]
  )
]

#slide[
  == 有界線形作用素

  #set align(horizon + center)

  *有界で線形な作用素*

  #image("img/hatena.png", height: 50%)

  $&arrow.r.stroked$ 順番に考えてみる

]

#slide[
  == 有界線形作用素
  #set align(horizon + center)

  #cetz.canvas({
    import cetz.draw: *
    circle((0,0),radius:(4,1))
    circle((0,0.5),radius:(6,2))
    circle((0,1),radius:(8,3))

    content((0,0),[有界線形作用素])
    content((0,1.7),[線形作用素])
    content((0,3.2),[作用素])
  })

  \ 図 有界で線形な作用素のベン図
]

#slide[
  == 有界線形作用素 - 作用素

  #figure(image("img/operator.png", height:85%))
]

#slide[
  == 有界線形作用素 - 線形作用素

  #set align(horizon)

    #emph(text(font: "Harano Aji Gothic")[*定義*])　線形作用素

  - $X,Y$ : 線形空間，$A$ : $X &arrow.r Y$への作用素
  - 任意の$u,v &in cal(D)(A) &subset X$と$alpha &in bb(K)$， に対し
    - $cal(D)(A)$が$X$の線形部分空間
    - $A(u+v) = A u + A v$
    - $A( alpha u) = alpha A u$

  を満たす，$A$のこと．
]


#slide[
  == 有界線形作用素 - 線形作用素

  #set align(center+horizon)

  #image("img/linear_venn.png", height: 85%)
]

#slide[
  == 有界線形作用素 - 線形作用素

  #set align(horizon)

    #emph(text(font: "Harano Aji Gothic")[*定義*])　線形作用素

  - #emph(text(gray)[$X,Y$ : 線形空間，$A$ : $X &arrow.r Y$への作用素])
  - #emph(text(gray)[任意の$u,v &in cal(D)(A) &subset X$と$alpha &in bb(K)$， に対し])
    - #emph(text(gray)[$cal(D)(A)$が$X$の線形部分空間])
    - $A(u+v) = A u + A v$
    - $A( alpha u) = alpha A u$
    #emph(text(gray)[を満たす，$A$のこと．])

  #sym.arrow.r.stroked つまり，#underline[足し算と定数倍が保存される]
]

#slide[
  == 有界線形作用素 - 線形作用素

  #set align(center+horizon)

  #image("img/linear_tree.png", height: 85%)
]

#slide[
  == 有界線形作用素 - 有界線形作用素
  #set align(horizon)

    #emph(text(font: "Harano Aji Gothic")[*定義*])　有界線形作用素

  - $X,Y$を線形空間，$A$ : $X &arrow.r Y$への作用素
  - $u &in cal(D)(A) $に対し,
    - $norm(A u)_Y &lt.eq K norm(u)_X$
    を満たす，正の定数$K$が存在するとき，\
    $A$を有界な作用素と呼ぶ．

]

#slide[
  == 有界線形作用素 - 線形作用素

  #set align(center+horizon)

  #image("img/linear_eg.png", height: 85%)
]

#slide[
  == 有界線形作用素

  #set align(horizon)

  つまり，有界線形作用素$A$とは，

  - *有界で*（上限と下限があって）
    - $norm(A u)_Y &lt.eq K norm(u)_X$
  - *線形な*
    - $A(u+v) = A u + A v$
    - $A( alpha u) = alpha A u$
  - *作用素*（写像みたいなやつ）
]

#slide[
  == おしながき
  #set align(horizon)

  #grid(
    columns: 2,
    [
      - #emph(text(gray)[自己紹介])
      - 知識
        - #emph(text(gray)[有界線形作用素])
        - #fc()微分
      - #emph(text(gray)[研究テーマの紹介])
    ],[
      #align(right + bottom)[
        #image("img/mathsha.png", width: 50%)
      ]
    ]
  )
]

#slide[
  == #fc()微分
  #set align(horizon)


    #emph(text(font: "Harano Aji Gothic")[*定義*])　#fc()微分

  - $X,Y$をBanach空間，開部分集合$ U subset X$
  - 定義域を$cal(D)(f) = U$とする，$U$から$Y$への作用素$f$は$U$上で連続
  - ある点$v in U$に対し，$v+h$となる任意の$h in X$について
  $
    frac(norm(f(v+h)-f(v)-f'[v]h)_Y,norm(h)_X) arrow.r 0, (h arrow.r 0)
  $
    を満たす線形作用素$f'[v] in cal(B)(X,Y)$を$f$の点$v$における#fc()微分と呼ぶ．
]

#slide[
  == #fc()微分

  #set align(center+horizon)

  #image("img/hatena.png", height: 50%)

  #sym.arrow.r.stroked 順番に考えてみる
]


#slide[
  == #fc()微分

  #set align(center+horizon)

  #fc()微分は，全微分の拡張

  #sym.arrow.b.stroked


  全微分は？

  #sym.arrow.b.stroked


  微分は？
]

#slide[
  == #fc()微分 - 微分
  #set align(horizon)

    #emph(text(font: "Harano Aji Gothic")[*定義*])　微分可能

  関数$f(x)$が微分可能であるとは，

  $
    frac(f(x+h)-f(x),h) arrow.r c,  (h arrow.r 0)
  $

  となる $c$ が存在すること．
]

#slide[
  == #fc()微分 - 微分
  #set align(horizon)

  $f(x) = x^2$のときの $c$ は？

  $
    c &= frac(f(x+h)-f(x),h) \
    &= frac((x^2 + 2x h + h^2)-x^2,h) \
    &= 2x + h \
    &arrow.r 2x quad (h arrow.r 0)
  $

  これより，$f'(x)=c=2x$
]

#slide[
  == #fc()微分 - 微分
  #set align(horizon)

  $
    frac(f(x+h)-f(x),h) arrow.r c, (h arrow.r 0)
  $

  より，

  $
    frac(f(x+h)-f(x)-c h,h) arrow.r 0, (h arrow.r 0)
  $
]

#slide[
  == #fc()微分 - 微分
  #set align(horizon)

  #align(center,image("img/diff0.png", width: 61%))
]

#slide[
  == #fc()微分 - 全微分
  #set align(horizon)

    #emph(text(font: "Harano Aji Gothic")[*定義*])　全微分可能

  関数$f(x)$が微分可能であるとは，

  $
    frac(f(x+h)-f(x)- c dot h, norm(h)) arrow.r 0, (h arrow.r 0)
  $

  となる $c$ が存在すること．\
  ($x,c,h$はベクトル)
]

#slide[
  == #fc()微分 - 全微分
  #set align(horizon)

  #align(center,image("img/diff1.png", width: 61%))
]

#slide[
  == #fc()微分 - #fc()微分
  #set align(horizon)


    #emph(text(font: "Harano Aji Gothic")[*定義*])　#fc()微分（再掲）

  - $X,Y$をBanach空間，開部分集合$ U subset X$
  - 定義域を$cal(D)(f) = U$とする，$U$から$Y$への作用素$f$は$U$上で連続
  - ある点$v in U$に対し，$v+h$となる任意の$h in X$について
  $
    frac(norm(f(v+h)-f(v)-f'[v]h)_Y,norm(h)_X) arrow.r 0, (h arrow.r 0)
  $
    を満たす線形作用素$f'[v] in cal(B)(X,Y)$を$f$の点$v$における#fc()微分と呼ぶ．
]

#slide[
  == #fc()微分 - #fc()微分
  #set align(horizon)

  - $X,Y$をBanach空間，開部分集合$ U subset X$
  - 定義域を$cal(D)(f) = U$とする，$U$から$Y$への作用素$f$は$U$上で連続

  #set align(center)
  #sym.arrow.b.stroked


  #rect(
    radius: 7pt,
    inset: 15pt,
    [
      定義域の指定 \
      関数を作用素に拡張\
      $f$が連続 $arrow.r$ 微分できることの明示\
    ],
  )
]

#slide[
  == #fc()微分 - #fc()微分
  #set align(horizon)

  - ある点$v in U$に対し，$v+h$となる任意の$h in X$について
  $
    frac(norm(f(v+h)-f(v)-f'[v]h)_Y,norm(h)_X) arrow.r 0, (h arrow.r 0)
  $
    を満たす線形作用素$f'[v] in cal(B)(X,Y)$を...

  #set align(center)
  #sym.arrow.b.stroked

  #rect(
    radius: 7pt,
    inset: 15pt,
    [
      ほぼ全微分と同じ
    ],
  )
]

#slide[
  == #fc()微分 - #fc()微分
  #set align(horizon)

  （さっきの式）を満たす線形作用素$f'[v] in cal(B)(X,Y)$を$f$の点$v$における#fc()微分と呼ぶ．

  #set align(center)
  #sym.arrow.b.stroked

  #rect(
    radius: 7pt,
    inset: 15pt,
    [
      微分の結果を$f'[v]$とするよ（$c$と意味が同じ）
    ],
  )

  #line(length: 100%)

  $cal(B)(X,Y) dots.c$定義域が$X$の全体となる有界な線形作用素全体の集合 \
  #sym.arrow.r.curve 定義域が$X$となる有界線形作用素を，すべて集めたやつ
]

#slide[
  == #fc()微分 - #fc()微分
  #set align(horizon)

  #grid(columns: 2,
    gutter: 14pt,
    align(center,image("img/diff2.png", width: 110%)),
    [
      #rect(
        inset: 15pt,
        [
          線形作用素$f'[v]$を\
          $f$ \の点$v$ における\
          #fc()微分と呼ぶ．
        ],
      )
    ]
  )

]

#slide[
  == #fc()微分（再掲）
  #set align(horizon)

  #emph(text(font: "Harano Aji Gothic")[*定義*])　#fc()微分

  - $X,Y$をBanach空間，開部分集合$ U subset X$
  - 定義域を$cal(D)(f) = U$とする，$U$から$Y$への作用素$f$は$U$上で連続
  - ある点$v in U$に対し，$v+h$となる任意の$h in X$について
  $
    frac(norm(f(v+h)-f(v)-f'[v]h)_Y,norm(h)_X) arrow.r 0, (h arrow.r 0)
  $
    を満たす線形作用素$f'[v] in cal(B)(X,Y)$を$f$の点$v$における#fc()微分と呼ぶ．
]

#slide[
  == おしながき
  #set align(horizon)

  #grid(
    columns: 2,
    [
      - #emph(text(gray)[自己紹介])
      - #emph(text(gray)[#fc()知識])
        - #emph(text(gray)[#fc()有界線形作用素])
        - #emph(text(gray)[#fc()微分])
      - 研究テーマの紹介
    ],[
      #align(right + bottom)[
        #image("img/mathsha.png", width: 50%)
      ]
    ]
  )
]

#slide[
  == 研究テーマの紹介 - 目的

  #set align(horizon + center)

  *#rp()を改善する*
]

#slide[
  == #rp()とは
  #set align(horizon)

  #underline([非線形方程式の近似解の精度保証を行う定理])

  #v(15pt)

  #emph(text(font: "Harano Aji Gothic")[*定理*])　#rp()
  - $X,Y$を，Banach空間
  - $cal(L)(X,Y)$を，有界線形作用素全体の集合
      - 定義域が$X$となる有界線形作用素を，すべて集めたやつ
  - 作用素$F:X arrow.r Y$が，$C^1$-#fc()可能
    - $C^1$-#fc()可能  ...  1回微分可能　かつ　$F'$が連続
]

#slide[
  == #rp()とは
  #set align(horizon)

  $tilde(x) in X$に対して、\ 正定数$Y_0,Z_0,Z_1$および、非減少関数$Z_2(r)(r>0)$が存在して
  $
  ||A F [tilde(x)]||_X &lt.eq Y_0 \
  ||I-A A^dagger||_(cal(L)(X)) &lt.eq Z_0 \
  ||A (F'[tilde(x)]-A^dagger)||_(cal(L)(X)) &lt.eq Z_1 \
  ||A (F'[b]-F' (tilde(x)))||_(cal(L)(X)) &lt.eq Z_2(r), quad forall b &in overline(B(tilde(x),r))
  $
  を満たすとする．
]

#slide[
  == #rp()とは
  #set align(horizon)

  このとき，radii polynomialを以下で定義する．
  $
  p(r) := Z_2(r)r^2 - (1-Z_1-Z_0)r + Y_0
  $

  これに対し、$p(r_0)<0$となる$r_0>0$が存在するならば、$F(x)=0$を満たす解 $tilde(x)$ が $overline(B(x,r))$ 内に一意に存在する．
]

#slide[
  == 研究テーマの紹介 - 背景
  // table

  例えば...　$||A F (tilde(x))||_X &lt.eq Y_0$を計算したい．

  #set align(center + horizon)

  従来手法では，$A$を$F'[macron(x)]^(-1)$の#underline([近似])として，計算する．

  #sym.arrow.r.curve 精度があまり良くない

  #sym.arrow.b.stroked

  #underline([無限次元のガウスの消去法])を使って$F'[macron(x)]^(-1)$を計算し，\
  精度を良くしよう
]

//------------------


#slide[
  == 研究テーマの紹介 - 目的・手法
  #set align(center + horizon)
  #set text(size: 28pt)

  #underline[*近似部分に無限次元ガウスの消去法を用いる*]\ #sym.arrow.r.curve 精度の改善
]
