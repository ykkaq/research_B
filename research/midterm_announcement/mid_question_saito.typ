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
  //numbering: "1.",
  numbering: none,
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

//---------------------//
// --- main content ---

#align(center, text([
  #set text(size:20pt)
  中間発表　質問まとめ
]))

#align(right, text([
  齋藤　悠希
]))

#v(20pt)

中間発表で問われた質問とその想定回答を以下に簡単に記述する．

= 質問

== 精度
- 提案手法を用いることでどの程度精度が変わるのか\
  西窪さんの研究/*(#rp()における零点探索手順の削除, 2024/01/24, p.52-54)\*/より，従来手法と比較して，どの程度精度が良くなるかは一概には言えなく，現時点ではわからない．

== julia
- juliaを用いた検証の手法\
  現時点で未定．

- なぜjuliaを使うのか\
  Juliaは，数値計算に関するライブラリが豊富であるため．また，先行研究で用いていたため．

/*
- Juliaとは何か\
  Pythonと同じインタプリタ型の言語．数値計算が得意．

== その他，上記に当てはまらないもの
- 内容の確認
*/
