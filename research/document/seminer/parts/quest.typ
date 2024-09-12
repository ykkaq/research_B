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
  #set text(size:15pt)
  質問 #datetime.today().display()
]))

/*
#align(right, text([
  齋藤　悠希
]))

*/

#line(length:100%)
#v(20pt)

#let dfx() = $F' [macron(x)]$
#let fx() = $F \(macron(x)\)$

$
#dfx() phi.alt = #fx()
$

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