# c++やつ

## オブジェクト指向
- カプセル化
  - private - class外で使えない
  - public  - どこでも使える
  - protected - 1つ下の子クラスまで使える
- 継承
  - 親クラスのメソッド関数・変数を引き継いだ新しいクラス
  - 小さい単位のクラスにする．
- polymorphism（多態性）
  - 抽象化．

### template \<typename _T\>
```c++
template <typename _T>
```
マクロでの文字置換のように，コンパイラが型を推測して文字置換してくれる．

## 2024-04-15
クラス継承とか

### ファイル構造
- main.cpp
  - メインファイル．いろいろ呼び出して具体的計算を指示．
  - ←interval.cpp
- interval.hpp
  - ヘッダファイルその１．


### memo
ひし形継承
