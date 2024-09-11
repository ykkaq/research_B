# ニュートン法の実装

以下の3手法で実装せよ
- 関数ポインタ　（）
- 関数オブジェクト　（）
- 仮想関数（）


## 仮想関数
```cpp
// 親クラス
class Animal{
  public:
    virtual void sound(){
      std::println("any sound");
    }
}

// 子クラス．ここで継承
class cat : public Animal{
  public:
  // オーバーライド := 仮想関数の上書き != overload
    void sound() override{
      std::println("にゃー");
    }
}

```