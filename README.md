# bayes-mapper

ML系統樹にベイズの事後確率をマップします。

- [bayes-mapper](#bayes-mapper)
  - [Installation](#installation)
    - [Singularityを使う場合](#singularityを使う場合)
    - [Singularityを使わない場合](#singularityを使わない場合)
  - [Usage](#usage)
    - [Singularityの場合](#singularityの場合)
    - [それ以外](#それ以外)
    - [使い方](#使い方)

## Installation

### Singularityを使う場合

`scripts/build-singularity.sh` を実行すると `package/` の中にコンテナ（`bayes-mapper.sif`）が生成されます。

### Singularityを使わない場合

要求：Python3.11

最初に `scripts/init.sh` を実行し，仮想環境の設定を行います。
その後に `scripts/bayes-mapper.sh` を編集します。
5行目の `cd "bayes-mapper.pyのあるディレクトリ"` の部分を適宜変更して下さい。
そして `bayes-mapper.sh` をパスの通った場所に配置すれば完成です。

## Usage

### Singularityの場合

```bash
singularity run bayes-mapper.sif 引数
```

### それ以外

```bash
bayes-mapper.sh 引数
```

### 使い方

ヘルプの表示
```bash
bayes-mapper.sh -h
```

バージョンの表示
```bash
bayes-mapper.sh -v
```

基本形
```bash
bayes-mapper.sh -m <MLツリーのファイル> -b <ベイズのツリーファイル> -o <出力ツリーファイル>
```

BP値やPP値でフィルタリングする場合  
指定した値に満たない場合はサポート値が表示されません。
```bash
bayes-mapper.sh -m <MLツリーのファイル> -b <ベイズのツリーファイル> -o <出力ツリーファイル> --min-bp <BP (0-100)> --min-pp <PP (0-1)>
```
