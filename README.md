# Orb Transform Library (OTL) — Ver. 1.5

OTL は、星表データをもとにプラネタリウム投影機用の原板（星の穴を開けるための版下）を生成するPython製のツールです。東京大学地文研究会天文部が2007年にJavaで開発した OTL をもとに、Pythonに移植し機能の追加や使いやすさの改善を進めています。

## これは何をするものか

自作の投影機式プラネタリウム(正十二面体などにレンズユニットを複数配置して天球全体を投影するタイプ)を作るとき、各投影機ユニットに取り付ける「原板」に星の穴をどこに開けるかを計算し、SVG または印刷用 PDF として出力します。

1. `SphereReader` が Hipparcos/Tycho の星表、RC3 銀河カタログ、IAU88 星座線カタログを読み込み、天球上の星・星座を取得
2. `Transformer`(または銀河専用の `GalaxyTransformer`)が、指定したドーム半径・投影機配置・レンズ焦点距離などの幾何パラメータに基づいて、天球上の位置を各ユニットの原板上の位置(mm単位のXY座標)に変換
3. `PlateWriterSVG` / `PlateWriterPDF` が変換結果を版下データとして書き出す

次の機能にも対応しています。

- 星座線を構成する恒星だけを、指定した倍率で拡大
- 各原盤が担当する天球領域をポリゴンSVG / PDFとして自動生成

## 構成

| ファイル | 役割 |
|---|---|
| `transformer.py` | 星・星座用の原板を生成するメインスクリプト(エントリポイント) |
| `galaxy_transformer.py` | 天の川(銀河)専用の原板を生成するスクリプト(エントリポイント) |
| `basic_transformer.py` | 上記2つの変換処理に共通する抽象基底クラス |
| `sphere_reader.py` | `hip_main.dat` / `tyc_main.dat` / `rc3.dat` / `IAU88.hlc` を読み込むリーダー |
| `unit_arrangement.py` | 正十二面体をもとにした投影機ユニットの配置計算 |
| `geometry.py` | 錐体(GeneralizedCone)・多面体(Polyhedron)などの幾何計算 |
| `mathvector.py` | 3次元ベクトル演算 |
| `models.py` | 恒星・星座・原板上の位置などのデータクラス |
| `plate_writer.py` | SVG / 印刷用 PDF 出力 |
| `plate_polygon.py` | 原盤ごとの担当星域ポリゴン計算 |
| `config.py` | Java プロパティ形式(`.properties`)の設定ファイル読み込み |
| `starconfig.properties` | `transformer.py` 用の設定サンプル |
| `galaxyconfig.properties` | `galaxy_transformer.py` 用の設定サンプル |

## 必要なデータファイル

`sphere_reader.py` は、スクリプトと同じディレクトリに置かれた以下のカタログファイルを読み込みます(パスはスクリプトの場所basisで解決されるため、実行時のカレントディレクトリには依存しません)。

- `hip_main.dat` — Hipparcos 星表
- `tyc_main.dat` — Tycho 星表
- `rc3.dat` — RC3 銀河カタログ(天の川の疑似的な描画に使用)
- `IAU88.hlc` — IAU88 星座線カタログ

## セットアップ

`tyc_main.dat` は Git LFS で管理しています。初回のみ Git LFS を有効化し、星表データを取得してください。

```bash
git lfs install
git lfs pull
```

続いて Python の依存ライブラリをインストールします。

```bash
pip install -r requirements.txt
```

依存ライブラリは `numpy`、`Pillow`(ユニット配置図の出力に使用)、`ReportLab`(印刷用 PDF の出力に使用)、`pypdf`(PDF 出力の自動検証に使用)です。

## 使い方

対話モードと、設定ファイル(`.properties`)を使う非対話モードの2種類があります。

```bash
# 対話モードで実行(質問に答えながらパラメータを決める)
python3 transformer.py

# 設定ファイルを指定して実行(SVG出力)
python3 transformer.py -f starconfig.properties

# 印刷用PDF形式で出力する場合
python3 transformer.py -PDF -f starconfig.properties

# 設定にかかわらず担当星域ポリゴンも生成する場合
python3 transformer.py --polygons -f starconfig.properties

# 天の川専用原板の生成
python3 galaxy_transformer.py -f galaxyconfig.properties
```

`transformer.py` の対話モードでは、`star_SVG`、`star_pdf`、`polygon_SVG`、`polygon_pdf` をそれぞれ出力するか選択できます。複数の形式を同時に選択することもできます。未入力時は従来どおり `star_SVG` のみを出力します。

主なオプション:

- `-f [configFile]` — プラネタリウムの設定を `.properties` ファイルから読み込む
- `-PDF` — 原板データとポリゴンを印刷用 PDF 形式で出力する(未指定時は SVG)
- `--polygons` — 各原盤の担当星域ポリゴンを原板データと同じ形式で別ファイルに出力する
- `-h`, `-help` — 使い方を表示

設定ファイルの各項目(ドーム半径、投影機とレンズの距離、正十二面体上のユニット配置箇所、星の等級の上限・下限など)は `starconfig.properties` / `galaxyconfig.properties` にコメント付きで記載しています。

設定ファイルの `output.directory` は出力ルートフォルダです。未指定時は `output` が使われ、生成物は種類ごとに次のサブフォルダへ出力されます。

```text
output/
├── star_SVG/      # 恒星原盤の SVG
├── star_pdf/      # 恒星原盤の PDF
├── polygon_SVG/   # 担当星域ポリゴンの SVG
├── polygon_pdf/   # 担当星域ポリゴンの PDF
├── galaxy/        # 天の川専用原盤
└── unit_position/ # ユニット配置画像
```

印刷用 PDF と SVG の用紙サイズは ISO A4（幅 210 mm × 高さ 297 mm）です。印刷時は「用紙に合わせる」を無効にし、100%(実寸)で出力してください。旧 `-PS` オプションは廃止され、使用するとエラーになります。

SVG / PDF ともに、円形原盤の左上の枠外に黒い×印（線幅 0.1 mm）を出力します。印は原盤とともに配置され、北天・南天や白黒反転の設定によらず同じ位置に描かれます。

`plate.frame-radius` は恒星用の円形原盤の枠半径です。星穴の中心座標 `(x, y)` が `x² + y² <= plate.frame-radius²` を満たす場合だけ、SVG / PDF に書き出します。`0` を指定した場合は、用紙サイズと原盤の配置数から枠半径を自動計算します。

### 星座線構成星の拡大

`star.enlarge-rate` に0以上2以下の値を指定します。`1.0` を指定すると、4等星相当の穴は1等級分大きくなります。4等星より明るい星では、固定幅 `ω = 2` のガウス型テーパー `T = exp(-((4 - 等級) / 2)²)` により、実際の等級変化量を `star.enlarge-rate × T` へ滑らかに抑えます。4等星およびそれより暗い星では指定値をそのまま適用し、`0.0` では拡大しません。旧Java版の設定名 `star-EnlargeRate` も使用できます。

拡大する星座は `IAU88.hlc` の `Name:` で指定できます。非対話モードでは `star.enlarge-constellations = Orion, Canis Major` のようにカンマ区切りで指定します。対話モードでは、`Name` を1行に1つ記載したUTF-8テキストファイルのパスを入力します。いずれも空欄の場合は全星座が対象です。存在しない名前は設定エラーになります。

全星座を対象とする場合、対象星は `IAU88.hlc` のHIP番号と、旧版の `hip_constellation_line_star.csv` の補助データを併合した893個です。星座を指定した場合は、星座名を持たない補助データを除き、選択した `Name:` ブロックのHIP番号だけを使用します。HIP番号による照合なので、旧実装の原盤座標の丸め比較より安定しています。

### 担当星域ポリゴン

設定ファイルで `polygon.enabled = yes` とするか、`--polygons` を指定すると生成します。`polygon.samples` は1原盤あたりの頂点数（12以上、既定値180）、`polygon.file-prefix` は出力ファイル名の接頭辞です。

ポリゴンは各原盤座標から天球へ光線を逆投影し、通常の星配置と同じユニット選択規則で境界を二分探索して求めます。枠の円でクリップ済みで、星原盤とは別ファイルに同じ用紙配置・黒い原盤・白い投影領域で出力されます。通常は `polygon_SVG` にSVG出力し、`-PDF` 使用時は `polygon_pdf` にベクトルPDF出力します。

### 天の川原盤の担当領域

`GalaxyTransformer` は星表のICRS赤道座標を銀河座標へ変換し、**銀緯 −20度以上 +20度以下**の星を次の4枚に割り当てます。銀経は0度以上360度未満に正規化し、共有する境界上の星が重複しないよう、各区間の下端を含み上端を含めません。

| 原盤番号 | 銀経の範囲 | 中心銀経 |
|---|---|---|
| N0 | 0度以上60度未満 | 30度 |
| N1 | 60度以上120度未満 | 90度 |
| S0 | 300度以上360度未満 | 330度 |
| S1 | 180度以上240度未満 | 210度 |

銀経120度以上180度未満、240度以上300度未満、および銀緯±20度の外は出力対象外です。N/Sは用紙の上段/下段の識別子で、赤緯・銀緯の正負を意味しません。`plate.column = 1`、`plate.row = 1` の標準設定では、`galaxy-0.svg` にN0/S0、`galaxy-1.svg` にN1/S1を配置し、計4枚の原盤を2用紙に出力します。PDFも同じ配置です。星がない領域も枠を出力します。

各原盤では中心銀経・銀緯0度への光線を原盤中央に合わせ、円筒投影で位置を求めます。`projector-horizontal` と `projector-vertical` は従来どおり赤道座標での距離で、投影機位置は `(0, ±projector-horizontal, ±projector-vertical)` です。横成分の符号は原盤番号0で正、1で負、縦成分はNで正、Sで負です。星座名は領域内だけを出力し、星座線は銀経・銀緯で線形補間して各領域の境界で切り分けます。

座標変換の回転行列はHipparcosのICRS定義に基づきます（[ERFAの参照実装](https://github.com/liberfa/erfa/blob/master/src/icrs2g.c)）。

## 元のJava版との関係・互換性

本プロジェクトは Java 版 OTL (Ver. 1.4) を出発点に、Python版として機能追加・改良を続けています。元のJava版との完全な動作一致を目的とするものではありません。既存の実装との関係で、次の点に留意してください。

- 乱数生成は Python 標準の `random` モジュールを使用しています。Java 版の `java.util.Random` とはアルゴリズムが異なるため、同じシード値(0)でも生成される乱数列そのものは一致しません(統計的性質は同等です)。
- `models.py` の `PlatePosition.__eq__` など、Java 版の挙動をそのまま踏襲した(一見不自然に見える)比較ロジックが一部にあります。詳細はコード中の `NOTE` コメントを参照してください。

## ライセンス

GNU General Public License v2 (or later)。詳細は各ソースファイル冒頭のヘッダーを参照してください。

元のJava版 OTL の著作権表記：

Copyright (C) 2007, 2012-2014 東京大学地文研究会天文部

Python版の2026年の変更分の著作権表記：

Copyright (C) 2026 東京大学地文研究会天文部

本Python版では、元の OTL の著作権表記とライセンスを引き継ぎ、継続的な開発・拡張を行っています。
