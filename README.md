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
| `galaxy_profile.py` | RC3銀河の指数関数の表面輝度モデルの尺度を計算 |
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
- `rc3.dat` — RC3 銀河カタログ(系外銀河を点の集合で描画するために使用)
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

### 暗い恒星とRC3銀河の描画

`transformer.py` では、次の2項目を独立して指定できます。対話モードでも、それぞれ別の質問で選択します。

| 設定 | `yes` の動作 | 省略時 |
|---|---|---|
| `star.under-minimum` | `star.minimum` より暗い恒星を確率的に採用し、最微等級の星として描画 | `no` |
| `rc3.enabled` | RC3カタログの系外銀河を、最微等級の点の集合で描画 | `no` |

例えば、暗い恒星だけを追加する場合は `star.under-minimum = yes` と `rc3.enabled = no`、RC3銀河だけを追加する場合はその逆を指定します。銀河の点の等級には `star.minimum` を使い、銀河の明るさに応じて個数を決めます。

以前は `star.under-minimum = yes` で両方の処理が有効になっていました。既存の設定ファイルで同じ描画を続けるには、`rc3.enabled = yes` を追加してください。同梱の `starconfig.properties` には追加済みです。天の川専用の `galaxy_transformer.py` は従来どおり両方の処理を有効にします。

#### RC3銀河の点群の広がり

銀河の表面輝度を `I(R) = I0 exp(-R/h)` とする楕円状のモデルを使います。`q` は長径/短径比、`R = sqrt(x² + (q y)²)`、`h` は長軸方向の尺度です。面積要素を含めた半径分布 `p(R) = R exp(-R/h) / h²`（形状母数2のガンマ分布）から点を選び、方向角は一様に選びます。従来の対数正規分布で生じていた中心の密度のへこみを解消します。

尺度は次の優先順で求めます。

1. **Aeがある場合:** 総B光量の半分を含む円形開口の直径Aeを読み、円形開口内のモデル光量が半分になるよう、軸比を考慮した数値積分と二分法でhを求めます。Aeを楕円の半長径としては扱いません。円形銀河では `h = (Ae/2) / 1.67834699` です。この経路ではD25は不要です。
2. **Aeがなく、D25と総B等級がある場合:** 総B光量と、D25の輪郭でのB表面輝度25等級/平方秒角を合わせます。D25の半長径を秒角単位で `a25` とすると、`t = a25/h` について `t² exp(-t) = (2π a25²/q) × 10^(0.4(BT-25))` を解きます。二解がある場合は、中心に集中する `t >= 2` の解を選びます。これはモデル上の選択であり、観測から一意に決まるものではありません。
3. **D25からの解がない場合、またはB等級が不明な場合:** `h = a25/3` とする近似に切り替えます。この近似ではD25楕円内にモデル光量の約80.1%が入りますが、輪郭の表面輝度25との一致は保証しません。AeもD25もなければ、その銀河を読み飛ばします。

RC3のAeは177–180列、軸比の常用対数は162–165列を使います。D25とAeは0.1分角単位の直径の常用対数なので、`3 × 10^格納値` で秒角単位の半径に戻します。軸比や位置角が欠損する銀河は引き続き読み飛ばします。

点の個数は従来どおり総V等級から `floor(2.512^(m_min - V))` で決めます。BT欄がV等級を表す場合は、色指数B−VがあればB等級を復元して尺度の計算に使い、なければAeまたは上記の近似を使います。BとVで同じ形状を仮定し、銀河内の色の変化は再現しません。指数関数の円盤はSérsic分布のn=1に相当する簡略モデルで、渦巻腕やバルジなどの個別構造は表現しません。

D25で点を打ち切らず、外側の淡い光も生成します。位置角と天球への変換、原盤への投影は共通の処理を使います。乱数のシードは従来どおり0ですが、分布を変更したため旧版とは点配置が変わります。

仕様の出典: [RC3のカタログ欄とAeの定義](https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/VII/155?format=html&tex=true)、[指数関数・Sérsic表面輝度モデル](https://ned.ipac.caltech.edu/level5/March05/Graham/Graham2.html)。

### 星座線構成星の拡大

`star.enlarge-rate` に0以上2以下の値を指定します。`1.0` を指定すると、4等星相当の穴は1等級分大きくなります。4等星より明るい星では、固定幅 `ω = 2` のガウス型テーパー `T = exp(-((4 - 等級) / 2)²)` により、実際の等級変化量を `star.enlarge-rate × T` へ滑らかに抑えます。4等星およびそれより暗い星では指定値をそのまま適用し、`0.0` では拡大しません。旧Java版の設定名 `star-EnlargeRate` も使用できます。

拡大する星座は `IAU88.hlc` の `Name:` で指定できます。非対話モードでは `star.enlarge-constellations = Orion, Canis Major` のようにカンマ区切りで指定します。対話モードでは、`Name` を1行に1つ記載したUTF-8テキストファイルのパスを入力します。いずれも空欄の場合は全星座が対象です。存在しない名前は設定エラーになります。

全星座を対象とする場合、対象星は `IAU88.hlc` のHIP番号と、旧版の `hip_constellation_line_star.csv` の補助データを併合した893個です。星座を指定した場合は、星座名を持たない補助データを除き、選択した `Name:` ブロックのHIP番号だけを使用します。HIP番号による照合なので、旧実装の原盤座標の丸め比較より安定しています。

### 担当星域ポリゴン

設定ファイルで `polygon.enabled = yes` とするか、`--polygons` を指定すると生成します。`polygon.samples` は1原盤あたりの頂点数（12以上、既定値180）、`polygon.file-prefix` は出力ファイル名の接頭辞です。

ポリゴンは各原盤座標から天球へ光線を逆投影し、通常の星配置と同じユニット選択規則で境界を二分探索して求めます。枠の円でクリップ済みで、星原盤とは別ファイルに同じ用紙配置・黒い原盤・白い投影領域で出力されます。通常は `polygon_SVG` にSVG出力し、`-PDF` 使用時は `polygon_pdf` にベクトルPDF出力します。

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
