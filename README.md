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
- Gaia DR3の微光星の光量を集計し、指定の最小穴径・金属幅を満たす天の川原盤を生成

## 構成

| ファイル | 役割 |
|---|---|
| `transformer.py` | 星・星座用の原板を生成するメインスクリプト(エントリポイント) |
| `galaxy_transformer.py` | 天の川(銀河)専用の原板を生成するスクリプト(エントリポイント) |
| `gaia_catalog.py` | Gaia DR3の光量集計・取得・キャッシュ検証・CSV読み込み |
| `galaxy_etching.py` | 局所光量を加工可能な円形穴へ変換し、光量誤差を集計 |
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
- `rc3.dat` — RC3 銀河カタログ(系外銀河の光を疑似星群として描画)
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

# Gaia天の川原盤: 初回は光量マップを取得してSVGを生成
python3 galaxy_transformer.py --download-gaia -f galaxyconfig.properties

# 2回目以降は保存済みのGaiaデータで再生成（加工条件の調整など）
python3 galaxy_transformer.py -f galaxyconfig.properties

# 天の川専用原板を印刷用PDFで生成
python3 galaxy_transformer.py -PDF -f galaxyconfig.properties
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

`galaxyconfig.properties` の `plate.frame-size` で、正方形原盤の一辺を出力上のmmで指定できます。恒星用の `starconfig.properties` とは独立した設定です。例えば `100` なら100 × 100 mmの枠になります。`0` または省略時は用紙サイズと配置数から自動計算し、標準配置では138.5 × 138.5 mmです。`scale` は枠サイズには掛かりません。SVG / PDFおよびlegacy / Gaiaの両モードに適用され、Gaiaの穴配置も指定した枠内に制限されます。

`GalaxyTransformer` は星表のICRS赤道座標を銀河座標へ変換し、**銀緯 −20度以上 +20度以下**の星を次の4枚に割り当てます。銀経は0度以上360度未満に正規化し、共有する境界上の星が重複しないよう、各区間の下端を含み上端を含めません。

| 原盤番号 | 銀経の範囲 | 中心銀経 |
|---|---|---|
| N0 | 0度以上60度未満 | 30度 |
| N1 | 60度以上120度未満 | 90度 |
| S0 | 300度以上360度未満 | 330度 |
| S1 | 240度以上300度未満 | 270度 |

銀経120度以上240度未満、および銀緯±20度の外は出力対象外です。N/Sは用紙の上段/下段の識別子で、赤緯・銀緯の正負を意味しません。`plate.column = 1`、`plate.row = 1` の標準設定では、`galaxy-0.svg` にN0/S0、`galaxy-1.svg` にN1/S1を配置し、計4枚の原盤を2用紙に出力します。PDFも同じ配置です。星がない領域も枠を出力します。

各原盤では中心銀経・銀緯0度への光線を原盤中央に合わせ、円筒投影で位置を求めます。`projector-horizontal` と `projector-vertical` は従来どおり赤道座標での距離で、投影機位置は `(0, ±projector-horizontal, ±projector-vertical)` です。横成分の符号は原盤番号0で正、1で負、縦成分はNで正、Sで負です。星座名は領域内だけを出力し、星座線は銀経・銀緯で線形補間して各領域の境界で切り分けます。

座標変換の回転行列はHipparcosのICRS定義に基づきます（[ERFAの参照実装](https://github.com/liberfa/erfa/blob/master/src/icrs2g.c)）。

### Gaia-DR3からぎんとう原盤を作成

付属の `galaxyconfig.properties` は `galaxy.mode = gaia-etch` を使用します。
従来のHipparcos・Tycho・RC3による出力は `galaxy.mode = legacy` で利用できます。
モード指定のない古い設定ファイルと対話モードは従来方式です。

Gaiaモードは **Gaia DR3のG帯の観測光量だけ**を使います。Hipparcos・Tycho・RC3を重ねず、
V等級への変換や星間減光の除去もしません。暗黒帯の効果を残した、地球から見える恒星光を扱います。
Gaiaの混雑領域での欠測や限界等級より暗い恒星、散光星雲の光は補完しません。
これはG帯に基づく投影原稿で、肉眼の暗所視や写真の完全な再現ではありません。

初回の `--download-gaia` はESA Gaia Archiveに非同期ADQLクエリを送り、4領域それぞれの
`7.5 < G <= 17` の天体を銀経・銀緯0.1度のセルへ集計します。星をランダムに間引かず、
各セルで `SUM(10**(-0.4*(G - reference_magnitude)))` と天体数を取得します。
取得量を抑えるため個々の星をダウンロードせず、セル中心に合計光量を置いて再投影します。
`--gaia-query` で使用する4本のクエリを通信せず確認できます。
取得には時間がかかる場合があり、ジョブURLと状態を表示します。結果の上限超過（OVERFLOW）は拒否します。
途中で失敗した場合、同じコマンドを再実行すると、検証済みの領域ファイルは再利用します。

`galaxy.gaia.input` の既定値は設定ファイルからの相対パス `data/gaia` です。
クエリ・等級・領域・集計幅・チェックサムを保存し、通常実行時は通信せず検証して読み込みます。
Gaiaの等級範囲・取得ビン幅を変えた場合は `--download-gaia` で再取得してください。
個別星のCSV/CSV.gzをこの設定に指定することもできます。列は
`source_id,ra,dec,phot_g_mean_mag` が必須（赤経・赤緯はICRS、度）で、追加列は無視します。
不正値・重複IDはエラーにします。個別CSVは網羅性を確認できないため、その旨をレポートに記録します。
個別CSVの重複チェックはメモリを使用するため、大規模データでは集計キャッシュを推奨します。

| 設定 | 初期値 | 意味 |
|---|---:|---|
| `galaxy.etch.min-hole-diameter-mm` | 0.02 | 最終原盤での最小穴直径。すべての穴にこの径を使用 |
| `galaxy.etch.min-web-mm` | 0.02 | 穴の縁と縁の間に残す金属の最小幅 |
| `galaxy.etch.cell-size-mm` | 0.2 | 原盤で光量を合計するセル幅。穴径＋金属幅の整数倍に切り下げ |
| `galaxy.etch.flux-gain` | 0.05 | 全領域共通の光量倍率（穴面積の合計に掛ける）。微光星増加による飽和を抑える初期値 |
| `galaxy.etch.seed` | 0 | セル内の穴の選択を再現する乱数種 |
| `galaxy.gaia.bright-limit` | 7.5 | G等級の最輝側境界。この等級自体は含まない |
| `galaxy.gaia.faint-limit` | 17 | G等級の最微側境界。この等級自体を含む |
| `galaxy.gaia.sky-bin-deg` | 0.1 | 取得時の天球上の集計幅。0.05〜1度、60度を等分する値 |
| `galaxy.gaia.reference-magnitude` | 7.5 | `star-radius-7.5` の穴面積に対応させるG等級 |

加工寸法は未実測の仮値です。寸法は最終金属原盤でのmmで指定し、出力時は `scale` を掛けます。
`star-radius-7.5` は従来どおり**縮小前の製図上の半径**です。`scale` を変えるときはこの半径も
同じ倍率にしてください。SVG/PDFは「用紙に合わせる」を使わず、指定倍率で製版します。

光量から必要な穴面積を計算し、原盤上のセルごとに合計します。穴の候補点は原盤全体で共通の格子上に
置き、ピッチを穴径＋金属幅とするため、セル境界をまたいでも金属幅を確保します。
各セルの必要穴数の小数部分は蛇行順に次のセルへ繰り越し、容量制限後の原盤全体の量子化誤差を
最小穴0.5個分以内に抑えます。ゼロ光量セルに穴を足すことはありません。
配置可能な穴数を超えた光は別の場所へ移さず、飽和損失として記録します。
穴は**局所光量の標本で、明るい星も含め一対一の実際の恒星位置を表すものではありません**。
実効的な細かさは取得ビン・原盤セル・投影光学系で制限されます。

Gaiaモードは星座線を加工用原盤へ追加せず、入力待ちなしで生成します。
`output/galaxy/gaia-etch-0.svg` と `gaia-etch-1.svg`（`-PDF`指定時はPDF）に計4原盤を出力し、
既存の `galaxy-*.svg/pdf` は残します。接頭辞は `galaxy.gaia.file.prefix` で変更できます。
`gaia-etch-report.json` に入力の出典・設定・天体数・穴数・目標面積・出力面積・飽和損失・量子化誤差を記録します。
穴面積と投影光量の比例は均一照明・薄板を仮定しており、板厚、光源、エッチング精度、
投影時のぼけによる影響は別途試作・測定してください。最小寸法を満たすことは加工成功の保証ではありません。

GaiaモードのG等級境界は、恒星投影機のV等級境界と同一ではありません。
恒星投影機と同時に使う際は、色による違いも含めて重複・明るさを確認してください。

出典: [Gaia DR3](https://www.cosmos.esa.int/web/gaia/dr3)、
[ESA Archiveのプログラムアクセス](https://www.cosmos.esa.int/web/gaia-users/archive/programmatic-access)。
GaiaデータはESA/Gaia/DPACによります。発表等では[Gaiaの謝辞・引用案内](https://www.cosmos.esa.int/web/gaia-users/credits)に従ってください。

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
