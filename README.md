# Orb Transform Library (OTL) — Python版

星表データから、ドーム型プラネタリウム投影機用の「原板」データ(星の位置を打ち抜く版下)を生成するツールです。東京大学地文研究会天文部が開発した Java 版 OTL (Ver. 1.4) の Python 移植版で、挙動を変えないことを優先した忠実な移植になっています。

## これは何をするものか

自作の投影機式プラネタリウム(正十二面体などにレンズユニットを複数配置して天球全体を投影するタイプ)を作るとき、各投影機ユニットに取り付ける「原板」に星の穴をどこに開けるかを計算し、SVG または印刷用 PDF として出力します。

1. `SphereReader` が Hipparcos/Tycho の星表、RC3 銀河カタログ、IAU88 星座線カタログを読み込み、天球上の星・星座を取得
2. `Transformer`(または銀河専用の `GalaxyTransformer`)が、指定したドーム半径・投影機配置・レンズ焦点距離などの幾何パラメータに基づいて、天球上の位置を各ユニットの原板上の位置(mm単位のXY座標)に変換
3. `PlateWriterSVG` / `PlateWriterPDF` が変換結果を版下データとして書き出す

元の `OTL_prog` に後から追加された次の機能にも対応しています。

- 星座線を構成する恒星だけを、指定した倍率で拡大
- 各原盤が担当する天球領域をポリゴンSVGとして自動生成

## 構成

| ファイル | 役割 |
|---|---|
| `transformer.py` | 星・星座用の原板を生成するメインスクリプト(エントリポイント) |
| `galaxy_transformer.py` | 天の川(銀河)専用の原板を生成するスクリプト(エントリポイント) |
| `basic_transformer.py` | 上記2つの変換処理に共通する抽象基底クラス |
| `sphere_reader.py` | `hip_main.dat` / `tyc_main.dat` / `rc3.dat` / `IAU88.hlc` を読み込むリーダー |
| `unit_arrangement.py` | 正十二面体をもとにした投影機ユニットの配置計算 |
| `geometry.py` | 錐体(GeneralizedCone)・多面体(Polyhedron)などの幾何計算 |
| `mathvector.py` | 3次元ベクトル演算(`MathVector.java` の移植) |
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

いずれも公開の天文カタログで、サイズが大きいため(`tyc_main.dat` は数百MB)Git管理下に置く場合は Git LFS 等を検討してください。

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

主なオプション:

- `-f [configFile]` — プラネタリウムの設定を `.properties` ファイルから読み込む
- `-PDF` — 原板データを印刷用 PDF 形式で出力する(未指定時は SVG)
- `--polygons` — 各原盤の担当星域ポリゴンを別の SVG に出力する
- `-h`, `-help` — 使い方を表示

設定ファイルの各項目(ドーム半径、投影機とレンズの距離、正十二面体上のユニット配置箇所、星の等級の上限・下限など)は `starconfig.properties` / `galaxyconfig.properties` にコメント付きで記載しています。

生成した SVG / PDF とユニット配置画像は、設定ファイルの `output.directory` で指定したフォルダに出力されます。未指定時は `output` フォルダが自動作成されます。

印刷用 PDF は幅 191 mm × 高さ 277 mm の1ページ・ベクトルデータです。原板枠と星穴以外の表示は `N0` / `S0` などの原盤番号だけで、SVG の向きマーク、星座名、星座線は含みません。印刷時は「用紙に合わせる」を無効にし、100%(実寸)で出力してください。旧 `-PS` オプションは廃止され、使用するとエラーになります。

### 星座線構成星の拡大

`star.enlarge-rate` に0以上の値を指定します。`1.0` を指定すると、4等星相当の穴は1等級分大きくなり、明るい星ほど拡大量が滑らかに抑えられます。`0.0` では拡大しません。旧Java版の設定名 `star-EnlargeRate` も使用できます。

対象星は `IAU88.hlc` のHIP番号と、旧版の `hip_constellation_line_star.csv` の補助データを併合した893個です。HIP番号による照合なので、旧実装の原盤座標の丸め比較より安定しています。

### 担当星域ポリゴン

設定ファイルで `polygon.enabled = yes` とするか、`--polygons` を指定すると生成します。`polygon.samples` は1原盤あたりの頂点数（12以上、既定値180）、`polygon.file-prefix` は出力ファイル名の接頭辞です。

ポリゴンは各原盤座標から天球へ光線を逆投影し、通常の星配置と同じユニット選択規則で境界を二分探索して求めます。枠の円でクリップ済みで、星原盤とは別のSVGに同じ用紙配置・黒い原盤・白い投影領域で出力されます。PDFで星を出力する場合も、ポリゴンはSVGです。

## 移植上の注意

- 乱数生成は Python 標準の `random` モジュールを使用しています。Java 版の `java.util.Random` とはアルゴリズムが異なるため、同じシード値(0)でも生成される乱数列そのものは一致しません(統計的性質は同等です)。
- `models.py` の `PlatePosition.__eq__` など、Java 版の挙動をそのまま踏襲した(一見不自然に見える)比較ロジックが一部にあります。詳細はコード中の `NOTE` コメントを参照してください。

## ライセンス

GNU General Public License v2 (or later)。詳細は各ソースファイル冒頭のヘッダーを参照してください。

Copyright (C) 2007, 2012-2014 東京大学地文研究会天文部
