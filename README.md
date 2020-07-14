# pathway_vari

Pythonで動く、変異遺伝子解析ツールです。変異がどの代謝経路に多く蓄積しているかを可視化することができます。代謝経路はKEGGpathway(url:https://www.genome.jp/kegg/pathway.html) を使用しています。多世代に対応しています。
rekgはテスト用のgeneIDが入っています。

## 動かす前に

ツールを動かす前にこちらのライブラリをインストールして下さい。
インストールは以下の方法で行ってください。

```
pip install -U git+https://github.com/CompBio-TDU-Japan/nsdm
```

## 入力

pathway_vari.pyは以下のコマンドで動きます。

```
python payhway_vari.py [geneIDfile1]  [geneIDfile2] ・・・ [gfffile] 
```

[geneIDfile] は1世代分のgeneIDを一つにまとめたファイルです。左から世代順に入力してください。

[gfffile] は対象生物のgffファイルです

### オプション

コマンドに-pまたは-ppiを加えるとタンパク質間相互作用のデータを統合した代謝経路で解析を行います。

コマンドは以下のとおりです。

```
python payhway_vari.py [geneIDfile1] [geneIDfile2] ・・・ [gfffile] -p [ppifile]
```

ppifileはオプションのi後に付けてください。
ppifileは STRINGのデータベース(url:https://string-db.org/cgi/download.pl?sessionId=l3jvlSNOTnAs) から対象の生物のfullバージョンを使用してください。

ppiはコマンド上でevidence typeの指定ができます。
-pまたは--ppiのみで行うと全てのevidence_typeの情報を付与します。
-pの後にevidence typeの名称か下記の番号を指定して下さい。

```
	  evidence_type
		0.  neighborhood
		1.  neighborhood_transferred
		2.  fusion
		3.  cooccurence
		4.  homology
		5.  coexpression
		6.  coexpression_transferred
		7.  experiments
		8.  experiments_transferred
		9.  database
		10. database_transferred
		11. textmining
		12. textmining_transferred
```


evidence typeはSTRINGに準じています。 http://version10.string-db.org/help/getting_started/

例えば、 

```
python pathway_vari.py [geneIDfile1] [gfffile] -p [ppifile] 1  database
```

は```neighborhood_transferred```と```database```のevidence typeを持っているppiの情報をpathwayに付与します。

## 出力

出力されるファイルは以下のとおりです

*  q-value.html  
*  q-value.tsv
*  pathway_images
*  gene_result


### q-value.html
各pathwayのQ valueを世代ごとに表示したグラフ
qvalueはstorey[1]の方法を使用しています。
テストデータで行うと以下の画像のようになります。
【画像】

[1]John D. Storey .et al:Statistical significance for genomewide studies. Proc Natl Acad Sci U S A. 100: 9440–9445,2003

### q-value.tsv
q-value.htmlのQvalueをテキストにまとめたもの

### pathway_images
KEGGpathwayのノードに変異の蓄積量を色付けし、可視化した図が保存されています。

変異の蓄積量は赤の濃さで表しており、色が濃い程変異蓄積量が高いことを示しています。
変異が全くない場合、ノードは青くなります。
白いノードは遺伝子がアノテーションされていないノードです。

### rekg_result
代謝経路に存在する遺伝子と、その内の変異遺伝子が書き込まれています。

### ターミナル上
以下のものが表示されます。
*  pathway number
*  pathway name
*  number of variant proteins
*  The number of proteins in each pathway
*  p-value(Qvalueではない)
