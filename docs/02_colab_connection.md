# ② IBM Quantum実習用アカウントとGoogle Colabの接続確認

先に① 利用登録を完了し、実習用インスタンスが見えることを確認してください。

今回は認証とバックエンド一覧の取得までを確認します。
Googleアカウントと、Google Colabを利用できる環境をご用意ください。

## 1. 練習用のAPIキーを作成する

1. ご自身のIBMidでIBM Quantum Platform（[https://quantum.cloud.ibm.com/](https://quantum.cloud.ibm.com/)）にログインします。
2. 今回招待された東京大学のアカウントを選びます。
3. ホーム画面の **API key → Create** からキーを作成します。
4. 名前を求められたら、例えば `FTQC-HandsOn-test` とします。
5. 作成されたキーをその場でコピーし、安全に保管します。同じキーを後から再表示することはできません。

**APIキーは本人の接続に使う秘密情報です。主催者への送付は不要です。**
メール、チャット、共有ノートブック、GitHubへ貼り付けないでください。管理者のキーや他人のキーは使いません。

## 2. 実習用インスタンスのCRNをコピーする

1. 画面上部のリージョンで **us-east** を選び、**Instances** を開きます。
2. **FTQC-HandsOn-20261109** のCRN欄のコピーアイコンをクリックします。

CRNは `crn:` で始まる長い識別子です。今回使う実習用インスタンスを指定するもので、APIキーとは異なります。

## 3. Google Colabを準備する

Google Colab（[https://colab.research.google.com/](https://colab.research.google.com/)）を開き、新しいノートブックを作成します。
最初のコードセルに次を貼り付け、実行ボタンを押してください。

```python
%pip -q install -U qiskit-ibm-runtime
```

インストールが完了してから次へ進みます。Colabからランタイム再起動を求められた場合は、再起動してから次のセルを実行してください。

## 4. 接続し、利用可能なバックエンドを表示する

新しいコードセルに、次をそのまま貼り付けて実行します。

```python
from getpass import getpass
from qiskit_ibm_runtime import QiskitRuntimeService

api_key = getpass("本人のAPIキーを貼り付け、Enter：").strip()
instance_crn = input("実習用インスタンスのCRNを貼り付け、Enter：").strip()

service = QiskitRuntimeService(
    channel="ibm_quantum_platform",
    token=api_key,
    instance=instance_crn,
)

del api_key

backends = service.backends()
print(f"利用可能なバックエンド数：{len(backends)}")
for backend in backends:
    print(backend.name)

if backends:
    print("接続確認完了です。")
else:
    print("一覧が空です。主催者へお知らせください。")
```

入力欄が順番に出ます。手順1のAPIキー、手順2のCRNを入力し、それぞれEnterを押します。
APIキーの入力内容は画面に表示されません。コード本文を書き換えてキーを埋め込む必要はありません。

バックエンド名（`ibm_...` など）が表示されれば、本人のAPIキーで指定インスタンスへ接続できています。

## 5. 確認後の操作

1. Colabの **ランタイム** メニューから、実行環境を切断・削除します。
2. IBM Quantum PlatformのAPIキー管理画面（**View all** など）で、今回作った練習用キーを削除します。
3. 当日は新しいAPIキーを作成して使用します。

削除するのは **練習用のAPIキーだけ** です。IBMid、アカウントへの参加登録、実習用インスタンス、CRNはそのままで構いません。

## うまくいかない場合

| 状況 | 確認すること |
|---|---|
| 認証エラー・401 | キーのコピー漏れ、削除済みキー、別の文字列を貼り付けていないか |
| 権限エラー・403 | 招待を承諾した本人のキーか、実習用Collaboratorsグループに登録されているか |
| インスタンスが見つからない | 今回の実習用CRNを正しくコピーしたか |
| 一覧が空 | 主催者へ連絡し、インスタンスのQPU割当を確認する |
| 通信エラー | 社内ネットワークでColabやIBM Quantumへの接続が制限されていないか |

## IBM公式資料

- APIキー・CRNの取得と接続（[https://quantum.cloud.ibm.com/docs/en/guides/initialize-account](https://quantum.cloud.ibm.com/docs/en/guides/initialize-account)）
- 外部環境での接続とAPIキーの扱い（[https://quantum.cloud.ibm.com/docs/en/guides/cloud-setup-untrusted](https://quantum.cloud.ibm.com/docs/en/guides/cloud-setup-untrusted)）

最終更新：2026年11月6日
