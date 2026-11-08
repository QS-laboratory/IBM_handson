"""Participant wording when importing the immutable first-edition notebook."""
REPLACEMENTS = {'## 振り返りと講師用メモ': '## 振り返り', '**参加者への問い**': '**考えてみましょう**', '**講師向け回答の要点**：1. 積状態なので各切断のランクは1。2. 1つの入力での一致は全入力での一致を保証しない。\n3. 誤ったユニタリも存在する。4. 消えるのは統計変動で、系統的な近似誤差は残る。\n': '', '**学生・企業人向け／目安60–90分＋発展20分**  \n': ''}

def clean_legacy_markdown(source):
    for before, after in REPLACEMENTS.items():
        source = source.replace(before, after)
    return source.strip()
