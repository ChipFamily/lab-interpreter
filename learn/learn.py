#!/usr/bin/env python3
"""
learn.py — обучение lab-interpreter на нераспознанных показателях (офлайн).

Читает learned/unmatched.log (его пишет движок с флагом --learn), группирует
нераспознанные имена по частоте и fuzzy-сопоставляет их с известными
показателями. Подтверждённые синонимы складывает в references/user_overrides.json
(базовый справочник не трогается). Перед сохранением прогоняет self-test — если
обучение что-то сломало, изменения откатываются.

Команды:
  python3 learn/learn.py review              # показать нераспознанное + подсказки
  python3 learn/learn.py apply KEY=НАЗВАНИЕ  # добавить синоним(ы) вручную
  python3 learn/learn.py auto                # автодобавить уверенные совпадения (≥0.88)
  python3 learn/learn.py auto --threshold 0.9
"""
import os, sys, io, json, argparse, difflib, shutil, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "engine"))
import interpret as E  # noqa: E402

REFS = os.path.join(ROOT, "references", "reference_ranges.json")
OVERLAY = os.path.join(ROOT, "references", "user_overrides.json")
LOG = os.path.join(ROOT, "learned", "unmatched.log")
SELFTEST = os.path.join(ROOT, "tests", "selftest.py")


def load_log() -> dict:
    """Нераспознанные имена → {count, examples:[(value,unit)...]}."""
    groups: dict[str, dict] = {}
    if not os.path.exists(LOG):
        return groups
    for line in io.open(LOG, encoding="utf-8"):
        line = line.strip()
        if not line:
            continue
        try:
            rec = json.loads(line)
        except Exception:
            continue
        name = (rec.get("name") or "").strip()
        if not name:
            continue
        g = groups.setdefault(name, {"count": 0, "examples": []})
        g["count"] += 1
        if len(g["examples"]) < 3:
            g["examples"].append((rec.get("value"), rec.get("unit")))
    return groups


def load_overlay() -> dict:
    if os.path.exists(OVERLAY):
        try:
            return json.load(io.open(OVERLAY, encoding="utf-8"))
        except Exception:
            pass
    return {"add_synonyms": {}, "analytes": {}}


def save_overlay(ov: dict):
    os.makedirs(os.path.dirname(OVERLAY), exist_ok=True)
    json.dump(ov, io.open(OVERLAY, "w", encoding="utf-8"), ensure_ascii=False, indent=2)


def run_selftest() -> bool:
    r = subprocess.run([sys.executable, SELFTEST], capture_output=True, text=True)
    if r.returncode != 0:
        sys.stderr.write(r.stdout + r.stderr + "\n")
    return r.returncode == 0


def suggest_matches(name: str, rb: E.RefBook, n=3):
    """Fuzzy-кандидаты (ключ аналита, его отображаемое имя, оценка сходства)."""
    nn = E.norm_name(name)
    keys = list(rb.index.keys())
    close = difflib.get_close_matches(nn, keys, n=n, cutoff=0.6)
    out = []
    for syn in close:
        key = rb.index[syn]
        ratio = difflib.SequenceMatcher(None, nn, syn).ratio()
        ru = rb.analytes[key].get("ru", key)
        out.append((key, ru, round(ratio, 2)))
    return out


def cmd_review(args):
    rb = E.RefBook(REFS)
    groups = load_log()
    if not groups:
        print("Лог нераспознанного пуст (learned/unmatched.log). Нечего разбирать.")
        return 0
    print(f"Нераспознано разных показателей: {len(groups)}\n")
    for name, g in sorted(groups.items(), key=lambda kv: -kv[1]["count"]):
        ex = ", ".join(f"{v} {u}".strip() for v, u in g["examples"])
        print(f"• «{name}»  ×{g['count']}   примеры: {ex}")
        for key, ru, ratio in suggest_matches(name, rb):
            mark = "↳ скорее всего" if ratio >= 0.88 else "↳ возможно"
            print(f"      {mark} «{ru}»  (ключ: {key}, сходство {ratio})  "
                  f"→ apply {key}=\"{name}\"")
        print()
    print("Добавить синоним:  python3 learn/learn.py apply <ключ>=\"<название>\"")
    print("Автодобавление:    python3 learn/learn.py auto")
    return 0


def _add_synonyms(pairs: list[tuple[str, str]]) -> int:
    """pairs: [(key, name)]. Возвращает число реально добавленных."""
    rb = E.RefBook(REFS)  # базовый, без overlay — чтобы знать валидные ключи
    ov = load_overlay()
    ov.setdefault("add_synonyms", {})
    added = 0
    for key, name in pairs:
        if key not in rb.analytes:
            sys.stderr.write(f"✗ ключ '{key}' не найден в справочнике — пропуск\n")
            continue
        lst = ov["add_synonyms"].setdefault(key, [])
        if name not in lst:
            lst.append(name)
            added += 1
            print(f"+ «{name}» → {key} ({rb.analytes[key].get('ru', key)})")
    if not added:
        return 0
    # self-test как ворота: пишем во временный backup, проверяем, откатываем при провале
    backup = OVERLAY + ".bak" if os.path.exists(OVERLAY) else None
    if backup:
        shutil.copy(OVERLAY, backup)
    save_overlay(ov)
    if run_selftest():
        print(f"\n✓ self-test пройден, добавлено синонимов: {added}")
        if backup and os.path.exists(backup):
            os.remove(backup)
        return added
    else:
        print("\n✗ self-test ПРОВАЛЕН — изменения откачены")
        if backup:
            shutil.move(backup, OVERLAY)
        else:
            os.remove(OVERLAY)
        return 0


def cmd_apply(args):
    pairs = []
    for item in args.pairs:
        if "=" not in item:
            sys.stderr.write(f"✗ неверный формат '{item}', нужно KEY=НАЗВАНИЕ\n")
            continue
        key, name = item.split("=", 1)
        pairs.append((key.strip(), name.strip().strip('"')))
    if not pairs:
        return 1
    return 0 if _add_synonyms(pairs) >= 0 else 1


def cmd_auto(args):
    rb = E.RefBook(REFS)
    groups = load_log()
    pairs = []
    skipped = []
    for name in groups:
        sug = suggest_matches(name, rb, n=1)
        if sug and sug[0][2] >= args.threshold:
            pairs.append((sug[0][0], name))
        else:
            skipped.append(name)
    if not pairs:
        print(f"Нет совпадений с сходством ≥{args.threshold}. "
              f"Разбери вручную: python3 learn/learn.py review")
        return 0
    print(f"Автодобавление {len(pairs)} синонимов (порог {args.threshold}):")
    _add_synonyms(pairs)
    if skipped:
        print(f"\nОсталось на ручной разбор ({len(skipped)}): "
              + ", ".join(f'«{s}»' for s in skipped[:10]) + ("…" if len(skipped) > 10 else ""))
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description="Обучение lab-interpreter на нераспознанном (офлайн).")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("review", help="показать нераспознанное и подсказки")
    pa = sub.add_parser("apply", help="добавить синоним(ы): KEY=НАЗВАНИЕ")
    pa.add_argument("pairs", nargs="+")
    au = sub.add_parser("auto", help="автодобавить уверенные совпадения")
    au.add_argument("--threshold", type=float, default=0.88)
    args = ap.parse_args(argv)
    return {"review": cmd_review, "apply": cmd_apply, "auto": cmd_auto}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
