#!/usr/bin/env python3
"""
update.py — безопасное обновление базового справочника lab-interpreter.

Берёт новую версию reference_ranges.json (из файла, по URL или git pull),
прогоняет на ней self-test и применяет ТОЛЬКО если тесты зелёные. Иначе
откатывает. Старый справочник сохраняется в references/*.bak.

Это единственное место, где возможна сеть, и запускается оно вручную/по задаче.

Команды:
  python3 update/update.py --from-file НОВЫЙ.json
  python3 update/update.py --from-url https://.../reference_ranges.json
  python3 update/update.py --git            # git pull в каталоге пакета
  python3 update/update.py --check          # только прогнать self-test на текущем
"""
import os, sys, json, argparse, shutil, subprocess, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
REFS = os.path.join(ROOT, "references", "reference_ranges.json")
SELFTEST = os.path.join(ROOT, "tests", "selftest.py")


def selftest_ok() -> bool:
    r = subprocess.run([sys.executable, SELFTEST], capture_output=True, text=True)
    sys.stdout.write(r.stdout)
    if r.returncode != 0:
        sys.stderr.write(r.stderr)
    return r.returncode == 0


def validate_json(path: str) -> bool:
    try:
        d = json.load(open(path, encoding="utf-8"))
    except Exception as e:
        sys.stderr.write(f"✗ не JSON: {e}\n")
        return False
    if not isinstance(d.get("analytes"), dict) or not d["analytes"]:
        sys.stderr.write("✗ нет непустого блока 'analytes'\n")
        return False
    return True


def apply_candidate(cand_path: str) -> int:
    """Проверяет cand_path и, если он валиден и self-test зелёный, ставит на место REFS."""
    if not validate_json(cand_path):
        return 1
    # число аналитов не должно внезапно обрушиться (защита от битого обновления)
    old_n = len(json.load(open(REFS, encoding="utf-8"))["analytes"]) if os.path.exists(REFS) else 0
    new_n = len(json.load(open(cand_path, encoding="utf-8"))["analytes"])
    if old_n and new_n < old_n * 0.8:
        sys.stderr.write(f"✗ новый справочник резко меньше ({new_n} vs {old_n}) — отказ, проверьте источник\n")
        return 1
    backup = REFS + ".bak"
    if os.path.exists(REFS):
        shutil.copy(REFS, backup)
    shutil.copy(cand_path, REFS)
    if selftest_ok():
        print(f"\n✓ обновлено: {old_n} → {new_n} аналитов. Прежняя версия в {os.path.basename(backup)}")
        return 0
    else:
        print("\n✗ self-test ПРОВАЛЕН на новом справочнике — откат к прежней версии")
        if os.path.exists(backup):
            shutil.move(backup, REFS)
        return 1


def cmd_from_file(path: str) -> int:
    if not os.path.exists(path):
        sys.stderr.write(f"✗ файл не найден: {path}\n")
        return 1
    return apply_candidate(path)


def cmd_from_url(url: str) -> int:
    import urllib.request
    print(f"Загрузка {url} …")
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            data = resp.read()
    except Exception as e:
        sys.stderr.write(f"✗ не удалось скачать: {e}\n")
        return 1
    with tempfile.NamedTemporaryFile("wb", suffix=".json", delete=False) as tmp:
        tmp.write(data)
        tmp_path = tmp.name
    try:
        return apply_candidate(tmp_path)
    finally:
        os.unlink(tmp_path)


def cmd_git() -> int:
    if not os.path.isdir(os.path.join(ROOT, ".git")):
        sys.stderr.write("✗ это не git-репозиторий — используйте --from-url или --from-file\n")
        return 1
    backup = REFS + ".bak"
    if os.path.exists(REFS):
        shutil.copy(REFS, backup)
    r = subprocess.run(["git", "-C", ROOT, "pull", "--ff-only"], capture_output=True, text=True)
    sys.stdout.write(r.stdout)
    if r.returncode != 0:
        sys.stderr.write(r.stderr)
        return 1
    if selftest_ok():
        print("\n✓ git pull применён, self-test пройден")
        if os.path.exists(backup):
            os.remove(backup)
        return 0
    else:
        print("\n✗ self-test провален после pull — откат справочника")
        if os.path.exists(backup):
            shutil.move(backup, REFS)
        return 1


def main(argv=None):
    ap = argparse.ArgumentParser(description="Безопасное обновление справочника (self-test как ворота).")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--from-file", metavar="PATH")
    g.add_argument("--from-url", metavar="URL")
    g.add_argument("--git", action="store_true")
    g.add_argument("--check", action="store_true", help="только прогнать self-test на текущем справочнике")
    args = ap.parse_args(argv)
    if args.check:
        return 0 if selftest_ok() else 1
    if args.from_file:
        return cmd_from_file(args.from_file)
    if args.from_url:
        return cmd_from_url(args.from_url)
    if args.git:
        return cmd_git()


if __name__ == "__main__":
    sys.exit(main())
