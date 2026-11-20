#!/usr/bin/env python3
"""
Self-test пакета lab-interpreter.

Проверяет:
  1. Целостность справочника: у каждого аналита валидные поля, диапазоны,
     множители единиц приводят canonical к canonical (×1), нет дублей синонимов
     между разными аналитами, critical не противоречит ranges.
  2. Поведение движка на эталонных кейсах: статусы, конверсии единиц,
     производные формулы, guard на ошибки единиц, качественные значения.

Запуск:  python3 tests/selftest.py
Код возврата 0 — все тесты прошли, 1 — есть провалы.
"""
import os, sys, io, json

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "engine"))
import interpret as E  # noqa: E402

REFS = os.path.join(ROOT, "references", "reference_ranges.json")

fails = []
def check(cond, msg):
    if not cond:
        fails.append(msg)

# ───────────────────── 1. целостность справочника ─────────────────────
def validate_refs():
    data = json.load(io.open(REFS, encoding="utf-8"))
    analytes = data["analytes"]
    syn_owner = {}
    for key, a in analytes.items():
        check(isinstance(a.get("names", []), list) and a["names"], f"{key}: пустой names")
        check("panel" in a, f"{key}: нет panel")
        if a.get("qualitative"):
            continue  # у качественных нет единиц/диапазонов
        cu = a.get("canonical_unit")
        check(cu is not None, f"{key}: нет canonical_unit")
        acc = a.get("units_accepted", {})
        # canonical-единица, если присутствует в accepted, должна иметь множитель 1
        ncu = E.norm_unit(cu)
        for u, mult in acc.items():
            if E.norm_unit(u) == ncu:
                check(mult == 1, f"{key}: множитель canonical-единицы '{u}' = {mult}, ожидался 1")
            check(isinstance(mult, (int, float)) and mult > 0, f"{key}: множитель '{u}' не положителен")
        # диапазоны: low <= high
        for grp, r in a.get("ranges", {}).items():
            lo, hi = r.get("low"), r.get("high")
            if lo is not None and hi is not None:
                check(lo <= hi, f"{key}.{grp}: low {lo} > high {hi}")
        # critical не внутри нормы (crit_low <= range low, crit_high >= range high)
        crit = a.get("critical") or {}
        dflt = a.get("ranges", {}).get("default", {})
        if crit.get("high") is not None and dflt.get("high") is not None:
            check(crit["high"] >= dflt["high"], f"{key}: critical.high {crit['high']} < range.high {dflt['high']}")
        if crit.get("low") is not None and dflt.get("low") is not None:
            check(crit["low"] <= dflt["low"], f"{key}: critical.low {crit['low']} > range.low {dflt['low']}")
        # синонимы-дубли между разными аналитами
        for n in a["names"]:
            nn = E.norm_name(n)
            if nn in syn_owner and syn_owner[nn] != key:
                fails.append(f"дубль синонима '{n}': {syn_owner[nn]} и {key}")
            syn_owner.setdefault(nn, key)
    return analytes

# ───────────────────── 2. поведение движка ─────────────────────
def run(text, **ctx):
    rb = E.RefBook(REFS)
    c, rows = E.parse_any(text)
    c.update(ctx)
    return E.evaluate(rb, c, rows)

def status_of(rep, key):
    for r in rep["results"]:
        if r["key"] == key:
            return r["status"], r["value"], r["unit"]
    return None, None, None

def derived_val(rep, name_part):
    for x in rep["derived"]:
        if name_part.lower() in x["name"].lower():
            return x["value"]
    return None

def behaviour_tests():
    # конверсия единиц: Hb g/dl → g/L
    r = run("Hb 13.5 g/dl", sex="male")
    st, val, _ = status_of(r, "hemoglobin")
    check(val == 135.0, f"Hb g/dl→g/L: ожидал 135, получил {val}")

    # критическое значение
    r = run("Калий 6.3 ммоль/л")
    st, _, _ = status_of(r, "potassium")
    check(st == "CRITICAL_HIGH", f"калий 6.3: ожидал CRITICAL_HIGH, получил {st}")

    # guard на ошибку единицы (Hb 118 g/dl → 1180)
    r = run("Hb 118 g/dl", sex="male")
    hb = [x for x in r["results"] if x["key"] == "hemoglobin"][0]
    check(hb.get("implausible") is True, "guard не сработал на Hb 118 g/dl")

    # возрастной порог: ЩФ 300 норма у ребёнка, выше у взрослого
    child = run("ЩФ 300 ед/л", sex="male", age=5)
    adult = run("ЩФ 300 ед/л", sex="male", age=40)
    check(status_of(child, "alp")[0] == "NORMAL", "ЩФ 300 у ребёнка должна быть нормой")
    check(status_of(adult, "alp")[0] == "HIGH", "ЩФ 300 у взрослого должна быть выше нормы")

    # качественные: реактивный → POSITIVE, нереактивный → NORMAL
    r = run("HBsAg реактивный\nAnti-HCV нереактивный")
    check(status_of(r, "hbsag")[0] == "POSITIVE", "HBsAg реактивный → POSITIVE")
    check(status_of(r, "anti_hcv")[0] == "NORMAL", "Anti-HCV нереактивный → NORMAL")

    # производные формулы
    r = run("Креатинин 140 мкмоль/л", sex="male", age=50)
    check(derived_val(r, "CKD-EPI") is not None, "CKD-EPI не рассчитан")
    r = run("Креатинин 45 мкмоль/л", sex="male", age=8, height=128)
    check(derived_val(r, "Шварц") is not None, "Шварц не рассчитан у ребёнка")
    check(derived_val(r, "CKD-EPI") is None, "CKD-EPI не должен считаться у ребёнка")

    r = run("Инсулин 18 мкЕд/мл\nГлюкоза 5.8 ммоль/л")
    check(derived_val(r, "HOMA") == 4.64, f"HOMA-IR: ожидал 4.64, получил {derived_val(r,'HOMA')}")

    # ГТТ-статус собирает несколько критериев
    r = run("Глюкоза 6.5 ммоль/л\nГлюкоза через 2 часа 9.2 ммоль/л\nHbA1c 6.1 %")
    gs = derived_val(r, "Гликемический")
    check(gs and gs.count(";") >= 2, f"гликемический статус: ожидал ≥3 критерия, получил «{gs}»")

    # парсер: диапазон в скобках и 10^9 не путаются со значением
    r = run("Гемоглобин 118 г/л (130-170)\nЛейкоциты 14.2 10^9/л", sex="male")
    hb = [x for x in r["results"] if x["key"] == "hemoglobin"][0]
    wbc = [x for x in r["results"] if x["key"] == "wbc"][0]
    check(hb["value"] == 118.0 and hb["ref_low"] == 130.0, f"диапазон в скобках: {hb['value']} [{hb['ref_low']}]")
    check(wbc["value"] == 14.2, f"10^9 парсинг: ожидал 14.2, получил {wbc['value']}")

    # CK-MB/CK индекс
    r = run("КФК 1000 ед/л\nКФК-МВ 80 ед/л", sex="male")
    check(derived_val(r, "CK-MB/CK") == 8.0, f"индекс CK-MB/CK: ожидал 8.0, получил {derived_val(r,'CK-MB/CK')}")

    # delta-check: значимый сдвиг показывается, шум — нет
    hist = [
        {"date": "2024-01-01", "key": "creatinine", "value": 80.0, "unit": "umol/L", "status": "NORMAL"},
        {"date": "2024-01-01", "key": "potassium", "value": 4.0, "unit": "mmol/L", "status": "NORMAL"},
    ]
    cur = [
        {"key": "creatinine", "name_ru": "Креатинин", "value": 125.0, "unit": "umol/L", "status": "HIGH", "unit_known": True},
        {"key": "potassium", "name_ru": "Калий", "value": 4.1, "unit": "mmol/L", "status": "NORMAL", "unit_known": True},
    ]
    dyn = E.compute_dynamics(cur, hist, "2024-02-01")
    keys = {d["key"] for d in dyn}
    check("creatinine" in keys, "delta-check: значимый рост креатинина должен попасть в динамику")
    check("potassium" not in keys, "delta-check: малый сдвиг калия не должен считаться значимым")

    # польские названия и качественные формы
    r = run("Glukoza 88 mg/dl\nŻelazo 111 ug/dl\nBiałko całkowite 74 g/l\nKwas moczowy 4,0 mg/dl")
    for k in ("glucose", "iron", "total_protein", "uric_acid"):
        check(any(x["key"] == k for x in r["results"]), f"польский синоним не распознан: {k}")
    r = run("HIV Ag/Ab NIEREAKTYWNY\nAntygen HBs NIEREAKTYWNY")
    check(status_of(r, "hiv_agab")[0] == "NORMAL", "польский 'niereaktywny' → NORMAL")

    # OCR-искажения единиц: Ul→U/l, ulU/ml→µIU/ml, g/l с мусором, min/ul→10^6/ul
    r = run("AST 19 Ul\nTSH 2,1 ulU/ml\nErytrocyty 5,15 min/ul", sex="male")
    check(status_of(r, "ast")[1] == 19.0 and [x for x in r["results"] if x["key"] == "ast"][0]["unit_known"],
          "OCR-единица 'Ul' у АСТ не распозналась")
    check([x for x in r["results"] if x["key"] == "tsh"][0]["unit_known"], "OCR 'ulU/ml' у ТТГ не распозналась")

    # голый референс без скобок + флаг H/L не должны стать значением
    r = run("Hemoglobina 15,0 g/dl 13,7 - 17,5", sex="male")
    hb = [x for x in r["results"] if x["key"] == "hemoglobin"][0]
    check(hb["value"] == 150.0, f"голый референс: значение должно быть 150, получено {hb['value']}")
    check(hb["status"] == "NORMAL", f"Hb 15 г/дл при референсе из бланка должен быть норма, а не {hb['status']}")

    # качественная серология по числовому порогу assay
    r = run("Herpes simplex virus IgM 0,01")
    check(status_of(r, "hsv_igm")[0] == "NORMAL", "HSV IgM 0.01 ниже порога 2.0 → норма, не POSITIVE")

    # boilerplate-строки не становятся результатами
    r = run("Glukoza 88 mg/dl\ncukrzyca: > 126 mg/dl\nnieprawidłowa glikemia na czczo: 100-125 mg/dl")
    glu = [x for x in r["results"] if x["key"] == "glucose"]
    check(len(glu) == 1, f"boilerplate-строки глюкозы не отфильтрованы: найдено {len(glu)} вместо 1")

    # №6: сверка с флагом лаборатории — расхождение помечается
    r = run("Limfocyty 38 % 20 - 40 H")
    lymp = [x for x in r["results"] if x["key"] == "lymphocytes_pct"][0]
    check(lymp.get("flag_mismatch") is True, "расхождение с флагом лаборатории (H vs норма) не помечено")

    # №5: пограничная зона серологии
    r = run("Chlamydia trachomatis IgA 0,9")
    check(status_of(r, "chlamydia_iga")[0] == "BORDERLINE", "хламидия IgA 0.9 должна быть пограничной")
    r = run("Chlamydia trachomatis IgA 1,5")
    check(status_of(r, "chlamydia_iga")[0] == "POSITIVE", "хламидия IgA 1.5 должна быть положительной")

    # сверка %/абс: процент вне нормы при нормальном абсолюте → пометка относительного сдвига
    r = run("Eozynofile 0,35 tys/ul 0,05 - 0,50\nEozynofile 6,9 % 0,0 - 5,0 H", sex="male", age=26)
    eos_pct = [x for x in r["results"] if x["key"] == "eosinophils_pct"]
    check(eos_pct and eos_pct[0].get("relative_shift_only") is True,
          "эозинофилы: процент ▲ при норм. абсолюте должен помечаться относительным сдвигом")

    # №7: тренд строится по истории
    hist7 = [{"date": "2024-01-01", "key": "tsh", "value": 5.0, "unit": "mIU/L", "status": "HIGH"}]
    cur7 = [{"key": "tsh", "name_ru": "ТТГ", "value": 2.1, "unit": "mIU/L", "status": "NORMAL", "unit_known": True}]
    tr = E.build_trend(cur7, hist7, "2024-03-01", {"tsh": "ТТГ"})
    check(any(row["key"] == "tsh" and len(row["values"]) == 2 for row in tr["rows"]),
          "тренд ТТГ по двум точкам не построен")

def main():
    analytes = validate_refs()
    behaviour_tests()
    total = len(analytes)
    if fails:
        print(f"ПРОВАЛЕНО ({len(fails)}):")
        for f in fails:
            print("  ✗", f)
        print(f"\nсправочник: {total} аналитов")
        return 1
    print(f"ВСЕ ТЕСТЫ ПРОШЛИ ✓   справочник: {total} аналитов, "
          f"{len(set(a['panel'] for a in analytes.values()))} панелей")
    return 0

if __name__ == "__main__":
    sys.exit(main())
