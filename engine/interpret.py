#!/usr/bin/env python3
"""
Lab Interpreter — детерминированное ядро расшифровки анализов.

Принимает бланк в любом из форматов:
  - JSON:  {"sex":"male","age":45,"results":[{"name":"Hemoglobin","value":135,"unit":"g/L"}, ...]}
  - CSV:   name,value,unit[,ref_low,ref_high]
  - текст: свободный «Гемоглобин 135 г/л», «ALT: 52 U/L (0-41)», по строке на аналит

Делает:
  - матч аналита по словарю синонимов (RU/EN, регистронезависимо, с нормализацией)
  - приведение единицы входа к canonical через множители из reference_ranges.json
  - выбор диапазона по полу; приоритет у диапазона, указанного в самом бланке
  - статус: LOW / HIGH / NORMAL / CRITICAL_LOW / CRITICAL_HIGH / UNKNOWN
  - производные показатели и флаги-паттерны (De Ritis, не-ЛПВП, стадия ХБП и т.п.)
  - машинный (JSON) и человекочитаемый отчёт

Запуск:
  python3 interpret.py BLANK_FILE [--sex male|female] [--age N] [--json] [--refs PATH]
  cat blank.txt | python3 interpret.py - --sex female
"""
from __future__ import annotations
import sys, os, re, json, argparse, csv, io, math
from typing import Any, Optional

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(HERE)
DEFAULT_REFS = os.path.join(HERE, "..", "references", "reference_ranges.json")

# ───────────────────────────── нормализация ──────────────────────────────

_UNIT_FIX = {
    "µ": "u", "μ": "u", "×": "x", "·": ".", "−": "-", "–": "-",
}
# Кириллические обозначения единиц → латиница. Порядок важен: длинные подстроки
# раньше коротких, иначе «л» внутри «ммоль» сломает замену.
_CYR_UNIT = [
    ("мкмоль", "umol"), ("ммоль", "mmol"), ("пмоль", "pmol"),
    ("нмоль", "nmol"), ("моль", "mol"),
    ("мкме", "uiu"), ("мме", "miu"), ("ме", "iu"),
    ("мкл", "ul"), ("мкг", "ug"),
    ("тыс", "10^3"), ("млн", "10^6"),
    ("мг", "mg"), ("нг", "ng"), ("пг", "pg"), ("дл", "dl"), ("мл", "ml"),
    ("мм", "mm"), ("ед", "u"), ("мин", "min"),
    ("г", "g"), ("л", "l"), ("ч", "hr"), ("сек", "s"), ("с", "s"), ("х", "x"),
]
# Типичные искажения единиц после OCR (наблюдались на реальных польских сканах):
# pe/ml→pg/ml, ulU/ml/LIU→µIU/ml, 5/1→g/l, Ul→u/l, min/ul→mln/ul, mmol/I→mmol/l, ps/mi→pg/ml
_OCR_UNIT_FIX = [
    ("pe/ml", "pg/ml"), ("pg/mi", "pg/ml"), ("ps/mi", "pg/ml"), ("ps/ml", "pg/ml"),
    ("ulu/ml", "uiu/ml"), ("liu/ml", "uiu/ml"), ("liu", "uiu"), ("ulu", "uiu"),
    ("miu/mi", "miu/ml"), ("iu/mi", "iu/ml"),
    ("min/ul", "10^6/ul"), ("mln/ul", "10^6/ul"),
    ("mmol/i", "mmol/l"), ("umol/i", "umol/l"), ("nmol/i", "nmol/l"), ("mmol/i*", "mmol/l"),
    ("g/i", "g/l"), ("u/i", "u/l"), ("iu/i", "iu/l"), ("mg/i", "mg/l"), ("ng/mi", "ng/ml"),
    ("tys/ul", "10^3/ul"), ("5/1", "g/l"), ("5/i", "g/l"),
    ("ul", "u/l"), ("ui", "u/l"),  # bare «Ul»/«UI» у ферментов = U/l (счётные клетки всегда со слешем)
    ("tu/ml", "iu/ml"), ("tu/l", "iu/l"), ("lu/ml", "iu/ml"),  # OCR: I→T/L у IU
]
def norm_unit(u: str) -> str:
    if u is None:
        return ""
    u = u.strip().lower()
    for a, b in _UNIT_FIX.items():
        u = u.replace(a, b)
    u = u.replace(" ", "")
    u = u.rstrip("*")  # хвостовая звёздочка-пометка из бланка
    u = u.replace("10e", "10^").replace("10*", "10^")
    for a, b in _CYR_UNIT:
        u = u.replace(a, b)
    for a, b in _OCR_UNIT_FIX:
        if u == a:
            u = b
            break
    return u

def norm_name(s: str) -> str:
    if s is None:
        return ""
    s = s.strip().lower()
    for a, b in _UNIT_FIX.items():
        s = s.replace(a, b)
    s = re.sub(r"[\t]+", " ", s)
    s = re.sub(r"\s+", " ", s)
    return s.strip(" :;.- ")

def to_float(v: Any) -> Optional[float]:
    if v is None:
        return None
    if isinstance(v, (int, float)):
        return float(v)
    s = str(v).strip()
    if not s:
        return None
    # отбросить операторы сравнения и единицы-хвосты
    s = s.replace(",", ".")
    m = re.search(r"[-+]?\d*\.?\d+(?:e[-+]?\d+)?", s, re.I)
    return float(m.group(0)) if m else None

# ───────────────────────────── справочник ────────────────────────────────

class RefBook:
    def __init__(self, path: str, overlay: Optional[str] = None):
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        self.analytes: dict[str, dict] = data["analytes"]

        # overlay (user_overrides.json) — всё «выученное»: новые синонимы к
        # существующим аналитам и новые аналиты. Базовый справочник не меняется.
        if overlay is None:
            cand_overlay = os.path.join(os.path.dirname(path), "user_overrides.json")
            overlay = cand_overlay if os.path.exists(cand_overlay) else None
        if overlay and os.path.exists(overlay):
            try:
                with open(overlay, "r", encoding="utf-8") as f:
                    ov = json.load(f)
                # add_synonyms: {existing_key: [names...]} — расширяет names аналита
                for key, syns in (ov.get("add_synonyms") or {}).items():
                    if key in self.analytes:
                        names = self.analytes[key].setdefault("names", [])
                        for s in syns:
                            if s not in names:
                                names.append(s)
                # analytes: новые аналиты (или переопределения) из overlay
                for key, a in (ov.get("analytes") or {}).items():
                    self.analytes[key] = a
            except Exception:
                pass  # overlay best-effort: повреждён — работаем на базовом справочнике

        # обратный индекс: нормализованный синоним → ключ аналита
        self.index: dict[str, str] = {}
        for key, a in self.analytes.items():
            cand = set(a.get("names", [])) | {a.get("ru", ""), a.get("en", ""), key}
            for n in cand:
                nn = norm_name(n)
                if nn:
                    self.index.setdefault(nn, key)

    def match(self, name: str) -> Optional[str]:
        nn = norm_name(name)
        # снять хвостовые сноски-маркеры (одиночные цифры/номера после названия: «... 2»)
        nn = re.sub(r"\s+\d{1,3}$", "", nn).strip()
        nn = re.sub(r"[\"'’`]+", "", nn)
        if nn in self.index:
            return self.index[nn]
        # приоритет многословному совпадению: самый длинный индекс-ключ, целиком входящий
        # в имя как фраза (чтобы «hemoglobina glikowana» → hba1c, а не hemoglobin)
        phrase = [(syn, key) for syn, key in self.index.items()
                  if " " in syn and len(syn) >= 6 and syn in nn]
        if phrase:
            phrase.sort(key=lambda sk: -len(sk[0]))
            return phrase[0][1]
        # смешанные кирилло-латинские имена: «анти-dsDNA» ↔ «anti-dsdna», «ат к …» ↔ «anti-…»
        for a, b in (("анти-", "anti-"), ("анти ", "anti "), ("anti-", "анти-"),
                     ("anti ", "анти "), ("ат к ", "anti-"), ("ат-", "anti-")):
            if nn.startswith(a):
                alt = b + nn[len(a):]
                if alt in self.index:
                    return self.index[alt]
        # мягкий матч: убрать %/абс пометки и единичные скобки, попробовать подстроку
        stripped = re.sub(r"[(\[].*?[)\]]", "", nn).strip()
        stripped = re.sub(r"\b(абс|abs|absolute|count|conc)\b", "", stripped).strip()
        if stripped in self.index:
            return self.index[stripped]
        for syn, key in self.index.items():
            if len(syn) >= 4 and (syn == stripped or syn in nn.split()):
                return key
        # нечёткий матч: мелкие OCR-опечатки в названии (напр. «креатенин»→«креатинин»,
        # « bialko calkowite » с пропущенными диакритиками). Порог высокий, чтобы не путать.
        import difflib
        cand = difflib.get_close_matches(stripped or nn, list(self.index.keys()), n=1, cutoff=0.9)
        if cand:
            return self.index[cand[0]]
        return None

    def convert(self, key: str, value: float, unit: str) -> tuple[Optional[float], str, bool]:
        """→ (значение в canonical, canonical_unit, конверсия_известна)"""
        a = self.analytes[key]
        canon = a["canonical_unit"]
        u = norm_unit(unit)
        if not u:
            return value, canon, True  # единицы нет — считаем, что уже canonical
        accepted = {norm_unit(k): v for k, v in a.get("units_accepted", {}).items()}
        cu = norm_unit(canon)
        if u == cu or u in ("", "ratio") and norm_unit(canon) in ("ratio", ""):
            return value, canon, True
        if u in accepted:
            return value * accepted[u], canon, True
        if u == cu:
            return value, canon, True
        if a.get("unit_optional"):
            return value, canon, True  # счётный показатель — единица несущественна
        return value, canon, False  # неизвестная единица — не домысливаем

# ───────────────────────────── парсеры входа ─────────────────────────────

def parse_json(text: str) -> tuple[dict, list[dict]]:
    obj = json.loads(text)
    ctx = {k: obj.get(k) for k in ("sex", "age", "weight", "height", "pregnant", "on_anticoagulant") if k in obj}
    # height уже включён выше
    rows = obj.get("results") or obj.get("analytes") or obj.get("tests") or []
    out = []
    for r in rows:
        out.append({
            "name": r.get("name") or r.get("analyte") or r.get("test") or "",
            "value": r.get("value", r.get("result")),
            "unit": r.get("unit", r.get("units", "")),
            "ref_low": r.get("ref_low", r.get("low")),
            "ref_high": r.get("ref_high", r.get("high")),
        })
    return ctx, out

def parse_csv(text: str) -> tuple[dict, list[dict]]:
    out = []
    reader = csv.DictReader(io.StringIO(text))
    if reader.fieldnames and any(norm_name(h) in ("name", "test", "analyte", "показатель") for h in reader.fieldnames):
        fmap = {norm_name(h): h for h in reader.fieldnames}
        def g(row, *keys):
            for k in keys:
                if k in fmap:
                    return row[fmap[k]]
            return None
        for row in reader:
            out.append({
                "name": g(row, "name", "test", "analyte", "показатель"),
                "value": g(row, "value", "result", "значение"),
                "unit": g(row, "unit", "units", "единица", "ед"),
                "ref_low": g(row, "ref_low", "low", "мин"),
                "ref_high": g(row, "ref_high", "high", "макс"),
            })
        return {}, out
    # без заголовка: name,value,unit[,low,high]
    for parts in csv.reader(io.StringIO(text)):
        if not parts:
            continue
        rec = {"name": parts[0]}
        if len(parts) > 1: rec["value"] = parts[1]
        if len(parts) > 2: rec["unit"] = parts[2]
        if len(parts) > 3: rec["ref_low"] = parts[3]
        if len(parts) > 4: rec["ref_high"] = parts[4]
        out.append(rec)
    return {}, out

# референс в конце строки: в скобках (low-high) ИЛИ голый «low - high» (как в бланках)
_REF_TAIL = re.compile(
    r"[(\[]\s*(?P<low>[-+]?\d+(?:[.,]\d+)?)\s*[-–—]\s*(?P<high>[-+]?\d+(?:[.,]\d+)?)\s*[)\]]\s*$")
_REF_TAIL_BARE = re.compile(
    r"(?<=\s)(?P<low>\d+(?:[.,]\d+)?)\s*[-–—]\s*(?P<high>\d+(?:[.,]\d+)?)\s*$")
# хвостовой флаг лаборатории (LIW): одиночная буква H/L/N (лат/кир) после результата
_LIW_TAIL = re.compile(r"(?<=\s)[HLNНЛ]\s*$")
_NUM = re.compile(r"[-+]?\d+(?:[.,]\d+)?")

def _split_measurement(line: str):
    """Разбор «<имя> <значение> <единица> [реф.диапазон] [H/L]». Значение — последнее
    число ПОСЛЕ отсечения хвостового референса и флага; не входящее в 10^9."""
    low = high = None
    lab_flag = None
    core = line.rstrip()
    # 1) отсечь хвостовой флаг LIW (H/L/N) — но СОХРАНИТЬ для сверки с вердиктом движка
    mflag = _LIW_TAIL.search(core)
    if mflag:
        ch = mflag.group(0).strip().upper()
        lab_flag = {"H": "H", "L": "L", "N": "N", "Н": "H", "Л": "L"}.get(ch)
        core = core[:mflag.start()].rstrip()
    # 2) отсечь хвостовой референс — сперва в скобках, затем голый «low - high»
    mref = _REF_TAIL.search(core)
    if mref:
        low, high = mref.group("low"), mref.group("high")
        core = core[:mref.start()].rstrip()
    else:
        mb = _REF_TAIL_BARE.search(core)
        # голый диапазон засчитываем только если перед ним есть ещё хотя бы одно число
        # (само значение), иначе это может быть сам результат
        if mb and len(_NUM.findall(core[:mb.start()])) >= 1:
            low, high = mb.group("low"), mb.group("high")
            core = core[:mb.start()].rstrip(" \t-–—")
    cands = [m for m in _NUM.finditer(core)
             if (m.start() == 0 or core[m.start() - 1] in " \t")        # пробел/начало слева
             and core[m.end():m.end() + 1] != "^"                        # не часть 10^9
             and not re.match(r"\s*[-(]?\s*(oh|он)", core[m.end():m.end() + 5], re.I)]  # не «25(OH)»
    if not cands:
        return None
    mv = cands[-1]
    name = core[:mv.start()].strip(" :=\t-–")
    unit = core[mv.end():].strip()
    # вычистить хвостовой мусор у единицы: « -», «—», «:», «*», скобки-пометки
    unit = re.sub(r"[\s:*\-–—]+$", "", unit).strip()
    unit = re.sub(r"^[\s:*\-–—]+", "", unit)
    if not name:
        return None
    return {"name": name, "value": mv.group(0), "unit": unit,
            "ref_low": low, "ref_high": high, "lab_flag": lab_flag}
_QUAL_LINE = re.compile(
    r"^\s*(?P<name>.+?)[\s:=\-–]+"
    r"(?P<val>\+{1,4}|отриц\w*|отр|положит\w*|полож|не\s*обнаружен\w*|обнаружен\w*|"
    r"нереактивн\w*|реактивн\w*|non-?reactive|reactive|"
    r"niereaktywn\w*|reaktywn\w*|ujemn\w*|dodatni\w*|nie\s*wykryt\w*|wykryt\w*|"
    r"negatywny|pozytywny|prawidłowy|prawidlowy|obecny|nieobecny|"
    r"следы|след|слабоположит\w*|negative|positive|neg|pos|trace|нет)\s*$",
    re.I)
_CA_MARKER = re.compile(r"\b(ca|са)[\s.-]*(125|19[\s.-]?9|15[\s.-]?3|72[\s.-]?4|242|50|199|153|724)\b", re.I)
def _glue_ca_markers(text: str) -> str:
    # «CA 125» / «CA 19-9» / «СА 15-3» → «CA125» и т.п., иначе парсер примет число за значение
    def r(m):
        return m.group(1) + re.sub(r"[\s.-]", "", m.group(2))
    return _CA_MARKER.sub(r, text)

# пояснительные/пороговые строки бланка — не результаты (иначе ложные показатели из OCR)
_BOILERPLATE = re.compile(
    r"(glikemia|cukrzyca|deficyt|stężenie|stezenie|zalecane|prawidłow|prawidlow|nieprawidłow|"
    r"ryzyk|poniżej|ponizej|powyżej|powyzej|wskazane|dla osób|dla osob|wg |według|norma[:\s]|"
    r"kryterium|przedcukrzycow|suboptymaln|optymaln|toksyczne|wysokie|w przypadku|"
    r"interpret|oznacza|metod|aparacie|wzór|wzor|przed godzin|po godzinie|OGTT|"
    r"дефицит|оптимальн|рекоменд|норма\s|при концентрации|порог)", re.I)

def parse_text(text: str) -> tuple[dict, list[dict]]:
    ctx: dict = {}
    out = []
    text = _glue_ca_markers(text)
    for raw in text.splitlines():
        line = raw.strip()
        # убрать коды «(ICD-9: 091)» / «(ICD-10 ...)» — иначе их цифры читаются как значение
        line = re.sub(r"\(ICD[^)]*\)", "", line, flags=re.I).strip()
        if not line or line.startswith("#"):
            continue
        if _BOILERPLATE.search(line):   # пояснительная строка, не результат
            continue
        low = line.lower()
        if re.search(r"\b(пол|sex)\b", low):
            if re.search(r"\b(муж|male|m|м)\b", low): ctx["sex"] = "male"
            elif re.search(r"\b(жен|female|f|ж)\b", low): ctx["sex"] = "female"
            continue
        if re.search(r"\b(возраст|age)\b", low):
            n = to_float(line)
            if n is not None: ctx["age"] = n
            continue
        if re.search(r"\b(вес|масса|weight)\b", low):
            n = to_float(line)
            if n is not None: ctx["weight"] = n
            continue
        if re.search(r"\b(рост|height)\b", low):
            n = to_float(line)
            if n is not None: ctx["height"] = n
            continue
        # качественная строка: «<название> +/++/отриц/следы/реактивный»
        qm = _QUAL_LINE.match(line)
        if qm:
            out.append({"name": qm.group("name"), "value": qm.group("val"), "unit": ""})
            continue
        rec = _split_measurement(line)
        if rec:
            # дифференциал ОАК: «Neutrofile 40,2 %» vs «Neutrofile 2,6 tys/ul» — по единице
            # решаем, это процент или абсолют, чтобы оба не ушли в один аналит
            nm = norm_name(rec["name"])
            diff = ("neutrofile", "limfocyty", "monocyty", "eozynofile", "bazofile",
                    "нейтрофилы", "лимфоциты", "моноциты", "эозинофилы", "базофилы")
            if nm in diff:
                u = norm_unit(rec.get("unit", ""))
                if u == "%":
                    rec["name"] = rec["name"] + " %"
                elif u and u != "%":
                    rec["name"] = rec["name"] + " абс"
            out.append(rec)
    return ctx, out

def parse_any(text: str) -> tuple[dict, list[dict]]:
    t = text.strip()
    if t.startswith("{") or t.startswith("["):
        try:
            return parse_json(t)
        except Exception:
            pass
    first = t.splitlines()[0] if t.splitlines() else ""
    if "," in first and len(first.split(",")) >= 2 and _split_measurement(first) is None:
        try:
            return parse_csv(t)
        except Exception:
            pass
    return parse_text(t)

# ───────────────────────────── оценка ────────────────────────────────────

_QUAL_NEG = {"neg", "negative", "отриц", "отрицательно", "отр", "нет", "необнаружено",
             "notdetected", "abs", "0", "-", "—", "none", "норма", "normal",
             # польский
             "ujemny", "ujemnie", "niewykryto", "niewykryta", "prawidłowy", "prawidlowy",
             "negatywny", "nieobecny", "brak"}
_QUAL_TRACE = {"trace", "следы", "след", "±", "+-", "-+", "слабоположит", "слабо", "sladowe", "ślad"}

def eval_qualitative(a: dict, key: str, raw: Any) -> dict:
    s = norm_name(str(raw)) if raw is not None else ""
    s_compact = s.replace(" ", "")
    normal_includes_trace = bool(a.get("normal_includes_trace"))
    # числовое значение у качественной серологии с порогом (HSV AU/ml, Chlamydia RU/ml и т.п.)
    num = to_float(raw)
    cut = a.get("pos_cutoff")
    is_wordy = bool(re.search(r"[a-zа-яąćęłńóśźż]", s))  # есть буквы → это словесный результат
    if num is not None and not is_wordy:
        if cut is not None:
            bl = a.get("borderline_low")
            if num >= cut:
                st, val_disp = "POSITIVE", f"{raw} (порог ≥{cut})"
            elif bl is not None and num >= bl:
                st, val_disp = "BORDERLINE", f"{raw} (погранично {bl}–{cut})"
            else:
                st, val_disp = "NORMAL", f"{raw} (порог ≥{cut})"
            notes = list(a.get("notes", []))
            if st == "BORDERLINE":
                notes.append("Пограничный результат — рекомендуется повтор через несколько недель.")
            return {
                "key": key, "name_ru": a.get("ru", key), "name_en": a.get("en", key),
                "panel": a.get("panel", "serology"),
                "input_value": raw, "input_unit": "", "value": val_disp, "unit": "",
                "ref_low": None, "ref_high": None, "range_source": "порог assay",
                "status": st, "unit_known": True, "notes": notes,
            }
        else:
            # число без известного порога — не гадаем (иначе ложный +)
            return {
                "key": key, "name_ru": a.get("ru", key), "name_en": a.get("en", key),
                "panel": a.get("panel", "serology"),
                "input_value": raw, "input_unit": "", "value": str(raw), "unit": "",
                "ref_low": None, "ref_high": None, "range_source": "нет порога",
                "status": "UNKNOWN", "unit_known": False,
                "notes": list(a.get("notes", [])) + ["Числовой результат без порога assay — сверьте с референсом бланка."],
            }
    # нереактивный / реактивный (RU/EN/PL) — порядок важен: «не...» проверяем первым
    if (s_compact.startswith("нереактив") or s_compact.startswith("nonreactive")
            or s_compact.startswith("non-reactive") or s_compact.startswith("niereaktyw")):
        status, val_disp = "NEG", "нереактивный"
    elif s_compact.startswith("реактив") or s_compact.startswith("reactive") or s_compact.startswith("reaktyw"):
        status, val_disp = "POS", "реактивный"
    elif s_compact.startswith("dodatni") or s_compact.startswith("pozytywny") or s_compact.startswith("wykryto"):
        status, val_disp = "POS", "положительно"
    elif s_compact in _QUAL_NEG or s_compact == "":
        status, val_disp = "NEG", "отриц"
    elif s_compact in _QUAL_TRACE or s_compact.startswith("след"):
        status, val_disp = ("NORMAL" if normal_includes_trace else "TRACE"), "следы"
    else:
        # «+», «++», «+++», «положит», «обнаружено», число КОЕ и т.п. → положительно
        status, val_disp = "POS", str(raw).strip()
    return {
        "key": key, "name_ru": a.get("ru", key), "name_en": a.get("en", key),
        "panel": a.get("panel", "urinalysis"),
        "input_value": raw, "input_unit": "", "value": val_disp, "unit": "",
        "ref_low": None, "ref_high": None, "range_source": "качеств.",
        "status": {"NEG": "NORMAL", "TRACE": "TRACE", "NORMAL": "NORMAL", "POS": "POSITIVE"}[status],
        "unit_known": True, "notes": list(a.get("notes", [])),
    }

def pick_range(a: dict, sex: Optional[str], age: Optional[float]) -> Optional[dict]:
    # Возрастные полосы приоритетнее, если возраст известен и полосы заданы.
    # Возраст — в годах (дробный для младенцев: 1 мес ≈ 0.083, 28 дней ≈ 0.077).
    if age is not None and a.get("age_ranges"):
        for band in a["age_ranges"]:
            if band.get("min", 0) <= age < band.get("max", 1e9):
                inner = band.get(sex) if sex and isinstance(band.get(sex), dict) else band
                return {"low": inner.get("low"), "high": inner.get("high"), "_band": band.get("label")}
    r = a.get("ranges", {})
    if sex in r:
        return r[sex]
    return r.get("default")

def evaluate(refbook: RefBook, ctx: dict, rows: list[dict]) -> dict:
    sex = (ctx.get("sex") or "").lower() or None
    results = []
    unmatched = []
    index_by_key: dict[str, dict] = {}

    for row in rows:
        name = row.get("name", "")
        if not name:
            continue
        raw_val = row.get("value")
        val = to_float(raw_val)
        key = refbook.match(name)
        if key is None:
            if raw_val not in (None, ""):
                unmatched.append({"name": name, "value": raw_val, "unit": row.get("unit", "")})
            continue
        a = refbook.analytes[key]

        # качественные показатели (ОАМ: +/++/отриц/следы) — обрабатываем отдельно
        if a.get("qualitative"):
            rec = eval_qualitative(a, key, raw_val)
            results.append(rec)
            index_by_key[key] = rec
            continue
        if val is None:
            continue
        conv_val, canon_unit, known = refbook.convert(key, val, row.get("unit", ""))

        # диапазон: приоритет — из бланка, затем возрастной/половой, затем взрослый.
        # Референс из бланка указан в ЕДИНИЦЕ БЛАНКА — конвертируем тем же множителем,
        # что и значение, иначе сравнение value(canonical) с ref(input) ложно.
        bl, bh = to_float(row.get("ref_low")), to_float(row.get("ref_high"))
        if known and val not in (None, 0):
            factor = conv_val / val
            if bl is not None: bl *= factor
            if bh is not None: bh *= factor
        age = ctx.get("age")
        rng = pick_range(a, sex, age)
        band_label = rng.get("_band") if rng else None
        if bl is not None or bh is not None:
            lo = bl if bl is not None else (rng.get("low") if rng else None)
            hi = bh if bh is not None else (rng.get("high") if rng else None)
            rng_source = "бланк"
        elif rng and (rng.get("low") is not None or rng.get("high") is not None):
            lo, hi = rng.get("low"), rng.get("high")
            rng_source = f"возраст: {band_label}" if band_label else "справочник"
        else:
            lo = hi = None
            rng_source = "нет"

        crit = a.get("critical") or {}
        cl, ch = crit.get("low"), crit.get("high")

        status = "UNKNOWN"
        if known:
            if cl is not None and conv_val <= cl:
                status = "CRITICAL_LOW"
            elif ch is not None and conv_val >= ch:
                status = "CRITICAL_HIGH"
            elif lo is not None and conv_val < lo:
                status = "LOW"
            elif hi is not None and conv_val > hi:
                status = "HIGH"
            elif lo is not None or hi is not None:
                status = "NORMAL"

        rec = {
            "key": key,
            "name_ru": a.get("ru", key),
            "name_en": a.get("en", key),
            "panel": a.get("panel", "other"),
            "input_value": val,
            "input_unit": row.get("unit", ""),
            "value": round(conv_val, 4) if known else val,
            "unit": canon_unit,
            "ref_low": lo,
            "ref_high": hi,
            "range_source": rng_source,
            "status": status,
            "unit_known": known,
            "lab_flag": row.get("lab_flag"),
            "notes": [],
        }
        # №6: сверка с флагом лаборатории (LIW: H/L). Расхождение = сигнал об
        # OCR-ошибке, разнице референсов или неверно распознанной единице.
        lf = row.get("lab_flag")
        if lf in ("H", "L") and known and status != "UNKNOWN":
            engine_dir = {"HIGH": "H", "CRITICAL_HIGH": "H", "LOW": "L", "CRITICAL_LOW": "L",
                          "NORMAL": "N"}.get(status, "N")
            if engine_dir != lf:
                rec["flag_mismatch"] = True
                rec["notes"].append(
                    f"⚠ Лаборатория пометила результат как '{lf}' (вне нормы), а расчёт дал "
                    f"'{'норма' if engine_dir=='N' else engine_dir}'. Вероятна разница референсов "
                    f"или ошибка распознавания — сверьте бланк.")
        if not known:
            rec["notes"].append(f"Единица входа '{row.get('unit','')}' не распознана — конверсия не выполнена, статус по исходному числу невозможен.")
        # sanity-guard: физиологически неправдоподобное значение после конверсии —
        # почти всегда опечатка в единице или вводе (напр. Hb «118 g/dl» → 1180 g/L).
        # Привязка к КРИТИЧЕСКИМ порогам: именно там кратное превышение означает
        # не патологию, а почти наверняка ошибку единицы/ввода (ложная «критика»).
        if known and status != "UNKNOWN":
            if ch is not None and ch > 0 and conv_val > ch * 4:
                rec["implausible"] = True
                rec["notes"].append(
                    f"⚠ Значение {round(conv_val,2)} {canon_unit} многократно превышает даже критический порог "
                    f"({ch}) — вероятна ошибка единицы или ввода (входные: {val} '{row.get('unit','')}'). Сверьте бланк.")
            elif cl is not None and cl > 0 and conv_val < cl / 4:
                rec["implausible"] = True
                rec["notes"].append(
                    f"⚠ Значение {round(conv_val,4)} {canon_unit} многократно ниже критического порога "
                    f"({cl}) — вероятна ошибка единицы или ввода (входные: {val} '{row.get('unit','')}'). Сверьте бланк.")
        if a.get("_note"):
            rec["notes"].append(a["_note"])
        if a.get("_note_infant") and age is not None and age < 0.25:
            rec["notes"].append(a["_note_infant"])
        if ctx.get("on_anticoagulant") and key == "inr":
            rec["notes"].append("Пациент на антикоагулянте — высокий МНО может быть целевым, а не патологией.")
        results.append(rec)
        index_by_key[key] = rec

    reconcile_differential(index_by_key)
    derived = compute_derived(index_by_key, sex, ctx.get("age"),
                              to_float(ctx.get("weight")), to_float(ctx.get("height")))
    summary = summarize(results)
    return {
        "context": {"sex": sex, "age": ctx.get("age")},
        "summary": summary,
        "results": results,
        "derived": derived,
        "unmatched": unmatched,
    }

# пары «процент ↔ абсолют» в лейкоформуле
_DIFF_PAIRS = [
    ("neutrophils_pct", "neutrophils_abs"),
    ("lymphocytes_pct", "lymphocytes_abs"),
    ("monocytes_pct", "monocytes_abs"),
    ("eosinophils_pct", "eosinophils_abs"),
    ("basophils_pct", "basophils_abs"),
]
def reconcile_differential(idx: dict[str, dict]) -> None:
    """Если процент в лейкоформуле вне нормы, а абсолютное число — в норме, помечает
    это как относительный (клинически обычно незначимый) сдвиг. Проценты связаны
    суммой до 100%, поэтому дефицит одной доли механически поднимает другие."""
    for pct_key, abs_key in _DIFF_PAIRS:
        p = idx.get(pct_key)
        a = idx.get(abs_key)
        if not p or not a:
            continue
        if p.get("status") in ("HIGH", "LOW") and a.get("status") == "NORMAL":
            direction = "повышен" if p["status"] == "HIGH" else "снижен"
            p.setdefault("notes", []).append(
                f"Отклонён только ПРОЦЕНТ ({direction}), абсолютное число в норме "
                f"({a['value']} {a['unit']}) — относительный сдвиг формулы, истинного "
                f"избытка/дефицита клеток нет; клинически обычно незначим.")
            p["relative_shift_only"] = True


def compute_derived(idx: dict[str, dict], sex: Optional[str],
                    age: Optional[float] = None, weight: Optional[float] = None,
                    height: Optional[float] = None) -> list[dict]:
    out = []
    def v(k):
        r = idx.get(k)
        val = r["value"] if r and r.get("unit_known") else None
        return val if isinstance(val, (int, float)) else None

    # ── СКФ из креатинина: CKD-EPI 2021 (race-free) и Кокрофт-Голт ──
    scr = v("creatinine")  # umol/L
    if scr and scr > 0 and age is not None and age >= 18 and sex in ("male", "female"):
        scr_mgdl = scr / 88.42
        k = 0.7 if sex == "female" else 0.9
        al = -0.241 if sex == "female" else -0.302
        fem = 1.012 if sex == "female" else 1.0
        ratio = scr_mgdl / k
        egfr = 142 * (min(ratio, 1) ** al) * (max(ratio, 1) ** -1.200) * (0.9938 ** age) * fem
        egfr = round(egfr, 1)
        if egfr >= 90: st = "G1"
        elif egfr >= 60: st = "G2"
        elif egfr >= 45: st = "G3a"
        elif egfr >= 30: st = "G3b"
        elif egfr >= 15: st = "G4"
        else: st = "G5"
        out.append({"name": "СКФ расч. (CKD-EPI 2021)", "value": egfr, "unit": "мл/мин/1.73м²",
                    "note": f"Стадия ХБП {st}. Для детей формула не валидна — используй Шварца."})
        if weight and weight > 0:
            crcl = (140 - age) * weight / (72 * scr_mgdl)
            if sex == "female":
                crcl *= 0.85
            out.append({"name": "Клиренс креатинина (Кокрофт-Голт)", "value": round(crcl, 1),
                        "unit": "мл/мин",
                        "note": "Для подбора доз препаратов; по фактической массе тела."})

    # ── СКФ у детей: формула Шварца (bedside, 2009) ──
    if scr and scr > 0 and age is not None and age < 18 and height and height > 0:
        scr_mgdl = scr / 88.42
        egfr_s = round(0.413 * height / scr_mgdl, 1)
        out.append({"name": "СКФ у детей (Шварц, bedside)", "value": egfr_s, "unit": "мл/мин/1.73м²",
                    "note": "Педиатрическая формула (0.413 × рост[см] / креатинин[мг/дл]); для детей точнее CKD-EPI."})

    # ── Скорректированный кальций по альбумину ──
    ca, alb = v("calcium"), v("albumin")
    if ca and alb:
        ca_corr = round(ca + 0.02 * (40 - alb), 2)
        out.append({"name": "Кальций скорр. по альбумину", "value": ca_corr, "unit": "mmol/L",
                    "note": "При альбумине <40 г/л истинный кальций выше измеренного; норма 2.15–2.55."})

    # ── Индексы из формулы крови: NLR / PLR / SII ──
    neu, lym, plt, wbc = v("neutrophils_abs"), v("lymphocytes_abs"), v("platelets"), v("wbc")
    neu_pct, lym_pct = v("neutrophils_pct"), v("lymphocytes_pct")
    if neu is None and wbc and neu_pct:   # оценка абс. из % × WBC, если абс. не дан
        neu = wbc * neu_pct / 100
    if lym is None and wbc and lym_pct:
        lym = wbc * lym_pct / 100
    if neu and lym and lym > 0:
        nlr = round(neu / lym, 2)
        out.append({"name": "НЛИ (нейтрофил/лимфоцитарный индекс, NLR)", "value": nlr,
                    "note": "Норма ~1–2; >3 — системное воспаление/стресс; высокие значения — неблагоприятный прогноз."})
    if plt and lym and lym > 0:
        out.append({"name": "ТЛИ (тромбоцит/лимфоцитарный, PLR)", "value": round(plt / lym, 1),
                    "note": "Повышение — воспаление/тромбовоспалительная активность."})
    if plt and neu and lym and lym > 0:
        out.append({"name": "Системный иммуно-воспалит. индекс (SII)", "value": round(plt * neu / lym, 0),
                    "note": "SII = тромбоциты × нейтрофилы / лимфоциты; рост — системное воспаление."})

    # ── Триглицерид-глюкозный индекс (TyG) — маркёр инсулинорезистентности ──
    tg, glu = v("triglycerides"), v("glucose")
    if tg and glu and tg > 0 and glu > 0:
        tg_mgdl = tg / 0.01129
        glu_mgdl = glu / 0.0555
        tyg = round(math.log(tg_mgdl * glu_mgdl / 2), 2)
        out.append({"name": "Индекс TyG (триглицерид-глюкозный)", "value": tyg,
                    "note": "Маркёр инсулинорезистентности; ориентир >4.49 (≈8.5 по части шкал) — вероятная ИР. Нужны ТГ и глюкоза натощак."})

    # ── Гликемический статус (ADA/ВОЗ): натощак + ГТТ 2ч + HbA1c ──
    glu_f = v("glucose")
    glu_2h = v("glucose_2h")
    a1c = v("hba1c")
    gverdict = []
    if glu_f is not None:
        if glu_f >= 7.0: gverdict.append("глюкоза натощак ≥7.0 — критерий диабета")
        elif glu_f >= 6.1: gverdict.append("нарушенная гликемия натощак (6.1–6.9)")
    if glu_2h is not None:
        if glu_2h >= 11.1: gverdict.append("ГТТ 2ч ≥11.1 — критерий диабета")
        elif glu_2h >= 7.8: gverdict.append("нарушенная толерантность к глюкозе (7.8–11.0)")
    if a1c is not None:
        if a1c >= 6.5: gverdict.append("HbA1c ≥6.5% — критерий диабета")
        elif a1c >= 5.7: gverdict.append("предиабет по HbA1c (5.7–6.4%)")
    if gverdict:
        out.append({"name": "Гликемический статус", "value": "; ".join(gverdict),
                    "note": "Диагноз диабета требует подтверждения (повтор или два разных критерия). Пороги ADA/ВОЗ."})

    # ── Индекс свободных андрогенов (FAI) ──
    testo, shbg = v("testosterone_total"), v("shbg")
    if testo and shbg and shbg > 0:
        fai = round(testo / shbg * 100, 1)
        hint = ("норма муж ~30–150" if sex == "male"
                else ("норма жен <5; выше — гиперандрогения (СПКЯ)" if sex == "female"
                      else "муж ~30–150, жен <5"))
        out.append({"name": "Индекс свободных андрогенов (FAI)", "value": fai,
                    "note": f"FAI = тестостерон общий / ГСПГ × 100. {hint}."})

    # ── Соотношение Apo B/A1 — атерогенный риск ──
    apob, apoa1 = v("apo_b"), v("apo_a1")
    if apob and apoa1 and apoa1 > 0:
        out.append({"name": "Соотношение Apo B/A1", "value": round(apob / apoa1, 2),
                    "note": "Атерогенный индекс; ориентир риска: муж <0.9, жен <0.8 — выше значит выше риск ССЗ."})

    # ── HOMA-IR (инсулинорезистентность) ──
    ins = v("insulin")  # mU/L
    if ins and glu and ins > 0 and glu > 0:
        homa = round(ins * glu / 22.5, 2)
        out.append({"name": "HOMA-IR", "value": homa,
                    "note": "Инсулинорезистентность: ориентир <2.0–2.7 норма, выше — вероятная ИР. Инсулин и глюкоза — строго натощак."})

    # ── Дефицит железа по Ганзони (мг для восполнения) ──
    hb = v("hemoglobin")  # g/L
    if hb and weight and weight > 0 and hb < 150:
        target = 150  # целевой Hb, g/L (≈15 g/dL)
        deficit = weight * (target - hb) / 10 * 2.4 + 500  # +500 мг депо (взрослый)
        if deficit > 0:
            out.append({"name": "Дефицит железа по Ганзони", "value": round(deficit), "unit": "мг",
                        "note": "Расчётная доза элементарного железа для восполнения (целевой Hb 150 г/л + депо 500 мг). Ориентир для терапии, не диагноз дефицита."})

    # ── Коэффициент де Ритиса, липидные индексы, анионный провал, трансферрин, глобулины ──
    ast, alt = v("ast"), v("alt")
    if ast and alt and alt > 0:
        ratio = round(ast / alt, 2)
        note = "De Ritis >2 — ориентир на алкогольное/митохондриальное поражение; <1 — чаще вирусное/метаболическое." if ratio else ""
        out.append({"name": "Коэффициент де Ритиса (АСТ/АЛТ)", "value": ratio, "note": note})

    chol, hdl = v("cholesterol_total"), v("hdl")
    if chol and hdl:
        non_hdl = round(chol - hdl, 2)
        out.append({"name": "Не-ЛПВП холестерин", "value": non_hdl, "unit": "mmol/L",
                    "note": "Цель обычно <3.4 ммоль/л; жёстче при высоком риске."})
        if hdl > 0:
            out.append({"name": "Индекс атерогенности (Chol/HDL)", "value": round(chol / hdl, 2),
                        "note": "Ориентир <4.0."})

    na, cl, hco3b = v("sodium"), v("chloride"), v("bicarbonate")
    if na and cl and hco3b:
        ag = round(na - (cl + hco3b), 1)
        out.append({"name": "Анионный провал", "value": ag, "unit": "mmol/L",
                    "note": "Норма ~8–12; рост — метаболический ацидоз с высоким AG."})

    egfr = v("egfr")
    if egfr is not None:
        if egfr >= 90: stage = "G1 (норма/высокая)"
        elif egfr >= 60: stage = "G2 (незначительно снижена)"
        elif egfr >= 45: stage = "G3a"
        elif egfr >= 30: stage = "G3b"
        elif egfr >= 15: stage = "G4 (выраженно снижена)"
        else: stage = "G5 (почечная недостаточность)"
        out.append({"name": "Стадия ХБП по СКФ (из бланка)", "value": stage})

    iron, tibc = v("iron"), v("tibc")
    if iron and tibc and tibc > 0:
        out.append({"name": "Насыщение трансферрина (расч.)", "value": round(iron / tibc * 100, 1), "unit": "%",
                    "note": "<20% — дефицит железа; >45% — перегрузка."})

    tp, alb2 = v("total_protein"), v("albumin")
    if tp and alb2 and (tp - alb2) > 0:
        out.append({"name": "Глобулины (расч.)", "value": round(tp - alb2, 1), "unit": "g/L",
                    "note": f"A/G = {round(alb2/(tp-alb2),2)}"})

    # ── Пунктаты: критерии Лайта (экссудат/транссудат) и SAAG ──
    ef_prot, ef_ldh = v("effusion_protein"), v("effusion_ldh")
    s_prot, s_ldh, s_alb = v("total_protein"), v("ldh"), v("albumin")
    ef_alb = v("effusion_albumin")
    light = []
    if ef_prot and s_prot and s_prot > 0:
        light.append(ef_prot / s_prot > 0.5)
    if ef_ldh and s_ldh and s_ldh > 0:
        light.append(ef_ldh / s_ldh > 0.6)
        light.append(ef_ldh > (2 / 3) * 280)  # >2/3 верхней нормы ЛДГ сыворотки (~280)
    if light:
        verdict = "ЭКССУДАТ (хотя бы один критерий Лайта +)" if any(light) else "транссудат (все критерии Лайта −)"
        out.append({"name": "Критерии Лайта (пунктат)", "value": verdict,
                    "note": "Экссудат — воспаление/инфекция/опухоль; транссудат — СН/цирроз/нефроз. Нужны парные белок и ЛДГ выпота и сыворотки."})
    if ef_alb is not None and s_alb is not None:
        saag = round(s_alb - ef_alb, 1)
        v2 = "≥11 г/л — портальная гипертензия (цирроз, СН)" if saag >= 11 else "<11 г/л — непортальная причина (карциноматоз, ТБ, панкреатит)"
        out.append({"name": "SAAG (сывороточно-асцитический градиент)", "value": saag, "unit": "g/L",
                    "note": v2})

    # ── Осмолярный провал ──
    osm_meas = v("osmolality_measured")
    na, glu_o, urea_o = v("sodium"), v("glucose"), v("urea")
    if osm_meas and na:
        calc = 2 * na + (glu_o or 0) + (urea_o or 0)
        gap = round(osm_meas - calc, 1)
        if gap > 10:
            out.append({"name": "Осмолярный провал", "value": gap, "unit": "mOsm/kg",
                        "note": ">10 — осмотически активный токсин (этиленгликоль, метанол, этанол); расчётная осмолярность = 2×Na+глюкоза+мочевина."})

    # ── Коэффициент глюкозы ликвор/кровь (менингит) ──
    csf_glu = v("csf_glucose")
    if csf_glu and glu and glu > 0:
        ratio = round(csf_glu / glu, 2)
        verdict = "норма (≥0.6)" if ratio >= 0.6 else ("пограничный" if ratio >= 0.4 else "снижен — бактериальный/ТБ/грибковый менингит")
        out.append({"name": "Коэффициент глюкозы ликвор/кровь", "value": ratio,
                    "note": f"{verdict}. Сравнение с одновременной глюкозой крови; изолированное значение ликвора менее информативно."})

    # ── Индекс CK-MB/CK (миокардиальное происхождение) ──
    ckmb, ck = v("ck_mb"), v("ck")
    if ckmb and ck and ck > 0:
        mb_idx = round(ckmb / ck * 100, 1)
        hint = "≥5% — в пользу поражения миокарда" if mb_idx >= 5 else "<5% — скорее скелетная мышца (рабдомиолиз/травма)"
        out.append({"name": "Относительный индекс CK-MB/CK", "value": mb_idx, "unit": "%",
                    "note": f"{hint}. Тропонин специфичнее для миокарда."})

    # ── Ретикулоцитарный индекс (RPI) — адекватность эритропоэза ──
    retic = v("reticulocytes_pct")
    hct = v("hematocrit")
    if retic and hct and hct > 0:
        corrected = retic * hct / 45.0
        if hct >= 35: mf = 1.0
        elif hct >= 25: mf = 1.5
        elif hct >= 20: mf = 2.0
        else: mf = 2.5
        rpi = round(corrected / mf, 2)
        verdict = ("адекватный костномозговой ответ (гемолиз/кровопотеря)" if rpi >= 2
                   else "неадекватный ответ — гипопролиферативная анемия (дефицит/ХБП/аплазия)")
        out.append({"name": "Ретикулоцитарный индекс (RPI)", "value": rpi,
                    "note": f"Порог ~2: {verdict}. Корректирован на Hct и созревание — истиннее «сырого» % ретикулоцитов."})

    # ── Кислотно-щелочное состояние (КЩС) ──
    ab = acid_base(v("ph_abg"), v("pco2"), v("bicarbonate"))
    if ab:
        out.append(ab)
    return out


# ─────────────────── delta-check (динамика во времени) ───────────────────
# RCV (reference change value) — порог клинически значимого изменения между
# заборами, % (аналитическая + биологическая вариация). Для ключевых аналитов —
# литературные значения; для прочих действует DEFAULT_RCV.
_RCV = {
    "creatinine": 14, "egfr": 15, "urea": 25, "uric_acid": 18,
    "alt": 45, "ast": 30, "ggt": 23, "alp": 12, "bilirubin_total": 25,
    "tsh": 50, "ft4": 12, "ft3": 15,
    "hemoglobin": 8, "hematocrit": 9, "rbc": 7, "platelets": 11, "wbc": 15,
    "neutrophils_abs": 30, "lymphocytes_abs": 25,
    "potassium": 14, "sodium": 1.5, "chloride": 3, "calcium": 5, "magnesium": 10,
    "glucose": 15, "hba1c": 5, "insulin": 35,
    "crp": 50, "ferritin": 20, "iron": 30, "procalcitonin": 50,
    "cholesterol_total": 15, "ldl": 20, "hdl": 12, "triglycerides": 25,
    "psa_total": 30, "d_dimer": 30, "troponin": 30, "nt_probnp": 40, "bnp": 40,
    "inr": 12, "fibrinogen": 15, "albumin": 8, "total_protein": 6,
    "vitamin_d": 25, "vitamin_b12": 25, "folate": 30,
    "anti_tpo": 40, "cortisol": 40, "testosterone_total": 25,
}
DEFAULT_RCV = 25.0


def compute_dynamics(results: list[dict], history: list[dict], current_date: str) -> list[dict]:
    """Сравнивает текущие показатели с последними прежними того же пациента.
    history: записи {date, key, value, unit, status}. Возвращает значимые сдвиги."""
    # последняя прежняя запись по каждому key (по дате, строго раньше текущей)
    prev: dict[str, dict] = {}
    for rec in sorted(history, key=lambda r: r.get("date", "")):
        if rec.get("date", "") >= current_date:
            continue
        prev[rec["key"]] = rec
    out = []
    for r in results:
        if not r.get("unit_known") or not isinstance(r.get("value"), (int, float)):
            continue
        p = prev.get(r["key"])
        if not p or not isinstance(p.get("value"), (int, float)):
            continue
        old, new = p["value"], r["value"]
        if old == 0:
            continue
        delta_pct = round((new - old) / abs(old) * 100, 1)
        rcv = _RCV.get(r["key"], DEFAULT_RCV)
        status_change = p.get("status") != r["status"] and \
            {p.get("status"), r["status"]} != {"NORMAL"}  # смена норма↔отклонение/крит
        significant = abs(delta_pct) >= rcv or status_change
        if not significant:
            continue
        direction = "↑" if new > old else "↓"
        out.append({
            "key": r["key"], "name_ru": r["name_ru"],
            "prev_value": old, "prev_date": p.get("date"), "value": new, "unit": r["unit"],
            "delta_pct": delta_pct, "rcv": rcv, "direction": direction,
            "prev_status": p.get("status"), "status": r["status"],
            "status_change": status_change,
        })
    return out


def _history_path(patient_id: str) -> str:
    safe = re.sub(r"[^\w.-]", "_", patient_id)[:64]
    return os.path.join(ROOT_DIR, "history", f"{safe}.jsonl")


def load_history(patient_id: str) -> list[dict]:
    path = _history_path(patient_id)
    out = []
    if os.path.exists(path):
        for line in open(path, encoding="utf-8"):
            line = line.strip()
            if line:
                try:
                    out.append(json.loads(line))
                except Exception:
                    pass
    return out


def save_run(patient_id: str, date: str, results: list[dict]):
    path = _history_path(patient_id)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "a", encoding="utf-8") as f:
        for r in results:
            if r.get("unit_known") and isinstance(r.get("value"), (int, float)):
                f.write(json.dumps({"date": date, "key": r["key"], "value": r["value"],
                                    "unit": r["unit"], "status": r["status"]},
                                   ensure_ascii=False) + "\n")


_SPARK = "▁▂▃▄▅▆▇█"
def _sparkline(vals):
    nums = [x for x in vals if isinstance(x, (int, float))]
    if len(nums) < 2:
        return ""
    lo, hi = min(nums), max(nums)
    if hi == lo:
        return _SPARK[0] * len(vals)
    out = []
    for x in vals:
        if isinstance(x, (int, float)):
            out.append(_SPARK[int((x - lo) / (hi - lo) * (len(_SPARK) - 1))])
        else:
            out.append(" ")
    return "".join(out)


def build_trend(results: list[dict], history: list[dict], current_date: str,
                name_by_key: dict) -> dict:
    """Таблица динамики по датам: для каждого показателя, встречавшегося более чем
    в один день, собираем значения по датам (+ спарклайн)."""
    series: dict[str, dict] = {}  # key -> {date: (value, unit, status)}
    for rec in history:
        if rec.get("date", "") >= current_date:
            continue
        series.setdefault(rec["key"], {})[rec["date"]] = (rec.get("value"), rec.get("unit"), rec.get("status"))
    for r in results:
        if r.get("unit_known") and isinstance(r.get("value"), (int, float)):
            series.setdefault(r["key"], {})[current_date] = (r["value"], r["unit"], r["status"])
    dates = sorted({d for s in series.values() for d in s})
    rows = []
    for key, byd in series.items():
        if len(byd) < 2:
            continue
        vals = [byd.get(d, (None,))[0] for d in dates]
        units = [byd[d][1] for d in byd]
        rows.append({
            "key": key, "name_ru": name_by_key.get(key, key),
            "unit": units[-1] if units else "",
            "values": vals, "spark": _sparkline(vals),
            "last_status": byd[sorted(byd)[-1]][2],
        })
    rows.sort(key=lambda r: r["name_ru"])
    return {"dates": dates, "rows": rows}


def acid_base(ph: Optional[float], pco2: Optional[float], hco3: Optional[float]) -> Optional[dict]:
    """Упрощённая интерпретация КЩС по pH, pCO2 (mmHg), HCO3 (mmol/L)."""
    if ph is None or (pco2 is None and hco3 is None):
        return None
    parts = []


def acid_base(ph: Optional[float], pco2: Optional[float], hco3: Optional[float]) -> Optional[dict]:
    """Упрощённая интерпретация КЩС по pH, pCO2 (mmHg), HCO3 (mmol/L)."""
    if ph is None or (pco2 is None and hco3 is None):
        return None
    parts = []
    acidemia = ph < 7.35
    alkalemia = ph > 7.45
    if acidemia:
        if pco2 is not None and pco2 > 45:
            parts.append("респираторный ацидоз")
        if hco3 is not None and hco3 < 22:
            parts.append("метаболический ацидоз")
        if not parts:
            parts.append("ацидемия (уточни первичное звено)")
    elif alkalemia:
        if pco2 is not None and pco2 < 35:
            parts.append("респираторный алкалоз")
        if hco3 is not None and hco3 > 26:
            parts.append("метаболический алкалоз")
        if not parts:
            parts.append("алкалемия (уточни первичное звено)")
    else:
        if pco2 is not None and hco3 is not None and (pco2 > 45 or pco2 < 35 or hco3 > 26 or hco3 < 22):
            parts.append("pH в норме при сдвиге pCO₂/HCO₃ — компенсированное или смешанное нарушение")
        else:
            parts.append("КЩС в пределах нормы")
    comp = ""
    if acidemia and hco3 is not None and pco2 is not None:
        if "респираторный" in parts[0] and hco3 > 26:
            comp = " с метаболической компенсацией"
        elif "метаболический" in parts[0] and pco2 < 35:
            comp = " с респираторной компенсацией"
    if alkalemia and hco3 is not None and pco2 is not None:
        if "респираторный" in parts[0] and hco3 < 22:
            comp = " с метаболической компенсацией"
        elif "метаболический" in parts[0] and pco2 > 45:
            comp = " с респираторной компенсацией"
    return {"name": "Интерпретация КЩС", "value": " + ".join(parts) + comp,
            "note": f"pH={ph}, pCO₂={pco2}, HCO₃={hco3}. При метаболическом ацидозе смотри анионный провал."}


def summarize(results: list[dict]) -> dict:
    order = {"CRITICAL_HIGH": 0, "CRITICAL_LOW": 0, "HIGH": 0, "LOW": 0,
             "POSITIVE": 0, "BORDERLINE": 0, "TRACE": 0, "NORMAL": 0, "UNKNOWN": 0}
    for r in results:
        order[r["status"]] = order.get(r["status"], 0) + 1
    crit = [r for r in results if r["status"].startswith("CRITICAL")]
    abn = [r for r in results if r["status"] in ("HIGH", "LOW", "POSITIVE")]
    return {
        "total": len(results),
        "critical": len(crit),
        "abnormal": len(abn),
        "normal": order["NORMAL"],
        "counts": order,
        "critical_items": [r["name_ru"] for r in crit],
    }

# ───────────────────────────── вывод ─────────────────────────────────────

_FLAG = {
    "CRITICAL_HIGH": "⛔↑ КРИТ", "CRITICAL_LOW": "⛔↓ КРИТ",
    "HIGH": "▲ выше", "LOW": "▼ ниже", "NORMAL": "✓ норма", "UNKNOWN": "? н/д",
    "POSITIVE": "⊕ полож", "TRACE": "± следы", "BORDERLINE": "~ погранич",
}
_PANEL_RU = {
    "cbc": "Общий анализ крови", "metabolic": "Углеводный/метаболизм",
    "renal": "Почки", "electrolytes": "Электролиты", "liver": "Печень/ЖКТ",
    "lipid": "Липиды", "thyroid": "Щитовидная железа", "iron": "Обмен железа",
    "vitamins": "Витамины", "inflammation": "Воспаление", "cardiac": "Сердце",
    "coagulation": "Коагулограмма", "urinalysis": "Общий анализ мочи",
    "hormones": "Гормоны", "blood_gas": "Газы крови / КЩС",
    "tumor_markers": "Онкомаркеры", "immunology": "Иммунология / ревмопанель",
    "bone": "Костный метаболизм", "serology": "Серология / инфекции",
    "protein_electrophoresis": "Белковые фракции (электрофорез)",
    "catecholamines": "Катехоламины / метанефрины", "coag_factors": "Факторы свёртывания",
    "micronutrients": "Витамины / микроэлементы", "autoantibodies": "Аутоантитела",
    "allergy": "Аллергология", "csf": "Ликвор (СМЖ)",
    "drug_monitoring": "Лекарственный мониторинг", "toxicology": "Токсикология",
    "effusion": "Пунктаты серозных полостей", "metabolites": "Специальные метаболиты",
    "other": "Прочее",
}

def render_text(report: dict) -> str:
    L = []
    ctx = report["context"]
    s = report["summary"]
    sx = {"male": "муж", "female": "жен"}.get(ctx.get("sex"), "—")
    L.append(f"РАСШИФРОВКА АНАЛИЗОВ   пол: {sx}   возраст: {ctx.get('age') or '—'}")
    L.append("=" * 64)
    L.append(f"Показателей: {s['total']}   критических: {s['critical']}   "
             f"отклонений: {s['abnormal']}   норма: {s['normal']}")
    if s["critical_items"]:
        L.append("⛔ КРИТИЧЕСКИЕ ЗНАЧЕНИЯ: " + ", ".join(s["critical_items"]))
    L.append("")

    by_panel: dict[str, list] = {}
    for r in report["results"]:
        by_panel.setdefault(r["panel"], []).append(r)
    for panel, items in by_panel.items():
        L.append(f"── {_PANEL_RU.get(panel, panel)} " + "─" * max(0, 48 - len(_PANEL_RU.get(panel, panel))))
        for r in items:
            rng = ""
            if r["ref_low"] is not None or r["ref_high"] is not None:
                rng = f"[{r['ref_low'] if r['ref_low'] is not None else '—'}–{r['ref_high'] if r['ref_high'] is not None else '—'}]"
            L.append(f"  {_FLAG[r['status']]:<10} {r['name_ru']:<34} "
                     f"{r['value']} {r['unit']} {rng}")
            for n in r["notes"]:
                L.append(f"             ↳ {n}")
        L.append("")

    if report["derived"]:
        L.append("── Производные показатели " + "─" * 30)
        for d in report["derived"]:
            unit = d.get("unit", "")
            L.append(f"  • {d['name']}: {d['value']} {unit}".rstrip())
            if d.get("note"):
                L.append(f"      {d['note']}")
        L.append("")

    dyn = report.get("dynamics")
    if dyn:
        pat = report.get("patient", {})
        L.append(f"── Динамика (delta-check) с прошлых заборов " + "─" * 13)
        for d in dyn:
            flag = "⚠ смена статуса" if d["status_change"] else "значимо"
            sc = (f"  [{_FLAG.get(d['prev_status'],'?').split()[-1]}→{_FLAG.get(d['status'],'?').split()[-1]}]"
                  if d["status_change"] else "")
            L.append(f"  {d['direction']} {d['name_ru']}: {d['prev_value']} → {d['value']} {d['unit']}  "
                     f"({'+' if d['delta_pct']>=0 else ''}{d['delta_pct']}%, порог RCV {d['rcv']}%, {flag}){sc}")
            L.append(f"      относительно {d['prev_date']}")
        L.append("")
    elif report.get("patient"):
        pr = report["patient"].get("prior_runs", 0)
        L.append(f"── Динамика: {'значимых изменений нет' if pr else 'это первый сохранённый забор — не с чем сравнивать'} ──")
        L.append("")

    tr = report.get("trend")
    if tr and tr["rows"]:
        dts = [d[2:] if len(d) == 10 else d for d in tr["dates"]]  # YY-MM-DD коротко
        L.append("── Тренд по датам " + "─" * 39)
        L.append("  " + " " * 26 + "  ".join(f"{d:>8}" for d in dts) + "   тренд")
        for row in tr["rows"]:
            cells = []
            for val in row["values"]:
                cells.append(f"{val:>8}" if isinstance(val, (int, float)) else f"{'·':>8}")
            nm = (row["name_ru"][:24] + "…") if len(row["name_ru"]) > 25 else row["name_ru"]
            L.append(f"  {nm:<26}" + "  ".join(cells) + f"   {row['spark']}")
        L.append("")

    if report["unmatched"]:
        L.append("── Не распознаны (нет в справочнике) " + "─" * 20)
        for u in report["unmatched"]:
            L.append(f"  ? {u['name']} = {u['value']} {u['unit']}")
        L.append("")
    return "\n".join(L)

# ───────────────────────────── CLI ───────────────────────────────────────

def main(argv=None):
    ap = argparse.ArgumentParser(description="Детерминированная расшифровка анализов.")
    ap.add_argument("input", help="файл с бланком, либо '-' для stdin")
    ap.add_argument("--sex", choices=["male", "female"], help="пол (переопределяет бланк)")
    ap.add_argument("--age", type=float, help="возраст в годах (дробный для младенцев: 0.077 ≈ 28 дней)")
    ap.add_argument("--weight", type=float, help="вес, кг (для клиренса по Кокрофт-Голту)")
    ap.add_argument("--height", type=float, help="рост, см (для педиатрической СКФ по Шварцу)")
    ap.add_argument("--refs", default=DEFAULT_REFS, help="путь к reference_ranges.json")
    ap.add_argument("--json", action="store_true", help="вывести машинный JSON вместо текста")
    ap.add_argument("--learn", nargs="?", const="__default__", metavar="LOGPATH",
                    help="логировать нераспознанные показатели для последующего обучения (по умолчанию learned/unmatched.log)")
    ap.add_argument("--patient", metavar="ID",
                    help="идентификатор пациента: включает delta-check (сравнение с прошлыми заборами) и сохраняет результат в историю")
    ap.add_argument("--date", metavar="YYYY-MM-DD", help="дата забора (по умолчанию сегодня)")
    ap.add_argument("--no-save", action="store_true", help="не сохранять этот прогон в историю пациента")
    args = ap.parse_args(argv)

    text = sys.stdin.read() if args.input == "-" else open(args.input, "r", encoding="utf-8").read()
    refbook = RefBook(args.refs)
    ctx, rows = parse_any(text)
    if args.sex: ctx["sex"] = args.sex
    if args.age is not None: ctx["age"] = args.age
    if args.weight is not None: ctx["weight"] = args.weight
    if args.height is not None: ctx["height"] = args.height

    report = evaluate(refbook, ctx, rows)

    # delta-check: динамика относительно прошлых заборов того же пациента
    if args.patient:
        import datetime
        run_date = args.date or datetime.date.today().isoformat()
        history = load_history(args.patient)
        report["dynamics"] = compute_dynamics(report["results"], history, run_date)
        name_by_key = {r["key"]: r["name_ru"] for r in report["results"]}
        name_by_key.update({k: refbook.analytes[k].get("ru", k) for k in refbook.analytes})
        report["trend"] = build_trend(report["results"], history, run_date, name_by_key)
        report["patient"] = {"id": args.patient, "date": run_date, "prior_runs": len({h.get("date") for h in history})}
        if not args.no_save:
            save_run(args.patient, run_date, report["results"])

    # обучение: логируем нераспознанные показатели (append JSONL) для learn.py
    if args.learn is not None and report["unmatched"]:
        import datetime
        logpath = (os.path.join(ROOT_DIR, "learned", "unmatched.log")
                   if args.learn == "__default__" else args.learn)
        os.makedirs(os.path.dirname(logpath), exist_ok=True)
        with open(logpath, "a", encoding="utf-8") as f:
            for u in report["unmatched"]:
                f.write(json.dumps({"ts": datetime.datetime.now().isoformat(timespec="seconds"),
                                    "name": u["name"], "value": u["value"], "unit": u["unit"]},
                                   ensure_ascii=False) + "\n")

    if args.json:
        print(json.dumps(report, ensure_ascii=False, indent=2))
    else:
        print(render_text(report))
    # код возврата: 2 если есть критические, 1 если есть отклонения, 0 если всё в норме
    s = report["summary"]
    return 2 if s["critical"] else (1 if s["abnormal"] else 0)

if __name__ == "__main__":
    sys.exit(main())
