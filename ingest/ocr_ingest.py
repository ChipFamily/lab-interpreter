#!/usr/bin/env python3
"""
ocr_ingest.py — вход PDF/изображение → текст показателей + авто-метаданные.

Полностью ОФЛАЙН (приватность). Многоуровневый OCR с автоопределением лучшего
доступного движка и graceful-откатом:
    PaddleOCR  →  EasyOCR  →  Tesseract(rus+eng+pol)
Каждый кадр предобрабатывается (оттенки серого, апскейл, контраст, бинаризация),
что заметно поднимает качество даже у Tesseract на сканах бланков.

Дополнительно извлекает из бланка (№4): дату забора и идентификатор пациента —
чтобы delta-check по папке PDF запускался автоматически.

Использование:
    python3 ingest/ocr_ingest.py BLANK.pdf                 # печатает текст
    python3 ingest/ocr_ingest.py BLANK.pdf --meta          # + дата/пациент (JSON)
    python3 ingest/ocr_ingest.py BLANK.pdf --engine tesseract
    python3 ingest/ocr_ingest.py scan.png                  # изображение тоже можно

Связка с движком:
    python3 ingest/ocr_ingest.py v.pdf > /tmp/v.txt
    python3 engine/interpret.py /tmp/v.txt --patient <id> --date <date>
Или всё сразу:
    python3 ingest/ocr_ingest.py v.pdf --run               # OCR → interpret, авто дата/пациент
"""
import sys, os, re, json, argparse, subprocess, tempfile, datetime

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

OCR_LANGS = "rus+eng+pol"


# ─────────────────────────── предобработка кадра ───────────────────────────
def _preprocess(pil_img):
    """Оттенки серого + апскейл мелких сканов + автоконтраст. Без тяжёлых зависимостей."""
    try:
        from PIL import ImageOps
        img = pil_img.convert("L")
        w, h = img.size
        if max(w, h) < 2200:                      # мелкий скан → апскейл
            k = 2200 / max(w, h)
            img = img.resize((int(w * k), int(h * k)))
        img = ImageOps.autocontrast(img)
        return img
    except Exception:
        return pil_img


def _pages(path):
    """PDF → список PIL-изображений; изображение → один кадр."""
    ext = os.path.splitext(path)[1].lower()
    if ext in (".pdf",):
        from pdf2image import convert_from_path
        return convert_from_path(path, dpi=300)
    from PIL import Image
    return [Image.open(path)]


# ─────────────────────────── движки OCR (по убыванию качества) ──────────────
def _try_paddle(images):
    from paddleocr import PaddleOCR  # noqa
    ocr = PaddleOCR(use_angle_cls=True, lang="ru", show_log=False)
    out = []
    for im in images:
        import numpy as np
        res = ocr.ocr(np.array(_preprocess(im)), cls=True)
        for block in (res or []):
            for line in (block or []):
                out.append(line[1][0])
    return "\n".join(out)


def _try_easyocr(images):
    import easyocr, numpy as np  # noqa
    reader = easyocr.Reader(["ru", "en", "pl"], gpu=False, verbose=False)
    out = []
    for im in images:
        for det in reader.readtext(np.array(_preprocess(im)), detail=0, paragraph=True):
            out.append(det)
    return "\n".join(out)


def _try_tesseract(images):
    import pytesseract
    out = []
    for im in images:
        out.append(pytesseract.image_to_string(_preprocess(im), lang=OCR_LANGS))
    return "\n".join(out)


_ENGINES = [("paddle", _try_paddle), ("easyocr", _try_easyocr), ("tesseract", _try_tesseract)]


def run_ocr(path, engine=None):
    """Пробует движки по порядку (или заданный). Возвращает (текст, имя_движка)."""
    images = _pages(path)
    order = [(n, f) for n, f in _ENGINES if engine in (None, n)]
    errs = []
    for name, fn in order:
        try:
            txt = fn(images)
            if txt and txt.strip():
                return txt, name
        except Exception as e:
            errs.append(f"{name}: {type(e).__name__}")
    raise RuntimeError("Ни один OCR-движок не доступен/не сработал. "
                       "Установите tesseract (+язык pol) или easyocr/paddleocr. Детали: " + "; ".join(errs))


# ─────────────────────────── авто-метаданные (№4) ──────────────────────────
_DATE_PATTERNS = [
    r"pobrani[ae][:\s]*?(\d{4})[.\-/](\d{2})[.\-/](\d{2})",   # Data pobrania: 2026-11-03
    r"pobrani[ae][:\s]*?(\d{2})[.\-/](\d{2})[.\-/](\d{4})",   # 03.11.2026
    r"(\d{4})[.\-/](\d{2})[.\-/](\d{2})",                     # любой ISO как запасной
]
def extract_date(text):
    low = text.lower()
    for pat in _DATE_PATTERNS:
        m = re.search(pat, low)
        if m:
            g = m.groups()
            y, a, b = (g if len(g[0]) == 4 else (g[2], g[1], g[0]))
            try:
                return datetime.date(int(y), int(a), int(b)).isoformat()
            except ValueError:
                continue
    return None

def extract_patient(text):
    """Идентификатор пациента: PESEL (11 цифр) приоритетно, иначе имя после 'Pacjent:'."""
    m = re.search(r"pesel[:\s]*?(\d{11})", text, re.I)
    if m:
        return "pesel_" + m.group(1)
    m = re.search(r"pacjent[:\s]*([A-ZА-ЯŁŚŻŹĆĄĘÓŃ][\w\-]+(?:\s+[A-ZА-ЯŁŚŻŹĆĄĘÓŃ][\w\-]+)*)", text)
    if m:
        return re.sub(r"\s+", "_", m.group(1).strip()).lower()[:48]
    return None


def main(argv=None):
    ap = argparse.ArgumentParser(description="OCR-ингест бланка (офлайн, многоуровневый).")
    ap.add_argument("input", help="PDF или изображение бланка")
    ap.add_argument("--engine", choices=["paddle", "easyocr", "tesseract"], help="форсировать движок")
    ap.add_argument("--meta", action="store_true", help="вывести JSON с датой/пациентом вместо текста")
    ap.add_argument("--run", action="store_true", help="сразу прогнать через interpret.py с авто-датой/пациентом")
    ap.add_argument("--sex", choices=["male", "female"])
    ap.add_argument("--age", type=float)
    args = ap.parse_args(argv)

    text, used = run_ocr(args.input, args.engine)
    date = extract_date(text)
    patient = extract_patient(text)
    sys.stderr.write(f"[ocr_ingest] движок={used}  дата={date}  пациент={patient}\n")

    if args.meta:
        print(json.dumps({"engine": used, "date": date, "patient": patient, "text": text},
                         ensure_ascii=False, indent=2))
        return 0
    if args.run:
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False, encoding="utf-8") as tf:
            tf.write(text); tmp = tf.name
        cmd = [sys.executable, os.path.join(ROOT, "engine", "interpret.py"), tmp]
        if patient: cmd += ["--patient", patient]
        if date: cmd += ["--date", date]
        if args.sex: cmd += ["--sex", args.sex]
        if args.age is not None: cmd += ["--age", str(args.age)]
        try:
            return subprocess.call(cmd)
        finally:
            os.unlink(tmp)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
