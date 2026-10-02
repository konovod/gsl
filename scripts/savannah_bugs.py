#!/usr/bin/env python3
"""Inventory the open Savannah bugs of a GNU project and collect their patches.

Fetches the bug list and then each individual bug page, downloads the
attached patch-ish files, checks whether each patch still applies to the
working tree, and writes a TSV inventory.  Nothing in the repository is
modified; the inventory and the downloaded patches go to an output
directory outside the tree.

  * the bug list carries id, summary and submission date for every item in
    a single request, which is all the ordering information we need;
  * an individual bug page carries the full discussion (patches are often
    pasted inline rather than attached) plus the list of attached files;
  * attachments live on file.savannah.gnu.org and are downloadable without
    authentication.

A bug page also lists its field values, which include a Status and a
Category; the summary text is additionally prefixed by the submitter with
words such as "Bug:" or "Feature:".  Both are only hints, but they are
cheap and good enough to sort the inventory.

Usage:
    python scripts/savannah_bugs.py
    python scripts/savannah_bugs.py --group gsl --out /tmp/inv --delay 0.7
    python scripts/savannah_bugs.py --max-bugs 20 --no-download
"""

import argparse
import html
import json
import os
import re
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

USER_AGENT = "gsl-savannah-inventory/1.0 (read-only bug tracker survey)"

LIST_URL = ("https://savannah.gnu.org/bugs/?group={group}&func=&set=open"
            "&msort=0&advsrch=0&morder=date&order=desc&max_rows={rows}")
BUG_URL = "https://savannah.gnu.org/bugs/?{bug_id}"
FILE_URL = "https://file.savannah.gnu.org/file/{name}?file_id={file_id}"

# A list row looks like
#   <tr class="priore"><td><a href="?68704">#68704</a></td>
#   <td><a href="?68704">summary</a></td> ... <td>2026-09-18</td></tr>
ROW_RE = re.compile(
    r"<a href=\"\?(?P<id>\d+)\">#\d+</a></td>\s*"
    r"<td><a href=\"\?\d+\">(?P<summary>.*?)</a></td>.*?"
    r"<td[^>]*>(?P<date>\d{4}-\d{2}-\d{2})</td>", re.S)

ATTACH_RE = re.compile(r"file/(?P<name>[^\"?]+)\?file_id=(?P<file_id>\d+)")

FIELD_RE = re.compile(
    r">(?P<label>[A-Za-z/ ]{2,25}):</span>(?:</a>)?&nbsp;</td>\s*"
    r"<td[^>]*>(?P<value>.*?)</td>", re.S)

# Markers of a patch pasted into the bug text rather than attached.
INLINE_MARKERS = ("```", "diff --git", "--- a/", "@@ ", "+++ b/")

# Suffixes we can do something useful with.  .log and .pdf attachments are
# build output and reports, not patches, so they are listed but not
# downloaded.
PATCH_SUFFIXES = (".patch", ".diff", ".c", ".h", ".rst", ".am", ".in", ".el")
SKIP_SUFFIXES = (".log", ".pdf", ".tar.gz", ".zip", ".gz", ".png", ".jpg")


def fetch(url, delay, retries=3):
    """GET a URL, returning the decoded body, or None if it cannot be had."""
    for attempt in range(retries):
        if delay:
            time.sleep(delay)
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=60) as resp:
                return resp.read().decode("utf-8", "replace")
        except (urllib.error.URLError, OSError, ValueError) as err:
            if attempt == retries - 1:
                print("  ! %s: %s" % (url, err), file=sys.stderr)
                return None
    return None


def clean(text):
    """Strip tags and entities from a chunk of bug page HTML."""
    text = re.sub(r"<[^>]+>", " ", text)
    return " ".join(html.unescape(text).split())


def parse_list(page):
    """Extract (id, date, summary) for every row of a bug list page."""
    rows = []
    for m in ROW_RE.finditer(page):
        rows.append((m.group("id"), m.group("date"),
                     clean(m.group("summary"))))
    return rows


def parse_bug(page):
    """Extract the fields, attachments and patch markers of a bug page."""
    fields = {}
    for m in FIELD_RE.finditer(page):
        label = m.group("label").strip()
        if label not in fields:
            fields[label] = clean(m.group("value"))

    # The attached file links appear twice: once in the "Attached Files"
    # section and once per entry in the history table.
    files = []
    for m in ATTACH_RE.finditer(page):
        entry = (urllib.parse.unquote(m.group("name")), m.group("file_id"))
        if entry not in files:
            files.append(entry)

    # The field table and the navigation come before the messages, and the
    # order of the section anchors is not reliable, so look for the markers
    # in the whole page: none of them occurs outside a pasted patch.
    inline = any(marker in page for marker in INLINE_MARKERS)

    return fields, files, inline


def git_apply_check(path, numstat=False):
    """Run git apply on a patch; return (verdict, detail).

    verdict is "clean", "dirty" (parses but does not apply) or an error
    string.  detail is the +/- line counts or the first line of git's
    complaint.
    """
    args = ["git", "apply", "--numstat" if numstat else "--check", "-p1",
            "--whitespace=nowarn", path]
    proc = subprocess.run(args, cwd=ROOT, capture_output=True, text=True)
    if proc.returncode == 0:
        if not numstat:
            return "clean", ""
        added = removed = files = 0
        for line in proc.stdout.splitlines():
            parts = line.split("\t")
            if len(parts) == 3:
                files += 1
                try:
                    added += int(parts[0])
                    removed += int(parts[1])
                except ValueError:
                    pass
        return "clean", "+%d/-%d in %d file(s)" % (added, removed, files)

    err = proc.stderr.strip().splitlines()
    detail = err[-1] if err else "git apply failed"
    if "does not exist" in detail or "No such file" in detail:
        return "missing", detail
    return "dirty", detail


def download(name, file_id, dest_dir, delay):
    """Download one attachment; return the local path or None."""
    url = FILE_URL.format(name=urllib.parse.quote(name), file_id=file_id)
    body = fetch(url, delay)
    if body is None:
        return None
    safe = re.sub(r"[^\w.-]", "_", name)
    path = os.path.join(dest_dir, safe)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write(body)
    return path


def classify(summary):
    """Guess a kind from the submitter's summary prefix and wording."""
    m = re.match(r"([A-Za-z]+):", summary)
    tag = m.group(1).lower() if m else ""
    if tag in ("bug", "bugfix", "fix", "error", "incorrectness"):
        return "bug"
    if tag in ("doc", "documentation"):
        return "doc"
    if tag in ("feature", "breaking", "test", "tests"):
        return tag
    low = summary.lower()
    if re.search(r"\b(documentation|docs? for|man page)\b", low):
        return "doc"
    if re.search(r"\b(test|testcase|test case)s? (for|of|case)\b", low):
        return "test"
    if re.search(r"\b(inaccurate|incorrect|wrong|bug|uninitial|undefined|nan|"
                 r"does not|fails|failure|overflow|underflow|error|leak)\b", low):
        return "bug"
    if re.search(r"\b(feature|add|implement|support|extend)\b", low):
        return "feature"
    return "-"


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--group", default="gsl",
                        help="Savannah group (default: gsl)")
    parser.add_argument("--out", default=None,
                        help="output directory (default: a temp directory)")
    parser.add_argument("--max-bugs", type=int, default=0,
                        help="only look at the N most recent bugs (0: all)")
    parser.add_argument("--delay", type=float, default=0.5,
                        help="seconds between requests (default: 0.5)")
    parser.add_argument("--no-download", action="store_true",
                        help="do not download or apply-check the patches")
    parser.add_argument("--skip-known", default="FORKNEWS",
                        help="comma separated files whose '#NNNNN' references "
                             "mark bugs that are already handled")
    parser.add_argument("--refresh", action="store_true",
                        help="ignore the cached bug pages")
    args = parser.parse_args()

    # Bug summaries are UTF-8 and the console may not be, so do not let a
    # stray character abort the report at the end of a long run.
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError):
            pass

    out_dir = args.out or os.path.join(tempfile.gettempdir(),
                                       "savannah-%s" % args.group)
    patch_dir = os.path.join(out_dir, "patches")
    os.makedirs(patch_dir, exist_ok=True)
    cache_path = os.path.join(out_dir, "cache.json")

    known = set()
    for name in filter(None, args.skip_known.split(",")):
        path = os.path.join(ROOT, name)
        if os.path.exists(path):
            with open(path, encoding="utf-8", errors="replace") as handle:
                known.update(re.findall(r"#(\d{4,6})", handle.read()))

    cache = {}
    if os.path.exists(cache_path) and not args.refresh:
        with open(cache_path, encoding="utf-8") as handle:
            cache = json.load(handle)

    print("group %s, known (already handled): %s" %
          (args.group, " ".join(sorted(known)) or "none"))

    rows = parse_list(fetch(LIST_URL.format(group=args.group, rows=500),
                            args.delay))
    if args.max_bugs:
        rows = rows[:args.max_bugs]
    print("found %d open items" % len(rows))

    inventory = []
    for n, (bug_id, date, summary) in enumerate(rows, 1):
        entry = cache.get(bug_id)
        if entry is None:
            page = fetch(BUG_URL.format(bug_id=bug_id), args.delay)
            if page is None:
                print("[%3d/%d] #%s skipped" % (n, len(rows), bug_id))
                continue
            fields, files, inline = parse_bug(page)
            entry = {"id": bug_id, "date": date, "summary": summary,
                     "category": fields.get("Category", ""),
                     "status": fields.get("Status", ""),
                     "assigned": fields.get("Assigned to", ""),
                     "inline": inline, "files": files}
            cache[bug_id] = entry
            with open(cache_path, "w", encoding="utf-8") as handle:
                json.dump(cache, handle, indent=1)
            time.sleep(args.delay)

        files = entry["files"]
        patchish = [(name, fid) for name, fid in files
                    if name.lower().endswith(PATCH_SUFFIXES)]
        other = [(name, fid) for name, fid in files
                 if name.lower().endswith(SKIP_SUFFIXES)]

        verdict, detail = ("inline" if entry["inline"] else "none"), ""
        applied = []
        if patchish and not args.no_download:
            bug_dir = os.path.join(patch_dir, bug_id)
            os.makedirs(bug_dir, exist_ok=True)
            for name, fid in patchish:
                local = os.path.join(bug_dir, "%s-%s" % (bug_id,
                                 re.sub(r"[^\w.-]", "_", name)))
                if not os.path.exists(local):
                    local = download(name, fid, bug_dir, args.delay) or ""
                if not local:
                    continue
                if not name.lower().endswith((".patch", ".diff")):
                    # a whole replacement file, not a diff: it tells us the
                    # fix exists but git apply has nothing to say about it
                    applied.append((name, "source",
                                    "%d bytes" % os.path.getsize(local)))
                    continue
                stat, stat_detail = git_apply_check(local, numstat=True)
                applied.append((name, stat, stat_detail))
            if applied:
                verdicts = [a[1] for a in applied]
                if "clean" in verdicts:
                    verdict = "clean" if all(v == "clean" for v in verdicts
                                             ) else "partial"
                elif verdicts[0] in ("source", "dirty"):
                    verdict = verdicts[0]
                detail = "; ".join("%s %s %s" % a for a in applied)

        inventory.append({
            "id": bug_id,
            "date": date,
            "kind": classify(summary),
            "status": entry["status"],
            "category": entry["category"],
            "known": bug_id in known,
            "summary": summary,
            "patches": ",".join(name for name, _ in patchish),
            "other": ",".join(name for name, _ in other),
            "inline": "yes" if entry["inline"] else "",
            "verdict": verdict,
            "detail": detail,
        })

    tsv = os.path.join(out_dir, "inventory.tsv")
    cols = ["id", "date", "kind", "status", "category", "known", "summary",
            "patches", "inline", "verdict", "detail"]
    with open(tsv, "w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(cols) + "\n")
        for row in inventory:
            handle.write("\t".join(
                str(row[c]).replace("\t", " ").replace("\n", " ")
                for c in cols) + "\n")

    with_patches = [r for r in inventory
                    if r["patches"] or r["inline"] == "yes"]
    print("\n%d/%d bugs carry a patch, %d already handled" %
          (len(with_patches), len(inventory),
           sum(1 for r in with_patches if r["known"])))
    print("inventory written to %s" % tsv)
    print("patches written to %s" % patch_dir)

    print("\n%-7s %-11s %-8s %-8s %-5s %s" %
          ("bug", "date", "kind", "applies", "seen", "summary"))
    for row in sorted(with_patches, key=lambda r: r["date"], reverse=True):
        if row["known"]:
            continue
        print("%-7s %-11s %-8s %-8s %-5s %s" %
              ("#" + row["id"], row["date"], row["kind"], row["verdict"],
               "yes" if row["known"] else "", row["summary"][:66]))


if __name__ == "__main__":
    main()