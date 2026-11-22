"""Download a reviewed asset with a hard project disk limit and provenance."""
import argparse
import hashlib
import json
import os
import shutil
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LIMIT = 1_500_000_000_000

def open_response(url):
    """Retry transient HTTP response errors before opening any partial file."""
    delays=(2,5,15)
    for attempt in range(len(delays)+1):
        try:
            return urllib.request.urlopen(url,timeout=90)
        except urllib.error.HTTPError as exc:
            if exc.code not in (429,500,502,503,504) or attempt==len(delays):raise
            print(f'Transient HTTP{exc.code}; retry{attempt+1} in{delays[attempt]}s',flush=True)
            time.sleep(delays[attempt])

def directory_footprint(root):
    """Count every regular file using cached directory-entry metadata."""
    total = 0
    pending = [root]
    while pending:
        directory = pending.pop()
        try:
            with os.scandir(directory) as entries:
                for entry in entries:
                    try:
                        if entry.is_symlink():
                            continue
                        if entry.is_dir(follow_symlinks=False):
                            pending.append(entry.path)
                        elif entry.is_file(follow_symlinks=False):
                            total += entry.stat(follow_symlinks=False).st_size
                    except FileNotFoundError:
                        # Validated cleanup can remove an entry during traversal.
                        continue
        except (FileNotFoundError, PermissionError):
            # An unrelated protected temporary directory can coexist under the
            # workspace. The free-space check below still bounds downloads.
            continue
    return total


def footprint():
    return directory_footprint(ROOT)

def download(url, destination, expected_bytes=None):
    destination = Path(destination).resolve()
    temporary_root = (ROOT / 'data' / 'downloads').resolve()
    if not destination.is_relative_to(temporary_root):
        raise ValueError('Downloads must stay inside project data/downloads')
    destination.parent.mkdir(parents=True, exist_ok=True)
    manifest = destination.with_suffix(destination.suffix + '.manifest.json')
    previous_manifest = json.loads(manifest.read_text()) if manifest.exists() else None
    if destination.exists() and manifest.exists():
        saved = json.loads(manifest.read_text())
        if saved['url'] == url and saved['sha256'] == hashlib.file_digest(destination.open('rb'), 'sha256').hexdigest():
            return destination
        raise ValueError('Existing download fails provenance verification')
    used = footprint()
    free = shutil.disk_usage(ROOT).free
    # Reserve room for compact outputs and keep under the user limit even if a response is unexpectedly large.
    allowance = min(LIMIT - used - 2_000_000_000, free - 2_000_000_000)
    if allowance <= 0 or (expected_bytes and expected_bytes > allowance):
        raise RuntimeError('Disk budget would be exceeded')
    part = destination.with_suffix(destination.suffix + '.part')
    if part.exists() and part.stat().st_size:
        raise RuntimeError('A partial download exists; inspect it before restarting')
    # A zero-byte placeholder contains no source data. The write below can
    # safely reuse it after a connection failed before sending any bytes.
    digest = hashlib.sha256(); size = 0
    with open_response(url) as response, part.open('wb') as out:
        headers = dict(response.headers)
        length = response.headers.get('Content-Length')
        if length and int(length) > allowance:
            raise RuntimeError('Response exceeds disk budget')
        while block := response.read(8 * 1024 * 1024):
            if size + len(block) > allowance:
                raise RuntimeError('Streaming disk budget exceeded')
            out.write(block); digest.update(block); size += len(block)
    if expected_bytes and size != expected_bytes:
        raise ValueError(f'Expected {expected_bytes} bytes, received {size}')
    if previous_manifest and previous_manifest['url'] == url and previous_manifest['sha256'] != digest.hexdigest():
        raise ValueError('Source checksum changed since prior retrieval; preserve partial file for inspection')
    part.replace(destination)
    manifest.write_text(json.dumps(dict(url=url, retrieved_at=datetime.now(timezone.utc).isoformat(),
        bytes=size, sha256=digest.hexdigest(), response_headers=headers), indent=2), encoding='utf-8')
    print(json.dumps(dict(path=str(destination), bytes=size, project_bytes_before=used)), flush=True)
    return destination

if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('url'); parser.add_argument('destination'); parser.add_argument('--expected-bytes', type=int)
    args = parser.parse_args(); download(args.url, args.destination, args.expected_bytes)
