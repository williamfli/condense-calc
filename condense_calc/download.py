"""Download a URL into the current folder with several connections, then check its md5 (used by get_reference.sh).

    python download.py URL MD5 [CONNECTIONS]      (or python -m condense_calc.download ...)

Some servers give little per connection (~0.35 MB/s from encode-public on one cluster), so the file is fetched
as 32 MB HTTP ranges, CONNECTIONS (default 16) at a time, and written in order to <name>.part. An interrupted
download resumes from the whole ranges already in <name>.part. On success the file is renamed to <name>.
Exit code: 0 ok; 2 the md5 does not match (the file is deleted); 1 another failure (run again to resume).
Only the Python standard library: the same on Linux and macOS.
"""
import hashlib
import os
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

CHUNK = 32 << 20     # bytes per range request
AHEAD = 2            # ranges in flight per connection (memory: CONNECTIONS * AHEAD * CHUNK = 1 GB at 16)
TRIES = 8
FATAL_HTTP = {400, 401, 403, 404, 405, 410, 416}   # retrying will not help
AGENT = {'User-Agent': 'condense-calc'}


class Fatal(Exception):
    pass


class Mismatch(Exception):
    pass


def _open(url, headers, tries=TRIES):
    """urlopen with retries on network errors; Fatal for errors that retrying does not fix."""
    for attempt in range(1, tries + 1):
        try:
            return urllib.request.urlopen(urllib.request.Request(url, headers={**AGENT, **headers}), timeout=300)
        except urllib.error.HTTPError as e:
            if e.code in FATAL_HTTP:
                raise Fatal(f'HTTP {e.code} {e.reason}: {url}')
            err = e
        except Exception as e:
            if 'CERTIFICATE_VERIFY_FAILED' in repr(e):
                raise Fatal(f'{e!r}: check SSL_CERT_FILE / SSL_CERT_DIR')
            err = e
        print(f'download: {url}: {err!r}; retry {attempt}', file=sys.stderr, flush=True)
        time.sleep(min(300, 10 * attempt ** 2))
    raise IOError(f'{url}: failed {tries} times')


def remote_size(url):
    """(bytes, True if the server answers range requests)."""
    with _open(url, {'Range': 'bytes=0-0'}) as r:
        if r.status == 206 and '/' in r.headers.get('Content-Range', ''):
            total = r.headers['Content-Range'].rsplit('/', 1)[1]
            if total.isdigit():
                return int(total), True
        return int(r.headers.get('Content-Length') or -1), False


def get_range(url, a, b):
    for attempt in range(1, TRIES + 1):
        try:
            with _open(url, {'Range': f'bytes={a}-{b}'}) as r:
                if r.status != 206:
                    raise Fatal(f'the server ignored a range request (HTTP {r.status})')
                data = r.read(b - a + 2)
            if len(data) == b - a + 1:
                return data
            err = IOError(f'{len(data)} bytes for a range of {b - a + 1}')
        except Fatal:
            raise
        except Exception as e:   # a connection lost while reading
            err = e
        print(f'download: bytes {a}-{b}: {err!r}; retry {attempt}', file=sys.stderr, flush=True)
        time.sleep(min(300, 10 * attempt ** 2))
    raise IOError(f'bytes {a}-{b} failed {TRIES} times')


def md5_of(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(8 << 20), b''):
            h.update(block)
    return h.hexdigest()


def download(url, md5, connections=16, folder='.'):
    out = Path(folder) / url.rstrip('/').rsplit('/', 1)[-1]
    part = out.with_name(out.name + '.part')
    if out.exists():
        if md5_of(out) == md5:
            return out
        out.unlink()
    size, ranges = remote_size(url)
    h = hashlib.md5()
    if not ranges:   # one connection, from the start
        part.unlink(missing_ok=True)
        with _open(url, {}) as r, open(part, 'wb') as f:
            for block in iter(lambda: r.read(8 << 20), b''):
                f.write(block)
                h.update(block)
    else:
        done = 0
        if part.exists():   # resume: keep the whole ranges already written
            done = min(part.stat().st_size, size) // CHUNK * CHUNK
            with open(part, 'r+b') as f:
                f.truncate(done)
                for _ in range(done // CHUNK):
                    h.update(f.read(CHUNK))
            if done:
                print(f'download: {out.name}: resuming at {done / 1e9:.1f} of {size / 1e9:.1f} GB', flush=True)
        starts = list(range(done, size, CHUNK))
        pool = ThreadPoolExecutor(connections)
        try:
            with open(part, 'ab') as f:
                pending, nxt = {}, 0
                for k in range(len(starts)):
                    while nxt < len(starts) and nxt < k + connections * AHEAD:
                        a = starts[nxt]
                        pending[nxt] = pool.submit(get_range, url, a, min(a + CHUNK, size) - 1)
                        nxt += 1
                    data = pending.pop(k).result()
                    f.write(data)
                    h.update(data)
        finally:
            pool.shutdown(wait=False, cancel_futures=True)
        if part.stat().st_size != size:
            raise IOError(f'{part.name} has {part.stat().st_size} bytes, expected {size}')
    if h.hexdigest() != md5:
        part.unlink()
        raise Mismatch(f'md5 mismatch, deleted: {out.name} (md5 {h.hexdigest()}, expected {md5})')
    part.rename(out)
    return out


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) not in (2, 3):
        print(__doc__, file=sys.stderr)
        return 1
    try:
        download(argv[0], argv[1], int(argv[2]) if len(argv) == 3 else 16)
    except Mismatch as e:
        print(f'download: {e}', file=sys.stderr)
        return 2
    except Exception as e:
        print(f'download: {argv[0]}: {e}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    rc = main()
    sys.stdout.flush(); sys.stderr.flush()
    os._exit(rc)   # do not wait for ranges still being retried in other threads
