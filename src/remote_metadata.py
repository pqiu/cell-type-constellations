"""Inspect pinned H5AD specimen metadata via verified HTTP byte ranges."""
import io
import urllib.request
from collections import OrderedDict

class RangeReader(io.RawIOBase):
    def __init__(self, url, size, block_size=1024 * 1024, cache_blocks=64):
        self.url = url; self.size = int(size); self.position = 0
        self.block_size = block_size; self.cache_blocks = cache_blocks
        self.cache = OrderedDict(); self.etag = None; self.bytes_retrieved = 0

    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.position

    def seek(self, offset, whence=0):
        new = offset if whence == 0 else (self.position + offset if whence == 1 else self.size + offset)
        if new < 0: raise ValueError('Negative seek')
        self.position = new; return new

    def block(self, index):
        if index in self.cache:
            self.cache.move_to_end(index); return self.cache[index]
        start = index * self.block_size; end = min(self.size, start + self.block_size) - 1
        request = urllib.request.Request(self.url, headers={'Range': f'bytes={start}-{end}'})
        with urllib.request.urlopen(request, timeout=90) as response:
            if response.status != 206 or response.headers.get('Content-Range') != f'bytes {start}-{end}/{self.size}':
                raise ValueError('Server did not honor the exact byte range')
            etag = response.headers.get('ETag')
            if self.etag is not None and etag != self.etag: raise ValueError('Source changed during inspection')
            self.etag = etag
            data = response.read(end - start + 2)
            if len(data) != end - start + 1: raise ValueError('Unexpected range length')
        self.bytes_retrieved += len(data); self.cache[index] = data
        if len(self.cache) > self.cache_blocks: self.cache.popitem(last=False)
        return data

    def read(self, size=-1):
        if size < 0: raise ValueError('Whole-file reads are disabled for metadata inspection')
        remaining = min(size, self.size - self.position); blocks = []
        while remaining > 0:
            index, offset = divmod(self.position, self.block_size)
            block = self.block(index); take = min(remaining, len(block) - offset)
            blocks.append(block[offset:offset + take]); self.position += take; remaining -= take
        return b''.join(blocks)

    def readinto(self, buffer):
        data = self.read(len(buffer)); buffer[:len(data)] = data; return len(data)
