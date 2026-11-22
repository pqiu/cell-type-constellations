"""Recover only provable LF-to-CRLF expansion of an otherwise valid gzip stream.

Original downloads are never modified. Recovery requires an exact forward byte
roundtrip, gzip decompression with CRC/length validation, and independent trailer
CRC32/ISIZE agreement with the complete expanded payload. No CRC errors are ignored.
"""
import gzip
import hashlib
import struct
import zlib


def read_verified_gzip(path):
    raw=path.read_bytes()
    try:
        expanded=gzip.decompress(raw)
        candidate=raw
        mode='original_gzip'
        raw_valid=True
    except (gzip.BadGzipFile,EOFError,zlib.error):
        assert raw.count(b'\n')>0 and raw.count(b'\n')==raw.count(b'\r\n'), 'Not a complete LF-to-CRLF expansion'
        candidate=raw.replace(b'\r\n',b'\n')
        assert len(candidate)<len(raw) and candidate.replace(b'\n',b'\r\n')==raw
        expanded=gzip.decompress(candidate)
        mode='reverse_complete_LF_to_CRLF_expansion'
        raw_valid=False
    # This study supplies single-member gzip objects. Reject other layouts rather
    # than infer what a source CRC was intended to cover.
    crc32,isize=struct.unpack('<II',candidate[-8:])
    assert zlib.crc32(expanded)&0xffffffff==crc32
    assert len(expanded)&0xffffffff==isize
    return expanded,dict(raw_bytes=len(raw),raw_sha256=hashlib.sha256(raw).hexdigest(),raw_gzip_valid=raw_valid,
        transport_operation=mode,recovered_gzip_bytes=len(candidate),recovered_gzip_sha256=hashlib.sha256(candidate).hexdigest(),
        exact_forward_byte_roundtrip=True,expanded_bytes=len(expanded),expanded_sha256=hashlib.sha256(expanded).hexdigest(),
        source_trailer_CRC32=crc32,source_trailer_ISIZE=isize,full_payload_CRC32_and_ISIZE_verified=True)
