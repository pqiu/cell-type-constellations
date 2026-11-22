"""Stream XDR RDS compression and support R's documented long-vector lengths."""
from contextlib import ExitStack
import gzip
import enum
import struct
import numpy as np
from rdata.parser._xdr import ParserXDR
from rdata.parser._parser import DEFAULT_ALTREP_MAP, RObject, RObjectInfo, RObjectType, BYTECODE_SPECIAL_SET


class SentinelType(enum.IntEnum):
    MISSINGARG = 251
    UNBOUNDVALUE = 252


class PushbackStream:
    def __init__(self, stream):
        self.stream = stream
        self.pending = b''

    def read(self, size):
        if size < 0:
            raise ValueError('Unbounded reads are not supported')
        prefix = self.pending[:size]
        self.pending = self.pending[len(prefix):]
        return prefix + self.stream.read(size-len(prefix))

    def unread(self, data):
        self.pending = data + self.pending

    def tell(self):
        return self.stream.tell()-len(self.pending)


class LongVectorParser(ParserXDR):
    def parse_R_object(self, reference_list=None, bytecode_rep_list=None, info_int=None):
        if info_int is not None and info_int in {kind.value for kind in BYTECODE_SPECIAL_SET}:
            return super().parse_R_object(reference_list, bytecode_rep_list, info_int)
        # R's ReadBCLang default ignores its non-special type marker and calls
        # ReadItem on the next serialized flags. ReadItem is ordinary recursive
        # parsing, so it must not propagate bytecode CAR/CDR marker semantics.
        # https://github.com/wch/r-source/blob/trunk/src/main/serialize.c
        flag = self.parse_int()
        if flag not in (251, 252):
            self.file.unread(struct.pack('>i', flag))
            return super().parse_R_object(reference_list, bytecode_rep_list, None)
        sentinel = SentinelType(flag)
        # Preserve distinct unsupported-conversion sentinels. A requested counts
        # or metadata slot containing one fails conversion rather than becoming NULL.
        return RObject(RObjectInfo(sentinel, False, False, False, 0, 0), sentinel.name, None)

    def _parse_array(self, dtype):
        length = self.parse_int()
        if length == -1:
            high = self.parse_int() & 0xffffffff
            low = self.parse_int() & 0xffffffff
            length = (high << 32) | low
            if length < 2**31 or length > 2**52:
                raise ValueError('Invalid R long-vector length')
        elif length < 0:
            raise ValueError('Invalid R array length')
        return self._parse_array_values(dtype, length)

    def _parse_array_values(self, dtype, length):
        dtype = np.dtype(dtype)
        size = length * dtype.itemsize
        buffer = self.file.read(size)
        if len(buffer) != size:
            raise ValueError('Truncated R array')
        return np.frombuffer(buffer, dtype=dtype.newbyteorder('>')).astype(dtype, copy=True)


def parse_file(path, compression_layers):
    if compression_layers not in (0, 1, 2):
        raise ValueError('Explicit reviewed compression depth required')
    with ExitStack() as stack:
        stream = stack.enter_context(open(path, 'rb'))
        for _ in range(compression_layers):
            stream = stack.enter_context(gzip.GzipFile(fileobj=stream, mode='rb'))
        if stream.read(2) != b'X\n':
            raise ValueError('Reviewed XDR RDS required')
        parser = LongVectorParser(memoryview(b''), expand_altrep=True, altrep_constructor_dict=DEFAULT_ALTREP_MAP)
        parser.file = PushbackStream(stream)
        result = parser.parse_all()
        parser.check_complete()
        return result
