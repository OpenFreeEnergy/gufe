import zstandard as zstd
from zstandard import ZstdError

def zst_compress(data: bytes):
    compressor = zstd.ZstdCompressor()
    return compressor.compress(data)

def zst_decompress(data: bytes):
    try:
        decompressor = zstd.ZstdDecompressor()
        return decompressor.decompress(data)
    # need to ensure backwards compatibility for noncompressed artifacts until gufe 2.0
    except ZstdError:
        return data
