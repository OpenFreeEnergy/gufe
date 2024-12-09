import zstandard as zstd


def zst_compress(data: bytes):
    compressor = zstd.ZstdCompressor()
    return compressor.compress(data)


def zst_decompress(data: bytes):
    # need to include backwards compatibility for noncompressed artifacts until gufe 2.0
    try:
        decompressor = zstd.ZstdDeCompressor()
        return decompressor.decompress(data)
    # TODO: find out what exception is thrown if we try to decompress data that isn't compressed
    except ...:
        return data
