import numpy as np

from storage import (
    iter_legacy_pair_chunks,
    packed_symmetric_get,
    packed_upper_index,
    packed_upper_size,
    pair_count,
    unpack_upper_triangle,
)


def test_packed_index_and_unpack():
    n = 5
    packed = np.arange(packed_upper_size(n), dtype=np.float64)
    full = unpack_upper_triangle(packed, n)
    for i in range(n):
        for j in range(n):
            expected = packed[packed_upper_index(i, j, n)]
            assert full[i, j] == expected
            assert packed_symmetric_get(packed, i, j, n) == expected


def test_legacy_pair_chunks_preserve_old_order():
    n = 7
    observed = []
    packed_indices = []
    starts = []
    stops = []
    for chunk in iter_legacy_pair_chunks(n, chunk_size=4):
        starts.append(chunk.start)
        stops.append(chunk.stop)
        observed.extend(zip(chunk.row_i.tolist(), chunk.row_j.tolist()))
        packed_indices.extend(chunk.packed_indices.tolist())

    expected = [(i, j) for i in range(n) for j in range(i)]
    assert observed == expected
    assert len(observed) == pair_count(n)
    assert starts[0] == 0
    assert stops[-1] == pair_count(n)
    for (i, j), packed_index in zip(observed, packed_indices):
        assert packed_index == packed_upper_index(i, j, n)


def test_zero_or_one_row_has_no_pair_chunks():
    assert list(iter_legacy_pair_chunks(0, 10)) == []
    assert list(iter_legacy_pair_chunks(1, 10)) == []
