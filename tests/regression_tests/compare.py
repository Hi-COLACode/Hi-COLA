"""Strict artifact readers and baseline comparisons; no implicit regridding."""
import struct
from pathlib import Path
import numpy as np

# These are migration gates for identical inputs/toolchains, not physics accuracy.
TOLERANCES = {
    'tables': (2e-5, 1e-8),
    'force': (1e-8, 1e-11),
    'power': (1e-4, 1e-6),
    'positions': (0., 2e-4),  # Mpc/h; GADGET stores float32 coordinates.
    'velocities': (2e-5, 2e-3),  # km/s
    'lightcone': (2e-5, 2e-3),  # Mpc/h; ASCII output has six significant digits.
}


def numeric(path, columns=None):
    data = np.loadtxt(path, ndmin=2)
    if data.size == 0 or not np.isfinite(data).all():
        raise ValueError(f'Empty or nonfinite numeric output: {path}')
    if columns is not None and data.shape[1] != columns:
        raise ValueError(f'Expected {columns} columns: {path}: {data.shape}')
    return data


def close(actual, reference, kind):
    if actual.shape != reference.shape or actual.size == 0:
        raise ValueError(f'{kind}: incompatible or empty shapes {actual.shape}, {reference.shape}')
    if not np.isfinite(actual).all() or not np.isfinite(reference).all():
        raise ValueError(f'{kind}: nonfinite output')
    rtol, atol = TOLERANCES[kind]
    np.testing.assert_allclose(actual, reference, rtol=rtol, atol=atol, err_msg=kind)
    return dict(max_absolute_error=float(np.max(np.abs(actual-reference))), rtol=rtol, atol=atol)


def record(stream):
    size_bytes = stream.read(4)
    if len(size_bytes) != 4:
        raise ValueError('Truncated GADGET record')
    size = struct.unpack('<I', size_bytes)[0]
    payload = stream.read(size)
    if len(payload) != size or stream.read(4) != size_bytes:
        raise ValueError('Invalid GADGET record markers')
    return payload


def particles(paths):
    """Read GADGET format-1 snapshots and align particles by stable ID."""
    ids, positions, velocities = [], [], []
    for path in paths:
        with Path(path).open('rb') as stream:
            header = record(stream)
            if len(header) != 256:
                raise ValueError('Expected a 256-byte GADGET format-1 header')
            count = sum(struct.unpack_from('<6I', header))
            time, redshift = struct.unpack_from('<2d', header, 72)
            if abs(time-1) > 1e-10 or abs(redshift) > 1e-10:
                raise ValueError('Snapshot did not reach z=0')
            positions.append(np.frombuffer(record(stream), dtype='<f4').reshape(count,3))
            velocities.append(np.frombuffer(record(stream), dtype='<f4').reshape(count,3))
            raw_ids = record(stream)
            if count == 0:
                ids.append(np.array([],dtype=np.uint64))
            elif len(raw_ids) == count*4:
                ids.append(np.frombuffer(raw_ids,dtype='<u4'))
            elif len(raw_ids) == count*8:
                ids.append(np.frombuffer(raw_ids,dtype='<u8'))
            else:
                raise ValueError('Invalid particle ID record')
    ids, positions, velocities = map(np.concatenate,(ids,positions,velocities))
    if len(ids) != 16**3 or len(np.unique(ids)) != len(ids):
        raise ValueError(f'Expected {16**3} unique particles, got {len(ids)}')
    order = np.argsort(ids)
    return ids[order], positions[order], velocities[order]


def force(paths):
    data = np.concatenate([numeric(path,6) for path in paths])
    if len(data) != 16**3:
        raise ValueError('Incomplete force grid')
    order = np.lexsort(data[:,:3].T[::-1])
    data = data[order]
    if len(np.unique(data[:,:3],axis=0)) != len(data):
        raise ValueError('Duplicate force grid coordinates')
    return data


def lightcone(paths):
    # Empty rank/step files are legitimate; the complete cone must contain crossings.
    rows = []
    for path in paths:
        if any(line.strip() and not line.lstrip().startswith('#') for line in Path(path).read_text().splitlines()):
            rows.append(numeric(path,3))
    if not rows:
        raise ValueError('No lightcone particle crossings were output')
    data = np.concatenate(rows)
    return data[np.lexsort(data.T[::-1])]
