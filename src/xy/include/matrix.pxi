cdef mat_repr(mat):
    data = ', '.join([repr(list(x)) for x in mat])
    return f'matrix({data})'


cdef mat_pretty(mat, fmt_str='%.3f'):
    data = [fmt(x, fmt_str=fmt_str) for x in mat.flat]
    cols = [data[i::mat.nrows] for i in range(mat.ncols)]
    sizes = [max(map(len, col)) for col in cols]
    cols = [[x.rjust(k) for x in col] for (col, k) in zip(cols, sizes)]
    lines = ['|%s|' % '  '.join(line) for line in zip(*cols)]
    return '\n'.join(lines)

cdef fmt(x, fmt_str):
    if isinstance(x, float):
        return (fmt_str % x).rstrip('0').rstrip('.')
    else:
        return repr(x)


