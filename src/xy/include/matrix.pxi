cpdef mat_repr(mat):
    data = ', '.join([repr(list(x)) for x in mat])
    return f'matrix({data})'


cpdef mat_pretty(mat, fmt_str='%.3f'):
    data = [fmt_number(x, fmt_str) for x in mat.flat]
    cols = [data[i::mat.nrows] for i in range(mat.ncols)]
    sizes = [max(map(len, col)) for col in cols]
    cols = [[x.rjust(k) for x in col] for (col, k) in zip(cols, sizes)]
    lines = ['|%s|' % '  '.join(line) for line in zip(*cols)]
    return '\n'.join(lines)


cpdef affine_pretty(mat, vec, fmt_str='%.3f'):
    data = [fmt_number(x, fmt_str) for x in mat.flat]
    cols = [data[i::mat.nrows] for i in range(mat.ncols)]
    cols.append([fmt_number(x, fmt_str) for x in vec])
    sizes = [max(map(len, col)) for col in cols]
    cols = [[x.rjust(k) for x in col] for (col, k) in zip(cols, sizes)]
    return '\n'.join(affine_lines(cols))


def affine_lines(cols):
    for line in zip(*cols):
        *rest, last = line
        data = ' '.join(rest)
        yield f'|{data} : {last}|'


cdef fmt_number(x, fmt_str):
    if isinstance(x, float):
        return (fmt_str % x).rstrip('0').rstrip('.')
    else:
        return repr(x)
