"""Minimal S-expression reader for KiCad files (.kicad_pcb, .kicad_sch).

Nodes are plain Python lists: ["footprint", "Lib:Name", ["layer", "F.Cu"], ...].
Atoms are str; quoted strings keep no quotes. Numbers are left as strings and
converted on demand with num().
"""

import re

_TOKEN = re.compile(r'\s*(?:(\()|(\))|"((?:[^"\\]|\\.)*)"|([^\s()"]+))', re.S)


def parse(text):
    stack = [[]]
    pos = 0
    n = len(text)
    while pos < n:
        m = _TOKEN.match(text, pos)
        if not m:
            if text[pos:].strip() == "":
                break
            raise ValueError(f"S-expression parse error at offset {pos}")
        pos = m.end()
        if m.group(1):
            stack.append([])
        elif m.group(2):
            node = stack.pop()
            stack[-1].append(node)
        elif m.group(3) is not None:
            stack[-1].append(m.group(3).replace('\\"', '"').replace("\\\\", "\\"))
        elif m.group(4) is not None:
            stack[-1].append(m.group(4))
    if len(stack) != 1:
        raise ValueError("Unbalanced parentheses")
    return stack[0][0] if len(stack[0]) == 1 else stack[0]


def load(path):
    with open(path, encoding="utf-8") as f:
        return parse(f.read())


def num(x):
    return float(x)


def children(node, tag):
    """All direct child lists of `node` whose head is `tag`."""
    return [c for c in node[1:] if isinstance(c, list) and c and c[0] == tag]


def child(node, tag, default=None):
    for c in node[1:]:
        if isinstance(c, list) and c and c[0] == tag:
            return c
    return default


def value(node, tag, default=None):
    """First atom of child `tag`, e.g. value(fp, 'layer') -> 'F.Cu'."""
    c = child(node, tag)
    if c is None or len(c) < 2:
        return default
    return c[1]


def xy(node, tag):
    c = child(node, tag)
    if c is None:
        return None
    return (num(c[1]), num(c[2]))
