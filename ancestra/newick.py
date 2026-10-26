"""Minimal, dependency-free Newick tree utilities.

Only what the ``ancestra`` package needs: a light :class:`Node`, an iterative
parser (no recursion limit for deep gene trees), a writer, and
:class:`TreeIndex` for O(1) ancestor tests and fast MRCA queries on species
trees.

Supported Newick subset: names (unquoted, or single-quoted), branch lengths,
internal node names, and bracketed comments (e.g. NHX ``[&&NHX:...]``) which
are skipped.
"""
from __future__ import annotations

from typing import Iterable, Iterator

__all__ = ["Node", "parse_newick", "read_newick", "TreeIndex"]


class Node:
    """A rooted tree node."""

    __slots__ = ("name", "length", "children", "parent", "data")

    def __init__(self, name: str = "", length: float | None = None):
        self.name = name
        self.length = length
        self.children: list[Node] = []
        self.parent: Node | None = None
        self.data: dict = {}

    def add_child(self, child: "Node") -> "Node":
        child.parent = self
        self.children.append(child)
        return child

    def is_leaf(self) -> bool:
        return not self.children

    def preorder(self) -> Iterator["Node"]:
        stack = [self]
        while stack:
            n = stack.pop()
            yield n
            stack.extend(reversed(n.children))

    def postorder(self) -> Iterator["Node"]:
        out = list(self.preorder())
        return reversed(out)

    def leaves(self) -> list["Node"]:
        return [n for n in self.preorder() if not n.children]

    def leaf_names(self) -> list[str]:
        return [n.name for n in self.leaves()]

    def find(self, name: str) -> "Node | None":
        for n in self.preorder():
            if n.name == name:
                return n
        return None

    def to_newick(self, lengths: bool = True, internal_names: bool = True,
                  fmt: str = "{:.6g}") -> str:
        """Serialise to a Newick string (with trailing ';')."""
        parts: list[str] = []
        # iterative serialisation: (node, state) where state 0=enter, 1=exit
        stack: list[tuple[Node, int]] = [(self, 0)]
        while stack:
            n, st = stack.pop()
            if st == 0:
                if n.children:
                    parts.append("(")
                    stack.append((n, 1))
                    for i, c in enumerate(reversed(n.children)):
                        stack.append((c, 0))
                        if i < len(n.children) - 1:
                            stack.append((None, 2))  # type: ignore[arg-type]
                else:
                    parts.append(_label(n, lengths, True, fmt))
            elif st == 1:
                parts.append(")")
                parts.append(_label(n, lengths, internal_names, fmt))
            else:
                parts.append(",")
        return "".join(parts) + ";"

    def __repr__(self) -> str:  # pragma: no cover
        return f"Node({self.name!r}, n_children={len(self.children)})"


def _label(n: Node, lengths: bool, names: bool, fmt: str) -> str:
    s = _quote(n.name) if (names and n.name) else ""
    if lengths and n.length is not None:
        s += ":" + fmt.format(n.length)
    return s


_SPECIAL = set("(),:;[]' \t\n")


def _quote(name: str) -> str:
    if any(c in _SPECIAL for c in name):
        return "'" + name.replace("'", "''") + "'"
    return name


def parse_newick(s: str) -> Node:
    """Parse a Newick string into a :class:`Node` tree (iterative)."""
    s = s.strip()
    if not s.endswith(";"):
        s += ";"
    root = Node()
    cur = root
    i, n = 0, len(s)
    while i < n:
        c = s[i]
        if c == "(":
            cur = cur.add_child(Node())
            i += 1
        elif c == ",":
            if cur.parent is None:
                raise ValueError(f"Malformed Newick at position {i}")
            cur = cur.parent.add_child(Node())
            i += 1
        elif c == ")":
            if cur.parent is None:
                raise ValueError(f"Unbalanced ')' at position {i}")
            cur = cur.parent
            i += 1
        elif c == ";":
            break
        elif c == "[":
            j = s.find("]", i)
            if j < 0:
                raise ValueError("Unterminated comment")
            i = j + 1
        elif c == ":":
            j = i + 1
            while j < n and s[j] not in ",);[":
                j += 1
            txt = s[i + 1:j].strip()
            cur.length = float(txt) if txt else None
            i = j
        elif c.isspace():
            i += 1
        else:
            if c == "'":
                j = i + 1
                buf = []
                while j < n:
                    if s[j] == "'":
                        if j + 1 < n and s[j + 1] == "'":
                            buf.append("'")
                            j += 2
                            continue
                        break
                    buf.append(s[j])
                    j += 1
                cur.name = "".join(buf)
                i = j + 1
            else:
                j = i
                while j < n and s[j] not in "(),:;[":
                    j += 1
                cur.name = s[i:j].strip()
                i = j
    return root


def read_newick(path: str) -> Node:
    with open(path) as fh:
        return parse_newick(fh.read())


class TreeIndex:
    """Index on a (species) tree for fast clade queries.

    * ``is_ancestor(a, b)``: True if ``a`` is an ancestor of, or equal to, ``b``.
    * ``mrca(names)``: name of the most recent common ancestor of a set of node
      names (leaves or internal).
    * ``clade(name)``: set of leaf names below a node.
    """

    def __init__(self, root: Node):
        self.root = root
        self.nodes: dict[str, Node] = {}
        self.pre: dict[str, int] = {}
        self.post_end: dict[str, int] = {}
        self.depth: dict[str, int] = {}
        order = list(root.preorder())
        for k, nd in enumerate(order):
            if not nd.name:
                raise ValueError("TreeIndex requires all nodes to be named")
            if nd.name in self.nodes:
                raise ValueError(f"Duplicate node name {nd.name!r}")
            self.nodes[nd.name] = nd
            self.pre[nd.name] = k
            self.depth[nd.name] = 0 if nd.parent is None else self.depth[nd.parent.name] + 1
        for nd in root.postorder():
            if nd.children:
                self.post_end[nd.name] = max(self.post_end[c.name] for c in nd.children)
            else:
                self.post_end[nd.name] = self.pre[nd.name]
        self._clade_cache: dict[str, frozenset] = {}
        self.leaves = [n.name for n in order if not n.children]
        self.internal = [n.name for n in order if n.children]

    def __contains__(self, name: str) -> bool:
        return name in self.nodes

    def is_ancestor(self, a: str, b: str) -> bool:
        """a is an ancestor of b or a == b."""
        pa = self.pre[a]
        return pa <= self.pre[b] <= self.post_end[a]

    def is_strict_descendant(self, a: str, b: str) -> bool:
        """a is a strict descendant of b."""
        return a != b and self.is_ancestor(b, a)

    def clade(self, name: str) -> frozenset:
        c = self._clade_cache.get(name)
        if c is None:
            c = frozenset(n.name for n in self.nodes[name].leaves())
            self._clade_cache[name] = c
        return c

    def mrca(self, names: Iterable[str]) -> str:
        it = iter(set(names))
        try:
            first = next(it)
        except StopIteration:
            raise ValueError("mrca of empty set") from None
        lo, hi = self.pre[first], self.post_end[first]
        for nm in it:
            p = self.pre[nm]
            e = self.post_end[nm]
            lo = min(lo, p)
            hi = max(hi, e)
        node = self.nodes[first]
        while not (self.pre[node.name] <= lo and self.post_end[node.name] >= hi):
            node = node.parent
            assert node is not None
        return node.name
