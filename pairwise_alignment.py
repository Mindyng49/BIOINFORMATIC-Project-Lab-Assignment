from __future__ import annotations

from dataclasses import dataclass


# Trace-back directions

DIAG, UP, LEFT, STOP = "diag", "up", "left", "stop"

ARROW = {DIAG: "\u2196", UP: "\u2191", LEFT: "\u2190", STOP: "\u2022"}  # ↖ ↑ ← •

GAP_CHAR = "-"

# Refuse absurdly large problems instead of freezing the browser.
MAX_CELLS = 4_000_000



# Input handling

def clean_sequence(raw: str) -> str:
    """Normalise user input: drop FASTA headers, digits, spaces and punctuation."""
    if raw is None:
        return ""
    residues = []
    for line in str(raw).splitlines():
        text = line.strip()
        if not text:
            continue
        if text[0] in ">;":          # FASTA/comment header
            continue
        if "|" in text:              # header pasted without the leading '>'
            continue
        residues.append("".join(ch for ch in text if ch.isascii() and ch.isalpha()))
    return "".join(residues).upper()



# Result container

@dataclass
class AlignmentResult:
    mode: str                      # "global" | "local"
    seq1: str
    seq2: str
    match_score: float
    mismatch_score: float
    gap_score: float
    matrix: list                   # (n+1) x (m+1) integers
    traceback: list                # (n+1) x (m+1) tuple of directions
    score: float                   # final alignment score
    path: list                     # cells from start -> end of the traceback
    aligned1: str
    aligned2: str
    operations: list               # per alignment column, see OPERATIONS below

    # ---- geometry ---------------------------------------------------------
    @property
    def n_rows(self) -> int:
        return len(self.seq1)

    @property
    def n_cols(self) -> int:
        return len(self.seq2)

    @property
    def is_local(self) -> bool:
        return self.mode == "local"

    @property
    def start_cell(self) -> tuple:
        return self.path[0] if self.path else (0, 0)

    @property
    def end_cell(self) -> tuple:
        return self.path[-1] if self.path else (0, 0)

    # ---- statistics -------------------------------------------------------
    @property
    def n_matches(self) -> int:
        return self.operations.count("match")

    @property
    def n_mismatches(self) -> int:
        return self.operations.count("mismatch")

    @property
    def n_gaps(self) -> int:
        return sum(1 for op in self.operations if op.startswith("gap"))

    @property
    def alignment_length(self) -> int:
        return len(self.aligned1)

    @property
    def identity(self) -> float:
        if not self.alignment_length:
            return 0.0
        return 100.0 * self.n_matches / self.alignment_length

    @property
    def path_set(self) -> set:
        return set(self.path)



# The algorithm

def align(
    seq1: str,
    seq2: str,
    match: float = 2,
    mismatch: float = -1,
    gap: float = -2,
    mode: str = "global",
) -> AlignmentResult:
    """Align two sequences with Needleman-Wunsch (global) or Smith-Waterman (local)."""
    a, b = clean_sequence(seq1), clean_sequence(seq2)
    if not a or not b:
        raise ValueError("Both sequences must contain at least one residue.")

    mode = mode.lower()
    if mode not in ("global", "local"):
        raise ValueError("mode must be 'global' or 'local'.")

    n, m = len(a), len(b)
    if (n + 1) * (m + 1) > MAX_CELLS:
        raise ValueError(
            f"Sequences are too long for this demo ({n} x {m} cells). "
            f"Keep the product below {MAX_CELLS:,}."
        )

    local = mode == "local"

    # ---- 1. fill the matrix ------------------------------------------------
    matrix = [[0] * (m + 1) for _ in range(n + 1)]
    traceback = [[()] * (m + 1) for _ in range(n + 1)]

    if not local:
        # first row / column: only gaps are possible
        for i in range(1, n + 1):
            matrix[i][0] = matrix[i - 1][0] + gap
            traceback[i][0] = (UP,)
        for j in range(1, m + 1):
            matrix[0][j] = matrix[0][j - 1] + gap
            traceback[0][j] = (LEFT,)
    else:
        # Smith-Waterman: edges are zero and a cell may restart (STOP)
        for j in range(m + 1):
            traceback[0][j] = (STOP,)
        for i in range(n + 1):
            traceback[i][0] = (STOP,)

    best_value, best_cell = 0, (0, 0)
    lookups = {}

    for i in range(1, n + 1):
        row = matrix[i]
        prev = matrix[i - 1]
        tb_row = traceback[i]
        residue = a[i - 1]
        for j in range(1, m + 1):
            other = b[j - 1]
            sub = match if residue == other else mismatch
            d = prev[j - 1] + sub          # diagonal
            u = prev[j] + gap              # up
            l = row[j - 1] + gap           # left
            value = d if d >= u else u
            if l > value:
                value = l
            if local and value < 0:
                value = 0

            if local and value == 0:
                tb_row[j] = (STOP,)
            else:
                dirs = []
                if d == value:
                    dirs.append(DIAG)
                if u == value:
                    dirs.append(UP)
                if l == value:
                    dirs.append(LEFT)
                tb_row[j] = tuple(dirs) or (STOP,)

            row[j] = value

            if local and value > best_value:
                best_value, best_cell = value, (i, j)

    if local:
        score = best_value
        end = best_cell
    else:
        score = matrix[n][m]
        end = (n, m)

    # ---- 2. trace-back -----------------------------------------------------
    aligned1, aligned2 = [], []
    operations, path = [], []

    i, j = end
    path.append((i, j))
    while True:
        dirs = traceback[i][j]
        if not dirs:
            break
        if local and (dirs[0] == STOP or matrix[i][j] == 0):
            break
        if DIAG in dirs and i > 0 and j > 0:
            move = DIAG
        elif UP in dirs and i > 0:
            move = UP
        elif LEFT in dirs and j > 0:
            move = LEFT
        else:                                    # pragma: no cover - defensive
            break

        if move == DIAG:
            aligned1.append(a[i - 1])
            aligned2.append(b[j - 1])
            operations.append("match" if a[i - 1] == b[j - 1] else "mismatch")
            i -= 1
            j -= 1
        elif move == UP:
            aligned1.append(a[i - 1])
            aligned2.append(GAP_CHAR)
            operations.append("gap_seq2")
            i -= 1
        else:
            aligned1.append(GAP_CHAR)
            aligned2.append(b[j - 1])
            operations.append("gap_seq1")
            j -= 1
        path.append((i, j))

    path.reverse()
    aligned1.reverse()
    aligned2.reverse()
    operations.reverse()

    return AlignmentResult(
        mode=mode,
        seq1=a,
        seq2=b,
        match_score=match,
        mismatch_score=mismatch,
        gap_score=gap,
        matrix=matrix,
        traceback=traceback,
        score=score,
        path=path,
        aligned1="".join(aligned1),
        aligned2="".join(aligned2),
        operations=operations,
    )



# Text helpers

def match_line(aligned1: str, aligned2: str) -> str:
    """'|' identical, '.' mismatch, ' ' gap."""
    out = []
    for x, y in zip(aligned1, aligned2):
        if x == GAP_CHAR or y == GAP_CHAR:
            out.append(" ")
        elif x == y:
            out.append("|")
        else:
            out.append(".")
    return "".join(out)


def alignment_blocks(aligned1: str, aligned2: str, width: int = 60) -> list:
    """Split an alignment into (seq1, markers, seq2) blocks for monospace display."""
    mid = match_line(aligned1, aligned2)
    blocks = []
    for s in range(0, max(len(aligned1), 1), width):
        blocks.append((aligned1[s:s + width], mid[s:s + width], aligned2[s:s + width]))
    return blocks


def as_text(result: AlignmentResult) -> str:
    """Plain-text report, used for the download button."""
    kind = "Global (Needleman-Wunsch)" if not result.is_local else "Local (Smith-Waterman)"
    lines = [
        "Pairwise sequence alignment",
        f"Algorithm : {kind}",
        f"Scoring   : match {result.match_score:+g}  "
        f"mismatch {result.mismatch_score:+g}  gap {result.gap_score:+g}",
        f"Score     : {result.score:g}",
        f"Identity  : {result.identity:.1f}%  "
        f"({result.n_matches}/{result.alignment_length} columns)",
        f"Gaps      : {result.n_gaps}",
        "",
        "Aligned sequences",
        "-----------------",
    ]
    for s1, mid, s2 in alignment_blocks(result.aligned1, result.aligned2):
        lines.extend(["Seq1 " + s1, "     " + mid, "Seq2 " + s2, ""])
    lines.append("Trace-back path (row, col) values")
    lines.append("--------------------------------")
    for k, (i, j) in enumerate(result.path, start=1):
        lines.append(f"{k:>4}  ({i:>4}, {j:>4})  {result.matrix[i][j]}")
    return "\n".join(lines)
