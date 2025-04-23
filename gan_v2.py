"""Room placement in a square with three rectangular corner cut‑outs
------------------------------------------------------------------
Given ≤10 axis‑aligned rectangular rooms (integer sizes) plus an
adjacency list, pack them on an integer grid inside an S×S square with
up to three rectangular corners removed.  Rooms must not overlap and
must lie completely within the boundary.  We maximise the number of
adjacency edges that become shared walls; normally the optimum equals
the full set, so all required adjacencies are met.

We model the problem as a MILP and solve it with OR‑Tools CP‑SAT.
The model uses `NoOverlap2D` for non‑overlap and boolean reified
constraints for adjacency.

This script demonstrates a full workflow: build model → solve → plot.

Usage (example at bottom): simply run the file.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple, Set

import matplotlib.pyplot as plt
import ortools.sat.python.cp_model as cp

# ----------------------------------------------------------------------
#  Data structures
# ----------------------------------------------------------------------
@dataclass
class RoomSpec:
    id: str
    w: int  # default width  (integer)
    h: int  # default height (integer)
    rot: bool = True  # allow 90° rotation


@dataclass
class BoundarySpec:
    S: int  # side length of outer square (grid units)
    chopped: Dict[str, Tuple[int, int]]  # {corner: (cw, ch)}, corners in {"tl","tr","bl","br"}


# ----------------------------------------------------------------------
#  CP‑SAT model builder
# ----------------------------------------------------------------------

def build_model(rooms: List[RoomSpec], edges: Set[Tuple[str, str]], bound: BoundarySpec, maximise=True):
    m = cp.CpModel()
    n = len(rooms)
    # index mapping for convenience
    idx = {r.id: i for i, r in enumerate(rooms)}

    # Decision vars
    x = [m.NewIntVar(0, bound.S, f"x_{r.id}") for r in rooms]  # bottom‑left X
    y = [m.NewIntVar(0, bound.S, f"y_{r.id}") for r in rooms]  # bottom‑left Y

    # rotation selector 0=no rot,1=rot (swap w/h)
    rot = [m.NewBoolVar(f"rot_{r.id}") if r.rot and r.w != r.h else None for r in rooms]

    w = []  # actual width var per room
    h = []  # actual height var per room
    for i, r in enumerate(rooms):
        if rot[i] is None:
            w.append(m.NewConstant(r.w))
            h.append(m.NewConstant(r.h))
        else:
            # w_i = rot * h0 + (1‑rot) * w0
            wi = m.NewIntVar(min(r.w, r.h), max(r.w, r.h), f"w_{r.id}")
            hi = m.NewIntVar(min(r.w, r.h), max(r.w, r.h), f"h_{r.id}")
            m.Add(wi == r.w).OnlyEnforceIf(rot[i].Not())
            m.Add(hi == r.h).OnlyEnforceIf(rot[i].Not())
            m.Add(wi == r.h).OnlyEnforceIf(rot[i])
            m.Add(hi == r.w).OnlyEnforceIf(rot[i])
            w.append(wi)
            h.append(hi)

    # ------------------ boundary constraints ------------------
    for i in range(n):
        m.Add(x[i] + w[i] <= bound.S)
        m.Add(y[i] + h[i] <= bound.S)

    # chopped rectangles exclusion
    # helper: create booleans and Big‑M constraints for each corner cut
    bigM = bound.S
    def exclude_corner(i: int, corner: str, cw: int, ch: int):
        """Ensure room i stays outside rectangular cut at *corner*."""
        bx, by = x[i], y[i]
        wi, hi = w[i], h[i]
        b = m.NewBoolVar(f"outside_{corner}_{i}")
        if corner == "tl":
            # strip x <= cw‑1  intersects?  if in strip then y+hi <= S‑ch
            m.Add(bx <= cw - 1).OnlyEnforceIf(b)
            m.Add(bx >= cw).OnlyEnforceIf(b.Not())
            m.Add(by + hi <= bound.S - ch).OnlyEnforceIf(b)
        elif corner == "tr":
            # x+wi >= S‑cw+1?  if in strip then y+hi <= S‑ch
            m.Add(bx + wi >= bound.S - cw + 1).OnlyEnforceIf(b)
            m.Add(bx + wi <= bound.S - cw).OnlyEnforceIf(b.Not())
            m.Add(by + hi <= bound.S - ch).OnlyEnforceIf(b)
        elif corner == "bl":
            # x <= cw‑1 strip -> by >= ch ? actually bottom left remove rectangle width cw height ch
            m.Add(bx <= cw - 1).OnlyEnforceIf(b)
            m.Add(bx >= cw).OnlyEnforceIf(b.Not())
            m.Add(by >= ch).OnlyEnforceIf(b)
        elif corner == "br":
            m.Add(bx + wi >= bound.S - cw + 1).OnlyEnforceIf(b)
            m.Add(bx + wi <= bound.S - cw).OnlyEnforceIf(b.Not())
            m.Add(by >= ch).OnlyEnforceIf(b)
        else:
            raise ValueError("bad corner")

    for i in range(n):
        for c, (cw, ch) in bound.chopped.items():
            exclude_corner(i, c, cw, ch)

    # ------------------ No overlap ------------------
    intervals_x = []
    intervals_y = []
    for i in range(n):
        # Build explicit end variables because CP‑SAT requires an IntVar, not an expression
        end_x = m.NewIntVar(0, bound.S, f"ex_{i}")
        m.Add(end_x == x[i] + w[i])
        end_y = m.NewIntVar(0, bound.S, f"ey_{i}")
        m.Add(end_y == y[i] + h[i])

        ix = m.NewIntervalVar(x[i], w[i], end_x, f"ix_{i}")
        iy = m.NewIntervalVar(y[i], h[i], end_y, f"iy_{i}")
        intervals_x.append(ix)
        intervals_y.append(iy)
    m.AddNoOverlap2D(intervals_x, intervals_y)

    # ------------------ adjacency variables & constraints ------------------ & constraints ------------------
    adj_bools = {}
    for (a, b) in edges:
        i, j = idx[a], idx[b]
        adj = m.NewBoolVar(f"adj_{a}_{b}")
        adj_bools[(i, j)] = adj

        # Four possible touching orientations encoded with helper bools
        touch = []
        # i right of j
        t1 = m.NewBoolVar(f"t1_{a}_{b}")
        m.Add(x[i] == x[j] + w[j]).OnlyEnforceIf(t1)
        m.Add(y[i] < y[j] + h[j]).OnlyEnforceIf(t1)
        m.Add(y[j] < y[i] + h[i]).OnlyEnforceIf(t1)
        touch.append(t1)
        # j right of i
        t2 = m.NewBoolVar(f"t2_{a}_{b}")
        m.Add(x[j] == x[i] + w[i]).OnlyEnforceIf(t2)
        m.Add(y[i] < y[j] + h[j]).OnlyEnforceIf(t2)
        m.Add(y[j] < y[i] + h[i]).OnlyEnforceIf(t2)
        touch.append(t2)
        # i above j
        t3 = m.NewBoolVar(f"t3_{a}_{b}")
        m.Add(y[i] == y[j] + h[j]).OnlyEnforceIf(t3)
        m.Add(x[i] < x[j] + w[j]).OnlyEnforceIf(t3)
        m.Add(x[j] < x[i] + w[i]).OnlyEnforceIf(t3)
        touch.append(t3)
        # j above i
        t4 = m.NewBoolVar(f"t4_{a}_{b}")
        m.Add(y[j] == y[i] + h[i]).OnlyEnforceIf(t4)
        m.Add(x[i] < x[j] + w[j]).OnlyEnforceIf(t4)
        m.Add(x[j] < x[i] + w[i]).OnlyEnforceIf(t4)
        touch.append(t4)

        # At most one orientation can be true (non‑overlap would forbid >1 anyway)
        m.AddBoolOr(touch + [adj.Not()])  # if adj false, all touch may be false
        # adj iff any touch true
        m.AddBoolOr([t1, t2, t3, t4]).OnlyEnforceIf(adj)
        m.Add(adj == 1).OnlyEnforceIf([t1, t2, t3, t4])  # if any orientation then adj true

    # Objective – maximise satisfied adjacencies
    if maximise and edges:
        m.Maximize(sum(adj_bools.values()))
    return m, (x, y, w, h, rot, adj_bools, idx)


# ----------------------------------------------------------------------
#  Solve & extract layout
# ----------------------------------------------------------------------

def solve_layout(rooms: List[RoomSpec], edges: Set[Tuple[str, str]], bound: BoundarySpec):
    model, aux = build_model(rooms, edges, bound)
    solver = cp.CpSolver()
    solver.parameters.max_time_in_seconds = 30
    solver.parameters.num_search_workers = 8
    res = solver.Solve(model)
    assert res in (cp.OPTIMAL, cp.FEASIBLE), "No placement found!"
    x, y, w, h, rot, _, idx = aux
    placements = {}
    for i, r in enumerate(rooms):
        placements[r.id] = {
            "x": solver.Value(x[i]),
            "y": solver.Value(y[i]),
            "w": solver.Value(w[i]),
            "h": solver.Value(h[i]),
        }
    satisfied = {(a, b) for (a, b), var in aux[5].items() if solver.Value(var) == 1}
    return placements, satisfied


# ----------------------------------------------------------------------
#  Visualisation helper
# ----------------------------------------------------------------------

def plot_layout(placements: Dict[str, Dict[str, int]], edges: Set[Tuple[str, str]], bound: BoundarySpec):
    S = bound.S
    fig, ax = plt.subplots(figsize=(6, 6))
    # draw outer square
    ax.plot([0, S, S, 0, 0], [0, 0, S, S, 0], "k-")
    # mask chopped corners
    for c, (cw, ch) in bound.chopped.items():
        if c == "tl":
            ax.fill([0, cw, cw, 0], [S - ch, S - ch, S, S], color="lightgrey")
        elif c == "tr":
            ax.fill([S - cw, S, S, S - cw], [S - ch, S - ch, S, S], color="lightgrey")
        elif c == "bl":
            ax.fill([0, cw, cw, 0], [0, 0, ch, ch], color="lightgrey")
        elif c == "br":
            ax.fill([S - cw, S, S, S - cw], [0, 0, ch, ch], color="lightgrey")
    ax.set_aspect("equal")
    # grid
    ax.set_xticks(range(S + 1))
    ax.set_yticks(range(S + 1))
    ax.grid(True, which="both", linestyle=":", lw=0.5)

    # draw rooms
    centers = {}
    for rid, rec in placements.items():
        ax.add_patch(plt.Rectangle((rec["x"], rec["y"]), rec["w"], rec["h"], fc="tab:blue", alpha=0.4, ec="k"))
        ax.text(rec["x"] + rec["w"] / 2, rec["y"] + rec["h"] / 2, rid, ha="center", va="center", color="k")
        centers[rid] = (rec["x"] + rec["w"] / 2, rec["y"] + rec["h"] / 2)

    for a, b in edges:
        xa, ya = centers[a]
        xb, yb = centers[b]
        ax.plot([xa, xb], [ya, yb], "r--", lw=0.8)

    ax.set_xlim(-0.5, S + 0.5)
    ax.set_ylim(-0.5, S + 0.5)
    ax.set_xlabel("X (grid)")
    ax.set_ylabel("Y (grid)")
    ax.set_title("Grid‑snapped room layout with chopped corners")
    plt.tight_layout()
    plt.show()

def main():
    main_menu = """1. Enter custom problem
2. Show examples
3. Exit\n>> """

    while (1):
        print(main_menu)
        c = int(input())   
        if (c == 1):
            rooms = [] # blocks
            edges = set() # adjacent pairs
            chopped = {} # dictionary
            print("Enter side length of grid boundary:")
            s = int(input())
            print("Three corners will be cut (except top right) to obtain proposed grid boundary. Enter dimensions of the rectangles to be cut from the 3 corners:")
            print("Top-left:-")
            print("Width:")
            w1 = int(input())
            print("Height:")
            h1 = int(input())
            chopped["tl"] = (w1, h1)

            print("Bottom-left:")
            print("Width:")
            w3 = int(input())
            print("Height:")
            h3 = int(input())
            chopped["bl"] = (w3, h3)

            print("Bottom-right:")
            print("Width:")
            w4 = int(input())
            print("Height:")
            h4 = int(input())
            chopped["br"] = (w4, h4)
            boundary = BoundarySpec(s, chopped)

            print("Enter number of blocks:")
            n = int(input())
            for i in range(n):
                print(f"Enter block {chr(65+i)} details:")
                print("Width:")
                w = int(input())
                print("Height:")
                h = int(input())
                rooms.append(RoomSpec(chr(65+i), w, h, True))

            print("Enter number of edges in adjacency graph:")
            ne = int(input())
            for i in range(ne):
                print(f"Enter edge in space separated format{i+1} (A B):")
                a, b = [x.upper() for x in input().split()]
                edges.add((a,b))
            
            placement, ok_edges = solve_layout(rooms, edges, boundary)
            print("Satisfied adjacencies:", ok_edges)
            plot_layout(placement, edges, boundary)

        elif (c == 2):
            print("Please enter the example number (1-4)\n>> ")
            ex = int(input())
            rooms = [
                    [  
                        RoomSpec("A", 2, 1, True),
                        RoomSpec("B", 3, 1, False),
                        RoomSpec("C", 2, 1, True),
                    ],
                    [
                        RoomSpec("A", 3, 1, True),
                        RoomSpec("B", 3, 2, True),
                        RoomSpec("C", 1, 1, True),
                        RoomSpec("D", 1, 1, True),
                        RoomSpec("E", 2, 1, True),
                    ],
                    [
                        RoomSpec("A", 2, 1, True),
                        RoomSpec("B", 2, 1, True),
                        RoomSpec("C", 2, 11, True),
                        RoomSpec("D", 17, 4, True),
                        RoomSpec("E", 15, 11, True),
                        RoomSpec("F", 2, 9, True),
                    ],
                    [
                        RoomSpec("A", 5, 1, True),
                        RoomSpec("B", 4, 1, True),
                        RoomSpec("C", 1, 1, True),
                        RoomSpec("D", 1, 3, True),
                        RoomSpec("E", 1, 1, True),
                        RoomSpec("F", 1, 1, True),
                        RoomSpec("G", 3, 1, True),
                    ]
                ]
            edges = [
                {("A","B"), ("B","C"),("A","C")},
                {("A","B"), ("B","C"),("A","C"), ("B", "D"), ("B", "E"), ("C", "D")},
                {("A","B"), ("B","C"),("A","C"), ("B", "D"), ("D", "E"), ("C", "D"), ("E","F"),("F","D")},
                {("A","B"), ("B","C"), ("B", "E"), ("B", "D"), ("E", "D"), ("C","F"),("C","D"),("D","F"),("A","G"), ("D","G")}
            ]

            boundaries = [
                BoundarySpec(S=4, chopped={"tl": (1, 1), "tr": (1, 1), "bl": (1, 1)}),
                BoundarySpec(S=4, chopped={"tl": (1, 1), "tr": (1, 1), "bl": (1, 1)}),
                BoundarySpec(S=17, chopped={"tl": (2, 2), "tr": (2, 2), "bl": (2, 2)}),
                BoundarySpec(S=5, chopped={"tl": (2, 2), "tr": (1, 1), "br": (2, 1)})
            ]
            if (ex > 4 or ex < 1):
                print("Invalid input")
            else:
                print('-'*50)
                print(f"Example - {ex}")
                ex -= 1
                for room in rooms[ex]:
                    print(f"{room.id} - ({room.h}X{room.w})") 
                print(f"Edges - {edges[ex]}")
                print(f"Boundary specifications ({boundaries[ex].S}x{boundaries[ex].S}) with {boundaries[ex].chopped} corners removed")

                placement, ok_edges = solve_layout(rooms[ex], edges[ex], boundaries[ex])
                print("Satisfied adjacencies:", ok_edges)
                print('-'*50)
                plot_layout(placement, edges[ex], boundaries[ex])

        elif (c == 3):
            break
        else:
            print("Invalid input please try again")

if __name__ == "__main__":
    main()
    