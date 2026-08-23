"""
sj_session.py -- an interactive Scientific Jigsaw session.

Realises the four persistent actions of the paper (TEST, WHY, COMPARE, UPDATE)
as a command line, so a user can interrogate an assembly model without writing
Python. Standard library only (uses cmd, which gives history and tab
completion).

    python3 sj_session.py                # starts with the eIF3 encoding
    python3 sj_session.py proteasome     # or another built-in board

Typical session:

    load eif3
    ask exact Tif35 Tif34         -> impossible, with the reason named
    mode merger                   -> switch representation
    ask exact Tif35 Tif34         -> 5/14 = 0.357
    see before Tif34 Tif35        -> condition; space narrows
    walk Nip1                     -> stand in the tree and look around
    do drop-p Tif34 Tif35         -> intervene on the encoding
    compare                       -> counterfactual report vs the factual model
    back                          -> undo the last step

Every answer that is a refusal names the relation or the representational limit
responsible; nothing is returned as a bare score.
"""

from __future__ import annotations

import cmd
import shlex
import sys
from typing import Dict, List, Optional, Tuple

from sj_prob import TypedModel, merger_F, merger_support, drop_relations
from sj_navigate import (Navigator, before, last, at_step, AND, OR, NOT)

# --------------------------------------------------------------- boards
def _fs(pairs):
    return frozenset(frozenset(p) for p in pairs)


BOARDS: Dict[str, Dict] = {
    "eif3": {
        "desc": "yeast eIF3 core",
        "V": ("Tif32", "Prt1", "Nip1", "Tif35", "Tif34"),
        "C": _fs([("Prt1", "Tif34"), ("Tif32", "Prt1"),
                  ("Tif35", "Tif34"), ("Tif32", "Nip1")]),
        "P": (("Tif34", "Tif35"),),
        "seed": "Prt1",
        "P_seeded": (("Prt1", "Tif34"), ("Prt1", "Tif32"),
                     ("Tif34", "Tif35"), ("Tif32", "Nip1")),
        "queries": {"g:i": ("Tif35", "Tif34"), "b:i": ("Prt1", "Tif34"),
                    "b:g:i": ("Prt1", "Tif35", "Tif34")},
    },
    "eif3-af3": {
        "desc": "eIF3 with AlphaFold 3 recurrent contacts (10/10)",
        "V": ("Tif32", "Prt1", "Nip1", "Tif35", "Tif34"),
        "C": _fs([("Tif32", "Prt1"), ("Tif32", "Nip1"), ("Tif35", "Tif34"),
                  ("Prt1", "Tif34"), ("Prt1", "Tif35"), ("Prt1", "Nip1"),
                  ("Tif32", "Tif35")]),
        "P": (("Tif34", "Tif35"),),
        "seed": "Prt1",
        "queries": {"g:i": ("Tif35", "Tif34")},
    },
    "proteasome": {
        "desc": "20S late stage; b1/b5/b6 onto the 13S precursor",
        "V": ("13S", "b1", "b5", "b6"),
        "C": _fs([("13S", "b1"), ("13S", "b5"), ("13S", "b6")]),
        "P": (),
        "seed": "13S",
        "P_seeded": (),
        "queries": {"13S+b1": ("13S", "b1"), "13S+b5+b6": ("13S", "b5", "b6"),
                    "13S+b1+b5": ("13S", "b1", "b5")},
    },
    "mtor": {
        "desc": "mTORC1 with competing FRB occupants (rapamycin example)",
        "V": ("mTOR", "mLST8", "RAPTOR", "S6K1", "FKBP12rap"),
        "C": _fs([("mTOR", "mLST8"), ("mTOR", "RAPTOR"),
                  ("mTOR", "S6K1"), ("mTOR", "FKBP12rap")]),
        "P": (),
        "X": _fs([("S6K1", "FKBP12rap")]),
        "seed": "mTOR",
        "queries": {"both-FRB": ("mTOR", "S6K1", "FKBP12rap")},
    },
}


# ------------------------------------------------------------- the session
class Session(cmd.Cmd):
    intro = ("Scientific Jigsaw interactive session.  'help' for commands, "
             "'boards' to list encodings, 'quit' to leave.\n")
    prompt = "jigsaw> "

    def __init__(self, board: str = "eif3"):
        super().__init__()
        self.board_name = ""
        self.mode = "merger"
        self.history: List[Tuple[str, object]] = []
        self.factual: Optional[TypedModel] = None
        self.model: Optional[TypedModel] = None
        self.nav: Optional[Navigator] = None
        self.seed = ""
        self.last_answer: Optional[Dict] = None
        self.do_load(board)

    # ------------------------------------------------------------ helpers
    def _sync(self) -> None:
        if self.mode == "seeded":
            self.nav = Navigator(self.model, self.seed)
        else:
            self.nav = None

    def _board(self) -> Dict:
        return BOARDS[self.board_name]

    def _count(self) -> int:
        if self.mode == "seeded":
            return self.nav.n_histories() if self.nav else 0
        return merger_F(self.model)(frozenset(self.model.V))

    def _count_of_snapshot(self, snap) -> int:
        if snap["mode"] == "seeded":
            return snap["nav"].n_histories() if snap["nav"] else 0
        return merger_F(snap["model"])(frozenset(snap["model"].V))

    def _parse_hypothesis(self, args: List[str]):
        """exact A B | contains A B | before A B | last A | step k A"""
        if not args:
            return None, "no hypothesis given"
        kind, rest = args[0].lower(), args[1:]
        try:
            if kind in ("exact", "forms"):
                if not rest:
                    return None, "usage: ask exact <components...>"
                return ("subset", frozenset(rest)), None
            if kind == "before":
                return before(rest[0], rest[1]), None
            if kind == "last":
                return last(rest[0]), None
            if kind == "step":
                return at_step(int(rest[0]), rest[1]), None
        except (IndexError, ValueError):
            return None, f"bad arguments for '{kind}'"
        return None, f"unknown hypothesis type '{kind}'"

    def _check_components(self, names) -> Optional[str]:
        unknown = [n for n in names if n not in self.model.V]
        if unknown:
            return (f"unknown component(s): {', '.join(unknown)}. "
                    f"this board has: {', '.join(self.model.V)}")
        return None

    def _push(self, label: str) -> None:
        """Store the complete session state, not only the model: several
        successive 'see' operations must undo one at a time."""
        snapshot = {"model": self.model,
                    "nav": self.nav,
                    "mode": self.mode,
                    "factual": self.factual}
        self.history.append((label, snapshot))

    # -------------------------------------------------------------- setup
    def do_boards(self, _):
        """List the built-in encodings."""
        for k, v in BOARDS.items():
            print(f"  {k:<12} {v['desc']}")

    def do_load(self, arg):
        """load <board>  -- load a built-in encoding (see 'boards')."""
        name = (arg or "eif3").strip()
        if name not in BOARDS:
            print(f"unknown board '{name}'. try: {', '.join(BOARDS)}")
            return
        b = BOARDS[name]
        self.board_name = name
        self.seed = b["seed"]
        self.model = TypedModel(V=b["V"], C=b["C"], P=b.get("P", ()),
                                X=b.get("X", frozenset()))
        self.factual = self.model
        self.history = []
        self._sync()
        print(f"loaded '{name}': {b['desc']}")
        self.do_status("")

    def do_mode(self, arg):
        """mode seeded|merger  -- switch representation."""
        m = arg.strip().lower()
        if m not in ("seeded", "merger"):
            print("usage: mode seeded|merger"); return
        b = self._board()
        if m == "seeded":
            self.model = TypedModel(V=b["V"], C=b["C"],
                                    P=b.get("P_seeded", b.get("P", ())),
                                    X=b.get("X", frozenset()))
        else:
            self.model = TypedModel(V=b["V"], C=b["C"], P=b.get("P", ()),
                                    X=b.get("X", frozenset()))
        self.mode = m
        self.factual = self.model        # the factual model is mode-specific
        self.history = []
        self._sync()
        print(f"mode: {m}")
        self.do_status("")

    def do_status(self, _):
        """Current board, mode, relations and mechanism count."""
        m = self.model
        print(f"  board {self.board_name} [{self.mode}]  seed={self.seed}")
        print(f"  components: {', '.join(m.V)}")
        print(f"  contacts:   {', '.join('--'.join(sorted(c)) for c in sorted(m.C, key=lambda s: sorted(s)))}")
        print(f"  prereqs:    {', '.join(f'{a}->{b}' for a, b in m.P) or '(none)'}")
        if m.X:
            print(f"  exclusions: {', '.join('--'.join(sorted(x)) for x in m.X)}")
        print(f"  admissible mechanisms: {self._count()}")
        if self.nav is not None and self.nav.trail:
            print(f"  conditioning: {'; '.join(self.nav.trail)}")
        if self.history:
            print(f"  steps taken: {len(self.history)} ('back' to undo)")

    # ---------------------------------------------------------- TEST / WHY
    def do_ask(self, arg):
        """ask exact A B... | before A B | last A | step k A   -- test a hypothesis."""
        args = shlex.split(arg)
        h, err = self._parse_hypothesis(args)
        if err:
            print(err); return
        if isinstance(h, tuple) and h[0] == "subset":
            E = h[1]
            bad = self._check_components(E)
            if bad:
                print(bad); return
            if self.mode == "merger":
                g, f = merger_support(self.model, E)
                if f == 0:
                    print("  no admissible mechanisms at all under this model"); return
                verdict = ("impossible" if g == 0 else
                           "necessary" if g == f else f"possible")
                print(f"  {sorted(E)}: {verdict}   {g}/{f} = {g/f:.3f}")
                if g == 0:
                    print(f"      cause: {self.model.blocking_reason(E) or 'no admissible tree contains it'}")
                self.last_answer = {"E": E, "n": g, "tot": f}
            else:
                r = self.nav.ask(self.nav.forms_exact(E))
                self._print_ask(r)
            return
        if self.mode == "merger":
            print("  ordering hypotheses need the seeded mode: 'mode seeded'"); return
        r = self.nav.ask(h)
        self._print_ask(r)

    def _print_ask(self, r: Dict) -> None:
        print(f"  {r['hypothesis']}: {r['verdict']}   p = {r['prob']:.3f}")
        if r.get("cause"):
            print(f"      cause: {r['cause']}")
        if r.get("given"):
            print(f"      given: {'; '.join(r['given'])}")
        self.last_answer = r

    def do_why(self, _):
        """Explain the last answer in full."""
        if not self.last_answer:
            print("nothing asked yet"); return
        r = self.last_answer
        if "E" in r:
            E = r["E"]
            print(f"  query {sorted(E)} occurs in {r['n']} of {r['tot']} mechanisms")
            if r["n"] == 0:
                print(f"  cause: {self.model.blocking_reason(E) or 'unreachable under this representation'}")
                print("  try: 'mode merger' if seeded, or 'do drop-p ...' to relax a prerequisite")
        else:
            print(f"  {r['hypothesis']}: {r['verdict']}")
            if r.get("cause"):
                print(f"  cause: {r['cause']}")

    def do_queries(self, _):
        """Test all the recorded queries of this board at once."""
        qs = self._board().get("queries", {})
        if not qs:
            print("no recorded queries for this board"); return
        for label, comps in qs.items():
            E = frozenset(comps)
            if self.mode == "merger":
                g, f = merger_support(self.model, E)
                tag = "IMPOSSIBLE" if g == 0 else f"{g}/{f} = {g/f:.3f}"
            else:
                r = self.nav.ask(self.nav.forms_exact(E))
                tag = "IMPOSSIBLE" if r["prob"] == 0 else f"{r['prob']:.3f}"
            print(f"  {label:<10} {sorted(E)}  {tag}")

    # -------------------------------------------------------------- narrow
    def do_see(self, arg):
        """see before A B | last A | step k A   -- suppose an observation."""
        if self.mode != "seeded":
            print("  conditioning needs the seeded mode: 'mode seeded'"); return
        h, err = self._parse_hypothesis(shlex.split(arg))
        if err or isinstance(h, tuple):
            print(err or "use an ordering hypothesis here"); return
        try:
            new = self.nav.see(h)
        except ValueError as e:
            print(f"  {e}"); return
        self._push(f"see {arg}")
        self.nav = new
        print(f"  conditioned on {getattr(h, 'label', arg)}: "
              f"{self.nav.n_histories()} mechanisms remain "
              f"(N_eff {self.nav.n_eff():.2f})")

    # -------------------------------------------------------------- UPDATE
    def do_do(self, arg):
        """do drop-p A B | drop-c A B | drop-x A B   -- edit the declared model."""
        args = shlex.split(arg)
        if len(args) != 3:
            print("usage: do drop-p|drop-c|drop-x A B"); return
        kind, a, b = args[0].lower(), args[1], args[2]
        bad = self._check_components([a, b])
        if bad:
            print(bad); return
        self._push(f"do {arg}")
        if kind == "drop-p":
            self.model = drop_relations(self.model, drop_P=[(a, b)])
        elif kind == "drop-c":
            self.model = drop_relations(self.model, drop_C=[(a, b)])
        elif kind == "drop-x":
            self.model = drop_relations(self.model, drop_X=[frozenset({a, b})])
        else:
            print("unknown edit"); self.history.pop(); return
        self._sync()
        print(f"  edited model: {self._count()} admissible mechanisms "
              f"(was {self._count_of_snapshot(self.history[-1][1])})")

    def _count_of(self, model) -> int:
        if self.mode == "seeded":
            return Navigator(model, self.seed).n_histories()
        return merger_F(model)(frozenset(model.V))

    def do_back(self, _):
        """Undo the last see/do step, restoring model AND conditioning."""
        if not self.history:
            print("nothing to undo"); return
        label, snap = self.history.pop()
        self.model = snap["model"]
        self.mode = snap["mode"]
        self.factual = snap["factual"]
        self.nav = snap["nav"]            # restores the conditioning trail
        n = self._count()
        trail = "; ".join(self.nav.trail) if (self.nav and self.nav.trail) else "none"
        print(f"  undid '{label}': {n} mechanisms   conditioning: {trail}")

    def do_reset(self, _):
        """Return to the factual model as loaded."""
        self.do_load(self.board_name)

    # ------------------------------------------------------------- COMPARE
    def do_compare(self, _):
        """Counterfactual report: current model versus the model as loaded."""
        if self.model is self.factual:
            print("  current model is the factual one; make an edit first"); return
        if self.mode != "seeded":
            n0 = merger_F(self.factual)(frozenset(self.factual.V))
            n1 = self._count()
            print(f"  mechanisms {n0} -> {n1}")
            for label, comps in self._board().get("queries", {}).items():
                E = frozenset(comps)
                g0, f0 = merger_support(self.factual, E)
                g1, f1 = merger_support(self.model, E)
                s0 = g0 / f0 if f0 else 0.0
                s1 = g1 / f1 if f1 else 0.0
                print(f"    {label:<10} support {s0:.3f} -> {s1:.3f}  ({s1-s0:+.3f})")
            return
        base = Navigator(self.factual, self.seed)
        r = base.compare_with(Navigator(self.model, self.seed))
        print(f"  mechanisms {r['n_factual']} -> {r['n_alternative']}")
        print(f"  killed {len(r['killed'])} (mass {r['killed_mass_factual']:.3f}), "
              f"resurrected {len(r['resurrected'])} (mass {r['resurrected_mass_alt']:.3f})")
        print(f"  dH2 {r['delta_H2_bits']:+.3f} bits   dN_eff {r['delta_N_eff']:+.2f}")
        ed = {k: v for k, v in r["edited_relations"].items() if v}
        print(f"  cause: {ed or 'no relation changed'}")

    # -------------------------------------------------------------- NAVIGATE
    def do_walk(self, arg):
        """walk [A B ...]  -- stand after these additions and look around."""
        if self.mode != "seeded":
            print("  walking the tree needs the seeded mode: 'mode seeded'"); return
        added = shlex.split(arg)
        bad = self._check_components(added)
        if bad:
            print(bad); return
        w = self.nav.walk(added)
        if not w["reachable"]:
            print("  that state is not reachable under the current model"); return
        print(f"  at {'->'.join(w['state'])}   mass here {w['mass_here']:.3f}")
        print("  next:  " + ", ".join(f"{k} {v:.3f}" for k, v in w["branches"].items()))
        print(f"  branch entropy {w['branch_entropy_bits']:.3f} bits; "
              f"remaining N_eff {w['remaining_n_eff']:.2f}")

    def do_branches(self, arg):
        """branches [k]  -- the k most uncertain decision points."""
        if self.mode != "seeded":
            print("  needs the seeded mode: 'mode seeded'"); return
        k = int(arg) if arg.strip().isdigit() else 3
        for r in self.nav.uncertain_branches(k):
            print(f"  {'->'.join(r['state']):<28} mass {r['mass_here']:.3f}  "
                  f"H2 {r['branch_entropy_bits']:.3f} bits  -> "
                  + ", ".join(f"{c}:{p:.2f}" for c, p in r["branches"].items()))

    def do_top(self, arg):
        """top [k]  -- most probable mechanisms."""
        if self.mode != "seeded":
            print("  needs the seeded mode: 'mode seeded'"); return
        k = int(arg) if arg.strip().isdigit() else 5
        for h, p in self.nav.top(k):
            print(f"  {p:.4f}  {self.seed}->" + "->".join(h))

    # ----------------------------------------------------------------- misc
    def do_quit(self, _):
        """Leave the session."""
        print("bye")
        return True

    do_EOF = do_quit
    do_exit = do_quit

    def emptyline(self):
        pass

    def default(self, line):
        print(f"unknown command: {line.split()[0]}. type 'help'.")


if __name__ == "__main__":
    board = sys.argv[1] if len(sys.argv) > 1 else "eif3"
    try:
        Session(board).cmdloop()
    except KeyboardInterrupt:
        print("\nbye")
