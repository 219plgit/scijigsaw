"""
sj_visual.py -- generate a standalone, interactive HTML view of an assembly model.

Enumeration happens once in Python; the resulting mechanism space (and one
variant per removable relation) is embedded in a single HTML file. The browser
then does conditioning, hypothesis testing and counterfactual comparison
instantly, with no server and no dependencies.

    python3 sj_visual.py                # writes jigsaw_eif3.html
    python3 sj_visual.py proteasome     # any board defined in sj_session.py
    open jigsaw_eif3.html               # macOS; or double-click the file

In the page:
  * click a node to condition on reaching that state; click the root to reset
  * tick components to ask whether that exact subcomplex can form
  * choose a relation to remove and see which mechanisms die or revive
"""

from __future__ import annotations

import json
import sys
from typing import Dict, List, Tuple

from sj_prob import TypedModel, drop_relations
from sj_navigate import Navigator
from sj_session import BOARDS


def merger_trees(model: TypedModel) -> List:
    """All admissible merger trees as nested lists, deduplicated."""
    import itertools

    def rec(S):
        S = frozenset(S)
        if len(S) == 1:
            return [next(iter(S))]
        if not model.valid_intermediate(S):
            return []
        out, anchor = [], min(S)
        for r in range(len(S)):
            for combo in itertools.combinations(sorted(S - {anchor}), r):
                A = frozenset({anchor}) | frozenset(combo)
                B = S - A
                if not B or not (model.valid_intermediate(A)
                                 and model.valid_intermediate(B)):
                    continue
                if not any(model.has_contact(u, v) for u in A for v in B):
                    continue
                for ta in rec(A):
                    for tb in rec(B):
                        out.append([ta, tb])
        return out

    seen, uniq = set(), []
    for t in rec(frozenset(model.V)):
        k = json.dumps(canon(t))
        if k not in seen:
            seen.add(k); uniq.append(canon(t))
    return uniq


def canon(t):
    if isinstance(t, str):
        return t
    a, b = canon(t[0]), canon(t[1])
    return [a, b] if json.dumps(a) <= json.dumps(b) else [b, a]


def space(model: TypedModel, seed: str) -> Dict:
    """Enumerate one model into a JSON-ready description."""
    nav = Navigator(model, seed)
    hists = [list(h) for h, _ in nav.masses]
    return {
        "histories": hists,
        "n": len(hists),
        "relations": {
            "P": [list(p) for p in model.P],
            "C": [sorted(c) for c in model.C],
            "X": [sorted(x) for x in model.X],
        },
    }


def build(board_name: str) -> Tuple[str, Dict]:
    b = BOARDS[board_name]
    seed = b["seed"]
    model = TypedModel(V=b["V"], C=b["C"],
                       P=b.get("P_seeded", b.get("P", ())),
                       X=b.get("X", frozenset()))
    merger_model = TypedModel(V=b["V"], C=b["C"], P=b.get("P", ()),
                              X=b.get("X", frozenset()))
    data: Dict = {
        "merger": merger_trees(merger_model),
        "board": board_name,
        "desc": b["desc"],
        "seed": seed,
        "components": [c for c in b["V"] if c != seed],
        "prereqs": [list(p) for p in b.get("P_seeded", b.get("P", ()))],
        "all_components": list(b["V"]),
        "factual": space(model, seed),
        "variants": {},
        "queries": {k: sorted(v) for k, v in b.get("queries", {}).items()},
    }
    for (a, c) in model.P:
        data["variants"][f"prereq {a}->{c}"] = space(
            drop_relations(model, drop_P=[(a, c)]), seed)
    for e in sorted(model.C, key=lambda s: sorted(s)):
        u, v = sorted(e)
        data["variants"][f"contact {u}--{v}"] = space(
            drop_relations(model, drop_C=[(u, v)]), seed)
    for e in sorted(model.X, key=lambda s: sorted(s)):
        u, v = sorted(e)
        data["variants"][f"exclusion {u}--{v}"] = space(
            drop_relations(model, drop_X=[e]), seed)
    return seed, data


HTML = r"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>Scientific Jigsaw &mdash; __BOARD__</title>
<style>
 body{font:14px/1.45 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;
      margin:0;color:#111;background:#fff}
 header{padding:14px 20px;border-bottom:1px solid #ddd}
 h1{font-size:17px;margin:0 0 3px} .sub{color:#666;font-size:12.5px}
 main{display:flex;gap:0;align-items:stretch;min-height:calc(100vh - 62px)}
 #treewrap{flex:1 1 auto;overflow:auto;padding:14px 8px}
 aside{width:330px;flex:0 0 330px;border-left:1px solid #ddd;padding:14px 16px;
       overflow:auto;background:#fafafa}
 .panel{margin-bottom:20px}
 .panel h2{font-size:12px;letter-spacing:.06em;text-transform:uppercase;
           color:#555;margin:0 0 8px;font-weight:600}
 .stat{display:flex;justify-content:space-between;padding:2px 0;font-variant-numeric:tabular-nums}
 .stat span:last-child{font-weight:600}
 .chips{display:flex;flex-wrap:wrap;gap:5px}
 .chip{border:1px solid #bbb;border-radius:13px;padding:3px 9px;cursor:pointer;
       background:#fff;font-size:12.5px;user-select:none}
 .chip.on{background:#111;color:#fff;border-color:#111}
 .verdict{margin-top:9px;padding:9px 11px;border-radius:6px;background:#fff;
          border:1px solid #ddd}
 .verdict b{display:block;margin-bottom:2px}
 .imp{border-color:#c00;background:#fff5f5} .nec{border-color:#070;background:#f4fbf4}
 select,button{font:inherit;padding:5px 8px;border:1px solid #bbb;border-radius:5px;
               background:#fff;width:100%;margin-top:4px}
 button{cursor:pointer;width:auto;padding:5px 12px}
 .cause{color:#a00;font-size:12.5px;margin-top:5px}
 .trail{font-size:12.5px;color:#555;margin-top:6px;min-height:1.2em}
 .mech{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12px;
       padding:1px 0;color:#333}
 node{cursor:pointer}
 .nd{fill:#fff;stroke:#333;stroke-width:1.4}
 .nd.dim{opacity:.25} .nd.sel{fill:#111} .nd.hit{fill:#ffe08a;stroke:#a67c00}
 .lk{stroke:#999;fill:none} .lk.dim{opacity:.15} .lk.on{stroke:#111;stroke-width:2}
 .lbl{font-size:11.5px;fill:#111} .lbl.dim{opacity:.3}
 .plb{font-size:10px;fill:#666} .plb.dim{opacity:.25}
 .legend{font-size:12px;color:#666;padding:0 20px 12px}
 .modes{padding:9px 20px 0;display:flex;align-items:center;gap:0}
 .mode{border:1px solid #bbb;background:#fff;padding:5px 14px;cursor:pointer;
       font:inherit;width:auto;margin:0}
 .mode:first-child{border-radius:6px 0 0 6px} .mode:last-child{border-radius:0 6px 6px 0;border-left:none}
 .mode.on{background:#111;color:#fff;border-color:#111}
 #modenote{margin-left:14px;color:#666;font-size:12.5px}
 #grid{display:flex;flex-wrap:wrap;gap:10px;padding:6px 4px}
 .tcard{border:1px solid #ddd;border-radius:7px;padding:6px 8px 2px;background:#fff}
 .tcard.hit{border-color:#a67c00;background:#fffdf3;box-shadow:0 0 0 2px #ffe08a inset}
 .tcard .cap{font-size:10.5px;color:#888;text-align:center;margin-top:1px}
</style></head><body>
<header>
 <h1>Scientific Jigsaw &mdash; __DESC__</h1>
 <div class="sub">Click a node to suppose the assembly reached it. Tick components to test whether that exact subcomplex can form.</div>
</header>
<div class="modes">
 <button id="m-seeded" class="mode on">seed-anchored</button
 ><button id="m-merger" class="mode">merger trees</button>
 <span id="modenote"></span>
</div>
<main>
 <div id="treewrap"><svg id="tree"></svg><div id="grid"></div></div>
 <aside>
  <div class="panel"><h2>Current space</h2>
    <div class="stat"><span>mechanisms</span><span id="s-n"></span></div>
    <div class="stat"><span>entropy (bits)</span><span id="s-h"></span></div>
    <div class="stat"><span>N<sub>eff</sub></span><span id="s-e"></span></div>
    <div class="trail" id="s-trail"></div>
    <button id="reset">reset conditioning</button>
  </div>
  <div class="panel"><h2>Can this subcomplex form?</h2>
    <div class="chips" id="chips"></div>
    <div class="verdict" id="verdict">tick one or more components</div>
  </div>
  <div class="panel"><h2>Counterfactual: remove a relation</h2>
    <select id="variant"></select>
    <div class="verdict" id="cf">factual model</div>
  </div>
  <div class="panel"><h2>Mechanisms</h2><div id="mechs"></div></div>
 </aside>
</main>
<div class="legend">Node = assembled state. Edge label = probability of that next
 addition under the current space. Yellow = state matching the ticked subcomplex.</div>
<script>
const DATA = __DATA__;
let cur = DATA.factual, cond = [], picked = new Set(), variantName = "", mode = "seeded";

const H2 = n => n>0 ? Math.log2(n) : 0;

function filtered(){
  return cur.histories.filter(h => cond.every((c,i)=>h[i]===c));
}
function buildTree(hs){
  const root = {name:DATA.seed, path:[], kids:{}, n:0};
  for(const h of hs){
    let node = root; node.n++;
    h.forEach((p,i)=>{
      node.kids[p] = node.kids[p] || {name:p, path:h.slice(0,i+1), kids:{}, n:0};
      node = node.kids[p]; node.n++;
    });
  }
  return root;
}
function layout(root){
  const levels=[]; let maxw=0;
  (function walk(n,d){ (levels[d]=levels[d]||[]).push(n);
    Object.values(n.kids).forEach(k=>walk(k,d+1)); })(root,0);
  levels.forEach(l=>maxw=Math.max(maxw,l.length));
  const DX=Math.max(120, 780/Math.max(1,maxw)), DY=78;
  levels.forEach((l,d)=>l.forEach((n,i)=>{
    n.x = (i+0.5)*DX + 20; n.y = d*DY + 34;
  }));
  return {levels, w: maxw*DX+60, h: levels.length*DY+40};
}
function draw(){
  if(mode === "merger"){ drawMerger(); return; }
  document.getElementById('grid').innerHTML = '';
  const hs = filtered(), root = buildTree(hs), L = layout(root);
  const svg = document.getElementById('tree');
  svg.setAttribute('width', L.w); svg.setAttribute('height', L.h);
  let out = '';
  (function walk(n){
    for(const k of Object.values(n.kids)){
      const p = (k.n/n.n);
      out += `<path class="lk" d="M${n.x},${n.y+11} C${n.x},${n.y+45} ${k.x},${k.y-45} ${k.x},${k.y-13}"/>`;
      out += `<text class="plb" x="${(n.x+k.x)/2}" y="${(n.y+k.y)/2+3}" text-anchor="middle">${p.toFixed(2)}</text>`;
      walk(k);
    }
  })(root);
  const target = picked.size ? new Set([...picked, DATA.seed]) : null;
  (function walk2(n){
    const state = new Set([DATA.seed, ...n.path]);
    let cls = 'nd';
    if(target && state.size===target.size && [...target].every(x=>state.has(x))) cls += ' hit';
    if(cond.length && n.path.length<=cond.length &&
       n.path.every((p,i)=>p===cond[i])) cls += ' sel';
    out += `<circle class="${cls}" cx="${n.x}" cy="${n.y}" r="11" `
        +  `onclick="pick(${JSON.stringify(n.path).replace(/"/g,'&quot;')})"/>`;
    out += `<text class="lbl" x="${n.x}" y="${n.y-17}" text-anchor="middle">${n.name}</text>`;
    Object.values(n.kids).forEach(walk2);
  })(root);
  svg.innerHTML = out;

  document.getElementById('s-n').textContent = hs.length;
  document.getElementById('s-h').textContent = H2(hs.length).toFixed(3);
  document.getElementById('s-e').textContent = hs.length.toFixed(2);
  document.getElementById('s-trail').textContent =
     cond.length ? 'given: ' + DATA.seed + '\u2192' + cond.join('\u2192') : 'no conditioning';
  document.getElementById('mechs').innerHTML =
     hs.slice(0,40).map(h=>`<div class="mech">${DATA.seed}\u2192${h.join('\u2192')}</div>`).join('')
     + (hs.length>40 ? `<div class="mech">\u2026 ${hs.length-40} more</div>` : '');
  verdict();
}
function pick(path){ cond = path; draw(); }

/* ---------------- merger-tree view ---------------- */
function tsize(t){ return (typeof t === "string") ? 1 : tsize(t[0]) + tsize(t[1]); }
function leaves(t){ return (typeof t === "string") ? [t] : leaves(t[0]).concat(leaves(t[1])); }
function subcomplexes(t, acc){
  /* every node of the tree is a species that exists during assembly:
     leaves are free monomers, internal nodes are assembled subcomplexes. */
  acc = acc || [];
  if(typeof t === "string"){ acc.push(t); return acc; }
  acc.push(leaves(t).slice().sort().join("|"));
  subcomplexes(t[0], acc); subcomplexes(t[1], acc);
  return acc;
}
function drawTree(t, W, H){
  const LH = 26, order = leaves(t), gap = Math.max(34, W/Math.max(1,order.length));
  const xs = {}; order.forEach((n,i)=> xs[n] = i*gap + 18);
  let out = '', maxd = 0;
  function place(node, depth){
    maxd = Math.max(maxd, depth);
    if(typeof node === "string") return {x: xs[node], y: H - 16, leaf:true, name:node};
    const a = place(node[0], depth+1), b = place(node[1], depth+1);
    const x = (a.x + b.x)/2, y = Math.min(a.y, b.y) - LH;
    out += `<path class="lk" d="M${a.x},${a.y-4} L${a.x},${y} L${b.x},${y} L${b.x},${b.y-4}"/>`;
    return {x, y, leaf:false};
  }
  place(t, 0);
  order.forEach(n=>{ out += `<text class="lbl" x="${xs[n]}" y="${H-4}" `
      + `text-anchor="middle" transform="rotate(-32 ${xs[n]},${H-4})">${n}</text>`; });
  return {svg: out, w: order.length*gap + 24};
}
function fmt(t){ return (typeof t==="string") ? t : "(" + fmt(t[0]) + "+" + fmt(t[1]) + ")"; }
function drawMerger(){
  document.getElementById('tree').innerHTML = '';
  document.getElementById('tree').setAttribute('height', 0);
  const want = picked.size ? [...picked].sort().join("|") : null;
  const g = document.getElementById('grid');
  let hits = 0, html = '';
  DATA.merger.forEach((t,i)=>{
    const subs = new Set(subcomplexes(t));
    const hit = want && subs.has(want);
    if(hit) hits++;
    const H = 118, r = drawTree(t, 150, H);
    html += `<div class="tcard${hit?' hit':''}">`
         +  `<svg width="${r.w}" height="${H}">${r.svg}</svg>`
         +  `<div class="cap">tree ${i+1}</div></div>`;
  });
  g.innerHTML = html;
  document.getElementById('s-n').textContent = DATA.merger.length;
  document.getElementById('s-h').textContent = H2(DATA.merger.length).toFixed(3);
  document.getElementById('s-e').textContent = DATA.merger.length.toFixed(2);
  document.getElementById('s-trail').textContent =
     'merger-tree representation (independent subcomplexes may form and join)';
  document.getElementById('mechs').innerHTML =
     DATA.merger.slice(0,14).map((t,i)=>`<div class="mech">${i+1}: ${fmt(t)}</div>`).join('');
  const el = document.getElementById('verdict');
  if(!want){ el.className='verdict'; el.textContent='tick components to test a subcomplex'; }
  else {
    const name = [...picked].sort().join(' + ');
    if(hits===0){ el.className='verdict imp';
      el.innerHTML = `<b>impossible</b>${name}<div class="cause">no admissible merger tree `
        + `contains this subcomplex</div>`; }
    else if(hits===DATA.merger.length){ el.className='verdict nec';
      el.innerHTML = `<b>necessary</b>${name}<div>all ${hits} trees contain it</div>`; }
    else { el.className='verdict';
      el.innerHTML = `<b>possible &mdash; ${(hits/DATA.merger.length).toFixed(3)}</b>${name}`
        + `<div>${hits} of ${DATA.merger.length} trees</div>`; }
  }
}
function verdict(){
  const el = document.getElementById('verdict');
  if(!picked.size){ el.className='verdict'; el.textContent='tick one or more components'; return; }
  const want = new Set([...picked, DATA.seed]);
  const hs = filtered();
  let n = 0;
  for(const h of hs){
    let cur = new Set([DATA.seed]); let ok = (cur.size===want.size);
    for(const p of h){ cur.add(p);
      if(cur.size===want.size && [...want].every(x=>cur.has(x))){ ok=true; break; } }
    if(ok) n++;
  }
  const p = hs.length ? n/hs.length : 0;
  const name = [...want].sort().join(' + ')
     + (mode==="seeded" && !picked.has(DATA.seed)
        ? ` <span style="color:#888">(seed ${DATA.seed} included automatically)</span>` : '');
  if(n===0){
    el.className='verdict imp';
    el.innerHTML = `<b>impossible</b>${name}<div class="cause">${causeFor(want)}</div>`;
  } else if(n===hs.length){
    el.className='verdict nec';
    el.innerHTML = `<b>necessary</b>${name}<div>every one of the ${hs.length} mechanisms forms it</div>`;
  } else {
    el.className='verdict';
    el.innerHTML = `<b>possible &mdash; ${p.toFixed(3)}</b>${name}`
      + `<div>${n} of ${hs.length} mechanisms</div>`;
  }
}
function causeFor(want){
  /* name the decisive relation rather than returning a bare zero */
  const blocking = (DATA.prereqs||[]).filter(([p,q]) => want.has(q) && !want.has(p));
  if(blocking.length){
    const list = blocking.map(([p,q])=>`${p}\u2192${q}`).join(", ");
    return `blocked by declared prerequisite ${list}: the target cannot be present `
         + `before its prerequisite, so this exact state never occurs.`;
  }
  if(cond.length){
    return `excluded by the current conditioning (${DATA.seed}\u2192${cond.join("\u2192")}); `
         + `reset conditioning to test it against the full space.`;
  }
  return `not expressible in the seed-anchored representation: every intermediate `
       + `contains ${DATA.seed}, so an off-scaffold species cannot appear. `
       + `Switch to merger trees to express it.`;
}
function setVariant(name){
  variantName = name;
  cur = name ? DATA.variants[name] : DATA.factual;
  cond = [];
  const el = document.getElementById('cf');
  if(!name){ el.className='verdict'; el.textContent='factual model'; }
  else {
    const f = new Set(DATA.factual.histories.map(h=>h.join('>')));
    const a = new Set(cur.histories.map(h=>h.join('>')));
    const killed = [...f].filter(x=>!a.has(x)).length;
    const res = [...a].filter(x=>!f.has(x)).length;
    el.className='verdict';
    el.innerHTML = `<b>removed: ${name}</b>`
      + `<div>mechanisms ${DATA.factual.n} \u2192 ${cur.n}</div>`
      + `<div>killed ${killed}, resurrected ${res}</div>`
      + `<div>\u0394 entropy ${(H2(cur.n)-H2(DATA.factual.n)>=0?'+':'')}`
      + `${(H2(cur.n)-H2(DATA.factual.n)).toFixed(3)} bits</div>`;
  }
  draw();
}
// build controls
const chips = document.getElementById('chips');
DATA.all_components.forEach(c=>{
  const d=document.createElement('div'); d.className='chip'; d.textContent=c;
  d.onclick=()=>{ d.classList.toggle('on');
    picked.has(c)?picked.delete(c):picked.add(c); draw(); };
  chips.appendChild(d);
});
const sel = document.getElementById('variant');
sel.innerHTML = '<option value="">(factual model)</option>' +
  Object.keys(DATA.variants).map(k=>`<option>${k}</option>`).join('');
sel.onchange = e => setVariant(e.target.value);
document.getElementById('reset').onclick = () => { cond=[]; draw(); };
function setMode(m){
  mode = m;
  document.getElementById('m-seeded').classList.toggle('on', m==="seeded");
  document.getElementById('m-merger').classList.toggle('on', m==="merger");
  document.getElementById('modenote').textContent = m==="seeded"
    ? "every intermediate contains " + DATA.seed
    : "independent subcomplexes may form and then join";
  document.getElementById('variant').disabled = (m==="merger");
  cond = []; draw();
}
document.getElementById('m-seeded').onclick = ()=>setMode('seeded');
document.getElementById('m-merger').onclick = ()=>setMode('merger');
setMode('seeded');
</script></body></html>
"""


def main() -> None:
    board = sys.argv[1] if len(sys.argv) > 1 else "eif3"
    if board not in BOARDS:
        print(f"unknown board '{board}'. available: {', '.join(BOARDS)}")
        return
    seed, data = build(board)
    html = (HTML.replace("__DATA__", json.dumps(data))
                .replace("__BOARD__", board)
                .replace("__DESC__", data["desc"]))
    out = f"jigsaw_{board}.html"
    with open(out, "w") as fh:
        fh.write(html)
    print(f"wrote {out}  ({data['factual']['n']} seeded, "
          f"{len(data['merger'])} merger trees, "
          f"{len(data['variants'])} counterfactual variants precomputed)")
    print(f"open it with:  open {out}")


if __name__ == "__main__":
    main()
