# Request for a hostile review: break-even sizes of the Appendix D route

Please try to refute the computation below. We are not asking for encouragement; we want the weakest points found.

## What this is about

The paper "The Exponential Time Hypothesis is False" (OpenAI, 2026) gives a 2^o(n) algorithm for counting
satisfying assignments of 3-CNF formulas. Being subexponential, it beats every c^n algorithm for large enough n.
We asked from which n it actually does, for the tensor-contraction route of its Appendix D, on 3-CNF in which every
variable occurs at most three times (the algorithm's best case). Over eleven rounds a loop (Claude as head writing
ideas, proofs and a frozen referee; a local model writing code; several ideas from GPT-6 Astra and corrections from
Fable) lowered the break-even sizes. After round 11:

| competitor | log2 n from which the route stays ahead | n about |
|---|---|---|
| 2^n (brute force) | 103 | 10^31 |
| 1.637^n (#3-SAT record cited by the paper, all formulas) | 159 | 10^48 |
| 2^(n/3) (counting for this class: path decomposition, Fomin and Hoie 2006) | 533 | 10^160 |
| 1.307^n (PPSZ, decision only, all formulas) | 413 | 10^124 |
| 1.1199^n (decision for this class, arXiv 2411.07389) | 2,034 | 10^612 |

As written, the route's cost formula beats 2^(n/3) only from log2 n = 1.11 x 10^11 (and the paper's own invariant
u >= lambda^(4/5) holds only from about 3.9 x 10^25); a tighter accounting of the paper's own construction starts at
496,177. We do not claim practical relevance: at a break-even both algorithms take about 2^(n/3) steps.

## How the numbers are computed

- appd_assembly.py (attached, frozen) is the cost model. For parameters (levels, groups r, tuple size u, packet
  size a, balance tolerance delta, chunk length Lc, and flags that switch on our modifications of the paper's
  construction) it returns an exponent e with cost at most 2^(e n) poly(n), or rejects the parameters when a size
  condition fails (precision headroom, label balance, padding; lines 270-310). The exponent is
  e = max(worst, cum + leaf), with one term per level of the form tup(u, eps, ...) * M_d / sqrt(u_d).
- certificates.json (attached) gives, for every cell, the parameters and the exponent the referee found.
- results_page.pdf (attached) is the results page: every round, its tables and its history.
- claims.md (attached) is the ledger: every flag is a named step (appd-tight, appd-bernstein, appd-memo,
  appd-modedp, appd-chunks, appd-perclass, appd-leafvars, appd-labelfree, appd-newton, appd-ldu, appd-shellmm,
  appd-guarded, appd-bincross, appd-profile) with a proof sketch, its toy checks and the outside reviews so far.
- A frozen referee verified the code of every step against exact counts on small formulas and hidden test nodes.
  It checks identities and mechanics, not the asymptotic bounds; those rest on the proof sketches.
- A break-even is the smallest log2 n from which e stays below the competitor's exponent. In the model e does not
  depend on log2 n and every size condition only relaxes as log2 n grows, so the cells are monotone; after round 10
  we also scanned the search from log2 n = 64 up to 10^9.
- Status of each ingredient, kept apart: (i) code checked by the referee on small instances; (ii) the cost
  formulas' constants and certificates checked in exact integer arithmetic; (iii) the bounds themselves resting on
  head-level proof sketches, not peer reviewed; (iv) the paper's own theorem assumed throughout.

What fixes the sizes: in every cell of round 11 the binding condition is a label balance condition, at the top level
in four cells and at depth 2 against 1.1199^n (the exponent e itself does not depend on log2 n in the model). Against
2^(n/3) the top level has r = 8 groups and u = 79, so 8^158 = 2^474 label types, and the
balancing argument asks for 2^58.2 wires of each type: n must exceed 2^532.2, hence 533
(0.83 bits to spare). Against 2^n: 7^22 = 2^61.8 types and 40.8 bits give 103.

## What we ask

1. The steps. For each appd-* step in claims.md: is the proof sketch correct? Give a gap or a small counterexample
   where you find one, and rank the three weakest steps. The newest steps had the least scrutiny: appd-profile
   (round 11: the cheapest of a menu of contraction orders per shell profile, with an averaging certificate),
   appd-bincross (round 10), appd-guarded (round 8), appd-shellmm (round 7).
2. Upper bound. Is the model an upper bound on the algorithm's work at these sizes? Look for terms it undercounts:
   the number of profiles, block padding of the matrix products, precision growth across levels, memoisation of
   the child tables, memory traffic.
3. Size conditions. Is the root label balance (Bernstein bound plus union bound, slack formula in lines 288-301)
   correct, and are the precision headroom and padding conditions sufficient for the paper's recovery step?
4. Competitors. Are 2^(n/3) and 1.1199^n the right fair competitors for 3-CNF with at most three occurrences per
   variable? Is a faster counting algorithm known for this class? Are 1.637^n and 1.307^n the right records for
   general formulas?
5. The paper. Do our modifications keep the paper's correctness conditions (Definition 3.1, Lemmas 4.5 and 4.7,
   Lemma D.5)? Where does a modification rely on a property the paper does not state?
6. General formulas. With m clauses the root has about 2m - n wires instead of n; at m = 4.26n the route needs
   6,063 / 20,415 / 3,089 against 1.637^n / 1.307^n / 2^n, at m = 5n 8,660 / 29,219 / 4,401 (n about 10^1325
   against 2^n). Anthropic's own paper is reported to need n > 10^10000 against 2^n at m = 5n for its algorithm.
   Is our density model (chains of copies above 1,024 occurrences, wire count, balance cost) right, and is that
   comparison fair?
7. Anything else that would move a cell by more than 1 percent, in either direction.

For each point: the step or line, the severity (breaks the result / changes a cell / cosmetic), the argument, and a
small computation or counterexample if you can. "Checked, no problem found" for a step is useful too.

A note on independence: three of the steps came from GPT-6 Astra in earlier conversations, appd-guarded (round 8),
appd-bincross (round 10) and appd-profile (round 11, which Claude extended to menus of three and four orders).
Please review them as if someone else had proposed them, and re-derive their key inequalities instead of relying
on any earlier discussion.
