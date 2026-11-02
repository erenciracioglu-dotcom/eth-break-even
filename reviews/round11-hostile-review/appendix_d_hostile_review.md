**Independent hostile review of the Appendix D break-even model**  
2 November 2026

**Verdict.** I reproduced the five supplied certificates and all thirteen profile-menu constants. I did not find a counterexample to the algebra of `appd-profile`, `appd-bincross`, `appd-shellmm`, or `appd-guarded`. I did find a direct counterexample to the claim that the reported thresholds are the smallest thresholds of the frozen model: **the unchanged model accepts 100 / 157 / 532 / 411 / 2,033**, in the order brute force / general counting / class counting / general decision / class decision. These are feasible witnesses, not new claims of optimality.

The general-density argument has a separate, genuine gap: its global `2m - n + n/Lc` wire inequality omits small connected components. The bounded-occurrence argument does not have this particular defect. The numerical general-density cells cannot be independently reproduced from this package, which contains neither their model patch nor their parameter certificates.

The report also gives a simpler deterministic balancing construction. It proves a sufficient condition linear in the maximum chunk degree, in place of the ledger's quadratic dependence on chunk size. This supports the existing conservative root balance test and gives further, explicitly separated, model improvements.

The five original cells should therefore be described as **conditional feasible leading-exponent certificates found by a particular search**, not as proved least break-even sizes or fully audited running-time comparisons.

**What was actually checked.** I read the attached brief, ledger, model, JSON certificates and results PDF; inspected the current source of the OpenAI contraction proof and the relevant competitor papers; ran the frozen exponent function; independently reconstructed the profile forms and their exact integer certificates; checked finite LDU majorization, Hermite losses, and binary-crossing factorizations. I did not have `appd_core.py`, `appd_gates.py`, the frozen referee, the hidden cases, or the density patch. The reported historical referee passes are evidence supplied by the ledger, not tests independently rerun in this review. This is a mathematical/code-model audit, not a verification of either ETH paper's theorem.

**1. A concrete refutation: the stated thresholds are not minimal.**

Severity: **changes cells**. Locations: the brief's definition of a break-even; ledger lines 965–981 and 1049–1055; model lines 230–328. No mathematical flag was added or removed to obtain the following witnesses.

| Target | Supplied log2 n | New feasible log2 n | New e | Target minus e |
|---|---:|---:|---:|---:|
| 2^n | 103 | **100** | 0.999534479794686 | 0.000465520205314 |
| 1.637^n | 159 | **157** | 0.710254458410584 | 0.000799863401483 |
| 2^(n/3) | 533 | **532** | 0.333319277422392 | 0.000014055910941 |
| 1.307031578^n | 413 | **411** | 0.386265410589750 | 0.000028586606076 |
| 1.1199^n | 2,034 | **2,033** | 0.163366065628536 | 0.000003848846240 |

The first two reductions exceed 1% in log2 n: 2.91% and 1.26%. The new parameters are deliberately rounded away from the exponent boundary.

| Target | r | u | a | delta | Lc |
|---|---|---|---|---:|---:|
| Brute | [7] | [11] | [29] | 0.0362 | 32 |
| General counting | [8] | [19] | [38] | 0.0200 | 64 |
| Class counting | [8,9,11] | [79,63,38] | [101,166,206] | 0.0029 | 512 |
| General decision | [8,9,11] | [59,47,29] | [81,130,161] | 0.0021 | 512 |
| Class decision | [8,9,16] | [323,257,173] | [349,608,783] | 0.0004 | 1,048,576 |

All use kappa = 10 at every level and the existing flags. `leafvars` is true exactly for the one-level witnesses. The complete dictionaries are in `additional_certificates.json` under `model: frozen-unchanged`.

The smallest reproduction is:

```python
import copy, json
import appd_assembly as A
c = json.load(open("certificates.json"))["certificates"]["brute"]
p = copy.deepcopy(c["params"])
p.update(delta=0.0362, Lc=32)
v = A.exponent_of(p, 100)
assert v["valid"] and v["e"] < 1
print(v["e"])  # 0.9995344797946857
```

For the 1.637 target, just change the supplied certificate's delta from 0.01 to 0.02 and evaluate it at 157.

For fixed parameters the model's monotonicity argument is correct. Every size condition relaxes as lambda increases, and e stays constant. Thus these witnesses remain feasible thereafter. A scan of a heuristic parameter selector, however dense, proves neither global optimality nor the absence of an earlier feasible dictionary. The search code was not supplied, and these witnesses independently establish the issue.

**2. Step-by-step findings.**

“No problem found” means the stated identity or inequality survived this review; it is not a claim to have verified the missing implementation.

| Step and ledger location | Finding | Severity / consequence |
|---|---|---|
| `appd-tight`, 309–323 | The strict multiplicity choice must be floor(degree / total taper) + 1, as now written. The exact taper deficit can absorb inherited precision and jet evaluation loss. Exact rational checks pass for every retained valuation level of all five supplied certificates. | No problem found at the supplied certificates. A finite-size proof should state the exact inequalities instead of replacing h by an asymptotic equivalent. |
| `appd-bernstein`, 324–330 | The bipartite two-stage argument is credible with the code's 8 Delta constant. Its extension to correlated chunk endpoints is not justified merely by naming Bernstein and allocating a per-color deviation budget. A concrete replacement proof is supplied below. | Proof-completeness issue; the replacement supports the existing root condition and improves it. |
| `appd-memo`, 331–340 | A child depends on incident modes, not the full tuple. Reuse across shell profiles is valid when the node, input tensor descriptions, Taylor/evaluation data and precision are fixed. | No identity problem found. Entry-oracle and lookup costs still need accounting; see point 4. |
| `appd-modedp`, 341–351 | Every original mode has one valuation-bearing target. Partial contractions remain sums of terms with the required divisibility. Terms omitted by a weight cutoff cannot become nonzero later because no division by 2 occurs. The ordered-pair count agrees with `eps_dp`. | No problem found. |
| `appd-chunks`, 352–370 | The connected-chunk saving and M_eff <= (1+1/Lc)n are sound component by component. Nonbipartite root tensors still factor over vertices/groups. Children recover the old-vertex/packet-vertex bipartition. | No contraction problem found. Root balancing and the cost of evaluating a chunk tensor must be included. |
| `appd-perclass`, 387–394 | Padding once per unordered pattern is valid. The per-coordinate pair rounding is small enough for the 1+2 delta mode-load factor; see derivation below. | No problem found. |
| `appd-leafvars`, 395–403 | The two-product expansion of an unsatisfied designated clause with at most two external literals is valid. The code's `(1+delta)(1+M_root/2)/r` follows from c2 <= M_eff/2. | Needs simultaneous balance of the vertex weight “variables + two-external-wire clauses.” Port balance alone does not imply this. The required extra events are affordable at these sizes. |
| `appd-labelfree`, 446–455 | Nodes with nonzero multiplicity read only valuation coordinates below u. Making all higher symbols equal preserves packing and node locality. | No problem found. |
| `appd-newton`, 544–584 | The Newton identity and factorial valuation gain are correct. There is a simple finite bound using the sharp weights, stronger than the asymptotic justification in the original paragraph. | No problem found when sharp weights are used; explicit finite argument below. |
| `appd-ldu`, 585–622; finite bound 684–695 | The residue ordering makes the triangular factors integral. The full virtual block is essential. The finite weight majorization supports cK=1/2. | No problem found. Do not replace the full virtual block with an arbitrary occupied submatrix. |
| `appd-shellmm`, 647–681 | The negative exponent on the smallest dimension must be treated as a maximum of three positive-coefficient log-cost forms. Doing so gives exactly the stated squared discount `(omega-2)^2`. Rectangular padding has only constant-factor cost per product. | No problem found. |
| `appd-guarded`, 747–818 | Central real slots, complete node residue grids, odd column units and the full-block LDU fit together. Degree falls to M(Q-1). The exact integer multiplicities and inherited precision pass the supplied certificates. | No problem found. It is a new recovery/packing lemma, not a literal application of the old lemma with unchanged hypotheses. |
| `appd-bincross`, 868–964 | The integral support factorizations and fast digit transforms survive re-derivation. The channel bound is sqrt(12)^m and theta must keep the original two shell weights. | No problem found. The polynomial transform factors must not be mistaken for a valuation gain. |
| `appd-profile`, 996–1055 | The min/average/max argument is a valid upper bound. All thirteen PF_EPS constants independently reproduce exactly. Counting complete profiles is polynomial for the fixed menus used here. | No problem found. This is one of the better-supported steps, not one of the weakest. |

Other `appd-*` entries are not all active in the five scores. `appd-levels` and `appd-twolevel` require immutable creation-time representatives and inherited precision, which the ledger explicitly describes. `appd-reset` and `appd-multbal` are marked unadopted and are not flags in the frozen assembly; their reported gains or tests cannot explain any scored cell. The CRT reset requires the stated entry-height bounds and polynomially many reconstruction calls; I found no contradiction in the sketch but did not reconstruct its omitted program. The multiplicity-sensitive balance sketch does not supply its named (L) and (R) inequalities in this package, so that lemma cannot be fully audited here. `appd-density` is addressed in point 6.

**The three weakest links.** Among the fourteen scored steps, I rank (1) the chunk-root extension of `appd-bernstein`, (2) the complete work accounting surrounding `appd-memo` and recursive packet-entry evaluation, and (3) the simultaneous weighted balancing needed by `appd-leafvars`. These are weaker proof interfaces, not discovered false identities. Outside those scored steps, `appd-density` has the most concrete mathematical defect. None of the newest four steps belongs in this top three on the evidence checked here.

There is also a useful dependency distinction: at every level of every supplied round-11 certificate, PF_EPS is smaller than the binary-crossing alternative. Turning `bincross` off leaves all five exponents unchanged. Likewise, with `ldu` enabled the separate `newton` flag does not control cK. The round-11 numerical claims therefore do not require an additional runtime saving from binary crossings or the standalone Newton route.

**3. Re-derived inequalities and size conditions.**

**Profile menus.** Write E for the unordered group pairs, ell_e = log2 N_e, and alpha = log2(7)-2. A blocked Strassen step has logarithmic cost

    X + Y + Z - (1-alpha) min(X,Y,Z)
      = max(alpha X+Y+Z, X+alpha Y+Z, X+Y+alpha Z).

The coefficients are nonnegative. Therefore replacing alpha by 101/125 is safe. For an order j, let C_j(ell) = max_{a in F_j} a dot ell. With s orders,

    min_j C_j(ell)
      <= (1/s) sum_j C_j(ell)
      = max_{a_j in F_j} ((1/s) sum_j a_j) dot ell.

Using ell_e <= beta sqrt(2 w_e) and sum w_e < P, Cauchy–Schwarz gives

    min_j C_j(ell)
      <= beta sqrt(2P) max_{a_j in F_j} ||(1/s) sum_j a_j||_2.

This is the direction needed for an upper bound. There is no illicit interchange of min and max. The exact resulting constants include

    r=7: 651717/1531250
    r=8: 102163/250000
    r=9: 338143/843750
    r=10: 2440461/6250000
    r=16: 6840623/18000000.

All r=4,...,16 agree with the attachment. I also checked the integer inequality 7^125 < 2^351. Rounding each edge log-size upwards adds at most |E| to an order's envelope. Empty shells should be skipped, and the unique zero-weight tuple has alphabet size one. Neither case invalidates the bound.

**Shell matrix multiplication.** For dimensions nx, ny, nz, take S the largest power of two at most their minimum. The ceilings in the number of S-blocks each cost less than a factor two, while S is within a factor two of the minimum. Thus padding changes a product's cost by an absolute constant. If an edge family has entropy coefficient A_j, optimizing the positive-coefficient forms over a shared weight budget gives

    A_x + A_y + A_z - [1-alpha^2] min(A_x,A_y,A_z).

The square on alpha is necessary and present in the code. Reading two inputs and writing the output is also covered: for example X+Y <= X+Y+Z-(1-alpha)min(X,Y,Z).

**Binary crossing.** With j=(0,1,-1,0), I-jj^T is exactly SWAP. After choosing an integral basis whose first column is j, the coefficient of degree i has row support on words with at least i distinguished digits. The alternate flattening has support on words with at least m-i distinguished digits. These support restrictions, not an inverse of 2, give the two integral factorizations.

For the chosen smaller support,

    sum_{k <= floor(m/2)} binom(m,k) 3^k
      <= 3^(m/2) sum_{k=0}^m binom(m,k)
      = sqrt(12)^m.

Hence, including unmatched bits, log2|theta| <= max(a,b)+(c-1)min(a,b), with c=log2(12)/2. Applying Cauchy–Schwarz to coefficients 1 and c-1 gives `u_theta = 2[1+(c-1)^2]`. The index is graded by wx+wy; its piece number is not an additional 2-adic valuation. I checked both integral support factorizations for m=0,...,4 and the channel inequality for m=0,...,100. This does not rerun the missing contraction implementation.

**Finite mode counts.** Let k=a-i-1 and 0<=b<2^k. For Newton's sharp weights,

    w_b = k+2+2b-s2(b) >= 2b+2,

because s2(b)<=k. The decreasing generating-function summand can therefore be bounded by the integral on a grid of spacing two. This proves the finite density constant cK=1, without dropping an s2(b) correction asymptotically.

For LDU, the relevant weights are `k+2+4b-2s2(b)`. Their partial sums dominate those of `4b+2`, with equal totals. Applying the convex function `-log(1-exp(-s x))` and then the midpoint integral gives pi^2/(24s) per chain. I checked the partial-sum statement through k=16. Together with the argument, this supports the finite cK=1/2 used in the model.

**Guarded recovery and precision.** Put Q=2^a, H=2^(a-kappa), b=a-kappa, and

    S = sum_{i=0}^{u-1} (H/2^i)(u-i)
      = 2H(u-1+2^(-u)),
    h = floor(M(Q-1)/S)+1.

Then hS>M(Q-1), including the strict endpoint. For a retained level i, the exact cardinal-factor loss is

    loss_i = h (H/2^i) sum_{j=1}^{b-i} j min(j,u-i)/2^j.

The infinite sum is `6-(u-i+3)/2^(u-i-1)`. Consequently

    6hH/2^i - loss_i >= 4hH/2^(u-1).

Add inherited precision and the jet-evaluation loss `(h(u-i)-1)(a+1-i)`, then compare against the chosen precision. This explains the headroom test: it uses the deficit below six, rather than simply forgetting the inherited precision. The checker includes the guarded degree/capacity ratio and performs exact rational comparisons for every retained i at all five original certificates. They pass. It also checks 5,658 small-grid cardinal-factor losses directly against the sum formula.

The code's root `lp_in=0` represents n rather than n+1. The exact checks use n+1. No certificate fails. Integer wire ceilings and floor/ceiling multiplicities should nevertheless be written explicitly in a theorem; they should not be justified by floating-point equality.

**Per-pattern padding.** Put S_pat=r(r+1)/2. With per-label load at most (1+delta)M/r^(2u), the padded slots belonging to a specified unequal coordinate pair are at most

    (1+delta) 2M/r^2 + Q S_pat^(u-1).

The code requires `Q S_pat^u <= delta M`, so the second term is at most `2 delta M/[r(r+1)]`, which is below `2 delta M/r^2`. This proves the 1+2 delta pair-load allowance. Guarded packets use capacity N=Q-2H; their per-pattern rounding is at most N<=Q, while the Q/N increase in the number of packets is separately charged by gamma. Real child ports need not include dummy positions. Padding and child-port accounting should not be conflated.

**A replacement proof for root label balance.** This also supplies a useful improvement. Regard the oriented crossing wires as the edges of a loopless undirected multigraph for coloring purposes. Parallel and oppositely oriented wires are allowed. The chunk-size bound gives maximum incident degree Delta<=63Lc. Greedy edge coloring partitions the wires into q<=2Delta-1<=126Lc-1 matchings.

Assign each vertex an independent uniform label in a set of size K=r^u, and let D=K^2. Fix a matching j and an ordered label pair gamma. Its load X_j is binomial with mean mu_j=|E_j|/D: different matching edges have disjoint endpoints. Put mu=M/D and

    t_j = (delta/2)(mu_j+mu/q).

If all these events satisfy |X_j-mu_j|<=t_j, then summing over matchings gives |X_gamma-mu|<=delta mu. A Bernstein calculation already gives a tail bounded by `2 exp(-delta^2 mu/(3q))` for delta<=1/4. For an explicitly rational pessimistic estimator one can use z=1+delta/2 and the products for z^X_j and z^(-X_j). The weaker bound

    Pr[bad matching/label event] <= 2 exp(-delta^2 mu/(5q))

suffices. To see the upper-tail estimate, expand its logarithmic moment; the coefficient of mu_j is nonpositive because `(1+x)log(1+x)>=x`, where x=delta/2. The remaining exponent is at most `-delta mu log(1+delta/2)/(2q) <= -delta^2 mu/(5q)`. The lower tail is bounded the same way, using `exp(-t)-1+t<=t^2/2` and t=log(1+x)<=x.

Thus the explicit sufficient condition

    M/D >= (8q/delta^2) log(8qD)

leaves the union bound over 2qD tails below 1/4. These estimates can be derandomized: within a matching the conditional moments factor over edges, and their sum is an explicitly computable pessimistic estimator as vertex labels are fixed. Integer event thresholds make the powers of the rational z rational. The number of events, candidate labels, and arithmetic bit lengths is polynomial in the network and label-space sizes.

One also includes the weighted vertex events required by `leafvars`, with weight “variables plus designated clauses with two external wires.” Each weight is at most 28Lc. Their standard weighted Chernoff bound has ample slack under the displayed root condition; row/column port bounds already follow by summing the pair-label upper bounds. This is why wire-load balance alone is not the entire `leafvars` proof, but the extra requirement is repairable.

For Lc>=1 and the model's r,u ranges, the original `8192 Lc^2 log(172LcD)/delta^2` condition dominates this matching condition. Therefore I do **not** conclude that the existing root balance inequality is false. The original color-stage sketch is incomplete, but a different explicit argument supports the inequality and allows it to be improved.

Replacing only the root balance line with this sufficient matching condition gives the following additional witnesses, retaining the other cost and size formulas:

| Target | New feasible log2 n with matching balance | e |
|---|---:|---:|
| 2^n | **92** | 0.998878510359880 |
| 1.637^n | **148** | 0.710254458410584 |
| 2^(n/3) | **526** | 0.333312175935314 |
| 1.307031578^n | **403** | 0.386239290207910 |
| 1.1199^n | **2,033** | 0.163366221426409 |

These are **proposed balance-extension witnesses**, not outputs accepted by the unchanged frozen model at those sizes. They have a separate model tag in the JSON and separate size checks in the checker. They are not claimed optimal, and the matching-label construction has not been integrated into the omitted SAT implementation. In particular, the profile-order certificates and all contraction algebra are unchanged.

**4. Is the exponent an upper bound on the algorithm's whole work?**

Severity: **a gap in the certification of the full finite-size runtime**, not a discovered contradiction to the five leading exponents.

Complete profiles are not missing an exponential factor when r is fixed. At a node they number at most `(p+1)^(r(r-1)/2)`. At the five supplied witnesses, r<=16 and log2 p is only a few thousand. A per-profile bound must be multiplied by this polynomial, not by an exponential in n. The word “polynomial” must be qualified if r is later allowed to grow with n; the frozen profile menus only cover r<=16.

The blocked matrix-product cost is also safe up to constants, and its input/output traffic is covered in ring-operation count. Indexed access to exponentially many tables requires an explicit machine model: addresses have O(n) bits, and a suitable indexed-memory or sorting implementation contributes polynomial overhead. An implementation that scans an entire memo table for each lookup would not have the stated exponent. The attachments do not include the implementation to check which occurs.

The more significant omitted term is **implicit tensor-entry evaluation**. The paper's D.4.6 construction describes an entry algorithm with `exp(O(Q log s))` work, not a uniformly bounded-degree polynomial in n when Q grows. Chunk tensor entries also have a `2^O(Lc)` construction cost. For a fixed certificate these are constants or parameter-dependent polynomial factors, but a uniform finite-size theorem must explicitly bound them as parameters change. The assembly does not do that.

For the supplied witnesses this omission does not look large enough to overturn a cell: their largest packet exponents are 29, 38, 218, 170 and 783, versus lambda 103, 159, 533, 413 and 2,034, respectively. A packet-entry exponent of order `2^a (lambda+a)` is vastly smaller than the available `n(target-e)` margin. That is a quantitative scale check, not an independently counted implementation bound. The report does not assert an unspecified polynomial “cannot move a cell” merely because n is large.

A clean statement should bound the cost as

    polynomial overhead × entry-evaluation overhead
      × 2^(n max(worst, cum+leaf)),

with explicit parameter dependence. It should also specify how unusually small child networks are handled when they do not themselves meet a lower bound on M: an early direct-contraction rule or a justified padding construction. Testing a worst-case upper bound on M does not alone establish a lower-bound size condition for every smaller child.

The memo recurrence itself is sound. Child tables are shared among profiles at one fixed node; they are not asserted to be the same across unrelated parent tensor descriptions or higher precisions. The cumulative child-table exponent pays for that distinction. Canonical lifts must remain attached to the tensor's creation precision; reducing and then reevaluating a different formula at a higher precision is not interchangeable with preserving those lifts.

Finally, the claim that the root fixes the *size* is different from saying its arithmetic is the largest *work* term. In the supplied `third` certificate the three arithmetic terms are approximately 0.3296866444, **0.3330366656**, 0.3284031091, and the leaf total is 0.3323297786. In `occ3sat`, the leaf total **0.1633438831** is largest. The root balance is the binding size condition in four cells, as reported; it is not the largest arithmetic term in every cell.

**5. Competitors and the paper's correctness interfaces.**

The 1.637^n deterministic general-counting benchmark is supported by Kutzkov's IJCAI 2026 paper. The general-decision number actually used in the code, 1.307031578^n, agrees with Jiang and Cai's current PPSZ preprint. The 1.1199^n bound is a decision algorithm for SAT with at most three occurrences, and applies to the width-at-most-three subclass. I found no primary-source improvement that invalidates those particular exponential benchmarks. They should be called the conventional exponential benchmarks: if the assumed subexponential theorems are accepted, they are no longer literally the fastest asymptotic algorithms.

For the counting benchmark, the incidence graph has maximum degree three. Degree-three variable vertices number at most n and degree-three clause vertices at most n. After handling low-degree attachments/subdivisions, the cubic pathwidth bound gives the leading rate n/3; a path DP stores variable values and clause-satisfaction bits. The safe notation is **2^(n/3+o(n))**, or a specified epsilon version, rather than claiming the cited theorem alone gives a sharp finite `2^(n/3) poly(n)` runtime at an identified n. Fomin–Høie's Theorem 5 is stated in terms of all graph vertices; the degree-three-vertex refinement needs the usual reduction argument, not that theorem number alone. I did not establish that this counting benchmark is an optimal or exhaustive current record for the intersection of the two restrictions. The general #3-SAT clause-count bound 1.4423^m does not replace it with a uniformly smaller rate for this class.

Define “3-CNF” explicitly as **clauses of width at most three**, after removing repeated literals. If instead every clause has exactly three distinct variables and each variable occurs at most three times, all such formulas are satisfiable. Indeed any set of s clauses has 3s incidences and therefore at least s neighboring variables; Hall's theorem matches every clause to a different variable, which can satisfy its matched literal. On that stricter class the decision comparison is with a polynomial algorithm. Counting remains nontrivial. This is a class-definition issue, not a criticism of the intended width-at-most-three class.

Definition 3.1 and Lemmas 4.5 and 4.7 govern the paper's main replacement-form construction. They are **not the interface that every Appendix D modification must preserve**. The correct obligations here are packet placement and real-port balance; F(0) equalling the desired contraction; integral jets at the selected nodes; strict interpolation degree coverage; sufficient recovery precision; integral, local source/target modes; and polynomially controlled immutable tensor descriptions. The report's independent derivations address those obligations. Reusing the paper's theorem name cannot replace proofs of the modified obligations.

The guarded deformation does not satisfy the old literal formula for W, and the proposed r,u,a do not satisfy the paper's original fixed parameter schedule or `u>=lambda^(4/5)` condition. That is intentional. The latter condition was a sufficient route to precision domination in the original analysis, not an intrinsic condition for Hermite correctness. Explicit finite recovery inequalities can replace it, as checked above.

**6. General formulas: a missing component term and a limited comparison.**

Severity: **breaks the stated global density inequality**; whether and how much it changes the listed high-density cells requires their missing patch and certificates. Locations: ledger lines 727–735 and 1091–1097.

For a given connected chunk partition, the actual combinatorial statement is

    M_eff <= sum_f (|f|-1) - n + U,

where U is the number of chunks, including one chunk for each small connected component. One may replace U by n/Lc only when every retained component is large enough for all its chunks to have at least Lc variables. Small components contribute zero crossing wires, but do not satisfy U<=n/Lc.

A minimal counterexample to the globally written bound is one clause on three variables with Lc=4. It forms one small component, so M_eff=0, while

    2m - n + n/Lc = 2 - 3 + 3/4 = -1/4.

Thus `0 <= -1/4` is false. No isolated variables or repeated clauses are needed. This counterexample does not by itself show that the m=4.26n or m=5n numerical witnesses fail; it shows that the universal inequality used to obtain them needs a correction. At fixed high total density, mixtures of small sparse components and a dense large component can also invalidate simply replacing the active component's variable count by the original n.

Taking any number of disjoint copies gives the same failure at arbitrarily large n, so this is not an objection restricted to inputs below the model's lambda=64 cutoff.

A safe repair is to solve small components first, track n_large and the clauses/factors in the large components, and use

    M_eff <= sum_{f in large components} (|f|-1)
             - n_large + n_large/Lc.

Alternatively retain U explicitly. The comparator exponent is still expressed in the original input n; do not silently renormalize it to n_large. After this repair, rerun and publish the density certificates.

The copying idea itself is sound **when an equality is represented as one arity-two tensor constraint**. Splitting a high-degree variable into k copies adds k-1 variables and k-1 equality constraints, so `sum(|f|-1)-n` is unchanged. It is parsimonious. A chain with internal copies using at most B-2 original incidences gives a number of added variables below the total original occurrence count divided by B-2. If equality is instead encoded by two CNF clauses, the bookkeeping must be redone; the asserted cancellation does not follow from counting that encoding as one clause.

The B^4 scaling of the old chunk-balancing constant is plausible from unit-size and conflict-degree growth, but a scaling slogan is not an exact constant certificate. The logarithmic event-count term and simultaneous weighted events also need updating. The attached assembly has no density parameter or occurrence-cap parameter, so it cannot verify the six general-density numbers by itself.

I checked the Anthropic statement in its own PDF, rather than relying on lilting.ch. Its introduction does state the `n>10^10000` requirement for its proved bound at m=5n. It is phrased as a **necessary threshold for that displayed bound to improve on exhaustive search**, not a certified exact crossover and not a lower bound on every algorithm in that family. Comparing it to a conditionally sufficient witness for a heavily modified different route is meaningful if those statuses are stated. It is not a fair basis for saying one implementation or one underlying algorithm is intrinsically that much faster. The more recent crossover-analysis preprint also makes clear that the original Anthropic constants are not the only analysis available.

**7. Reproduction, remaining limits, and recommended wording.**

The audit archive contains the unchanged assembly and supplied certificates in `upload/`, the new certificate dictionaries, `review_checks.py`, and its output JSON. Run:

```bash
python3 review_checks.py
```

NumPy is needed for the batched exact integer profile checks and the small binary matrices. The core exponent evaluations use the frozen attachment. The checker does not claim to execute the missing SAT implementation. `local_search.py` records the exploratory search and additionally requires SciPy; the certificate checks do not need SciPy or a search rerun.

Completed checks include all five original numerical certificates; all thirteen rational PF_EPS constants; 5,658 exact small-grid Hermite-loss checks; every retained precision level in the original and additional witnesses; 17 finite LDU-majorization sizes; 15 binary-factorization pieces; and 101 channel-bound cases. No failed original certificate or failed profile constant was found.

I would replace “the smallest size from which the route stays ahead” with: **“A feasible threshold for this conditional leading-exponent model, found by the stated parameter search; not proved minimal. Polynomial and tensor-entry costs require explicit accounting for a full finite-runtime comparison.”** Keep the implementation-test status, cost-certificate status, extension-proof status and assumption of the paper's theorem separate, as the brief already tries to do.

The two immediate corrections are to accept the smaller unchanged-model witnesses and to repair the component accounting before presenting general-density results as certified. The strongest additional mathematical opportunity found in this review is the matching-based label-balance proof, with its separately reported witnesses.

**Primary sources consulted.** The external links support the literature/interface statements; the mathematical audits of the proposed extensions above are independent derivations from the supplied ledger.

- OpenAI, *The Exponential Time Hypothesis is False*, 25 October 2026 revision: https://github.com/openai/math/tree/main/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026 ; contraction source `build/quantitative-contraction.tex`.
- Fomin and Høie, *Pathwidth of cubic graphs and exact algorithms*, IPL 97 (2006), 191–196, especially Theorem 5: https://fedorvf.github.io/articles/2006/2006b.pdf .
- Jain, Neoh and Stephan, *An Improved Algorithm for Sparse Instances of SAT*: https://arxiv.org/abs/2411.07389 .
- Kutzkov, *Complementary Branching for #SAT: Bounded-Degree Bounds and Faster #3-SAT*, IJCAI 2026: https://ijcai.org/proceedings/2026/254 .
- Jiang and Cai, *A Better Analysis For PPSZ For 3-SAT*, current v2: https://arxiv.org/abs/2607.10697 .
- Peng, Sheng and Xiao, *New Algorithms for #2-SAT and #3-SAT*, IJCAI 2025: https://ijcai.org/proceedings/2025/297 .
- Anthropic/Claude, *A 2^{o(n)} Time Algorithm for 3-SAT*: https://github.com/anthropics/formal-math/blob/main/eth/paper.pdf , introduction and Appendix A.
- Majdoddin, *Lowered crossover points of Anthropic's subexponential 3-SAT algorithm*, 30 October 2026 preprint: https://arxiv.org/abs/2610.41622 . This review does not independently validate that preprint's improvements.
