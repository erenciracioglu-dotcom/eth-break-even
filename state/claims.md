# Claims ledger

PROVISIONAL = accepted by the referee under a head extension that is not yet audited.
AUDITED = the head checked the argument. REJECTED = the argument failed and results built on
it are reverted.

## Regimes

- paper [BASE]: Lemma 4.5 of the paper as written: 2(b+w) D^2 <= M, h <= M/D + sqrt(2(b+w) M).
  Accounting refinements used by the frozen assembly (rows touch <= 2 groups, restricted row
  abs sum <= 2^(Q-1), replacement-row degree <= R-1, n/r groups) are checked on real designs by
  the measure gate on every referee run.
- freedman [PROVISIONAL, head 2026-10-31]: D ln D <= M / (4(b+w)) with h <= 2M/D, from a
  Freedman/Bernstein bound per label, a union bound over the D labels and derandomisation by
  pessimistic estimators. Needs a written audit against the exact dependency structure of
  Lemma 4.5 before any result that uses it is reported as more than provisional.

## Results that depend on a provisional regime
