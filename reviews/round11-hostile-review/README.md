# Hostile review of round 11 (GPT-6 Astra, 2 November 2026)

After round 11 the owner asked GPT-6 Astra to try to refute the computation. This folder holds what was sent and what came back.

| File | What it is |
| --- | --- |
| `review_prompt.md` | Our request and the questions asked |
| `upload/` | The cost model as reviewed (version 7, SHA-256 `cf422a0e…6748a8c7`) and the round-11 certificates |
| `appendix_d_hostile_review.md` | Her report |
| `review_checks.py`, `review_checks.json` | Her independent arithmetic checks and their output |
| `local_search.py`, `search_results.json`, `additional_certificates.json` | Her parameter search and the smaller valid parameters it found |
| `audit_manifest.json` | SHA-256 of her files |

Run her checker from this folder (it needs numpy; `local_search.py` also needs scipy):

```
python review_checks.py --assembly upload/appd_assembly.py --certificates upload/certificates.json --output my_checks.json
```

What changed because of the review is in round 12: the matching-based root balance, the finer parameter search that took in her smaller parameters, and the repaired general-density entries. Our answer to each point is in `state/claims.md`.
