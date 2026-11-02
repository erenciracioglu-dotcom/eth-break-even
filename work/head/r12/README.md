# Round 12 (head)

Round 12 changed only the cost model (the matching-based root balance, `appd-matchbal`, as referee version 8) and the parameter search; it had no coding task.

Runnable from the repository root:

| File | What it does |
| --- | --- |
| `round12_witnesses.json` | The five certificates behind the cells of the results table; `python check_certificates.py` re-checks them |
| `whatif_dens_v8.py` | The general-formula entries: `python work/head/r12/whatif_dens_v8.py 4.26 out.json brute,sharp3sat,ppsz` (or 5.0) |
| `scan12.py` | The referee's threshold search on a code folder with a given model: `python work/head/r12/scan12.py algo_appd referee_appd/appd_assembly.py out.json` |

Records of how round 12 was built, validated and installed; they refer to the owner's folders (`D:/ETH-loop`) and to staging copies that are not in this repository: `patch_v8.py` (the model patch), `patch_params_v8.py` and `patch_params_v8b.py` (the search), `validate_v8.sh` (the validation runs: head solution, untouched round-11 code, a deliberately broken variant), `install_r12.py` (the installation) and `check_astra_witnesses.py` (the check of the review's smaller parameters).
