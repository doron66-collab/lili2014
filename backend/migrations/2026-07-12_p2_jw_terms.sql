-- p2_jw_terms — the full Jordan-Wigner Hamiltonian (every Pauli term + coefficient)
-- as computed at source. Closes a real gap: P2 ("compilation lineage") already
-- promises "source Hamiltonian" in the P1-P9 schema, but the term list was only
-- ever written to a local jw_*.json file on the compute cluster and discarded by
-- the backend after a one-time consistency check — never part of the permanent,
-- LEON-notarized record. Because this column is named with the p2_ prefix, LEON's
-- existing seal (build_p8_payload filters by p1_..p9_) automatically covers it —
-- no leon.py change needed; tampering with the Hamiltonian is now detectable too.
-- Run once in the Supabase SQL editor.

alter table public.simulation_runs
  add column if not exists p2_jw_terms jsonb;
