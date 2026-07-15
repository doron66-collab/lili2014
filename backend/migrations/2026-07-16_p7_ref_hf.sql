-- p7_ref_hf_ha — the exact classical <H> on the SAME fixed HF reference state a
-- Phase 3B QPU run measured. For a fixed-state QPU smoke test, the honest measure
-- of hardware fidelity is |measured - exact_HF|, NOT |measured - CASSCF_ground|
-- (the latter conflates real hardware noise with the intended, deliberate
-- HF-vs-ground gap). Storing this lets the HPC Ladder show the true hardware
-- noise for QPU rows instead of a misleadingly large ΔE.
-- Named with the p7_ prefix so LEON's existing P8 seal already covers it.
-- Run once in the Supabase SQL editor.

alter table public.simulation_runs
  add column if not exists p7_ref_hf_ha double precision;
