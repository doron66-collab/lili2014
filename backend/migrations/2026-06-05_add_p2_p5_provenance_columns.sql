-- ─────────────────────────────────────────────────────────────────────────────
-- Migration: add 4 missing provenance columns to simulation_runs
-- Date: 2026-06-05
-- Reason: the P1–P9 record assembled in backend/routes/simulate.py grew to include
--         model-compound identity (P2) and the energy decomposition (P5: core,
--         active-space, and CASSCF reference energies). These columns were absent
--         from the live table, causing every insert to fail with PGRST204 and
--         silently dropping all new runs from the Live Results table.
--
-- Run this once in the Supabase SQL Editor (service_role / dashboard).
-- Safe to re-run: every statement uses IF NOT EXISTS.
-- ─────────────────────────────────────────────────────────────────────────────

alter table public.simulation_runs add column if not exists p2_model_compound   text;
alter table public.simulation_runs add column if not exists p5_ecore_ha          float;
alter table public.simulation_runs add column if not exists p5_active_energy_ha  float;
alter table public.simulation_runs add column if not exists p5_casscf_ref_ha     float;

-- Force PostgREST to reload its schema cache so the new columns are visible
-- to the API immediately (otherwise inserts keep failing until the next reload).
notify pgrst, 'reload schema';
