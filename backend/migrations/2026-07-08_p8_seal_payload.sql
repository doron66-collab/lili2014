-- Store the exact JSON string that the P8 seal hashes over, so re-verification is
-- ROBUST: text survives the Postgres round-trip byte-identically, unlike floats and
-- timestamps (which get reformatted and made every "Verify" spuriously FAIL).
-- Run once in the Supabase SQL editor.

alter table public.simulation_runs
  add column if not exists p8_seal_payload text;
