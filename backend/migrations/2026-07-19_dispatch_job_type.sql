-- job_type on the dispatch queue — routes a queued job to the right pull-agent:
--   'hpc' (default) → the classical Laguna agent (solange_hpc.py --agent)
--   'qpu'           → the real-hardware agent (solange_qpu.py --agent)
-- Existing rows backfill to 'hpc' so the classical agent keeps claiming them.
-- The QPU backend/shots are agent-level choices (set when you start the QPU
-- agent), so the queue row itself needs only key/side + this type.
-- Run once in the Supabase SQL editor.

alter table public.hpc_dispatch
  add column if not exists job_type text not null default 'hpc';

update public.hpc_dispatch set job_type = 'hpc' where job_type is null;
