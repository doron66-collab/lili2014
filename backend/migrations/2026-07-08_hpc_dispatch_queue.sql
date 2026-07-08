-- HPC outbound dispatch queue: SOLANGE → (pull agent on the cluster) → SOLANGE.
-- SOLANGE queues a job here; a solange_hpc.py --agent running inside the user's
-- Laguna session pulls it, runs it, and --submits the result back.
-- Run this once in the Supabase SQL editor.

create table if not exists public.hpc_dispatch (
  id            uuid primary key default gen_random_uuid(),
  created_at    timestamptz not null default now(),
  requested_by  uuid,
  status        text not null default 'queued',   -- queued | running | done | failed
  key           text,          -- e.g. ARID2_LOF
  side          text default 'native',
  compound      text,          -- model compound (GEOM key), e.g. acetamide
  basis         text default '6-31g',
  ncas          int,
  nelecas       int,
  run_vqe       boolean default false,
  residue       text,
  note          text,
  claimed_at    timestamptz,
  finished_at   timestamptz,
  run_id        text           -- id of the resulting simulation_runs row, once submitted
);

create index if not exists hpc_dispatch_status_idx on public.hpc_dispatch (status, created_at);

-- The backend uses the service_role key, which bypasses RLS. Enable RLS so the
-- table is not world-readable via the anon key; the backend endpoints mediate access.
alter table public.hpc_dispatch enable row level security;
