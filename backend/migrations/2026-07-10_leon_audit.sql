-- LEON audit trail — a persistent, append-only, time-stamped ledger of every
-- notarization, re-verification, and rejection LEON performs. This is the running-code
-- backing for the 21 CFR §11.10(e) audit-trail claim (§06.iii): audit records are
-- computer-generated, time-stamped, and never mutated in place.
-- Run once in the Supabase SQL editor.

create table if not exists public.leon_audit (
  id            uuid primary key default gen_random_uuid(),
  created_at    timestamptz not null default now(),
  event         text not null,              -- notarize | reverify | reject
  run_id        text,                       -- simulation_runs.id (when known)
  integrity     text,                       -- PASS | FAIL | REJECTED | LEGACY-UNVERIFIABLE
  seal_ok       boolean,
  consistency_ok boolean,
  method        text,                       -- sealed-payload | legacy-reconstruction | ingestion
  stored_hash   text,
  recomputed_hash text,
  actor         text,                       -- provenance_source / requested_by, when available
  note          text
);

create index if not exists leon_audit_created_idx on public.leon_audit (created_at desc);
create index if not exists leon_audit_run_idx     on public.leon_audit (run_id);

alter table public.leon_audit enable row level security;

-- Append-only by policy: allow inserts and reads; no update/delete policy is created,
-- so even the anon/service paths cannot rewrite history through this table's policies.
do $$
begin
  if not exists (select 1 from pg_policies where tablename='leon_audit' and policyname='leon_audit_insert') then
    create policy leon_audit_insert on public.leon_audit for insert with check (true);
  end if;
  if not exists (select 1 from pg_policies where tablename='leon_audit' and policyname='leon_audit_read') then
    create policy leon_audit_read on public.leon_audit for select using (true);
  end if;
end $$;
