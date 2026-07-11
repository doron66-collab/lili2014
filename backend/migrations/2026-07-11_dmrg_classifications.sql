-- DMRG A/B/C classification records — LEON-notarized the same way as HPC/CASSCF
-- runs, but a different shape (classification + entanglement diagnostic, not a
-- P1-P9 energy record). Storing these here means a classification is safe in
-- SOLANGE the instant solange_dmrg.py --submit succeeds, independent of any later
-- loss of connectivity to the HPC cluster it ran on.
-- Run once in the Supabase SQL editor.

create table if not exists public.dmrg_classifications (
  id                   uuid primary key default gen_random_uuid(),
  created_at           timestamptz not null default now(),
  key                  text,
  compound             text,
  basis                text,
  ncas                 integer,
  nelecas              integer,
  e_casscf             double precision,
  dmrg_energies        jsonb,
  s_max                double precision,
  bqp_class            text,          -- 'A' | 'B' | 'C'
  class_rationale      text,
  time_budget_hit      boolean default false,
  bond_dims_requested  jsonb,
  method               text,
  dmrg_seal_payload    text,          -- exact canonical JSON LEON's seal hashes over
  dmrg_hash            text,          -- SHA-256 seal, computed at source, re-verified by LEON
  provenance_source    text
);

create index if not exists dmrg_classifications_created_idx on public.dmrg_classifications (created_at desc);
create index if not exists dmrg_classifications_key_idx     on public.dmrg_classifications (key);

alter table public.dmrg_classifications enable row level security;

do $$
begin
  if not exists (select 1 from pg_policies where tablename='dmrg_classifications' and policyname='dmrg_class_read') then
    create policy dmrg_class_read on public.dmrg_classifications for select using (true);
  end if;
  if not exists (select 1 from pg_policies where tablename='dmrg_classifications' and policyname='dmrg_class_insert') then
    create policy dmrg_class_insert on public.dmrg_classifications for insert with check (true);
  end if;
end $$;
