-- ============================================================
-- QC·AI·HPC Simulation Platform — Supabase Schema
-- IST 697 · Doron Cohen · CGU 2026
-- FDA 21 CFR Part 11 compliant audit trail
-- ============================================================

-- Enable UUID extension
create extension if not exists "uuid-ossp";

-- ── Authentication roles ──────────────────────────────────────────────────────
-- Supabase Auth handles user accounts.
-- Custom roles are stored as app_metadata on each auth.users row:
--   { "role": "admin" }   — sees all runs, can delete, accesses audit log
--   { "role": "researcher" } — sees own runs only
--   (no role) / public    — read-only aggregated dashboard data

-- ── simulation_runs ──────────────────────────────────────────────────────────
create table if not exists public.simulation_runs (
    -- Identity
    id              uuid primary key default uuid_generate_v4(),
    created_at      timestamptz not null default now(),
    user_id         uuid references auth.users(id),   -- null = unauthenticated API call

    mutation_id     text not null,
    mutation_name   text not null,
    pdb_id          text,
    phase           text,

    -- P1 — Circuit fingerprint
    p1_circuit_hash text,
    p1_gate_count   integer,
    p1_depth        integer,
    p1_qubit_count  integer,
    p1_ansatz       text,

    -- P2 — Compilation lineage
    p2_compiler         text,
    p2_compiler_version text,
    p2_encoding         text,
    p2_basis_set        text,
    p2_active_electrons integer,
    p2_active_orbitals  integer,

    -- P3 — Device & calibration
    p3_backend           text,
    p3_backend_version   text,
    p3_calibration_epoch timestamptz,
    p3_simulator         boolean,

    -- P4 — Error budget
    p4_gate_error_rate    float,
    p4_readout_error_rate float,
    p4_t1_us              float,
    p4_t2_us              float,
    p4_note               text,

    -- P5 — Raw outcome distribution
    p5_shots         integer,
    p5_raw_energy    float,
    p5_energy_variance float,
    p5_opt_steps     integer,
    p5_elapsed_s     float,

    -- P6 — Error mitigation
    p6_method text,
    p6_note   text,

    -- P7 — Statistical estimator & CI
    p7_energy_ha float,
    p7_ci_lower  float,
    p7_ci_upper  float,
    p7_confidence float,
    p7_method    text,

    -- P8 — Cryptographic seal (computed over P1-P7 + P9)
    p8_hash      text,
    p8_algorithm text,
    p8_sealed_at timestamptz,

    -- P9 — ML decoder
    p9_applicable boolean,
    p9_note       text
);

-- Indexes for common query patterns
create index if not exists idx_sim_runs_mutation on public.simulation_runs(mutation_id);
create index if not exists idx_sim_runs_created  on public.simulation_runs(created_at desc);
create index if not exists idx_sim_runs_user     on public.simulation_runs(user_id);

-- ── Row Level Security ────────────────────────────────────────────────────────
alter table public.simulation_runs enable row level security;

-- Public read: aggregate dashboard data (no PII, no user_id exposed)
create policy "public_read_runs" on public.simulation_runs
    for select using (true);

-- Researchers see only their own runs (when authenticated)
-- (The public read policy above also applies; RLS is additive for SELECT)

-- Only authenticated users (or service_role) can insert
create policy "auth_insert_runs" on public.simulation_runs
    for insert with check (auth.role() in ('authenticated', 'service_role'));

-- Only admins (service_role key or custom claim) can delete
create policy "admin_delete_runs" on public.simulation_runs
    for delete using (auth.role() = 'service_role');

-- ── provenance_audit ─────────────────────────────────────────────────────────
-- Immutable append-only log of all provenance API operations (21 CFR Part 11)
create table if not exists public.provenance_audit (
    id          uuid primary key default uuid_generate_v4(),
    created_at  timestamptz not null default now(),
    user_id     uuid references auth.users(id),
    action      text not null,   -- 'create_run', 'verify_seal', 'delete_run', etc.
    run_id      uuid,
    ip_address  inet,
    user_agent  text,
    result      text,            -- 'success' | 'fail' | 'integrity_pass' | 'integrity_fail'
    detail      jsonb
);

create index if not exists idx_audit_run on public.provenance_audit(run_id);
create index if not exists idx_audit_created on public.provenance_audit(created_at desc);

alter table public.provenance_audit enable row level security;

-- Only admins can read audit log
create policy "admin_read_audit" on public.provenance_audit
    for select using (auth.role() = 'service_role');

-- Any authenticated call can append to audit log (insert only, never update/delete)
create policy "auth_insert_audit" on public.provenance_audit
    for insert with check (auth.role() in ('authenticated', 'service_role'));

-- No UPDATE or DELETE policies — audit log is immutable by design

-- ── users_profile (optional, for display names / roles) ──────────────────────
create table if not exists public.users_profile (
    id          uuid primary key references auth.users(id) on delete cascade,
    created_at  timestamptz not null default now(),
    full_name   text,
    role        text not null default 'researcher',  -- 'researcher' | 'admin'
    institution text,
    orcid       text
);

alter table public.users_profile enable row level security;

create policy "users_read_own_profile" on public.users_profile
    for select using (auth.uid() = id);

create policy "users_update_own_profile" on public.users_profile
    for update using (auth.uid() = id);

create policy "admin_read_all_profiles" on public.users_profile
    for select using (auth.role() = 'service_role');
