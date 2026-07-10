-- Agent liveness — a single-row table the cluster pull-agent upserts on every poll,
-- so SOLANGE can show whether an execution agent is currently active (online/offline).
-- Run once in the Supabase SQL editor.

create table if not exists public.agent_heartbeat (
  id         text primary key default 'default',
  last_seen  timestamptz not null default now(),
  agent      text,
  note       text
);

alter table public.agent_heartbeat enable row level security;

do $$
begin
  if not exists (select 1 from pg_policies where tablename='agent_heartbeat' and policyname='agent_hb_all') then
    create policy agent_hb_all on public.agent_heartbeat for all using (true) with check (true);
  end if;
end $$;
