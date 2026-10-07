`fpocket` — precompiled x86_64 Linux binary, GPL-licensed, source:
https://github.com/Discngine/fpocket

Vendored here (rather than relying solely on compiling it during Render's
build) because Render's Python build image has no confirmed gcc/make/git
toolchain, and that's not something this environment can verify directly
against a live Render deploy. `render.yaml`'s buildCommand still attempts a
from-source compile as a fallback/refresh path if this binary is ever
missing or non-executable; this committed copy is what actually ships on
the first deploy regardless of whether that compile step works.

Built from a clean `git clone --depth 1` + `make` of the upstream repo,
then `strip`ped (debug symbols only — GPL §3 source-availability obligation
is met by the upstream link above, not by anything removed here).
Reproduce with:
    git clone --depth 1 https://github.com/Discngine/fpocket.git
    cd fpocket && make && strip bin/fpocket
