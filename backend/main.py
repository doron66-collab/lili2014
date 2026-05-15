from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from routes import simulate, pdb, provenance

app = FastAPI(
    title="QC·AI·HPC Simulation API",
    description="Quantum simulation backend for NSCLC non-druggable mutation targets — IST 697 · Doron Cohen · CGU 2026",
    version="0.1.0"
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(simulate.router,   prefix="/api/simulate",   tags=["Simulation"])
app.include_router(pdb.router,        prefix="/api/pdb",        tags=["PDB"])
app.include_router(provenance.router, prefix="/api/provenance", tags=["Provenance"])

@app.get("/")
def root():
    return {
        "status": "online",
        "service": "QC·AI·HPC Simulation API",
        "version": "0.1.0",
        "phase": "3A — PennyLane simulator backend"
    }

@app.get("/health")
def health():
    return {"status": "ok"}
