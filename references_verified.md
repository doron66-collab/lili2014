# Verified Literature Review — Classical↔Quantum Architectural Foundations
### Deep-research output · for dissertation § 03 (Literature Review) and § 04 (Theoretical Framework)
*Compiled 2026-06-05 · 30 sources, organized under five headings · all identifiers cross-corroborated*

> **Verification caveat (applies to every entry).** DOIs/identifiers were confirmed by
> agreement across multiple independent indexes (publisher pages, NASA ADS, PubMed/PMC,
> APS/AIP exports, OSTI, W3C, govinfo/eCFR), **not** by live byte-level resolution — the
> research environment blocked outbound fetches (uniform HTTP 403). Every DOI follows the
> correct registrant pattern and matched title/authors/year/volume across sources. Before
> final submission, one-click verify each with `curl -LI https://doi.org/<DOI>`. Items with
> residual uncertainty are explicitly **FLAGGED**.

---

## 3.1 · Physical foundations of quantum simulation

1. **Feynman, R. P. (1982).** "Simulating Physics with Computers." *Int. J. Theor. Phys.* 21(6/7), 467–488. DOI: `10.1007/BF02650179`.
   — The originating proposal: classical simulation of quantum systems is intractable; a quantum machine can do it. Conceptual root of the classical→quantum thesis.

2. **Lloyd, S. (1996).** "Universal Quantum Simulators." *Science* 273(5278), 1073–1078. DOI: `10.1126/science.273.5278.1073`.
   — Constructive proof that a universal quantum computer can efficiently simulate any local quantum system via Trotterized evolution. Legitimizes the Trotter/sqDRIFT framing.

3. **Abrams, D. S., & Lloyd, S. (1999).** "Quantum Algorithm Providing Exponential Speed Increase for Finding Eigenvalues and Eigenvectors." *Phys. Rev. Lett.* 83(24), 5162–5165. DOI: `10.1103/PhysRevLett.83.5162`.
   — Phase-estimation eigenvalue algorithm for fermionic Hamiltonians; the bridge from "simulate dynamics" to "extract ground-state energy." Estimates 50–100 qubits — directly germane to the C275F 88-qubit target.

4. **Aspuru-Guzik, A., Dutoi, A. D., Love, P. J., & Head-Gordon, M. (2005).** "Simulated Quantum Computation of Molecular Energies." *Science* 309(5741), 1704–1707. DOI: `10.1126/science.1113479`.
   — Foundational chemistry application: maps phase estimation onto real molecules. Establishes the polynomial-scaling quantum route the classical proxy emulates.

5. **Cao, Y., Romero, J., Olson, J. P., … Aspuru-Guzik, A. (2019).** "Quantum Chemistry in the Age of Quantum Computing." *Chem. Rev.* 119(19), 10856–10915. DOI: `10.1021/acs.chemrev.8b00803`.
   — Comprehensive modern review tying JW mappings, VQE, and phase estimation together. Best single citation to situate the 4-qubit JW/VQE proxy.

6. **McArdle, S., Endo, S., Aspuru-Guzik, A., Benjamin, S. C., & Yuan, X. (2020).** "Quantum Computational Chemistry." *Rev. Mod. Phys.* 92(1), 015003. DOI: `10.1103/RevModPhys.92.015003`.
   — Authoritative RMP review of NISQ-era quantum chemistry; supports the Phase 3A→3B narrative (classical proxy → hardware with noise characterization).

---

## 3.2 · Hybrid classical–quantum algorithms

7. **Peruzzo, A., McClean, J., Shadbolt, P., … O'Brien, J. L. (2014).** "A variational eigenvalue solver on a photonic quantum processor." *Nat. Commun.* 5, 4213. DOI: `10.1038/ncomms5213`.
   — The founding VQE paper: the hybrid loop (quantum expectation value + classical optimizer) that SOLANGE's Phase 3A proxy is modeled on.

8. **McClean, J. R., Romero, J., Babbush, R., & Aspuru-Guzik, A. (2016).** "The theory of variational hybrid quantum-classical algorithms." *New J. Phys.* 18(2), 023023. DOI: `10.1088/1367-2630/18/2/023023`.
   — Formal theoretical framework for VQE (variational ansätze, UCC connection, error suppression). Why the hybrid scheme is robust on noisy hardware.

9. **Preskill, J. (2018).** "Quantum Computing in the NISQ era and beyond." *Quantum* 2, 79. DOI: `10.22331/q-2018-08-06-79`.
   — Defines the NISQ paradigm motivating hybrid variational methods as the near-term route before fault tolerance.

10. **Kandala, A., Mezzacapo, A., Temme, K., … Gambetta, J. M. (2017).** "Hardware-efficient variational quantum eigensolver for small molecules and quantum magnets." *Nature* 549(7671), 242–246. DOI: `10.1038/nature23879`.
    — First hardware-efficient VQE on a superconducting processor (up to BeH₂). The experimental template for IBM-class hardware execution.

11. **Campbell, E. (2019).** "Random Compiler for Fast Hamiltonian Simulation." *Phys. Rev. Lett.* 123(7), 070503. DOI: `10.1103/PhysRevLett.123.070503`.
    — Introduces **qDRIFT**, the randomized Trotter compiler whose gate cost is independent of the number of Hamiltonian terms. The algorithmic basis for the **sqDRIFT** step in Phase 3B.

12. **Kiss, O., Grossi, M., & Roggero, A. (2023).** "Importance sampling for stochastic quantum simulations." *Quantum* 7, 977. DOI: `10.22331/q-2023-04-13-977`.
    — A stochastic qDRIFT-style variant unifying qDRIFT with importance sampling. Directly supports the stochastic-sampling refinement underlying "sqDRIFT."
    — **NOTE:** "sqDRIFT" is not a single canonical paper title but a descriptor for stochastic qDRIFT variants; this is the strongest peer-reviewed primary instance.

13. **Tilly, J., Chen, H., Cao, S., … Tennyson, J. (2022).** "The Variational Quantum Eigensolver: A review of methods and best practices." *Phys. Rep.* 986, 1–128. DOI: `10.1016/j.physrep.2022.08.003`.
    — The authoritative VQE review (ansätze, measurement, optimizers, error mitigation). Best citation for justifying design choices in the hybrid pipeline.

---

## 3.3 · The classical ceiling & HPC↔quantum architectural evolution

14. **Helgaker, T., Jørgensen, P., & Olsen, J. (2000).** *Molecular Electronic-Structure Theory.* Wiley. DOI: `10.1002/9781119019572` (chapter-level DOIs resolve, e.g. Ch. 11 `…ch11`).
    — Canonical reference for the formal scaling of CI/CASSCF/coupled-cluster and the combinatorial growth of the FCI determinant space. The authoritative anchor for *why* exact methods hit a wall.

15. **Levine, B. G., Durden, A. S., Esch, M. P., Liang, F., & Shu, Y. (2021).** "CAS without SCF—Why to use CASCI and where to get the orbitals." *J. Chem. Phys.* 154(9), 090902. DOI: `10.1063/5.0042147`.
    — Concise peer-reviewed statement of the exponential/factorial scaling of CASSCF/CASCI (~16–20 electron/orbital wall). Modern articulation of the classical-intractability boundary.

16. **Reiher, M., Wiebe, N., Svore, K. M., Wecker, D., & Troyer, M. (2017).** "Elucidating reaction mechanisms on quantum computers." *PNAS* 114(29), 7555–7560. DOI: `10.1073/pnas.1619152114`.
    — Landmark resource-estimation paper (nitrogenase FeMoco). The prototypical "classical breaks down → quantum takes over" argument the platform instantiates.

17. **Bauer, B., Bravyi, S., Motta, M., & Chan, G. K.-L. (2020).** "Quantum Algorithms for Quantum Chemistry and Quantum Materials Science." *Chem. Rev.* 120(22), 12685–12717. DOI: `10.1021/acs.chemrev.9b00829`.
    — Definitive survey mapping classical-method limits onto quantum resource requirements. One-stop citation to situate the VQE proxy in the quantum-advantage landscape.

18. **Alexeev, Y., Amsler, M., Barroca, M. A., et al. (2024).** "Quantum-centric supercomputing for materials science: A perspective on challenges and future directions." *Future Gener. Comput. Syst.* 160, 666–710. DOI: `10.1016/j.future.2024.04.060`.
    — The flagship "quantum-centric supercomputing" perspective on HPC↔quantum integration. Direct citation for the classical→quantum migration design and Phase 3B HPC-quantum coupling.

19. **Fowler, A. G., Mariantoni, M., Martinis, J. M., & Cleland, A. N. (2012).** "Surface codes: Towards practical large-scale quantum computation." *Phys. Rev. A* 86, 032324. DOI: `10.1103/PhysRevA.86.032324`.
    — Standard source on surface-code error-correction thresholds and qubit overhead toward logical qubits. Grounds the roadmap from Heron to fault tolerance.
    — **Optional newer complement (FLAG — DOI unverified this session):** Google Quantum AI (2024), "Quantum error correction below the surface code threshold," *Nature*, DOI `10.1038/s41586-024-08449-y`. Confirm before citing.

---

## 3.4 · Enterprise architecture as a design lens

20. **The Open Group (2022).** *TOGAF Standard, 10th Edition* — Architecture Development Method (ADM). The Open Group Series. ISBN 978-94-018-0859-0 (Core Concepts); ADM module ISBN 978-94-018-0862-0. Catalog: publications.opengroup.org/i220.
    — Provides the ADM: a repeatable, phase-driven method letting the classical-proxy and quantum-hardware phases be expressed as governed, traceable architecture transitions.
    — **NOTE:** 10th Edition is modular — cite the specific volume (Core Concepts vs. ADM) matching the claim.

21. **The Open Group (2022).** *ArchiMate 3.2 Specification.* Document No. C226. ISBN 978-94-018-0955-9. Standard: pubs.opengroup.org/architecture/archimate32-doc/.
    — Formal modeling language (Business/Application/Technology layers + motivation/migration) for visualizing how compute services, provenance, and compliance relate across layers.

22. **Lankhorst, M. (2017).** *Enterprise Architecture at Work: Modelling, Communication and Analysis,* 4th ed. Springer (The Enterprise Engineering Series). DOI: `10.1007/978-3-662-53933-0`.
    — Standard scholarly reference connecting ArchiMate modeling to architecture communication and impact analysis.
    — **NOTE:** 3rd ed. (2013) has a distinct DOI `10.1007/978-3-642-29651-2`; cite the edition used.

23. **Zachman, J. A. (1987).** "A Framework for Information Systems Architecture." *IBM Systems Journal* 26(3), 276–292. DOI: `10.1147/sj.263.0276`.
    — The foundational paper establishing enterprise architecture as a discipline of structured perspectives. Historical legitimacy for treating a scientific platform as a multi-stakeholder architecture.

24. **Pierantoni, G., Kiss, T., Bolotov, A., et al. (2023).** "Toward a reference architecture based science gateway framework with embedded e-learning support." *Concurrency Computat. Pract. Exper.* 35(18), e6872. DOI: `10.1002/cpe.6872`.
    — Peer-reviewed application of reference-architecture thinking to science gateways over clouds/grids/HPC. The closest analog for formally architecting a scientific platform spanning classical and quantum back-ends.

---

## 3.5 · Provenance, reproducibility & compliance architectures

25. **Wilkinson, M. D., Dumontier, M., Aalbersberg, I. J., et al. (2016).** "The FAIR Guiding Principles for scientific data management and stewardship." *Scientific Data* 3, 160018. DOI: `10.1038/sdata.2016.18`.
    — Establishes Findable/Accessible/Interoperable/Reusable criteria a provenance-complete pipeline must satisfy to be machine-actionable and auditable.

26. **Lebo, T., Sahoo, S., McGuinness, D. (eds.) / W3C (2013).** *PROV-O: The PROV Ontology.* W3C Recommendation, 30 April 2013. Stable URI: w3.org/TR/2013/REC-prov-o-20130430/.
    — Standard data model (entities/activities/agents) for serializing the P1–P9 chain in an interoperable, queryable form rather than ad hoc logs.

27. **Sandve, G. K., Nekrutenko, A., Taylor, J., & Hovig, E. (2013).** "Ten Simple Rules for Reproducible Computational Research." *PLOS Comput. Biol.* 9(10), e1003285. DOI: `10.1371/journal.pcbi.1003285`.
    — Actionable rules (track provenance, record exact versions, archive raw inputs) mapping directly to making a VQE/PySCF pipeline reproducible end-to-end.

28. **Stodden, V., McNutt, M., Bailey, D. H., et al. (2016).** "Enhancing reproducibility for computational methods." *Science* 354(6317), 1240–1241. DOI: `10.1126/science.aah6168`.
    — "Reproducibility Enhancement Principles" for disclosing workflows, code, and data. Policy-level complement to Sandve's operational rules.

29. **U.S. FDA (1997).** "Electronic Records; Electronic Signatures" (Final Rule), codified at **21 CFR Part 11.** *Federal Register* 62(54), 13430. Primary: govinfo.gov/content/pkg/FR-1997-03-20/pdf/97-6833.pdf; current text: ecfr.gov/current/title-21/.../part-11.
    — Authoritative legal basis for the §11.10(e) audit-trail, access-control, and signature requirements the pipeline satisfies. Cite alongside FDA's 2003 "Part 11 — Scope and Application" guidance.

30. **Greene, J. E. (2008).** "Implementation of a low-cost Interim 21CFR11 compliance solution for laboratory environments." *JALA: J. Assoc. Lab. Autom.* PMID 18924616; PMC2548384.
    — Worked example of operationalizing Part 11 controls (audit trail, backup/restore, access security) in an automated lab.
    — **⚠ FLAGGED — PARTIALLY UNVERIFIED:** title/PMID/PMCID confirmed, but journal name/volume/issue/pages/DOI not independently confirmed. Retrieve from the PMC record before citing; treat as secondary to #29.

---

## Synthesis & research gap (3.6 — for prose)

The thirty sources establish that **each circle of the focus Venn is mature in isolation**:
the physical theory of quantum simulation (3.1), the hybrid algorithmic layer (3.2), the
classical ceiling and emerging HPC↔quantum integration (3.3), enterprise-architecture method
(3.4), and provenance/compliance standards (3.5). **What is absent in the literature** is a
single system that conjugates all five: a hybrid classical→quantum pipeline (3.1–3.3) that is
(i) framed by a formal enterprise-architecture method — TOGAF/ADM (3.4); (ii) governed by a
complete, interoperable provenance chain under regulatory-grade controls — PROV-O + 21 CFR
Part 11 + FAIR (3.5); and (iii) aimed specifically at the *non-druggable* tumour-suppressor
class. Notably, the EA literature (3.4) has been applied to science gateways and HPC
(Pierantoni 2023) but **never to a classical→quantum scientific pipeline**, and the
quantum-chemistry literature (3.1–3.3) is essentially silent on enterprise-architecture
governance and Part-11 provenance. SOLANGE™ occupies that intersection.

---

### Caveats & gaps to disclose
- **ALCOA+ data-integrity:** no peer-reviewed primary paper coins it; it originates in FDA/MHRA/EMA/WHO *guidance*. Cite the guidance documents (FDA 2018 "Data Integrity and Compliance With Drug CGMP"; MHRA 2018 "GxP Data Integrity Guidance"; WHO TRS 996 Annex 5, 2016) — not a journal DOI.
- **Live DOI resolution** was blocked in-environment; all identifiers cross-corroborated across ≥2 indexes. One-click verify before submission.
- Source **#30 (Greene)** and the optional **Google 2024** error-correction paper are the only entries with incomplete verification — both flagged inline.
