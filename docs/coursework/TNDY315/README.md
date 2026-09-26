# TNDY 315: Principles of Project Management in a Complex World

Coursework for TNDY 315 at Claremont Graduate University (Dr. Matthew Muga), Fall 2026.
SOLANGE is the applied case throughout.

## final-project/ (due September 27, 2026)

Team: Doron Cohen, Manuel Blanco, Amit Mekel.

| File | What it is |
|---|---|
| `TNDY315_SOLANGE_Final_Paper.docx` / `.pdf` | Group paper, APA 7. Eight body pages: text, landscape timeline (Figure 1), artifacts table (Table 1), PM model. References and an appendix of supporting screenshots follow. |
| `TNDY315_SOLANGE_Final_Presentation.pptx` | 18-slide deck. Speaker notes on every slide name the presenter. |
| `SOLANGE_Speaker_Script.docx` | The speaker notes as one document, ordered by slide. Doron: 1-5 and 16-17. Manuel: 6-10. Amit: 11-15. Slide 18 (references) is not presented. |
| `src/` | Generators. `paper_content.js` holds the paper text; `paper.js` builds the .docx; `deck.js` and `deck_notes.js` build the deck; `script.js` builds the speaker script; `gantt.py` draws the timeline. |
| `assets/` | Timeline image, platform screenshots, poster and manuscript thumbnails used by the paper and the deck. |

Rebuild (from `src/`, with the assets copied next to the scripts): `python3 gantt.py`, `node paper.js`, `node deck.js`, `node script.js`.

## module5-essay/

Module 5 individual essay: leveraging a PMBOK Chapter 4 model (Cynefin, with the cross-cultural communication model as a complement). `build.js` and `content.json` generate the .docx.

## Related

- `arxiv/main.tex` and `arxiv/main.pdf`: manuscript prepared for arXiv submission, September 2026 working draft (not yet submitted). Cited in the final paper as Cohen (2026b).
- `public/dissertation_revised.html`: the dissertation proposal, cited as Cohen (2026a).
