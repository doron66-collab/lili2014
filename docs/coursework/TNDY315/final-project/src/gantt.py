import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import date, timedelta
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

plt.rcParams["font.family"] = "DejaVu Sans"

C_DONE = "#8a96a3"
C_ACTIVE = "#b8441f"
C_PLAN = "#34475e"
C_MGMT = "#2e7d6b"
C_FUND = "#c28a1b"
C_DIP = "#6b4fa0"

d = date
rows = [
    # label, owner, start, end, color, style, bar text
    ("PHASES", None, None, None, None, "header", ""),
    ("P1  Foundation and literature review", "Doron", d(2026,3,1), d(2026,3,20), C_DONE, "solid", "Complete"),
    ("P2  Platform prototype deployed", "Doron", d(2026,3,1), d(2026,4,30), C_DONE, "solid", "Complete"),
    ("P3A Architecture validation on simulator", "Doron", d(2026,5,1), d(2026,8,31), C_ACTIVE, "solid", "20e/40q run complete Aug 4"),
    ("      Early test jobs on real IBM Heron hardware (4-8 qubits)", "Doron", d(2026,7,15), d(2026,8,11), C_DONE, "solid", "Pipeline validated"),
    ("      P3A intermediate 24e local-site run", "Doron", d(2026,9,1), d(2026,11,30), C_ACTIVE, "dashed", "Pending"),
    ("P3B Scaling run on IBM Heron r3 (48e)", "Doron", d(2026,9,1), d(2026,12,31), C_ACTIVE, "dashed", "Contingent on IBM access"),
    ("P4  Triangulated DSR evaluation", "Doron", d(2027,1,1), d(2027,4,30), C_PLAN, "solid", "Jan-Apr 2027"),
    ("      4A Technical harness  |  4B Benchmark", "Doron", d(2027,1,4), d(2027,2,28), C_PLAN, "light", "4A 4 wk, 4B 3 wk"),
    ("      4C Expert-panel recruitment", "Doron + Amit", d(2027,1,4), d(2027,1,25), C_DIP, "light", "3 wk"),
    ("      4D Panel sessions  |  4E Walkthroughs", "Doron", d(2027,2,1), d(2027,3,31), C_PLAN, "light", "4D 4 wk, 4E 2 wk"),
    ("      4F Triangulation and evaluation report", "Doron", d(2027,4,1), d(2027,4,30), C_PLAN, "light", "3 wk"),
    ("P5  Dissertation writing (ch. 1-8)", "Doron", d(2027,5,1), d(2027,9,30), C_PLAN, "solid", "May-Sep 2027"),
    ("P6  Pre-defense committee revisions", "Doron", d(2027,10,1), d(2027,10,31), C_PLAN, "solid", "Oct"),
    ("P7  Defense", "Doron", d(2027,11,1), d(2027,12,31), C_PLAN, "solid", "Nov-Dec"),
    ("Post-defense full-scale Heron r3 runs", "Doron", d(2028,1,1), d(2028,2,28), C_PLAN, "dashed", "2028+"),
    ("MANAGEMENT TRACKS", None, None, None, None, "header", ""),
    ("Quarterly re-baseline against measured capability", "All", d(2026,10,1), d(2027,12,31), C_MGMT, "ticks", ""),
    ("Funding: grants, IBM collaboration, investor briefings", "Manuel", d(2026,10,1), d(2027,6,30), C_FUND, "solid", "BSF, IBM Research Israel, microgrant"),
    ("International partners: hospitals, health diplomacy", "Amit", d(2026,11,1), d(2027,3,31), C_DIP, "solid", "Partner MOUs, expert sourcing"),
    ("Results dissemination to partner health systems", "Amit + Manuel", d(2027,9,1), d(2028,2,28), C_DIP, "dashed", "After OSF-frozen results"),
]

milestones = [
    ("Proposal accepted", d(2026,5,1)),
    ("P3B access decision", d(2026,9,1)),
    ("OSF pre-registration freeze", d(2027,1,1)),
    ("Defense", d(2027,12,31)),
]
today = d(2026,9,25)

fig, ax = plt.subplots(figsize=(13.5, 8.2), dpi=200)
n = len(rows)
y = list(range(n))[::-1]
xmin, xmax = d(2026,3,1), d(2028,3,1)

for yi, (label, owner, s, e, col, style, txt) in zip(y, rows):
    if style == "header":
        ax.axhspan(yi-0.5, yi+0.5, color="#eef1f4", zorder=0)
        ax.text(xmin - timedelta(days=5), yi, label, ha="right", va="center", fontsize=9.5, fontweight="bold", color="#1a2028")
        continue
    bold = not label.startswith("      ")
    ax.text(xmin - timedelta(days=5), yi, label.strip(), ha="right", va="center",
            fontsize=8.6 if bold else 8.0, color="#1a2028" if bold else "#4a5563",
            fontweight="bold" if bold and label[:2] in ("P1","P2","P3","P4","P5","P6","P7") else "normal")
    ax.text(xmax + timedelta(days=5), yi, owner, ha="left", va="center", fontsize=8.2, color="#1a2028")
    width = (e - s).days
    if style == "solid":
        ax.barh(yi, width, left=s, height=0.56, color=col, zorder=2)
        tc = "white"
    elif style == "light":
        ax.barh(yi, width, left=s, height=0.42, color=col, alpha=0.55, zorder=2)
        tc = "#1a2028"
    elif style == "dashed":
        ax.barh(yi, width, left=s, height=0.56, color="white", edgecolor=col, linestyle="--", linewidth=1.4, zorder=2)
        tc = col
    elif style == "ticks":
        ax.plot([s, e], [yi, yi], color=col, linewidth=1.2, zorder=1)
        q = s
        while q <= e:
            ax.plot(q, yi, marker="s", color=col, markersize=6, zorder=3)
            m = q.month + 3
            q = d(q.year + (m-1)//12, (m-1)%12 + 1, 1)
        tc = None
    if txt and tc:
        mid = s + (e - s)/2
        if width > 55:
            ax.text(mid, yi, txt, ha="center", va="center", fontsize=7.2, color=tc, zorder=4)
        else:
            ax.text(e + timedelta(days=4), yi, txt, ha="left", va="center", fontsize=7.2, color="#1a2028", zorder=4)

top = n - 0.5
for i, (name, md) in enumerate(milestones):
    ax.plot(md, top + 0.35, marker="D", color="#5a3a2a", markersize=7, clip_on=False, zorder=5)
    ax.axvline(md, color="#5a3a2a", linewidth=0.6, linestyle=":", zorder=1)
    off, ha = [(0.95,"center"),(1.75,"right"),(1.2,"center"),(0.95,"center")][i]
    dx = {"right": -timedelta(days=4), "left": timedelta(days=4), "center": timedelta(0)}[ha]
    ax.text(md + dx, top + off, name, ha=ha, va="bottom", fontsize=7.4, color="#5a3a2a", zorder=6, bbox=dict(facecolor="white", edgecolor="none", pad=0.6))
ax.axvline(today, color=C_ACTIVE, linewidth=1.3, zorder=4)
ax.text(today, -1.05, "Today, Sep 25, 2026", ha="center", va="top", fontsize=7.6, color=C_ACTIVE, fontweight="bold")

ax.set_xlim(xmin, xmax)
ax.set_ylim(-1.6, n + 3.1)
ax.set_yticks([])
ax.xaxis.set_major_locator(mdates.MonthLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%y"))
ax.tick_params(axis="x", labelsize=7)
for sp in ["top", "right", "left"]:
    ax.spines[sp].set_visible(False)
ax.grid(axis="x", color="#e3e6ea", linewidth=0.5, zorder=0)
ax.text(xmax + timedelta(days=5), n - 0.5 + 0.35, "Owner", ha="left", va="center", fontsize=8.4, fontweight="bold")

legend = [
    Patch(color=C_DONE, label="Complete"),
    Patch(color=C_ACTIVE, label="Active phase"),
    Patch(facecolor="white", edgecolor=C_ACTIVE, linestyle="--", label="Pending or contingent"),
    Patch(color=C_PLAN, label="Planned phase"),
    Patch(color=C_FUND, label="Funding track"),
    Patch(color=C_DIP, label="International partnership track"),
    Line2D([0], [0], marker="s", color=C_MGMT, label="Quarterly re-baseline", markersize=6),
    Line2D([0], [0], marker="D", color="#5a3a2a", linestyle="none", label="Milestone", markersize=6),
]
ax.legend(handles=legend, loc="lower center", bbox_to_anchor=(0.45, -0.17), ncol=4, fontsize=7.2, frameon=False)
plt.subplots_adjust(left=0.27, right=0.91, top=0.9, bottom=0.17)
plt.savefig("gantt.png", dpi=200)
print("saved")
