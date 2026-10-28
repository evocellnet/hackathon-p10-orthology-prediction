// Hackathon slide deck for the Ancestra project.
// Build: NODE_PATH=<dir with pptxgenjs> node technical_report/slides/build_deck.js  (run from the repository root)
// Numbers are copied from REPORT.md; pictures come from make_assets.py.
const pptxgen = require("pptxgenjs");
const path = require("path");
const A = (f) => path.join(__dirname, "assets", f);
const OUT = path.join(__dirname, "ancestra_hackathon_slides.pptx");

const THEME = {
  name: "Ancestra",
  headFontFace: "Cambria",
  bodyFontFace: "Calibri",
  colors: {
    dk1: "22303F", lt1: "FFFFFF", dk2: "1F3A5F", lt2: "E9EEF5",
    accent1: "3E6FA3", accent2: "D9822B", accent3: "4E8A5B", accent4: "B5473A",
    accent5: "8E6BB0", accent6: "8A929C", hlink: "3E6FA3", folHlink: "8E6BB0",
  },
};
const HEX = THEME.colors;

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE"; // 13.33 x 7.5
pres.theme = { headFontFace: THEME.headFontFace, bodyFontFace: THEME.bodyFontFace };
pres.title = "Ancestra: an AI-led orthology project";
pres.author = "PI agent (Claude Opus 5.5)";
const C = pres.SchemeColor;
const NAVY = C.text2, INK = C.text1, LIGHT = C.background2, WHITE = C.background1;
const BLUE = C.accent1, ORANGE = C.accent2, GREEN = C.accent3, RED = C.accent4, GREY = C.accent6;

// ---------------------------------------------------------------- layouts
pres.defineSlideMaster({
  title: "TITLE_DARK",
  background: { color: HEX.dk2 },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.8, y: 2.0, w: 11.7, h: 1.6, fontSize: 54, bold: true,
      color: WHITE, fontFace: "Cambria", valign: "bottom", align: "left" }, text: "" } },
    { placeholder: { options: { name: "body", type: "body", x: 0.8, y: 3.75, w: 11.7, h: 1.4, fontSize: 24,
      color: LIGHT, valign: "top" }, text: "" } },
  ],
});
pres.defineSlideMaster({
  title: "SECTION",
  background: { color: HEX.dk2 },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.8, y: 2.9, w: 11.7, h: 1.2, fontSize: 44, bold: true,
      color: WHITE, fontFace: "Cambria", align: "left" }, text: "" } },
    { placeholder: { options: { name: "body", type: "body", x: 0.8, y: 4.15, w: 11.7, h: 0.8, fontSize: 20,
      color: LIGHT }, text: "" } },
  ],
  slideNumber: { x: 12.4, y: 7.0, w: 0.6, h: 0.3, fontSize: 11, color: LIGHT },
});
pres.defineSlideMaster({
  title: "CONTENT",
  background: { color: HEX.lt1 },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.6, y: 0.3, w: 12.1, h: 0.85, fontSize: 32, bold: true,
      color: NAVY, fontFace: "Cambria", valign: "middle", align: "left" }, text: "" } },
    { text: { text: "Ancestra · AI hackathon · agent-generated, not yet verified by a human",
      options: { x: 0.6, y: 7.0, w: 8, h: 0.3, fontSize: 10, color: GREY } } },
  ],
  slideNumber: { x: 12.4, y: 7.0, w: 0.6, h: 0.3, fontSize: 11, color: GREY },
});

let section = "";
function slide(master, sec) {
  if (sec && sec !== section) { pres.addSection({ title: sec }); section = sec; }
  return pres.addSlide({ masterName: master, sectionTitle: section });
}
function title(s, t) { s.addText(t, { placeholder: "title" }); }

// ---------------------------------------------------------------- helpers
function txt(s, text, o) {
  s.addText(text, Object.assign({ isTextBox: true, fontSize: 16, color: INK, margin: 0, valign: "top" }, o));
}
function card(s, x, y, w, h, fill, o = {}) {
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h, rectRadius: 0.12, fill: { color: fill, transparency: o.t || 0 },
    line: { color: o.line || fill, width: o.lw || 0 }, objectName: o.name });
}
function chip(s, x, y, w, h, label, fill, color = WHITE, fs = 14) {
  card(s, x, y, w, h, fill);
  s.addText(label, { x, y, w, h, isTextBox: true, fontSize: fs, color, bold: true, align: "center", valign: "middle", margin: 0.04 });
}
function dot(s, x, y, d, fill, label, fs = 16) {
  s.addShape(pres.shapes.OVAL, { x, y, w: d, h: d, fill: { color: fill }, line: { color: fill, width: 0 } });
  if (label) s.addText(label, { x, y, w: d, h: d, isTextBox: true, fontSize: fs, bold: true, color: WHITE, align: "center", valign: "middle", margin: 0 });
}
function arrowLine(s, x1, y1, x2, y2, color = GREY, dash) {
  const o = { x: Math.min(x1, x2), y: Math.min(y1, y2), w: Math.abs(x2 - x1) || 0.001, h: Math.abs(y2 - y1) || 0.001,
    line: { color, width: 2, endArrowType: "triangle", dashType: dash || "solid" } };
  if (x2 < x1) o.flipH = true;
  if (y2 < y1) o.flipV = true;
  s.addShape(pres.shapes.LINE, o);
}
function stat(s, x, y, w, big, small, color = NAVY) {
  txt(s, big, { x, y, w, h: 0.95, fontSize: 44, bold: true, color, fontFace: "Cambria", valign: "bottom" });
  txt(s, small, { x, y: y + 1.0, w, h: 0.8, fontSize: 14, color: INK });
}
function img(s, file, x, y, w, h, alt) { s.addImage({ path: A(file), x, y, w, h, altText: alt }); }

// ================================================================ 1. title
{
  const s = slide("TITLE_DARK", "Introduction");
  s.addText("Ancestra", { placeholder: "title" });
  s.addText("An AI agent leads an orthology research project in 48 hours", { placeholder: "body" });
  txt(s, "AI hackathon · Beltrao lab, ETH Zurich · [date to be added]", { x: 0.8, y: 5.4, w: 11, h: 0.4, fontSize: 18, color: WHITE });
  txt(s, "Slides generated by the PI agent (Claude Opus 5.5). Results are agent-generated and not yet verified by a human.",
    { x: 0.8, y: 6.5, w: 11.5, h: 0.4, fontSize: 13, color: LIGHT, italic: true });
  s.addNotes("Title. This deck documents the orthology project of the AI hackathon. Everything after the title page, including these slides, was produced by the agent.");
}

// ================================================================ 2. hackathon
{
  const s = slide("CONTENT");
  title(s, "The hackathon: can AI agents run bioinformatics projects?");
  img(s, "hackathon_flow.png", 0.6, 1.35, 7.7, 4.44, "Hackathon project flow: phase 1 predict X, phase 2 models and datasets, phase 3 discovery");
  txt(s, "Project flow shared by the organisers. This project: orthology. Its output feeds phase 3 (protein-universe similarity, complex evolution).",
    { x: 0.6, y: 5.9, w: 7.7, h: 0.8, fontSize: 13, color: GREY, italic: true });
  const rows = [
    ["Organisers", "Beltrao (ETH Zurich), Ryan (UCD), Gonçalves (IST), Dessimoz (UNIL) labs", BLUE],
    ["Projects", "12 projects, each led by Claude Opus 5.5 with Opus/Sonnet subagents", GREEN],
    ["Rule", "AI builds, benchmarks and documents; humans only verify", ORANGE],
    ["Dates · people", "[to be added]", GREY],
    ["Code", "github.com/evocellnet/ai-hackathon", NAVY],
  ];
  rows.forEach(([h, b, c], i) => {
    const y = 1.4 + i * 1.08;
    dot(s, 8.7, y + 0.12, 0.42, c);
    txt(s, h, { x: 9.3, y, w: 3.5, h: 0.35, fontSize: 15, bold: true, color: NAVY });
    txt(s, b, { x: 9.3, y: y + 0.36, w: 3.6, h: 0.65, fontSize: 13 });
  });
  s.addNotes("Information from the hackathon repository README and the project-flow figure. Dates, participant list and the final hackathon link are left blank on purpose: the agent has no access to them.");
}

// ================================================================ 3. set-up
{
  const s = slide("CONTENT");
  title(s, "One PI agent, many subagents, one repository");
  img(s, "agents.png", 0.5, 1.4, 8.1, 3.96, "PI agent delegating to Opus and Sonnet subagents");
  stat(s, 9.0, 1.3, 3.8, "48 h", "budget, checked against START.txt at every milestone");
  stat(s, 9.0, 3.1, 3.8, "24 → 48", "CPUs; an A100 GPU was added after 19.5 h");
  stat(s, 9.0, 4.9, 3.8, "7", "human instructions after the first prompt");
  s.addNotes("The PI agent never wrote method code itself. It split the work, launched subagents in the background (often several at once), reviewed their reports and committed decisions to git. Opus subagents: literature, design, algorithms, verification. Sonnet subagents: environment, competitor and scaling runs, bookkeeping, reference checks. Human instructions: Foldcomp databases, include an OMA resource in every comparison, pause and move to a GPU node, resume, soft deadline plus README header, write the report, second round of data layers.");
}

// ================================================================ 4. prompt
{
  const s = slide("CONTENT");
  title(s, "The prompt: half research protocol, half orthology task");
  card(s, 0.6, 1.3, 12.1, 0.85, LIGHT);
  txt(s, "“Obtain a new orthology prediction method that is conceptually different from the state of the art, with comparable if not better accuracy and/or compute time.”",
    { x: 0.8, y: 1.38, w: 11.7, h: 0.7, fontSize: 16, italic: true, color: NAVY, valign: "middle" });
  const cols = [
    ["General: any AI-led project", BLUE, ["PI role: delegate, review, decide", "Sub-tasks, each with its own benchmark",
      "Literature + design critiqued by agents", "Reproducible: seeds, DATA.md, no metric changes",
      "Results checked by a different agent", "Research tree + lessons file", "48 h budget with milestones", "Journal-style REPORT"]],
    ["Specific: orthology", ORANGE, ["Conceptually different method", "Dozens → thousands of proteomes",
      "Hierarchical groups in OrthoXML", "Alternatives to one species tree", "Extra data: structure, synteny, omics, isoforms, DNA",
      "Benchmarks: QfO, OrthoBench, SwissTree", "Separate benchmark repository", "Consider GPU speed-ups"]],
  ];
  cols.forEach(([h, c, items], k) => {
    const x0 = 0.6 + k * 6.15;
    txt(s, h, { x: x0, y: 2.4, w: 5.9, h: 0.45, fontSize: 20, bold: true, color: c, fontFace: "Cambria" });
    items.forEach((t, i) => {
      const x = x0 + (i % 2) * 3.0, y = 3.0 + Math.floor(i / 2) * 0.95;
      card(s, x, y, 2.85, 0.8, c, { t: 82 });
      txt(s, t, { x: x + 0.12, y, w: 2.65, h: 0.8, fontSize: 13, valign: "middle" });
    });
  });
  s.addNotes("The full prompt is in PROMPT.md (verbatim). Left: rules that would apply to any agent-led project and could be reused by the other hackathon projects. Right: what is specific to orthology.");
}

// ================================================================ 5. timeline
{
  const s = slide("CONTENT");
  title(s, "48 hours at a glance");
  const x0 = 0.8, x1 = 12.5, y = 3.6, sc = (x1 - x0) / 48;
  s.addShape(pres.shapes.LINE, { x: x0, y, w: x1 - x0, h: 0, line: { color: NAVY, width: 3 } });
  s.addShape(pres.shapes.RECTANGLE, { x: x0 + 15.4 * sc, y: y - 0.06, w: (48 - 15.4) * sc, h: 0.12, fill: { color: ORANGE, transparency: 40 }, line: { color: ORANGE, width: 0 } });
  const ev = [
    [0, "Start", "prompt + START.txt", 1],
    [1.1, "Design v2", "8 agent reviews", -1],
    [5.2, "Core works", "OrthoBench dev", 1],
    [8.0, "QfO peek 1", "recall too low", -1],
    [15.4, "Freeze v1", "two presets", 1],
    [17.5, "Test scored", "OrthoBench, once", -1],
    [19.3, "QfO final", "then GPU node", 1],
    [22.5, "Report", "writing starts", -1],
    [27.3, "Layers round 2", "5 agents in parallel", 1],
    [39.6, "1,000 proteomes", "scaling done", -1],
    [48, "Deadline", "soft", 1],
  ];
  ev.forEach(([t, h, b, side]) => {
    const x = x0 + t * sc;
    dot(s, x - 0.11, y - 0.11, 0.22, t <= 15.4 ? NAVY : ORANGE);
    const ty = side > 0 ? y - 1.55 : y + 0.35;
    s.addShape(pres.shapes.LINE, { x, y: side > 0 ? y - 0.45 : y + 0.12, w: 0, h: 0.33, line: { color: GREY, width: 1 } });
    txt(s, h, { x: x - 0.8, y: side > 0 ? ty + 0.35 : ty + 0.12, w: 1.6, h: 0.35, fontSize: 14, bold: true, color: NAVY, align: "center" });
    txt(s, b, { x: x - 0.8, y: side > 0 ? ty + 0.68 : ty + 0.45, w: 1.6, h: 0.45, fontSize: 12, color: GREY, align: "center" });
  });
  card(s, 0.8, 5.6, 3.6, 0.55, NAVY, { t: 85 });
  txt(s, "Build and tune on development data", { x: 0.9, y: 5.6, w: 3.5, h: 0.55, fontSize: 13, color: NAVY, valign: "middle" });
  card(s, 4.7, 5.6, 4.6, 0.55, ORANGE, { t: 75 });
  txt(s, "After the freeze: evaluation, GPU, layers, scaling", { x: 4.8, y: 5.6, w: 4.5, h: 0.55, fontSize: 13, color: INK, valign: "middle" });
  s.addNotes("Times are hours on the budget clock (the one-hour pause to move nodes is excluded). Freeze at T+15.4 h. Everything after the freeze is evaluation or post-freeze exploration on development data only. An API session limit at T+11.7 h stopped four agents; their background jobs survived.");
}

// ================================================================ 6. section
{
  const s = slide("SECTION", "How the method was chosen");
  s.addText("How the method was chosen", { placeholder: "title" });
  s.addText("Literature, agent critique, cost estimates and pre-registered rules", { placeholder: "body" });
}

// ================================================================ 7. decision path
{
  const s = slide("CONTENT", "How the method was chosen");
  title(s, "From idea to frozen method");
  const top = [
    ["Literature", "SYNERGY, Hieranoid, OMA/FastOMA, OrthoFinder", LIGHT, INK],
    ["Design v0", "merge whole ancestral proteomes up the tree", RED, WHITE],
    ["Design v2", "families first, merge inside each family", GREEN, WHITE],
    ["Dev tuning", "OrthoBench, 34 reference groups", BLUE, WHITE],
  ];
  const bot = [
    ["QfO peek 1", "recall too low: 78% of misses are merge failures", ORANGE, WHITE],
    ["Fixes", "small FastTree trees + deep-node rule", BLUE, WHITE],
    ["Freeze v1", "accurate + fast presets", NAVY, WHITE],
    ["Test, once", "OrthoBench + QfO", GREEN, WHITE],
  ];
  const w = 2.6, gap = 0.55, xs = (i) => 0.7 + i * (w + gap);
  top.forEach(([h, b, f, c], i) => {
    card(s, xs(i), 1.6, w, 1.5, f);
    txt(s, h, { x: xs(i) + 0.15, y: 1.7, w: w - 0.3, h: 0.4, fontSize: 17, bold: true, color: c });
    txt(s, b, { x: xs(i) + 0.15, y: 2.12, w: w - 0.3, h: 0.9, fontSize: 13, color: c });
    if (i < 3) arrowLine(s, xs(i) + w + 0.05, 2.35, xs(i + 1) - 0.05, 2.35);
  });
  bot.forEach(([h, b, f, c], i) => {
    card(s, xs(i), 4.3, w, 1.5, f);
    txt(s, h, { x: xs(i) + 0.15, y: 4.4, w: w - 0.3, h: 0.4, fontSize: 17, bold: true, color: c });
    txt(s, b, { x: xs(i) + 0.15, y: 4.82, w: w - 0.3, h: 0.9, fontSize: 13, color: c });
    if (i < 3) arrowLine(s, xs(i) + w + 0.05, 5.05, xs(i + 1) - 0.05, 5.05);
  });
  arrowLine(s, xs(3) + w / 2, 3.15, xs(0) + w / 2, 4.25, GREY, "dash");
  txt(s, "✗ only ~2.5× cheaper than all-vs-all (reviewer's cost estimate)", { x: xs(1) - 0.1, y: 3.2, w: 3.4, h: 0.6, fontSize: 12, color: RED, italic: true });
  txt(s, "Dropped after review: “novel” parts that SYNERGY, InParanoid/Hieranoid and OrthoFinder/SHOOT already had.",
    { x: 0.7, y: 6.1, w: 11.9, h: 0.5, fontSize: 13, color: GREY, italic: true });
  s.addNotes("v0 was rejected before any code was written: a reviewer agent estimated 2.0e11 pairs vs 4.9e11 for all-vs-all on QfO, because ancestral proteomes stay large and 60% of deep ancestral genes have no partner. v2 puts all proteins into families first and never searches ancestral genes again; their scores are aggregated from member proteins (as in SYNERGY). QfO peek 1 exposed low recall; a dev-only diagnosis found 78% of missed TreeFam-A pairs were merge failures at the pair's last common ancestor.");
}

// ================================================================ 8. rules
{
  const s = slide("CONTENT");
  title(s, "Four rules steered every decision");
  const r = [
    ["1", "Pre-registration", "Hypothesis, arms and keep rule written before results", BLUE],
    ["2", "Dev / test split", "OrthoBench split by homology: 34 dev, 36 test groups", GREEN],
    ["3", "Shuffled controls", "Every data layer vs a shuffled copy (seeds 11, 23, 47)", ORANGE],
    ["4", "Independent checks", "Another agent recomputes each result before it is used", NAVY],
  ];
  r.forEach(([n, h, b, c], i) => {
    const x = 0.6 + (i % 2) * 6.15, y = 1.45 + Math.floor(i / 2) * 2.25;
    card(s, x, y, 5.9, 1.95, LIGHT);
    dot(s, x + 0.3, y + 0.35, 0.8, c, n, 24);
    txt(s, h, { x: x + 1.35, y: y + 0.3, w: 4.4, h: 0.45, fontSize: 20, bold: true, color: NAVY });
    txt(s, b, { x: x + 1.35, y: y + 0.85, w: 4.3, h: 0.9, fontSize: 15 });
  });
  txt(s, "83", { x: 0.6, y: 5.75, w: 1.3, h: 0.8, fontSize: 40, bold: true, color: NAVY, fontFace: "Cambria", valign: "middle" });
  txt(s, "claims in the evidence ledger", { x: 1.95, y: 5.75, w: 1.9, h: 0.8, fontSize: 14, valign: "middle" });
  txt(s, "19", { x: 4.0, y: 5.75, w: 1.0, h: 0.8, fontSize: 40, bold: true, color: NAVY, fontFace: "Cambria", valign: "middle" });
  txt(s, "deviations disclosed", { x: 5.0, y: 5.75, w: 1.9, h: 0.8, fontSize: 14, valign: "middle" });
  txt(s, "Test sets were scored once, after the freeze. Keep rule: lower 95% CI of the gain above 0, beating the shuffled control, no guard metric worse.",
    { x: 7.4, y: 5.75, w: 5.3, h: 0.9, fontSize: 13, color: GREY, italic: true });
  s.addNotes("The non-inferiority margin was -0.03 (3 F points on OrthoBench). Paired bootstrap B=2000, seed 11. Verification found real errors, for example omitted pre-registered comparators and a preset mismatch in the scaling claims.");
}

// ================================================================ 9. section method
{
  const s = slide("SECTION", "The method");
  s.addText("The method: Ancestra", { placeholder: "title" });
  s.addText("Rebuild the genes of every ancestor, one tree level at a time, inside protein families", { placeholder: "body" });
}

// ================================================================ 10. method schema
{
  const s = slide("CONTENT", "The method");
  title(s, "Ancestra in three steps");
  img(s, "method_schema.png", 0.6, 1.5, 12.1, 4.22, "Pipeline: proteomes to families, scores inside families, merge up the species tree, OrthoXML output");
  const items = [["No reference database", BLUE], ["Ancestral genes never searched again", GREEN], ["Every tree level gives HOGs", ORANGE], ["Optional: gene order, 3Di", GREY]];
  items.forEach(([t, c], i) => chip(s, 0.6 + i * 3.05, 6.05, 2.9, 0.6, t, c, WHITE, 13));
  s.addNotes("S0: DIAMOND clustering at 30% identity, all genes vs cluster representatives, connected components, Leiden split of components above 7,500 genes. S2: DIAMOND inside families plus MMseqs2 re-alignment of missed pairs. S3: walk the NCBI taxonomy bottom-up; at each node merge child ancestral genes by mutual best matches, absorb in-paralogs, and resolve ambiguous clusters with a small FastTree (accurate) or NJ (fast) tree labelled by species overlap. Scores of an ancestral gene are aggregated from its member proteins.");
}

// ================================================================ 11. comparison
{
  const s = slide("CONTENT");
  title(s, "Same goal as FastOMA, a different route");
  const Y = "✓", N = "✗";
  const hdr = ["", "Ancestra", "FastOMA", "OrthoFinder", "SonicParanoid2"].map((t, i) => ({ text: t, options: { bold: true, color: i === 1 ? WHITE : NAVY, fill: { color: i === 1 ? HEX.dk2 : HEX.lt2 }, align: "center" } }));
  const rows = [
    ["Reference database needed", N, "OMA (OMAmer)", N, N],
    ["All-vs-all search", N, N, Y, Y + " (reduced)"],
    ["HOGs level by level up the species tree", Y, Y, "from gene trees", N],
    ["Families found first", "de novo", "placed in reference", "orthogroups", N],
    ["Gene trees", "small, ambiguous cases", "per HOG and level", "every orthogroup", N],
    ["Gene order used", "optional", N, N, N],
  ];
  const cell = (v, j) => ({ text: v, options: { align: j === 0 ? "left" : "center", bold: j === 1,
    color: v === Y ? HEX.accent3 : v === N ? HEX.accent4 : HEX.dk1, fill: { color: j === 1 ? "EEF3F9" : "FFFFFF" } } });
  s.addTable([hdr, ...rows.map((r) => r.map(cell))], { x: 0.6, y: 1.4, w: 12.1, colW: [3.9, 2.05, 2.05, 2.05, 2.05],
    fontSize: 15, rowH: 0.62, border: { type: "solid", pt: 0.75, color: "D5DCE6" }, valign: "middle" });
  card(s, 0.6, 5.95, 12.1, 0.75, LIGHT);
  txt(s, "Ideas credited: bottom-up merging with synteny (SYNERGY) · in-paralog absorption (InParanoid, Hieranoid) · HOGs per level (OMA, FastOMA) · representatives (OrthoFinder/SHOOT)",
    { x: 0.8, y: 5.95, w: 11.8, h: 0.75, fontSize: 13, color: NAVY, valign: "middle" });
  s.addNotes("Commonality with FastOMA: both build hierarchical orthologous groups level by level along a species tree, family by family. Differences: FastOMA places proteins into a precomputed OMA reference with OMAmer; Ancestra finds families de novo and uses score rules with small local trees only where candidates are ambiguous. SonicParanoid2 reduces its all-vs-all searches with machine learning but remains pairwise.");
}

// ================================================================ 12. section results
{
  const s = slide("SECTION", "Results");
  s.addText("Results", { placeholder: "title" });
  s.addText("OrthoBench, Quest for Orthologs 2022 and scaling to 1,000 proteomes", { placeholder: "body" });
}

// ================================================================ 13. OrthoBench
{
  const s = slide("CONTENT", "Results");
  title(s, "OrthoBench: highest F, lead not significant");
  const names = ["FastOMA 0.5.1*", "ProteinOrtho", "OMA HOGs*", "OrthoFinder 3.2", "SonicParanoid2", "Ancestra fast", "Ancestra accurate"];
  const f = [33.1, 47.5, 60.4, 69.7, 70.4, 74.8, 76.2];
  const isA = (n) => n.startsWith("Ancestra");
  s.addChart(pres.charts.BAR, [
    { name: "Ancestra", labels: names, values: names.map((n, i) => (isA(n) ? f[i] : 0)) },
    { name: "Others", labels: names, values: names.map((n, i) => (isA(n) ? 0 : f[i])) },
  ], { x: 0.5, y: 1.35, w: 7.6, h: 5.4, barDir: "bar", barGrouping: "stacked", chartColors: [HEX.accent2, HEX.accent1],
    showValue: true, dataLabelPosition: "inEnd", dataLabelFormatCode: "0.0;;;", dataLabelColor: "FFFFFF", dataLabelFontSize: 12,
    dataLabelFontFace: "+mn-lt", showLegend: false, showTitle: true, title: "Test F (36 held-out groups)", titleFontSize: 14,
    titleColor: HEX.dk2, titleFontFace: "+mn-lt", catAxisLabelColor: HEX.dk1, valAxisLabelColor: HEX.accent6, catAxisLabelFontSize: 13,
    catAxisLabelFontFace: "+mn-lt", valAxisLabelFontFace: "+mn-lt", valAxisMaxVal: 100, valAxisMinVal: 0,
    valGridLine: { color: "E2E7EE", size: 0.5 }, catGridLine: { style: "none" }, barGapWidthPct: 45 });
  stat(s, 8.6, 1.4, 4.3, "+6.5 F", "vs OrthoFinder, 95% CI [−6.5, +18.1]: non-inferiority not shown", ORANGE);
  stat(s, 8.6, 3.2, 4.3, "+15.7 F", "vs OMA HOGs [+1.8, +29.2], significant", GREEN);
  stat(s, 8.6, 5.0, 4.3, "3.9 vs 8.2", "CPU-hours, Ancestra vs OrthoFinder", NAVY);
  txt(s, "* OMA contains all 12 benchmark species", { x: 0.6, y: 6.7, w: 6, h: 0.3, fontSize: 10, color: GREY });
  s.addNotes("Official Open OrthoBench scoring, test half only (36 reference orthogroups split by homology component), scored once after the freeze. Paired bootstrap over groups. With 36 groups the standard error of the difference is 5-6 F points. FastOMA was patched for a gene-dropping bug; F did not change (33.5 to 33.1).");
}

// ================================================================ 14. QfO scatter
{
  const s = slide("CONTENT");
  title(s, "QfO 2022: precise, but recall stays below the best methods");
  const pts = [
    ["Ancestra accurate", 0.583, 0.921], ["Ancestra fast", 0.608, 0.916], ["OrthoFinder", 0.737, 0.918],
    ["SonicParanoid2-sens-g", 0.667, 0.936], ["Hieranoid 2", 0.621, 0.930], ["InParanoid 5", 0.627, 0.933],
    ["FastOMA", 0.546, 0.919], ["OMA HOGs", 0.580, 0.923], ["OMA Pairs", 0.565, 0.920],
  ];
  const anc = pts.filter((p) => p[0].startsWith("Ancestra")), oth = pts.filter((p) => !p[0].startsWith("Ancestra"));
  const xs = pts.map((p) => p[1]);
  const series = [{ name: "TPR", values: xs },
    { name: "Ancestra", values: pts.map((p) => (p[0].startsWith("Ancestra") ? p[2] : null)), labels: pts.map((p) => (p[0].startsWith("Ancestra") ? p[0] : "")) },
    { name: "Other methods", values: pts.map((p) => (p[0].startsWith("Ancestra") ? null : p[2])), labels: pts.map((p) => (p[0].startsWith("Ancestra") ? "" : p[0])) }];
  s.addChart(pres.charts.SCATTER, series, { x: 0.5, y: 1.35, w: 8.3, h: 5.45, lineSize: 0, lineDataSymbolSize: 12,
    chartColors: [HEX.accent2, HEX.accent1], showLabel: true, dataLabelFormatScatter: "custom", dataLabelPosition: "r",
    dataLabelFontSize: 11, dataLabelFontFace: "+mn-lt", dataLabelColor: HEX.dk1,
    showTitle: true, title: "TreeFam-A test: recall (TPR) vs precision (PPV)", titleFontSize: 14, titleColor: HEX.dk2, titleFontFace: "+mn-lt",
    valAxisMinVal: 0.91, valAxisMaxVal: 0.94, catAxisMinVal: 0.52, catAxisMaxVal: 0.78, valAxisLabelFormatCode: "0.000",
    catAxisLabelFormatCode: "0.00", showValAxisTitle: true, valAxisTitle: "precision (PPV)", showCatAxisTitle: true, catAxisTitle: "recall (TPR)",
    valAxisTitleFontSize: 12, catAxisTitleFontSize: 12, valAxisTitleColor: HEX.accent6, catAxisTitleColor: HEX.accent6,
    catAxisLabelColor: HEX.accent6, valAxisLabelColor: HEX.accent6, valAxisLabelFontFace: "+mn-lt", catAxisLabelFontFace: "+mn-lt",
    valGridLine: { color: "E2E7EE", size: 0.5 }, catGridLine: { color: "E2E7EE", size: 0.5 }, showLegend: false });
  stat(s, 9.2, 1.4, 3.7, "−0.16", "TreeFam-A recall vs OrthoFinder", RED);
  stat(s, 9.2, 3.2, 3.7, "≈ OMA", "recall on par with FastOMA and OMA", NAVY);
  txt(s, "Misses sit at deep nodes: 34% Bilateria, 21% Eukaryota, 13% Opisthokonta, where an outgroup is missing or distant.",
    { x: 9.2, y: 5.0, w: 3.7, h: 1.6, fontSize: 14, color: INK });
  s.addNotes("Full official QfO 2022 suite run locally (78 species, 988,778 proteins). Local scoring reproduced 60 of 60 published TPR/PPV values. Competitor rows are the published OpenEBench submissions. Ancestra also has the worst G-STD2 Eukaryota species-tree discordance of the main comparators.");
}

// ================================================================ 15. QfO forest
{
  const s = slide("CONTENT");
  title(s, "QfO test: non-inferior only to FastOMA and OMA Pairs");
  img(s, "qfo_forest.png", 0.5, 1.35, 8.4, 5.02, "Forest plot of TreeFam-A recall differences with 95% confidence intervals");
  const v = [["non-inferior", "FastOMA, OMA Pairs", GREEN], ["not shown", "OMA HOGs, Hieranoid 2", ORANGE], ["fails", "InParanoid 5, SonicParanoid2, OrthoFinder", RED]];
  v.forEach(([h, b, c], i) => {
    const y = 1.6 + i * 1.5;
    chip(s, 9.3, y, 1.9, 0.5, h, c, WHITE, 14);
    txt(s, b, { x: 9.3, y: y + 0.6, w: 3.6, h: 0.7, fontSize: 15 });
  });
  txt(s, "Margin −0.03 recall, family bootstrap over 749 test families.", { x: 9.3, y: 6.1, w: 3.6, h: 0.6, fontSize: 12, color: GREY, italic: true });
  s.addNotes("Ancestra accurate minus each method, TreeFam-A test recall, 95% CI. Five primary comparators were pre-registered (OrthoFinder, FastOMA, Hieranoid 2, InParanoid 5, SonicParanoid2-sens-g); the OMA resources were added after the human instruction to include an OMA resource in every comparison.");
}

// ================================================================ 16. scaling
{
  const s = slide("CONTENT");
  title(s, "Scaling: 1,000 proteomes on one server node");
  img(s, "fig3_scaling.png", 0.4, 1.3, 9.4, 4.09, "CPU hours against number of proteins for QfO subsets and 125 to 1,000 UniProt proteomes");
  stat(s, 10.0, 1.3, 3.0, "10.6 h", "wall time, 1,000 proteomes (4.75 M proteins)");
  stat(s, 10.0, 3.1, 3.0, "221 GiB", "peak memory, the likely first limit", RED);
  stat(s, 10.0, 4.9, 3.0, "16 vs 37", "CPU-h at 64 species: Ancestra vs FastOMA", GREEN);
  txt(s, "Cost is superlinear (CPU exponent 1.16 accurate on QfO; 1.35 fast at 0.5–4.75 M proteins). The family graph step grows about quadratically: the next thing to fix.",
    { x: 0.6, y: 5.6, w: 9.0, h: 1.0, fontSize: 14, color: INK });
  s.addNotes("OrthoFinder at 64 species was stopped at the 6 h cap (more than 45.9 CPU-h); SonicParanoid2 failed at 64 species with an internal error. Single runs on shared nodes; the QfO series spans two nodes. No competitor was run on the 1,000-proteome set and accuracy cannot be measured there; all outputs passed OrthoXML validation.");
}

// ================================================================ 17. take-home
{
  const s = slide("CONTENT");
  title(s, "Take-home messages");
  const t = [
    ["76.2", "OrthoBench F (OrthoFinder 69.7); lead not significant", ORANGE],
    ["−0.16", "QfO recall vs OrthoFinder: deep nodes are the weak spot", RED],
    ["1,000", "proteomes in 10.6 h; memory, not time, limits next", GREEN],
    ["0 / 13", "extra data layers improved held-out results", NAVY],
  ];
  t.forEach(([b, sm, c], i) => {
    const x = 0.6 + i * 3.1;
    card(s, x, 1.55, 2.85, 3.3, LIGHT);
    txt(s, b, { x: x + 0.2, y: 1.85, w: 2.5, h: 1.1, fontSize: 48, bold: true, color: c, fontFace: "Cambria", align: "center" });
    txt(s, sm, { x: x + 0.2, y: 3.1, w: 2.45, h: 1.5, fontSize: 15, align: "center" });
  });
  card(s, 0.6, 5.3, 12.1, 1.2, NAVY);
  txt(s, "Building HOGs ancestor by ancestor inside families is viable and cheap, but it does not yet beat the best methods on distant orthologs.",
    { x: 0.9, y: 5.3, w: 11.5, h: 1.2, fontSize: 20, color: WHITE, valign: "middle", fontFace: "Cambria" });
}

// ================================================================ 18. section layers
{
  const s = slide("SECTION", "Extra data layers");
  s.addText("Extra data layers", { placeholder: "title" });
  s.addText("13 layers in two rounds, each against a shuffled control", { placeholder: "body" });
}

// ================================================================ 19. where layers plug in
{
  const s = slide("CONTENT", "Extra data layers");
  title(s, "Where each layer entered the pipeline");
  const stages = [["Input", 0.6], ["S0 families", 3.2], ["S2 scores", 5.8], ["S3 merge", 8.4], ["Quality control", 11.0]];
  stages.forEach(([t, x], i) => {
    chip(s, x, 1.45, 1.75, 0.65, t, NAVY, WHITE, 14);
    if (i < 4) arrowLine(s, x + 1.8, 1.78, stages[i + 1][1] - 0.05, 1.78);
  });
  const L = (t, c) => [t, c];
  const col = {
    0.6: [L("Isoforms", RED)],
    3.2: [L("3D structure rescue", ORANGE), L("pLM neighbours", RED)],
    5.8: [L("GPU alignment", GREEN)],
    8.4: [L("Gene order", ORANGE), L("Species tree", GREEN), L("Domains", RED), L("Expression", RED), L("Networks", RED), L("Abundance", RED), L("Protein LMs", RED), L("Introns", RED)],
    11.0: [L("Genomic DNA", GREEN), L("Synteny repair", ORANGE)],
  };
  Object.entries(col).forEach(([x, arr]) => arr.forEach(([t, c], i) => {
    const xx = Number(x) + (Number(x) === 8.4 && i >= 4 ? 0 : 0);
    const cx = Number(x) === 8.4 ? 8.4 + (i % 2) * 1.32 - 0.45 : xx;
    const cy = Number(x) === 8.4 ? 2.45 + Math.floor(i / 2) * 0.72 : 2.45 + i * 0.72;
    card(s, cx, cy, Number(x) === 8.4 ? 1.25 : 1.75, 0.58, c, { t: 15 });
    txt(s, t, { x: cx, y: cy, w: Number(x) === 8.4 ? 1.25 : 1.75, h: 0.58, fontSize: 12, bold: true, color: WHITE, align: "center", valign: "middle" });
  }));
  const leg = [["kept / useful", GREEN], ["option, no held-out gain", ORANGE], ["rejected", RED]];
  leg.forEach(([t, c], i) => { dot(s, 0.6 + i * 3.6, 6.25, 0.3, c); txt(s, t, { x: 1.0 + i * 3.6, y: 6.22, w: 3.1, h: 0.35, fontSize: 14 }); });
  s.addNotes("Round 1 (before the freeze): gene order, 3D structure, expression, isoforms, genomic DNA. Round 2 (after the freeze, development data only): networks, abundance, domain architecture, protein language models including a model trained on our own synteny labels, intron positions, synteny-driven fission repair, expression revisited, structure rescue at full coverage. The species tree and GPU slides are not data layers but were explicit asks of the prompt.");
}

// ================================================================ 20-31. one slide per layer
function layerSlide(o) {
  const s = slide("CONTENT");
  title(s, o.title);
  img(s, o.img, 0.5, 1.45, o.w || 6.4, o.h || 4.0, o.alt);
  chip(s, 0.6, 5.9, 2.2, 0.48, o.round, GREY, WHITE, 12);
  txt(s, o.what, { x: 2.95, y: 5.85, w: 4.1, h: 0.95, fontSize: 13, color: GREY, italic: true });
  const x = 7.35, w = 5.4;
  txt(s, "✓  Pros", { x, y: 1.45, w, h: 0.4, fontSize: 18, bold: true, color: GREEN });
  s.addText(o.pros.map((t, i) => ({ text: t, options: { bullet: true, breakLine: i < o.pros.length - 1 } })),
    { x, y: 1.9, w, h: 1.65, isTextBox: true, fontSize: 15, color: INK, margin: 0, valign: "top", paraSpaceAfter: 6 });
  txt(s, "✗  Cons", { x, y: 3.6, w, h: 0.4, fontSize: 18, bold: true, color: RED });
  s.addText(o.cons.map((t, i) => ({ text: t, options: { bullet: true, breakLine: i < o.cons.length - 1 } })),
    { x, y: 4.05, w, h: 1.65, isTextBox: true, fontSize: 15, color: INK, margin: 0, valign: "top", paraSpaceAfter: 6 });
  card(s, x, 5.85, w, 0.75, o.vc);
  txt(s, o.verdict, { x: x + 0.2, y: 5.85, w: w - 0.3, h: 0.75, fontSize: 15, bold: true, color: WHITE, valign: "middle" });
  if (o.notes) s.addNotes(o.notes);
}

layerSlide({ title: "Gene order (synteny)", img: "layer_synteny.png", alt: "Conserved neighbours link orthologs between human and mouse",
  round: "Rounds 1 + 2", what: "Neighbouring genes stay together; the true ortholog shares them.",
  pros: ["Strong alone: AUC 0.92 human–mouse (shuffled 0.50)", "Hard cases +0.16 on development data", "Repair rule cut fission candidates by 29–40%"],
  cons: ["No signal at deep splits (human–fly 0.49)", "Sequence already solves the same cases", "Test: −0.1 F and −0.012 recall when removed"],
  verdict: "Kept as an option (default on); no measurable held-out gain", vc: ORANGE,
  notes: "Gene order may only re-rank partners inside the group the sequence-only pass chose (constrained mode); an unconstrained bonus cost 4-8 recall points. Round 2 tried gated synteny repair of fission errors: fission candidates fell from 706 to 426-499 with no change under shuffled order, but only 7 executed joins touched a labelled pair (5 correct), so its accuracy is unmeasured." });

layerSlide({ title: "3D structure (3Di)", img: "layer_structure.png", alt: "Two proteins with different sequences but the same 3Di structure string",
  round: "Rounds 1 + 2", what: "AlphaFold structures written as 3Di letters (Foldseek); ProstT5 fills gaps on the GPU.",
  pros: ["Family completeness 0.794 → 0.853, above the sequence ceiling", "ProstT5 on GPU: 40–60× faster than CPU"],
  cons: ["Paralogs share folds: no orthology signal", "OrthoBench test precision 69.5 → 60.1", "Full 3Di coverage added nothing more"],
  verdict: "Option, off by default", vc: ORANGE,
  notes: "3Di was tested in family building (rescue of fragments missed by sequence) and in merge scoring (deep AUC 0.694 vs 0.751 for sequence, rejected). Coverage: 94.5% of QfO proteins from AlphaFold DB/Foldcomp." });

layerSlide({ title: "Expression", img: "layer_expression.png", alt: "Tissue expression heatmap for an ortholog and a paralog",
  round: "Rounds 1 + 2", what: "Tissue expression of each gene copy (Bgee v15.2).",
  pros: ["Some signal overall (AUC 0.667)", "Open data (CC0)"],
  cons: ["Hard cases −0.019 [−0.20, +0.15]", "Specificity (τ) tie-break: −0.076", "Only 10 of 12 benchmark species covered"],
  verdict: "Rejected by the pre-registered rule", vc: RED });

layerSlide({ title: "Isoforms", img: "layer_isoforms.png", alt: "One gene with several alternative transcripts",
  round: "Round 1", what: "Search all splice isoforms, not one protein per gene.",
  pros: ["Could pair orthologous isoforms", "Data in every Ensembl release"],
  cons: ["Hard cases +0.013 [−0.018, +0.044]", "2–5× more search time"],
  verdict: "Rejected", vc: RED });

layerSlide({ title: "Genomic DNA: are lost genes really lost?", img: "layer_genomic.png", alt: "Searching the syntenic locus in the genome finds a pseudogene",
  round: "Round 1", what: "Search the genome at the syntenic locus for each apparent loss (miniprot).",
  pros: ["92% of planted controls found, 0 / 512 false hits", "53% of 2,523 apparent mammal losses were not real", "Exposed our own fission errors (32%)"],
  cons: ["Needs genomes and annotation", "Validates results, does not infer them"],
  verdict: "Kept as a quality-control add-on", vc: GREEN,
  notes: "25% were confirmed pseudogenes or deletions, for example the UOX pseudogene in human; GULO was confirmed by a separate clade-level search." });

layerSlide({ title: "Species tree: one tree for all families?", img: "layer_speciestree.png", w: 6.4, h: 3.46, alt: "NCBI, data-derived and random guide trees",
  round: "Round 1", what: "The guide tree decides the order of merges.",
  pros: ["A random tree costs 16.4 F: the tree does real work", "Small local gene trees correct it per family"],
  cons: ["Data-derived trees: +0.02–0.03 recall (n.s.), EC −0.12", "Long-branch errors, 25–30% conflict in bacteria"],
  verdict: "NCBI taxonomy + per-family local trees", vc: GREEN,
  notes: "Two further policies were planned but not run: delaying merges at uncertain nodes, and per-family guide trees." });

layerSlide({ title: "Interaction networks and co-fractionation", img: "layer_network.png", alt: "Two interaction networks with shared partners",
  round: "Round 2", what: "STRING v12 (own-species evidence only) and hu.MAP 3.0 complexes.",
  pros: ["Some signal over all pairs (AUC 0.63)", "Cheap: 0.4 CPU-hours"],
  cons: ["Covers 3.4% of candidates", "Hard cases at chance (AUC 0.50)", "Rescued pairs were orthologs only 48–54% of the time"],
  verdict: "Fails the first gate, mainly on coverage", vc: RED,
  notes: "Per-species co-fractionation data (Wan et al. 2015, McWhite et al. 2020) could not be obtained quickly and were not used; hu.MAP 3.0 is human only." });

layerSlide({ title: "Protein abundance", img: "layer_abundance.png", alt: "Abundance of an ortholog and paralog in proteomics data",
  round: "Round 2", what: "Whole-organism protein abundance from PaxDb 5.0.",
  pros: ["AUC 0.68 over all pairs, 0.69 on hard pairs"],
  cons: ["Right in only 39 of 116 hard decisions", "As a re-ranker: −0.131"],
  verdict: "Fails the first gate", vc: RED });

layerSlide({ title: "Domain architecture", img: "layer_domains.png", alt: "Domain architectures of orthologs and paralogs",
  round: "Round 2", what: "Pfam-A 38.2 domains and clans (TED/CATH not needed: Pfam covered >96%).",
  pros: ["Separates different groups (AUC 0.98)", "Best gate: TreeFam-A F1 +0.018 [+0.003, +0.034]"],
  cons: ["Inside a group no better than shuffled (0.54)", "OrthoBench −2.9 OG-F, precision −0.025", "48 CPU-hours of annotation"],
  verdict: "Fails the guard metrics", vc: RED });

layerSlide({ title: "Protein language models", img: "layer_plm.png", alt: "Embedding space where ortholog and paralog are equally close",
  round: "Round 2", what: "ESM-2 (35M–3B), ProtT5, and a head trained on our own synteny pairs.",
  pros: ["Above shuffled on all triplets (+0.23–0.33)", "Own head: +0.040 on all triplets", "Cheap on GPU: 3 GPU-hours"],
  cons: ["Below sequence (deep 0.61–0.72 vs 0.757)", "Fails where sequence fails; hard −0.007", "Neighbour rescue chained 229k of 251k genes"],
  verdict: "Fails the first gate", vc: RED,
  notes: "ProtT5 licence was unclear on the model hub; it was used locally only and no weights or embeddings were committed." });

layerSlide({ title: "Intron positions", img: "layer_introns.png", alt: "Intron positions mapped onto aligned proteins",
  round: "Round 2", what: "Shared intron positions and phases from Ensembl gene models.",
  pros: ["Some signal (AUC 0.578)", "Free with any genome annotation"],
  cons: ["Below sequence (0.654)", "Tie-break −0.071 [−0.149, +0.007]"],
  verdict: "Not integrated (a logged judgement call)", vc: RED });

layerSlide({ title: "GPU acceleration", img: "layer_gpu.png", alt: "Small families on the GPU, large families on DIAMOND",
  round: "Post-freeze", what: "An A100 GPU arrived after 19.5 h; all options are opt-in on branch v2-gpu.",
  pros: ["Own CUDA Smith–Waterman: 600–640 GCUPS", "1.18× faster end to end, same accuracy", "ProstT5 3Di 40–60× faster"],
  cons: ["All-GPU scoring changes the score graph: −5 to −8 OG-F", "Local trees are not GPU-bound", "MMseqs2-GPU slower at these sizes"],
  verdict: "Safe GPU option kept; faster variants not equivalent", vc: GREEN,
  notes: "78-85% of DIAMOND time goes to a few families above 1,000 genes. A rank-ratio sparsified all-GPU variant needed about 15x fewer CPU-hours on QfO with similar TreeFam-A dev accuracy, but lost 8.4 OrthoBench dev OG-F. A 1% jitter of CPU scores alone cost 3-7 OG-F, showing how sensitive the small benchmark is." });

// ================================================================ 32. summary matrix
{
  const s = slide("CONTENT");
  title(s, "Informative alone is not useful together");
  const G = "strong", W = "weak", Nn = "none", NA = "—";
  const colorOf = (v) => (v === G || v.startsWith("kept") || v.startsWith("QC") ? HEX.accent3 : v === W || v.startsWith("option") ? HEX.accent2 : v === Nn || v.startsWith("loss") || v.startsWith("rejected") || v.startsWith("fails") ? HEX.accent4 : HEX.dk1);
  const hdr = ["Layer", "Signal, all pairs", "Signal, hard cases", "Inside the pipeline", "Verdict"].map((t) => ({ text: t, options: { bold: true, color: HEX.dk2, fill: { color: HEX.lt2 } } }));
  const rows = [
    ["Gene order", G, "shallow only", "no test gain", "option, on"],
    ["3D structure", G, Nn, "loss of precision", "option, off"],
    ["Expression", W, Nn, "not run", "rejected"],
    ["Isoforms", NA, Nn, "not run", "rejected"],
    ["Genomic DNA", NA, NA, "finds errors", "QC add-on"],
    ["Networks", W, Nn, "not run", "fails"],
    ["Abundance", W, W, "not run", "fails"],
    ["Domains", G, W, "loss of precision", "fails"],
    ["Protein LMs", G, Nn, "not run", "fails"],
    ["Introns", W, W, "not run", "fails"],
  ];
  s.addTable([hdr, ...rows.map((r) => r.map((v, j) => ({ text: v, options: { color: j === 0 ? HEX.dk1 : colorOf(v), bold: j === 0 || j === 4 } })))],
    { x: 0.6, y: 1.35, w: 8.6, colW: [1.9, 1.6, 1.7, 1.9, 1.5], fontSize: 13, rowH: 0.45, border: { type: "solid", pt: 0.75, color: "D5DCE6" }, valign: "middle" });
  card(s, 9.6, 1.35, 3.1, 4.95, NAVY);
  txt(s, "Most layers tell orthologs from paralogs in general.\n\nFew help in the hard cases, the near-equal paralogs where sequence fails.\n\nThat is where a layer must add something new.",
    { x: 9.85, y: 1.6, w: 2.65, h: 4.5, fontSize: 16, color: WHITE, fontFace: "Cambria" });
}

// ================================================================ 33. section decisions
{
  const s = slide("SECTION", "Decision making");
  s.addText("Decision making", { placeholder: "title" });
  s.addText("How a layer earned its place, and what to reuse next time", { placeholder: "body" });
}

// ================================================================ 34. funnel
{
  const s = slide("CONTENT", "Decision making");
  title(s, "Every layer went through the same funnel");
  const steps = [
    ["Pre-register", "hypothesis, arms, keep rule", 12.1, BLUE],
    ["Signal alone", "AUC vs shuffled control", 10.4, BLUE],
    ["Hard cases", "near-equal paralogs", 8.7, ORANGE],
    ["Inside the pipeline (dev)", "lower CI > 0, guards hold", 7.0, ORANGE],
    ["Held-out test, once", "separate agent", 5.3, GREEN],
  ];
  steps.forEach(([h, b, w, c], i) => {
    const x = 0.6 + (12.1 - w) / 2, y = 1.4 + i * 0.98;
    card(s, x, y, w, 0.82, c);
    txt(s, h + "  ·  " + b, { x: x + 0.2, y, w: w - 0.4, h: 0.82, fontSize: 16, bold: false, color: WHITE, align: "center", valign: "middle" });
  });
  txt(s, "13 layers entered · 3 helped on development data (gene order, structure rescue, domain gate) · 0 improved held-out results",
    { x: 0.6, y: 6.35, w: 12.1, h: 0.45, fontSize: 15, color: NAVY, bold: true, align: "center" });
  s.addNotes("Gene order and the structure rescue reached the held-out test before the freeze; neither replicated. In round 2 no layer passed the keep rule, so nothing reached test. The domain gate passed its own metric but failed the guard metrics.");
}

// ================================================================ 35. here vs future
{
  const s = slide("CONTENT");
  title(s, "Decisions here, and how to carry them forward");
  const cols = [
    ["In this project", NAVY, ["Pre-registered keep rule + guard metrics", "Dev set small: ±2–4 F is noise", "Test scored once, after the freeze", "Every decision and reason in the research tree"]],
    ["Next time", GREEN, ["Build a hard-case benchmark first", "Same funnel, run by parallel agents", "Fresh, unseen benchmark for version 2", "Share the protocol with other projects"]],
  ];
  cols.forEach(([h, c, items], k) => {
    const x = 0.6 + k * 6.3;
    card(s, x, 1.4, 5.8, 0.7, c);
    txt(s, h, { x: x + 0.25, y: 1.4, w: 5.4, h: 0.7, fontSize: 20, bold: true, color: WHITE, valign: "middle", fontFace: "Cambria" });
    items.forEach((t, i) => {
      const y = 2.35 + i * 0.98;
      dot(s, x + 0.1, y + 0.15, 0.5, c, String(i + 1), 15);
      txt(s, t, { x: x + 0.8, y, w: 4.9, h: 0.8, fontSize: 17, valign: "middle" });
    });
  });
  s.addNotes("The protocol file docs/L2_protocol.md is general: worktrees per agent, flags default off, data rules, pre-registration, two-step test, resources and deadlines. Five agents ran eight layers in about 14 hours with it.");
}

// ================================================================ 36. AI-led research
{
  const s = slide("CONTENT", "Conclusion");
  title(s, "AI-led research: what worked and what did not");
  const good = ["~30 approaches tried and dropped, all documented", "Parallel agents: 5 layers explored overnight", "Independent checks caught real errors", "Full record in git: tree, claims, deviations"];
  const bad = ["Tuning pressure on a small dev set", "Dev gains that did not replicate on test", "Overstated claims needed correction", "No human verification yet"];
  [[good, "Worked", GREEN, "✓"], [bad, "Did not", RED, "✗"]].forEach(([items, h, c, m], k) => {
    const x = 0.6 + k * 6.3;
    txt(s, h, { x, y: 1.4, w: 5.8, h: 0.5, fontSize: 22, bold: true, color: c, fontFace: "Cambria" });
    items.forEach((t, i) => {
      const y = 2.05 + i * 1.0;
      card(s, x, y, 5.8, 0.85, c, { t: 85 });
      dot(s, x + 0.15, y + 0.16, 0.52, c, m, 16);
      txt(s, t, { x: x + 0.85, y, w: 4.8, h: 0.85, fontSize: 16, valign: "middle" });
    });
  });
  txt(s, "Examples of caught errors: omitted pre-registered comparators, too-narrow intervals, a preset mismatch in the scaling claims.",
    { x: 0.6, y: 6.25, w: 12.1, h: 0.45, fontSize: 13, color: GREY, italic: true });
  s.addNotes("Other friction: an API session limit stopped all agents at T+11.7 h (jobs survived); two QfO runs were lost to script edits during a run; the project was paused once to move to a GPU node.");
}

// ================================================================ 37. conclusion
{
  const s = slide("TITLE_DARK", "Conclusion");
  s.addText("Conclusion", { placeholder: "title" });
  txt(s, "Full record: REPORT.md · technical_report/ · minimal repository with code and predictions",
    { x: 0.8, y: 6.7, w: 11.7, h: 0.35, fontSize: 13, color: LIGHT, italic: true });
  const cols = [
    ["Outcome", "A working HOG method in 48 h: strong on OrthoBench, below the best on QfO recall, runs on 1,000 proteomes."],
    ["AI-led research", "Agents ran the whole cycle, from literature to report, with 7 human instructions. Verification by other agents was essential."],
    ["Next", "Outgroup tests at deep nodes · linear family step and streamed scores · hard-case benchmark · human verification · reuse in phase 3."],
  ];
  cols.forEach(([h, b], i) => {
    const x = 0.8 + i * 4.0;
    txt(s, h, { x, y: 3.95, w: 3.7, h: 0.5, fontSize: 24, bold: true, color: ORANGE, fontFace: "Cambria" });
    txt(s, b, { x, y: 4.55, w: 3.6, h: 2.0, fontSize: 18, color: WHITE });
  });
}

pres.writeFile({ fileName: OUT }).then(async () => {
  const { applyTheme } = require(process.env.PPTX_SKILL + "/scripts/apply_theme.js");
  await applyTheme(OUT, THEME);
  console.log("wrote", OUT);
});
