# Hackathon slide deck

`ancestra_hackathon_slides.pptx`: 37 slides documenting the whole experiment (hackathon context, prompt, decision path,
method, results, the 13 extra data layers one by one, decision making, AI-led research). Agent-generated; not yet
verified by a human. Blanks marked "[to be added]" are information the agent had no access to (dates, participants).

Rebuild from the repository root:

```bash
python technical_report/slides/make_assets.py                    # schematics and plots -> assets/ (matplotlib)
NODE_PATH=<dir containing pptxgenjs> node technical_report/slides/build_deck.js   # deck (pptxgenjs)
```

`build_deck.js` expects `PPTX_SKILL` to point to a directory with `scripts/apply_theme.js`, which writes the theme
colours into the file; without it the deck still builds, but the colour step fails. All numbers are copied from
`REPORT.md`. `assets/hackathon_flow.png` is the organisers' project-flow figure.
