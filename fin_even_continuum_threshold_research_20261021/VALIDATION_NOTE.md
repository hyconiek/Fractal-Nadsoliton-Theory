# Validation note

Date: 2026-10-21.

The first metadata run failed before writing RESULTS/INPUTS/MANIFEST.
Its Markdown-link scanner interpreted the indented mathematical formula
with a square-bracket prefactor followed by parentheses as a local link
to `1+k^2`. This was not a scientific algebra failure or a broken prose link.

The scanner was corrected to exclude CommonMark indented code blocks,
in addition to fenced blocks and inline code, before checking prose links.
No formula, input, predecessor file or scientific condition was changed.
The exact algebra checks had already passed; the corrected full metadata
run and subsequent read-only replay are recorded in the final results.
