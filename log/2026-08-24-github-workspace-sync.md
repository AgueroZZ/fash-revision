# GitHub workspace synchronization — 2026-08-24

## Scope

Archive the current working-tree changes in `master`, including simulation and
analysis source code, workflowr analysis pages, rendered documentation, and
their generated figure assets.

## Notes

- This is a repository synchronization commit; it does not rerun analyses or
  alter the scientific results.
- `git diff --check` reports pre-existing and generated trailing whitespace in
  several R Markdown and HTML artifacts. The files are retained unchanged so
  this commit preserves the current rendered outputs.
