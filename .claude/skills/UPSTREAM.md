# Vendored skills: upstream record

| item | value |
|---|---|
| upstream | https://github.com/BootLoops-ai/skills |
| commit | ca892277dcf0468d995f0036f3bd6d753a8afe7d (HEAD of `main` when fetched) |
| fetched | 2026-11-20, `git clone` + `git checkout` of the commit above |
| copied | `skills/planted-truth/`, `skills/acceptance-gate/`, `skills/independence-bookkeeping/`, `skills/timing-discipline/` |
| local paths | `.claude/skills/<name>/` (verbatim; `diff -r` against the clone empty) |
| licences | skill texts (`SKILL.md`): CC BY 4.0, `LICENSE-bootloops-CC-BY-4.0` (upstream `LICENSE-CONTENT`); scripts and manifests: MIT, `LICENSE-bootloops-MIT` (upstream `LICENSE`). The four skills copied contain no scripts. |

SHA-256 of the copied files:

```
72e4895a2fa148a31eb25ad977c6505e2eb77b55387d3dfc70f5e8e06da7ebdc  acceptance-gate/SKILL.md
e655aaf792063af7f51031d94b209c0157f72ce3ca7f6497e61bd36b06bd3909  independence-bookkeeping/SKILL.md
9ad58a3c5bbf0b1766968cc7deeeca1911b7f0c5d6383b564273b4692900b5a7  planted-truth/SKILL.md
7e0ae5519c466773940607ee6f018a039032c49174de7bdfa1efc11524db5da5  timing-discipline/SKILL.md
```

Do not edit these files. To update, fetch a new upstream commit, copy verbatim, and update this record and `NOTICE`.
