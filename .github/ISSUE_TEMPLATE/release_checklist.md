---
name: Release Checklist
about: Checklist for publishing a new release
title: "Release v"
labels: release
assignees: ''
---

## Release Checklist

**Target version:** vX.Y.Z

### Pre-release
- [ ] Version number updated in `CMakeLists.txt` (`project(FEM_2d VERSION X.Y.Z)`)
- [ ] `CHANGELOG.md` updated with new version section
- [ ] `RELEASE_NOTES.md` updated
- [ ] CI passes on all platforms (Linux, macOS, Windows)
- [ ] Documentation updated if needed

### Release
- [ ] Tag created: `git tag -a vX.Y.Z -m "Release vX.Y.Z"`
- [ ] Tag pushed: `git push origin vX.Y.Z`
- [ ] Release workflow completed successfully
- [ ] GitHub Release page shows correct assets (Linux/macOS/Windows binaries)

### Post-release
- [ ] Verify download links work
- [ ] Announce release (if applicable)
