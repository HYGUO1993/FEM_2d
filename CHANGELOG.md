# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [1.0.1] - 2026-08-11

### Fixed
- Removed accidentally committed `build/` directory containing hardcoded Linux paths, which caused macOS/Windows CI build failures
- Added `.gitignore` to prevent build artifacts from being committed again

### Added
- `--version` command-line flag for displaying program version
- `CHANGELOG.md` for tracking changes
- Release checklist issue template
- `workflow_dispatch` trigger for release workflow (allows manual testing)

## [1.0.0] - 2026-04-18

### Added
- Initial release
- 2D truss and frame structural analysis (FEM core in C++)
- Support for multiple material and cross-section types
- Nodal concentrated force loading
- Displacement, internal force, and reaction force computation
- Python visualization scripts (deformation, internal forces, reactions)
- Basic GUI framework (PySide6/PyQt6)
- Unit tests with CTest
- CI workflow for Linux

[1.0.1]: https://github.com/HYGUO1993/FEM_2d/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/HYGUO1993/FEM_2d/releases/tag/v1.0.0
